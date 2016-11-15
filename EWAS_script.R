#  Copyright (C) 2014-2016 Roby Joehanes

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, a copy is available at
#  https://www.gnu.org/licenses/gpl-3.0.en.html

# Data loading
# Assumption: rows = probes, cols = num samples
# The file is assumed to be big, so let's use fread to speed up the loading
library(data.table);
mdata <- data.frame(fread("my-methylation-data.txt"), check.names=FALSE, stringsAsFactors=FALSE);
rownames(mdata) <- mdata[,1];
mdata <- mdata[,-1];

# Phenotype data loading
# Assumption: comma delimited, rows = num samples, cols = num phenotypes, has a header with column names
pdata <- read.csv("my-phenotype-data.csv", as.is=TRUE, header=T);

#Make sure chip, row, and column effects are factors
pdata[,"ChipNo"] <- as.factor(pdata[,"ChipNo"]);
pdata[,"RowNo"] <- as.factor(pdata[,"RowNo"]);
pdata[,"ColNo"] <- as.factor(pdata[,"ColNo"]);

# IMPORTANT: pdata file and mdata file must be in the same order (by ID) before proceeding!

#Formula (you can customize this formula as long as the phenotype of interest is placed first)
fm <- y ~ cIMT + Age + Sex + WBC + LymPct + MonPct + EosPct + BasPct + PLT + (1|ChipNo) + (1|RowNo) + (1|ColNo);

#Fixed effect version of the formula above
fm_fixed <- y ~ cIMT + Age + Sex + WBC + LymPct + MonPct + EosPct + BasPct + PLT + ChipNo + RowNo + ColNo;

# Compute DFE (useful for p-value computation) using R's lm
# Since we only cares about the DFE, we can use any dummy values for y
pdata[,"y"] <- rnorm(NROW(pdata));
DFE <- lm(fm_fixed, data=pdata)$df.residual;

method = "LMM"; # Set to RLM if you want robust linear regression.

if (method == "LMM") {
	# LMM
	library(lme4);
	doOne <- function(i) {
		pdata$y <- as.numeric(mdata[i,]);
		result <- lmer(fm, data=pdata);
		ssr <- as.vector(getME(result, "RX") %*% getME(result, "beta"));
		reduced_y <- result@frame[, "y"];
		sst <- var(reduced_y) * (length(reduced_y) - 1);
		rsq <- ssr / sst;
		tbl <- summary(result);
		# Automatic detection between new and old lme4
		if (any(names(attributes(summary(tbl))) == "coefs")) {
			tbl <- tbl@coefs;
		} else {
			tbl <- tbl$coef;
		}
		# Returns P-value, Effect size, Standard Error, T-statistics, R^2
		result <- c(2*pt(abs(tbl[2,3]), df=DFE, lower.tail=FALSE), tbl[2,1], tbl[2,2], tbl[2,3], rsq[2]);
		return (result);
	}
} else if (method == "RLM") {
	# RLM
	library(robust);
	# Set iteration number to 1000 to overcome convergence problem.
	ctl <- lmRob.control(mxr=1000,mxf=1000,mxs=1000);
	doOne <- function(i) {
		pdata$y <- as.numeric(mdata[i,]);
		result <- lmRob(fm_fixed, data=pdata, x=TRUE, y=TRUE, control=ctl);
		sst <- var(result$y) * (length(result$y) - 1);
		tbl <- summary(result)$coef;
		ssr <- sum(as.vector(result$x[,2]) * tbl[2,1]^2);
		rsq <- ssr / sst;
		# Returns P-value, Effect size, Standard Error, T-statistics, R^2
		result <- c(tbl[2,4], tbl[2,1], tbl[2,2], tbl[2,3], rsq);
		return (result);
	}
} else {
	stop("Unknown method!");
}

library(parallel);

# NOTE: If you have a big machine with big RAM, use the parallel version.
# If not, use the serial version.

# Show how many cores your machine has:
num_cores <- detectCores(logical=TRUE);
print(num_cores);

# Show the size of your dataset in Mb
(used_mem <- gc(reset=TRUE)[2,2]);

# Make sure you have enough free RAM. The line below utilizes all the cores your machine has.
# However, you need to make sure that you have (num_cores * dataset) free RAM or else
# the parallel code will quit with out of memory error!
# Example: If you have 24 cores and 9 GB data size, make sure you have 216GB of free RAM.
# If you don't have that much free RAM, you can opt for fewer cores to run, e.g.: options(mc.cores = 10);

# Parallel version
options(mc.cores = num_cores);
result_all <- do.call(rbind, mclapply(1:NROW(mdata), doOne));

# Serial version. Uncomment the code below if you want to use it.

#result_all <- do.call(rbind, lapply(1:NROW(mdata), doOne));

rownames(result_all) <- rownames(mdata);
colnames(result_all) <- c("P", "Fx", "SE", "T", "RSq");

result_all <- as.data.frame(result_all)

write.csv(result_all, "outputfile.txt", row.names=TRUE, quote=FALSE, na="");


