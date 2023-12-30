# Script to run discriminant analysis of principal components
# With environmental variables superimposed

library(adegenet)
library(vegan)
library(tidyverse)


obj1 <- read.genepop(file = "./analysis/8.populations/ALL_noSMLL_50k_p20_r60/")
obj2 <- read.genepop(file = "./analysis/8.populations/EUR_50k_p13_r70/populations_snps.gen")
obj3 <- read.genepop(file = "./analysis/8.populations/USA_noSMLL_50k_p7_r70/populations_snps.gen")


obj1_noRMCC <- obj1[!obj1$pop %in% c("RM-40", "CC-38")]
obj2_noRM <- obj2[obj2$pop != "RM-40"]
obj3_noCC <- obj3[obj3$pop != "CC-38"]


pop(obj1)
pop(obj1_noRMCC)
pop(obj2)
pop(obj2_noRM)
pop(obj3)
pop(obj3_noCC)

# DAPC cross validation

mat <- tab(obj1_noRMCC, NA.method = "mean")
grp <- pop(obj1_noRMCC)
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

mat <- tab(obj2_noRM, NA.method = "mean")
grp <- pop(obj2_noRM)
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

mat <- tab(obj3_noCC, NA.method = "mean")
grp <- pop(obj3_noCC)
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)


dapc.obj1_noRMCC <- dapc(obj1_noRMCC, var.contrib = TRUE, scale = FALSE, n.pca = 100, n.da = nPop(obj1_noRMCC) - 1)
temp <- optim.a.score(dapc.obj1_noRMCC)

dapc.obj1_noRMCC <- dapc(obj1_noRMCC, var.contrib = TRUE, scale = FALSE, n.pca = 14, n.da = nPop(obj1_noRMCC) - 1)
pdf(file = "./analysis/manuscript_figures/dapc_scatter_all_noRMCC.pdf")
  scatter(dapc.obj1_noRMCC, scree.da = FALSE, cell = 0, cstar = 1, mstree = FALSE, leg = TRUE,
          col = c("darkblue", "darkblue", "darkblue", "darkblue", "darkblue", "darkblue", "darkblue",
                  "darkblue", "darkblue", "darkblue", "darkblue", "darkblue", "darkblue", "darkblue",
                  "orange", "orange", "orange", "orange", "orange", "orange", "orange"),
          pch = 19, cex = 2, lwd = 2, lty = 2)
dev.off()



dapc.obj2_noRM <- dapc(obj2_noRM, var.contrib = TRUE, scale = FALSE, n.pca = 100, n.da = nPop(obj2_noRM) - 1)
temp <- optim.a.score(dapc.obj2_noRM)

dapc.obj2_noRM <- dapc(obj2_noRM, var.contrib = TRUE, scale = FALSE, n.pca = 10, n.da = nPop(obj2_noRM) - 1)
pdf(file = "./analysis/manuscript_figures/dapc_scatter_eur_noRM.pdf")
  scatter(dapc.obj2_noRM, scree.da = FALSE, cell = 0, cstar = 1, mstree = FALSE, leg = TRUE,
          col = c("darkblue", "blue", "red", "orange", "orange", "blue", "red", "darkblue", "darkgreen", "lightblue", "darkgreen", "purple", "brown", "black"),
          pch = 19, cex = 2, lwd = 2, lty = 2)
dev.off()


# Superimpose environment EUR

env_data <- read.csv("./info/SAG-RAD_gonyPop_lakedata_EUR.csv", header = TRUE)

selection <- c("abr","wc", "pH", "TP", "Ca", "Fe", "DOC", "TN")
env_analyse <- env_data[selection]

indivs <- read.csv("./info/populations_popmaps/popmap50k_eur", sep = "\t", header = FALSE)
colnames(indivs) <- c("ind", "abr")
indivs <- filter(indivs, !abr %in% c("RM"))
env_vars <- left_join(indivs, env_analyse, by = "abr")

env_vars_ord <- env_vars[match(dimnames(dapc.obj2_noRM$ind.coord)[[1]], env_vars$ind), ]

env <- envfit(dapc.obj2_noRM$ind.coord, env_vars_ord[-c(1, 2)], permutations = 999, na.rm = TRUE)

pdf(file = "./analysis/manuscript_figures/dapc_scatter_eur_noRM_ENV.pdf")
  scatter(dapc.obj2_noRM, scree.da = FALSE, cell = 0, cstar = 1, mstree = FALSE, leg = TRUE,
        col = c("darkblue", "blue", "red", "orange", "orange", "blue", "red", "darkblue", "darkgreen", "lightblue", "darkgreen", "purple", "brown", "black"),
        pch = 19, cex = 2, lwd = 2, lty = 2)
  plot(env, col = "black", lwd = 100)
dev.off()


dapc.obj3_noCC <- dapc(obj3_noCC, var.contrib = TRUE, scale = FALSE, n.pca = 100, n.da = nPop(obj3_noCC) - 1)
temp <- optim.a.score(dapc.obj3_noCC)


dapc.obj3_noCC <- dapc(obj3_noCC, var.contrib = TRUE, scale = FALSE, n.pca = 13, n.da = nPop(obj3_noCC) - 1)
pdf(file = "./analysis/manuscript_figures/dapc_scatter_usa_noCC.pdf")
  scatter(dapc.obj3_noCC, scree.da = FALSE, cell = 0, cstar = 1, mstree = FALSE, leg = TRUE,
          col = c("darkblue", "purple", "darkgreen", "orange", "darkgreen", "darkblue", "darkgreen"),
          pch = 19, cex = 2, lwd = 2, lty = 2)
  plot(env)
dev.off()


# Superimpose environment USA

env_data <- read.csv("./info/SAG-RAD_gonyPop_lakedata_USA.csv", header = TRUE)

selection <- c("abr", "wc", "pH", "TP", "Ca", "Fe", "DOC", "TN")
env_analyse <- env_data[selection]

indivs <- read.csv("./info/populations_popmaps/popmap50k_usa", sep = "\t", header = FALSE)
colnames(indivs) <- c("ind", "abr")
indivs <- filter(indivs, !abr %in% c("SM", "LL", "CC"))
env_vars <- left_join(indivs, env_analyse, by = "abr")

env_vars_ord <- env_vars[match(dimnames(dapc.obj3_noCC$ind.coord)[[1]], env_vars$ind), ]

env <- envfit(dapc.obj3_noCC$ind.coord, env_vars_ord[-c(1, 2)], permutations = 999, na.rm = TRUE)

pdf(file = "./analysis/manuscript_figures/dapc_scatter_usa_noCC_ENV.pdf")
  scatter(dapc.obj3_noCC, scree.da = FALSE, cell = 0, cstar = 1, mstree = FALSE, leg = TRUE,
        col = c("darkblue", "purple", "darkgreen", "orange", "darkgreen", "darkblue", "darkgreen"),
        pch = 19, cex = 2, lwd = 2, lty = 2)
  plot(env, col = "black")
dev.off()