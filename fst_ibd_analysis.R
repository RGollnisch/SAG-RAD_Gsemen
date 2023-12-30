# Script to plot Fst matrix and isolation by distance (IBD) analysis

library(dplyr)
library(ggplot2)
library(ggnewscale)
library(plotly)
library(ade4)
library(vegan)
library(ggpubr)


# Function to order distance matrix
reord_distmat <- function(distmat){
  dd <- as.dist(distmat)
  hc <- hclust(dd, method = "complete")
  neworder <- dimnames(distmat)[[1]][hc$order]
  distmat_reord <- distmat[neworder, neworder]
  return(distmat_reord)
}


## Read geographic distance data from genodive
eur_geo <- read.csv("./analysis/10.genodive/eur_50k_p13_r70_genodive_geo.csv", row.names = 1)
eur_geo <- eur_geo[!dimnames(eur_geo)[[1]] %in% c("RM"), !dimnames(eur_geo)[[2]] %in%  c("RM")]
dimnames(eur_geo)[[1]] <- c("EE/VP", "EE/PS", "LT/NT", "LT/PB", "CZ/KO", "PL/SE", "DE/RO", "DE/GF", "NL/BR", "NL/SI", "DK/HU", "SE/HE", "PT/PE", "PT/PI")
dimnames(eur_geo)[[2]] <- c("EE/VP", "EE/PS", "LT/NT", "LT/PB", "CZ/KO", "PL/SE", "DE/RO", "DE/GF", "NL/BR", "NL/SI", "DK/HU", "SE/HE", "PT/PE", "PT/PI")

usa_geo <- read.csv("./analysis/10.genodive/usa_50k_noSMLL_p7_r70_genodive_geo.csv", row.names = 1)
usa_geo <- usa_geo[!dimnames(usa_geo)[[1]] %in% c("CC"), !dimnames(usa_geo)[[2]] %in%  c("CC")]
dimnames(usa_geo)[[1]] <- c("WA/VO", "MI/RE", "ME/NO", "ME/WN", "ME/PQ", "MA/PT", "MA/CS")
dimnames(usa_geo)[[2]] <- c("WA/VO", "MI/RE", "ME/NO", "ME/WN", "ME/PQ", "MA/PT", "MA/CS")

## Read genetic distance (Fst) data from 
eur_fst <- read.csv("./analysis/8.populations/EUR_50k_p13_r70/populations.fst_summary.tsv", sep = "\t", row.names = 1)
eur_fst <- rbind(eur_fst, RM = c(rep(NA, 14), 0))
diag(eur_fst) <- 0
eur_fst[lower.tri(eur_fst)]  <- t(eur_fst)[lower.tri(eur_fst)]
eur_fst <- eur_fst[!dimnames(eur_fst)[[1]] %in% c("RM"), !dimnames(eur_fst)[[2]] %in%  c("RM")]
dimnames(eur_fst)[[1]] <- c("EE/VP", "LT/NT", "DE/GF", "NL/BR", "NL/SI", "LT/PB", "DE/RO", "EE/PS", "PT/PE", "PL/SE", "PT/PI", "SE/HE", "DK/HU", "CZ/KO")
dimnames(eur_fst)[[2]] <- c("EE/VP", "LT/NT", "DE/GF", "NL/BR", "NL/SI", "LT/PB", "DE/RO", "EE/PS", "PT/PE", "PL/SE", "PT/PI", "SE/HE", "DK/HU", "CZ/KO")

usa_fst <- read.csv("./analysis/8.populations/USA_noSMLLCC_50k_p7_r70/populations.fst_summary.tsv", sep = "\t", row.names = 1)
usa_fst <- rbind(usa_fst, NO = c(rep(NA, 6), 0))
diag(usa_fst) <- 0
usa_fst[lower.tri(usa_fst)]  <- t(usa_fst)[lower.tri(usa_fst)]
dimnames(usa_fst)[[1]] <- c("MA/PT", "MI/RE", "ME/PQ", "WA/VO", "ME/WN", "MA/CS", "ME/NO")
dimnames(usa_fst)[[2]] <- c("MA/PT", "MI/RE", "ME/PQ", "WA/VO", "ME/WN", "MA/CS", "ME/NO")


## Convert to matrix and arrange order
eur_fst_matrix <- as.matrix(eur_fst)
dimnames(eur_fst_matrix)[[1]]
eur_fst_matrix <- reord_distmat(eur_fst_matrix)
eur_order <- dimnames(eur_fst_matrix)[[1]]
dimnames(eur_fst_matrix)[[1]]

hc <- hclust(as.dist(eur_fst_matrix), method = "complete")
hc <- hclust(as.dist(eur_fst_matrix), method = "median")
plot(hc, ylab = "Fst")
hcd <- as.dendrogram(hc)

pdf(file = "./analysis/manuscript_figures/eur_fst_tree.pdf", height = 5, width = 10)
  plot(hcd, type = "triangle", ylab = "Fst", )
dev.off()

eur_geo_matrix <- as.matrix(eur_geo)
eur_geo_matrix <- eur_geo_matrix[eur_order, eur_order]


usa_fst_matrix <- as.matrix(usa_fst)
dimnames(usa_fst_matrix)[[1]]
usa_fst_matrix <- reord_distmat(usa_fst_matrix)
usa_order <- dimnames(usa_fst_matrix)[[1]]
dimnames(usa_fst_matrix)[[1]]

hc <- hclust(as.dist(usa_fst_matrix), method = "complete")
hcd <- as.dendrogram(hc)

pdf(file = "./analysis/manuscript_figures/usa_fst_tree.pdf", height = 5, width = 10)
  plot(hcd, type = "triangle", ylab = "Fst")
dev.off()

usa_geo_matrix <- as.matrix(usa_geo)
usa_geo_matrix <- usa_geo_matrix[usa_order, usa_order]


## Convert geo dist/matrix object to data frame
eur_geo_ind = which(upper.tri(eur_geo_matrix, diag = TRUE), arr.ind = TRUE)
usa_geo_ind = which(upper.tri(usa_geo_matrix, diag = TRUE), arr.ind = TRUE)


eur_geo.df = data.frame(loc1_f = factor(dimnames(eur_geo_matrix)[[2]][eur_geo_ind[,2]], levels = eur_order),
                    loc2_f = factor(dimnames(eur_geo_matrix)[[1]][eur_geo_ind[,1]], levels = eur_order),
                    dist = eur_geo_matrix[eur_geo_ind] %>% round(digits = 0))

usa_geo.df = data.frame(loc1_f = factor(dimnames(usa_geo_matrix)[[2]][usa_geo_ind[,2]], levels = usa_order),
                        loc2_f = factor(dimnames(usa_geo_matrix)[[1]][usa_geo_ind[,1]], levels = usa_order),
                        dist = usa_geo_matrix[usa_geo_ind] %>% round(digits = 0))

eur_geo.df <- eur_geo.df %>%
  mutate(pair1 = paste(loc1_f, loc2_f, sep = "-"))
eur_geo.df <- eur_geo.df %>%
  mutate(pair2 = paste(loc2_f, loc1_f, sep = "-"))

usa_geo.df <- usa_geo.df %>%
  mutate(pair1 = paste(loc1_f, loc2_f, sep = "-"))
usa_geo.df <- usa_geo.df %>%
  mutate(pair2 = paste(loc2_f, loc1_f, sep = "-"))


## Convert Fst dist/matrix object to data frame
eur_fst_ind = which(upper.tri(eur_fst_matrix, diag = TRUE), arr.ind = TRUE)
usa_fst_ind = which(upper.tri(usa_fst_matrix, diag = TRUE), arr.ind = TRUE)

eur_fst.df = data.frame(loc1_f = factor(dimnames(eur_fst_matrix)[[2]][eur_fst_ind[,2]], levels = eur_order),
                        loc2_f = factor(dimnames(eur_fst_matrix)[[1]][eur_fst_ind[,1]], levels = eur_order),
                        fst = eur_fst_matrix[eur_fst_ind] %>% round(digits = 3))

usa_fst.df = data.frame(loc1_f = factor(dimnames(usa_fst_matrix)[[2]][usa_fst_ind[,2]], levels = usa_order),
                        loc2_f = factor(dimnames(usa_fst_matrix)[[1]][usa_fst_ind[,1]], levels = usa_order),
                        fst = usa_fst_matrix[usa_fst_ind] %>% round(digits = 3))


eur_fst.df <- eur_fst.df %>%
  mutate(pair = paste(loc1_f, loc2_f, sep = "-"))
usa_fst.df <- usa_fst.df %>%
  mutate(pair = paste(loc1_f, loc2_f, sep = "-"))

# Convert zero values to NA
eur_fst.df$fst[eur_fst.df$fst == 0] = NA
usa_fst.df$fst[usa_fst.df$fst == 0] = NA

eur_geo.df$dist[eur_geo.df$dist == 0] = NA
usa_geo.df$dist[usa_geo.df$dist == 0] = NA

# Convert negative values to zero
eur_fst.df$fst[eur_fst.df$fst < 0] = 0
usa_fst.df$fst[usa_fst.df$fst < 0] = 0

# Fst italic label
fst.label = expression(italic("F")[ST])
geo.label = expression("d (km)")


## plot EUR Fst and geo dist heatmap
eur_heatmap <- ggplot(data = cbind(eur_fst.df[1:3], eur_geo.df[3])) +
  geom_tile(aes(x = loc1_f, y = loc2_f, fill = fst), color = "black") +
  geom_text(aes(x = loc1_f, y = loc2_f, label = fst), color = "black", size = 3, fontface = "bold") +
  scale_fill_viridis_c(option = "D", alpha = 0.75, begin = 1, end = 0, na.value = "white", name = fst.label, limits = c(min(eur_fst.df$fst), max(eur_fst.df$fst))) +
  new_scale_fill() +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0), position = "right", limits = rev(levels(eur_fst.df$loc2_f))) +
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

eur_heatmap
ggsave("./analysis/manuscript_figures/figure_fst-eur.pdf", device = cairo_pdf, width = 210, height = 105, units = "mm")
ggsave("./analysis/manuscript_figures/figure_fst-eur.png", width = 210, height = 105, units = "mm", dpi = 300)


# plot USA Fst heatmap
usa_heatmap <- ggplot(data = cbind(usa_fst.df[1:3], usa_geo.df[3])) +
  geom_tile(aes(x = loc1_f, y = loc2_f, fill = fst), color = "black") +
  geom_text(aes(x = loc1_f, y = loc2_f, label = fst), color = "black", size = 3, fontface = "bold") +
  scale_fill_viridis_c(option = "D", alpha = 0.75, begin = 1, end = 0, na.value = "white", name = fst.label, limits = c(min(usa_fst.df$fst), max(usa_fst.df$fst))) +
  new_scale_fill() +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0), position = "right", limits = rev(levels(usa_fst.df$loc2_f))) +
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

usa_heatmap
ggsave("./analysis/manuscript_figures/figure_fst-usa.pdf", device = cairo_pdf, width = 210, height = 105, units = "mm")
ggsave("./analysis/manuscript_figures/figure_fst-usa.png", width = 210, height = 105, units = "mm", dpi = 300)


# Combine both possible pair combinations for IBD plot
eur_combined1 <- inner_join(eur_geo.df, eur_fst.df[,3:4], by = c("pair1" = "pair"), keep = TRUE)
eur_combined2 <- inner_join(eur_geo.df, eur_fst.df[,3:4], by = c("pair2" = "pair"), keep = TRUE)

usa_combined1 <- inner_join(usa_geo.df, usa_fst.df[,3:4], by = c("pair1" = "pair"), keep = TRUE)
usa_combined2 <- inner_join(usa_geo.df, usa_fst.df[,3:4], by = c("pair2" = "pair"), keep = TRUE)

eur_ibd <- rbind(eur_combined1, eur_combined2)
usa_ibd <- rbind(usa_combined1, usa_combined2)

eur_ibd$fst[eur_ibd$fst == 0] = NA
eur_ibd$dist[eur_ibd$dist == 0] = NA

usa_ibd$fst[usa_ibd$fst == 0] = NA
usa_ibd$dist[usa_ibd$dist == 0] = NA


## Plot IBD in Europe
## RM excluded
eur_ibd_noRM <- filter(eur_ibd, !loc1_f %in% c("CZ/RM") & !loc2_f %in% c("CZ/RM"))
eur_ibd_east <- filter(eur_ibd, loc1_f %in% c("EE/VP", "EE/PS", "SE/HE", "LT/NT", "LT/PB", "PL/SE", "DE/GF") & loc2_f %in% c("EE/VP", "EE/PS", "SE/HE", "LT/NT", "LT/PB", "PL/SE", "DE/GF"))
eur_ibd_west <- filter(eur_ibd, loc1_f %in% c("DK/HU", "NL/SI", "NL/BR", "DE/RO", "CZ/KO", "PT/PE", "PT/PI") & loc2_f %in% c("DK/HU", "NL/SI", "NL/BR", "DE/RO", "CZ/KO", "PT/PE", "PT/PI"))

eur_ibd.plot <- ggplot(mapping = aes(x = dist, y = fst)) +
  geom_point(data = eur_ibd_noRM, size = 3, alpha = 0.75) +
  geom_point(data = eur_ibd_east, size = 4, color = "firebrick3") +
  geom_smooth(data = eur_ibd_east, method = "lm", se = FALSE, color = "firebrick3") +
  geom_point(data = eur_ibd_west, size = 4, color = "royalblue3") +
  geom_smooth(data = eur_ibd_west, method = "lm", se = FALSE, color = "royalblue3") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x  = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y  = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
eur_ibd.plot

ggplotly(eur_ibd.plot)

eur_ibd.plot <- eur_ibd.plot +
  labs(x = "Geographic distance (km)", y = fst.label)

eur_ibd.plot
ggsave("./analysis/manuscript_figures/eur_ibd.pdf", device = cairo_pdf, width = 210, height = 105, units = "mm")
ggsave("./analysis/manuscript_figures/eur_ibd.png", width = 210, height = 105, units = "mm", dpi = 300)


lm.eur_ibd_noRM <- lm(fst ~ dist, data = eur_ibd_noRM)
summary(lm.eur_ibd_noRM)

lm.eur_ibd_east <- lm(fst ~ dist, data = eur_ibd_east)
summary(lm.eur_ibd_east)

lm.eur_ibd_west <- lm(fst ~ dist, data = eur_ibd_west)
summary(lm.eur_ibd_west)


eur_geo_matrix
as.dist(eur_geo_matrix)
eur_fst_matrix
as.dist(eur_fst_matrix)

# arrange matrices noRM, east, and west
eur_fst_matrix_noRM <- eur_fst_matrix[!dimnames(eur_fst_matrix)[[1]] %in% c("CZ/RM"),
                                      !dimnames(eur_fst_matrix)[[2]] %in%  c("CZ/RM")]
eur_geo_matrix_noRM <- eur_geo_matrix[!dimnames(eur_geo_matrix)[[1]] %in% c("CZ/RM"),
                                      !dimnames(eur_geo_matrix)[[2]] %in%  c("CZ/RM")]
eur_fst_matrix_east <- eur_fst_matrix[dimnames(eur_fst_matrix)[[1]] %in% c("EE/VP", "EE/PS", "SE/HE", "LT/NT", "LT/PB", "PL/SE", "DE/GF"),
                                      dimnames(eur_fst_matrix)[[2]] %in% c("EE/VP", "EE/PS", "SE/HE", "LT/NT", "LT/PB", "PL/SE", "DE/GF")]
eur_geo_matrix_east <- eur_geo_matrix[dimnames(eur_geo_matrix)[[1]] %in% c("EE/VP", "EE/PS", "SE/HE", "LT/NT", "LT/PB", "PL/SE", "DE/GF"),
                                      dimnames(eur_geo_matrix)[[2]] %in% c("EE/VP", "EE/PS", "SE/HE", "LT/NT", "LT/PB", "PL/SE", "DE/GF")]
eur_fst_matrix_west <- eur_fst_matrix[dimnames(eur_fst_matrix)[[1]] %in% c("DK/HU", "NL/SI", "NL/BR", "DE/RO", "CZ/KO", "PT/PE", "PT/PI"),
                                      dimnames(eur_fst_matrix)[[2]] %in% c("DK/HU", "NL/SI", "NL/BR", "DE/RO", "CZ/KO", "PT/PE", "PT/PI")]
eur_geo_matrix_west <- eur_geo_matrix[dimnames(eur_geo_matrix)[[1]] %in% c("DK/HU", "NL/SI", "NL/BR", "DE/RO", "CZ/KO", "PT/PE", "PT/PI"),
                                      dimnames(eur_geo_matrix)[[2]] %in% c("DK/HU", "NL/SI", "NL/BR", "DE/RO", "CZ/KO", "PT/PE", "PT/PI")]

eur_geo_matrix


# Mantel test
vegan::mantel(eur_fst_matrix_noRM, eur_geo_matrix_noRM, permutations = 999, na.rm = TRUE)
ibd_noRM <- mantel.randtest(as.dist(eur_fst_matrix_noRM), as.dist(eur_geo_matrix_noRM), nrepet = 999)
ibd_noRM

vegan::mantel(eur_fst_matrix_east, eur_geo_matrix_east, permutations = 999, na.rm = TRUE)
ibd_east <- mantel.randtest(as.dist(eur_fst_matrix_east), as.dist(eur_geo_matrix_east), nrepet = 999)
ibd_east

vegan::mantel(eur_fst_matrix_west, eur_geo_matrix_west, permutations = 999, na.rm = TRUE)
ibd_west <- mantel.randtest(as.dist(eur_fst_matrix_west), as.dist(eur_geo_matrix_west), nrepet = 999)
ibd_west


## Plot IBD in United States
usa_ibd.plot <- ggplot(mapping = aes(x = dist, y = fst)) +
  geom_point(data = usa_ibd, size = 3, alpha = 0.75) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x  = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y  = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
usa_ibd.plot

ggplotly(usa_ibd.plot)

usa_ibd.plot <- usa_ibd.plot +
  labs(x = "Geographic distance (km)", y = fst.label)

usa_ibd.plot
ggsave("./analysis/manuscript_figures/usa_ibd.pdf", device = cairo_pdf, width = 210, height = 105, units = "mm")
ggsave("./analysis/manuscript_figures/usa_ibd.png", width = 210, height = 105, units = "mm", dpi = 300)


lm.usa_ibd <- lm(fst ~ dist, data = usa_ibd)
summary(lm.usa_ibd)

usa_geo_matrix
as.dist(usa_geo_matrix)
usa_fst_matrix
as.dist(usa_fst_matrix)

# Mantel test
vegan::mantel(usa_fst_matrix, usa_geo_matrix, permutations = 999, na.rm = TRUE)
ibd_usa <- mantel.randtest(as.dist(usa_fst_matrix), as.dist(usa_geo_matrix), nrepet = 999)
ibd_usa


ggarrange(eur_ibd.plot, usa_ibd.plot, ncol = 1, nrow = 2, labels = c("A", "B"))
ggsave("./analysis/manuscript_figures/figure6_ibd.pdf", device = cairo_pdf, width = 210, height = 210, units = "mm")
ggsave("./analysis/manuscript_figures/figure6_ibd.png", width = 210, height = 210, units = "mm", dpi = 300)
ggsave("./analysis/manuscript_figures/figure6_ibd.eps", device = cairo_ps, width = 210, height = 210, units = "mm")
