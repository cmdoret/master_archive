# This script performs visualisations to make sense out of the multiQC and samtools stats output.

# Set to CL arg 

in_folder <- "data/wgs_wild/qc_output/multiqc_data/"

sam <- read.table(paste0(in_folder, "multiqc_samtools_stats.txt"), header=T)

clean_sam <- sam %>% 
  select_if(is.numeric) %>%
  select(-is_sorted) %>%
  select_if(function(x) sd(x) > 0) %>%
  select(contains("percent")) %>%
  summarise_all(scale)

clean_names <- str_extract(sam$Sample, "[^_]*")

pairs(clean_sam, col=1:15, pch=".")
cormat <- cor(t(clean_sam))

heatmap(cormat, symm = T, labRow = clean_names, labCol = clean_names)

cormat[upper.tri(cormat)] <- NA
heatmap(cormat, labRow = clean_names, labCol = clean_names, Rowv = NA, Colv = NA, verbose = T)

PC <- prcomp(t(clean_sam))
plot(PC$rotation[,1:2], pch=" ");text(PC$rotation[,1:2], labels = clean_names)

pairs(PC$rotation[,1:4], col=1:15, pch=)
PCpos <- apply(PC$rotation, MARGIN = 2, function(x) x-min(x))
#scatterplot3d(PC$rotation[,1:3], main="3D Scatterplot",highlight.3d = T, type = 'h')
plot3d(PCpos[,1:3], col="red", size=3, type='h', pch='.') ;text3d(PCpos[,1:3],  texts = clean_names, add = T, 
       colkey = FALSE, cex = 1)       


# PCA stuff


res.pca <- PCA(clean_sam, graph = FALSE)

# Proportion of variance explained by each PC
fviz_screeplot(res.pca, ncp=10)


# Which variables correlate together ? which are best described by the 2 first PC
# The sum of the cos2 for variables on the principal components is equal to one.
# If a variable is perfectly represented by only two components, the sum of the 
# cos2 is equal to one. In this case the variables will be positioned on the circle of correlations.

fviz_pca_var(res.pca, col.var="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.5) + theme_minimal()

# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",axes =c(1,2)) +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=5) + theme_minimal()

# Contribution of each variable to variation between species on PC1 and PC2
fviz_pca_contrib(res.pca, choice = "var", axes = 1)
fviz_pca_contrib(res.pca, choice = "var", axes = 2)
fviz_pca_contrib(res.pca, choice = "var", axes = 3)
fviz_pca_contrib(res.pca, choice = "var", axes = 1:2) # Contribution on PC 1 and 2 together
fviz_pca_ind(res.pca,axes=c(1,2),geom="point")
fviz_pca_ind(res.pca,axes=c(2,3),geom="point")



