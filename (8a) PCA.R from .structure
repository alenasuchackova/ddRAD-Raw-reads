#### PCA from SNP table###################################################################
#This script produces PCA from populations.structure file produced by stacks populations##
#and customizes the population colors#####################################################

###set your working directory in your local computer to the folder where you genotype structure file is
### download populations.structure file produced by stacks populations and change the suffix to .stru

library(adegenet)
library(ade4)
library(ggplot2)
library(dplyr)
library (svglite) # Save the plot as SVG

# Read the structure file and save it as a genind object called data
# we need the number of individuals and the number of loci (SNPs)= number of lines in populations.plink.map -1
##camilla
#data <- read.structure("populations.stru", n.ind=113, n.loc=77885, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=2, NA.char=0, ask=FALSE)
#rivularis
#data <- read.structure("populations.stru", n.ind=45, n.loc=3855, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=2, NA.char=0, ask=FALSE)
#paphia
data <- read.structure("populations.stru", n.ind=87, n.loc=143416, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=2, NA.char=0, ask=FALSE)


sample_names <- indNames(data)

# Fill missing data with average
X <- scaleGen(data, NA.method="mean")

# Perform PCA
pca1 <- dudi.pca(X, cent=FALSE, scale=FALSE, scannf=FALSE, nf=3)

# Estimate contribution of the first PC's to print them in the final plot
PC1 <- (100 * pca1$eig / sum(pca1$eig))[1]
PC2 <- (100 * pca1$eig / sum(pca1$eig))[2]

# Create a data frame for ggplot2
pca_df <- data.frame(PC1 = pca1$li[,1], PC2 = pca1$li[,2], Population = pop(data))
pca_df$Sample <- indNames(data)


# Define colors for populations
##rivularis
#pop_colors <- c("Czechia" = "#f08080", "China" = "green", "Hungary" = "blue", "Romania" = "darkred",
#                "Kazakhstan" = "#ffd700", "Kyrgyzstan" = "yellow", "Dagestan" = "black", #"FarEast" = "#228b22")
##camilla
#pop_colors <- c(Belgium = "#1f78b4", Bulgaria = "#33a02c", Caucasus = "#e31a1c", Czechia = "#6a3d9a", Estonia = "#ff7f00", FarEast = "#b15928", France = "#a6cee3", Germany = "#b2df8a", Hungary = "#fb9a99", Italy = "#fdbf6f", Japan = "#cab2d6", Kaliningrad = "#ffff99", Poland = "#8dd3c7", Rom_Sa = "#fb8072", Romania = "#80b1d3", Serbia = "#f4cae4", Slovakia = "#d9d9d9", Spain = "#bc80bd", Switzerland = "#ccebc5")
##paphia
pop_colors <- c(Algeria = "#f08080", Armenia = "green", Corsica = "blue", Czechia = "darkred", Czechia_M = "#ffd700", Dagestan = "yellow", England = "black", FarEast = "#228b22", France = "#bebada", Georgia = "#fb8072", Hungary = "#80b1d3", China = "#fdb462", China_South = "#b3de69", Italy_central = "#fccde5", Italy_south = "#d9d9d9", Japan_Hokk = "#bc80bd", Japan_Hon = "#ccebc5", Kazakhstan = "#ffed6f", Kyrgyzstan = "#d95f02", Poland = "#e7298a", Romania = "#66a61e", Sardinia = "#e6ab02", Sicily = "#a6761d", Slovenia = "#666666", Spain_central = "#1b9e77", Spain_north = "#d53e4f", Sweden = "#5e4fa2", Ukraine = "#f46d43", Ural = "#fee08b")




# Plot PCA using ggplot2 - with sample labels
p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Population)) +
  geom_point(size = 5, alpha = 0.7) +  # Increase point size
  scale_color_manual(values = pop_colors) +
  geom_text (size = 3) +
  labs(title = "PCA of Population Structure",
       x = paste("PC1 (", round(PC1, 2), "%)", sep=""),
       y = paste("PC2 (", round(PC2, 2), "%)", sep=""),
       color = "Locality") +
  theme_minimal() +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),  # Increase legend text size
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Increase title size
        axis.title = element_text(size = 16),  # Increase axis title size
        axis.text = element_text(size = 14),  # Increase axis text size
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # Use linewidth instead of size
        axis.line = element_line(linewidth = 0.5, color = "black"))  # Use linewidth instead of size



#plot PCA without labels
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 5, alpha = 0.7) +  # Increase point size
  scale_color_manual(values = pop_colors) +
  labs(title = "PCA of Population Structure",
       x = paste("PC1 (", round(PC1, 2), "%)", sep=""),
       y = paste("PC2 (", round(PC2, 2), "%)", sep=""),
       color = "Locality") +
  theme_minimal() +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),  # Increase legend text size
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Increase title size
        axis.title = element_text(size = 16),  # Increase axis title size
        axis.text = element_text(size = 14),  # Increase axis text size
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # Use linewidth instead of size
        axis.line = element_line(linewidth = 0.5, color = "black"))  # Use linewidth instead of size


# Display the plot
print(p)


# Save the plot to a file
ggsave("paphia_PCA.png", plot = p, width = 10, height = 8, dpi = 150)
ggsave("paphia_PCA.svg", plot = p, width = 10, height = 8, dpi = 300, device = "svg")
ggsave("paphia_PCA_nolabs.pdf", plot = p, width = 10, height = 8, dpi = 300, device = "pdf")




