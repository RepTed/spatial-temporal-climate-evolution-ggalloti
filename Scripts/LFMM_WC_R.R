#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("LEA")
library(LEA)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)
library(tidyr)

#to run the lfmm it requires the files in specific formats, so to prep them accordingly use the write.env and write.lfmm functions to re-write and re-read the files if needed.

#Generate a single PC for environmental variables
data <- read_table("env_WC_sst.env", col_names = FALSE)

#scale the data
data_scaled <- scale(data)

#run pca
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

#extrac PC1
first_pc <- pca_result$x[, 1]
write.csv(first_pc, "1PC_wc_subset.csv", row.names = FALSE)

#PC1 = 85% explained variance
#loadings:
#ann.mean.temp 0.44 positive strong
#max.warm.month 0.4 positive strong
#ann.precip -0.43 negative strong
#precip.season 0.35 positive strong
#ELEV -0.43 negative strong
#SRAD.MEAN -0.35 negative strong

#re-write .env file with 1st PC
write.env(first_pc, "PC1_WC_sst.env")

#1st PC represents gradient where low elevation and high mean temperature, along with sky radiant temp


#---------------optional to investigate PC------------------
# Calculate the proportion of variance and loading
explained_variance <- summary(pca_result)$importance[2, ]
loadings <- pca_result$rotation[, 1]
#----------------------------------------------------------

#run lfmm
#K-value determined prior, see "finding K" script
lfmm_project <- lfmm("geno_clim_sst.lfmm", "PC1_WC_sst.env", 
                     K = 3, 
                     repetitions = 5, 
                     project = "new")

#change d = for each environmental variable
# compute adjusted p-values
p1 = lfmm.pvalues(lfmm_project, K = K, d = 1)
pvalues = p1$pvalues

# GEA significance test
par(mfrow = c(2,1))
hist(pvalues, col = "lightblue")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)

# Convert p-values to -log10 scale for visualization
log_pvalues <- -log10(pvalues)

# Create a data frame for ggplot
results_df <- data.frame(
  SNP = seq_along(log_pvalues),
  log_pvalues = log_pvalues
)
snp_pos<-colnames(AllFreq)
snp_pos<-as.data.frame(snp_pos)
results_df <- results_df %>%
  separate(pos, into = c("gene", "pos"), sep = "_", convert = TRUE)
# Simulate significant SNPs (adjust the threshold as needed)
results_df$significant <- results_df$log_pvalues > -log10(0.05 / nrow(results_df))


# Create a Manhattan Plot with contrasting colors and significant SNP labels
manhattan_plot<-ggplot(results_df, aes(x = SNP, y = log_pvalues, color = gene)) +
  geom_point(aes(size = significant), alpha = 0.8) +
  theme_classic(base_size = 16) +
  labs(title = "Manhattan Plot of LFMM Results (PC1, 85%)",
       x = "SNP Index",
       y = "-log10(p-value)",
       color = "gene") +
  geom_hline(yintercept = -log10(0.05 / nrow(results_df)), col = "red", linetype = "dashed") +
  theme(legend.position = "right") +
  scale_x_discrete(labels = function(x) gsub("SNP", "", x)) +  
  geom_text_repel(data = subset(results_df, significant), max.overlaps = Inf,
                  aes(label = paste(gene, pos, sep = ", "),
                      color = gene),
                  size = 4, box.padding = 0.5)