#Code by Edward Gilbert, 2024

###--------------------------------------------Plotting-ADMIXTURE------------------------------------------------###

##ASSUMES YOU HAVE ALREADY GENERATED ADMIXTURE Q FILES
##For a tutorial to run ADMIXTURE look here: https://speciationgenomics.github.io/ADMIXTURE/ Copyright Mark Ravient & Joana Meier.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Directory containing the ADMIXTURE .Q files
admix_dir <- "Admixture/"

# Load the sample populations file
sample_pops_file <- file.path(admix_dir, "sample_pops.csv")
sample_pops_df <- read_csv(sample_pops_file)

# Ensure the Population column is a factor
sample_pops_df$Population <- factor(sample_pops_df$Population)

# Function to read .Q files
read_Q_file <- function(file_path) {
  read_delim(file_path, delim = " ", col_names = FALSE, trim_ws = TRUE)
}

# Find all .Q files in the directory
Q_files <- list.files(admix_dir, pattern = "spatial_control\\.[0-9]+\\.Q", full.names = TRUE)

# Sort files by K value
Q_files <- Q_files[order(as.numeric(gsub(".+\\.([0-9]+)\\.Q", "\\1", Q_files)))]

# Loop through each Q file and create a plot
for (Q_file in Q_files) {
  # Extract K value from the file name
  K <- as.numeric(gsub(".+\\.([0-9]+)\\.Q", "\\1", Q_file))
  
  # Read the Q file
  Q_df <- read_Q_file(Q_file)
  
  # Add sample and population information to the dataframe
  Q_df <- Q_df %>%
    mutate(Sample = sample_pops_df$Sample, Population = sample_pops_df$Population) %>%
    gather(key = "Cluster", value = "Ancestry_Proportion", -Sample, -Population) %>%
    mutate(Cluster = factor(Cluster, levels = paste0("X", 1:K)))
  
  # Create the plot
  p <- ggplot(Q_df, aes(x = Sample, y = Ancestry_Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    facet_grid(~ Population, scales = "free_x", space = "free_x") +
    labs(x = "Localities", y = "Ancestry Proportion", title = paste("ADMIXTURE Plot for K=", K)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top") +
    guides(fill = guide_legend(nrow = 1))
  
  # Print the plot to the RStudio plot viewer
  print(p)
  
  # Save the plot as a PNG file
  ggsave(filename = file.path(admix_dir, paste0("admixture_plot_K", K, ".png")), plot = p, width = 10, height = 6)
}


###--------------------------------------------DAPC----------------------------------------------------###

vcf <- read.vcfR("VCFs//control.vcf.gz")
metadata<-read.csv("sample_site_info.csv")

genlight<- vcfR2genlight(vcf)
num_clust <- find.clusters(genlight)
#PCs retained = 40
#number of clusters = 4


dapc <- dapc(genlight, num_clust$grp)
#PCs retained = 12
#too many PCs can result in overfitting, while too few can result in overdispersion
#discriminant functions retained = 3


scatter(dapc, posi.da="bottomleft")

dapc_data_df <-
  as_tibble(dapc$ind.coord, rownames = "individual") %>%
  mutate(population = metadata$Location)

dapc_data_df$env<-metadata$Env
dapc_data_df$pop<-metadata$Pop

dapc_plot <-
  ggplot(dapc_data_df, aes(
    x = LD1,
    y = LD2,
    fill = population
  )) +
  geom_point(shape = 21, size = 3) + 
  theme_bw(base_size = 16)

dapc_plot

#k=3
###-------------------------------------------------------snmf-------------------------------------------------------###

#NOTE: LEA requires the data in specific formats. The best way to do this is to save your genotype matrix using the function:
write.geno(control_geno, "control-loci.geno")

library(reshape)

genotype_data <- read.geno("control-loci.geno")
snmf <- snmf("control-loci.geno", K = 3:8, 
                repetitions = 5, entropy = TRUE, 
                project = "new", alpha = 100)


plot(snmf, cex = 1.2, col = "lightblue", pch = 19)
ce <-  cross.entropy(snmf, K = 3)
ce
best_run <- which.min(ce)
best_run
q_mat <- LEA::Q(snmf, K = 3, run = best_run) 
colnames(q_mat) <- paste0("P", 1:3)
head(q_mat)

q_df <- as.data.frame(q_mat)
q_df$Individual <- factor(1:nrow(q_df))
metadata <- data.frame(
  Individual = 1:236,  
  Locality = locality
)
q_df$Individual <- as.integer(as.character(q_df$Individual))
q_df <- merge(q_df, metadata, by = "Individual")
# Melt the dataframe for ggplot2 plotting
q_melt <- melt(q_df, id.vars = c("Individual", "Locality"))

# Order q_melt dataframe by Locality and then by Individual (if needed)
q_melt <- q_melt[order(q_melt$Locality, q_melt$Individual), ]

ggplot(q_melt, aes(x = Locality, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("skyblue", "orange", "lightgreen", "red", "darkblue", "pink", "brown", "purple"), name = "Ancestry") +
  theme_classic(base_size = 16) +
  labs(x = "Locality", y = "Ancestry Proportion", 
       title = "Admixture Plot by Locality") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot(q_melt, aes(x = Individual, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("skyblue", "orange", "lightgreen", "red", "darkblue", "pink", "brown", "purple"), name = "Ancestry") +
  theme_minimal() +
  labs(x = "Individual", y = "Ancestry Proportion", 
       title = "Admixture Plot with Individuals Ordered by Locality") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###------------------------------------screeplot----------------------------------------------###
library(vegan)

# Extract the genotype matrix
control_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Function to convert genotypes to numeric format
geno_to_numeric <- function(geno) {
  if (is.na(geno)) {
    return(9)  # NA values are converted to 9
  } else if (geno == "0/0") {
    return(0)
  } else if (geno == "0/1" || geno == "1/0") {
    return(1)
  } else if (geno == "1/1") {
    return(2)
  } else {
    return(9)  # Any other unexpected value is also treated as missing
  }
}

# Apply the conversion function to the genotype matrix
numeric_control_matrix <- apply(control_matrix, 2, function(x) sapply(x, geno_to_numeric))

# Write the numeric genotype matrix to a .geno file
write.table(numeric_control_matrix, "spatial_control.geno", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "")
write.table(numeric_control_matrix, "spatial_control.txt")
gen.pca <- rda(numeric_control_matrix, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")
#The broken stick stopping rule states that principal components should be retained as long as 
#observed eigenvalues are higher than corresponding random broken stick components.

#K = 5-6
