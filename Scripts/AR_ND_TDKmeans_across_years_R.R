#Code by Edward Gilbert, 2024

###------------------------------ALLELIC-RICHNESS-ACROSS-YEARS-FOR-EACH-ENVIRONMENT--------------------------------###

# Required libraries
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(hierfstat)
library(tibble)

envA.geno<-read.csv("AR plots//all_geno_envA.csv")
envB.geno<-read.csv("AR plots//all_geno_envB.csv")
envC.geno<-read.csv("AR plots//all_geno_envC.csv")
envD.geno<-read.csv("AR plots//all_geno_envD.csv")

# Filter out Year-Env combinations with fewer than 3 samples. Change Env type and re-run. 
filtered_dataC <- envC.geno %>%
  group_by(Year) %>%
  filter(n() >= 3)

filtered_dataC <- as.data.frame(filtered_dataC)

# Calculate allelic richness using hierfstat::allelic.richness
allelic_richnessA <- allelic.richness(filtered_dataA)
allelic_richnessB <- allelic.richness(filtered_dataB)
allelic_richnessC <- allelic.richness(filtered_dataC)
allelic_richnessD <- allelic.richness(filtered_dataD)

allelic_richnessA <- as.data.frame(allelic_richnessA)
allelic_richnessB <- as.data.frame(allelic_richnessB)
allelic_richnessC <- as.data.frame(allelic_richnessC)
allelic_richnessD <- as.data.frame(allelic_richnessD)

allelic_richnessA <- allelic_richnessA %>% mutate(Env = "A")
allelic_richnessB <- allelic_richnessB %>% mutate(Env = "B")
allelic_richnessC <- allelic_richnessC %>% mutate(Env = "C")
allelic_richnessD <- allelic_richnessD %>% mutate(Env = "D")


# Filter the data for each environment by Year, ensuring there are at least 3 samples per year
filtered_dataA <- envA.geno %>%
  group_by(Year) %>%
  filter(n() >= 3) %>%
  as.data.frame()

filtered_dataB <- envB.geno %>%
  group_by(Year) %>%
  filter(n() >= 3) %>%
  as.data.frame()

filtered_dataC <- envC.geno %>%
  group_by(Year) %>%
  filter(n() >= 3) %>%
  as.data.frame()

filtered_dataD <- envD.geno %>%
  group_by(Year) %>%
  filter(n() >= 3) %>%
  as.data.frame()

# Calculate allelic richness for each filtered dataset
ar_A <- allelic.richness(filtered_dataA)
ar_B <- allelic.richness(filtered_dataB)
ar_C <- allelic.richness(filtered_dataC)
ar_D <- allelic.richness(filtered_dataD)

# Convert each allelic richness result to a dataframe
ar_df_A <- as.data.frame(ar_A)
ar_df_B <- as.data.frame(ar_B)
ar_df_C <- as.data.frame(ar_C)
ar_df_D <- as.data.frame(ar_D)

# Reshape each dataframe from wide to long format and add Env column to indicate environment type
ar_long_A <- ar_df_A %>%
  rownames_to_column(var = "SNP") %>%
  pivot_longer(cols = starts_with("Ar"), names_to = "Year", values_to = "AllelicRichness") %>%
  mutate(Year = as.numeric(sub("Ar\\.", "", Year)), Env = "A")  # Extract the year and add Env type

ar_long_B <- ar_df_B %>%
  rownames_to_column(var = "SNP") %>%
  pivot_longer(cols = starts_with("Ar"), names_to = "Year", values_to = "AllelicRichness") %>%
  mutate(Year = as.numeric(sub("Ar\\.", "", Year)), Env = "B")

ar_long_C <- ar_df_C %>%
  rownames_to_column(var = "SNP") %>%
  pivot_longer(cols = starts_with("Ar"), names_to = "Year", values_to = "AllelicRichness") %>%
  mutate(Year = as.numeric(sub("Ar\\.", "", Year)), Env = "C")

ar_long_D <- ar_df_D %>%
  rownames_to_column(var = "SNP") %>%
  pivot_longer(cols = starts_with("Ar"), names_to = "Year", values_to = "AllelicRichness") %>%
  mutate(Year = as.numeric(sub("Ar\\.", "", Year)), Env = "D")


# Combine them into one dataframe for further analysis.
combined_allelic_richness <- bind_rows(ar_long_A, ar_long_B, ar_long_C, ar_long_D)


#optional for summary data--------------------------------
# Add 'Env' column to each dataset
filtered_dataA <- filtered_dataA %>% mutate(Env = "A")
filtered_dataB <- filtered_dataB %>% mutate(Env = "B")
filtered_dataC <- filtered_dataC %>% mutate(Env = "C")
filtered_dataD <- filtered_dataD %>% mutate(Env = "D")

# Combine all filtered datasets
combined_filtered_data <- bind_rows(
  filtered_dataA,
  filtered_dataB,
  filtered_dataC,
  filtered_dataD
)

# Summarize the combined data
filtered_data_summary <- combined_filtered_data %>%
  group_by(Year, Env) %>%
  summarize(SampleCount = n())
#---------------------------------------------------------

#Assuming model assumptions are met (not described here)
# Fit a GLM model
glm_model <- glm(AllelicRichness ~ Year * Env, data = combined_allelic_richness, family = Gamma(link = "log"))

# Summary of the model
summary(glm_model)
#Statistical test of the model
#If assumptions are met for glm, assumptions are met for this test
anova(glm_model, test = "Chisq")

# Create a new data frame for predictions
new_data <- expand.grid(Year = unique(combined_allelic_richness$Year), 
                        Env = unique(combined_allelic_richness$Env))

# Get predicted values
new_data$Predicted <- predict(glm_model, newdata = new_data, type = "response", se.fit = TRUE)$fit

# Calculate confidence intervals
new_data$CI <- predict(glm_model, newdata = new_data, type = "response", se.fit = TRUE)
new_data$CI_lower <- new_data$Predicted - 1.96 * new_data$CI$se.fit
new_data$CI_upper <- new_data$Predicted + 1.96 * new_data$CI$se.fit

#Plot
ggplot() +
  geom_line(data = new_data, aes(x = Year, y = Predicted, color = Env)) +
  geom_ribbon(data = new_data, aes(x = Year, ymin = CI_lower, ymax = CI_upper, fill = Env), alpha = 0.2) +
  labs(title = "Predicted Allelic Richness Over Time by Environment",
       x = "Year", y = "Predicted Allelic Richness") +
  theme_classic(base_size = 16) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("darkorange", "coral3", "chartreuse3", "antiquewhite4")) +
  theme(legend.position = "right")

###-------------------------------------NUCLEOTIDE-DIVERSITY-ACROSS-YEARS-FOR-EACH-ENVIRONMENT--------------------------------###

library(dplyr)
library(ape)
library(ggplot2)

# Set the directory containing the fasta files
fasta_dir <- "###"
metadata <- read.csv("###Sequence_info.csv")

# Initialize an empty list to store results
results_list <- list()

# List all fasta files in the directory
fasta_files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Function to extract data and calculate nucleotide diversity
process_fasta_file <- function(file_path) {
  # Extract file name and components
  file_name <- basename(file_path)
  parts <- strsplit(file_name, "_")[[1]]
  
  # Assuming the format is consistent
  gene_name <- parts[1]
  year <- parts[length(parts) - 1]
  env_type <- parts[length(parts)]
  
  # Read the fasta file
  fasta <- read.dna(file_path, format = "fasta")
  
  # Calculate nucleotide diversity
  diversity <- nuc.div(fasta)
  
  # Get the number of samples
  num_samples <- nrow(fasta)
  
  # Create a data frame with the results
  result <- data.frame(
    Gene = gene_name,
    Year = year,
    Env = env_type,
    NumSamples = num_samples,
    NucleotideDiversity = mean(diversity, na.rm = TRUE)
  )
  
  return(result)
}

# Loop through each fasta file and process it
for (fasta_file in fasta_files) {
  result <- process_fasta_file(fasta_file)
  results_list[[length(results_list) + 1]] <- result
}

# Combine all results into a single data frame
final_results <- bind_rows(results_list)

# Clean the Gene and Env columns and filter out low sample sizes
final_results_tidy <- final_results %>%
  mutate(
    Gene = gsub("\\.fa$", "", Gene),  # Remove ".fa" from Gene
    Env = gsub("\\.fasta$", "", Env)   # Remove ".fasta" from Env
  ) %>%
  filter(NumSamples >= 6)  # Keep only rows with NumSamples >= 6


# Calculate mean nucleotide diversity for each year and environment
mean_values <- final_results_tidy %>%
  group_by(Year, Env) %>%
  summarise(MeanNucleotideDiversity = mean(NucleotideDiversity, na.rm = TRUE), .groups = 'drop')

# Create the boxplot with mean trendline
ggplot(final_results_tidy, aes(x = factor(Year), y = NucleotideDiversity)) +
  geom_boxplot(aes(fill = Env), outlier.shape = NA) +  # Boxplot for Nucleotide Diversity
  geom_line(data = mean_values, aes(x = factor(Year), y = MeanNucleotideDiversity, group = Env, color = Env), size = 1, show.legend = TRUE) +  # Trendline for mean
  geom_point(data = mean_values, aes(x = factor(Year), y = MeanNucleotideDiversity, color = Env), size = 3) +  # Points for mean values
  labs(title = "Nucleotide Diversity Across Years by Gene",
       x = "Year",
       y = "Nucleotide Diversity") +
  scale_fill_brewer(palette = "Set1") +  # Color palette for boxplots
  scale_color_brewer(palette = "Set1") +  # Color palette for trendline
  theme_minimal(base_size = 14) +  # Minimal theme for clarity
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

#assuming model assumptions are met as before (not displayed in this script)
# Fit a GLM model for Nucleotide Diversity
glm_model_nd <- glm(NucleotideDiversity ~ Year * Env, data = final_results_tidy, family = Gamma(link = "log"))

# Check summary of models
summary(glm_model_nd)
anova(glm_model_nd, test = "Chisq")


# Create a new data frame for predictions
new_data_nd <- expand.grid(Year = unique(final_results_tidy$Year), 
                           Env = unique(final_results_tidy$Env))

# Ensure Year is treated as a factor in predictions
new_data_nd$Year <- factor(new_data_nd$Year, levels = unique(final_results_tidy$Year))

# Get predicted values and standard errors
predictions <- predict(glm_model_nd, newdata = new_data_nd, type = "response", se.fit = TRUE)
new_data_nd$Predicted <- predictions$fit
new_data_nd$CI <- predictions$se.fit

# Calculate confidence intervals
new_data_nd$CI_lower <- new_data_nd$Predicted - 1.96 * new_data_nd$CI
new_data_nd$CI_upper <- new_data_nd$Predicted + 1.96 * new_data_nd$CI

# Create plot
ggplot() +
  geom_line(data = new_data_nd, aes(x = Year, y = Predicted, color = Env, group = Env), size = 1) +
  geom_ribbon(data = new_data_nd, aes(x = Year, ymin = CI_lower, ymax = CI_upper, fill = Env, group = Env), alpha = 0.2) +
  labs(title = "Nucleotide Diversity Over Time by Environment",
       x = "Year", y = "Nucleotide Diversity") +
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("blue", "green", "orange", "red")) +
  theme(legend.position = "right")


###---------------------------------------TAJIMA'S-D-CLUSTERING-AND-PLOTS------------------------------------------###
filtered_results<-read.csv("metrics_by-year-env.csv")

filtered_results_no_2014 <- filtered_results %>%
  filter(Year != 2014)

#Group by Gene and Year, calculate mean and SD
cluster_data <- filtered_results_no_2014 %>%
  group_by(Gene, Year) %>%
  summarise(
    TajimasD.D_mean = mean(TajimasD.D, na.rm = TRUE),    
    TajimasD.D_sd = sd(TajimasD.D, na.rm = TRUE)      
  ) %>%
  ungroup()

#Reshape
cluster_data_wide <- cluster_data %>%
  pivot_wider(names_from = Year, values_from = c(TajimasD.D_mean, TajimasD.D_sd))

#Impute missing values
cluster_data_imputed <- cluster_data_wide %>%
  mutate(across(starts_with("TajimasD.D_mean_20"), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

#Perform k-means clustering
set.seed(123)
k <- 3 
kmeans_result <- kmeans(cluster_data_imputed %>% select(-Gene), centers = k)

#Cluster assignment
cluster_data_imputed$Cluster <- as.factor(kmeans_result$cluster)

#Reshape
cluster_data_long <- cluster_data_imputed %>%
  pivot_longer(cols = starts_with("TajimasD.D_mean_20"), 
               names_to = "Year", 
               names_prefix = "TajimasD.D_mean_", 
               values_to = "TajimasD.D_mean") %>%
  left_join(
    cluster_data_imputed %>%
      pivot_longer(cols = starts_with("TajimasD.D_sd_20"), 
                   names_to = "Year", 
                   names_prefix = "TajimasD.D_sd_", 
                   values_to = "TajimasD.D_sd"),
    by = c("Gene", "Year", "Cluster")
  )

#Remove control genes
filtered_cluster_data <- cluster_data_long %>%
  filter(!Gene %in% c('CR.fa', 'CytB.fa', 'RAG1.fa'))

#Plot function
plot_cluster <- function(cluster_num) {
  ggplot(filtered_cluster_data %>% filter(Cluster == cluster_num), 
         aes(x = Year, y = TajimasD.D_mean, group = Gene, color = Gene)) +
    geom_line(size = 1) +  # Line for each gene
    geom_point(size = 3) +  # Points showing mean Tajima's D for each gene
    geom_ribbon(aes(ymin = TajimasD.D_mean - TajimasD.D_sd, ymax = TajimasD.D_mean + TajimasD.D_sd, fill = Gene),
                alpha = 0.2, linetype = "blank") +  # Shaded area for standard deviation
    labs(title = paste("Cluster", cluster_num, "- Tajima's D Trends per Gene (with SD)"),
         x = "Year", y = "Mean Tajima's D ± SD") +
    theme_classic(base_size = 16) +
    theme(legend.position = "right")
}

#display plots
cluster_plots <- lapply(1:k, plot_cluster)

for (i in 1:k) {
  plot <- plot_cluster(i)
  print(plot)
  Sys.sleep(2)
}

#Separate data frame for control genes
control_genes_data <- cluster_data_long %>%
  filter(Gene %in% c('CR.fa', 'CytB.fa', 'RAG1.fa'))

#plot function
plot_control_genes <- function() {
  ggplot(control_genes_data, aes(x = Year, y = TajimasD.D_mean, group = Gene, color = Gene)) +
    geom_line(size = 1) +  # Line for each control gene
    geom_point(size = 3) +  # Points showing mean Tajima's D for each control gene
    geom_ribbon(aes(ymin = TajimasD.D_mean - TajimasD.D_sd, ymax = TajimasD.D_mean + TajimasD.D_sd, fill = Gene),
                alpha = 0.2, linetype = "blank") +  # Shaded area for standard deviation
    labs(title = "Control Genes (CR.fa, CytB.fa, RAG1.fa): Tajima's D Trends",
         x = "Year", y = "Mean Tajimas D ± SD") +
    theme_classic(base_size = 16) +
    theme(legend.position = "right")
}

#Display
control_gene_plot <- plot_control_genes()
print(control_gene_plot)

