#Code by Edward Gilbert, 2024

#Generate Tajima's D across sequences in a sliding window approach

library(ape)
library(pegas)

# Define function to calculate Tajima's D in sliding windows
sliding_window_tajimas_d <- function(dna, window_size, step_size) {
  num_sites <- ncol(dna)
  window_starts <- seq(1, num_sites - window_size + 1, by = step_size)
  results <- data.frame(Window_Start = numeric(), Tajima_D = numeric())
  
  for (start in window_starts) {
    end <- start + window_size - 1
    segment <- dna[, start:end]
    taj_d <- tajima.test(segment)$D
    results <- rbind(results, data.frame(Window_Start = start, Tajima_D = taj_d))
  }
  
  return(results)
}

# Get list of FASTA files in the current working directory
fasta_files <- list.files(pattern = "*master.fasta")
getwd()
# Initialize an empty data frame to store all results
all_results <- data.frame(Gene = character(), Window_Start = numeric(), Tajima_D = numeric(), stringsAsFactors = FALSE)

# Loop through each FASTA file
for (file in fasta_files) {
  # Extract gene name from file name (assuming it's the part before '.fa_aligned.fasta')
  gene_name <- sub("\\.fa_aligned\\.fasta", "", file)
  
  # Read the DNA sequence from the FASTA file
  dna <- read.dna(file, format = "fasta")
  
  # Define window size and step size
  window_size <- 500
  step_size <- 50
  
  # Calculate Tajima's D
  tajimas_d_values <- sliding_window_tajimas_d(dna, window_size, step_size)
  
  # Add gene name to the results
  tajimas_d_values$Gene <- gene_name
  
  # Append the results to the main data frame
  all_results <- rbind(all_results, tajimas_d_values)
}

# Optionally save the results to a CSV file with write.csv

###----------------------------------plotting tajD-----------------------------------------------###

##--------------------heatmap------------------------##
tajcont<-read.csv("Tables//tajimas_d_results_all_500-50_control.csv")
taj<-read.csv("Tables//tajimas_d_results_all_500-50_no_control.csv")
taj<-rbind(taj, tajcont)
# Convert relevant columns to numeric if needed
taj$Position <- as.numeric(taj$Position)
taj$D <- as.numeric(taj$D)
theme_set(theme_minimal(base_size = 16))


# Average D across gene
# Calculate average D for each gene
average_D <- taj %>%
  group_by(Gene) %>%
  summarize(avg_D = mean(D, na.rm = TRUE))

# Reorder genes based on average D
taj$Gene <- factor(taj$Gene, levels = average_D$Gene[order(average_D$avg_D)])

# Create a heatmap with ordered genes
td_heat<-ggplot(taj, aes(x = Position, y = Gene, fill = D)) +
  geom_tile() +
  scale_fill_gradient2(low = "#377eb8", mid = "#999999", high = "#e41a1c", midpoint = 0) +
  labs(title = "Tajima's D across Gene Regions",
       x = "Genomic Position",
       y = "Gene") +  theme_bw(base_size = 16) +
  scale_x_continuous(breaks = seq(-2000, 1000, by = 100))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.2)
td_heat
##------------------------kerneldensity---------------------##

td_density<-ggplot(taj, aes(x = Position, y = D)) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 9) +
  scale_fill_distiller(palette = "Blues", direction = 1)+
  scale_x_continuous(breaks = seq(-2200, 500, by = 100), limits = c(-2200, 500),
                     expand = c(0, 0))+
  scale_y_continuous(limits = c(-3.5, 3.5)) +
  labs(title = "Climate genes",
       x = "Genomic Position",
       y = "Tajima's D") +
  theme_bw(base_size = 16) + 
  geom_point(alpha = 0.1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.2) +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.2)

td_density

###----------------------control genes KD----------------------###
td_density<-ggplot(taj, aes(x = Position, y = D)) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 6) +
  scale_fill_distiller(palette = "Blues", direction = 1)+
  scale_x_continuous(breaks = seq(-1500, 1500, by = 100), limits = c(-1500, 1500),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.5, 3)) +
  labs(title = "Control genes",
       x = "Genomic Position",
       y = "Tajima's D") +
  theme_bw(base_size = 16) + 
  geom_point(alpha = 0.1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.2) +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.2)
td_density
