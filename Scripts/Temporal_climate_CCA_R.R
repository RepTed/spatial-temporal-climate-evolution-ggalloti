library(dplyr)
library(tidyr)
library(readr)
library(scales)

climate<-read.csv("monthly_AEMET_Climate_data_2012-2023.csv")

# Filter df
climate_filtered <- climate %>%
  filter(Climate_Variable %in% c("TM_MAX", "TM_MES", "TM_MIN", "P_MES"))

# Calculate monthly thresholds
monthly_quartiles <- climate_filtered %>%
  group_by(Indicativo, Month, Climate_Variable) %>%
  summarize(
    heatwave_threshold = if_else(Climate_Variable %in% c("TM_MAX", "TM_MES", "TM_MIN"), 
                                 quantile(Value, 0.9, na.rm = TRUE), 
                                 NA_real_),
    drought_threshold = if_else(Climate_Variable == "P_MES", 
                                quantile(Value, 0.1, na.rm = TRUE), 
                                NA_real_),
    .groups = "drop"
  ) %>%
  distinct(Indicativo, Month, Climate_Variable, .keep_all = TRUE)

# Join
climate_highlight <- climate_filtered %>%
  left_join(monthly_quartiles, by = c("Indicativo", "Month", "Climate_Variable")) %>%
  mutate(
    Heatwave = Climate_Variable %in% c("TM_MAX", "TM_MES") & Value >= heatwave_threshold,
    Drought = Climate_Variable == "P_MES" & Value <= drought_threshold
  )

#date type
climate_highlight <- climate_highlight %>%
  mutate(Year_Month = as.Date(paste(Year, Month, "01", sep = "-")))

# Plot
temp_plot <- ggplot(climate_highlight %>% filter(Climate_Variable %in% c("TM_MAX", "TM_MES", "TM_MIN")),
                    aes(x = Year_Month, y = Value, color = Climate_Variable)) +
  geom_line() +
  geom_point() +
  geom_rect(data = subset(climate_highlight, Heatwave),
            aes(xmin = Year_Month - 15, xmax = Year_Month + 15, ymin = -Inf, ymax = Inf), 
            fill = "red", alpha = 0.2, inherit.aes = FALSE) +
  facet_wrap(~Indicativo) +
  labs(title = "Monthly Temperature Trends with Heatwave Highlights by Station",
       x = "Year",
       y = "Temperature (°C)",
       color = "Temperature Type") +
  scale_x_date(breaks = date_breaks("6 months"), labels = date_format("%Y-%m")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot
precip_plot <- ggplot(climate_highlight %>% filter(Climate_Variable == "P_MES"),
                      aes(x = Year_Month, y = Value, color = Indicativo)) +
  geom_line() +
  geom_point() +
  geom_rect(data = subset(climate_highlight, Drought),
            aes(xmin = Year_Month - 15, xmax = Year_Month + 15, ymin = -Inf, ymax = Inf), 
            fill = "blue", alpha = 0.2, inherit.aes = FALSE) +  # Adjusted transparency
  labs(title = "Monthly Precipitation Trends with Drought Highlights",
       x = "Year",
       y = "Total Precipitation (mm)",
       color = "Weather Station") +
  scale_x_date(breaks = date_breaks("6 months"), labels = date_format("%Y-%m")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plots
temp_plot
precip_plot


###--------------CLIMATE AND ADMIXTURE PROPORTIONS-------------------------###
library(dplyr)
library(purrr)
library(tidyr)
library(zoo)
library(ggplot2)

lag<-read.csv("admixture-climate.csv")

#Interpolation
interpolate_all_proportions <- function(df) {
  df$Average_Proportion <- na.approx(df$Average_Proportion, na.rm = FALSE)
  total <- sum(df$Average_Proportion, na.rm = TRUE)
  
  if (total > 0) {
    df <- df %>%
      mutate(Average_Proportion = Average_Proportion / total)
  }
  
  return(df)
}

#Prep
lag <- lag %>%
  filter(!is.na(pop) & pop != "") %>% 
  complete(year, pop) %>% 
  group_by(pop) %>% 
  do(interpolate_all_proportions(.)) %>% 
  ungroup()

lag <- lag %>%
  mutate(year = as.integer(year))

#Plot
ggplot(lag, aes(x = year)) +
  geom_line(aes(y = Average_Proportion, color = pop), size = 1.5) +
  geom_point(aes(y = Average_Proportion, color = pop), size = 2) +
  geom_line(aes(y = Total_Heatwaves / 35), color = "red", linetype = "dashed") +
  geom_line(aes(y = Total_Droughts / 35), color = "blue", linetype = "dashed") +
  scale_y_continuous(
    name = "Average Admixture Proportion",
    sec.axis = sec_axis(~ . * 10, name = "Count of Heatwaves / Droughts")
  ) +
  labs(
    title = "Admixture Proportion vs. Climate Events Over Time",
    x = "Year",
    color = "Population"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title.y.right = element_text(color = "red"),
    axis.title.y = element_text(color = "black"),
    legend.position = "top"
  ) +
  scale_x_continuous(breaks = seq(min(lag$year), max(lag$year), by = 1))

# Min-max normalization across all years and populations for Average_Proportion
min_proportion <- min(lag$Average_Proportion, na.rm = TRUE)
max_proportion <- max(lag$Average_Proportion, na.rm = TRUE)
lag <- lag %>%
  mutate(
    Normalized_Proportion = (Average_Proportion - min_proportion) / (max_proportion - min_proportion)
  )


###---------------------------------CROSS-CORRELATION-ANALYSIS--------------------------------###
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

#Prep
lag <- lag %>%
  mutate(
    Total_Heatwaves = replace_na(Total_Heatwaves, 0),
    Total_Droughts = replace_na(Total_Droughts, 0),
    Average_Proportion = replace_na(Average_Proportion, 0)
  )

#CCA by population function
cross_correlation_analysis <- function(data, pop) {
  pop_data <- data %>% filter(pop == !!pop)
  
  max_lag <- 6
  ccf_heatwave <- ccf(pop_data$Total_Heatwaves, pop_data$Average_Proportion, lag.max = max_lag, plot = FALSE)
  ccf_drought <- ccf(pop_data$Total_Droughts, pop_data$Average_Proportion, lag.max = max_lag, plot = FALSE)
  
  heatwave_corr <- ccf_heatwave$acf
  drought_corr <- ccf_drought$acf
  lags <- seq(-max_lag, max_lag)
  
  if (length(heatwave_corr) < length(lags)) {
    heatwave_corr <- c(heatwave_corr, rep(NA, length(lags) - length(heatwave_corr)))
  }
  if (length(drought_corr) < length(lags)) {
    drought_corr <- c(drought_corr, rep(NA, length(lags) - length(drought_corr)))
  }
  data.frame(
    Lag = lags,
    Heatwave_Correlation = heatwave_corr,
    Drought_Correlation = drought_corr,
    Population = pop
  )
}

#Apply function and combine
ccf_results <- do.call(rbind, lapply(unique(lag$pop), function(p) cross_correlation_analysis(lag, p)))

#Remove negative lag
ccf_results_negative <- ccf_results %>%
  filter(Lag <= 0)

#Add age structure
age<-read.csv("predicted_age.csv")
median_age <- -median(age$predicted.age, na.rm = TRUE)
sd_age <- -sd(age$predicted.age, na.rm = TRUE)

#Plot
ggplot(ccf_results_negative, aes(x = Lag)) +
  geom_line(aes(y = Heatwave_Correlation, color = "Heatwave Correlation"), size = 1) +
  geom_line(aes(y = Drought_Correlation, color = "Drought Correlation"), size = 1) +
  geom_rect(aes(xmin = median_age - sd_age, xmax = median_age + sd_age, 
                ymin = -Inf, ymax = Inf), fill = "gray80", alpha = 0.4) +
  geom_vline(xintercept = median_age) +
  facet_wrap(~Population) +
  labs(title = "Cross-Correlation of Admixture Proportion with Climate Events",
       x = "Lag (years)",
       y = "Correlation Coefficient",
       color = "Climate Event") +
  theme_minimal(base_size = 14)




###-------------------------------------PERMUTATION-TESTING---------------------------------------###

permutation_test_ccf <- function(data, pop, num_permutations = 1000, max_lag = 6) {
  pop_data <- data %>% filter(pop == !!pop)
  
  #observed max correlation
  observed_ccf_heatwave <- ccf(pop_data$Average_Proportion, pop_data$Total_Heatwaves, lag.max = max_lag, plot = FALSE)$acf
  observed_ccf_drought <- ccf(pop_data$Average_Proportion, pop_data$Total_Droughts, lag.max = max_lag, plot = FALSE)$acf
  observed_stat_heatwave <- max(abs(observed_ccf_heatwave), na.rm = TRUE)
  observed_stat_drought <- max(abs(observed_ccf_drought), na.rm = TRUE)
  
  #prep vectors
  permuted_stats_heatwave <- numeric(num_permutations)
  permuted_stats_drought <- numeric(num_permutations)
  
  for (i in 1:num_permutations) {
    permuted_heatwaves <- sample(pop_data$Total_Heatwaves, replace = TRUE)
    permuted_droughts <- sample(pop_data$Total_Droughts, replace = TRUE)
    
    #Calculate
    permuted_ccf_heatwave <- ccf(pop_data$Average_Proportion, permuted_heatwaves, lag.max = max_lag, plot = FALSE)$acf
    permuted_ccf_drought <- ccf(pop_data$Average_Proportion, permuted_droughts, lag.max = max_lag, plot = FALSE)$acf
    
    #Store
    permuted_stats_heatwave[i] <- max(abs(permuted_ccf_heatwave), na.rm = TRUE)
    permuted_stats_drought[i] <- max(abs(permuted_ccf_drought), na.rm = TRUE)
  }
  
  #p-values
  p_value_heatwave <- mean(permuted_stats_heatwave >= observed_stat_heatwave)
  p_value_drought <- mean(permuted_stats_drought >= observed_stat_drought)
  
  #stats
  mean_heatwave <- mean(permuted_stats_heatwave)
  sd_heatwave <- sd(permuted_stats_heatwave)
  ci_heatwave <- quantile(permuted_stats_heatwave, c(0.025, 0.975))
  
  mean_drought <- mean(permuted_stats_drought)
  sd_drought <- sd(permuted_stats_drought)
  ci_drought <- quantile(permuted_stats_drought, c(0.025, 0.975))
  
  #Z-scores
  z_score_heatwave <- (observed_stat_heatwave - mean_heatwave) / sd_heatwave
  z_score_drought <- (observed_stat_drought - mean_drought) / sd_drought
  
  return(list(
    observed_stat_heatwave = observed_stat_heatwave,
    observed_stat_drought = observed_stat_drought,
    p_value_heatwave = p_value_heatwave,
    p_value_drought = p_value_drought,
    mean_heatwave = mean_heatwave,
    sd_heatwave = sd_heatwave,
    ci_heatwave = ci_heatwave,
    mean_drought = mean_drought,
    sd_drought = sd_drought,
    ci_drought = ci_drought,
    z_score_heatwave = z_score_heatwave,
    z_score_drought = z_score_drought
  ))
}

#Apply test
permutation_results <- do.call(rbind, lapply(unique(lag$pop), function(p) {
  result <- permutation_test_ccf(lag, p)
  data.frame(
    Population = p,
    Observed_Heatwave = result$observed_stat_heatwave,
    Observed_Drought = result$observed_stat_drought,
    P_Value_Heatwave = result$p_value_heatwave,
    P_Value_Drought = result$p_value_drought,
    Mean_Heatwave = result$mean_heatwave,
    SD_Heatwave = result$sd_heatwave,
    CI_Lower_Heatwave = result$ci_heatwave[1],
    CI_Upper_Heatwave = result$ci_heatwave[2],
    Z_Score_Heatwave = result$z_score_heatwave,
    Mean_Drought = result$mean_drought,
    SD_Drought = result$sd_drought,
    CI_Lower_Drought = result$ci_drought[1],
    CI_Upper_Drought = result$ci_drought[2],
    Z_Score_Drought = result$z_score_drought
  )
}))
