
library(tidyverse)
library(readxl)
library(growthcurver)
setwd("C:/Users/et491/OneDrive/Documents/MbyRes/tecan assays/Knockout panel growth curves")

df <- read_excel("tecan sunrise Robin 3.4.25.xlsx", sheet = "RAW")
dictionary <- read.csv("tecan sunrise dictionary.csv")

#changing seconds to hours
df$Time <- as.numeric(as.character(df$Time))
df <- df %>% mutate(Time = Time/3600)

#changing to long format
df_long <- pivot_longer(df, !Time, names_to = "Well", values_to = "Absorbance")

#join the dictionary to the dataframe
df_dictionary <- df_long %>% left_join(dictionary)

#mutating the wells to being factors and more readable for the computer (no I don't really understand how this code works)
df_dictionary$Well= factor(df_dictionary$Well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))
df <- df_dictionary %>% mutate(row=gsub('^([A-H])(.*)', '\\1', Well), column=factor(gsub('^([A-H])(.*)', '\\2', Well), levels=1:12))

#creating a blanks dataframe to normalise your data to the blank 
blks_df <- df %>% filter(Phage == "BLK") %>% group_by(Time) %>% summarise(mean_blk_od = mean(Absorbance))

#normalising your data to the blank 
df <- df %>% left_join(blks_df) %>% mutate(norm_od = Absorbance - mean_blk_od)

#filtering so it is now just hsdR and wildtype
df_hsdR_WT <- df %>% 
  filter(Phage %in% c("PAO1_1","PAO1_2","PAO1_3","PAO1_4","PAO1_5","PAO1_6", "PAO1 hsdR_1", "PAO1 hsdR_2", "PAO1 hsdR_3", "PAO1 hsdR_4", "PAO1 hsdR_5", "PAO1 hsdR_6"))

#simplify
df_simplified <- df_hsdR_WT %>% 
  select(Phage, norm_od, Time)

#converting back to a wide format for growthcurver with columns for each strain and for time
df_wide <- df_simplified %>%
  pivot_wider(names_from = Phage, values_from = norm_od)


#fitting logistic growth models to each repeat of each strain
PAO1_1_fit <- SummarizeGrowth(df_wide$Time, df_wide$PAO1_1)
PAO1_2_fit <- SummarizeGrowth(df_wide$Time, df_wide$PAO1_2)
PAO1_3_fit <- SummarizeGrowth(df_wide$Time, df_wide$PAO1_3)
PAO1_4_fit <- SummarizeGrowth(df_wide$Time, df_wide$PAO1_4)
PAO1_5_fit <- SummarizeGrowth(df_wide$Time, df_wide$PAO1_5)
PAO1_6_fit <- SummarizeGrowth(df_wide$Time, df_wide$PAO1_6)
PAO1_hsdR_1_fit <- SummarizeGrowth(df_wide$Time, df_wide$`PAO1 hsdR_1`)
PAO1_hsdR_2_fit <- SummarizeGrowth(df_wide$Time, df_wide$`PAO1 hsdR_2`)
PAO1_hsdR_3_fit <- SummarizeGrowth(df_wide$Time, df_wide$`PAO1 hsdR_3`)
PAO1_hsdR_4_fit <- SummarizeGrowth(df_wide$Time, df_wide$`PAO1 hsdR_4`)
PAO1_hsdR_5_fit <- SummarizeGrowth(df_wide$Time, df_wide$`PAO1 hsdR_5`)
PAO1_hsdR_6_fit <- SummarizeGrowth(df_wide$Time, df_wide$`PAO1 hsdR_6`)



# Create a dataframe to store the results
results_df <- data.frame(
  sample = character(12),  # Adjust based on the number of samples
  N0 = numeric(12),
  K = numeric(12),
  r = numeric(12),
  DT = numeric(12),
  stringsAsFactors = FALSE
)


# Fit data and store the results
fit_list <- list(
  PAO1_1_fit, PAO1_2_fit, PAO1_3_fit, PAO1_4_fit, PAO1_5_fit, PAO1_6_fit,
  PAO1_hsdR_1_fit, PAO1_hsdR_2_fit, PAO1_hsdR_3_fit, PAO1_hsdR_4_fit, PAO1_hsdR_5_fit, PAO1_hsdR_6_fit
)

# Sample names (you can adjust this based on your data)
sample_names <- c("PAO1_1", "PAO1_2", "PAO1_3", "PAO1_4", "PAO1_5", "PAO1_6",
                  "PAO1_hsdR_1", "PAO1_hsdR_2", "PAO1_hsdR_3", "PAO1_hsdR_4", "PAO1_hsdR_5", "PAO1_hsdR_6")

# Loop through the fit results and extract the required parameters
for (i in 1:length(fit_list)) {
  # Extract relevant parameters from each fit
  fit <- fit_list[[i]]
  
  N0 <- fit$vals$n0
  K <- fit$vals$k
  r <- fit$vals$r
  DT <- log(2) / fit$vals$r  # Doubling time, derived from growth rate (r)
  
  # Store the values in the dataframe
  results_df[i, ] <- c(sample_names[i], N0, K, r, DT)
}

# View the results
print(results_df)


strains_df <- results_df %>% 
  mutate(strain = ifelse(sample %in% c("PAO1_1", "PAO1_2", "PAO1_3", "PAO1_4", "PAO1_5", "PAO1_6"), "PAO1", "hsdR"))

strains_df$N0 <- as.numeric(strains_df$N0)
strains_df$K <- as.numeric(strains_df$K)
strains_df$r <- as.numeric(strains_df$r)
strains_df$DT <- as.numeric(strains_df$DT)

means_df <- strains_df %>% group_by(strain) %>% summarise(N0 = mean(N0), K = mean(K), r = mean(r), DT = mean(DT))



# Subset the data based on the strain group (PAO1 vs hsdR)
PAO1_r_values <- strains_df$r[strains_df$strain == "PAO1"]
hsdr_r_values <- strains_df$r[strains_df$strain == "hsdR"]

# Perform a t-test to compare the growth rates (r) between PAO1 and hsdR
r_t_test_result <- t.test(PAO1_r_values, hsdr_r_values)
print(r_t_test_result)


# Subset the data based on the strain group (PAO1 vs hsdR)
PAO1_K_values <- strains_df$K[strains_df$strain == "PAO1"]
hsdr_K_values <- strains_df$K[strains_df$strain == "hsdR"]

# Perform a t-test to compare the K between PAO1 and hsdR
K_t_test_result <- t.test(PAO1_K_values, hsdr_K_values)
print(K_t_test_result)


# Subset the data based on the strain group (PAO1 vs hsdR)
PAO1_N0_values <- strains_df$N0[strains_df$strain == "PAO1"]
hsdr_N0_values <- strains_df$N0[strains_df$strain == "hsdR"]

# Perform a t-test to compare the N0 between PAO1 and hsdR
N0_t_test_result <- t.test(PAO1_N0_values, hsdr_N0_values)
print(N0_t_test_result)


# Subset the data based on the strain group (PAO1 vs hsdR)
PAO1_DT_values <- strains_df$DT[strains_df$strain == "PAO1"]
hsdr_DT_values <- strains_df$DT[strains_df$strain == "hsdR"]

# Perform a t-test to compare the DT between PAO1 and hsdR
DT_t_test_result <- t.test(PAO1_DT_values, hsdr_DT_values)
print(DT_t_test_result)
