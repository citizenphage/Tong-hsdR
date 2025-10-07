###A - PFU scatter plot

library(ggplot2)
library(tidyverse)
library(cowplot)
library(scales)

setwd("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/PAO1 hsdR data")
df = read.csv("all PAO1 hsdR PFU data.csv")

#calculating the EoP
df_wide <- df %>%
  pivot_wider(names_from = strain, values_from = PFU)

df_EoP <- df_wide %>%
  mutate(EoP = PAO1_hsdr/PAO1)

#calculating the mean EoP, but removing P163 and P284 because they had a PFU of 0 on one strain giving infinite EoP
df_avEoP <- df_EoP %>%
  filter(phage!="P163") %>%
  filter(phage!="P284")

mean(df_avEoP$EoP)


#defining the distance between each point and the line y=x
df_distance <- df_avEoP %>% mutate(distance2 = log10(df_avEoP$PAO1_hsdr) - log10(df_avEoP$PAO1))

#Plot the values on a log-log plot
#converting 0 to 1 such that the phages which could not infect WT PAO1 can be plotted on the log plot
df_EoP <- df_EoP %>%
  mutate(PAO1 = ifelse(PAO1 == 0, 1, PAO1),
         zero_flag = ifelse(PAO1 == 0, TRUE, FALSE))  # flag true zeros

df_EoP <- df_EoP %>%
  mutate(PAO1_plot = ifelse(PAO1 == 0, 1, PAO1),
         PAO1_hsdr_plot = ifelse(PAO1_hsdr == 0, 1, PAO1_hsdr))

P1 <- ggplot(df_EoP, aes(x = PAO1_plot, y = PAO1_hsdr_plot)) +
  geom_abline(intercept = 0, slope = 1, col = "#666666", lty = 1, size = 1.5) +
  geom_segment(aes(x = PAO1_plot, xend = PAO1_plot, y = PAO1_hsdr_plot, yend = PAO1_plot),
               col = "#666666", lty = 2, size = 1) +
  geom_point(aes(fill = colours), shape = 21, color = "black", size = 3, stroke = 0.8) +
  scale_fill_manual(name=NULL, values = c("bulk_up" = "#3B9AB2", "hsdR_isolation" = "#F21A00", "none" = "black"), 
                    labels = c("bulk_up"="Phage used in bulk-up \nexperiment", "hsdR_isolation"=expression("Phage isolated on " * Delta * italic('hsdR'))), 
                    breaks = c("bulk_up", "hsdR_isolation")) +
  geom_text_repel(data = df_EoP %>% filter(PAO1 == 1) %>%mutate(phage = case_when(phage == "P163" ~ "CPL00163",phage == "P284" ~ "CPL00284",TRUE ~ phage)),
    aes(label = phage),size = 4,box.padding = 1,point.padding = 0.5,force = 2,force_pull = 0.5,nudge_x = 0.1,max.overlaps = Inf,segment.colour = NA) +
  scale_x_continuous(trans = log_trans(base = 10),
                     breaks = 10^seq(0, 10, by = 2),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log_trans(base = 10),
                     breaks = 10^seq(0, 10, by = 2),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(y = expression("PFU/mL on " * Delta * italic('hsdR')),
       x = "PFU/mL on wildtype PAO1") +
  theme_bw(18) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5))




###B - Bulk up scatter plot and supplementary data. 

library(ggplot2)
library(tidyverse)
library(cowplot)
library(scales)
library(ggrepel)
library(broom)

setwd("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/PAO1 hsdR data/bulk up experiment")
df = read.csv("PFUs after bulk up all data.csv")

#filtering the spotting out host
df_WT<- df %>% filter(Spotting_out_host == "PAO1_WT")
df_hsdR<- df %>% filter(Spotting_out_host == "PAO1_hsdR")

#averaging the PFU
df_average_WT <- df_WT %>% group_by(phage, bulk_up_host) %>% summarise(mean_PFU = mean(PFU))
df_average_hsdR <- df_hsdR %>% group_by(phage, bulk_up_host) %>% summarise(mean_PFU = mean(PFU))

#determining phages with significant difference between hsdR and WT bulk ups with Welch two sample T-tests (two sided)
stats_WT <- df_WT %>%
  group_by(phage) %>%
  do(tidy(t.test(log10(PFU) ~ bulk_up_host, data = .))) %>%
  ungroup()

stats_hsdR <- df_hsdR %>%
  group_by(phage) %>%
  do(tidy(t.test(log10(PFU) ~ bulk_up_host, data = .))) %>%
  ungroup()

#joining the calculated means and the stats output
df_WT_calculated <- df_average_WT %>%
  left_join(stats_WT, by = "phage")

df_hsdR_calculated <- df_average_hsdR %>%
  left_join(stats_hsdR, by = "phage")

df_wide_WT <- df_WT_calculated %>%
  pivot_wider(names_from = bulk_up_host, values_from = mean_PFU)
df_wide_hsdR <- df_hsdR_calculated %>%
  pivot_wider(names_from = bulk_up_host, values_from = mean_PFU)


#defining which phages experience a significant difference in PFU when bulked up on hsdR compared to WT
df_wide_WT <- df_wide_WT %>%
  mutate(significant = ifelse(p.value < 0.05,"Significant (p < 0.05)", "Not Significant"))

df_wide_hsdR <- df_wide_hsdR %>%
  mutate(significant = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not Significant"))

#graph from spot assay on WT
S1 <- ggplot(df_wide_WT, aes(x = PAO1_WT, y = PAO1_hsdR)) + 
  geom_abline(intercept = 0, slope = 1, col = "#666666", lty = 1, size = 1.5) +
  geom_segment(aes(x = PAO1_WT, xend = PAO1_WT, y = PAO1_hsdR, yend = PAO1_WT), 
               col = "#666666", lty = 2, size = 1) +
  geom_point(aes(fill = significant), shape = 21, color = "black", size = 3, stroke = 0.8) +
  scale_fill_manual(name = NULL, values = c("Significant (p < 0.05)" = "#F21A00", "Not Significant" = "black")) +
  geom_text_repel(data = df_wide_WT %>% filter(significant == "Significant (p < 0.05)"), aes(label = phage), size = 4, box.padding = 1, point.padding = 0.5, force = 2,
                  force_pull = 0.5, nudge_y = 0.3, max.overlaps = Inf, segment.colour = NA ) +
  scale_x_continuous(trans = 'log10', labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(trans = 'log10', labels = trans_format("log10", math_format(10^.x))) +
  labs(y = expression("Mean PFU/mL from " * Delta * italic("hsdR") * " bulk-up (n=3)"), 
       x = "Mean PFU/mL from wildtype PAO1 bulk-up (n=3)") +
  theme_bw(18) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5))


#Graph from spot assay on hsdR
P2 <- ggplot(df_wide_hsdR, aes(x = PAO1_WT, y = PAO1_hsdR,)) + 
  geom_abline(intercept = 0, slope = 1, col = "#666666", lty = 1, size = 1.5) +
  geom_segment(aes(x = PAO1_WT, xend = PAO1_WT, y = PAO1_hsdR, yend = PAO1_WT), 
               col = "#666666", lty = 2, size = 1) +
  geom_point(aes(fill = significant), shape = 21, color = "black", size = 3, stroke = 0.8) +
  scale_fill_manual(name = NULL, values = c("Significant (p < 0.05)" = "#F21A00", "Not Significant" = "black")) +
  geom_text_repel(data = df_wide_hsdR %>% filter(significant == "Significant (p < 0.05)"), aes(label = phage), size = 4, box.padding = 0.5, point.padding = 0.3, force = 2,
                  force_pull = 0.5, nudge_y = 0.3, max.overlaps = Inf, segment.colour = NA) +
  scale_x_continuous(trans = 'log10', labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(trans = 'log10', labels = trans_format("log10", math_format(10^.x))) +
  labs(y = expression("Mean PFU/mL from " * Delta * italic("hsdR") * " bulk-up (n=3)"), 
       x = "Mean PFU/mL from wildtype PAO1 bulk-up (n=3)") +
  theme_bw(18) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5))


#bargraphs showing individual repeats for each phage on each host
df_signif_hsdR <- df_average_hsdR %>%left_join(df_hsdR_calculated %>% select(phage, bulk_up_host, p.value), by = c("phage","bulk_up_host"))
df_signif_hsdR  <- df_signif_hsdR %>% mutate(bulk_up_host=factor(bulk_up_host, levels = c("PAO1_WT", "PAO1_hsdR")))
df_hsdR  <- df_hsdR %>% mutate(bulk_up_host=factor(bulk_up_host, levels = c("PAO1_WT", "PAO1_hsdR")))

df_stars_hsdR <- df_signif_hsdR %>% 
  filter(p.value < 0.05) %>% 
  group_by(phage) %>% 
  summarise(y_position = max(log(mean_PFU)) + 2, .groups = "drop")

S2 <- ggplot(df_signif_hsdR, aes(x = factor(phage), y = log(mean_PFU), fill = bulk_up_host)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", alpha = 0.6) +
  geom_point(data= df_hsdR, aes(x = factor(phage), y = log(PFU), color = bulk_up_host),shape = 21,
             color = "black",
             stroke = 0.8,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 3) +
  scale_fill_manual(values = c("PAO1_WT"="#00A08A", "PAO1_hsdR"="#F98400"), labels = c("PAO1_WT"="PAO1", "PAO1_hsdR"=expression(Delta * italic("hsdR")))) +
  scale_color_manual(values = c("PAO1_WT"="#00A08A", "PAO1_hsdR"="#F98400")) +
  labs(x = "Phage", y = "log10(PFU)", fill = "Bulk-up host") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 18) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(data = df_stars_hsdR, aes(x = factor(phage), y = y_position, label = "*"),
            inherit.aes = FALSE, size = 12)


df_signif_WT <- df_average_WT %>%left_join(df_WT_calculated %>% select(phage, bulk_up_host, p.value), by = c("phage","bulk_up_host"))
df_signif_WT  <- df_signif_WT %>% mutate(bulk_up_host=factor(bulk_up_host, levels = c("PAO1_WT", "PAO1_hsdR")))
df_WT  <- df_WT %>% mutate(bulk_up_host=factor(bulk_up_host, levels = c("PAO1_WT", "PAO1_hsdR")))

df_stars_WT <- df_signif_WT %>% 
  filter(p.value <= 0.05) %>% 
  group_by(phage) %>% 
  summarise(y_position = max(log(mean_PFU)) + 2, .groups = "drop")

S3 <- ggplot(df_signif_WT, aes(x = factor(phage), y = log(mean_PFU), fill = bulk_up_host)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", alpha = 0.6) +
  geom_point(data= df_WT, aes(x = factor(phage), y = log(PFU), color = bulk_up_host),shape = 21,
             color = "black",
             stroke = 0.8,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 3) +
  scale_fill_manual(values = c("PAO1_WT"="#00A08A", "PAO1_hsdR"="#F98400"), labels = c("PAO1_WT"="PAO1", "PAO1_hsdR"=expression(Delta * italic("hsdR")))) +
  scale_color_manual(values = c("PAO1_WT"="#00A08A", "PAO1_hsdR"="#F98400")) +
  labs(x = "Phage", y = "log10(PFU)", fill = "Bulk-up host") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 18) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(data = df_stars_WT, aes(x = factor(phage), y = y_position, label = "*"),
            inherit.aes = FALSE, size = 12)



###C - isolation data bar graph and statistical test

library(tidyverse)
library(DescTools)
library(cowplot)
library(dplyr)

setwd("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/PAO1 hsdR data")
df = read.csv("phage isolation hsdR data.csv")


df_simple <- df %>%
  select(-by) %>%
  group_by(strain, source) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) 

df_proportions <- df_simple %>% mutate(proportion = number/total_samples*100)


P3 <- ggplot(df_proportions, aes(x = source, y = proportion, fill = strain)) + 
  geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.9)) +
  scale_y_continuous(limits = c(0, 12), expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 12, by = 3)) +
  scale_x_discrete(labels = c("Environmental", "Wastewater")) +
  scale_fill_manual(name = "Enrichment host", values = c("#00A08A", "#DD8D29"), 
                    labels = c("PAO1", expression(Delta * italic('hsdR')))) +
  labs(x = "Sample type", y = "Proportion of samples yielding phages (%)", fill = "Bacterial strain") +
  theme_bw(18) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = 11.8, yend = 11.8)) +  # horizontal line
  geom_segment(aes(x = 0.8, xend = 0.8, y = 11.8, yend = 11.5)) +  # left vertical
  geom_segment(aes(x = 1.2, xend = 1.2, y = 11.8, yend = 11.5)) +  # right vertical
  annotate("text", x = 1, y = 12, label = "*", size = 12)


#Fisher's exact tests for sewage samples and environmental water samples separately

environmental_data <- read.csv("water sample isolation data.csv")

contingency_table_environmental <- as.matrix(environmental_data[, c("phage", "no.phage")])
fisher_result_environmental <- fisher.test(contingency_table_environmental)
print(fisher_result_environmental)


sewage_data <- read.csv("sewage sample isolation data.csv")

contingency_table_sewage <- as.matrix(sewage_data[, c("phage", "no.phage")])
fisher_result_sewage <- fisher.test(contingency_table_sewage)
print(fisher_result_sewage)



###Putting together Figure 1
library(patchwork)
patchwork_1 <- P1 / P2
patchwork_1

patchwork_2 <- patchwork_1 | P3
patchwork_2 + plot_annotation(tag_levels = list(c("A", "B", "C")))


#Patching together suplementary bar graphs
supplementary_bargraphs <- S2/S3
supplementary_bargraphs + plot_annotation(tag_levels = 'A')




###supplementary figure 1 - growth curves

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


library(tidyverse)
library(readxl)
setwd("C:/Users/et491/OneDrive/Documents/MbyRes/tecan assays/Knockout panel growth curves")

df <- read_excel("tecan sunrise Robin 3.4.25.xlsx", sheet = "RAW")
dictionary <- read.csv("tecan sunrise dictionary 2.csv")

#changing to long format
df_long <- pivot_longer(df, !Time, names_to = "Well", values_to = "Absorbance")

#join the dataframes
df_dictionary <- df_long %>% left_join(dictionary)


#mutating the wells to being factors and more readable for the computer (no I don't really understand how this code works)
df_dictionary$Well= factor(df_dictionary$Well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))
df <- df_dictionary %>% mutate(row=gsub('^([A-H])(.*)', '\\1', Well), column=factor(gsub('^([A-H])(.*)', '\\2', Well), levels=1:12))


#creating a blanks dataframe to normalise your data to the blank 
blks_df <- df %>% filter(Phage == "BLK") %>% group_by(Time) %>% summarise(mean_blk_od = mean(Absorbance))

#normalising your data to the blank 
df <- df %>% left_join(blks_df) %>% mutate(norm_od = Absorbance - mean_blk_od)

df$Time <- as.numeric(as.character(df$Time))

#Creating a growth curve per well
overall_growth_curves <- ggplot(df, aes(x = Time, y = norm_od)) +
  geom_line(color='#0072B2') +
  scale_x_continuous('Time (Hours)') +
  scale_y_continuous('OD600') +
  ylab("Mean OD600") +
  facet_wrap(~Well, ncol = 12)
overall_growth_curves

#Mean Growth Curve with StDev
#You need to create a new dataframe here to group things and mutate it to be the mean 
df_mean <- df %>% group_by(Phage, Time) %>% summarise(Time = Time, Phage = Phage, OD600 = mean(norm_od), STDEV = sd(norm_od))
#removing the blk wells
df_mean <- df_mean %>% filter(!Phage %in% c("BLK", "Protector"))
#Graph time!
mean_growth_curve <- ggplot(data = df_mean, aes(x = Time, y = OD600, colour = Phage)) +
  geom_line() +
  geom_errorbar(aes(ymin = OD600 - STDEV, ymax = OD600 + STDEV), width = 0.2) +
  scale_y_continuous('Mean OD600 (n=10)') +
  scale_x_continuous('Time (Hours)') +
  ggtitle("XX") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
mean_growth_curve

## calculating the standard error and CI95
df_SE <- df_mean %>% mutate(SE = STDEV/sqrt(6))
df_CI <- df_SE %>% mutate(CI95 = 1.96*SE)
df_CI <- df_CI %>% mutate(CI95_min = OD600 - CI95) %>% mutate(CI95_max = OD600 + CI95)

#changing seconds to hours
df_hours <- df_CI %>% mutate(time_hours = Time/3600)

#custom labels for mutants

custom_labels <- as_labeller(c(
  "PAO1" = "PAO1",
  "PAO1 hsdR" = "Delta * italic('hsdR')"))

#just hsdR and WT growth curves 16 hours
df_16 <- df_hours %>% 
  filter(time_hours <= 16)

df_hsdR_WT <- df_16 %>% 
  filter(Phage == c("PAO1", "PAO1 hsdR"))

hsdR_WT_facetted_growth_curves_16 <- ggplot(data = df_hsdR_WT, aes(x = time_hours, y = OD600)) +
  geom_line(size=1.5, colour = "#56B4E9") +
  geom_ribbon(aes(ymin= CI95_min, ymax = CI95_max), alpha = 0.2, show.legend = FALSE, fill = "#56B4E9", colour = "#56B4E9")+
  scale_y_continuous('Mean OD600 (n=6)') +
  scale_x_continuous('Time (hours)',  breaks = seq(0,16, by=4), expand = expansion(mult = c(0.02, 0.02))) +
  theme(plot.margin = margin(10, 30, 10, 10)) +
  facet_wrap(~Phage, ncol = 1, labeller=custom_labels) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#aca9a9", size=0.4),
        panel.grid.minor = element_line(colour="#aca9a9", size=0.05),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        panel.spacing.x = unit(1.2, "lines"))
hsdR_WT_facetted_growth_curves_16

