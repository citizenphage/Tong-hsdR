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


ggplot(df_proportions, aes(x = source, y = proportion, fill = strain)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_y_continuous(limits = c(0, 12), expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 12, by = 3)) +
  scale_x_discrete(labels = c("Environmental", "Wastewater")) +
  scale_fill_manual(values = c("#009e73", "#0072b2"), 
                    labels = c("PAO1", expression(Delta * italic('hsdR')))) +
  labs(x = "Sample type", y = "Proportion of samples yielding phages (%)", fill = "Bacterial strain") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = 11.8, yend = 11.8)) +  # horizontal line
  geom_segment(aes(x = 0.8, xend = 0.8, y = 11.8, yend = 11.5)) +  # left vertical
  geom_segment(aes(x = 1.2, xend = 1.2, y = 11.8, yend = 11.5)) +  # right vertical
  annotate("text", x = 1, y = 12, label = "*", size = 10)


#Fisher's exact tests for sewage samples and environmental water samples separately

environmental_data <- read.csv("water sample isolation data.csv")

contingency_table_environmental <- as.matrix(environmental_data[, c("phage", "no.phage")])
fisher_result_environmental <- fisher.test(contingency_table_environmental)
print(fisher_result_environmental)


sewage_data <- read.csv("sewage sample isolation data.csv")

contingency_table_sewage <- as.matrix(sewage_data[, c("phage", "no.phage")])
fisher_result_sewage <- fisher.test(contingency_table_sewage)
print(fisher_result_sewage)
