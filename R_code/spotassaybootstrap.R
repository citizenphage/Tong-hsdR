library(ggplot2)
library(tidyverse)


setwd("C:/Users/et491/OneDrive/Documents/MbyRes/Paper/PAO1 hsdR data")
df = read.csv("all PAO1 hsdR PFU data.csv")


df_wide <- df %>%
  pivot_wider(names_from = strain, values_from = PFU)

df_noinf <- df_wide %>%
  filter(PAO1 != 0) %>%
  filter(PAO1_hsdr !=0)

df_EoP <- df_noinf %>%
  mutate(EoP = PAO1_hsdr/PAO1)
  

mean(df_EoP$EoP)

bootstrap_mean_eop <- function(data, func, size=nrow(data)){
  bs_sample <- sample(data, size=size, replace=TRUE)
  do.call(func, list(x=bs_sample))
}

generate_bootstrap_replicates <- function(data, func, n=10000, size=nrow(data)){
  replicate(n, do.call(bootstrap_mean_eop, list(data=data, func=func, size=size)))
  
}

bootstrap_means <- tibble(generate_bootstrap_replicates(df_EoP$EoP, mean, size=84, n=10000), .name_repair = ~ c("mean_EoP"))



ggplot(bootstrap_means, aes(x=mean_EoP)) + 
  geom_histogram(color="#000000", fill="#0072B2", bins=30) +
  scale_x_continuous(expression("Mean EoP"), breaks = scales::pretty_breaks(n = 11)) +
  scale_y_continuous(expression("Count"), limits = c(0, 1250), expand = c(0,0), breaks = scales::pretty_breaks(n = 7)) +
  geom_vline(xintercept=1, size=1, linetype = "longdash") +
theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))
