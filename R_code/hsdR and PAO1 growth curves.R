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


#graph for hsdR and WT growth curves 16 hours
df_16 <- df_hours %>% 
  filter(time_hours <= 16)
df_hsdR_WT <- df_16 %>% 
  filter(Phage == c("PAO1", "PAO1 hsdR"))

custom_labels <- as_labeller(c(
  "PAO1" = "PAO1",
  "PAO1 hsdR" = "Delta * italic('hsdR')"), label_parsed)

hsdR_WT_facetted_growth_curves_16 <- ggplot(data = df_hsdR_WT, aes(x = time_hours, y = OD600)) +
  geom_line(size=1.5, colour = "#56B4E9") +
  geom_ribbon(aes(ymin= CI95_min, ymax = CI95_max), alpha = 0.2, show.legend = FALSE, fill = "#56B4E9", colour = "#56B4E9")+
  scale_y_continuous('Mean OD600 (n=6)') +
  scale_x_continuous('Time (hours)',  breaks = seq(0,16, by=4), expand = expansion(mult = c(0.02, 0.02))) +
  theme(plot.margin = margin(10, 30, 10, 10)) +
  facet_wrap(~Phage, ncol = 1, labeller = custom_labels) +
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


