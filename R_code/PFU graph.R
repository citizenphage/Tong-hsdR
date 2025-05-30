library(ggplot2)
library(tidyverse)
library(cowplot)


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


#Plot the values on a log-log plot

library(scales)

ggplot(df_EoP, aes(x = PAO1, y = PAO1_hsdr)) + 
  geom_point(size = 3) + 
  scale_x_continuous(trans = 'log10', 
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(trans = 'log10', 
                     labels = trans_format("log10", math_format(10^.x))) +
  geom_abline(intercept = 0, slope = 1, col = "red", lty = 1, size = 1.5) +
  geom_segment(aes(x = PAO1, xend = PAO1, y = PAO1_hsdr, yend = PAO1), 
               col = "blue", lty = 2, size = 1) +
  labs(y = expression("PFU/mL on " * Delta * italic('hsdR')), 
       x = "PFU/mL on wildtype PAO1") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))



#The distance from a point to the identity in the y-axis (the blue line) is simply y - x
#Because we're dealing with powers of 10, we can log the numbers.
distance2<-log10(df_avEoP$PAO1_hsdr) - log10(df_avEoP$PAO1)

#so, a value of -4 means that host_1 produced 4 orders of magnitude fewer plaques per mL than host_2

mean(distance2)
print(distance2)











