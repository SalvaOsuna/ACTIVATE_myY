#genetic gain####
 #load data
 #BLUPs from SpATS package
library(readr)
myY <- read.csv("data/lentil_blups_wide_spats.csv")
myY_sub <- read.csv("data/lentil_sub_blups_wide_spats.csv")
#plot genetic gain
#install.packages("vctrs")
#install.packages("gander")
#library(gander)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr) # Essential for adding P-values/R2 easily

# 1. SETUP
df <- myY

# 2. DATA WRANGLING
trait_of_interest <- "DS"

plot_data <- df |>
  # Select ID, metadata, and columns containing the specific trait
  select(Genotype, dev_year, cot_color, ends_with(trait_of_interest)) |>
  
  # Pivot to Long Format
  # names_pattern regex explained: "(.*)_(.*)" captures (Env)_(Trait)
  pivot_longer(
    cols = ends_with(trait_of_interest),
    names_to = c("Environment", "Trait"),
    names_pattern = "(.*)_(.*)", 
    values_to = "BLUP_Value"
  )

# Preview to ensure structure is correct
head(plot_data)

# 3. VISUALIZATION

ggplot(plot_data, aes(x = dev_year, y = BLUP_Value, color = cot_color)) +
  
  # A. The Scatter Points
  geom_point(alpha = 0.6, size = 2) + 
  
  # B. The Regression Line (Genetic Gain Trend)
  geom_smooth(method = "lm", se = F, alpha = 0.1) +
  
  # C. Facet by Environment (Site-Year)
  facet_wrap(~ Environment, scales = "free_y") +
  
  # D. Add Statistics (R2 and P-value)
  # label.y.npc = "top": Places stats at the top of the plot area
  # label.x.npc = "left": Aligns stats to the left
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4,
    show.legend = FALSE 
  ) +
  
  # E. Custom Colors 
  scale_color_manual(values = c("yellow" = "#FFD700", 
                                "red" = "#D55E00", 
                                "green" = "#009E73")) +
  
  # F. Formatting
  labs(
    title = paste("Genetic Gain Analysis:", trait_of_interest),
    subtitle = "Regression of Predicted means on Release Year by cotyledon color",
    x = "Year of Development",
    y = paste("Predicted mean (", trait_of_interest, ")", sep=""),
    color = "Cotyledon Color"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

ggsave(paste("Predicted mean","_", trait_of_interest,".png", sep=""), width = 10, height = 5, bg = "white")
