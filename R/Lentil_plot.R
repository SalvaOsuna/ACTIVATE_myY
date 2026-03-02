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
df <- myY_sub

# 1. Define the trait you want to plot
target_trait <- "C.N.ratio.stover"

# 2. Prepare the data (Wide to Long)
plot_data <- df %>%
  # Select the genotype ID and only columns ending with your target trait
  select(Genotype, ends_with(target_trait)) %>%
  
  # Pivot to long format
  pivot_longer(
    cols = -Genotype, 
    names_to = c("Environment", "Trait"),
    names_pattern = "(.*)_(.*)", # Splits "Clavet2024_biomass.g" into Env and Trait
    values_to = "BLUP"
  ) %>%
  
  # Drop NAs just in case some plots were missing
  drop_na(BLUP) 

# 3. Calculate the overall mean per lentil line to sort the X-axis
plot_data <- plot_data %>%
  group_by(Genotype) %>%
  mutate(Mean_BLUP = mean(BLUP)) %>%
  ungroup() %>%
  # Reorder the Lentil factor levels by their Mean_BLUP
  mutate(Genotype = reorder(Genotype, Mean_BLUP))

# 4. Generate the Plot
ggplot(plot_data, aes(x = Genotype, y = BLUP, color = Environment)) +
  geom_point(alpha = 0.7, size = 2) +
  
  # Customizing the look
  theme_bw() +
  theme(
    # Rotate X-axis text 90 degrees and shrink it so 100 names can fit
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(), # Cleans up vertical grid lines
    legend.position = "top"
  ) +
  labs(
    title = paste("Genotype Performance across Environments:", target_trait),
    subtitle = "Genotypes ordered by overall mean performance",
    x = "Lentil Genotype",
    y = paste("BLUP Value (", target_trait, ")", sep = ""),
    color = "Site-Year"
  )

ggsave(paste("Figures/","Predicted mean","_", target_trait,".png", sep=""), width = 14, height = 5, bg = "white", dpi = 600)
