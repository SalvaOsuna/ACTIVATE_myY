#GWAS analysis####
 #install.packages("devtools")
 #devtools::install_github("jiabowang/GAPIT", force = TRUE)
 #library(GAPIT)
#loading the formulas from source
source("https://raw.githubusercontent.com/jiabowang/GAPIT/refs/heads/master/gapit_functions.txt", encoding = "UTF-8")

setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY/")
#Load myG for GWAS
myG <- read.table("data/myG_ACT100_filtered.hmp.txt", header = F) #this contains the 100 lines

#Load myY data
myY_SpATS <- as.data.frame(read.csv("data/ACTIVATE_lentils_myY_SpATS.csv"))

myY_SpATS <- myY_SpATS |>
  pivot_wider(names_from = c(ENV, Trait), values_from = PredictedMean)

myY_SpATS <- as.data.frame(myY_SpATS)

#Run GWAS
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY/GWAS_results")

myGAPIT <- GAPIT( Y = myY_SpATS[,-c(1)],
                  G = myG,
                  PCA.total = 4,
                  model = c("GLM","BLINK", "MLM", "FarmCPU", "CMLM", "MLMM", "SUPER"),
                  Random.model = F,  # Optional: use if GAPIT returns an error
                  Phenotype.View = F # Optional: use if GAPIT returns an error
)

#Check which traits and models have ran with is.ran()
 #devtools::install_github("derekmichaelwright/gwaspr")
 #install.packages("ggpubr")
library(ggpubr)
library(gwaspr)
is_ran(folder = "GWAS_results/")

#Check if the results files are ordered by p.value with is_ordered()
is_ordered(folder = "GWAS_results/")

order_GWAS_Results(folder = "GWAS_results/")

# Create summary of GWAS results
xx <- table_GWAS_Results(folder = "GWAS_results/", 
                         threshold = 6.8, 
                         sug.threshold = 5.3)

# Plot
mp <- gg_Manhattan_2(
  # Specify a folder with GWAS results
  folder = "GWAS_results/", 
  # Select a trait to plot
  trait = "Clavet_2025_lodging",
  threshold = 6.5
)
ggsave(path = "GWAS_results/","Fig_Clavet2025_lodging.png", mp, width = 12, height = 4, bg = "white")
