#GWAS analysis####
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT")
library(GAPIT)

#Load myG for GWAS
myG <- read.table("myG_ACT100_filtered.hmp.txt", header = F) #this contains the 100 lines

#Load myY data
setwd("~/SALVA/lentil")
myY_Hunter2024_SpATS <- read.delim("myY_Hunter2024_SpATS.txt")