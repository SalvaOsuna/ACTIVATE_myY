# ==============================================================================
# PIPELINE: LENTIL GxE, ENVIROMICS, AND GENETIC GAIN
# ==============================================================================

# 1. Load Required Libraries
# install.packages(c("metan", "dplyr", "tidyr", "ggplot2", "nasapower", "pls", "randomForest"))
library(metan)       # Core GxE, AMMI, GGE, and Stability indices
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping
library(ggplot2)     # Visualization
library(nasapower)   # Fetching NASA POWER climate data
library(pls)         # Partial Least Squares for Enviromics

# ==============================================================================
# STAGE 0: DATA PREPARATION
# ==============================================================================

# Assuming you have generated Stage 1 BLUPs for your 100 lines.
# Your dataframe 'lentil_blups' should look like this:
# ENV (Site_Year), GEN (Line), RELEASE_YEAR, YIELD, PHENO_DAYS, etc.

# Example structure based on your trial:
# lentil_blups <- read.csv("lentil_stage1_blups.csv")

# For the sake of the script, here is the expected structure:
# ENV: "Hunter_2024", "Hunter_2025", "Clavet_2024", "Clavet_2025"

lentil_blups <- all_blups_wide |>
  select(genotype, env, trait, SpATS) |>
  pivot_wider(names_from = trait, values_from = SpATS)

myY <- read.csv("data/lentil_blups_wide_spats.csv")

dev_year_df <- myY %>%
  select(Genotype, dev_year) %>%
  mutate(genotype = Genotype) %>%
  distinct()

lentil_blups <- lentil_blups %>%
  left_join(dev_year_df, by = "genotype") |>
  select(!Genotype) |>
  filter(!genotype %in% c("CL1", "CL2", "CL3", "Extra"))

# ==============================================================================
# STAGE 1: GxE STRUCTURE (AMMI, GGE, Finlay-Wilkinson)
# ==============================================================================

# 1.1 AMMI Analysis
# -----------------
# Fit the AMMI model for yield
dat <- dat |>
  filter(!Lentil %in% c("CL1", "CL2", "CL3", "Extra"))

ammi_model <- performs_ammi(dat, 
                            env = ENV, 
                            gen = Lentil, 
                            rep = Rep_gen,
                            resp = YLD,
                            verbose = FALSE)

# Generate AMMI1 Biplot (Main effects vs IPCA1)
p_ammi1 <- plot_scores(ammi_model, type = 1)
ggsave("Results/AMMI1_Biplot_YLD.png", p_ammi1, width = 8, height = 6)

# Generate AMMI2 Biplot (IPCA1 vs IPCA2 to see environmental clustering)
p_ammi2 <- plot_scores(ammi_model, type = 2)
ggsave("Results/AMMI2_Biplot_YLD.png", p_ammi2, width = 8, height = 6)

plot_scores(ammi_model, type = 4)

# 1.2 GGE Biplot Analysis
# -----------------------
# Fit GGE model (Genotype + Genotype-by-Environment)
gge_model <- gge(lentil_blups, env = env, gen = genotype, resp = YLD)

# Which-won-where pattern (Mega-environments)
p_gge_won <- plot(gge_model, type = 3) 
ggsave("Results/GGE_Which_Won_Where_YLD.png", p_gge_won, width = 8, height = 6)

# Mean performance vs stability
p_gge_stab <- plot(gge_model, type = 2)


# 1.3 Finlay-Wilkinson & Stability Statistics (ge_stats)
# ------------------------------------------------------
# metan's ge_stats calculates almost all the stability indices you mentioned at once
# including Wricke's, Eberhart-Russell, ASV, WAASB, and superiority indices.
stability_indices <- ge_stats(dat, 
                              env = ENV, 
                              gen = Lentil, 
                              rep = Rep_gen,
                              resp = YLD)
corr_stab_ind(stability_indices)
# Extract the compiled table of indices
stab_table <- stability_indices$YLD
write.csv(stab_table, "Results/Lentil_Stability_Indices_YLD.csv", row.names = FALSE)

# ==============================================================================
# STAGE 2: CLIMATE DATA & ENVIROMICS
# ==============================================================================

# 2.1 Fetch Climate Data
# ----------------------
# Example: Fetching daily weather for Clavet (approx coords) during the 2024 season.
# You will loop this for Hunter and for 2025 using precise GPS coordinates.
weather_clavet_24 <- get_power(
  community = "ag",
  lonlat = c(-106.229217, 52.115386), # Replace with exact Clavet/Hunter long/lat
  pars = c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN", "RH2M"),
  dates = c("2024-05-15", "2024-08-26"), # Approximate seeding to harvest
  temporal_api = "daily"
) |> mutate( ENV = "Clavet2024")
weather_clavet_25 <- get_power(
  community = "ag",
  lonlat = c(-106.221749, 52.135848), # Replace with exact Clavet/Hunter long/lat
  pars = c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN", "RH2M"),
  dates = c("2025-05-05", "2025-08-29"), # Approximate seeding to harvest
  temporal_api = "daily"
)|> mutate( ENV = "Clavet2025")
weather_hunter_24 <- get_power(
  community = "ag",
  lonlat = c(-106.51472222, 52.20522222), # Replace with exact Clavet/Hunter long/lat
  pars = c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN", "RH2M"),
  dates = c("2024-05-28", "2024-09-09"), # Approximate seeding to harvest
  temporal_api = "daily"
)|> mutate( ENV = "Hunter2024")
weather_hunter_25 <- get_power(
  community = "ag",
  lonlat = c(-106.506301, 52.221332), # Replace with exact Clavet/Hunter long/lat
  pars = c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN", "RH2M"),
  dates = c("2025-05-12", "2025-09-13"), # Approximate seeding to harvest
  temporal_api = "daily"
)|> mutate( ENV = "Hunter2025")

all_weather_data <- rbind(weather_clavet_24, weather_clavet_25, weather_hunter_24, weather_hunter_25)

# 2.2 Derive Crop-Stage Environmental Covariates (averages)
env_phenology <- data.frame(
  ENV = c("Hunter2024", "Clavet2024", "Hunter2025", "Clavet2025"),
  Sow_Date = as.Date(c("2024-05-28", "2024-05-15", "2025-05-12", "2025-05-05")),
  Avg_DTE = c(12.75, 14.77, 16.01, 11.11),  # Days to Emergence
  Avg_DTF = c(47.48, 53.37, 53.94, 55.34),  # Days to Flowering
  Avg_DTM = c(81.56, 87.07, 93.77, 95)   # Days to Maturity
)

# Process the climate data using your real phenology
env_covariates <- all_weather_data %>%
  # 1. Join the phenology averages to the weather data by Environment
  left_join(env_phenology, by = "ENV") %>%
  
  # 2. Calculate the exact days after sowing for every single weather record
  mutate(
    Date = as.Date(YYYYMMDD),
    Days_After_Sowing = as.numeric(Date - Sow_Date),
    
    # 3. Dynamically assign the growth phase using your DTE, DTF, and DTM variables
    Growth_Phase = case_when(
      Days_After_Sowing > 0 & Days_After_Sowing <= Avg_DTE ~ "1_Emergence",
      Days_After_Sowing > Avg_DTE & Days_After_Sowing <= Avg_DTF ~ "2_Vegetative",
      Days_After_Sowing > Avg_DTF & Days_After_Sowing <= Avg_DTM ~ "3_GrainFill",
      TRUE ~ NA_character_ # Ignore days before sowing or after harvest
    ),
    
    # Calculate daily GDD (assuming base temp of 5C)
    GDD = ((T2M_MAX + T2M_MIN) / 2) - 5
  ) %>%
  
  # 4. Remove any dates that fall outside your growing season
  filter(!is.na(Growth_Phase)) %>%
  
  # 5. Group by Environment AND the dynamic Growth Phase to summarize
  group_by(ENV, Growth_Phase) %>%
  summarise(
    Acc_GDD = sum(GDD[GDD > 0], na.rm = TRUE),
    Total_Precip = sum(PRECTOTCORR, na.rm = TRUE),
    Heat_Days = sum(T2M_MAX > 30, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  
  # 6. Pivot the table wide so each environment is exactly one row
  pivot_wider(
    names_from = Growth_Phase, 
    values_from = c(Acc_GDD, Total_Precip, Heat_Days)
  )

# 2.3 Linking Climate to GxE (PLS Regression)
# Extract the environmental IPCA scores from the AMMI model
env_scores <- ammi_model$YLD$MeansGxE %>% 
  select(ENV, envPC1) %>%
  distinct()

# Merge with your derived climate matrix (env_covariates)
enviromics_data <- left_join(env_scores, env_covariates, by = "ENV")

# Check the column names to see exactly what pivot_wider created
colnames(enviromics_data)

# Run Partial Least Squares (PLS) to see which climate variables drive envPC1
# We will use the exact column names generated by your previous step.
# Note: With only 4 environments, we remove validation = "CV" because 
# cross-validation requires a larger sample size.
pls_model <- plsr(envPC1 ~ Total_Precip_2_Vegetative + 
                    Total_Precip_3_GrainFill+
                    Acc_GDD_2_Vegetative + 
                    Acc_GDD_3_GrainFill + 
                    Heat_Days_3_GrainFill, 
                  data = enviromics_data, 
                  scale = TRUE)

summary(pls_model)

# Visualize how the climate vectors pull the environments along the GxE axis
biplot(pls_model)
validationplot(pls_model)
validationplot(pls_model, val.type="MSEP")
validationplot(pls_model, val.type="R2")

# ==============================================================================
# STAGE 3: GENETIC GAIN IN STABILITY
# ==============================================================================

# Merge stability indices with your line metadata (Release Year)
# metadata <- read.csv("lentil_line_metadata.csv") # Contains GEN and RELEASE_YEAR

myY <- read.csv("data/lentil_blups_wide_spats.csv")

dev_year_df <- myY %>%
  select(Genotype, dev_year) %>%
  mutate(GEN = Genotype) %>%
  distinct()

gain_data <- stab_table %>%
  left_join(dev_year_df, by = "GEN")

# 3.1 Plot WAASB (Weighted AMMI Stability) vs Release Year
# Note: Lower WAASB means higher stability.
p_gain_stab <- ggplot(gain_data, aes(x = dev_year, y = WAASB)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue") +
  theme_minimal() +
  labs(title = "Genetic Gain: Stability over Decades",
       x = "Year of Release",
       y = "WAASB (Lower = More Stable)")

ggsave("Genetic_Gain_WAASB.png", p_gain_stab, width = 7, height = 5)

# ==============================================================================
# 3.2a COMBINING YIELD AND STABILITY (WAASBY INDEX)
# ==============================================================================

# 1. Fit the WAASB model (Mixed-model version of AMMI)
waasb_model <- waasb(dat, 
                     env = ENV, 
                     gen = Lentil, 
                     rep = Rep_gen,
                     resp = YLD,
                     verbose = FALSE)

# 2. Calculate the WAASBY index
# wresp = 50 means we give 50% weight to mean Yield and 50% weight to Stability.
# You can adjust this (e.g., wresp = 65 for 65% yield / 35% stability)
waasby_index <- wsmp(waasb_model)

# 3. Generate the Yield x Stability Biplot (Highly recommended for the paper)
# Genotypes in the bottom-right quadrant (Q4) are highly productive and highly stable.
p_waasby_biplot <- plot(waasby_index, type = 1)
ggsave("WAASB_vs_Yield_Biplot.png", p_waasby_biplot, width = 8, height = 6)

# 4. Link WAASBY scores back to Release Year for the Genetic Gain regression
# Extract the calculated WAASBY scores
waasby_scores <- waasby_index$YLD$WAASY %>% select(Code, WAASBY)

# Merge with your metadata (Release Year)
# Assuming metadata object exists: metadata <- read.csv("lentil_line_metadata.csv")
myY <- read.csv("data/lentil_blups_wide_spats.csv")
dev_year_df <- myY %>%
  select(Genotype, dev_year) %>%
  mutate(Code = Genotype) %>%
  distinct()
gain_data_waasby <- waasby_scores %>%
  left_join(dev_year_df, by = "Code")

# Plot Genetic Gain of the combined index
# Note: For WAASBY, a HIGHER score is better (unlike raw WAASB where lower is better)
p_gain_combo <- ggplot(gain_data_waasby, aes(x = dev_year, y = WAASBY)) +
  geom_point(alpha = 0.6, color = "darkgray") +
  geom_smooth(method = "lm", color = "forestgreen", fill = "lightgreen") +
  theme_minimal() +
  labs(title = "Simultaneous Genetic Gain: Yield & Stability",
       subtitle = "WAASBY Index (50/50 Weighting)",
       x = "Year of Release",
       y = "WAASBY Score (Higher = Better)")

ggsave("Genetic_Gain_WAASBY.png", p_gain_combo, width = 7, height = 5)

# ==============================================================================
# 3.2b THE IDEOTYPE DISTANCE (MGIDI INDEX)
# ==============================================================================

# 1. Prepare a Genotype x Trait matrix
# We will use Mean Yield (Y) and Stability (WAASB) from your earlier stab_table,
# and pretend we are joining it with a dataframe of average Days to Maturity (DTM).
# (Higher Y is better, Lower WAASB is better, Lower DTM is better)

ideotype_df <- stab_table %>%
  select(GEN, Y, WAASB) %>%
  # left_join(mean_dtm_data, by = "GEN") %>% # Uncomment and join your DTM data here
  column_to_rownames("GEN") # MGIDI requires genotypes to be row names!

# 2. Run the MGIDI Index
# The ideotype argument tells the model what the "perfect" lentil looks like.
# "h" = higher value is better, "l" = lower value is better
# Assuming 3 traits: Yield (h), WAASB (l), DTM (l)
# If you only use Yield and WAASB, use ideotype = c("h", "l")
mgidi_model <- mgidi(ideotype_df, ideotype = c("h", "l")) 

# 3. Visualize the MGIDI index
# This plots the distance of each line from the ideotype. The selected lines are red.
p_mgidi_radar <- plot(mgidi_model, SI = 10)
ggsave("MGIDI_Trait_Contribution.png", p_mgidi_radar, width = 8, height = 6)

# 4. Link MGIDI distance to Release Year
# Extract the MGIDI index scores
mgidi_scores <- mgidi_model$MGIDI %>%
  select(Genotype, MGIDI) %>%
  rename(GEN = Genotype)
write.csv(mgidi_scores, "Results/mgidi_scores_YLD.csv" ,row.names = FALSE)
mgidi_scores <- read.csv("Results/mgidi_scores_YLD.csv") # I have change the genotype nomenclature manually

dev_year_df <- myY %>%
  select(Genotype, dev_year) %>%
  mutate(GEN = Genotype) %>%
  distinct()

gain_data_mgidi <- mgidi_scores %>%
  left_join(dev_year_df, by = "GEN")

# Plot Genetic Gain for the Ideotype
# Note: For MGIDI, a LOWER score is better (it means the line is "closer" to the perfect ideotype)
p_gain_ideotype <- ggplot(gain_data_mgidi, aes(x = dev_year, y = MGIDI)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "darkred") +
  theme_minimal() +
  labs(title = "Genetic Gain: Distance from Ideotype",
       subtitle = "MGIDI Index (Lower = Closer to Perfect Lentil)",
       x = "Year of Release",
       y = "MGIDI Distance")

ggsave("Genetic_Gain_MGIDI.png", p_gain_ideotype, width = 7, height = 5)