# =============================================================================
# Agronomic Performance, GxE Interaction & Enviromics Analysis
# ACTIVATE Lentil Rotation Trial — 100 lines × 4 environments
# =============================================================================
# Pipeline:
#   1. Load BLUPs from genetic gain script (Stage 1 outputs)
#   2. AMMI analysis (IPCA decomposition)
#   3. GGE biplot (which-won-where)
#   4. Finlay-Wilkinson stability regression
#   5. Stability indices (WAASB, ASV, ecovalence, MGIDI)
#   6. Climate data download & processing (NASA POWER + Environment Canada)
#   7. Crop-stage environmental covariates (GDD, precip, VPD, etc.)
#   8. Environmental kinship matrix (W matrix)
#   9. Enviromics: PLS linking AMMI IPCAs to climate variables
#  10. Correlation: stability vs. release year
#  11. Publication figures
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────

# Install any missing packages before loading:
# install.packages(c("metan","agricolae","nasapower","tidyverse",
#                    "factoextra","ggrepel","corrplot","pls",
#                    "randomForest","vegan","reshape2","scales"))

library(metan)        # GGE, AMMI, WAASB, MGIDI, stability indices
library(agricolae)    # AMMI (alternative), Finlay-Wilkinson
library(nasapower)    # NASA POWER climate data API
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggrepel)      # Non-overlapping labels
library(factoextra)   # PCA visualisation
library(corrplot)     # Correlation matrices
library(pls)          # Partial least squares
library(reshape2)     # Data reshaping
library(scales)       # Plot scales
library(viridis)      # Color palettes
library(lubridate)    # Date handling for climate data


# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
# ── 1a. Raw trial data (same file as genetic gain script) --------------------
dat <- read.csv("data/ACTIVATE_lentils_myYraw.csv", stringsAsFactors = FALSE)

#add a replication (gen_rep) col:
dat <- dat |>
  # 1. Group by Environment AND Genotype
  # This ensures the counter restarts for every new variety in every new site
  group_by(ENV, Lentil) |>
  # 2. Add the Rep column
  # row_number() simply assigns 1, 2, 3... to the rows in that group
  mutate(Rep_gen = row_number(),) |>
  # 3. Ungroup (Important! removes the grouping structure so it doesn't mess up future stats)
  ungroup() |>
  # Optional: Move 'Rep' to be near the Genotype column for easier viewing
  relocate(Rep_gen, .after = Lentil) |>
  #rename Rep by Rep_combo
  rename(Rep_combo = Rep)%>%
  mutate(
    genotype   = as.factor(Lentil),
    site       = as.factor(site),
    year       = as.factor(year),
    rep        = as.factor(Rep_gen),
    block      = as.factor(Block),
    env        = interaction(site, year, drop = TRUE)
  )

# ── 1b. Load Stage-1 BLUPs (output from genetic_gain_analysis.R) -------------
# These are environment-specific BLUPs (best model per trait × env)
blups_env <- read.csv("Results/combined_BLUPs_best_model.csv")  # genotype, BLUP_combined, trait
# Also load the per-environment BLUPs (needed for GxE matrix)
blups_SpATS <- read.csv("Results/blups_SpATS_all.csv", stringsAsFactors = FALSE)
blups_GIBD  <- read.csv("Results/blups_GIBD_all.csv",  stringsAsFactors = FALSE)
best_model  <- read.csv("Results/best_model_selection.csv",     stringsAsFactors = FALSE)

# ── 1c. Build per-environment BLUP matrix (genotype × environment) -----------
# Select best-model BLUPs per trait × env
blups_best_long <- bind_rows(
  blups_SpATS %>% semi_join(best_model %>% filter(best_model == "SpATS"),
                            by = c("trait", "env")),
  blups_GIBD  %>% semi_join(best_model %>% filter(best_model == "GIBD"),
                            by = c("trait", "env"))
)

# ── 1d. Release year lookup ---------------------------------------------------
myY <- read.csv("data/lentil_blups_wide_spats.csv")
release_yr <- myY %>%
  select(Genotype, dev_year) %>%
  distinct()

# ── 1e. Traits & environments ------------------------------------------------
traits <- unique(blups_best_long$trait)
envs   <- unique(blups_best_long$env)
cat("Traits :", paste(traits, collapse=", "), "\n")
cat("Envs   :", paste(envs,   collapse=", "), "\n")


# ── 2. AMMI ANALYSIS ──────────────────────────────────────────────────────────
# Uses metan::performs_ammi() — works on a two-way table (genotype × env)
# Requires: genotype, env, rep, trait columns in long format

ammi_results <- list()

for (tr in traits) {
  cat("\nAMMI —", tr, "\n")

  # metan needs the raw data (not BLUPs) for AMMI with reps
  d_tr <- dat %>%
    select(genotype, env, rep, value = all_of(tr)) %>%
    filter(!is.na(value))

  tryCatch({
    ammi_mod <- performs_ammi(
      .data      = d_tr,
      env        = env,
      gen        = genotype,
      rep        = rep,
      resp       = value,
      verbose    = FALSE
    )
    ammi_results[[tr]] <- ammi_mod

    # AMMI stability value (ASV) — distance from origin in IPCA1 × IPCA2 space
    ammi_stab <- ammi_indexes(ammi_mod)
    write.csv(ammi_stab$value, paste0("Results/AMMI_stability_", tr, ".csv"), row.names = FALSE)

  }, error = function(e) message("AMMI failed for ", tr, ": ", e$message))
}

# ── 3. GGE BIPLOT ─────────────────────────────────────────────────────────────
# gge() in metan — shows which-won-where and mega-environment structure

gge_results <- list()

for (tr in traits) {
  cat("\nGGE —", tr, "\n")

  d_tr <- dat %>%
    select(genotype, env, rep, value = all_of(tr)) %>%
    filter(!is.na(value))

  tryCatch({
    gge_mod <- gge(
      .data  = d_tr,
      env    = env,
      gen    = genotype,
      rep    = rep,
      resp   = value,
      svp       = "symmetrical"
    )
    gge_results[[tr]] <- gge_mod

    # Export GGE scores for custom plotting
    gge_scores <- gge_mod[["value"]][["coordgen"]]
    write.csv(gge_scores, paste0("Results/GGE_scores_", tr, ".csv"), row.names = FALSE)

  }, error = function(e) message("GGE failed for ", tr, ": ", e$message))
}


# ── 4. FINLAY-WILKINSON STABILITY REGRESSION ──────────────────────────────────
# Correct approach: metan does not have stability_measures().
# FW is implemented manually via OLS regression of each genotype's BLUPs
# on the environment index (env mean − grand mean), which is the exact
# Finlay & Wilkinson (1963) definition.
#
# Interpretation of slope b_i:
#   b_i = 1   → average stability (responds like the average genotype)
#   b_i < 1   → below-average response (stable; suits low-input environments)
#   b_i > 1   → above-average response (specifically adapted to good environments)
#
# s²d (deviation MS) measures unpredictable instability beyond the regression.

fw_results <- list()

finlay_wilkinson <- function(df, min_Ij_sd = 0.01, max_b = 10) {
  # df must have columns: genotype, env, BLUP, trait
  #
  # Guards against degenerate slopes:
  #   min_Ij_sd : minimum SD of environment indices required to attempt
  #               regression. If the four environments produce nearly
  #               identical trait means (SD of I_j < this threshold),
  #               the slope is mathematically undefined — skip the trait.
  #   max_b     : slopes with |b| > max_b are flagged as unreliable
  #               (numerical artefact of near-zero I_j variance) and
  #               set to NA so they don't pollute the plot or summaries.
  
  grand_mean <- mean(df$BLUP, na.rm = TRUE)
  
  # Environment index I_j = mean of all genotypes in env j minus grand mean
  env_index <- df %>%
    group_by(env) %>%
    summarise(I_j = mean(BLUP, na.rm = TRUE) - grand_mean, .groups = "drop")
  
  # Guard: if environments barely differ for this trait, FW is meaningless
  Ij_sd <- sd(env_index$I_j, na.rm = TRUE)
  if (is.na(Ij_sd) || Ij_sd < min_Ij_sd) {
    message("  [FW] Skipping — environment index SD = ", round(Ij_sd, 6),
            " (< ", min_Ij_sd, "). Trait has negligible E main effect; ",
            "FW regression is not estimable.")
    return(NULL)
  }
  
  df    <- df %>% left_join(env_index, by = "env")
  n_env <- length(unique(df$env))
  
  results <- df %>%
    group_by(genotype) %>%
    group_modify(function(g, key) {
      if (sum(!is.na(g$BLUP)) < 3) {
        return(data.frame(mean_BLUP = NA, b = NA, s2d = NA, r2 = NA,
                          p_slope = NA, reliable = FALSE,
                          stringsAsFactors = FALSE))
      }
      lm_fw <- lm(BLUP ~ I_j, data = g)
      sm    <- summary(lm_fw)
      b     <- coef(lm_fw)[["I_j"]]
      r2    <- sm$r.squared
      p_slp <- coef(sm)[2, 4]
      s2d   <- sum(resid(lm_fw)^2) / max(1, n_env - 2)
      
      # Flag slopes that exploded due to residual near-zero I_j variance
      reliable <- !is.na(b) && abs(b) <= max_b
      
      data.frame(mean_BLUP = mean(g$BLUP, na.rm = TRUE),
                 b         = ifelse(reliable, b, NA_real_),
                 s2d       = s2d,
                 r2        = ifelse(reliable, r2,    NA_real_),
                 p_slope   = ifelse(reliable, p_slp, NA_real_),
                 reliable  = reliable,
                 stringsAsFactors = FALSE)
    }) %>%
    ungroup() %>%
    mutate(
      trait   = unique(df$trait),
      b_class = case_when(
        is.na(b) ~ "Unreliable",
        b <  0.8 ~ "Stable / low-input adapted",
        b >  1.2 ~ "Responsive / high-input adapted",
        TRUE     ~ "Average stability"
      )
    )
  
  n_unreliable <- sum(!results$reliable, na.rm = TRUE)
  if (n_unreliable > 0)
    message("  [FW] ", n_unreliable, " genotype(s) had |b| > ", max_b,
            " and were set to NA.")
  
  return(results)
}

for (tr in traits) {
  cat("\nFinlay-Wilkinson —", tr, "\n")
  
  df_tr <- blups_best_long %>%
    filter(trait == tr) %>%
    select(genotype, env, BLUP, trait) %>%
    filter(!is.na(BLUP))
  
  tryCatch({
    fw_df <- finlay_wilkinson(df_tr)
    fw_results[[tr]] <- fw_df
    
    cat("  b range: [", round(min(fw_df$b, na.rm=TRUE), 3), ",",
        round(max(fw_df$b, na.rm=TRUE), 3), "]\n")
    cat("  Genotypes with b > 1.2 (responsive):",
        sum(fw_df$b > 1.2, na.rm=TRUE), "\n")
    cat("  Genotypes with b < 0.8 (stable)    :",
        sum(fw_df$b < 0.8, na.rm=TRUE), "\n")
    
    write.csv(fw_df, paste0("FW_stability_", tr, ".csv"), row.names = FALSE)
    
  }, error = function(e) message("FW failed for ", tr, ": ", e$message))
}

# Combine all FW results
fw_all <- bind_rows(fw_results)
write.csv(fw_all, "Results/FW_stability_all_traits.csv", row.names = FALSE)
cat("\nFinlay-Wilkinson complete for", length(fw_results), "traits\n")


# ── 5. COMPREHENSIVE STABILITY INDICES (metan) ────────────────────────────────
# Computes: WAASB, ASV, ecovalence (Wricke), s²d (Eberhart-Russell),
#           superiority index (Pi), CV, HMGV, RPGV

stab_all <- list()

for (tr in traits) {
  cat("\nStability indices —", tr, "\n")
  
  d_tr <- dat %>%
    select(genotype, env, rep, value = all_of(tr)) %>%
    filter(!is.na(value))
  
  tryCatch({
    # WAASB = Weighted Average of Absolute Scores from mixed model
    waasb_mod <- waasb(
      .data   = d_tr,
      env     = env,
      gen     = genotype,
      rep     = rep,
      resp    = value,
      verbose = FALSE
    )
    
    # All parametric and non-parametric stability statistics
    stab_mod <- ge_stats(
      .data   = d_tr,
      env     = env,
      gen     = genotype,
      rep     = rep,
      resp    = value,
      verbose = FALSE
    )
    
    stab_all[[tr]] <- list(waasb = waasb_mod, stats = stab_mod)
    write.csv(stab_mod$value, paste0("stability_indices_", tr, ".csv"),
              row.names = FALSE)
    
  }, error = function(e) message("Stability failed for ", tr, ": ", e$message))
}

# ── 5b. MGIDI — Multi-trait Genotype-Ideotype Distance Index -----------------
# mgidi() requires a single data frame where:
#   - First column  : genotype identifier (GEN), character or factor
#   - All other cols: numeric traits/stability indices (one column per variable)
# It does NOT accept a list of waasb objects directly.
#
# Strategy: extract the per-genotype mean BLUP and WAASB score for each trait,
# then join everything into one wide data frame with one row per genotype.
# Columns = trait means + WAASB scores (labelled e.g. yield_mean, yield_WAASB).
# For traits where *higher is better* leave as-is; for traits where
# *lower is better* (e.g. days to maturity, WAASB) pass them to the
# `ideotype` argument as minimisation targets.
mgidi_frames <-NULL
tryCatch({
  
  # ── Extract BOTH mean (Y) and WAASB score simultaneously -------------------
  mgidi_frames <- imap(stab_all, function(x, tr) {
    
    # 1. Skip if the trait failed in the previous loop
    if (is.null(x$waasb)) {
      cat("Skipping", tr, "- waasb object is missing\n")
      return(NULL)
    }
    
    # 2. Print success message to console
    cat("Successfully unpacked nested path for:", tr, "\n")
    
    # 3. Navigate the exact path you mapped out: stab_all[[tr]][["waasb"]][["value"]][["model"]]
    df <- x$waasb$value$model %>%
      as.data.frame() %>%
      select(Code, Y, WAASB) %>%
      #mutate(GEN = Code) |>
      rename(
        !!paste0(tr, "_mean") := Y,
        !!paste0(tr, "_WAASB") := WAASB,
        GEN = Code
      )
    
    return(df)
  })
  
  # Combine them into your wide dataframe for MGIDI
  all_frames <- Filter(Negate(is.null), mgidi_frames)
  mgidi_input <- reduce(all_frames, full_join, by = "GEN")
  
  cat("\nMGIDI input dimensions:", nrow(mgidi_input), "genotypes ×",
      ncol(mgidi_input) - 1, "variables\n")
  
  # ── Define ideotype direction per column -----------------------------------
  # "max" = higher is better (yield, seed weight)
  # "min" = lower is better (days to maturity, WAASB scores)
  ideotype_dir <- ifelse(
    grepl("days_|DTF|DTE|DTM|WAASB", names(mgidi_input)[-1], ignore.case = TRUE),
    "min", "max"
  )
  ideotype_str <- paste(ideotype_dir, collapse = ", ")
  
  # ── Fit MGIDI --------------------------------------------------------------
  mgidi_mod <- mgidi(
    .data    = mgidi_input,
    SI       = 15,            # 15% selection intensity
    ideotype = ideotype_str,
    verbose  = FALSE
  )
  
  # ── Extract and save results -----------------------------------------------
  mgidi_scores <- mgidi_mod$MGIDI %>%
    arrange(MGIDI) %>%
    inner_join(release_yr, by = "Genotype")
  
  write.csv(mgidi_scores, "Results/MGIDI_multi_trait_ranking.csv", row.names = FALSE)
  write.csv(mgidi_mod$sel_dif, "Results/MGIDI_selection_differential.csv", row.names = FALSE)
  
  cat("MGIDI ranking complete —", nrow(mgidi_scores), "genotypes ranked\n")
  cat("Top 10 selected genotypes:\n")
  print(head(mgidi_scores, 10))
  
}, error = function(e) {
  message("MGIDI failed: ", e$message)
})

# ── 6. CLIMATE DATA DOWNLOAD (NASA POWER) ─────────────────────────────────────
# Replace site coordinates with your actual site lat/lon values

site_info <- data.frame(
  site     = c("Clavet", "Hunter"),            # Replace with actual site names
  lat      = c(52.115386, 52.221332),              # Replace with actual latitudes
  lon      = c(-106.229217, -106.506301),            # Replace with actual longitudes
  stringsAsFactors = FALSE
)

trial_years <- c(2024, 2025)                 # Replace with actual trial years

# Variables to download from NASA POWER
# T2M = temp 2m, T2M_MAX, T2M_MIN, PRECTOTCORR = precip,
# ALLSKY_SFC_SW_DWN = solar radiation, RH2M = relative humidity,
# WS2M = wind speed
power_vars <- c("T2M", "T2M_MAX", "T2M_MIN",
                "PRECTOTCORR", "ALLSKY_SFC_SW_DWN",
                "RH2M", "WS2M")

climate_raw <- list()

for (s in site_info$site) {
  for (yr in trial_years) {
    key  <- paste(s, yr, sep = "_")
    lat  <- site_info$lat[site_info$site == s]
    lon  <- site_info$lon[site_info$site == s]
    cat("Downloading NASA POWER:", key, "\n")

    tryCatch({
      pw <- get_power(
        community  = "AG",                  # Agroclimatology
        lonlat     = c(lon, lat),
        pars       = power_vars,
        dates      = c(paste0(yr, "-01-01"), paste0(yr, "-12-31")),
        temporal_api = "daily"
      )
      pw$site <- s
      pw$trial_year <- yr
      climate_raw[[key]] <- pw
      Sys.sleep(1)   # polite pause between API calls
    }, error = function(e) message("NASA POWER download failed for ", key, ": ", e$message))
  }
}

# Combine all climate data
climate_daily <- bind_rows(climate_raw)
write.csv(climate_daily, "Results/climate_daily_NASA_POWER.csv", row.names = FALSE)
cat("Climate data saved:", nrow(climate_daily), "daily records\n")


# ── 7. CROP-STAGE ENVIRONMENTAL COVARIATES ────────────────────────────────────
# Define approximate phenological windows per environment
# Adjust these dates to match your actual sowing / harvest calendar

#dates
pheno_windows <- data.frame(
  site        = c("Clavet","Clavet","Hunter","Hunter"),
  trial_year  = c(2024,  2025,  2024,  2025),
  sowing      = as.Date(c("2024-05-15","2025-05-05","2024-05-28","2025-05-12")),
  harvest     = as.Date(c("2024-08-26","2025-08-29","2024-09-09","2025-09-13")),
  Avg_DTE = c(12.75, 14.77, 16.01, 11.11),  # Days to Emergence
  Avg_DTF = c(47.48, 53.37, 53.94, 55.34),  # Days to Flowering
  Avg_DTM = c(81.56, 87.07, 93.77, 95)   # Days to Maturity
) |>
  mutate(
    emergence   = as.Date(sowing+Avg_DTE),
    flowering   = as.Date(sowing+Avg_DTF),
    maturity    = as.Date(sowing+Avg_DTM),
  ) |>
  select(!c(Avg_DTE, Avg_DTF,Avg_DTM))

# ── 7a. Helper: compute covariates for a given growth window ------------------
compute_covariates <- function(clim_df, start_date, end_date, T_base = 5) {
  d <- clim_df %>%
    filter(YYYYMMDD >= start_date & YYYYMMDD <= end_date)

  if (nrow(d) == 0) return(NULL)

  # Growing degree days (GDD = ((Tmax + Tmin)/2 - Tbase), floor at 0)
  GDD <- sum(pmax((d$T2M_MAX + d$T2M_MIN) / 2 - T_base, 0), na.rm = TRUE)

  # Precipitation total (mm)
  precip_total <- sum(d$PRECTOTCORR, na.rm = TRUE)

  # Mean solar radiation (MJ m-2 d-1)
  mean_rad <- mean(d$ALLSKY_SFC_SW_DWN, na.rm = TRUE)

  # Mean VPD (kPa) — computed from Tmax, Tmin, RH2M
  # Tetens formula: es = 0.6108 * exp(17.27*T / (T+237.3))
  es_max <- 0.6108 * exp(17.27 * d$T2M_MAX / (d$T2M_MAX + 237.3))
  es_min <- 0.6108 * exp(17.27 * d$T2M_MIN / (d$T2M_MIN + 237.3))
  es     <- (es_max + es_min) / 2
  ea     <- es * (d$RH2M / 100)
  VPD    <- mean(es - ea, na.rm = TRUE)

  # Heat stress days (Tmax > 30°C)
  heat_days <- sum(d$T2M_MAX > 30, na.rm = TRUE)

  # Frost days (Tmin < 0°C)
  frost_days <- sum(d$T2M_MIN < 0, na.rm = TRUE)

  # Diurnal temperature range
  dtr <- mean(d$T2M_MAX - d$T2M_MIN, na.rm = TRUE)

  # Mean temperature
  tmean <- mean(d$T2M, na.rm = TRUE)

  data.frame(GDD, precip_total, mean_rad, VPD,
             heat_days, frost_days, dtr, tmean)
}

# ── 7b. Compute covariates for each growth phase × environment ----------------
covariate_list <- list()

for (i in seq_len(nrow(pheno_windows))) {
  pw_row  <- pheno_windows[i, ]
  env_key <- paste(pw_row$site, pw_row$trial_year, sep = ".")
  clim    <- climate_daily %>%
    filter(site == pw_row$site, trial_year == pw_row$trial_year) %>%
    mutate(YYYYMMDD = as.Date(as.character(YYYYMMDD), format = "%Y-%m-%d"))

  phases <- list(
    sowing_to_emergence = c(pw_row$sowing,    pw_row$emergence),
    emergence_to_flower = c(pw_row$emergence, pw_row$flowering),
    flowering_to_maturity = c(pw_row$flowering, pw_row$maturity),
    full_season         = c(pw_row$sowing,    pw_row$harvest)
  )

  env_covs <- map_dfr(names(phases), function(ph) {
    cov <- compute_covariates(clim, phases[[ph]][1], phases[[ph]][2])
    if (!is.null(cov)) {
      cov$phase <- ph
      cov$env   <- env_key
    }
    cov
  })

  covariate_list[[env_key]] <- env_covs
}

cov_long <- bind_rows(covariate_list)
write.csv(cov_long, "Results/climate_covariates_by_phase.csv", row.names = FALSE)

# ── 7c. Wide format: one row per environment, columns = phase_variable --------
cov_wide <- cov_long %>%
  pivot_wider(names_from  = phase,
              values_from = c(GDD, precip_total, mean_rad, VPD,
                               heat_days, frost_days, dtr, tmean),
              names_glue  = "{phase}_{.value}") %>%
  column_to_rownames("env")

write.csv(cov_wide, "Results/climate_covariates_wide.csv")
cat("Environmental covariates matrix:", nrow(cov_wide), "envs ×",
    ncol(cov_wide), "variables\n")

# ── 8. ENVIRONMENTAL KINSHIP MATRIX (W matrix) ────────────────────────────────
# Quantifies climatic similarity between environments
# Analogous to genomic relationship matrix but for environments

# Standardise covariates (z-scores)
cov_scaled <- scale(cov_wide)
cov_scaled[is.nan(cov_scaled)] <- 0

# W = (1/p) * Z %*% t(Z)  where Z = scaled covariate matrix, p = n covariates
W_matrix <- tcrossprod(cov_scaled) / ncol(cov_scaled)
rownames(W_matrix) <- colnames(W_matrix) <- rownames(cov_wide)

write.csv(as.data.frame(W_matrix), "Results/W_environmental_kinship.csv")

# Convert W to correlation-scale matrix (diagonal = 1)
D_inv_sqrt <- diag(1 / sqrt(diag(W_matrix)))
W_cor <- D_inv_sqrt %*% W_matrix %*% D_inv_sqrt
rownames(W_cor) <- colnames(W_cor) <- rownames(W_matrix)

write.csv(as.data.frame(W_cor), "Results/W_cor_environmental_kinship.csv")

# Visualise W matrix
png("Results/plot_W_matrix_heatmap.png", width = 1200, height = 1000, res = 180)
corrplot::corrplot(W_cor,
         method  = "color",
         type    = "upper",
         addCoef.col = "black",
         tl.col  = "black",
         col     = colorRampPalette(c("#D7191C","white","#2C7BB6"))(100),
         title   = "Environmental Kinship Matrix (W) — Climate Covariates",
         mar     = c(0,0,2,0))
dev.off()

# PCA of environments based on W
 #remove columns with 0 var
# Calculate the variance of every column
col_variances <- apply(cov_scaled, 2, var, na.rm = TRUE)
# Print the names of the columns where variance is exactly 0 (or NA)
names(col_variances[col_variances == 0 | is.na(col_variances)])
# Keep only columns where the variance is greater than 0
cov_filtered <- as.data.frame(cov_scaled) %>% 
  dplyr::select(where(~ var(., na.rm = TRUE) > 0))
# Now run your PCA on the filtered dataset
env_pca <- prcomp(cov_filtered, scale. = FALSE)
png("Results/plot_environment_PCA.png", width = 1400, height = 1000, res = 180)
fviz_pca_biplot(env_pca,
                repel      = TRUE,
                col.ind    = "#2C7BB6",
                col.var    = "#D7191C",
                title      = "Environment PCA — Climate Covariates") +
  theme_bw(base_size = 11)
dev.off()


# ── 9. ENVIROMICS: PLS LINKING AMMI IPCAs TO CLIMATE VARIABLES ────────────────
# For each trait: extract AMMI IPCA scores per environment,
# then use PLS to identify which climate covariates explain them

pls_results <- list()

for (tr in traits) {
  cat("\nPLS Enviromics —", tr, "\n")

  if (is.null(ammi_results[[tr]])) next

  tryCatch({
    # Extract IPCA scores per environment from AMMI
    ammi_env_scores <- ammi_results[[tr]]$value$model %>%
      filter(type == "ENV") %>%
      select(Code, PC1, PC2) %>%
      rename(env = Code)

    # Match with climate covariates
    cov_for_pls <- cov_wide[ammi_env_scores$env, , drop = FALSE]

    if (nrow(cov_for_pls) < 3) {
      message("  Not enough environments for PLS (n=", nrow(cov_for_pls), ")")
      next
    }

    Y_pls <- as.matrix(ammi_env_scores[, c("PC1","PC2")])
    X_pls <- as.matrix(cov_for_pls)

    # Remove zero-variance columns
    X_pls <- X_pls[, apply(X_pls, 2, var) > 0, drop = FALSE]

    pls_mod <- plsr(Y_pls ~ X_pls,
                    ncomp     = min(2, nrow(X_pls) - 1),
                    scale     = TRUE,
                    validation = "LOO")

    pls_results[[tr]] <- pls_mod

    # Variable importance (VIP scores proxy: loading weights)
    loadings_pls <- as.data.frame(loadings(pls_mod)[, 1:2])
    loadings_pls$variable <- rownames(loadings_pls)
    loadings_pls$trait    <- tr
    write.csv(loadings_pls, paste0("Results/PLS_loadings_", tr, ".csv"), row.names = FALSE)

    cat("  R² comp1:", round(R2(pls_mod)$val[1,1,1], 3),
        "| comp2:", round(R2(pls_mod)$val[1,1,2], 3), "\n")

  }, error = function(e) message("PLS failed for ", tr, ": ", e$message))
}

# ── 9b. Random Forest variable importance (complement to PLS) ----------------
library(randomForest)

rf_importance_list <- list()

for (tr in traits) {
  cat("\nRandom Forest importance —", tr, "\n")

  if (is.null(ammi_results[[tr]])) next

  tryCatch({
    ammi_env_scores <- ammi_results[[tr]]$value$model %>%
      filter(type == "ENV") %>%
      select(Code, PC1) %>%
      rename(env = Code)

    cov_rf <- cov_wide[ammi_env_scores$env, , drop = FALSE]
    cov_rf  <- cov_rf[, apply(cov_rf, 2, var) > 0, drop = FALSE]

    if (nrow(cov_rf) < 4) next

    df_rf <- data.frame(IPCA1 = ammi_env_scores$PC1, cov_rf)
    rf_mod <- randomForest(IPCA1 ~ ., data = df_rf,
                           ntree = 500, importance = TRUE)

    imp_df <- data.frame(
      variable  = rownames(importance(rf_mod)),
      IncMSE    = importance(rf_mod)[, "%IncMSE"],
      IncNodePur = importance(rf_mod)[, "IncNodePurity"],
      trait     = tr
    ) %>% arrange(desc(IncMSE))

    rf_importance_list[[tr]] <- imp_df
    write.csv(imp_df, paste0("RF_importance_", tr, ".csv"), row.names = FALSE)

  }, error = function(e) message("RF failed for ", tr, ": ", e$message))
}


# ── 10. STABILITY vs. RELEASE YEAR CORRELATION ────────────────────────────────
# Key question: Are newer lines more stable AND higher yielding?

stab_vs_time <- list()

for (tr in traits) {
  if (is.null(stab_all[[tr]])) next
  
  tryCatch({
    stab_df <- stab_all[[tr]]$stats$value %>%
      as.data.frame() %>%
      mutate(genotype = as.character(GEN)) %>%
      inner_join(release_yr, by = "genotype")
    
    # FIX 1: Use the exact column names output by metan ("s2di" instead of "S2")
    stab_vars <- c("WAASB", "ASV", "Wi_g", "s2di", "Pi_a")
    stab_vars_avail <- intersect(stab_vars, colnames(stab_df))
    
    cors <- map_dfr(stab_vars_avail, function(sv) {
      vals  <- as.numeric(stab_df[[sv]])
      
      # FIX 2: Safely extract the year column (checking for both common names)
      ryear <- if("dev_year" %in% colnames(stab_df)) stab_df$dev_year else stab_df$release_yr
      
      ok    <- !is.na(vals) & !is.na(ryear)
      
      if (sum(ok) < 5) return(NULL)
      
      ct <- cor.test(ryear[ok], vals[ok], method = "spearman")
      
      data.frame(trait = tr, 
                 stability_index = sv,
                 rho   = round(ct$estimate, 3),
                 p_val = round(ct$p.value,  4))
    })
    
    stab_vs_time[[tr]] <- list(data = stab_df, correlations = cors)
    cat("\n── Stability × Release Year correlations (", tr, ") ──\n")
    print(cors)
    
  }, error = function(e) message("Stab×Time failed for ", tr, ": ", e$message))
}

stab_time_cors <- bind_rows(lapply(stab_vs_time, function(x) x$correlations))
write.csv(stab_time_cors, "Results/stability_vs_release_year.csv", row.names = FALSE)


# ── 11. PUBLICATION FIGURES ───────────────────────────────────────────────────

# ── 11a. AMMI1 biplot (PC1 score vs. genotype mean) ─────────────────────────
for (tr in traits) {
  if (is.null(ammi_results[[tr]])) next
  tryCatch({
    p <- plot_scores(ammi_results[[tr]],
                     type = 1,        # AMMI1: mean vs IPCA1
                     #title = paste("AMMI1 Biplot —", tr)
                     )
    ggsave(paste0("Results/plot_AMMI1_", tr, ".png"), p,
           width = 9, height = 7, dpi = 300)
  }, error = function(e) NULL)
}

# ── 11b. AMMI2 biplot (IPCA1 vs IPCA2) ───────────────────────────────────────
for (tr in traits) {
  if (is.null(ammi_results[[tr]])) next
  tryCatch({
    p <- plot_scores(ammi_results[[tr]],
                     type = 2,        # AMMI2: IPCA1 vs IPCA2
                     #title = paste("AMMI2 Biplot —", tr)
                     )
    ggsave(paste0("Results/plot_AMMI2_", tr, ".png"), p,
           width = 9, height = 7, dpi = 300)
  }, error = function(e) NULL)
}

# ── 11c. GGE which-won-where biplot ──────────────────────────────────────────
for (tr in traits) {
  if (is.null(gge_results[[tr]])) next
  tryCatch({
    p <- plot(gge_results[[tr]],
              type  = 3,              # which-won-where polygon
              #title = paste("GGE Which-Won-Where —", tr)
              )
    ggsave(paste0("Results/plot_GGE_whichwonwhere_", tr, ".png"), p,
           width = 9, height = 7, dpi = 300)
  }, error = function(e) NULL)
}

# ── 11d. GxE heatmap ordered by release year ─────────────────────────────────
release_yr <- release_yr |>
  rename(genotype = Genotype)
for (tr in traits) {
  tryCatch({
    mat <- blups_best_long %>%
      filter(trait == tr) %>%
      inner_join(release_yr, by = "genotype") %>%
      arrange(dev_year) %>%
      select(genotype, env, BLUP) %>%
      pivot_wider(names_from = env, values_from = BLUP) %>%
      column_to_rownames("genotype")

    mat_scaled <- scale(mat)    # scale within env

    mat_long <- as.data.frame(mat_scaled) %>%
      rownames_to_column("genotype") %>%
      pivot_longer(-genotype, names_to = "env", values_to = "z_BLUP") %>%
      left_join(release_yr, by = "genotype") %>%
      mutate(genotype = reorder(genotype, release_yr))

    p_hm <- ggplot(mat_long%>%
                     mutate(name = fct_reorder(genotype, dev_year)), #ordered by release year
                   aes(x = env, y = name, fill = z_BLUP)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#D7191C", mid = "white", high = "#2C7BB6",
                           midpoint = 0, name = "Standardised\nBLUP") +
      labs(title    = paste("GxE Heatmap —", tr),
           subtitle = "Genotypes ordered by release year (bottom = oldest)",
           x = "Environment", y = "Genotype") +
      theme_bw(base_size = 9) +
      theme(axis.text.y = element_text(size = 5),
            axis.text.x = element_text(angle = 30, hjust = 1))

    ggsave(paste0("Results/plot_GxE_heatmap_", tr, ".png"), p_hm,
           width = 8, height = 14, dpi = 300)
  }, error = function(e) NULL)
}

# ── 11e. PLS loadings: climate drivers of AMMI IPCA1 ────────────────────────
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY/Results")
pls_loadings_all <- bind_rows(
  lapply(list.files(pattern = "^PLS_loadings_.*\\.csv$"), read.csv)
)
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY")
if (nrow(pls_loadings_all) > 0) {
  p_pls <- pls_loadings_all %>%
    rename(IPCA1_loading = `Comp.1`) %>%
    group_by(trait) %>%
    slice_max(abs(IPCA1_loading), n = 10) %>%
    ungroup() %>%
    mutate(variable = reorder(variable, IPCA1_loading)) %>%
    ggplot(aes(x = IPCA1_loading, y = variable, fill = IPCA1_loading > 0)) +
    geom_col(show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_fill_manual(values = c("#D7191C","#2C7BB6")) +
    labs(title    = "PLS Loadings: Climate Drivers of AMMI IPCA1",
         subtitle = "Top 10 climate covariates per trait",
         x = "Loading (Comp. 1)", y = "Climate Covariate") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "#F0F0F0"))

  ggsave("Results/plot_PLS_climate_loadings.png", p_pls,
         width = 14, height = 8, dpi = 300)
}

# ── 11f. Stability × Release Year scatter (WAASB example) ───────────────────
for (tr in traits) {
  if (is.null(stab_vs_time[[tr]])) next
  tryCatch({
    df_plot <- stab_vs_time[[tr]]$data %>%
      filter(!is.na(WAASB), !is.na(dev_year)) %>%
      mutate(genotype = as.character(GEN))

    rho_val <- stab_vs_time[[tr]]$correlations %>%
      filter(stability_index == "WAASB") %>%
      pull(rho)

    p_sv <- ggplot(df_plot, aes(x = dev_year, y = WAASB)) +
      geom_point(aes(color = WAASB), size = 2.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE,
                  color = "#2C7BB6", fill = "#AED6F1") +
      scale_color_viridis_c(direction = -1, name = "WAASB") +
      geom_text_repel(aes(label = genotype), size = 2.2,
                      max.overlaps = 15) +
      annotate("text", x = -Inf, y = Inf,
               label = paste0("Spearman \u03C1 = ", rho_val),
               hjust = -0.1, vjust = 1.5, size = 3.5, color = "firebrick") +
      labs(title    = paste("Stability (WAASB) vs. Release Year —", tr),
           subtitle = "Lower WAASB = more stable across environments",
           x = "Release Year", y = "WAASB") +
      theme_bw(base_size = 11)

    ggsave(paste0("plot_stability_vs_release_yr_", tr, ".png"), p_sv,
           width = 10, height = 7, dpi = 300)
  }, error = function(e) NULL)
}

# ── 11g. Mean performance vs. stability (WAASB) quadrant plot ────────────────
for (tr in traits) {
  if (is.null(stab_all[[tr]])) next
  tryCatch({
    p_quad <- plot_scores(stab_all[[tr]]$waasb,
                   type  = 3,
                   #title = paste("Mean × Stability (WAASB) —", tr)
                   )
    ggsave(paste0("Results/plot_mean_vs_WAASB_", tr, ".png"), p_quad,
           width = 9, height = 7, dpi = 300)
  }, error = function(e) NULL)
}

# ── 11h. RF variable importance heatmap across traits ────────────────────────
rf_all <- bind_rows(rf_importance_list)
if (nrow(rf_all) > 0) {
  top_vars <- rf_all %>%
    group_by(trait) %>%
    slice_max(IncMSE, n = 8) %>%
    pull(variable) %>% unique()

  p_rf <- rf_all %>%
    filter(variable %in% top_vars) %>%
    ggplot(aes(x = trait, y = variable, fill = IncMSE)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(IncMSE, 1)), size = 3) +
    scale_fill_viridis_c(name = "%IncMSE") +
    labs(title    = "Random Forest: Climate Variable Importance for AMMI IPCA1",
         x = "Trait", y = "Climate Covariate") +
    theme_bw(base_size = 10)

  ggsave("Results/plot_RF_climate_importance.png", p_rf,
         width = 10, height = 8, dpi = 300)
}


# ── 12. MASTER SUMMARY TABLE ──────────────────────────────────────────────────
# One row per genotype per trait: mean BLUP + all stability indices + release yr

summary_master <- map_dfr(traits, function(tr) {
  if (is.null(stab_all[[tr]])) return(NULL)
  stab_all[[tr]]$stats$value %>%
    as.data.frame() %>%
    rename(genotype = GEN) %>%
    mutate(genotype = as.character(genotype), trait = tr) %>%
    inner_join(release_yr, by = "genotype")
})

write.csv(summary_master, "Results/master_summary_stability.csv", row.names = FALSE)

cat("\n\n✔ GxE + Enviromics pipeline complete.\n")
cat("Output files include:\n")
cat("  • AMMI_stability_<trait>.csv\n")
cat("  • GGE_scores_<trait>.csv\n")
cat("  • FW_stability_<trait>.csv\n")
cat("  • stability_indices_<trait>.csv\n")
cat("  • MGIDI_multi_trait_ranking.csv\n")
cat("  • climate_daily_NASA_POWER.csv\n")
cat("  • climate_covariates_by_phase.csv\n")
cat("  • climate_covariates_wide.csv\n")
cat("  • W_environmental_kinship.csv\n")
cat("  • PLS_loadings_<trait>.csv\n")
cat("  • RF_importance_<trait>.csv\n")
cat("  • stability_vs_release_year.csv\n")
cat("  • master_summary_stability.csv\n")
cat("  • 10+ publication-quality figures\n")
