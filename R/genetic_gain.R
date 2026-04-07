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

#genetic gain analysis####
# Genetic Gain Analysis in Lentils — ACTIVATE Rotation Trial
# Methodology: SpATS (Spatial) vs. GIBD (lme4) Comparison
# Author: [Your Name]
# Date: 2026
# Trial structure:
#   - 100 lines x 20 reps x 2 sites x 2 years
#   - Trial dimensions per site: 76 rows x 30 columns
#   - Each line has a known release year

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────

library(SpATS)       # Spatial analysis with P-splines
library(lme4)        # Mixed models for GIBD
library(lmerTest)    # p-values for lme4
library(dplyr)       # Data wrangling
library(tidyr)       # Data reshaping
library(purrr)       # Functional programming (map)
library(ggplot2)     # Visualization
library(ggpubr)      # Publication-ready ggplot themes
library(cowplot)     # Multi-panel plots
library(viridis)     # Color palettes
library(broom.mixed) # Tidy model outputs
library(corrplot)    # Correlation visualization


# ── 1. DATA LOADING & PREPARATION ─────────────────────────────────────────────

# Load your dataset — adjust path and format as needed
# Expected columns:
#   genotype   : lentil line identifier (factor)
#   dev_year : year of line release (numeric, e.g. 1995, 2003, ...)
#   site       : location (factor, 2 levels)
#   year       : trial year (factor, 2 levels)
#   rep        : replicate (factor)
#   block      : incomplete block nested within rep (factor)
#   row        : field row position (integer, 1–76)
#   col        : field column position (integer, 1–30)
#   yield      : seed yield kg/ha (numeric) — add other traits below
#   [other traits: days_flower, days_mature, seed_weight, plant_height, etc.]

dat <- read.csv("data/ACTIVATE_lentils_myYraw.csv", stringsAsFactors = FALSE)

#add a replication (gen_rep) col:
dat <- dat |>
  # 1. Group by Environment AND Genotype
  # This ensures the counter restarts for every new variety in every new site
  group_by(ENV, Lentil) |>
  # 2. Add the Rep column
  # row_number() simply assigns 1, 2, 3... to the rows in that group
  mutate(Rep_gen = row_number()) |>
  # 3. Ungroup (Important! removes the grouping structure so it doesn't mess up future stats)
  ungroup() |>
  # Optional: Move 'Rep' to be near the Genotype column for easier viewing
  relocate(Rep_gen, .after = Lentil) |>
  #rename Rep by Rep_combo
  rename(Rep_combo = Rep)

# ── Factor conversions --------------------------------------------------------
dat <- dat %>%
  mutate(
    genotype   = as.factor(Lentil),
    site       = as.factor(site),
    year       = as.factor(year),
    rep        = as.factor(Rep_gen),
    block      = as.factor(Block),
    row        = as.integer(Row),
    col        = as.integer(Col),
    env        = interaction(site, year, drop = TRUE)  # 4 environments
  )

# ── Define traits to analyze --------------------------------------------------
traits <- c("DTE",	"DTF",	"VegP",	"DTM",	"RepP",	"lodging",	"YLD",	"PRT",	"DS")
# Add or remove traits as available in your dataset

# ── Create environment list for looping --------------------------------------
envs <- levels(dat$env)
cat("Environments detected:", paste(envs, collapse = ", "), "\n")
cat("Traits to analyze:    ", paste(traits, collapse = ", "), "\n")
cat("Number of genotypes:  ", nlevels(dat$genotype), "\n")

# ── 2. HELPER FUNCTIONS ───────────────────────────────────────────────────────
# ── 2a. Cullis et al. (2006) entry-mean heritability ─────────────────────────
# H² = 1 − (mean pairwise PEV) / (2 * σ²g)
# Works for both SpATS and lme4 BLUPs

cullis_H2 <- function(model_blups, sigma2g) {
  # model_blups: named numeric vector of BLUPs for genotypes
  # sigma2g    : estimated genetic variance component
  if (is.null(model_blups) || sigma2g <= 0) return(NA_real_)
  n   <- length(model_blups)
  # Approximate mean PEV from BLUP shrinkage:
  # PEV_i ≈ σ²g − var(BLUP_i) for balanced data;
  # for lme4 we use the conditional variance approach below.
  # For SpATS, the package provides BLUP standard errors directly.
  mean_pev <- mean(attr(model_blups, "se")^2, na.rm = TRUE)
  H2 <- 1 - mean_pev / (2 * sigma2g)
  return(max(0, min(1, H2)))   # bound to [0,1]
}

# ── 2b. Extract Cullis H² from lme4 object ───────────────────────────────────
cullis_H2_lme4 <- function(lme4_model, genotype_var = "genotype") {
  vc       <- as.data.frame(VarCorr(lme4_model))
  sigma2g  <- vc$vcov[vc$grp == genotype_var]
  if (length(sigma2g) == 0 || sigma2g <= 0) return(list(H2 = NA, sigma2g = NA))
  
  # Conditional variances (PEV) for each genotype BLUP
  re_vals  <- ranef(lme4_model, condVar = TRUE)[[genotype_var]]
  pev      <- attr(re_vals, "postVar")[1, 1, ]   # diagonal of posterior variance
  mean_pev <- mean(pev, na.rm = TRUE)
  H2       <- max(0, min(1, 1 - mean_pev / (2 * sigma2g)))
  return(list(H2 = H2, sigma2g = sigma2g, mean_pev = mean_pev))
}

# ── 2c. Extract Cullis H² from SpATS object ──────────────────────────────────
cullis_H2_SpATS <- function(spats_model) {
  # SpATS stores genotype predictions with standard errors
  preds    <- predict(spats_model, which = "genotype")
  sigma2g  <- spats_model$var.comp["genotype"]
  if (is.na(sigma2g) || sigma2g <= 0) return(list(H2 = NA, sigma2g = NA))
  
  mean_pev <- mean(preds$standard.errors^2, na.rm = TRUE)
  H2       <- max(0, min(1, 1 - mean_pev / (2 * sigma2g)))
  return(list(H2 = H2, sigma2g = sigma2g, mean_pev = mean_pev))
}

# ── 3. MODEL FITTING: SpATS PER ENVIRONMENT ───────────────────────────────────
# SpATS fits: y = μ + Genotype + f(row, col) + Rep + ε
# Genotype is fitted as random to obtain BLUPs

fit_SpATS_env <- function(data_env, trait) {
  df <- data_env %>%
    filter(!is.na(.data[[trait]])) %>%
    droplevels()
  
  nrow_val <- max(df$row)
  ncol_val <- max(df$col)
  
  # Number of knots: rule of thumb = min(round(n/4), 35)
  nseg_row <- min(round(nrow_val / 4), 35)
  nseg_col <- min(round(ncol_val / 4), 35)
  
  tryCatch({
    mod <- SpATS(
      response    = trait,
      spatial     = ~ PSANOVA(col, row,
                              nseg     = c(nseg_col, nseg_row),
                              degree   = c(3, 3),
                              pord     = c(2, 2)),
      genotype    = "genotype",
      genotype.as.random = TRUE,            # BLUPs for heritability
      fixed       = ~ rep,                  # Rep as fixed
      random      = ~ block,               # Block within rep as random (optional)
      data        = df,
      control     = list(maxit = 100, tolerance = 1e-06, monitoring = 0)
    )
    return(mod)
  }, error = function(e) {
    message("SpATS failed for trait ", trait, ": ", e$message)
    return(NULL)
  })
}

# ── 4. MODEL FITTING: GIBD (lme4) PER ENVIRONMENT ─────────────────────────────
# y_ijk = μ + Genotype_i + Rep_j + Block_k(Rep_j) + ε_ijk
# Genotype random → BLUPs; Rep fixed; Block(Rep) random

fit_GIBD_env <- function(data_env, trait) {
  df <- data_env %>%
    filter(!is.na(.data[[trait]])) %>%
    droplevels()
  
  formula_str <- paste0(trait,
                        " ~ rep + (1|genotype) + (1|rep:block)")
  f <- as.formula(formula_str)
  
  tryCatch({
    mod <- lmer(f, data = df,
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl   = list(maxfun = 2e5)))
    return(mod)
  }, error = function(e) {
    message("GIBD failed for trait ", trait, ": ", e$message)
    return(NULL)
  })
}


# ── 5. MAIN LOOP: FIT BOTH MODELS ACROSS ALL ENVIRONMENTS & TRAITS ────────────

results_list <- list()   # Store all results

for (tr in traits) {
  cat("\n══════════════════════════════════════════\n")
  cat(" Trait:", tr, "\n")
  cat("══════════════════════════════════════════\n")
  
  blups_SpATS_list <- list()
  blups_GIBD_list  <- list()
  H2_table_rows    <- list()
  
  for (env in envs) {
    cat("  → Environment:", env, "\n")
    d_env <- dat %>% filter(env == !!env) %>% droplevels()
    
    # ── SpATS ----------------------------------------------------------------
    mod_sp <- fit_SpATS_env(d_env, tr)
    if (!is.null(mod_sp)) {
      h2_sp  <- cullis_H2_SpATS(mod_sp)
      preds  <- predict(mod_sp, which = "genotype")
      pev_sp <- predict(mod_sp, which = "genotype")$standard.errors^2 # Extract PEVs from Stage 1 SpATS model
      blups_sp <- setNames(preds$predicted.values, preds$genotype)
      blups_SpATS_list[[env]] <- data.frame(
        genotype = names(blups_sp),
        BLUP     = as.numeric(blups_sp),
        env      = env,
        trait    = tr,
        pev      = pev_sp,
        model    = "SpATS"
      )
      ed_spatial <- summary(mod_sp)$p.table.random    # effective dimensions
    } else {
      h2_sp <- list(H2 = NA, sigma2g = NA)
      ed_spatial <- NA
    }
    
    # ── GIBD ----------------------------------------------------------------
    mod_gi <- fit_GIBD_env(d_env, tr)
    if (!is.null(mod_gi)) {
      h2_gi  <- cullis_H2_lme4(mod_gi, "genotype")
      re_gi  <- ranef(mod_gi)$genotype
      
      # Extract PEVs from Stage 1 lme4 model
      re_lme4  <- ranef(mod_gi, condVar = TRUE)$genotype
      pev_lme4 <- attr(re_lme4, "postVar")[1, 1, ]
      
      #Calculate the grand mean of the trait for this specific environment
      grand_mean <- mean(d_env[[tr]], na.rm = TRUE)
      
      #Add the grand mean to the random deviations to get predicted means
      blups_gi <- setNames(re_gi[, 1] + grand_mean, rownames(re_gi))
      
      blups_GIBD_list[[env]] <- data.frame(
        genotype = names(blups_gi),
        BLUP     = as.numeric(blups_gi),
        env      = env,
        trait    = tr,
        pev      = pev_lme4,
        model    = "GIBD"
      )
    } else {
      h2_gi <- list(H2 = NA, sigma2g = NA)
    }
    
    # ── Collect H² row -------------------------------------------------------
    H2_table_rows[[env]] <- data.frame(
      trait         = tr,
      env           = env,
      H2_SpATS      = round(h2_sp$H2, 3),
      sigma2g_SpATS = round(h2_sp$sigma2g, 4),
      H2_GIBD       = round(h2_gi$H2, 3),
      sigma2g_GIBD  = round(h2_gi$sigma2g, 4),
      delta_H2      = round(h2_sp$H2 - h2_gi$H2, 3)
    )
    
    # ── BLUP correlation between methods (Spearman) --------------------------
    if (!is.null(mod_sp) && !is.null(mod_gi)) {
      common_g <- intersect(names(blups_sp), names(blups_gi))
      rho <- cor(blups_sp[common_g], blups_gi[common_g],
                 method = "spearman", use = "complete.obs")
      H2_table_rows[[env]]$BLUP_rho_Spearman <- round(rho, 3)
      cat("    SpATS H²=", h2_sp$H2, "| GIBD H²=", h2_gi$H2,
          "| ΔH²=", h2_sp$H2 - h2_gi$H2,
          "| BLUP ρ=", round(rho, 3), "\n")
    }
  }
  
  results_list[[tr]] <- list(
    H2_table     = bind_rows(H2_table_rows),
    blups_SpATS  = bind_rows(blups_SpATS_list),
    blups_GIBD   = bind_rows(blups_GIBD_list)
  )
}
blups_SpATS_all <- bind_rows(lapply(results_list, function(x) x$blups_SpATS))
blups_GIBD_all <- bind_rows(lapply(results_list, function(x) x$blups_GIBD))

write.csv(blups_SpATS_all, "Results/blups_SpATS_all.csv", row.names = FALSE)
write.csv(blups_GIBD_all, "Results/blups_GIBD_all.csv", row.names = FALSE)

# ── Compile master H² summary table ------------------------------------------
H2_master <- bind_rows(lapply(results_list, function(x) x$H2_table))
cat("\n\n══ H² SUMMARY TABLE ══\n")
print(H2_master)
write.csv(H2_master, "Results/H2_comparison_SpATS_vs_GIBD.csv", row.names = FALSE)
H2_master <- read.csv("Results/H2_comparison_SpATS_vs_GIBD.csv")

# ── 6. SELECT BEST MODEL PER TRAIT × ENVIRONMENT ─────────────────────────────
# Criterion: higher Cullis H² → better model
# Tie-break: if |ΔH²| < 0.02 → prefer SpATS (explicit spatial correction)

best_model_selection <- H2_master %>%
  mutate(
    best_model = case_when(
      is.na(H2_SpATS) & is.na(H2_GIBD) ~ "none",
      is.na(H2_SpATS)  ~ "GIBD",
      is.na(H2_GIBD)   ~ "SpATS",
      delta_H2 >=  0.05 ~ "SpATS",
      delta_H2 <= -0.05 ~ "GIBD",
      TRUE              ~ "SpATS"   # tie → prefer SpATS
    )
  )

cat("\n══ BEST MODEL SELECTION ══\n")
print(best_model_selection %>% select(trait, env, H2_SpATS, H2_GIBD, delta_H2, best_model))
write.csv(best_model_selection, "Results/best_model_selection.csv", row.names = FALSE)
best_model_selection<- read.csv("Results/best_model_selection.csv")

# ── 7. STAGE 2: MULTI-ENVIRONMENT BLUPs (Best-model BLUPs) ───────────────────
# Combine best-model BLUPs across environments using a second-stage mixed model
# y_env = μ + Genotype_i + Environment_j + G×E_ij + ε

fit_stage2 <- function(blup_df, trait_name) {
  # blup_df: long-format data frame with columns: genotype, BLUP, env
  tryCatch({
    mod2 <- lmer(BLUP ~ (1 | genotype) + env,
                 data    = blup_df,
                 weights = 1 / pev, #inverse prediction error variance (PEV)
                 control = lmerControl(optimizer = "bobyqa"))
    
    # Extract the random deviations for the genotypes
    re   <- ranef(mod2)$genotype
    
    # [NEW] Calculate the grand mean of the trait across all environments
    grand_mean <- mean(blup_df$BLUP, na.rm = TRUE)
    
    # [NEW] Add the grand mean to the deviations to get the predicted mean BLUPs
    predicted_means <- re[, 1] + grand_mean
    
    blups_combined <- data.frame(
      genotype      = rownames(re),
      BLUP_combined = predicted_means,
      trait         = trait_name
    )
    return(blups_combined)
  }, error = function(e) {
    message("Stage 2 failed for ", trait_name, ": ", e$message)
    return(NULL)
  })
}

# Build best-model BLUP dataset for each trait
combined_blups_list <- list()

for (tr in traits) {
  blup_sp  <- results_list[[tr]]$blups_SpATS
  blup_gi  <- results_list[[tr]]$blups_GIBD
  best_sel <- best_model_selection %>% filter(trait == tr)
  
  blup_best <- bind_rows(
    blup_sp %>% inner_join(
      best_sel %>% filter(best_model == "SpATS") %>% select(env),
      by = "env"),
    blup_gi %>% inner_join(
      best_sel %>% filter(best_model == "GIBD") %>% select(env),
      by = "env")
  )
  
  if (nrow(blup_best) > 0) {
    combined_blups_list[[tr]] <- fit_stage2(blup_best, tr)
  }
}

combined_blups <- bind_rows(combined_blups_list)

# Attach release year to combined BLUPs
dev_year_df <- myY %>%
  select(Genotype, dev_year) %>%
  mutate(genotype = Genotype) %>%
  distinct()

combined_blups <- combined_blups %>%
  inner_join(dev_year_df, by = "genotype")

write.csv(combined_blups, "Results/combined_BLUPs_best_model.csv", row.names = FALSE)
combined_blups <- read.csv("Results/combined_BLUPs_best_model.csv")

# ── 8. GENETIC GAIN REGRESSION ────────────────────────────────────────────────
# Regress combined BLUPs on release year per trait
# Model: BLUP_combined ~ dev_year

genetic_gain_results <- combined_blups %>%
  group_by(trait) %>%
  summarise(
    n_lines    = n(),
    yr_range   = paste(min(dev_year, na.rm = TRUE),
                       max(dev_year, na.rm = TRUE), sep = "–"),
    slope      = coef(lm(BLUP_combined ~ dev_year))[2],
    slope_se   = summary(lm(BLUP_combined ~ dev_year))$coefficients[2, 2],
    intercept  = coef(lm(BLUP_combined ~ dev_year))[1],
    R2         = summary(lm(BLUP_combined ~ dev_year))$r.squared,
    p_value    = summary(lm(BLUP_combined ~ dev_year))$coefficients[2, 4],
    grand_mean = mean(BLUP_combined, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(
    gain_pct_yr = slope / grand_mean * 100,   # % gain per year relative to grand mean
    sig_label   = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

cat("\n══ GENETIC GAIN ESTIMATES ══\n")
print(genetic_gain_results)
write.csv(genetic_gain_results, "Results/genetic_gain_estimates.csv", row.names = FALSE)
genetic_gain_results <- read.csv("Results/genetic_gain_estimates.csv")

# ── 9. VISUALIZATION ──────────────────────────────────────────────────────────

# ── 9a. H² Comparison: SpATS vs GIBD (faceted by trait) ─────────────────────
p_H2 <- H2_master %>%
  pivot_longer(cols      = c(H2_SpATS, H2_GIBD),
               names_to  = "model",
               names_prefix = "H2_",
               values_to = "H2") %>%
  ggplot(aes(x = env, y = H2, fill = model)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "white") +
  geom_text(aes(label = round(H2, 2)),
            position = position_dodge(0.7), vjust = -0.4, size = 3) +
  facet_wrap(~ trait, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("SpATS" = "#2C7BB6", "GIBD" = "#D7191C"),
                    labels = c("SpATS", "GIBD (lme4)")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Cullis Entry-Mean Heritability: SpATS vs. GIBD",
    subtitle = "ACTIVATE lentil rotation trial — per trait × environment",
    x        = "Environment (Site × Year)",
    y        = expression(H^2 ~ "(Cullis)"),
    fill     = "Model"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#F0F0F0"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "top"
  )

ggsave("Results/plot_H2_comparison.png", p_H2, width = 12, height = 7, dpi = 300)

# ── 9b. ΔH² heatmap (SpATS − GIBD) ──────────────────────────────────────────
p_deltaH2 <- H2_master %>%
  ggplot(aes(x = env, y = trait, fill = delta_H2)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", delta_H2)), size = 3.5) +
  scale_fill_gradient2(
    low      = "#D7191C",
    mid      = "white",
    high     = "#2C7BB6",
    midpoint = 0,
    limits   = c(-0.2, 0.2),
    oob      = scales::squish,
    name     = "ΔH²\n(SpATS−GIBD)"
  ) +
  labs(
    title    = "ΔH² Heatmap (SpATS − GIBD)",
    subtitle = "Blue = SpATS superior | Red = GIBD superior",
    x        = "Environment", y = "Trait"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("Results/plot_deltaH2_heatmap.png", p_deltaH2, width = 9, height = 5, dpi = 300)

# ── 9c. BLUP Spearman correlation scatter (SpATS vs GIBD, per env × trait) ───
# Merge both BLUP sets for all traits/envs
all_blups_wide <- bind_rows(
  bind_rows(lapply(results_list, function(x) x$blups_SpATS)),
  bind_rows(lapply(results_list, function(x) x$blups_GIBD))
) %>%
  select(!pev) |>
  pivot_wider(names_from = model, values_from = BLUP)

p_rho <- all_blups_wide %>%
  ggplot(aes(x = SpATS, y = GIBD)) +
  geom_point(alpha = 0.5, size = 1.2, color = "#444444") +
  geom_smooth(method = "lm", se = TRUE, color = "#2C7BB6", linewidth = 0.8) +
  facet_wrap(trait ~ env, scales = "free") +
  stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.92,
           size = 3, color = "firebrick") +
  labs(
    title    = "BLUP Concordance: SpATS vs. GIBD",
    subtitle = "Pearson ρ per trait × environment",
    x        = "BLUP (SpATS)", y = "BLUP (GIBD)"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "#F0F0F0"))

ggsave("Results/plot_BLUP_concordance.png", p_rho, width = 14, height = 10, dpi = 300)

# ── 9d. Genetic Gain regression plots per trait ───────────────────────────────
p_gain <- combined_blups %>%
  ggplot(aes(x = dev_year, y = BLUP_combined)) +
  geom_point(alpha = 0.7, color = "#555555", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#2C7BB6", fill = "#AED6F1") +
  facet_wrap(~ trait, scales = "free_y", nrow = 2) +
  geom_text(
    data  = genetic_gain_results,
    aes(
      x     = -Inf, y = Inf,
      label = paste0("b = ", round(slope, 3), " ", sig_label,
                     "\nR² = ", round(R2, 2))
    ),
    hjust = -0.1, vjust = 1.3, size = 3.2, inherit.aes = FALSE
  ) +
  labs(
    title    = "Genetic Gain in Lentil — ACTIVATE Rotation Trial",
    subtitle = "Best-model combined BLUPs regressed on release year",
    x        = "Release Year",
    y        = "Combined BLUP"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F0F0F0"))

ggsave("Results/plot_genetic_gain.png", p_gain, width = 12, height = 7, dpi = 300)

# ── 9e. SpATS spatial trend map for each environment (first trait) ────────────
# Run separately for each environment — example for first trait & first env
env1    <- envs[2]
trait1  <- traits[7]
d_env1  <- dat %>% filter(env == env1) %>% droplevels()
mod_sp1 <- fit_SpATS_env(d_env1, trait1)

if (!is.null(mod_sp1)) {
  png("Results/plot_SpATS_spatial_trend_env1.png", width = 1800, height = 1200, res = 200)
  plot(mod_sp1, main = paste("SpATS Spatial Trend —", env1, "|", trait1))
  dev.off()
}


# ── 10. SUMMARY OUTPUT TABLE FOR PAPER ────────────────────────────────────────
paper_table <- genetic_gain_results %>%
  left_join(
    best_model_selection %>%
      group_by(trait) %>%
      summarise(
        H2_SpATS_mean          = round(mean(H2_SpATS, na.rm = TRUE), 3),
        H2_GIBD_mean           = round(mean(H2_GIBD,  na.rm = TRUE), 3),
        best_model_most_common = names(sort(table(best_model), decreasing = TRUE))[1]
      ),
    by = "trait"
  ) %>%
  select(
    trait, n_lines, yr_range,
    H2_SpATS_mean, H2_GIBD_mean, best_model_most_common,
    slope, slope_se, gain_pct_yr, R2, p_value, sig_label
  )

cat("\n══ PAPER TABLE: GENETIC GAIN SUMMARY ══\n")
print(paper_table)
write.csv(paper_table, "Results/paper_table_genetic_gain.csv", row.names = FALSE)

cat("\n\n✔ Analysis complete. Output files:\n")
cat("  • H2_comparison_SpATS_vs_GIBD.csv\n")
cat("  • best_model_selection.csv\n")
cat("  • combined_BLUPs_best_model.csv\n")
cat("  • genetic_gain_estimates.csv\n")
cat("  • paper_table_genetic_gain.csv\n")
cat("  • plot_H2_comparison.png\n")
cat("  • plot_deltaH2_heatmap.png\n")
cat("  • plot_BLUP_concordance.png\n")
cat("  • plot_genetic_gain.png\n")
cat("  • plot_SpATS_spatial_trend_env1.png\n")
