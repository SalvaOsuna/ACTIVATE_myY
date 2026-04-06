# =============================================================================
# Results Figures & Tables — Lentil Genetic Gain Paper
# Generates all outputs in section order (3.1 → 3.6)
# Assumes outputs from genetic_gain_analysis.R and gxe_enviromics_analysis.R
# are available in the working directory.
# =============================================================================

library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
library(viridis); library(scales); library(patchwork); library(tibble)
library(forcats); library(stringr)

dir.create("Results", showWarnings = FALSE)

# Traits and environments (adjust to match your actual data)
traits <- c("DTE",	"DTF",	"VegP",	"DTM",	"RepP",	"lodging",	"YLD",	"PRT",	"DS")
envs   <- c("Clavet.2024","Clavet.2025","Hunter.2024","Hunter.2025")
env_colors <- c("Clavet.2024"="#2C7BB6","Clavet.2025"="#74ADD1",
                "Hunter.2024"="#D7191C","Hunter.2025"="#F77964")

# =============================================================================
# SECTION 3.1 — ENVIRONMENTAL CHARACTERISATION & CLIMATE VARIATION
# =============================================================================

# ── Fig 3.1a: Climate covariate summary heatmap — actual growing season -------
# Covariates are computed for the sowing-to-harvest window per environment,
# using field-observed phenological dates rather than the calendar full_season.

# Actual phenological windows from field observations
pheno_windows <- data.frame(
  site       = c("Clavet", "Clavet", "Hunter", "Hunter"),
  trial_year = c(2024,     2025,     2024,     2025),
  sowing     = as.Date(c("2024-05-15","2025-05-05",
                         "2024-05-28","2025-05-12")),
  harvest    = as.Date(c("2024-08-26","2025-08-29",
                         "2024-09-09","2025-09-13")),
  Avg_DTE    = c(12.75, 14.77, 16.01, 11.11),
  Avg_DTF    = c(47.48, 53.37, 53.94, 55.34),
  Avg_DTM    = c(81.56, 87.07, 93.77, 95.00)
) %>%
  mutate(
    emergence = as.Date(sowing + Avg_DTE),
    flowering = as.Date(sowing + Avg_DTF),
    maturity  = as.Date(sowing + Avg_DTM),
    env       = paste(site, trial_year, sep = ".")
  ) %>%
  select(-Avg_DTE, -Avg_DTF, -Avg_DTM)

# Load daily climate data (output from gxe_enviromics_analysis.R)
climate_daily <- read.csv("Results/climate_daily_NASA_POWER.csv") %>%
  mutate(
    date = as.Date(as.character(YYYYMMDD)),
    env  = paste(site, trial_year, sep = ".")
  )

# Helper: compute covariates for a date window within one environment
compute_gs_covariates <- function(clim_env, start_date, end_date,
                                  T_base = 5) {
  d <- clim_env %>%
    filter(date >= start_date & date <= end_date)
  if (nrow(d) == 0) return(NULL)
  
  GDD          <- sum(pmax((d$T2M_MAX + d$T2M_MIN) / 2 - T_base, 0),
                      na.rm = TRUE)
  precip_total <- sum(d$PRECTOTCORR, na.rm = TRUE)
  mean_rad     <- mean(d$ALLSKY_SFC_SW_DWN, na.rm = TRUE)
  
  es_max <- 0.6108 * exp(17.27 * d$T2M_MAX / (d$T2M_MAX + 237.3))
  es_min <- 0.6108 * exp(17.27 * d$T2M_MIN / (d$T2M_MIN + 237.3))
  ea     <- ((es_max + es_min) / 2) * (d$RH2M / 100)
  VPD    <- mean((es_max + es_min) / 2 - ea, na.rm = TRUE)
  
  heat_days  <- sum(d$T2M_MAX > 30,  na.rm = TRUE)
  frost_days <- sum(d$T2M_MIN < 0,   na.rm = TRUE)
  dtr        <- mean(d$T2M_MAX - d$T2M_MIN, na.rm = TRUE)
  tmean      <- mean(d$T2M,     na.rm = TRUE)
  
  data.frame(GDD, precip_total, mean_rad, VPD,
             heat_days, frost_days, dtr, tmean)
}

# Compute growing-season covariates for each environment
gs_covariates <- purrr::map_dfr(seq_len(nrow(pheno_windows)), function(i) {
  pw  <- pheno_windows[i, ]
  cdf <- climate_daily %>% filter(env == pw$env)
  cov <- compute_gs_covariates(cdf, pw$sowing, pw$harvest)
  if (is.null(cov)) return(NULL)
  cov$env        <- pw$env
  cov$site       <- pw$site
  cov$trial_year <- pw$trial_year
  cov$gs_days    <- as.integer(pw$harvest - pw$sowing)   # season length
  cov
})

write.csv(gs_covariates,
          "Results/climate_covariates_growing_season.csv",
          row.names = FALSE)

# Pivot to long and scale for heatmap
covariate_labels <- c(
  "GDD"          = "GDD (\u00B0C\u00B7d)",
  "precip_total" = "Precipitation (mm)",
  "mean_rad"     = "Solar Rad. (MJ m\u207B\u00B2)",
  "VPD"          = "VPD (kPa)",
  "heat_days"    = "Heat stress days (>30\u00B0C)",
  "frost_days"   = "Frost days (<0\u00B0C)",
  "dtr"          = "Diurnal T range (\u00B0C)",
  "tmean"        = "Mean temperature (\u00B0C)",
  "gs_days"      = "Growing season length (d)"
)

cov_season <- gs_covariates %>%
  pivot_longer(cols      = names(covariate_labels),
               names_to  = "covariate",
               values_to = "value") %>%
  group_by(covariate) %>%
  mutate(z = as.numeric(scale(value))) %>%
  ungroup() %>%
  mutate(
    covariate = recode(covariate, !!!covariate_labels),
    covariate = factor(covariate, levels = covariate_labels),
    env       = factor(env, levels = envs),
    # Format label: integers for counts/days, 1 decimal for continuous
    label = ifelse(
      covariate %in% c("Heat stress days (>30\u00B0C)",
                       "Frost days (<0\u00B0C)",
                       "Growing season length (d)"),
      as.character(round(value, 0)),
      as.character(round(value, 1))
    )
  )

p_clim_heat <- ggplot(cov_season,
                      aes(x = env, y = fct_rev(covariate), fill = z)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = label,
                color  = abs(z) > 0.8),   # white text on saturated cells
            size = 3.2, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient2(
    low      = "#D7191C",
    mid      = "white",
    high     = "#2C7BB6",
    midpoint = 0,
    name     = "Z-score"
  ) +
  scale_color_manual(values = c("FALSE" = "grey20", "TRUE" = "grey20")) +
  scale_x_discrete(position = "top") +
  labs(
    #title    = "Growing-season climate covariates by environment",
    #subtitle = paste0("Sowing to harvest window per site \u00D7 year | ",
    #                  "Cell values = raw; colour = Z-score across environments"),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 0, face = "bold"),
    axis.text.y  = element_text(face = "bold"),
    panel.grid   = element_blank(),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("Results/Fig3.1a_climate_covariate_heatmap.png",
       p_clim_heat, width = 9, height = 6.5, dpi = 300)
cat("Fig 3.1a saved — growing-season climate heatmap\n")

# ── Fig 3.1b: W matrix (upper triangle, ggplot2) ------------------------------
W_matrix <- as.matrix(read.csv("Results/W_environmental_kinship.csv",
                               row.names = 1))
W_long <- as.data.frame(W_matrix) %>%
  rownames_to_column("env1") %>%
  pivot_longer(-env1, names_to = "env2", values_to = "w") %>%
  mutate(env1 = factor(env1, levels = envs),
         env2 = factor(env2, levels = envs),
         ri   = as.integer(env1),
         ci   = as.integer(env2)) %>%
  filter(ci >= ri)

p_W_2 <- ggplot(W_long, aes(x = env2, y = fct_rev(env1), fill = w)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(w, 3)), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "#D7191C", mid = "white", high = "#2C7BB6",
                       midpoint = 0, limits = c(-1, 1),
                       oob = squish, name = "Kinship (W)") +
  scale_x_discrete(position = "top") +
  labs(#title    = "Environmental kinship matrix (W)",
       #subtitle = "Derived from 28 standardised climate covariates",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 0, face = "bold"),
        axis.text.y  = element_text(face = "bold"),
        panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold"))
ggsave("Results/Fig3.1b_W_matrix.png", p_W,
       width = 6, height = 5, dpi = 300)

# ── Fig 3.1c: Environment PCA biplot -----------------------------------------
cov_wide   <- read.csv("Results/climate_covariates_wide.csv", row.names = 1)
cov_scaled <- scale(cov_wide)
cov_scaled[is.nan(cov_scaled)] <- 0
col_variances <- apply(cov_scaled, 2, var, na.rm = TRUE) # Calculate the variance of every column
names(col_variances[col_variances == 0 | is.na(col_variances)]) # Print the names of the columns where variance is exactly 0 (or NA)
cov_filtered <- as.data.frame(cov_scaled) %>% 
  dplyr::select(where(~ var(., na.rm = TRUE) > 0))# Keep only columns where the variance is greater than 0
env_pca    <- prcomp(cov_filtered, scale. = FALSE)# Now run your PCA on the filtered dataset
pct_var    <- round(100 * env_pca$sdev^2 / sum(env_pca$sdev^2), 1)

env_scores <- as.data.frame(env_pca$x[, 1:2]) %>%
  rownames_to_column("env")
scale_f <- max(abs(env_scores[, 2:3])) /
  max(abs(env_pca$rotation[, 1:2])) * 0.65
var_loads <- as.data.frame(env_pca$rotation[, 1:2]) %>%
  rownames_to_column("variable") %>%
  mutate(PC1s = PC1 * scale_f, PC2s = PC2 * scale_f,
         mag  = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(mag, n = 15) %>%
  mutate(label = str_replace_all(variable, "_", " "))

p_pca <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_segment(data = var_loads,
               aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
               arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
               color = "grey40", linewidth = 0.6) +
  geom_text_repel(data = var_loads,
                  aes(x = PC1s, y = PC2s, label = label),
                  color = "grey30", size = 2.8, box.padding = 0.3) +
  geom_point(data = env_scores,
             aes(x = PC1, y = PC2, color = env), size = 6) +
  geom_text_repel(data = env_scores,
                  aes(x = PC1, y = PC2, label = env, color = env),
                  fontface = "bold", size = 3.2, box.padding = 0.5) +
  scale_color_manual(values = env_colors, guide = "none") +
  labs(#title = "Environment PCA — climate covariates",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))
ggsave("Results/Fig3.1c_environment_PCA.png", p_pca,
       width = 8, height = 6, dpi = 300)

cat("Section 3.1 figures saved\n")

#combine figures together
Fig1 <- (p_clim_heat | (p_W / p_W_2)) / p_pca+ 
  plot_layout(tag_level = 'new', heights = unit(c(9, 7), c('cm', 'cm'))) +
  plot_annotation(tag_levels = list(c('a)', 'b)', 'c)', 'd)')))

ggsave("Figures/Fig1_ENVIRONMENTAL CHARACTERISATION & CLIMATE VARIATION.png", Fig1,
       width = 8.5, height = 8, dpi = 600)

# =============================================================================
# SECTION 3.2 — SPATIAL VARIATION & MODEL COMPARISON
# =============================================================================

# ── Table 3.2: H² summary per model × trait × environment -------------------
H2_master <- read.csv("Results/H2_comparison_SpATS_vs_GIBD.csv")
best_sel   <- read.csv("Results/best_model_selection.csv")

# Wide format for the paper table
H2_table <- H2_master %>%
  select(trait, env, H2_SpATS, H2_GIBD, delta_H2, BLUP_rho_Spearman) %>%
  arrange(trait, env)
write.csv(H2_table, "Results/Table3.2_H2_model_comparison.csv",
          row.names = FALSE)

# ── Fig 3.2a: ΔH² heatmap (SpATS − GIBD) ------------------------------------
p_dH2 <- H2_master %>%
  mutate(env   = factor(env, levels = envs),
         trait = factor(trait, levels = traits)) %>%
  ggplot(aes(x = env, y = fct_rev(trait), fill = delta_H2)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.3f", delta_H2)), size = 3.5,
            fontface = "bold") +
  scale_fill_gradient2(low = "#D7191C", mid = "white", high = "#2C7BB6",
                       midpoint = 0, limits = c(-0.25, 0.25),
                       oob = squish,
                       name = "\u0394H\u00B2\n(SpATS\u2212GIBD)") +
  labs(title    = "\u0394H\u00B2 — SpATS vs. GIBD",
       subtitle = "Blue = SpATS superior | Red = GIBD superior",
       x = "Environment", y = "Trait") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold"))
ggsave("Results/Fig3.2a_deltaH2_heatmap.png", p_dH2,
       width = 8, height = 5.5, dpi = 300)

# ── Fig 3.2b: H² grouped bar chart (SpATS vs GIBD, faceted by trait) ---------
p_H2bar <- H2_master %>%
  pivot_longer(cols = c(H2_SpATS, H2_GIBD),
               names_to = "model", names_prefix = "H2_",
               values_to = "H2") %>%
  mutate(env   = factor(env, levels = envs),
         trait = factor(trait, levels = traits),
         model = factor(model, levels = c("SpATS","GIBD"))) %>%
  ggplot(aes(x = env, y = H2, fill = model)) +
  geom_col(position = position_dodge(0.7), width = 0.65, color = "white") +
  geom_text(aes(label = round(H2, 2)),
            position = position_dodge(0.7),
            vjust = -0.3, size = 2.6) +
  facet_wrap(~ trait, nrow = 3) +
  scale_fill_manual(values = c("SpATS" = "#2C7BB6", "GIBD" = "#D7191C")) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(title = "Cullis entry-mean heritability: SpATS vs. GIBD",
       x = "Environment", y = expression(H^2~"(Cullis)"),
       fill = "Model") +
  theme_bw(base_size = 10) +
  theme(axis.text.x      = element_text(angle = 35, hjust = 1, size = 7),
        strip.background = element_rect(fill = "#F0F0F0"),
        legend.position  = "top",
        plot.title       = element_text(face = "bold"))
ggsave("Results/Fig3.2b_H2_barplot.png", p_H2bar,
       width = 12, height = 9, dpi = 300)

# ── Fig 3.2c: BLUP concordance scatter (SpATS vs GIBD) -----------------------
blups_SpATS <- read.csv("Results/blups_SpATS_all.csv")
blups_GIBD  <- read.csv("Results/blups_GIBD_all.csv")

blups_wide <- blups_SpATS %>%
  rename(BLUP_SpATS = BLUP) %>%
  inner_join(blups_GIBD %>% rename(BLUP_GIBD = BLUP),
             by = c("genotype","env","trait")) %>%
  mutate(env   = factor(env, levels = envs),
         trait = factor(trait, levels = traits))

rho_labels <- blups_wide %>%
  group_by(trait, env) %>%
  summarise(rho = round(cor(BLUP_SpATS, BLUP_GIBD,
                            method = "spearman",
                            use = "complete.obs"), 3),
            .groups = "drop") %>%
  mutate(label = paste0("\u03C1 = ", rho))

p_rho <- ggplot(blups_wide, aes(x = BLUP_SpATS, y = BLUP_GIBD)) +
  geom_point(alpha = 0.35, size = 1, color = "#444444") +
  geom_smooth(method = "lm", se = FALSE,
              color = "#2C7BB6", linewidth = 0.7) +
  geom_text(data = rho_labels,
            aes(label = label, x = -Inf, y = Inf),
            hjust = -0.1, vjust = 1.4, size = 2.6,
            color = "firebrick", inherit.aes = FALSE) +
  ggh4x::facet_grid2(trait ~ env, scales = "free", independent = "all") +
  labs(title = "BLUP concordance: SpATS vs. GIBD",
       x = "BLUP (SpATS)", y = "BLUP (GIBD)") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        plot.title       = element_text(face = "bold"))
ggsave("Results/Fig3.2c_BLUP_concordance.png", p_rho,
       width = 14, height = 18, dpi = 300)

cat("Section 3.2 figures saved\n")


# =============================================================================
# SECTION 3.3 — DESCRIPTIVE STATISTICS & TRAIT VARIATION
# =============================================================================

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

# ── Table 3.3: Descriptive statistics per trait × environment ----------------
library(moments)
desc_stats <- dat %>%
  pivot_longer(cols = all_of(traits), names_to = "trait",
               values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(trait, env) %>%
  summarise(
    n    = n(),
    mean = round(mean(value), 2),
    sd   = round(sd(value), 2),
    cv   = round(sd(value) / mean(value) * 100, 1),
    min  = round(min(value), 2),
    max  = round(max(value), 2),
    skewness = round(skewness(value),2),
    kurtosis= round(kurtosis(value), 2),
    .groups = "drop"
  )
write.csv(desc_stats, "Results/Table3.3_descriptive_stats.csv",
          row.names = FALSE)

# ── Fig 3.3a: Violin + boxplot per trait (all envs overlaid) -----------------
p_violin <- dat %>%
  pivot_longer(cols = all_of(traits), names_to = "trait",
               values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(env   = factor(as.character(env), levels = envs),
         trait = factor(trait, levels = traits)) %>%
  ggplot(aes(x = env, y = value, fill = env)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.6,
               color = "grey30", fill = "white", alpha = 0.8) +
  facet_wrap(~ trait, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = env_colors, name = "Environment") +
  labs(title = "Phenotypic distribution per trait and environment",
       x = NULL, y = "Phenotypic value") +
  theme_bw(base_size = 12) +
  theme(axis.text.x      = element_text(angle = 35, hjust = 1, size = 9),
        strip.background = element_rect(fill = "#F0F0F0"),
        legend.position  = "none",
        plot.title       = element_text(face = "bold"))
ggsave("Results/Fig3.3a_trait_distributions.png", p_violin,
       width = 13, height = 10, dpi = 600)

# ── Fig 3.3b: Phenotypic correlation heatmap (grand means across envs) -------
trait_means <- dat %>%
  group_by(genotype) %>%
  summarise(across(all_of(traits), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  select(all_of(traits))

cor_mat <- cor(trait_means, use = "pairwise.complete.obs",
               method = "pearson")
cor_long <- as.data.frame(cor_mat) %>%
  rownames_to_column("t1") %>%
  pivot_longer(-t1, names_to = "t2", values_to = "r") %>%
  mutate(t1 = factor(t1, levels = traits),
         t2 = factor(t2, levels = traits),
         ri = as.integer(t1), ci = as.integer(t2)) %>%
  filter(ci >= ri)

p_cor <- ggplot(cor_long, aes(x = t2, y = fct_rev(t1), fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r),
                color = abs(r) > 0.5), size = 3.2) +
  scale_fill_gradient2(low = "#D7191C", mid = "white", high = "#2C7BB6",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Pearson r") +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "black"),
                     guide = "none") +
  scale_x_discrete(position = "top") +
  labs(#title    = "Phenotypic correlations among traits",
       #subtitle = "Based on genotype means across all environments",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 0, face = "bold"),
        panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold"),legend.position = "none")
ggsave("Results/Fig3.3b_phenotypic_correlations.png", p_cor,
       width = 7, height = 6, dpi = 600)

cat("Section 3.3 figures saved\n")

#combine figures together
Fig3 <- p_violin + p_cor+ 
  plot_layout(tag_level = 'new',widths = c(1.5, 1)) +
  plot_annotation(tag_levels = list(c('a)', 'b)'))) & 
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))


ggsave("Figures/Fig3_trait variation.png", Fig3,
       width = 9, height = 5, dpi = 600)

# =============================================================================
# SECTION 3.4 — GENOTYPE × ENVIRONMENT INTERACTION
# =============================================================================

# ── Fig 3.4a: GxE heatmap ordered by release year (yield example) -----------
blups_best <- read.csv("Results/combined_BLUPs_best_model.csv")

release_yr <- read.csv("data/lentil_blups_wide_spats.csv") %>%
  rename(genotype = Genotype)%>%
  select(genotype, dev_year) %>% distinct()

# Per-env BLUPs (all traits)
blups_long <- read.csv("Results/blups_SpATS_all.csv") %>%
  bind_rows(read.csv("Results/blups_GIBD_all.csv")) %>%
  semi_join(read.csv("Results/best_model_selection.csv"),
            by = c("trait","env","model" = "best_model"))

for (tr in traits) {
  mat <- blups_long %>%
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
    mutate(genotype = reorder(genotype, dev_year))
  
  p_hm <- ggplot(mat_long,
                 aes(x = env, y = genotype, fill = z_BLUP)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#D7191C", mid = "white", high = "#2C7BB6",
                         midpoint = 0, name = "Std. BLUP") +
    labs(title    = paste0("GxE heatmap — ", tr),
         subtitle = "Genotypes ordered bottom-to-top by release year",
         x = "Environment", y = "Genotype (ordered by release year)") +
    theme_bw(base_size = 9) +
    theme(axis.text.y      = element_text(size = 4),
          axis.text.x      = element_text(angle = 30, hjust = 1),
          panel.grid       = element_blank(),
          plot.title       = element_text(face = "bold"))
  ggsave(paste0("Results/Fig3.4a_GxE_heatmap_", tr, ".png"), p_hm,
         width = 6, height = 14, dpi = 300)
}

# ── Fig 3.4b: Finlay-Wilkinson b distribution per trait ----------------------
# Traits where all environments produced nearly identical means (SD of I_j
# below the min_Ij_sd threshold) return NULL from finlay_wilkinson() and are
# excluded from the plot with a diagnostic message.

fw_all <- read.csv("Results/FW_stability_all_traits.csv") %>%
  filter(trait %in% traits,
         b_class != "Unreliable",    # drop numerical artefacts (|b| > 10)
         !is.na(b)) %>%
  mutate(trait = factor(trait, levels = traits))

# Report which traits were dropped
skipped <- setdiff(traits, unique(as.character(fw_all$trait)))
if (length(skipped) > 0)
  cat("FW plot — traits skipped (negligible environment index variance):\n  ",
      paste(skipped, collapse = ", "), "\n")

# Only plot traits that have valid slopes
fw_plot_traits <- intersect(traits, unique(as.character(fw_all$trait)))

if (length(fw_plot_traits) == 0) {
  cat("No traits with valid FW slopes — skipping Fig 3.4b\n")
} else {
  fw_plot <- fw_all %>%
    filter(trait %in% fw_plot_traits) %>%
    mutate(trait = factor(trait, levels = fw_plot_traits))
  
  # How many columns? Keep 3 unless fewer traits remain
  n_cols <- min(3, length(fw_plot_traits))
  
  p_fw <- ggplot(fw_plot, aes(x = b, fill = b_class)) +
    geom_histogram(bins = 20, color = "white", linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey30", linewidth = 0.7) +
    facet_wrap(~ trait, scales = "free_x",   # free x only — keep y comparable
               ncol = n_cols) +
    scale_fill_manual(
      values = c(
        "Stable / low-input adapted"     = "#2C7BB6",
        "Average stability"               = "#FFFFBF",
        "Responsive / high-input adapted" = "#D7191C"
      ),
      name = "Stability class"
    ) +
    labs(
      title    = "Finlay-Wilkinson regression slope distribution",
      subtitle = paste0("Dashed line = b = 1 (average stability) | ",
                        "Traits with negligible G\u00D7E excluded"),
      x = "Regression slope (b)",
      y = "Count"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "#F0F0F0"),
      legend.position  = "top",
      plot.title       = element_text(face = "bold")
    )
  
  ggsave("Results/Fig3.4b_FW_slope_distribution.png", p_fw,
         width = 4 * n_cols, height = 3 * ceiling(length(fw_plot_traits) / n_cols),
         dpi = 300)
  cat("Fig 3.4b saved —", length(fw_plot_traits), "traits plotted\n")
}

# ── Fig 3.4c: Mean performance vs. WAASB (yield, per env) --------------------
stab_df <- read.csv("Results/stability_indices_YLD.csv") %>%
  left_join(release_yr %>% mutate(GEN = as.character(genotype)) %>%
              select(GEN, dev_year), by = "GEN")

if ("WAASB" %in% names(stab_df) && "Y" %in% names(stab_df)) {
  grand_mean_yield <- mean(stab_df$Y, na.rm = TRUE)
  grand_waasb      <- mean(stab_df$WAASB, na.rm = TRUE)
  
  p_quad <- ggplot(stab_df,
                   aes(x = WAASB, y = Y, color = dev_year)) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_vline(xintercept = grand_waasb,
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = grand_mean_yield,
               linetype = "dashed", color = "grey50") +
    geom_text_repel(data = stab_df %>%
                      filter(Y > quantile(Y, 0.85, na.rm=TRUE) |
                               WAASB < quantile(WAASB, 0.15, na.rm=TRUE)),
                    aes(label = GEN), size = 2.5, max.overlaps = 15) +
    scale_color_viridis_c(name = "Release\nyear") +
    annotate("text", x = -Inf, y = Inf,
             label = "High mean\nHigh stability", hjust = -0.1,
             vjust = 1.3, size = 3, color = "#2C7BB6", fontface = "italic") +
    annotate("text", x = Inf, y = Inf,
             label = "High mean\nLow stability", hjust = 1.1,
             vjust = 1.3, size = 3, color = "#D7191C", fontface = "italic") +
    labs(title    = "Mean performance vs. stability (WAASB) — yield",
         subtitle = "Dashed lines = grand means | Colour = release year",
         x = "WAASB (lower = more stable)",
         y = "Mean yield (g plot\u207B\u00B9)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave("Results/Fig3.4c_mean_vs_WAASB_yield.png", p_quad,
         width = 9, height = 7, dpi = 300)
}

#??? ── Table 3.4: Stability indices summary (top 20 genotypes by MGIDI) ---------
mgidi <- read.csv("Results/MGIDI_multi_trait_ranking.csv") %>%
  inner_join(release_yr, by = "genotype")

write.csv(mgidi, "Results/Table3.4_MGIDI_ranking.csv", row.names = FALSE)

cat("Section 3.4 figures saved\n")

# =============================================================================
# SECTION 3.5 — CLIMATE DRIVERS OF GxE
# =============================================================================

# ── Fig 3.5a: PLS loadings — climate drivers of AMMI IPCA1 per trait --------
library(tidytext)
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY/Results")
pls_files <- list.files(pattern = "^PLS_loadings_.*\\.csv$")
if (length(pls_files) > 0) {
  pls_all <- bind_rows(lapply(pls_files, read.csv)) %>%
    filter(trait %in% traits) %>%
    rename(load_PC1 = matches("Comp.1|comp_1|Comp 1")) %>%
    group_by(trait) %>%
    slice_max(abs(load_PC1), n = 10) %>%
    ungroup() %>%
    mutate(
      direction = ifelse(load_PC1 > 0, "Positive", "Negative"),
      variable  = str_replace_all(variable, "_", " "),
      variable  = reorder_within(variable, load_PC1, trait)
    )
  
  p_pls <- ggplot(pls_all,
                  aes(x = load_PC1, y = variable, fill = direction)) +
    geom_col(color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    facet_wrap(~ trait, scales = "free_y", nrow = 3) +
    scale_y_reordered() +
    scale_fill_manual(values = c("Positive" = "#2C7BB6",
                                 "Negative" = "#D7191C"),
                      guide = "none") +
    labs(title    = "PLS loadings: climate drivers of AMMI IPCA1",
         subtitle = "Top 8 climate covariates per trait",
         x = "PLS loading (Component 1)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "#F0F0F0"),
          plot.title       = element_text(face = "bold"))
  ggsave("Fig3.5a_PLS_loadings.png", p_pls,
         width = 13, height = 10, dpi = 300)
}
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY")
# ── Fig 3.5b: Random Forest importance heatmap across traits -----------------
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY/Results")
rf_files <- list.files(pattern = "^RF_importance_.*\\.csv$")
if (length(rf_files) > 0) {
  rf_all  <- bind_rows(lapply(rf_files, read.csv)) %>%
    filter(trait %in% traits)
  top10   <- rf_all %>% group_by(trait) %>%
    slice_max(IncMSE, n = 10) %>% pull(variable) %>% unique()
  
  p_rf <- rf_all %>%
    filter(variable %in% top10) %>%
    mutate(trait    = factor(trait, levels = traits),
           variable = str_replace_all(variable, "_", " ")) %>%
    ggplot(aes(x = trait, y = variable, fill = IncMSE)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = round(IncMSE, 1)), size = 2.8) +
    scale_fill_viridis_c(name = "%IncMSE", option = "plasma") +
    labs(title    = "Random Forest variable importance — climate drivers of IPCA1",
         subtitle = "Top 8 variables per trait | % increase in MSE on permutation",
         x = "Trait", y = "Climate covariate") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
  ggsave("Fig3.5b_RF_importance.png", p_rf,
         width = 10, height = 7, dpi = 300)
}
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myY")
cat("Section 3.5 figures saved\n")


# =============================================================================
# SECTION 3.6 — GENETIC GAIN OVER TIME
# =============================================================================

combined_blups <- read.csv("Results/combined_BLUPs_best_model.csv") %>%
  filter(trait %in% traits) %>%
  mutate(trait = factor(trait, levels = traits))

# ── Table 3.6: Genetic gain estimates ----------------------------------------
gain_results <- read.csv("Results/genetic_gain_estimates.csv") %>%
  filter(trait %in% traits) %>%
  mutate(trait = factor(trait, levels = traits)) %>%
  arrange(trait) %>%
  select(trait, n_lines, yr_range, slope, slope_se,
         gain_pct_yr, R2, p_value, sig_label)
write.csv(gain_results, "Results/Table3.6_genetic_gain_estimates.csv",
          row.names = FALSE)

# ── Fig 3.6a: Genetic gain regression — all traits faceted ------------------
p_gain <- ggplot(combined_blups,
                 aes(x = dev_year, y = BLUP_combined)) +
  geom_point(aes(color = dev_year), size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2C7BB6", fill = "#AED6F1", linewidth = 0.9) +
  facet_wrap(~ trait, scales = "free_y", nrow = 3) +
  geom_text(
    data = gain_results,
    aes(x = -Inf, y = Inf,
        label = paste0("\u03B2 = ", round(slope, 3),
                       " (", sig_label, ")\n",
                       "R\u00B2 = ", round(R2, 2))),
    hjust = -0.07, vjust = 1.3, size = 3, inherit.aes = FALSE,
    color = "grey20"
  ) +
  scale_color_viridis_c(name = "Release\nyear", option = "plasma") +
  labs(title    = "Genetic gain in lentil — ACTIVATE rotation trial",
       subtitle = "Multi-environment BLUPs regressed on release year (2006\u20132027)",
       x = "Release year",
       y = "Combined BLUP") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        legend.position  = "right",
        plot.title       = element_text(face = "bold"))
ggsave("Results/Fig3.6a_genetic_gain_all_traits.png", p_gain,
       width = 13, height = 10, dpi = 300)

# ── Fig 3.6b: Stability (WAASB) vs. release year — key traits ----------------
stab_time <- read.csv("Results/stability_vs_release_year.csv")

p_stab_time <- stab_time %>%
  filter(stability_index == "WAASB") %>%
  inner_join(
    bind_rows(lapply(traits, function(tr) {
      f <- paste0("stability_indices_", tr, ".csv")
      if (file.exists(f))
        read.csv(f) %>% mutate(trait = tr) %>%
        rename(GEN = GEN) %>% select(GEN, WAASB, trait)
    })),
    by = "trait"
  ) %>%
  left_join(release_yr %>% mutate(GEN = as.character(genotype)) %>%
              select(GEN, dev_year), by = "GEN") %>%
  filter(!is.na(dev_year), !is.na(WAASB)) %>%
  mutate(trait = factor(trait, levels = traits)) %>%
  ggplot(aes(x = dev_year, y = WAASB, color = trait)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  facet_wrap(~ trait, scales = "free_y") +
  geom_text(
    data = stab_time %>%
      filter(stability_index == "WAASB",
             trait %in% c("yield","DTM","protein")),
    aes(x = -Inf, y = Inf,
        label = paste0("\u03C1 = ", rho, "  p = ", round(p_val, 3))),
    hjust = -0.1, vjust = 1.4, size = 3,
    color = "grey20", inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("yield"   = "#2C7BB6",
                                "DTM"     = "#D7191C",
                                "protein" = "#1A9641"),
                     guide = "none") +
  labs(title    = "Stability (WAASB) vs. release year",
       subtitle = "Spearman \u03C1 tests whether newer lines are more/less stable",
       x = "Release year", y = "WAASB (lower = more stable)") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        plot.title       = element_text(face = "bold"))
ggsave("Results/Fig3.6b_stability_vs_release_year.png", p_stab_time,
       width = 10, height = 7, dpi = 300)

cat("Section 3.6 figures saved\n")
cat("\n\u2714 All Results outputs written to ./Results/\n")
cat("\nFigures:\n")
invisible(lapply(sort(list.files("Results", pattern="\\.png$")),
                 function(f) cat("  \u2022", f, "\n")))
cat("\nTables:\n")
invisible(lapply(sort(list.files("Results", pattern="\\.csv$")),
                 function(f) cat("  \u2022", f, "\n")))