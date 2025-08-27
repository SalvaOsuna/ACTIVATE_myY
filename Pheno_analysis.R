#Spatial variation models #####
# Required libraries
library(readr)
library(dplyr)
library(SpATS)
library(reshape2)

# Load data
ACTIVATE_lentils_myYraw <- read_csv("data/ACTIVATE_lentils_myYraw.csv")

# Define environments and traits (adjust if needed)
envs <- c("Clavet_2024", "Hunter_2024", "Clavet_2025", "Hunter_2025")
traits <- c("DTE", "DTF", "DTM", "lodging")

# ---------- 1) quick diagnostics on classes ----------
message("Column classes (first 30 columns):")
print(sapply(ACTIVATE_lentils_myYraw[ , intersect(names(ACTIVATE_lentils_myYraw), c(traits, 'Row','Col','Block','Lentil','ENV'))], class))

# ---------- 2) safe type conversions ----------
# Convert genotype to factor, ENV to factor
ACTIVATE_lentils_myYraw <- ACTIVATE_lentils_myYraw %>%
  mutate(
    Lentil = as.factor(Lentil),
    ENV = as.factor(ENV),
    Block = as.factor(Block) # keep blocks as factor
  )

# Convert Row and Col to numeric coordinates if they are not already numeric
# (if they are labeled like "R1","R2" you may need as.numeric(gsub("[^0-9]","",Row)) instead)
ACTIVATE_lentils_myYraw <- ACTIVATE_lentils_myYraw %>%
  mutate(
    Row = suppressWarnings(as.numeric(as.character(Row))),
    Col = suppressWarnings(as.numeric(as.character(Col)))
  )

# Convert trait columns to numeric in a safe way and report any problematic values
for (t in traits) {
  if (t %in% names(ACTIVATE_lentils_myYraw)) {
    # try converting using as.numeric(as.character(.))
    orig <- ACTIVATE_lentils_myYraw[[t]]
    conv <- suppressWarnings(as.numeric(as.character(orig)))
    # show if conversion created NAs where orig was not NA
    newly_na <- sum(is.na(conv) & !is.na(orig))
    if (newly_na > 0) {
      warning(sprintf("Conversion created %d NA(s) in trait '%s' (non-numeric entries). Showing up to 10 unique problematic values:", newly_na, t))
      bad_vals <- unique(as.character(orig[is.na(conv) & !is.na(orig)]))
      print(head(bad_vals, 10))
    }
    ACTIVATE_lentils_myYraw[[t]] <- conv
  } else {
    stop(sprintf("Trait '%s' not found in the dataset.", t))
  }
}

# ---------- 3) containers and modeling loop ----------
spatial_results <- list()
geno_results <- list()

for (env in envs) {
  data_env <- ACTIVATE_lentils_myYraw %>% filter(ENV == env)
  if (nrow(data_env) == 0) {
    warning("No rows for environment: ", env)
    next
  }
  
  # Ensure Row/Col are numeric and Lentil/Block are factors in the subset
  data_env <- data_env %>%
    mutate(
      Row = as.numeric(Row),
      Col = as.numeric(Col),
      Lentil = as.factor(Lentil),
      Block = as.factor(Block)
    )
  
  for (trait in traits) {
    # skip trait if all NA in this env
    if (!(trait %in% names(data_env))) {
      warning("Trait ", trait, " missing in env ", env)
      next
    }
    if (all(is.na(data_env[[trait]]))) {
      warning("All NA for ", trait, " in env ", env)
      next
    }
    # If trait is not numeric (still), skip and warn
    if (!is.numeric(data_env[[trait]])) {
      warning("Trait ", trait, " in env ", env, " is not numeric after conversion. Skipping.")
      next
    }
    
    message("Fitting SpATS: ", env, " - ", trait)
    fit <- tryCatch(
      SpATS(
        response = trait,
        genotype.as.random = TRUE,
        genotype = "Lentil",
        spatial = ~ PSANOVA(Row, Col, nseg = c(15, 8)),
        data = data_env,
        fixed = ~ Block,
        control = list(tolerance = 1e-03)
      ),
      error = function(e) {
        warning("SpATS failed for ", env, " - ", trait, ": ", e$message)
        return(NULL)
      }
    )
    if (is.null(fit)) next
    
    # spatial trend
    sp_trend <- obtain.spatialtrend(fit, grid = c(30, 77))
    df_sp <- melt(sp_trend$fit, varnames = c("Row", "Col"), value.name = "Value")
    df_sp$Trait <- trait
    df_sp$ENV <- env
    spatial_results[[paste(env, trait, "spatial", sep = "_")]] <- df_sp
    
    # genotype predictions
    preds_geno <- tryCatch(
      predict(fit, which = "Lentil", predFixed = "marginal"),
      error = function(e) {
        warning("predict(..., which='geno') failed for ", env, " - ", trait, ": ", e$message)
        return(NULL)
      }
    )
    if (is.null(preds_geno)) next
    
    # detect predicted and SE columns robustly
    pred_col <- names(preds_geno)[grepl("pred", names(preds_geno), ignore.case = TRUE)][1]
    if (is.na(pred_col)) {
      numeric_cols <- names(preds_geno)[sapply(preds_geno, is.numeric)]
      if (length(numeric_cols) >= 1) pred_col <- numeric_cols[1] else stop("Cannot detect predicted column in predict output")
    }
    
    preds_df <- as.data.frame(preds_geno)
    # assume the genotype column is the first non-numeric or the first column if uncertain
    geno_col_name <- names(preds_df)[1]
    # build tidy df
    geno_df <- preds_df %>%
      dplyr::mutate(across(all_of(pred_col), as.numeric)) %>%
      dplyr::select(all_of(c(geno_col_name, pred_col))) %>%
      dplyr::rename(
        Lentil = !!geno_col_name,
        PredictedMean = !!pred_col
      )
    
    # compute grand mean and BLUP
    grand_mean <- mean(geno_df$PredictedMean, na.rm = TRUE)
    geno_df <- geno_df %>%
      mutate(
        BLUP = PredictedMean - grand_mean,
        Trait = trait,
        ENV = env,
        grand_mean = grand_mean
      )
    
    geno_results[[paste(env, trait, sep = "_")]] <- geno_df
  }
}

# Combine results and save
spatial_variation_lentil <- bind_rows(spatial_results)
genotype_means_and_blups <- bind_rows(geno_results)

write.csv(spatial_variation_lentil, "spatial_variation_lentil.csv", row.names = FALSE)
write.csv(genotype_means_and_blups, "genotype_means_and_blups.csv", row.names = FALSE)

# Quick look
head(spatial_variation_lentil)
head(genotype_means_and_blups)


#heat maps based on SpATS trends
# Load required libraries
library(ggplot2)
library(dplyr)

# Get unique combinations of Trait and Environment
trait_env_combos <- spatial_variation_lentil %>%
  distinct(Trait, ENV)

# Function to standardize a subset of data
standardize_subset <- function(df, trait, env) {
  df_subset <- df %>%
    filter(Trait == trait, ENV == env)
  if (nrow(df_subset) > 0 && sd(df_subset$Value, na.rm = TRUE) != 0) {
    df_subset <- df_subset %>%
      mutate(Value_std = (Value - mean(Value, na.rm = TRUE)) / sd(Value, na.rm = TRUE))
  } else {
    df_subset <- df_subset %>%
      mutate(Value_std = NA_real_)  # Handle cases with no variation or empty data
  }
  return(df_subset)
}

# Apply standardization to each Trait-Environment combination
standardized_dfs <- pmap_dfr(
  list(trait_env_combos$Trait, trait_env_combos$ENV),
  ~ standardize_subset(spatial_variation_lentil, ..1, ..2)
)

# Create heatmaps using standardized values
p <- ggplot(standardized_dfs |> filter(Trait %in% c("DTF","lodging")), aes(x = Row, y = Col, fill = Value_std)) +
  geom_tile() +  # Use tiles for heatmap
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,  # Midpoint at 0 for standardized data
                       name = "\nSpatial Trend\n(z-score)") +  # Diverging color gradient
  facet_grid(ENV ~ Trait, scales = "free") +  # Facet by Environment and Trait
  theme_minimal() +  # Clean theme
  theme(
    plot.title = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11),
    legend.position = "none",
    panel.grid = element_blank()  # Remove grid lines for clarity
  ) +
  labs(x = "Row", y = "Column", 
       title = "Spatial Trend Heatmaps by Trait and Location")


# Display the plot
print(p)

# Optionally, save the plot to a file
ggsave("spatial_trend_heatmaps.png", 
       plot = p, units = "in",width = 6, height = 3, dpi = 600, bg = "white")

#GET PREDICTED MEANS#



#linear mixed models####

#Augmented design models####