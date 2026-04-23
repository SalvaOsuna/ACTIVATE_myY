# ============================================================
# EM38 MK2 Survey — Hunter 2025 | Combined Heatmap + Raw Points
# Traits: Elevation, CV_1.0m, CV_0.5m, IV_1.0m, IV_0.5m
# ============================================================

# --- 1. Packages ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, sf, gstat, ggplot2, viridis, scales, patchwork)

# --- 2. Load & Clean Data ---
df <- read.csv("data/Rel_Electro_Hunter2025.csv", header = TRUE)

colnames(df) <- c("Reading", "Seconds_UTC", "Lon", "Lat", "Elevation",
                  "CV_1.0m", "CV_0.5m", "IV_1.0m", "IV_0.5m")

df_clean <- df %>%
  filter(!is.na(Lon), !is.na(Lat))

cat("Rows loaded:", nrow(df_clean), "\n")

# --- 2b. Remove Extreme Outliers (IQR method, per trait) ---
remove_outliers_iqr <- function(df, cols, k = 3) {
  df_out <- df
  removed_total <- 0
  
  for (col in cols) {
    x    <- df_out[[col]]
    Q1   <- quantile(x, 0.25, na.rm = TRUE)
    Q3   <- quantile(x, 0.75, na.rm = TRUE)
    IQR  <- Q3 - Q1
    lo   <- Q1 - k * IQR
    hi   <- Q3 + k * IQR
    
    outliers <- !is.na(x) & (x < lo | x > hi)
    n_out    <- sum(outliers)
    removed_total <- removed_total + n_out
    
    df_out[[col]][outliers] <- NA   # set to NA (row kept for other traits)
    
    cat(sprintf("  %-12s  IQR fence [%.3f, %.3f]  →  %d outlier(s) removed\n",
                col, lo, hi, n_out))
  }
  
  cat(sprintf("  Total outlier values nulled: %d\n\n", removed_total))
  return(df_out)
}

cat("--- Outlier removal (k = 3 × IQR) ---\n")
df_clean <- remove_outliers_iqr(df_clean, cols = traits, k = 3)

# --- 3. Shared settings ---
traits     <- c("Elevation", "CV_1.0m", "CV_0.5m", "IV_1.0m", "IV_0.5m")
grid_res   <- 0.000075                                      # ~5 m; increase to 0.0001 if slow
lat_mean   <- mean(df_clean$Lat)
asp        <- 1 / cos(lat_mean * pi / 180)                # lon/lat aspect correction (~52°N)

trait_labels <- c(
  Elevation = "Elevation (m)",
  CV_1.0m   = "CV 1.0 m (mS/m)",
  CV_0.5m   = "CV 0.5 m (mS/m)",
  IV_1.0m   = "IV 1.0 m (mS/m)",
  IV_0.5m   = "IV 0.5 m (mS/m)"
)

# --- 4. IDW interpolation grid (built once, reused) ---
lon_seq  <- seq(min(df_clean$Lon), max(df_clean$Lon), by = grid_res)
lat_seq  <- seq(min(df_clean$Lat), max(df_clean$Lat), by = grid_res)
grid_df  <- expand.grid(Lon = lon_seq, Lat = lat_seq)
grid_sf  <- st_as_sf(grid_df, coords = c("Lon", "Lat"), crs = 4326)

cat(sprintf("Interpolation grid: %d x %d = %d cells\n",
            length(lon_seq), length(lat_seq), nrow(grid_df)))

# ============================================================
# --- 5. Loop over traits ---
# ============================================================
for (trait in traits) {
  
  cat(sprintf("\nProcessing: %s ...\n", trait))
  
  # -- 5a. Drop NAs for this trait --
  df_trait <- df_clean %>% filter(!is.na(.data[[trait]]))
  
  pts_sf <- st_as_sf(df_trait, coords = c("Lon", "Lat"), crs = 4326)
  
  # -- 5b. IDW --
  idw_model  <- gstat(formula = as.formula(paste(trait, "~ 1")),
                      data = pts_sf, nmax = 12, set = list(idp = 2))
  idw_result <- predict(idw_model, grid_sf)
  
  interp_df <- grid_df %>%
    mutate(value = idw_result$var1.pred)
  
  # -- 5c. Shared colour scale limits (same for both panels) --
  zlim <- range(c(df_trait[[trait]], interp_df$value), na.rm = TRUE)
  
  fill_scale <- scale_fill_viridis_c(
    option  = "turbo",
    limits  = zlim,
    name    = trait_labels[[trait]],
    labels  = label_number(accuracy = 0.01)
  )
  colour_scale <- scale_colour_viridis_c(
    option  = "turbo",
    limits  = zlim,
    name    = trait_labels[[trait]],
    labels  = label_number(accuracy = 0.01),
    guide   = "none"          # legend shown only on heatmap panel
  )
  
  # ---------------------------------------------------------
  # Panel a — IDW Heatmap
  # ---------------------------------------------------------
  p_heat <- ggplot(interp_df, aes(x = Lon, y = Lat, fill = value)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = value), colour = "white", alpha = 0.35,
                 linewidth = 0.3,
                 breaks = pretty(interp_df$value, n = 8)) +
    fill_scale +
    coord_fixed(ratio = asp) +
    labs(
      tag      = "a",
      title    = "IDW Heatmap",
      x        = "Longitude (°)",
      y        = "Latitude (°)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.tag        = element_text(face = "bold", size = 13),
      plot.title      = element_text(face = "bold"),
      legend.position = "right",
      panel.grid      = element_line(colour = "grey88", linewidth = 0.25),
      axis.text       = element_text(size = 8),
      axis.text.x     = element_text(angle = 25, hjust = 1)
    )
  
  # ---------------------------------------------------------
  # Panel b — Raw GPS Points
  # ---------------------------------------------------------
  p_raw <- ggplot(df_trait, aes(x = Lon, y = Lat, colour = .data[[trait]])) +
    geom_point(size = 0.9, alpha = 0.8) +
    colour_scale +
    coord_fixed(ratio = asp) +
    labs(
      tag   = "b",
      title = "Raw GPS Points",
      x     = "Longitude (°)",
      y     = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.tag        = element_text(face = "bold", size = 13),
      plot.title      = element_text(face = "bold"),
      legend.position = "none",
      panel.grid      = element_line(colour = "grey88", linewidth = 0.25),
      axis.text       = element_text(size = 8),
      axis.text.x     = element_text(angle = 25, hjust = 1)
    )
  
  # ---------------------------------------------------------
  # Combine with patchwork
  # ---------------------------------------------------------
  fig <- p_heat + p_raw +
    plot_layout(guides = 'collect') +   
    plot_annotation(
      title    = sprintf("Hunter 2025 EM38 MK2 — %s", trait_labels[[trait]]),
      subtitle = "Fall survey | WGS84 (EPSG:4326) | IDW interpolation (idp = 2, nmax = 12)",
      caption  = "Left: IDW-interpolated surface with contours  |  Right: observed GPS track",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey40", size = 10),
        plot.caption  = element_text(colour = "grey55", size = 8)
      )
    )
  
  # ---------------------------------------------------------
  # Save
  # ---------------------------------------------------------
  fname <- sprintf("Results/Electroconductivity/Hunter2025_%s.png", gsub("\\.", "", trait))
  ggsave(fname, plot = fig, width = 13, height = 6, dpi = 600)
  cat(sprintf("  Saved: %s\n", fname))
  
}

cat("\nAll done! Files written to working directory.\n")


# Real image============================================================
# --- 1. Packages ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, sf, gstat, ggplot2, viridis, scales,
               patchwork, ggmap, ggnewscale, png)

# ── Google Maps API key ────────────────────────────────────
google_API <- read.delim("data/Google Map API key.txt")
register_google(key = google_API)

# --- 2. Load & Clean Data ---
df <- read.csv("data/Rel_Electro_Clavet2024.csv", header = TRUE)

colnames(df) <- c("Reading", "Seconds_UTC", "Lon", "Lat", "Elevation",
                  "CV_1.0m", "CV_0.5m", "IV_1.0m", "IV_0.5m")

df_clean <- df %>% filter(!is.na(Lon), !is.na(Lat))

# --- 3. Outlier removal ---
traits <- c("Elevation", "CV_1.0m", "CV_0.5m", "IV_1.0m", "IV_0.5m")

remove_outliers_iqr <- function(df, cols, k = 3) {
  df_out <- df
  cat("--- Outlier removal (k =", k, "× IQR) ---\n")
  for (col in cols) {
    x   <- df_out[[col]]
    Q1  <- quantile(x, 0.25, na.rm = TRUE)
    Q3  <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lo  <- Q1 - k * IQR;  hi <- Q3 + k * IQR
    out <- !is.na(x) & (x < lo | x > hi)
    df_out[[col]][out] <- NA
    cat(sprintf("  %-12s  [%.3f, %.3f]  → %d removed\n", col, lo, hi, sum(out)))
  }
  return(df_out)
}

df_clean <- remove_outliers_iqr(df_clean, traits, k = 3)

# --- 4. Fetch & SAVE Google satellite basemap ---
# Buffer: increase number to show more field context around GPS tracks
buf <- 0.0010    # ~100 m padding — increase to 0.002 for even more context

lon_c <- mean(df_clean$Lon)
lat_c <- mean(df_clean$Lat)

# zoom 17 = field scale (~1 m/px). Try 16 for wider view, 18 for sharper.
cat("\nFetching Google satellite basemap...\n")
basemap_gg <- get_googlemap(
  center   = c(lon = lon_c, lat = lat_c),
  zoom     = 17,
  maptype  = "satellite",
  # bbox approach: force extent to match your data + buffer
  style    = "feature:all|element:labels|visibility:off"  # clean: no labels
)

# Save the raw basemap as PNG for reuse (no re-fetching needed later)
basemap_path <- "Results/Electroconductivity/Clavet2024_satellite_basemap.png"
ggmap(basemap_gg) +
  theme_void()
ggsave(basemap_path, width = 8, height = 8, dpi = 300)
cat(sprintf("Basemap saved: %s\n", basemap_path))

# Extract basemap extent (lon/lat bbox from ggmap object)
bb      <- attr(basemap_gg, "bb")          # data.frame: ll.lat ll.lon ur.lat ur.lon
map_lon <- c(bb$ll.lon, bb$ur.lon)
map_lat <- c(bb$ll.lat, bb$ur.lat)

cat(sprintf("Basemap extent  Lon [%.6f, %.6f]  Lat [%.6f, %.6f]\n",
            map_lon[1], map_lon[2], map_lat[1], map_lat[2]))

# --- 5. IDW interpolation grid ---
grid_res <- 0.000075
lon_seq  <- seq(min(df_clean$Lon), max(df_clean$Lon), by = grid_res)
lat_seq  <- seq(min(df_clean$Lat), max(df_clean$Lat), by = grid_res)
grid_df  <- expand.grid(Lon = lon_seq, Lat = lat_seq)
grid_sf  <- st_as_sf(grid_df, coords = c("Lon", "Lat"), crs = 4326)

cat(sprintf("Grid: %d × %d = %d cells\n",
            length(lon_seq), length(lat_seq), nrow(grid_df)))

# --- 6. Shared settings ---
trait_labels <- c(
  Elevation = "Elevation (m)",
  CV_1.0m   = "CV 1.0 m (mS/m)",
  CV_0.5m   = "CV 0.5 m (mS/m)",
  IV_1.0m   = "IV 1.0 m (mS/m)",
  IV_0.5m   = "IV 0.5 m (mS/m)"
)

heatmap_alpha <- 0.70
rawpts_size   <- 0.8
rawpts_alpha  <- 0.85

base_theme <- theme_void(base_size = 11) +     # void = no grey background
  theme(
    plot.tag          = element_text(face = "bold", size = 13),
    plot.title        = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position   = "right",
    legend.key.height = unit(2, "cm"),
    axis.text         = element_text(size = 7.5, colour = "grey20"),
    axis.text.x       = element_text(angle = 25, hjust = 1),
    axis.title        = element_text(size = 9),
    panel.border      = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
  )

# ============================================================
# --- 7. Loop ---
# ============================================================
for (trait in traits) {
  
  cat(sprintf("\nProcessing: %s\n", trait))
  
  df_trait <- df_clean %>% filter(!is.na(.data[[trait]]))
  pts_sf   <- st_as_sf(df_trait, coords = c("Lon", "Lat"), crs = 4326)
  
  # IDW
  idw_model  <- gstat(formula = as.formula(paste(trait, "~ 1")),
                      data = pts_sf, nmax = 12, set = list(idp = 2))
  idw_result <- predict(idw_model, grid_sf)
  interp_df  <- grid_df %>% mutate(value = idw_result$var1.pred)
  
  zlim    <- range(c(df_trait[[trait]], interp_df$value), na.rm = TRUE)
  cbreaks <- pretty(interp_df$value, n = 8)
  
  # ── Panel a: Google satellite + IDW heatmap ───────────────
  p_heat <- ggmap(basemap_gg) +             # ggmap handles the satellite base
    
    # semi-transparent IDW surface
    geom_raster(data = interp_df,
                aes(x = Lon, y = Lat, fill = value),
                alpha = heatmap_alpha, interpolate = TRUE) +
    
    # contour lines
    geom_contour(data = interp_df,
                 aes(x = Lon, y = Lat, z = value),
                 colour = "white", alpha = 0.55,
                 linewidth = 0.35, breaks = cbreaks) +
    
    scale_fill_viridis_c(
      option = "turbo",
      limits = zlim,
      name   = trait_labels[[trait]],
      labels = label_number(accuracy = 0.01)
    ) +
    # Zoom extent to data + buffer (bigger than GPS tracks but inside tile)
    coord_cartesian(
      xlim = c(min(df_clean$Lon) - buf * 0.6,
               max(df_clean$Lon) + buf * 0.6),
      ylim = c(min(df_clean$Lat) - buf * 0.6,
               max(df_clean$Lat) + buf * 0.6)
    ) +
    labs(tag = "a", title = "IDW Heatmap",
         x = "Longitude (°)", y = "Latitude (°)") +
    base_theme
  
  # ── Panel b: Google satellite + raw GPS points ────────────
  p_raw <- ggmap(basemap_gg) +
    
    geom_point(data = df_trait,
               aes(x = Lon, y = Lat, colour = .data[[trait]]),
               size = rawpts_size, alpha = rawpts_alpha) +
    
    scale_colour_viridis_c(
      option = "turbo",
      limits = zlim,
      name   = trait_labels[[trait]],
      labels = label_number(accuracy = 0.01),
      guide  = "none"
    ) +
    coord_cartesian(
      xlim = c(min(df_clean$Lon) - buf * 0.6,
               max(df_clean$Lon) + buf * 0.6),
      ylim = c(min(df_clean$Lat) - buf * 0.6,
               max(df_clean$Lat) + buf * 0.6)
    ) +
    labs(tag = "b", title = "Raw GPS Points",
         x = "Longitude (°)", y = NULL) +
    base_theme + theme(legend.position = "none")
  
  # ── Combine ───────────────────────────────────────────────
  fig <- p_heat + p_raw +
    plot_layout(guides = 'collect') +
    plot_annotation(
      title    = sprintf("Clavet 2024 EM38 MK2 — %s", trait_labels[[trait]]),
      subtitle = "Google Maps Satellite | IDW (idp = 2, nmax = 12) | Outliers: k = 3 × IQR",
      caption  = "Map data © Google",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        plot.caption  = element_text(colour = "grey55", size = 7.5)
      )
    )
  
  fname <- sprintf("Results/Electroconductivity/Clavet2024_%s_googlesat.png", gsub("\\.", "", trait))
  ggsave(fname, plot = fig, width = 14, height = 7, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

cat("\nAll done!\n")
