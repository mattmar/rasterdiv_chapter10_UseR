# =============================================================================
# Title       : "Chapter: The rasterdiv package for measuring diversity from space"
# Author      : Matteo Marcantonio
# Date        : 2025-06-18
# Email       : marcantoniomatteo@gmail.com
# Description : rasterdiv applied to Stebbins Cold Canyon (CA, USA) and Macchiarvana (IT). 
# Description : Processes NAIP CIR and PlanetScope NDVI imagery:
#               - Computes NDVI for years 2016, 2020, 2022
#               - Computes Rényi entropy (α = 0,1,2)
#               - Computes Rao’s quadratic entropy (global & area‐based)
#               - Generates helical plots of seasonal NDVI
# Downloading : All the materials necessary to reproduce codes and figures in this chapter can be downloaded from:
#               https://osf.io/pgbmq/?view_only=73fed91d69b6409d95a72b39d6f24e5f 
# =============================================================================

# Load required packages
library(rasterdiv)
library(terra)
library(rasterVis)
library(RColorBrewer)
library(lattice)
library(tidyterra)
library(viridis)
library(ggplot2)
library(ggpubr)
library(gam)
library(mgcv)
library(ggnewscale)

# Set working directory
setwd("./cold_canyon_CIR/")

# Set number of parallel processes (e.g., detectCores()-1)
mynp=4

# -----------------------------------------------------------------------------
# Figure 03: NDVI
# -----------------------------------------------------------------------------

# Read and stretch NDVI rasters to 8-bit
ndvi       <- lapply(list.files(pattern = "*subset.tif$"), rast)
ndvi_8bit  <- lapply(ndvi, terra::stretch, 0, 255, smin = -1, smax = 1)
names(ndvi_8bit) <- paste("NDVI", c(2016, 2020, 2022))

# Plot
ggplot() +
  geom_spatraster(data = rast(ndvi_8bit)) +
  facet_wrap(~ lyr) +
  coord_sf(crs = 3857, expand = FALSE) +
  geom_sf_text(check_overlap = TRUE) +
  scale_fill_hypso_c(palette = "dem_screen", direction = -1) +
  theme(
    strip.text        = element_text(size = 14, face = "bold"),
    legend.position   = "top",
    axis.text.y       = element_text(size = 6),
    axis.text.x       = element_blank(),
    axis.title        = element_blank(),
    legend.text       = element_blank(),
    panel.spacing.x   = unit(1, "lines")
  ) +
  guides(fill = guide_colourbar(title = "NDVI"))

# Save Figure 03
ggsave(
  "../figures/figure03_NDVI_overview.png",
  units  = "cm",
  height = 12,
  width  = 15,
  dpi    = 600
)

# -----------------------------------------------------------------------------
# Figure 04: Rényi entropy (α = 0,1,2)
# -----------------------------------------------------------------------------

# Compute Rényi entropy (≈5 min with np = 100)
Renyi_ndvi <- lapply(
  ndvi_8bit,
  Renyi,
  window       = 17,
  alpha        = c(0, 1, 2),
  na.tolerance = 0.1,
  np           = mynp
)

# Stack and rename
Renyi_plot <- rast(unlist(lapply(Renyi_ndvi, function(x) lapply(x, identity))))
names(Renyi_plot) <- paste0(
  gsub("_.*", "", names(Renyi_plot)),
  ": Alpha = ",
  rep(c(0, 1, 2), 3)
)
# Reorder layers so α = 0,1,2 for each year
Renyi_plot <- select(Renyi_plot, c(1, 4, 7, 2, 5, 8, 3, 6, 9))

# Alternatively load the file already processed
Renyi_plot <- readRDS("../results/Renyi.plot.RDS")

# Plot
ggplot() +
  geom_spatraster(data = Renyi_plot) +
  facet_wrap(~ lyr) +
  coord_sf(crs = 3857, expand = FALSE) +
  geom_sf_text(check_overlap = TRUE) +
  scale_fill_viridis(option = "B", direction = -1, begin = 0.15) +
  theme(
    strip.text        = element_text(size = 8, face = "bold"),
    legend.position   = "top",
    axis.text         = element_blank(),
    axis.title        = element_blank(),
    legend.text       = element_blank(),
    legend.title      = element_text(size = 8, face = "bold"),
    panel.spacing     = unit(0.1, "lines"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.3, "cm")
  ) +
  guides(fill = guide_colourbar(title = "Rényi’s entropy"))

# Save Figure 04
ggsave(
  "../figures/figure04_Renyi_overview.png",
  units  = "cm",
  height = 12,
  width  = 15,
  dpi    = 600
)

# -----------------------------------------------------------------------------
# Figure 05: Rao’s quadratic entropy (global)
# -----------------------------------------------------------------------------

# paRao for each NDVI
ndvi_paRao <- lapply(
  ndvi_8bit,
  paRao,
  window       = 17,
  alpha        = 1,
  na.tolerance = 0.1,
  np           = mynp
)

paRao_plot <- rast(unlist(ndvi_paRao))
names(paRao_plot) <- paste0("Year ", c(2016, 2020, 2022))

# Alternatively load the file already processed
paRao_plot <- readRDS("../results/paRao.plot.RDS")

# Plot
ggplot() +
  geom_spatraster(data = paRao_plot) +
  facet_wrap(~ lyr) +
  coord_sf(crs = 3857, expand = FALSE) +
  geom_sf_text(check_overlap = TRUE) +
  scale_fill_viridis(option = "B", direction = -1, begin = 0) +
  theme(
    strip.text        = element_text(size = 8, face = "bold"),
    legend.position   = "top",
    axis.text.y       = element_text(size = 6),
    axis.text.x       = element_blank(),
    axis.title        = element_blank(),
    legend.text       = element_blank(),
    legend.title      = element_text(size = 8, face = "bold"),
    legend.margin     = margin(t = 0, b = -8, unit = "pt"),
    panel.spacing     = unit(0.1, "lines"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.3, "cm")
  ) +
  guides(fill = guide_colourbar(
    title          = "Rao’s quadratic entropy",
    title.position = "top"
  ))

# Save Figure 05
ggsave(
  "../figures/figure05_paRao_overview.png",
  units  = "cm",
  height = 12,
  width  = 15,
  dpi    = 600
)

# -----------------------------------------------------------------------------
# Figure 06: Rao’s quadratic entropy (area‐based)
# -----------------------------------------------------------------------------

# Load polygoal altitude layer for isoipses

iso <- vect("coldCanyon_usgs_iso30m.gpkg")
ndvi_paRaoArea <- lapply(
  ndvi_8bit,
  paRao,
  area         = iso,
  field        = "alt_min",
  alpha        = 1,
  na.tolerance = 0.1
)

paRao_plotArea       <- vect(unlist(ndvi_paRaoArea))
paRao_plotArea$lyr   <- rep(paste0("Year ", c(2016, 2020, 2022)), each = 6)
ndvi_plot            <- rast(ndvi_8bit)
names(ndvi_plot)     <- paste0("Year ", c(2016, 2020, 2022))

# Plot
ggplot() +
  geom_spatraster(data = ndvi_plot) +
  facet_wrap(~ lyr) +
  scale_fill_hypso_c(palette = "dem_screen", direction = -1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_spatvector(data = paRao_plotArea, aes(fill = alpha.1)) +
  scale_fill_viridis(option = "B", direction = -1, alpha = 0.7) +
  geom_spatvector_text(
    data      = paRao_plotArea,
    aes(label = alt_max, size = alt_max),
    fontface  = "bold",
    colour    = "red"
  ) +
  scale_size(range = c(2, 4)) +
  coord_sf(crs = 3857, expand = FALSE) +
  theme(
    strip.text        = element_text(size = 12, face = "bold"),
    legend.position   = "top",
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_blank(),
    axis.title        = element_blank(),
    legend.text       = element_blank(),
    legend.title      = element_text(size = 8, face = "bold"),
    legend.margin     = margin(t = 0, b = -8, unit = "pt"),
    panel.spacing     = unit(0.1, "lines"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.4, "cm")
  ) +
  guides(
    size = guide_legend(
      title          = "Altitude (m)",
      title.position = "top"
    ),
    fill = guide_colourbar(
      title          = "Rao’s entropy",
      title.position = "top"
    )
  )

# Save Figure 06
ggsave(
  "../figures/figure06_paRao_area.png",
  units  = "cm",
  height = 14,
  width  = 18,
  dpi    = 600,
  scale  = 1.2
)

# -----------------------------------------------------------------------------
# Figure 07: Multidimensional parametric Rao’s index
# -----------------------------------------------------------------------------

# Load a 3-band CIR image
CIR_bands <- rast("coldCanyon_subset_3Bands2016.tif", lyrs = 1:3)
names(CIR_bands) <- c("NIR", "Red", "Green")
CIR_crop <- crop(CIR_bands, ext(-13592400, -13592300, 4648700, 4648800))

# Plot bands 1 & 2
p1 <- ggplot() +
  geom_spatraster(data = select(CIR_crop, 1:2)) +
  coord_sf(crs = 3857, expand = FALSE) +
  facet_wrap(~ lyr) +
  scale_fill_gradientn(colors = gray.colors(100, start = 0.2, end = 0.8, rev = TRUE)) +
  theme(
    strip.text        = element_text(size = 8, face = "bold"),
    legend.position   = "top",
    axis.text         = element_blank(),
    axis.title        = element_blank(),
    legend.text       = element_blank(),
    legend.title      = element_text(size = 8, face = "bold"),
    panel.spacing     = unit(0.1, "lines"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.3, "cm")
  ) +
  guides(fill = guide_colourbar(title = "DN", title.position = "top"))

# Compute and plot multidimensional Rao
bands_paRaoMulti <- paRao(
  x           = select(CIR_crop, 1:2),
  window      = 9,
  alpha       = 1,
  np          = mynp,
  na.tolerance= 0.1,
  method      = "multidimension",
  rescale     = TRUE,
  dist_m      = "manhattan"
)

# Alternatively load the file already processed
bands_paRaoMulti <- readRDS("../results/paRao.Multi.RDS")

p2 <- ggplot() +
  geom_spatraster(data = bands_paRaoMulti[[1]][[1]]) +
  coord_sf(crs = 3857, expand = FALSE) +
  scale_fill_viridis(option = "B", direction = -1, begin = 0) +
  theme(
    strip.text        = element_text(size = 8, face = "bold"),
    legend.position   = "top",
    axis.text         = element_blank(),
    axis.title        = element_blank(),
    legend.text       = element_blank(),
    legend.title      = element_text(size = 8, face = "bold"),
    panel.spacing     = unit(0.1, "lines"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.3, "cm")
  ) +
  guides(fill = guide_colourbar(
    title          = "Multidimensional\nRao’s quadratic entropy",
    title.position = "top"
  ))

# Combine and save
ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(2, 1.08))
ggsave(
  "../figures/figure07_paRao_multi.png",
  units  = "cm",
  height = 14,
  width  = 18,
  dpi    = 600
)

# -----------------------------------------------------------------------------
# Example (no figure): AUC of Rao’s parametric function
# -----------------------------------------------------------------------------

# Takes a long time (16 mins with 100 cores, possibly hours with 4 cores)
ndvi_paRaoAUC <- lapply(
  ndvi_8bit,
  RaoAUC,
  alphas       = 1:5,
  window       = 17,
  dist_m       = "euclidean",
  na.tolerance = 0.1,
  np           = mynp,
  simplify     = FALSE
)

paRaoAUC.plot <- unlist(ndvi.paRaoAUC)
names(paRaoAUC.plot) <- paste0("Year ",c(2016,2020,2022))

# Alternatively load the file already processed
paRaoAUC.plot <- unwrap(readRDS("results/RaoAUC.alphas.RDS"))

ggplot() +
geom_spatraster(data = paRaoAUC.plot) +
facet_wrap(~lyr) +
coord_sf(crs=3857,expand = FALSE) +
geom_sf_text(check_overlap=TRUE) +
scale_fill_viridis(option = "B", direction=-1, begin=0) + 
theme(
  strip.text = element_text(size = 8, color = "black", face = "bold"),
  legend.position="top",
  axis.text.y=element_text(size = 6),
  axis.text.x=element_blank(),
  axis.title=element_blank(),
  legend.text=element_blank(),
  legend.title=element_text(size=8, color = "black", face = "bold"),
  legend.margin=margin(t = 0, r = 0, b = -8, l = 0, unit = "pt"),
  panel.spacing.x = unit(0.1, "lines"),
  panel.spacing.y = unit(0.1, "lines"),
  legend.key.height= unit(0.4, 'cm'),
  legend.key.width= unit(0.3, 'cm')) +
guides(fill = guide_colourbar(
  title = "Rao's quadratic entropy", title.position = "top"
  ))

# -----------------------------------------------------------------------------
# Figure 08: Helical plots of NDVI seasonality
# -----------------------------------------------------------------------------
# Unfortunately, I cannot share this dataset which is owned by Planet, I've only a Research licence. For more information contact me at: marcantoniomatteo@gmail.com

setwd("../macchiarvana_NDVI")
files   <- list.files(pattern = "*.tif")
dates   <- as.Date(sub("^(.{8}).*$", "\\1", files), "%Y%m%d")
week    <- as.numeric(format(dates, "%W"))
ndvi_r  <- lapply(files, terra::rast)
ndvi_avg <- as.numeric(sapply(ndvi_r, function(x) terra::global(x, mean, na.rm = TRUE)))

# Fit GAMs
outGAM1 <- gam(ndvi_avg ~ te(as.numeric(dates), week, bs = c("cr", "ps")),
               data = data.frame(dates, week, ndvi_avg))

# Predict over continuous dates
all_dates <- seq(min(dates), max(dates), by = 1)
all_week  <- as.numeric(format(all_dates, "%W"))
NDVI_pred1 <- predict(outGAM1, newdata = data.frame(dates = all_dates, week = all_week))

# Prepare for helical plot
NDVI_df <- data.frame(
  date =     all_dates,
  ndvi =     NDVI_pred1,
  year =     as.numeric(format(all_dates, "%Y")),
  week =     as.numeric(format(all_dates, "%W"))
)
NDVI_df <- subset(NDVI_df, year > 2017 & year < 2022)
NDVI_agg <- aggregate(ndvi ~ week + year, NDVI_df, mean)
NDVI_agg$date <- as.Date(paste(NDVI_agg$year, NDVI_agg$week, 1, sep = "-"),
                          "%Y-%W-%u")
# Plot time series
ggplot(NDVI_agg, aes(x = date, y = ndvi)) +
  geom_line() +
  labs(x = "Date", y = "NDVI weekly average") +
  theme_minimal()

# Helical plots
NDVI_prep <- heliPrep(NDVI_agg$date, NDVI_agg$ndvi, filterWidth = 3)
NDVI_prep$year <- format(NDVI_prep$date, "%Y")

heliPlot(NDVI_prep,
         facet        = TRUE,
         group        = "year",
         arrow        = FALSE,
         labelInterval= 90,
         sizeRange    = c(0.05, 2),
         facetScales  = "fixed")

heliPlot(NDVI_prep,
         facet        = FALSE,
         group        = "year",
         arrow        = FALSE,
         labelInterval= 24,
         sizeRange    = c(0.05, 3),
         ylab         = "NDVI weekly average",
         xlab         = "NDVI rate of weekly change") +
  theme(legend.position = c(1, 0.15)) +
  guides(color = guide_legend(title = "Year"))

ggsave(
  "../figures/figure08_NDVI_heliPlot.png",
  units  = "cm",
  height = 14,
  width  = 18,
  dpi    = 600
)
