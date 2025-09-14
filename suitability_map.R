# suitability_analysis_lahore.R
# ArcGIS-like spatial analysis in R: 3+ datasets, distance function, overlay, map export (PNG)
# Example study area: Lahore, Pakistan


# 0. Packages
# 1. Parameters
# 2. Study area bbox
# 3. Helper: safe OSM query
# 4. Download OSM layers
# 5. Project vectors to UTM
# 6. Raster template
# 7. Rasterize + distances
# 8. Distance -> scores
# 9. Weighted overlay + mask water
# 10. Top 5% polygons
# 11. Map
# 12. Save outputs



# 0.
pkgs <- c("sf","osmdata","terra","ggplot2","ggspatial","dplyr","scales")
for(p in pkgs) if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(sf); library(osmdata); library(terra); library(ggplot2); library(ggspatial); library(dplyr); library(scales)

# 1.
center <- c(lon = 74.3436, lat = 31.5204)   # Lahore center
buffer_m <- 15000                           # study radius (15 km)
target_crs <- 32643                         # UTM zone for Lahore
raster_res <- 50                            # raster resolution in meters
maxdist_for_scoring <- 2000                 # cap distance for suitability scoring

# 2.
pt <- st_sfc(st_point(c(center["lon"], center["lat"])), crs = 4326)
pt_prj <- st_transform(pt, crs = target_crs)
study_poly_prj <- st_buffer(pt_prj, dist = buffer_m)
study_poly_wgs84 <- st_transform(study_poly_prj, 4326)
bbox_vec <- st_bbox(study_poly_wgs84)
bbox_for_osm <- c(bbox_vec["xmin"], bbox_vec["ymin"], bbox_vec["xmax"], bbox_vec["ymax"])

# 3.
osm_safe <- function(bbox, key, value = NULL, geom_type = c("lines","polygons","points")) {
  geom_type <- match.arg(geom_type)
  q <- opq(bbox = bbox)
  if (is.null(value)) q <- add_osm_feature(q, key = key) else q <- add_osm_feature(q, key = key, value = value)
  res <- tryCatch(osmdata_sf(q), error = function(e) { message("OSM query failed: ", e$message); return(NULL) })
  if (is.null(res)) return(NULL)
  if (geom_type == "lines") return(res$osm_lines)
  if (geom_type == "polygons") {
    if (!is.null(res$osm_polygons) && nrow(res$osm_polygons) > 0) return(res$osm_polygons)
    if (!is.null(res$osm_multipolygons) && nrow(res$osm_multipolygons) > 0) return(res$osm_multipolygons)
    return(NULL)
  }
  if (geom_type == "points") return(res$osm_points)
}

# 4.
message("Downloading OSM data...")

roads_sf <- osm_safe(bbox_for_osm, key = "highway", value = NULL, geom_type = "lines")
if (!is.null(roads_sf) && "highway" %in% names(roads_sf)) {
  roads_sf <- roads_sf %>% filter(highway %in% c("motorway","trunk","primary","secondary","tertiary","residential","unclassified",
                                                 "primary_link","secondary_link","tertiary_link"))
}
parks_sf <- osm_safe(bbox_for_osm, key = "leisure", value = "park", geom_type = "polygons")
water_sf <- osm_safe(bbox_for_osm, key = "natural", value = "water", geom_type = "polygons")
if (is.null(water_sf)) water_sf <- osm_safe(bbox_for_osm, key = "waterway", value = NULL, geom_type = "lines")

# 5.
study_poly_prj <- st_transform(study_poly_prj, target_crs)
if (!is.null(roads_sf)) roads_prj <- st_transform(roads_sf, target_crs) else roads_prj <- NULL
if (!is.null(parks_sf)) parks_prj <- st_transform(parks_sf, target_crs) else parks_prj <- NULL
if (!is.null(water_sf)) water_prj <- st_transform(water_sf, target_crs) else water_prj <- NULL

clip_to_study <- function(x) {
  if (is.null(x)) return(NULL)
  st_intersection(st_make_valid(x), st_make_valid(study_poly_prj))
}
roads_prj <- clip_to_study(roads_prj)
parks_prj <- clip_to_study(parks_prj)
water_prj <- clip_to_study(water_prj)

# 6.
bbp <- st_bbox(study_poly_prj)
r <- rast(xmin = bbp["xmin"], xmax = bbp["xmax"], ymin = bbp["ymin"], ymax = bbp["ymax"],
          resolution = raster_res, crs = paste0("EPSG:", target_crs))
sf_to_vect <- function(sf_obj) { if(is.null(sf_obj)) return(NULL); terra::vect(sf_obj) }

# 7.
if(!is.null(roads_prj)) {
  v_roads <- sf_to_vect(roads_prj)
  roads_r <- rasterize(v_roads, r, field = 1, background = 0)
  d_roads <- distance(roads_r)
} else { d_roads <- rast(r); values(d_roads) <- NA }

if(!is.null(water_prj)) {
  v_water <- sf_to_vect(water_prj)
  water_r <- rasterize(v_water, r, field = 1, background = 0)
  d_water <- distance(water_r)
} else { water_r <- NULL; d_water <- rast(r); values(d_water) <- NA }

if(!is.null(parks_prj)) {
  v_parks <- sf_to_vect(parks_prj)
  parks_r <- rasterize(v_parks, r, field = 1, background = 0)
  d_parks <- distance(parks_r)
} else { d_parks <- rast(r); values(d_parks) <- NA }

# 8.
dist_to_score <- function(d_raster, maxd = maxdist_for_scoring) {
  s <- (maxd - d_raster) / maxd
  s_clamped <- terra::clamp(s, lower = 0, upper = 1, values = TRUE)
  s100 <- s_clamped * 100
  return(s100)
}
roads_score <- dist_to_score(d_roads)
water_score <- dist_to_score(d_water)
parks_score <- dist_to_score(d_parks)

# 9.
suit_raw <- 0.5 * roads_score + 0.3 * water_score + 0.2 * parks_score
if (!is.null(water_r)) {
  suitability <- mask(suit_raw, water_r, maskvalues = 1, inverse = TRUE)
} else { suitability <- suit_raw }

rng <- terra::minmax(suitability)
suitability <- (suitability - rng[1]) / (rng[2] - rng[1]) * 100

# 10.
vals <- values(suitability, mat = FALSE)
pct95 <- quantile(vals, probs = 0.95, na.rm = TRUE)
top_cells <- suitability >= pct95
top_polys <- as.polygons(top_cells, dissolve = TRUE)
top_polys_sf <- st_as_sf(top_polys)
top_polys_sf <- st_intersection(st_make_valid(top_polys_sf), st_make_valid(study_poly_prj))
top_polys_sf <- top_polys_sf %>% filter(as.numeric(st_area(.)) > 50*50)

# 11.
suit_df <- as.data.frame(suitability, xy = TRUE, na.rm = TRUE)
names(suit_df)[3] <- "suitability"
if (nrow(suit_df) > 50000) {
  suit_df <- suit_df[sample(nrow(suit_df), 50000), ]
}

p <- ggplot() +
  geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_viridis_c(option = "A", na.value = "transparent", name = "Suitability\n(0-100)") +
  geom_sf(data = study_poly_prj, fill = NA, color = "black", size = 0.6) +
  { if(!is.null(water_prj)) geom_sf(data = water_prj, fill = "lightblue", color = NA, alpha = 0.6) } +
  { if(!is.null(parks_prj)) geom_sf(data = parks_prj, fill = "darkgreen", color = NA, alpha = 0.6) } +
  { if(!is.null(roads_prj)) geom_sf(data = roads_prj, color = "gray20", size = 0.3) } +
  geom_sf(data = top_polys_sf, fill = NA, color = "red", size = 0.8) +
  coord_sf(crs = st_crs(32643)) +
  annotation_scale(location = "bl", width_hint = 0.35) +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering) +
  labs(title = "Suitability map: candidate locations for a new tourist facility",
       subtitle = "Study area: Lahore (15 km buffer)\nWeighted overlay of distance to roads, water, parks",
       caption = "Weights: roads 50%, water 30%, parks 20% | Masked water bodies") +
  theme_minimal() +
  theme(legend.position = c(0.88, 0.25),
        plot.title = element_text(face = "bold", size = 14))

# Save
out_png <- "suitability_map_lahore.png"
ggsave(out_png, plot = p, width = 10, height = 8, dpi = 300)
message("Map saved: ", normalizePath(out_png))

# 12.
writeRaster(suitability, filename = "suitability_lahore.tif", overwrite = TRUE)
if (nrow(top_polys_sf) > 0) st_write(top_polys_sf, "top_candidates_lahore.geojson", delete_dsn = TRUE)
message("Finished.")


