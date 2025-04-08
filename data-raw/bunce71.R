# Code for preparing pinewood survey data into format for jsdmstan
# Source of bunce data:
# * tree data: https://doi.org/10.5285/4d93f9ac-68e3-49cf-8a41-4d02a7ead81a
# * soil data: https://doi.org/10.5285/fb1e474d-456b-42a9-9a10-a02c35af10d2
# * site location: https://doi.org/10.5285/d6409d40-58fe-4fa7-b7c8-71a105b965b4
# Source of climate data: Met Office; Hollis, D.; McCarthy, M.; Kendon, M.; Legg, T. (2022): HadUK-Grid Gridded Climate Observations on a 1km grid over the UK, v1.1.0.0 (1836-2021). NERC EDS Centre for Environmental Data Analysis, 26 May 2022. doi:10.5285/bbca3267dc7d4219af484976734c9527. http://dx.doi.org/10.5285/bbca3267dc7d4219af484976734c9527

library(dplyr)
library(terra)
library(sf)

data.dir <- rstudioapi::askForSecret("Data directory")

tree <- read.csv(paste0(data.dir,"Woodlands_Survey_Tree_Diameter_Data_1971_2001.csv"))
site_info <- read.csv(paste0(data.dir,"Woodlands_Survey_Site_Information_1971_2001.csv"))
soil <- read.csv(paste0(data.dir,"Woodlands_Survey_Soil_Data_1971_2001.csv"))

# Convert to locations
site_loc <- site_info %>% select(SITE, EASTING, NORTHING) %>%
  filter(!is.na(EASTING)) %>%
  mutate(EASTING = EASTING*100, NORTHING = NORTHING*100)
site_locs_sf <- st_as_sf(site_loc, coords = c("EASTING","NORTHING"),
                         crs = 27700)

# Climate files extraction
clim_folder <- rstudioapi::askForSecret("Path of climate data folder")
clim_files <- list.files(clim_folder)

rainfall_files <- grep("rainfall", clim_files, value = TRUE)

rain <- sapply(rainfall_files,function(x){
  annual_rain <- rast(paste0(clim_folder,"/",x))
  annual_siterain <- extract(annual_rain, site_locs_sf)
  annual_siterain$rainfall
})

tasmax_files <- grep("tasmax", clim_files, value = TRUE)
summer_tasmax <- sapply(tasmax_files,function(x){
  tasmax <- rast(paste0(clim_folder,"/",x))
  annual_sitetasmax <- extract(tasmax, site_locs_sf)
  annual_sitetasmax$tasmax_3
})

tasmin_files <- grep("tasmin", clim_files, value = TRUE)
winter_tasmin <- sapply(tasmin_files,function(x){
  tasmin <- rast(paste0(clim_folder,"/",x))
  annual_sitetasmin <- extract(tasmin, site_locs_sf)
  annual_sitetasmin$tasmin_1
})

site_loc <- mutate(site_loc,
                   rainfall = rowMeans(rain),
                   tasmax_summer = rowMeans(summer_tasmax),
                   tasmin_winter = rowMeans(winter_tasmin))
str(site_loc)

# soil vars
site_env <- soil %>%
  mutate(SITEPLOT = paste0("S",Site,"P",Plot)) %>%
  select(SITE = Site, SITEPLOT, PH = pH1971, SOM = SOM1971) %>%
  na.omit() %>%
  group_by(SITE) %>%
  summarise(PH = mean(PH), SOM = mean(SOM), NPlots = length(unique(SITEPLOT))) %>%
  left_join(site_loc, by = "SITE") %>%
  arrange(SITE)
anyNA(site_env)

# tree data
tree_genus_site <- tree %>% filter(Amalgam_names != "" & Yr == 1 & Status == "Live") %>%
  mutate(Genus = sapply(strsplit(Amalgam_names," "), "[",1)) %>%
  group_by(Site, Genus, Plot) %>%
  summarise(Count = sum(Count)>0) %>%
  group_by(Site, Genus) %>%
  summarise(Count = sum(Count>0, na.rm = TRUE)) %>%
  tidyr::pivot_wider(names_from = Genus, values_from = Count,
              values_fill = 0) %>%
  arrange(Site) %>%
  tibble::column_to_rownames("Site")

all.equal(as.character(site_env$SITE), rownames(tree_genus_site))

bunce71 <- tree_genus_site
bunce71_env <- site_env %>%
  select(Site = SITE,
         Easting = EASTING,
         Northing = NORTHING,
         Nplots = NPlots,
         pH = PH,
         SOM = SOM,
         Rainfall = rainfall,
         SummerMaxTemp = tasmax_summer,
         WinterMinTemp = tasmin_winter)

all.equal(as.character(bunce71_env$Site), rownames(bunce71))

bunce71_abund <- bunce71
bunce71 <- list(abund = bunce71_abund,
                env = bunce71_env)

usethis::use_data(bunce71, overwrite = TRUE)
