# Code for preparing pinewood survey data into format for jsdmstan
# Source of pinewood data: https://doi.org/10.5285/56a48373-771c-4d4a-8b5a-45ef496c6e55
# Source of climate data: Met Office; Hollis, D.; McCarthy, M.; Kendon, M.; Legg, T. (2022): HadUK-Grid Gridded Climate Observations on a 1km grid over the UK, v1.1.0.0 (1836-2021). NERC EDS Centre for Environmental Data Analysis, 26 May 2022. doi:10.5285/bbca3267dc7d4219af484976734c9527. http://dx.doi.org/10.5285/bbca3267dc7d4219af484976734c9527

library(dplyr)
library(raster, exclude = "select")
library(sf)

site_locs <- read.csv("Scots_Pine_1971_Sites.csv")
str(site_locs)

# Convert to 1km and central point of 1 km
site_locs2 <- site_locs %>%
  mutate(OSGR = paste0(substring(OSGR, 1, 4),
                       substring(OSGR, 6, 7)),
         Easting = 1000*floor(POINT_X/1000)+500,
         Northing = 1000*floor(POINT_Y/1000)+500)
site_locs_sf <- st_as_sf(site_locs2, coords = c("Easting","Northing"),
                         crs = 27700)

# Climate files extraction
clim_folder <- rstudioapi::askForSecret("Path of climate data folder")
clim_files <- list.files(clim_folder)

rainfall_files <- grep("rainfall", clim_files, value = TRUE)

rain <- sapply(rainfall_files,function(x){
  annual_rain <- brick(paste0(clim_folder,x))
  annual_siterain <- extract(annual_rain, site_locs_sf)
  annual_siterain
})

tasmax_files <- grep("tasmax", clim_files, value = TRUE)
summer_tasmax <- sapply(tasmax_files,function(x){
  tasmax <- brick(paste0(clim_folder,x))
  annual_sitetasmax <- extract(tasmax, site_locs_sf)
  annual_sitetasmax[,3,drop=FALSE]
})

tasmin_files <- grep("tasmin", clim_files, value = TRUE)
winter_tasmin <- sapply(tasmin_files,function(x){
  tasmin <- brick(paste0(clim_folder,x))
  annual_sitetasmin <- extract(tasmin, site_locs_sf)
  annual_sitetasmin[,1,drop=FALSE]
})

site_locs2 <- mutate(site_locs2,
                     rainfall = rowMeans(rain),
                     tasmax_summer = rowMeans(summer_tasmax),
                     tasmin_winter = rowMeans(winter_tasmin))
str(site_locs2)

# soil vars
soils <- read.csv("PINEWOODS1971_SOIL_DATA.csv")
str(soils)
soils %>% mutate(LitterDepth = Litter_to - Litter_from,
                 OrganicDepth = Organic_to - Organic_from,
                 MixedDepth = Mixed_to - Mixed_from,
                 LeachedDepth = Leached_to - Leached_from,
                 WeatherDepth = Weather_to - Weather_from) %>%
  filter(!is.na(Soil_PH)) %>%
  select(Soil_PH, contains("Depth")) %>% psych::pairs.panels()

site_env <- select(site_locs2, Site_no = SITEID, OSGR, Easting, Northing,
                   rainfall, tasmax_summer, tasmin_winter) %>%
  full_join(soils %>%
              mutate(LitterDepth = Litter_to - Litter_from) %>%
              select(Site_no, Plot_no, Soil_PH, LitterDepth),
            by = "Site_no") %>%
  mutate(SitePlot = paste0(Site_no, "P",Plot_no))%>%
  filter(SitePlot != "26P16" | !is.na(Soil_PH))

# tree data
trees <- read.csv("PINEWOODS1971_TREE_DATA.csv")
str(trees)
tree_wide <- trees %>% mutate(SitePlot = paste0(Site_no, "P",Plot_no)) %>%
  filter(Dead == "") %>%
  filter(Species != "Unknown") %>%
  filter(Species != "Pinus sp.") %>%
  mutate(Species = recode(Species,
                          "Betula pendula" = "Betula sp.",
                          "Larix kaempferi" = "Larix sp.",
                          "Sorbus aucuparia" = "Sorbus sp.")) %>%
  count(Site_no, Plot_no, SitePlot, Species) %>%
  mutate(Species = gsub(" ","_",Species)) %>%
  mutate(Species = gsub("\\.", "", Species)) %>%
  tidyr::pivot_wider(names_from = Species,
                     values_from = n,
                     values_fill = 0)

common_flora <- c("Calluna vulgaris", "Vaccinium myrtillus", "Potentilla erecta",
                  "Deschampsia flexuosa", "Vaccinium vitis-idaea","Molinia caerulea",
                  "Blechnum spicant", "Betula seedling/sp", "Agrostis canina sens.lat.",
                  "Pteridium aquilinum", "Erica tetralix", "Carex echinata",
                  "Sorbus aucuparia", "Galium saxatile","Narthecium ossifragum",
                  "Luzula multiflora", "Erica cinerea", "Pinus sylvestris",
                  "Viola riviniana/reichenbiana", "Melampyrum pratense",
                  "Polygala serpyllifolia","Succisa pratensis", "Carex panicea",
                  "Carex binervis", "Eriophorum vaginatum")

flora <- read.csv("PINEWOODS1971_GROUND_FLORA.csv")
str(flora)
flora_wide <- flora %>% mutate(SitePlot = paste0(SITE_NO, "P",PLOT_NO)) %>%
  filter(BRC_NAME %in% common_flora) %>%
  rename(Species = BRC_NAME) %>%
  count(SITE_NO, PLOT_NO, SitePlot, Species) %>%
  mutate(Species = gsub(" ","_",Species)) %>%
  mutate(Species = gsub("\\.", "", Species)) %>%
  tidyr::pivot_wider(names_from = Species,
                     values_from = n,
                     values_fill = 0) %>%
  rename(Viola_sp = `Viola_riviniana/reichenbiana`,
         Pinus_sylvestris_seed = Pinus_sylvestris)

combined_data <- full_join(tree_wide, flora_wide, by = c("Site_no" = "SITE_NO",
                                                         "Plot_no" = "PLOT_NO",
                                                         "SitePlot")) %>%
  mutate(across(Pinus_sylvestris:Polygala_serpyllifolia, tidyr::replace_na, 0)) %>%
  mutate(Betula_sp = Betula_sp + `Betula_seedling/sp`,
         Pinus_sylvestris = Pinus_sylvestris + Pinus_sylvestris_seed) %>%
  select(-Pinus_sylvestris_seed, -`Betula_seedling/sp`) %>%
  arrange(SitePlot)
str(combined_data)
sum(combined_data$SitePlot %in% site_env$SitePlot)

site_env <- filter(site_env, SitePlot %in% combined_data$SitePlot) %>%
  arrange(SitePlot)

# Complete jsdmstan format
site_env <- filter(site_env, !is.na(Soil_PH)) %>%
  select(Site_no, Plot_no, SitePlot, OSGR,
         rainfall, tasmax_summer, tasmin_winter, Soil_PH)

pinewood <- filter(combined_data, SitePlot %in% site_env$SitePlot) %>%
  tibble::column_to_rownames("SitePlot") %>%
  select(-Site_no, -Plot_no)
pinewood <- Y[site_env$SitePlot,]
pinewood_env <- site_env

pinewood <- pinewood %>% mutate(Site = pinewood_env$Site_no) %>% group_by(Site) %>%
  summarise(across(everything(), function(x) sum(x > 0))) %>%
  tibble::column_to_rownames("Site") %>%
  select(-Prunus_padus, -Crataegus_monogyna, -Fraxinus_sp, -Populus_tremula,
         -Pseudotsuga_menziesii)
pinewood_env <- pinewood_env %>%
  select(Site_no, OSGR, rainfall:Soil_PH) %>%
  group_by(Site_no, OSGR) %>% summarise(across(.fns = mean),
                                        No_plots = n(), .groups = "drop") %>%
  mutate(rainfall = rainfall/1000 - 1.5,
         tasmax_summer = tasmax_summer - 16,
         tasmin_winter = tasmin_winter + 1,
         Soil_PH = Soil_PH - 4.3,
         No_plots = No_plots - 16)

usethis::use_data(pinewood, overwrite = TRUE)
usethis::use_data(pinewood_env, overwrite = TRUE)
