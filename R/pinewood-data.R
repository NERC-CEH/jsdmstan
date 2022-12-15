#' Vegetation data from native pinewoods in Scotland, 1971
#'
#' This is the tree data from a series of plots across native Scots Pine woodlands
#' in Scotland surveyed in 1971. It contains 16 tree species measured over 265 plots.
#' For further information on this data see Wood and Bunce (2016).
#'
#' This data is made available under the terms of the Open Government Licence.
#'
#' @format
#' A data frame with 265 rows and 16 columns. Each column is a
#'   species, or species aggregates and each row is a plot. The plot identity is marked
#'   in the rownames of the data frame. Every entry is given as the number of
#'   individuals of that species that is observed in within the plot.
#'
#' @references Wood, C. M., & Bunce, R. G. H. (2016). Ecological survey of the native
#'   pinewoods of Scotland 1971. Earth System Science Data, 8(1), 177-189.
#'   https://doi.org/10.5194/essd-8-177-2016
#'
#' @source <https://doi.org/10.5285/56a48373-771c-4d4a-8b5a-45ef496c6e55>
#'
#' @seealso [`pinewood_env`]
"pinewood"

#' Environmental data from native pinewoods in Scotland, 1971
#'
#' This is the environmental data from a series of plots across native Scots Pine
#' woodlands in Scotland surveyed in 1971. It contains the approximate location of
#' each site, the soil pH of each plot plus the average annual rainfall, maximum
#' summer temperature and minimum winter temperature from 1941-1970. All climate
#' variables were calculated based upon the HadUK-Grid climate data (Met Office et
#' al., 2022). For more information on the sites see Wood and Bunce (2016).  This
#' dataset has been limited to the 265 plots across 25 sites that have soil pH
#' measurements.
#'
#' This data is made available under the terms of the Open Government Licence.
#'
#' @format
#' A data frame with 265 rows and 10 columns.
#' \describe{
#'   \item{Site_no}{The site identifier (25 in total)}
#'   \item{Plot_no}{The plot identifier per site}
#'   \item{SitePlot}{The unique plot identifier}
#'   \item{OSGR}{The 1 km grid reference of each site (OSGB)}
#'   \item{Easting}{The Easting of each site (centre position of 1km square, British National Grid)}
#'   \item{Northing}{The Northing of each site (centre position of 1km square, British National Grid)}
#'   \item{rainfall}{The average annual rainfall in mm}
#'   \item{tasmax_summer}{The average maximum temperature in summer in deg Celsius}
#'   \item{tasmin_winter}{The average minimum temperature in winter in deg Celsius}
#'   \item{Soil_PH}{The soil pH of the top 10cm of soil}
#' }
#'
#' @references {
#' Wood, C. M., & Bunce, R. G. H. (2016). Ecological survey of the native pinewoods
#' of Scotland 1971. Earth System Science Data, 8(1), 177-189.
#' https://doi.org/10.5194/essd-8-177-2016
#'
#' Met Office; Hollis, D.; McCarthy, M.; Kendon, M.; Legg, T. (2022): HadUK-Grid
#' Gridded Climate Observations on a 1km grid over the UK, v1.1.0.0 (1836-2021). NERC
#' EDS Centre for Environmental Data Analysis, 26 May 2022.
#' http://dx.doi.org/10.5285/bbca3267dc7d4219af484976734c9527
#'   }
#'
#' @source <https://doi.org/10.5285/56a48373-771c-4d4a-8b5a-45ef496c6e55>
#'
#' @seealso [`pinewood`]
"pinewood_env"
