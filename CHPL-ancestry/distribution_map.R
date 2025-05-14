# This script contains the R code to create a distribution map
#
# Authors: Aguillon SM, et al.
# Year: 2025
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: Nature Ecology and Evolution
# DOI: https://doi.org/10.1038/s41559-025-02669-9
#
# Edited date: 22 Jan 2025
#
# Please cite the paper if you use these scripts
#

# load packages
library(tidyverse)
library(sf)
library(basemaps)
library(ggspatial)

### load data ###

# load KML distribution and CHPL point files
xbir_range <- read_sf(dsn="./data/map_data/birchmanni_range.kml")
xcor_range <- read_sf(dsn="./data/map_data/cortezi_range.kml")
CHPL_point <- read_sf(dsn="./data/map_data/CHPL.kml")
rios <- read_sf(dsn="./data/map_data/Rios.kml")

# plot colors
corcol <- "#359D6F"
bircol <- "#3A6BC6"



# transforming the coordinate reference system (CRS) to match the basemap
# https://stackoverflow.com/questions/78343424/plot-basemap-in-ggplot2-basemaps-gglayer-not-working
xbir_range <- xbir_range %>% 
  st_transform(crs=3857) 

xcor_range <- xcor_range %>% 
  st_transform(crs=3857) 

CHPL_point <- CHPL_point %>%
  st_transform(crs=3857)

rios <- rios %>%
  st_transform(crs=3857)

# intersect to get hybrid zone
HZ <- xbir_range %>%
  st_intersection(xcor_range)

# use the ranges to calculate a bounding box
bbox <- st_bbox(bind_rows(xbir_range,xcor_range))
bbox2 <- bbox*c(1.00085,.99,.99915,1.01)

# crop the rivers to the bounding box
bbp = st_as_sfc(st_bbox(bbox2), crs=3857)
rios <- rios %>%
  sf::st_intersection(bbp)



# basemap types
# https://jakob.schwalb-willmann.de/basemaps/index.html
# https://www.arcgis.com/home/group.html?sortField=title&sortOrder=asc&id=702026e41f6641fb85da88efe79dc166#content
#get_maptypes()

# set default map service and type
# using a base map that shows topographic relief
set_defaults(map_service = "esri", map_type = "world_hillshade_dark")

# build the map
# first three lines create the basemap
ggplot() +
  basemaps::basemap_gglayer(bbox2, dpi=600) +
  ggplot2::scale_fill_identity() + 
  ggplot2::coord_sf() +
  # next lines add the layers on top
  # distributions
  geom_sf(data=HZ, fill="#595959", alpha=0.75) +
  geom_sf(data=xbir_range, fill=bircol, alpha=0.5) +
  geom_sf(data=xbir_range, fill=NA, color=bircol, lwd=2) +
  geom_sf(data=xcor_range, fill=corcol, alpha=0.5) +
  geom_sf(data=xcor_range, fill=NA, color=corcol, lwd=2) +
  # rivers
  geom_sf(data=rios, size=1, color="cyan1")+
  # CHPL site location
  geom_sf(data=CHPL_point, size=3, color="black") +
  ggspatial::annotation_scale(width_hint = 0.3) +
  theme_void() # removes everything except the map itself






### plot inset figure of larger geography

inset_region <- ggplot2::map_data("world") %>%
  filter(region=="Mexico") %>%
  filter(is.na(subregion)==TRUE) #removes small islands around the border

inset_region2 <- ggplot2::map_data("state") %>%
  filter(region=="arizona" | region=="new mexico" | region=="texas" | region=="oklahoma" | 
           region=="arkansas" | region=="louisiana" | region=="mississippi" | 
           region=="alabama" | region=="florida" | region=="georgia")

inset_region3 <- ggplot2::map_data("world") %>%
  filter(region=="Guatemala" | region=="Honduras" | region=="El Salvador" | region=="Belize")

CHPL_lat_long <- c(-98.67022002063706,21.21083535036814)

ggplot() + 
  geom_polygon(data=inset_region, aes(x=long,y=lat,group=group), fill="darkgray", color="black", alpha=0.3) +
  geom_polygon(data=inset_region2, aes(x=long,y=lat,group=group), fill="darkgray", color="black", alpha=0.3) +
  geom_polygon(data=inset_region3, aes(x=long,y=lat,group=group), fill="darkgray", color="black", alpha=0.3) +
  geom_point(aes(x=CHPL_lat_long[1], y=CHPL_lat_long[2]), size=5, color="black") +
  theme_void()



# Both figures were saved as PDFs and then combined+edited in Illustrator to make the final panel.
# Needed to remove one river segment in Xbir's range to better fit the illustrations.
# Moved the gray hybrid zone polygon to be higher up in the stack. 
# Focused the inset of the larger map to be a horizontal rectangle centered on the CHPL point.
# Added in a north arrow and labels, moved scale bar onto the map.
