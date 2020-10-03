library(sf)
library(tidyverse)
library(nuwcru)
library(geosphere)


coords  <- read_csv("data/2012-2020_report.csv") %>% select(Lat_DD, Lon_DD)


sites <- st_as_sf(coords, coords = c("Lon_DD", "Lat_DD"))
road <- st_read("gis/bim_road.shp")


dist <- dist2Line(p = st_coordinates(sites), line = st_coordinates(road)[,1:2])
dist[,3] - dat$Lat_DD
dat <- read_csv("data/2012-2020_report.csv") %>% select(-X1)
dat$dist_disturb <- dist[,1]

write.csv(dat, "data/2012-2020_report.csv")
