library(cluster)
library(factoextra)
library(tidyverse)
library(sp)
library(rgdal)
library(geosphere)
library(dismo)
library(rgeos)



# load data ---------------------------------------------------------------

dat <- read_csv("data/2012-2020_report.csv") %>% select(-X1)

# site exploration

# RLHA ~~~
dat %>% select(y20v1:y20v3, y12v1:)


# spatial cluster ---------------------------------------------------------

# PEFA territory delineation

dat <- read.csv("Territory Cluster/pefa_territories.csv")
dat <- read.csv("Territory Cluster/rlha_territories.csv")

head(dat)
length(unique(dat$territory))	
# example data from the thread
x <- dat$long
y <- dat$lat

# convert data to a SpatialPointsDataFrame object
xy <- SpatialPointsDataFrame(
  matrix(c(x,y), ncol=2), data.frame(ID=seq(1:length(x))),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# use the distm function to generate a geodesic distance matrix in meters
mdist <- distm(xy)
hist(mdist)
dim(mdist)


# cluster all points using a hierarchical clustering approach
hc <- hclust(as.dist(mdist), method="complete")

# define the distance threshold, in this case 3500 m
d = 1000

# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
xy$clust <- cutree(hc, h=d)
length(unique(xy$clust)) # points are clustered into 50 territories
dim(dat)

dat$territory <- xy$clust # append territory number to all 144 historical nesting sites

head(dat)

write.csv(dat, file = "pefa_territories.csv")

xy@bbox[] <- as.matrix(extend(extent(xy),0.001))

# get the centroid coords for each cluster
cent <- matrix(ncol=2, nrow=max(xy$clust))
for (i in 1:max(xy$clust))
  # gCentroid from the rgeos package
  cent[i,] <- gCentroid(subset(xy, clust == i))@coords

# compute circles around the centroid coords using a 3500m radius
# from the dismo package
ci <- circles(cent, d=d, lonlat=T)

# plot
plot(ci@polygons, axes=T)
plot(xy, col=rainbow(length(unique(xy$clust)))[factor(xy$clust)], add=T)

# example cropping
range(dat$latitude)
range(dat$longitude)

plot(ci@polygons, axes=T, xlim = c(-117, -116), ylim = c(67.8,68))
plot(xy, col=rainbow(200)[factor(xy$clust)], add=T)



# verify territories / rename ---------------------------------------------


