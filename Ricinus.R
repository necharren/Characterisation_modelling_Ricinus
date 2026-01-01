
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


library(raster)
library(exactextractr)
library(readr)
library(sp)
library(tidyverse)
library(geodata)

getwd()
setwd("C:/Your_directory")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## ### ### ### ### PREPARING THE DATA FOR ANALYSIS ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


# We obtain the bioclimatic variables
worldclim <- worldclim_global(var = "bio", res = 2.5)

#We determine the occurrences of Ricinus
## For the Iberian Peninsula
iber_occ <- read.csv("C:/Your_directory/ric_iber_nic.csv", header = TRUE)
str(iber_occ)

## For the northeast of Africa
native_occ <- read.csv("C:/Your_directory/ric_nat.csv", header = TRUE)
str(native_occ)


# We are going to transform the raster information of the bioclimatic variables into a table so that we can link them to the occurrences
# We extract spatial points from Ricinus occurrences.

## For the Iberian Peninsula
iber_occ_spatial <- SpatialPoints(iber_occ[, c("decimalLon", "decimalLat")])

## For the northeast of Africa
native_occ_spatial <- SpatialPoints(native_occ[, c("decimalLongitude", "decimalLatitude")])

# Upload raster files
bio_rasters <- c("C:/Your_directory/wc2-5/bio1.bil", "C:/Your_directory/wc2-5/bio2.bil", "C:/Your_directory/wc2-5/bio3.bil", "C:/Your_directory/wc2-5/bio4.bil", "C:/Your_directory/wc2-5/bio5.bil","C:/Your_directory/wc2-5/bio6.bil", "C:/Your_directory/wc2-5/bio7.bil", "C:/Your_directory/wc2-5/bio8.bil","C:/Your_directory/wc2-5/bio9.bil", "C:/Your_directory/wc2-5/bio10.bil", "C:/Your_directory/wc2-5/bio11.bil","C:/Your_directory/wc2-5/bio12.bil", "C:/Your_directory/wc2-5/bio13.bil", "C:/Your_directory/wc2-5/bio14.bil","C:/Your_directory/wc2-5/bio15.bil", "C:/Your_directory/wc2-5/bio16.bil", "C:/Your_directory/wc2-5/bio17.bil","C:/Your_directory/wc2-5/bio18.bil", "C:/Your_directory/wc2-5/bio19.bil")

# Create a list of raster objects from file names
rasters <- lapply(bio_rasters, raster)

# Extract values from bioclimatic variables at castor bean locations
ricino_iber <- data.frame(iber_occ_spatial) # We create an empty dataframe to store the IP values.

detach("package:tidyr", unload = TRUE) #We deactivated it because there is a conflict with another function called "extract".

for (i in 1:length(rasters)) {
  values <- extract(rasters[[i]], iber_occ_spatial) # Extract the values of variable i at the castor locations
  col_name <- paste0("bio", i) # Name of the column corresponding to variable i
  ricino_iber[col_name] <- values # Add the values to the dataframe
}

write.csv(ricino_iber, file = "C:/Your_directory/ricino_iber.csv", row.names = FALSE)

library(dplyr)

# Repeat it for native 
ricino_nat <- data.frame(native_occ_spatial)

for (i in 1:length(rasters)) {
  values <- extract(rasters[[i]], native_occ_spatial) 
  col_name <- paste0("bio", i) 
  ricino_nat[col_name] <- values 
}

write.csv(ricino_nat, file = "C:/Your_directory/ricino_nat.csv", row.names = FALSE)


# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


library(maptools)
data(wrld_simpl)

#Checkpoint 
head(ricino_nat)
head(ricino_iber)


# plot bioclimatic variable 6 corresponding to the minimum temperature of the coldest month. 


bio6 <- raster("C:/Your_directory/22/wc2-5/bio6.bil")
plot(bio6) # plot bio6 worldwide
plot(wrld_simpl, add=TRUE) #add map limits


## We will exclude all castor plants that are less than 4º. 
# All temperature-related biovariables are above 10 (10 = 1ºC, 15 = 1.5º).

# Let's establish the temperature limits for growth. 
Ric_NoFrost_iber <- subset(ricino_iber, bio6 > 40) 
head(Ric_NoFrost_iber)


# Load map data for Spain and Portugal
data(wrld_simpl)
spain_portugal <- wrld_simpl[wrld_simpl$NAME %in% c("Spain", "Portugal"), ]


# Plot the bio6 map for Spain and Portugal only.
plot(bio6, xlim = c(-10, 5), ylim = c(35, 45))


# We plotted the points on the map, with black representing discarded occurrences (<4°C) and red representing selected occurrences (>4°C).
points(ricino_iber[,1:2], col='black', pch=20, cex=0.5)
points(Ric_NoFrost_iber[,1:2], col='red', pch=20, cex=0.5) #subset we retain


summary(Ric_NoFrost_iber)

str(Ric_NoFrost_iber)

write.table(Ric_NoFrost_iber[,1:21], file="Ric_NoFrost_iber.txt", col.names=T, row.names=F, sep="\t")
ric_iber <- read.table("Ric_NoFrost_iber.txt", header=TRUE, sep="\t") # all CLEAN occurrences



# Ricinus maps for the Iberian Peninsula and north-eastern Africa.

# Load parcels 

library(maptools)
library(ggplot2)
library(ggspatial)
library(rnaturalearth)
library(ggpubr)

# Load geographical data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Iberian Peninsula 
iberica<- ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitud", y = "Latitud") +
  coord_sf(xlim = c(-10.00, 5), ylim = c(36.00, 44.00), expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.4) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()

# Add a point layer to the existing chart
iberica_with_points <- iberica +
  geom_point(data = ric_iber, aes(x = decimalLon, y = decimalLat), color = "red", size = 0.5)
print(iberica_with_points)


# Native environment (northeast Africa)
Native <- ggplot(data = world) +
  geom_sf() +
  labs(x = "Longitud", y = "Latitud") +
  coord_sf(xlim = c(33.00, 51), ylim = c(-5.00, 18.00), expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.5) +  
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +  
  theme_bw()

# Add a point layer to the existing chart
Native_with_points <- Native +
  geom_point(data = ricino_nat, aes(x = decimalLongitude, y = decimalLatitude), color = "blue", size = 0.5)
print(Native_with_points)




ggarrange(iberica_with_points, Native_with_points, ncol = 2, nrow = 1)




# Combined map

global <- ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitud", y = "Latitud") +
  coord_sf(xlim = c(-180.00, 180.00), ylim = c(-90.00, 90.00), expand = FALSE) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()

combined_ric <- global +
  geom_point(data = ricino_nat, aes(x = decimalLongitude, y = decimalLatitude), color = "blue", size = 0.01) +
  geom_point(data = ric_iber, aes(x = decimalLon, y = decimalLat), color = "red", size = 0.01)

# Display the resulting graph
print(combined_ric)



# Now, let's obtain the background climate values for our areas of interest

# Native area
library(maptools)
data(wrld_simpl)

# We represent the geographical area of interest, in this case, the latitude and longitude coordinates that define the native range
native_range <- extent(33,51, -5, 18)
# create a native_shape object that contains the shape of the selected countries
native_shape<-subset(wrld_simpl,wrld_simpl@data$ISO3=='ETH' | wrld_simpl@data$ISO3=='SOM' | wrld_simpl@data$ISO3=='ERI' | wrld_simpl@data$ISO3=='KEN' | wrld_simpl@data$ISO3=='DJI')

# Here we create the background (bg): We cut the raster area of the bioclimatic variables that match the edges we have previously defined
Ricinus_native_bg<- crop(worldclim, native_shape) #crop raster layers
# Now, we discard all areas that have been left outside the area of interest (mask).
Ricinus_native_bg<- mask(Ricinus_native_bg, native_shape) #crop raster layers
# Let's see how it turned out (e.g. with bio1)
plot(Ricinus_native_bg$bio1)

# Now we are going to transform the information from each raster cell into a CSV file
Ricinus_native_bg <- rasterToPoints(Ricinus_native_bg) 
Ricinus_native_bg <- write.csv(Ricinus_native_bg, "Ricinus_native_bg.csv") #save on disk
Ricinus_native_bg <- read.csv("Ricinus_native_bg.csv") #load


# Repeat the process with the Iberian Peninsula

# Represents the geographical extent of interest, in this case, the latitude and longitude coordinates that delimit the native distribution area
invaded_range <- extent(-10, 5, 34, 46)
# Create a native_shape object that contains the shape of the selected countries 
invaded_shape<-subset(wrld_simpl,wrld_simpl@data$ISO3=='ESP' | wrld_simpl@data$ISO3=='PRT')

# Here we create the background (bg): We cut the raster area of the bioclimatic variables that match the edges we have previously defined
Ricinus_invaded_bg<- crop(worldclim, invaded_shape) #crop raster layers
# Now, we discard all areas that have been left outside the area of interest (mask)
Ricinus_invaded_bg<- mask(Ricinus_invaded_bg, invaded_shape) 
# Let's see how it turned out (e.g. with bio1)
plot(Ricinus_invaded_bg$bio1)
# We delimit the region of interest (Iberian Peninsula)
reg_de_interes <- extent(-10, 5, 34, 46)

# If there's any problem: Crop the layer to the region of interest 
# Ricinus_invaded_bg <- crop(Ricinus_invaded_bg, reg_de_interes)
# plot(Ricinus_invaded_bg$bio6)
# class(Ricinus_invaded_bg)









### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## ### ### ### ### CHARACTERIZATION ANALYSIS ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


library("ade4")
library(ecospat)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(RColorBrewer)
library(maptools)

# Follow the "Ecospat" script to repeat the ecological analyses.

Ricinus_native_bg <- read.csv("C:/Your_directory/Ricinus_native_bg.csv") #load
Ricinus_invaded_bg <- read.csv("C:/Your_directory/iberica_bg.csv") #load

ric_iber <- read.table("Ric_NoFrost_iber.txt", header=TRUE, sep="\t")
ricino_nat <- read.csv("C:/Your_directory/ricino_nat.csv") #load

write.csv(ric_iber, file = "C:/Your_directory/ECOREGIONES/ric_iber.csv", row.names = FALSE)

inv<-ric_iber
nat<-ricino_nat

clim_inv <- Ricinus_invaded_bg
clim_nat <- Ricinus_native_bg

## MERGE ALL THE OBSERVATIONS ##
# Add a "type" (in spanish, "tipo") column to each dataset to distinguish species.
nat$tipo <- "nat"
inv$tipo <- "inv"

#Verify!
print(names(nat))
print(names(inv))

# Rename the columns in the "nat" dataset
names(nat) <- c("lon", "lat", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "tipo")
names(inv) <- c("lon", "lat", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "tipo")

str(nat)
str(inv)


# Combine the data into a single dataset and...
ric_occ <- rbind(nat, inv)

# Save the combined dataset!
write.csv(ric_occ, "ric_occ.csv", row.names = FALSE)



## MERGE ENVIRONMENTAL SPACES ##
# Add a "type" column to each dataset to distinguish between species (we will reuse them for the boxplots later)
clim_nat$tipo <- "nat"
clim_inv$tipo <- "inv"

# Check the column names in nat and inv.
print(names(clim_nat))
print(names(clim_inv))


# Combine the data into a single dataset
clim_ric <- rbind(clim_nat, clim_inv)
write.csv(clim_ric, "clim_ric.csv", row.names = FALSE)


# Remove the extra columns

clim_ric_L <- subset(clim_ric, select = -c(X, x, y, tipo))
clim_nat_L <- subset(clim_nat, select = -c(X, x, y, tipo))
clim_inv_L <- subset(clim_inv, select = -c(X, x, y, tipo))

# Transform climate-related datasets into a matrix

clim_ric_M <- data.matrix(clim_ric_L, rownames.force = 19)
clim_nat_M <- data.matrix(clim_nat_L, rownames.force = 19)
clim_inv_M <- data.matrix(clim_inv_L, rownames.force = 19)

# And with the occs

nat_M <- data.matrix(nat, rownames.force = 22)
inv_M <- data.matrix(inv, rownames.force = 22)


# These analyses were performed using the available guide as a help: https://plantarum.ca/2023/07/28/ecospat-terra/
# PLANTARUM: PCA

pca.clim <- dudi.pca(clim_ric_M, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)
global.scores <- pca.clim$li

nativeRic.scores <- suprow(pca.clim, data.frame(nat_M)[, colnames(clim_ric_M)])$li   
invasiveRic.scores <- suprow(pca.clim, data.frame(inv_M)[, colnames(clim_ric_M)])$li

nativeEnv.scores <- suprow(pca.clim, clim_nat_M)$li
invasiveEnv.scores <- suprow(pca.clim, clim_inv_M)$li


plot(global.scores, pch = 20, asp = 1,
     col = adjustcolor(1, alpha.f = 0.2), cex = 2,
     xlab = "PC1", ylab = "PC2") 
points(nativeRic.scores, pch = 20, col = 2, cex = 2)
points(invasiveRic.scores, pch = 20, col = 4, cex = 2)



#PLANTARUM: Occurence Densities Grid

nativeGrid <- ecospat.grid.clim.dyn(global.scores,
                                    nativeEnv.scores,
                                    nativeRic.scores)

invasiveGrid <- ecospat.grid.clim.dyn(global.scores,
                                      invasiveEnv.scores, 
                                      invasiveRic.scores)

ecospat.plot.niche.dyn(nativeGrid, invasiveGrid,               #Niche overlap
                       quant = 0.05, name.axis1 = "PC1",
                       name.axis2 = "PC2" ) 

ecospat.shift.centroids(invasiveRic.scores, nativeRic.scores, invasiveEnv.scores, nativeEnv.scores,col= "purple")


?ecospat.plot.contrib
pca.clim$co
pca.clim$eig

ecospat.plot.contrib(contrib=pca.clim$co, eigen=pca.clim$eig)

ecospat.niche.overlap(nativeGrid,invasiveGrid ,cor = TRUE)

#$D
#[1] 0.001652451 # Schoener's D index (Overlap ratio)

#$I
#[1] 0.004110775

ecospat.niche.dyn.index(nativeGrid, invasiveGrid, intersection = 0.1)$dynamic.index.w


# PLANTARUM: Potential niche for each region with its own castor oil plant
data(wrld_simpl)

#  For the nat_M dataframe
nat_lon_lat <- nat_M[, c("lon", "lat")]

# For the inv_M dataframe
inv_lon_lat <- inv_M[, c("lon", "lat")]


geoGrid_inv <- expand.grid(longitude =
                         seq(-10, 5, length.out = 500),
                       latitude =
                         seq(36, 44, length.out = 500))

invGeoGrid <- ecospat.grid.clim.dyn(geoGrid_inv, geoGrid_inv, inv_lon_lat)

ecospat.plot.niche.dyn(invGeoGrid, invGeoGrid, quant = 0)
plot(wrld_simpl, add = TRUE)

naMask <- subset(wrld_simpl, NAME %in% c("Spain", "Portugal", "Andorra"))
invGeoGrid <- ecospat.grid.clim.dyn(geoGrid_inv, geoGrid_inv, inv_lon_lat, geomask = naMask)


ecospat.plot.niche.dyn(invGeoGrid, invGeoGrid, quant = 0)


geoGrid_nat <- expand.grid(longitude =
                             seq(33, 51, length.out = 500),
                           latitude =
                             seq(-5, 18, length.out = 500))

natGeoGrid <- ecospat.grid.clim.dyn(geoGrid_nat, geoGrid_nat, nat_lon_lat)

ecospat.plot.niche.dyn(natGeoGrid, natGeoGrid, quant = 0)
plot(wrld_simpl, add = TRUE)

naMask <- subset(wrld_simpl, NAME %in% c("Tanzania", "Djibouti", "Eritrea", "Sudan", "South Sudan", "Ethiopia", "Somalia", "Kenya", "Uganda"))
natGeoGrid <- ecospat.grid.clim.dyn(geoGrid_nat, geoGrid_nat, nat_lon_lat, geomask = naMask)


ecospat.plot.niche.dyn(natGeoGrid, natGeoGrid, quant = 0) 
plot(wrld_simpl, add = TRUE)


## More methods and tests #https://rpubs.com/janikh99/1044769

# Recalculate the niche overlap.
ecospat.niche.overlap(nativeGrid,invasiveGrid ,cor = TRUE)

# Perform niche equivalence testing

eq.test_cons <- ecospat.niche.equivalency.test(nativeGrid, invasiveGrid, rep = 100, 
                                        overlap.alternative = "higher", intersection = NA) #conservatism hypothesis

eq.test_diverg <- ecospat.niche.equivalency.test(nativeGrid, invasiveGrid, rep = 100, 
                                               overlap.alternative = "lower", intersection = NA) #niche divergence



# We perform the niche similarity test

sim.test_overlap.l <- ecospat.niche.similarity.test(nativeGrid, 
                                          invasiveGrid, rep = 100, overlap.alternative = "lower")
sim.test_overlap.h <- ecospat.niche.similarity.test(nativeGrid, 
                                                  invasiveGrid, rep = 100, overlap.alternative = "higher")


sim.test_expansion.h <- ecospat.niche.similarity.test(nativeGrid, 
                                                    invasiveGrid, rep = 100, expansion.alternative = "higher")
sim.test_expansion.l <- ecospat.niche.similarity.test(nativeGrid, 
                                                      invasiveGrid, rep = 100, expansion.alternative = "lower")


sim.test_stability.h <- ecospat.niche.similarity.test(nativeGrid, 
                                                      invasiveGrid, rep = 100, stability.alternative = "higher")
sim.test_stability.l <- ecospat.niche.similarity.test(nativeGrid, 
                                                      invasiveGrid, rep = 100, stability.alternative = "lower")


sim.test_unfilling.h <- ecospat.niche.similarity.test(nativeGrid, 
                                                      invasiveGrid, rep = 100, unfilling.alternative = "higher")
sim.test_unfilling.l <- ecospat.niche.similarity.test(nativeGrid, 
                                                      invasiveGrid, rep = 100, unfilling.alternative = "lower")



# plot the tests

par(mfrow=c(1,2))
ecospat.plot.overlap.test(eq.test_cons, "D", "Equivalency")
ecospat.plot.overlap.test(eq.test_diverg, "D", "Equivalency")


par(mfrow=c(1,2))
ecospat.plot.overlap.test(sim.test_overlap.h, "D", "Similarity")
ecospat.plot.overlap.test(sim.test_overlap.l, "D", "Similarity")


par(mfrow=c(1,2))
ecospat.plot.overlap.test(eq.test_cons, "D", "Equivalency")
ecospat.plot.overlap.test(sim.test_overlap.h, "D", "Similarity")



## Analysis not used

# Niche dynamics versus a gradient #https://rpubs.com/janikh99/1044769

# gridding the native niche

grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = clim_ric_M[,1],
                                         glob1 = data.frame(clim_nat_M[,1]),
                                         data.frame(nat_M)[,3], R = 1000, th.sp = 0)

# gridding the invasive niche

grid.clim.t.inv <- ecospat.grid.clim.dyn (glob = clim_ric_M[,1], 
                                          glob1 = data.frame(clim_inv_M[,1]), 
                                          data.frame(inv_M)[,3], R = 1000, th.sp = 0)

t.dyn <- ecospat.niche.dyn.index (grid.clim.t.nat, grid.clim.t.inv, intersection=0.1)

ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0.1, interest=2, title= "Niche Overlap", name.axis1="Average temperature")

# showing the shift of the niche centroid along the temperature gradient (compared to the shift of the available climate in the study area)

ecospat.shift.centroids(data.frame(nat_M)[,3],
                        data.frame(inv_M)[,3],
                        data.frame(clim_nat_M)[,1],
                        data.frame(clim_inv_M)[,1])





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## ### ### ### ### NICHE MODELIZATION ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


############# VIF analysis #############
library(raster)
library(rgdal)
library(rgeos)
library(usdm)
library(gtools)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(sf)
library(dplyr)

# So that the values are given in scientific notation and to take as characters
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999, stringsAsFactors = FALSE)

occ_iber <- read.csv("ric_iber_max.csv")
occ_nat <- read.csv("ric_nat_max.csv")
worldclim <- getData('worldclim', var='bio', res=2.5)


por <- raster::getData("GADM", country = "PORTUGAL", level = 1)
esp <- raster::getData("GADM", country = "SPAIN", level = 1)

plot(esp)
points(occ_iber$Longitud, occ_iber$Latitud, pch = 16, col = "red", cex = 0.4)

# Extract values from biovariables
occ_iber_swd <- raster::extract(worldclim, occ_iber[,3:2]) %>% cbind(occ_iber, .) %>% as.data.frame()
occ_nat_swd <- raster::extract(worldclim, occ_nat[,3:2]) %>% cbind(occ_nat, .) %>% as.data.frame()

occ_iber_swd <- occ_iber_swd[complete.cases(occ_iber_swd), ]
occ_nat_swd <- occ_nat_swd[complete.cases(occ_nat_swd), ]

# Extract matrix
mtx_iber <- occ_iber_swd[,4:ncol(occ_iber_swd)]
mtx_nat <- occ_nat_swd[,4:ncol(occ_nat_swd)]


# Correlation analysis
m_iber <- cor(mtx_iber)
m_iber
corrplot(m_iber, method = "ellipse")

m_nat <- cor(mtx_nat)
m_nat
corrplot(m_nat, method = "ellipse")


# VIF analysis

# Iberian Peninsula
vif.res_iber <- vif(x = occ_iber_swd[,4:ncol(occ_iber_swd)])
vif.step_iber <- vifstep(x = occ_iber_swd[,4:ncol(occ_iber_swd)], th = 10) # Our threshold will be 10, meaning we will eliminate all variables that have a correlation greater than 10.
vrs_iber <- vif.step_iber@results$Variables %>% as.character()

vif.step_iber
vrs_iber
#"bio4"  "bio7"  "bio8"  "bio14" "bio15" "bio18" are not highly correlated with the other variables

#BIO4 = Temperature Seasonality (standard deviation ×100)
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO18 = Precipitation of Warmest Quarter


# Northeast Africa
vif.res_nat <- vif(x = occ_nat_swd[,4:ncol(occ_nat_swd)])
vif.step_nat <- vifstep(x = occ_nat_swd[,4:ncol(occ_nat_swd)], th = 10) # Our threshold will be 10, meaning we will eliminate all variables that have a correlation greater than 10.

vrs_nat <- vif.step_nat@results$Variables %>% as.character()

vif.step_nat
vrs_nat

# "bio3"  "bio4"  "bio8"  "bio9"  "bio10" "bio12" "bio14"

# The message "essentially perfect fit: summary may be unreliable" appears. 
# This indicates that the fit of the linear regression model for calculating the VIF is almost perfect, 
# suggesting that there is a strong relationship between the predictor variables.

# In the previous VIF analysis (of the Iberian Peninsula), the values did not exceed 3.901351 (bio18), so we will exclude those that exceed this value.


# Northeast Africa
vif.res_nat <- vif(x = occ_nat_swd[,4:ncol(occ_nat_swd)])
vif.step_nat <- vifstep(x = occ_nat_swd[,4:ncol(occ_nat_swd)], th = 3.9) #Our threshold will be 10, meaning we will eliminate all variables that have a correlation greater than 10
vrs_nat <- vif.step_nat@results$Variables %>% as.character()

vif.step_nat
vrs_nat

# "bio3"  "bio4"  "bio8"  "bio12" "bio14"


# BIO3 = Isothermality (BIO2/BIO7) (×100)

# BIO4 = Temperature Seasonality (standard deviation ×100)

# BIO8 = Mean Temperature of Wettest Quarter

# BIO12 = Annual Precipitation

# BIO14 = Precipitation of Driest Month



####################################################################################
# Reissue of vector layers get it from MaxEnt and Worldclim for QGIS

library(raster)
library(rgdal)
library(leaflet)


bio6 <- raster("C:/Your_directory/22/wc2-5/bio6.bil")
bio6_nat_ecos <- as.data.frame(bio6, xy = TRUE)

bio6_nat_ecos_clean <- na.omit(bio6_nat_ecos)
detach(package:raster)
bio6_NoFrost <- subset.data.frame(bio6_nat_ecos_clean, bio6 > 40)
bio6_NoFrost_sf <- st_as_sf(bio6_NoFrost, coords = c("x", "y"), crs = 4326)

bio6_NoFrost_poly <- st_convex_hull(bio6_NoFrost_sf)
st_write(bio6_NoFrost_poly, "bio6_NoFrost_poly.shp")


plot(polys)

st_write(bio6_NoFrost_sf, "bio6_NoFrost_sf.shp")


# Iberian peninsula
iber_max_may  <- raster("C:/Your_directory/MAXENT/maxent/RESULTADOS/iberica_mayo/Ricinus_communis.bil")
nofrost  <- readOGR("C:/Your_directory/QGIS/maxent/bio6_nofrost_def.shp")
plot(nofrost)
plot(iber_max_may)

iber_max_masked <- mask(x = iber_max_may, mask = nofrost)
plot(iber_max_masked)

iber_max_clean <- crop(x = iber_max_masked, y = extent(nofrost))
plot(iber_max_clean)


# Northeast Africa
nat_max_may  <- raster("C:/Your_directory/MAXENT/maxent/RESULTADOS/native_mayo/Ricinus_communis.bil")
nofrost  <- readOGR("C:/Your_directory/QGIS/maxent/bio6_nofrost_def.shp")
plot(nofrost)
plot(nat_max_may)

nat_max_masked <- mask(x = nat_max_may, mask = nofrost)
plot(nat_max_masked)

nat_max_clean <- crop(x = nat_max_masked, y = extent(nofrost))
plot(nat_max_clean)


writeRaster(iber_max_clean, filename = "C:/Your_directory/QGIS/maxent/iber_max_clean.bil", overwrite = TRUE)
writeRaster(nat_max_clean, filename = "C:/Your_directory/QGIS/maxent/nat_max_clean.bil", overwrite = TRUE)






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## ### ### ### ### BOXPLOTS ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# library
library(ggplot2)
library(patchwork)
library(dplyr)

ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot()


boxplot1 <- ggplot(ric_occ, aes(x=tipo, y=bio1, fill=tipo)) + 
  geom_boxplot(notchwidth = 20)

boxplot2 <- ggplot(ric_occ, aes(x=tipo, y=bio2, fill=tipo)) + 
  geom_boxplot()

boxplot3 <- ggplot(ric_occ, aes(x=tipo, y=bio3, fill=tipo)) + 
  geom_boxplot()

boxplot4 <- ggplot(ric_occ, aes(x=tipo, y=bio4, fill=tipo)) + 
  geom_boxplot()

boxplot5 <- ggplot(ric_occ, aes(x=tipo, y=bio5, fill=tipo)) + 
  geom_boxplot()

boxplot6 <- ggplot(ric_occ, aes(x=tipo, y=bio6, fill=tipo)) + 
  geom_boxplot()

boxplot7 <- ggplot(ric_occ, aes(x=tipo, y=bio7, fill=tipo)) + 
  geom_boxplot()

boxplot8 <- ggplot(ric_occ, aes(x=tipo, y=bio8, fill=tipo)) + 
  geom_boxplot()

boxplot9 <- ggplot(ric_occ, aes(x=tipo, y=bio9, fill=tipo)) + 
  geom_boxplot()

boxplot10 <- ggplot(ric_occ, aes(x=tipo, y=bio10, fill=tipo)) + 
  geom_boxplot()

boxplot11 <- ggplot(ric_occ, aes(x=tipo, y=bio11, fill=tipo)) + 
  geom_boxplot()

boxplot12 <- ggplot(ric_occ, aes(x=tipo, y=bio12, fill=tipo)) + 
  geom_boxplot()

boxplot13 <- ggplot(ric_occ, aes(x=tipo, y=bio13, fill=tipo)) + 
  geom_boxplot()

boxplot14 <- ggplot(ric_occ, aes(x=tipo, y=bio14, fill=tipo)) + 
  geom_boxplot()

boxplot15 <- ggplot(ric_occ, aes(x=tipo, y=bio15, fill=tipo)) + 
  geom_boxplot()

boxplot16 <- ggplot(ric_occ, aes(x=tipo, y=bio16, fill=tipo)) + 
  geom_boxplot()

boxplot17 <- ggplot(ric_occ, aes(x=tipo, y=bio17, fill=tipo)) + 
  geom_boxplot()

boxplot18 <- ggplot(ric_occ, aes(x=tipo, y=bio18, fill=tipo)) + 
  geom_boxplot()

boxplot19 <- ggplot(ric_occ, aes(x=tipo, y=bio19, fill=tipo)) + 
  geom_boxplot()

boxplot1 + boxplot2 + boxplot3 + boxplot4 + boxplot5 + boxplot6 + boxplot7 + boxplot8 + 
  boxplot9 + boxplot10 + boxplot11 + boxplot12 + boxplot13 + boxplot14 + boxplot15 + boxplot16 +
  boxplot17 + boxplot18 + boxplot19 



# We selected the most independent variables with the highest percentage contribution from MaxEnt (bio3, 14, 8, 15, and 18).

boxplot3 + boxplot14 + boxplot8 + boxplot15 + boxplot18


# For bio2
lm.bio2 <- lm(bio2~tipo, ric_occ)
anova(lm.bio2)

stats_nat2 <- ric_occ %>% 
  filter(tipo == "nat") %>%
  summarise(
    minimo = min(bio2, na.rm = TRUE),
    maximo = max(bio2, na.rm = TRUE),
    media = mean(bio2, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio2), decreasing = TRUE)[1]))
  )

stats_inv2 <- ric_occ %>% 
  filter(tipo == "inv") %>%
  summarise(
    minimo = min(bio2, na.rm = TRUE),
    maximo = max(bio2, na.rm = TRUE),
    media = mean(bio2, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio2), decreasing = TRUE)[1]))
  )

stats_nat2
stats_inv2



# For bio3
lm.bio3 <- lm(bio3~tipo, ric_occ)
anova(lm.bio3)

stats_nat3 <- ric_occ %>% 
  filter(tipo == "nat") %>%
  summarise(
    minimo = min(bio3, na.rm = TRUE),
    maximo = max(bio3, na.rm = TRUE),
    media = mean(bio3, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio3), decreasing = TRUE)[1]))
  )

stats_inv3 <- ric_occ %>% 
  filter(tipo == "inv") %>%
  summarise(
    minimo = min(bio3, na.rm = TRUE),
    maximo = max(bio3, na.rm = TRUE),
    media = mean(bio3, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio3), decreasing = TRUE)[1]))
  )

stats_nat3
stats_inv3



# For bio8
lm.bio8 <- lm(bio8~tipo, ric_occ)
anova(lm.bio8)

stats_nat8 <- ric_occ %>% 
  filter(tipo == "nat") %>%
  summarise(
    minimo = min(bio8, na.rm = TRUE),
    maximo = max(bio8, na.rm = TRUE),
    media = mean(bio8, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio8), decreasing = TRUE)[1]))
  )

stats_inv8 <- ric_occ %>% 
  filter(tipo == "inv") %>%
  summarise(
    minimo = min(bio8, na.rm = TRUE),
    maximo = max(bio8, na.rm = TRUE),
    media = mean(bio8, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio8), decreasing = TRUE)[1]))
  )

stats_nat8
stats_inv8



# For bio14
lm.bio14 <- lm(bio14~tipo, ric_occ)
anova(lm.bio14)

stats_nat14 <- ric_occ %>% 
  filter(tipo == "nat") %>%
  summarise(
    minimo = min(bio14, na.rm = TRUE),
    maximo = max(bio14, na.rm = TRUE),
    media = mean(bio14, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio14), decreasing = TRUE)[1]))
  )

stats_inv14 <- ric_occ %>% 
  filter(tipo == "inv") %>%
  summarise(
    minimo = min(bio14, na.rm = TRUE),
    maximo = max(bio14, na.rm = TRUE),
    media = mean(bio14, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio14), decreasing = TRUE)[1]))
  )

stats_nat14
stats_inv14



# For bio15
lm.bio15 <- lm(bio15~tipo, ric_occ)
anova(lm.bio15)

stats_nat15 <- ric_occ %>% 
  filter(tipo == "nat") %>%
  summarise(
    minimo = min(bio15, na.rm = TRUE),
    maximo = max(bio15, na.rm = TRUE),
    media = mean(bio15, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio15), decreasing = TRUE)[1]))
  )

stats_inv15 <- ric_occ %>% 
  filter(tipo == "inv") %>%
  summarise(
    minimo = min(bio15, na.rm = TRUE),
    maximo = max(bio15, na.rm = TRUE),
    media = mean(bio15, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio15), decreasing = TRUE)[1]))
  )

stats_nat15
stats_inv15



# For bio18
lm.bio18 <- lm(bio18~tipo, ric_occ)
anova(lm.bio18)

stats_nat18 <- ric_occ %>% 
  filter(tipo == "nat") %>%
  summarise(
    minimo = min(bio18, na.rm = TRUE),
    maximo = max(bio18, na.rm = TRUE),
    media = mean(bio18, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio18), decreasing = TRUE)[1]))
  )

stats_inv18 <- ric_occ %>% 
  filter(tipo == "inv") %>%
  summarise(
    minimo = min(bio18, na.rm = TRUE),
    maximo = max(bio18, na.rm = TRUE),
    media = mean(bio18, na.rm = TRUE),
    moda = as.numeric(names(sort(table(bio18), decreasing = TRUE)[1]))
  )

stats_nat18
stats_inv18


# Bibliography

library(maptools)
library(ggspatial)
library(raster)
library(ade4)
library(ggplot2)
library(rnaturalearthdata)
library(sp)
library(ecospat)
library(usdm)


citation("maptools")      
citation("ggspatial")
citation("raster")
citation("ade4")
citation("ggplot2")
citation("rnaturalearthdata")
citation("sp")
citation("ecospat")
citation("usdm")


