## Mean protected area human footprint as a function of mean buffer human footprint
## Exploratory analysis and GLMM fitting
## August 2022

## Workspace ----
# Load libraries
library(tidyverse)
library(MuMIn)
library(DHARMa)
library(sp)
library(glmmTMB)

# Parameters
northernPaIds <- c('650', '4768', '100676_A', '61771', '1689', '209739', '200523')

# Load data
data <- read_csv("Data/master-dataset-netEBVs.csv")

## Prep dataset ----
# Create humanFootptint dataset
humanFootprintDf <- data[,c(2,3,4,5,6,9,10,11,21,22,23)]
# Remove NAs
humanFootprintDf <- na.omit(humanFootprintDf)
rm(data)

# Remove extreme northern PAs
# Based on assumption that human populations are lower in northern regions
humanFootprintDf <- humanFootprintDf[!(humanFootprintDf$WDPA_PID %in% northernPaIds),] 

# Set variables to factor
humanFootprintDf$BIOME25 <- factor(humanFootprintDf$BIOME25)
humanFootprintDf$IUCN_GROUP <- factor(humanFootprintDf$IUCN_GROUP)
humanFootprintDf$CONTINENT <- factor(humanFootprintDf$CONTINENT)
humanFootprintDf$WDPA_PID <- factor(humanFootprintDf$WDPA_PID)

# Identify duplicate PAs
# Duplicate coordinate points throws an error when calculating Moran's I
coords <- data.frame(LONG = c(humanFootprintDf$LONG), LAT = c(humanFootprintDf$LAT))
points <- SpatialPoints(coords)
duplicates <- as.data.frame(zerodist(points,zero = 0.0))
dup.coords <- data.frame(LONG1 = c(coords$LONG[401],coords$LONG[3196],coords$LONG[4939],coords$LONG[133],coords$LONG[11734],coords$LONG[11736],coords$LONG[11738]), 
                         LONG2 = c(coords$LONG[1447],coords$LONG[5502],coords$LONG[8640],coords$LONG[11735],coords$LONG[11740],coords$LONG[11741],coords$LONG[11742]), 
                         LAT1 = c(coords$LAT[401],coords$LAT[3196],coords$LAT[4939],coords$LAT[133],coords$LAT[11734],coords$LAT[11736],coords$LAT[11738]), 
                         LAT2 = c(coords$LAT[1447],coords$LAT[5502],coords$LAT[8640],coords$LAT[11735],coords$LAT[11740],coords$LAT[11741],coords$LAT[11742]))

# Omit duplicates
humanFootprintDf <- humanFootprintDf[!(humanFootprintDf$LONG %in% dup.coords$LONG1),]
rm(points,duplicates,dup.coords,coords)

# Filter by size limit of raster resolution
humanFootprintDf <- humanFootprintDf[!(humanFootprintDf$GIS_AREA < 1),]

## Exploratory analyses ----

# Check distribution of response variable
hist(humanFootprintDf$PA_HFP_MEAN) 
hist(humanFootprintDf$BUFF_HFP_MEAN) 

# Check for outliers
dotchart(humanFootprintDf$PA_HFP_MEAN, xlab = "PA Human Footprint", ylab = "Order of the data")
dotchart(humanFootprintDf$BUFF_HFP_MEAN, xlab = "Buffer Human Footprint", ylab = "Order of the data")
# No outliers

# Test fit of gaussian distribution
m <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP,
             data = humanFootprintDf,
             REML = TRUE)

# Check residuals
sims <- simulateResiduals(m)
plot(sims)

## Test random effects structures ----
m.all <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT),
                 data = humanFootprintDf,
                 REML = TRUE)

m.bio <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME),
             data = humanFootprintDf,
             REML = TRUE)

m.con <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|CONTINENT),
             data = humanFootprintDf,
             REML = TRUE)

model.sel(m,m.bio,m.con,m.all)
anova(m,m.bio,m.con,m.all)
# Full model (m.all) has best fit

## Test random slopes ----
m <- m.all

m.both <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1+BUFF_HFP_MEAN|BIOME) + (1+BUFF_HFP_MEAN|CONTINENT),
             data = humanFootprintDf,
             REML = TRUE)

m.bio <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1+BUFF_HFP_MEAN|BIOME) + (1|CONTINENT),
             data = humanFootprintDf,
             REML = TRUE)

m.con <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME) + (1+BUFF_HFP_MEAN|CONTINENT),
             data = humanFootprintDf,
             REML = TRUE)

model.sel(m,m.bio,m.con,m.both)
anova(m,m.bio,m.con,m.both)
# Full model again (m.both) has best fit
# Spatial model would not converge with random slopes
# Further model testing was done without random slopes (i.e. m)

## Check best fit model for spatial autocorrelation ----
sims <- simulateResiduals(m)
testSpatialAutocorrelation(sims,humanFootprintDf$LONG,humanFootprintDf$LAT) 
# Significant SAC is detected in residuals (0.06672287). 

## Test spatial correlation structures ----
# N.B: spatial models were tested on a subset of the data
# Computational demand of fitting spatial models to full dataset was not feasible

# Create a random subset of humanFootprintDf
humanFootprintDfSubset <- humanFootprintDf[sample(nrow(humanFootprintDf), 3970),] # Take a random sample

# Add position and group factor to fit spatial correlation parameter
humanFootprintDfSubset$pos <- numFactor(humanFootprintDfSubset$LONG,humanFootprintDfSubset$LAT)
humanFootprintDfSubset$group <- factor(1)

# Fit non-spatial model
m <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT),
             dispformula = ~0,  
             data = humanFootprintDfSubset,
             REML = TRUE)

# Fit spatial model - gaussian
# Did not converge
m.gaus <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT) + gau(0 + pos|group),
                  dispformula = ~0,  
                  data = humanFootprintDfSubset,
                  REML = TRUE)

# Fit spatial model - exponential
m.exp <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT) + exp(0 + pos|group),
                 dispformula = ~0,  
                 data = humanFootprintDfSubset,
                 REML = TRUE)

# Fit spatial model - matern
# Did not converge
m.mat <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT) + mat(0 + pos|group),
                 dispformula = ~0,  
                 data = humanFootprintDfSubset,
                 REML = TRUE)

# Compare model fits
anova(m,m.exp)
# exponential model has best fit

# Check residuals
sims <- simulateResiduals(m.exp)
plot(m.exp)
hist(resid(m.exp)) 

## Fit final model ----
# N.B: The final model was fit to multiple subsets of the full dataset
# The computational demand of fitting spatial models to the full dataset was not feasible
# Each subset was randomly sampled from the full dataset
# The outputs of this model fit to multiple subdatasets were averaged to obtain the final results
finalModel <- glmmTMB(PA_HFP_MEAN ~ BUFF_HFP_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT) + exp(0 + pos|group),
                      dispformula = ~0,  
                      data = humanFootprintDfSubset,
                      REML = TRUE)

# Save model object
saveRDS(finalModel, file = "human-footprint-final-model.rds")