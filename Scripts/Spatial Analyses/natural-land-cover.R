## Protected area natural land cover as a function of buffer natural land cover
## Exploratory analysis and GLMM fitting
## August 2022

## Workspace ----
# Load libraries
library(tidyverse)
library(MuMIn)
library(DHARMa)
library(sp)
library(glmmTMB)

# Load data
data <- read_csv("Data/master-dataset-netEBVs.csv")
sampledPAs <- read_csv("Data/natural-land-cover-corresponding-pa-IDs.csv")

## Prep dataset ----
# Create fpar dataset
landcoverDf <- data[,c(2,3,4,5,6,7,8,9,10,11,24,25,26)]
# Remove NAs
landcoverDf <- na.omit(landcoverDf)
rm(data)

# Calculate proportion of natural land cover change 
landcoverDf$PA_prop <- (landcoverDf$PA_NAT_AREA/landcoverDf$PA_Area) 
landcoverDf$BUFF_prop <- (landcoverDf$BUFF_NAT_AREA/landcoverDf$BUFF_Area) 

# Set variables to factor
landcoverDf$BIOME <- factor(landcoverDf$BIOME)
landcoverDf$IUCN_GROUP <- factor(landcoverDf$IUCN_GROUP)
landcoverDf$CONTINENT <- factor(landcoverDf$CONTINENT)
landcoverDf$WDPA_PID <- factor(landcoverDf$WDPA_PID)

# Identify duplicate PAs
# Duplicate coordinate points throws an error when calculating Moran's I
coords <- data.frame(LONG = c(landcoverDf$LONG), LAT = c(landcoverDf$LAT))
points <- SpatialPoints(coords)
duplicates <- as.data.frame(zerodist(points,zero = 0.0))
# No duplicates
rm(coords,points,duplicates)

# Filter by size limit of raster resolution
landcoverDf <- landcoverDf[(landcoverDf$GIS_AREA >= 1),]

# Remove protected areas that were not sampled time series analyses
landcoverDf <- landcoverDf[!(landcoverDf$WDPA_PID %in% sampledPAs$LandCover),]
rm(sampledPAs)

## Exploratory analyses ---- 

# Check distribution of response variable
hist(landcoverDf$PA_prop) 
hist(landcoverDf$BUFF_prop)

# Check for outliers
dotchart(landcoverDf$PA_prop, xlab = "PA AREA", ylab = "Order of the data")
dotchart(landcoverDf$BUFF_prop, xlab = "BUFF AREA", ylab = "Order of the data")
# Remove outliers
landcoverDf <- landcoverDf[(landcoverDf$BUFF_prop < 0.5),]
landcoverDf <- landcoverDf[(landcoverDf$BUFF_prop > -0.3),]

# Test fit of gaussian distribution
m <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP, 
                 data = landcoverDf, 
                 REML = TRUE)

# Check residuals
sims <- simulateResiduals(m)
plot(sims)

## Test random effects structures ----
m.all <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|BIOME) + (1|CONTINENT), 
              data = landcoverDf, 
              REML = TRUE,
              control = glmmTMBControl(optimizer=optim,optArgs=list(method="L-BFGS-B")))

m.bio <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|BIOME), 
              data = landcoverDf, 
              REML = TRUE,
              control = glmmTMBControl(optimizer=optim,optArgs=list(method="L-BFGS-B")))

m.cont <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|CONTINENT), 
              data = landcoverDf, 
              REML = TRUE)

anova(m.all,m.bio,m.con,m)
# m.bio has lowest AIC

## Test random slopes ----

m.bio.s <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (BUFF_prop|BIOME), 
                  data = landcoverDf, 
                  REML = TRUE,
                  control = glmmTMBControl(optimizer=optim,optArgs=list(method="L-BFGS-B")))

anova(m.bio,m.bio.s)
# Model with a random slope (m.bio.s) has the best fit. 
# Spatial models would not converge with random slope term
# Further model testing was done on simplified model without random slope (m.bio)

## Check best fit model for spatial autocorrelation ----
sims <- simulateResiduals(m.bio)
testSpatialAutocorrelation(sims, x = landcoverDf$LONG, y = landcoverDf$LAT)
# SAC was found to be significant (0.0226, p-value = 2.2e-16)

## Test spatial correlation structures ----
# N.B: spatial models were tested on a subset of the data
# Computational demand of fitting spatial models to full dataset was not feasible

# Create a random subset of landcoverDf
landcoverDfSubset <- landcoverDf[sample(nrow(landcoverDf), 1163),] # Take a random sample

# Add position and group factor to fit spatial correlation parameter
landcoverDfSubset$pos <- numFactor(landcoverDfSubset$LONG,landcoverDfSubset$LAT)
landcoverDfSubset$group <- factor(1)

# Fit non-spatial model
m <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|BIOME), 
             dispformula = ~0,  
             data = landcoverDfSubset, 
             REML = TRUE)

# Fit spatial model - gaussian
m.g <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|BIOME) + gau(0 + pos|group), 
                    dispformula = ~0,  
                    data = landcoverDfSubset, 
                    REML = TRUE)
                   
# Fit spatial model - exponential
m.e <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|BIOME) + exp(0 + pos|group), 
                    dispformula = ~0,  
                    data = landcoverDfSubset, 
                    REML = TRUE)
                    
# Did not converge
# Fit spatial model - matern
m.m <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|BIOME) + mat(0 + pos|group), 
                    dispformula = ~0,  
                    data = landcoverDfSubset, 
                    REML = TRUE)

# Compare model fits
anova(m,m.g,m.e)
# Exponential model has the best fit

# Check residuals
sims <- simulateResiduals(m.e)
plot(sims)

## Fit final model ----
# N.B: The final model was fit to multiple subsets of the full dataset
# The computational demand of fitting spatial models to the full dataset was not feasible
# Each subset was randomly sampled from the full dataset
# The outputs of this model fit to multiple subdatasets were averaged to obtain the final results
finalModel <- glmmTMB(PA_prop ~ BUFF_prop*IUCN_GROUP + (1|BIOME) + exp(0 + pos|group), 
                      dispformula = ~0,  
                      data = landcoverDfSubset, 
                      REML = TRUE)

# Save model object
saveRDS(finalModel, file = "natural-land-cover-final-model.rds")