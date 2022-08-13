## Mean protected area fPAR as a function of mean buffer fPAR
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
# Create fpar dataset
fparDf <- data[,c(2,3,4,5,6,9,10,11,12,13,14)]
# Remove NAs
fparDf <- na.omit(fparDf)
rm(data)

# Rescale values of fPAR [-1 and 1]
fparDf$PA_FPAR_MEAN <- (fparDf$PA_FPAR_MEAN/10000)
fparDf$BUFF_FPAR_MEAN <- (fparDf$BUFF_FPAR_MEAN/10000)

# Remove extreme northern PAs
# Based on assumption that fPAR values in northern PAs are extremely low
fparDf <- fparDf[!(fparDf$WDPA_PID %in% northernPaIds),] 

# Set variables to factor
fparDf$BIOME25 <- factor(fparDf$BIOME25)
fparDf$IUCN_GROUP <- factor(fparDf$IUCN_GROUP)
fparDf$CONTINENT <- factor(fparDf$CONTINENT)
fparDf$WDPA_PID <- factor(fparDf$WDPA_PID)

# Identify duplicate PAs
# Duplicate coordinate points throws an error when calculating Moran's I
coords <- data.frame(LONG = c(fparDf$LONG), LAT = c(fparDf$LAT))
points <- SpatialPoints(coords)
duplicates <- as.data.frame(zerodist(points,zero = 0.0))
dup.coords <- data.frame(LONG1 = c(coords$LONG[3567],coords$LONG[5315],coords$LONG[129],coords$LONG[11957],coords$LONG[11958],coords$LONG[11959]), 
                         LONG2 = c(coords$LONG[5885],coords$LONG[9115],coords$LONG[11951],coords$LONG[11961],coords$LONG[11962],coords$LONG[11963]), 
                         LAT1 = c(coords$LAT[3567],coords$LAT[5315],coords$LAT[129],coords$LAT[11957],coords$LAT[11958],coords$LAT[11959]), 
                         LAT2 = c(coords$LAT[5885],coords$LAT[9115],coords$LAT[11951],coords$LAT[11961],coords$LAT[11962],coords$LAT[11963]))

# Omit duplicates from fparDf 
fparDf <- fparDf[!(fparDf$LONG %in% dup.coords$LONG1),]
rm(points,duplicates,dup.coords,coords)

# Remove PAs with biome type 99 ('rock and ice')
fparDf <- fparDf[!(fparDf$BIOME == '99'),] 

# Found a PA that was sampled twice (WDPA_PID 67117)
fparDf <- fparDf[!(fparDf$WDPA_PID == '67117' & fparDf$GIS_AREA < 0.5),]

# Filter by size limit of raster resolution
fparDf <- fparDf[!(fparDf$GIS_AREA < 1),]

## Exploratory analyses ----

# Check distribution of response variable
hist(fparDf$PA_FPAR_MEAN) 

# Check for outliers
dotchart(fparDf$PA_FPAR_MEAN, xlab = "PA fPAR", ylab = "Order of the data")
dotchart(fparDf$BUFF_FPAR_MEAN, xlab = "Buffer fPAR", ylab = "Order of the data")
# No outliers

# Test fit of gaussian distribution
m <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP, 
              data = fparDf, 
              REML = TRUE)

# Check residuals
sims <- simulateResiduals(m)
plot(sims)

## Test random effects structures ----
m.all <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT), 
             data = fparDf, 
             REML = TRUE)

m.bio <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1|BIOME), 
                 data = fparDf, 
                 REML = TRUE)

m.con <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1|CONTINENT), 
                 data = fparDf, 
                 REML = TRUE)

m <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP, 
                 data = fparDf, 
                 REML = TRUE)

model.sel(m.all,m.bio,m.con,m)
anova(m.all,m.bio,m.con,m)
# m.all has the best fit model.

## Test random slopes ----
m <- m.all

m.all <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1+BUFF_FPAR_MEAN|BIOME) + (1+BUFF_FPAR_MEAN|CONTINENT), 
                 data = fparDf, 
                 REML = TRUE)

m.bio <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1+BUFF_FPAR_MEAN|BIOME) + (1|CONTINENT), 
                 data = fparDf, 
                 REML = TRUE)

m.con <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1|BIOME) + (1+BUFF_FPAR_MEAN|CONTINENT), 
                 data = fparDf, 
                 REML = TRUE)

model.sel(m,m.bio,m.con,m.all)
anova(m,m.bio,m.con,m.all)
# Full model has best fit (m.all)

## Check best fit model for spatial autocorrelation ----
sims <- simulateResiduals(m.all)
testSpatialAutocorrelation(sims,fparDf$LONG,fparDf$LAT) 
# Significant SAC detected (6.7846e-02; p-value < 2.2e-16)

## Test spatial correlation structures ----
# N.B: spatial models were tested on a subset of the data
# Computational demand of fitting spatial models to full dataset was not feasible

# Create a random subset of fparDf
fparDfSubset <- fparDf[sample(nrow(fparDf), 2326),] # Take a random sample

# Add position and group factor to fit spatial correlation parameter
fparDfSubset$pos <- numFactor(fparDfSubset$LONG,fparDfSubset$LAT)
fparDfSubset$group <- factor(1)

# Fit non-spatial model
m <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1|BIOME) + (1|CONTINENT),
             dispformula = ~0,
             data = fparDfSubset, 
             REML = TRUE)

# Fit spatial model - gaussian
m.gaus <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1+BUFF_FPAR_MEAN|BIOME) + (1+BUFF_FPAR_MEAN|CONTINENT) + gau(0 + pos|group), 
                  dispformula = ~0,  
                  data = fparDfSubset, 
                  REML = TRUE)

# Fit spatial model - exponential
m.exp <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1+BUFF_FPAR_MEAN|BIOME) + (1+BUFF_FPAR_MEAN|CONTINENT) + exp(0 + pos|group), 
                 dispformula = ~0,  
                 data = fparDfSubset, 
                 REML = TRUE)

# Fit spatial model - matern
# Did not converge
m.mat <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1+BUFF_FPAR_MEAN|BIOME) + (1+BUFF_FPAR_MEAN|CONTINENT) + mat(0 + pos|group), 
                 dispformula = ~0,  
                 data = fparDfSubset, 
                 REML = TRUE)

# Compare model fits
anova(m,m.gaus,m.exp)
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
finalModel <- glmmTMB(PA_FPAR_MEAN ~ BUFF_FPAR_MEAN*IUCN_GROUP + (1+BUFF_FPAR_MEAN|BIOME) + (1+BUFF_FPAR_MEAN|CONTINENT) + exp(0 + pos|group), 
                      dispformula = ~0,  
                      data = fparDfSubset, 
                      REML = TRUE)

# Save model object
saveRDS(finalModel, file = "fpar-final-model.rds")