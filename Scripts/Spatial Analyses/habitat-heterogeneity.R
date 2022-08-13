## Protected area habitat heterogeneity as a function of buffer habitat heterogeneity
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
sampledPAs <- read_csv("Data/habitat-heterogeneity-corresponding-pa-IDs.csv")

## Prep dataset ----
# Create habitat heterogeneity dataset
# Create HET dataset, remove NA values
heterogeneityDf <- data[,c(2,3,4,5,6,9,10,11,27,28,29)]
# Remove NAs
heterogeneityDf <- na.omit(heterogeneityDf)
rm(data)

# Set variables to factor
heterogeneityDf$BIOME <- factor(heterogeneityDf$BIOME)
heterogeneityDf$WDPA_PID <- factor(heterogeneityDf$WDPA_PID)
heterogeneityDf$IUCN_GROUP <- factor(heterogeneityDf$IUCN_GROUP)
heterogeneityDf$CONTINENT <- factor(heterogeneityDf$CONTINENT)

# Identify duplicate PAs in a single year
# Duplicate coordinate points thorws an error when calculating Moran's I
coords <- data.frame(LONG = c(heterogeneityDf$LONG), LAT = c(heterogeneityDf$LAT))
points <- SpatialPoints(coords)
duplicates <- as.data.frame(zerodist(points,zero = 0.0))
dup.coords <- data.frame(LONG1 = c(coords$LONG[429],coords$LONG[4617],coords$LONG[6590],coords$LONG[138],coords$LONG[14428],coords$LONG[14429],coords$LONG[14431]), 
                         LONG2 = c(coords$LONG[1790],coords$LONG[7224],coords$LONG[11012],coords$LONG[14422],coords$LONG[14433],coords$LONG[14434],coords$LONG[14435]), 
                         LAT1 = c(coords$LAT[429],coords$LAT[4617],coords$LAT[6590],coords$LAT[138],coords$LAT[14428],coords$LAT[14429],coords$LAT[14431]), 
                         LAT2 = c(coords$LAT[1790],coords$LAT[7224],coords$LAT[11012],coords$LAT[14422],coords$LAT[14433],coords$LAT[14434],coords$LAT[14435]))

# Omit duplicates 
heterogeneityDf <- heterogeneityDf[!(heterogeneityDf$LONG %in% dup.coords$LONG1),]
rm(points,duplicates,dup.coords,coords)

# Filter by size limit of raster resolution
heterogeneityDf <- heterogeneityDf[(heterogeneityDf$GIS_AREA >= 1),]

# Remove protected areas that were not sampled in temporal analyses
heterogeneityDf <- heterogeneityDf[!(heterogeneityDf$WDPA_PID %in% sampledPAs$Heterogeneity),]
rm(sampledPAs)

# Create a binary data set
# 0 <- no land cover heterogeneity change
# 1 <- land cover heterogeneity change
heterogeneityDf$PA_binary <- ifelse(heterogeneityDf$PA_HET != 0, 1, 0)
heterogeneityDf$BUFF_binary <- ifelse(heterogeneityDf$BUFF_HET != 0, 1, 0)

## Test random effects structures ----
m <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP, 
             family = binomial(),
             data = heterogeneityDf, 
             REML = TRUE)
              
m.all <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1|CONTINENT),
              family = binomial(),
              data = heterogeneityDf, 
              REML = TRUE)

m.bio <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME), 
              family = binomial(),
              data = heterogeneityDf, 
              REML = TRUE)

m.cont <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|CONTINENT), 
              family = binomial(),
              data = heterogeneityDf, 
              REML = TRUE)

anova(m,m.all,m.bio,m.cont)
# Model with both random effects has best fit (significant anova test, lowest AIC)

## Test random slopes ----
m <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1|CONTINENT),           
            family = binomial(),
            data = heterogeneityDf, 
            REML = TRUE)

m.bio <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1+BUFF_binary|BIOME) + (1|CONTINENT),           
             family = binomial(),
             data = heterogeneityDf, 
             REML = TRUE)


m.con <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1+BUFF_binary|CONTINENT),           
                 family = binomial(),
                 data = heterogeneityDf, 
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

m.both <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1+BUFF_binary|BIOME) + (1+BUFF_binary|CONTINENT),           
                 family = binomial(),
                 data = heterogeneityDf, 
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

anova(m,m.bio,m.con,m.both)
# m.con has best fit

## Check best fit model for spatial autocorrelation ----
sims <- simulateResiduals(m.con)
testSpatialAutocorrelation(sims,heterogeneityDf$LONG,heterogeneityDf$LAT)
# Significant SAC is detected (observed = 9.5051e-03, expected = -8.3886e-05, sd = 9.4698e-04, p-value < 2.2e-16)

## Test spatial correlation structures ----
# N.B: spatial models were tested on a subset of the data
# Computational demand of fitting spatial models to full dataset was not feasible

# Create a random subset of heterogeneityDf
heterogeneityDfSubset <- heterogeneityDf[sample(nrow(heterogeneityDf), 1163),] # Take a random sample

# Add position and group factor to fit spatial correlation parameter
heterogeneityDfSubset$pos <- numFactor(heterogeneityDfSubset$LONG,heterogeneityDfSubset$LAT)
heterogeneityDfSubset$group <- factor(1)

# Fit non-spatial model
m <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1+BUFF_binary|CONTINENT),           
                 dispformula = ~ 0,
                 family = binomial(),
                 data = heterogeneityDfSubset, 
                 REML = TRUE)

# Fit spatial model - gaussian
# Did not converge
m.g <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1+BUFF_binary|CONTINENT) + gau(0 + pos|group),
                   dispformula = ~ 0,
                   family = binomial(),
                   data = heterogeneityDfSubset, 
                   REML = TRUE)

# Fit spatial model - exponential
m.e <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1+BUFF_binary|CONTINENT) + exp(0 + pos|group),
                  dispformula = ~ 0,
                  family = binomial(),
                  data = heterogeneityDfSubset, 
                  REML = TRUE)

# Fit spatial model - matern
# Did not converge. 
m.m <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1+BUFF_binary|CONTINENT) + + mat(0 + pos|group),
                  dispformula = ~ 0,
                  family = binomial(),
                  data = heterogeneityDfSubset, 
                  REML = TRUE)

# Compare model fits
anova(m,m.e)
# Exponential model produced the best fit

# Check residuals
sims <- simulateResiduals(m.e)
plot(sims)

## Fit final model ----
# N.B: The final model was fit to multiple subsets of the full dataset
# The computational demand of fitting spatial models to the full dataset was not feasible
# Each subset was randomly sampled from the full dataset
# The outputs of this model fit to multiple subdatasets were averaged to obtain the final results
finalModel <- glmmTMB(PA_binary ~ BUFF_binary*IUCN_GROUP + (1|BIOME) + (1+BUFF_binary|CONTINENT) + exp(0 + pos|group),
                      dispformula = ~ 0,
                      family = binomial(),
                      data = heterogeneityDfSubset, 
                      REML = TRUE)

# Save model object
saveRDS(finalModel, file = "habitat-heterogeneity-final-model.rds")
