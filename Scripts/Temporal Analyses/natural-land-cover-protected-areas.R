## Protected area natural land cover as a function of time
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
data <- read_csv("Data/master-dataset-time(EBV).csv")
sampledPAs <- read_csv("Data/natural-land-cover-corresponding-pa-IDs.csv")

## Prep dataset ----
landcoverDf <- data[,c(2,3,4,5,6,7,8,9,10,11,12,28,29,30)]
# Remove NAs
landcoverDf <- na.omit(landcoverDf) 
rm(data)

# Calculate proportion of natural land cover
landcoverDf$PA_prop <- landcoverDf$PA_NAT_AREA/landcoverDf$PA_Area

# Remove observations where the proportion is > 1 (errors from ArcMap)
landcoverDf <- landcoverDf[(landcoverDf$PA_prop <= 1),]

# Set variables to factor
landcoverDf$BIOME25 <- factor(landcoverDf$BIOME25)
landcoverDf$IUCN_GROUP <- factor(landcoverDf$IUCN_GROUP)
landcoverDf$CONTINENT <- factor(landcoverDf$CONTINENT)
landcoverDf$WDPA_PID <- factor(landcoverDf$WDPA_PID)

# Identify duplicate PAs in a single year
# Duplicate coordinate points thorws an error when calculating Moran's I
df2003 <- landcoverDf[(landcoverDf$YEAR == '2003'),]
coords03 <- SpatialPoints(data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT)))
duplicates03 <- as.data.frame(zerodist(coords03,zero = 0.0))
# No duplicates
rm(df2003,coords03,duplicates03)

# Check if there are the same number of PAs sampled in each year
# Assuming the PAs in first and last year are representative of PAs sampled in all other years
df03 <- landcoverDf[(landcoverDf$YEAR == 2003),]
df15 <- landcoverDf[(landcoverDf$YEAR == 2015),]
diff <- df15[!(df15$WDPA_PID %in% df03$WDPA_PID),]

# Remove PAs that were not sampled across all years
landcoverDf <- landcoverDf[!(landcoverDf$WDPA_PID %in% diff$WDPA_PID),]
rm(diff,df03,df15)

# Repeat for excess PAs in 2003
df03 <- landcoverDf[(landcoverDf$YEAR == 2003),]
df15 <- landcoverDf[(landcoverDf$YEAR == 2015),]
diff <- df03[!(df03$WDPA_PID %in% df15$WDPA_PID),]

landcoverDf <- landcoverDf[!(landcoverDf$WDPA_PID %in% diff$WDPA_PID),]
rm(diff,df03,df15)

df03 <- landcoverDf[(landcoverDf$YEAR == 2003),]
df10 <- landcoverDf[(landcoverDf$YEAR == 2010),]
diff <- df03[!(df03$WDPA_PID %in% df10$WDPA_PID),]

landcoverDf <- landcoverDf[!(landcoverDf$WDPA_PID %in% diff$WDPA_PID),]
rm(diff,df03,df10)

# Filter by size limit of raster resolution
landcoverDf <- landcoverDf[(landcoverDf$GIS_AREA >= 1),]

# Remove protected areas that were not sampled in PA ~ Buffer analyses
# Remove protected areas that were not sampled in PA ~ Buffer analyses
landcoverDf <- landcoverDf[!(landcoverDf$WDPA_PID %in% sampledPAs$LandCover),]
rm(sampledPAs)

# Set values of 1 to just below 1 (to fit beta distribution)
landcoverDf$PA_prop[landcoverDf$PA_prop == 1] <- 0.9999999

# Scale year (set year 1 to 0)
landcoverDf$YEAR_s <- landcoverDf$YEAR - 2003

## Exploratory analyses ---- 

# Check distribution of response variable
hist(landcoverDf$PA_prop)
# These data fall between 0 and 1, suggesting best distribution is the beta family

# Check for outliers
dotchart(landcoverDf$PA_prop, xlab = "PA natural land cover", ylab = "Order of the data")
# Two outliers detected
# Remove outliers
landcoverDf <- landcoverDf[!(landcoverDf$WDPA_PID == 10070),]
landcoverDf <- landcoverDf[!(landcoverDf$WDPA_PID == 33617),]

# Test fit of gaussian vs. poisson and beta distributions
m <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP, 
              data = landcoverDf,
              REML = TRUE)

mp <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP, 
              family = poisson(),  
              data = landcoverDf,
              REML = TRUE)

mb <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP, 
              family = beta_family(),  
              data = landcoverDf,
              REML = TRUE)

# Compare model fits
AIC(m,mp,mb)
anova(m,mp,mb)
# Beta distribution has best fit

# Check residuals
sims <- simulateResiduals(mb)
plot(sims)

## Test random effects structures ----
m <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP, 
              data = landcoverDf,
              REML = TRUE)

m.all <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID),
              family = beta_family(),
              data = landcoverDf,
              REML = TRUE)

m.bw <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|WDPA_PID),
                  family = beta_family(),
                  data = landcoverDf,
                  REML = TRUE)

m.cw <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|CONTINENT) + (1|WDPA_PID),
                   family = beta_family(),
                   data = landcoverDf,
                   REML = TRUE)

m.w <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|WDPA_PID),
                   family = beta_family(),
                   data = landcoverDf,
                   REML = TRUE)

m.b <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|BIOME25),
                  family = beta_family(),
                  data = landcoverDf,
                  REML = TRUE)

m.c <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|CONTINENT),
                   family = beta_family(),
                   data = landcoverDf,
                   REML = TRUE)

m.bc <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT),
                 family = beta_family(),
                 data = landcoverDf,
                 REML = TRUE)

anova(m,m.all,m.b,m.c,m.w,m.bc,m.bw,m.cw)
# Full model has the best fit

## Test random slopes ----
m <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID),
             family = beta_family(),
             data = landcoverDf,
             REML = TRUE)

m.c <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (YEAR_s|CONTINENT) + (1|WDPA_PID),
                   family = beta_family(),
                   data = landcoverDf,
                   REML = TRUE)

m.b <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (YEAR_s|BIOME25) + (1|CONTINENT) + (1|WDPA_PID),
                  family = beta_family(),
                  data = landcoverDf,
                  REML = TRUE)

m.all <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (YEAR_s|BIOME25) + (YEAR_s|CONTINENT) + (1|WDPA_PID),
                   family = beta_family(),
                   data = landcoverDf,
                   REML = TRUE)

anova(m,m.b,m.c,m.all)
anova(m,m.all)
# Full model (m.all) has significantly best fit

## Check best fit model for spatial autocorrelation ----
df2003 <- landcoverDf[(landcoverDf$YEAR == 2003),]
coords <- data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT))

sims <- simulateResiduals(m0.both)
groupedSims <- recalculateResiduals(sims, group = landcoverDf$WDPA_PID)
testSpatialAutocorrelation(groupedSims,coords$LONG,coords$LAT) 
# Significant SAC detected (observed = 0.04028754, expected = -0.00016327, sd = 0.00120646, p-value < 2.2e-16)

## Test spatial correlation structures ----
# N.B: spatial models were tested on a subset of the data
# Computational demand of fitting spatial models to full dataset was not feasible

# Create a random subset of landcoverDf
df03 <- landcoverDf[(landcoverDf$YEAR == 2003),] # Get data for 2003
s1 <- df03[sample(nrow(df03), 1163),] # Take a random sample
df03 <- df03[!(df03$WDPA_PID %in% s1$WDPA_PID),] # Remove s1 from df03
landcoverDfSubset <- landcoverDf[(landcoverDf$WDPA_PID %in% s1$WDPA_PID),] # Get the full data set that matches the PA in s1

# Add position and group factor to fit spatial correlation parameter
landcoverDfSubset$pos <- numFactor(landcoverDfSubset$LONG,landcoverDfSubset$LAT)
landcoverDfSubset$group <- factor(1)

# Fit non-spatial model
m <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (YEAR_s|BIOME25) + (YEAR_s|CONTINENT) + (1|WDPA_PID),
              family = beta_family(),
              data = landcoverDfSubset,
              REML = TRUE)

# Fit spatial model - gaussian
# Did not converge
m.gaus<- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (YEAR_s|BIOME25) + (YEAR_s|CONTINENT) + (1|WDPA_PID) + gau(0 + pos|group),
                 family = beta_family(),
                 data = landcoverDfSubset,
                 REML = TRUE)

# Fit spatial model - exponential
m.exp <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (YEAR_s|BIOME25) + (YEAR_s|CONTINENT) + (1|WDPA_PID) + exp(0 + pos|group),
                 family = beta_family(),
                 data = landcoverDfSubset,
                 REML = TRUE)

# Fit spatial model - matern
# Did not converge
m.mat <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (YEAR_s|BIOME25) + (YEAR_s|CONTINENT) + (1|WDPA_PID) + mat(0 + pos|group),
                 family = beta_family(),
                 data = landcoverDfSubset,
                 REML = TRUE)

# Compare model fits
anova(m, m.exp)
AIC(m, m.exp)
# Non-spatial model has a slightly lower AIC value, but is not a significantly better fit than the spatial model
# Because spatial autocorrelation was detected, use spatial exponential model

# Check fit of exponential model
sims <- simulateResiduals(m.exp)
plot(sims) 

## Fit final model ----
# N.B: The final model was fit to multiple subsets of the full dataset
# The computational demand of fitting spatial models to the full dataset was not feasible
# Each subset was randomly sampled from the full dataset
# The outputs of this model fit to multiple subdatasets were averaged to obtain the final results
finalModel <- glmmTMB(PA_prop ~ YEAR_s*IUCN_GROUP + (YEAR_s|BIOME25) + (YEAR_s|CONTINENT) + (1|WDPA_PID) + exp(0 + pos|group),
                      family = beta_family(),
                      data = landcoverDfSubset,
                      REML = TRUE)

# Save model object
saveRDS(finalModel, file = "protected-areas-land-cover-final-model.rds")
