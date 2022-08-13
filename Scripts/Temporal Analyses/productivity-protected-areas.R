## Mean protected area fPAR as a function of time
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
data <- read_csv("Data/master-dataset-time(EBV).csv")

## Prep dataset ----
# Create fpar dataset
fparDf <- data[,c(2,3,4,5,6,7,10,11,12,13,14,15)]
# Remove NAs
fparDf <- na.omit(fparDf)
rm(data)

# Rescale values of fPAR [0,1]
fparDf$PA_FPAR_MEAN <- (fparDf$PA_FPAR_MEAN/10000)

# Remove extreme northern PAs
# Based on assumption that fPAR values in northern PAs are extremely low
fparDf <- fparDf[!(fparDf$WDPA_PID %in% northernPaIds),] 

# Set variables to factor
fparDf$BIOME25 <- factor(fparDf$BIOME25)
fparDf$IUCN_CAT <- factor(fparDf$IUCN_CAT)
fparDf$IUCN_GROUP <- factor(fparDf$IUCN_GROUP)
fparDf$CONTINENT <- factor(fparDf$CONTINENT)
fparDf$WDPA_PID <- factor(fparDf$WDPA_PID)

# Identify duplicate PAs in a single year
# Duplicate coordinate points throws an error when calculating Moran's I
df2003 <- fparDf[(fparDf$YEAR == '2003'),]
coords03 <- SpatialPoints(data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT)))
duplicates03 <- as.data.frame(zerodist(coords03,zero = 0.0))
dup.coords03 <- data.frame(LONG1 = c(coords03$LONG[3567],coords03$LONG[5315],coords03$LONG[129],coords03$LONG[11957],coords03$LONG[11958],coords03$LONG[11959]), 
                           LONG2 = c(coords03$LONG[5885],coords03$LONG[9115],coords03$LONG[11951],coords03$LONG[11961],coords03$LONG[11962],coords03$LONG[11963]), 
                           LAT1 = c(coords03$LAT[3567],coords03$LAT[5315],coords03$LAT[129],coords03$LAT[11957],coords03$LAT[11958],coords03$LAT[11959]), 
                           LAT2 = c(coords03$LAT[5885],coords03$LAT[9115],coords03$LAT[11951],coords03$LAT[11961],coords03$LAT[11962],coords03$LAT[11963]))

# Omit duplicates from fparDf 
fparDf <- fparDf[!(fparDf$LONG %in% dup.coords03$LONG1),]
rm(coords03,duplicates03,dup.coords03,df2003)

# Check if there are the same number of PAs sampled in each year
# Assuming the PAs in first and last year are representative of PAs sampled in all other years
df03 <- fparDf[(fparDf$YEAR == 2003),]
df15 <- fparDf[(fparDf$YEAR == 2015),]
diff <- df15[!(df15$WDPA_PID %in% df03$WDPA_PID),]

# Remove PAs that were not sampled across all years
fparDf <- fparDf[!(fparDf$WDPA_PID %in% diff$WDPA_PID),]
rm(diff,df03,df15)

# Remove PAs with biome type 99 ('rock and ice')
fparDf <- fparDf[!(fparDf$BIOME == '99'),] 

# Found a PA that was sampled twice in each year (WDPA_PID 67117)
fparDf <- fparDf[!(fparDf$WDPA_PID == '67117' & fparDf$GIS_AREA < 0.5),]

# Filter by size limit of raster resolution
fparDf <- fparDf[!(fparDf$GIS_AREA < 1),]

# Scale year (set year 1 to 0)
fparDf$YEAR_s <- fparDf$YEAR - 2003

# Convert values of 0 to 0.000001 to fit the beta family. 
fparDf$PA_FPAR_MEAN[fparDf$PA_FPAR_MEAN == 0] <- 0.000001

## Exploratory analyses ---- 

# Check distribution of response variable
hist(fparDf$PA_FPAR_MEAN)
# These data fall between 0 and 1, suggesting best distribution is the beta family

# Check for outliers
dotchart(fparDf$PA_FPAR_MEAN, xlab = "PA fPAR", ylab = "Order of the data")
# No outliers

# Test fit of gaussian vs. beta distribution
m0 <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP, 
             data = fparDf,
             REML = TRUE)

m1 <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP, 
              family = beta_family(),  
              data = fparDf,
              REML = TRUE)

# Check residuals
sims0 <- simulateResiduals(m0)
plot(sims0)
sims1 <- simulateResiduals(m1)
plot(sims1)

# Compare model fits
AIC(m0,m1)
anova(m0,m1)
# Gaussian model has a non-significantly better fit than the beta model based on AIC
# However, these are proportional data - a beta distribution is more appropriate

## Test random effects structures ----
m <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP, 
              family = beta_family(),  
              data = fparDf,
              REML = TRUE)

m.b <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25), 
             family = beta_family(),  
             data = fparDf,
             REML = TRUE)

m.c <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|CONTINENT), 
             family = beta_family(),  
             data = fparDf,
             REML = TRUE)

m.w <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|WDPA_PID), 
             family = beta_family(),  
             data = fparDf,
             REML = TRUE)

m.bw <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|WDPA_PID), 
             family = beta_family(),  
             data = fparDf,
             REML = TRUE)

m.cw <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|CONTINENT) + (1|WDPA_PID), 
             family = beta_family(),  
             data = fparDf,
             REML = TRUE)

m.bc <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT), 
             family = beta_family(),  
             data = fparDf,
             REML = TRUE)

m.all <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
                 family = beta_family(),  
                 data = fparDf,
                 REML = TRUE)

anova(m,m.all,m.b,m.c,m.w,m.bc,m.bw,m.cw)
anova(m,m.all)
anova(m.bw,m.all)
# Full model (m.all) has significantly best fit

## Test random slopes ----
m.all <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
                 family = beta_family(),  
                 data = fparDf,
                 REML = TRUE)

m.bio <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
                 family = beta_family(),  
                 data = fparDf,
                 REML = TRUE)

m.cont <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
                 family = beta_family(),  
                 data = fparDf,
                 REML = TRUE)

anova(m,m.all,m.bio,m.cont)
anova(m.all,m.cont)
# Full model (m.all) has significantly best fit

## Check best fit model for spatial autocorrelation ----
sims <- simulateResiduals(m.all)
df2003 <- fparDf[(fparDf$YEAR_s == 0),]
coords <- data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT))

groupedSims <- recalculateResiduals(sims, group = fparDf$WDPA_PID)
testSpatialAutocorrelation(groupedSims,coords$LONG,coords$LAT) 
# Significant SAC was detected (observed = 3.2501e-01; p-value < 2.2e-16).

## Test spatial correlation structures ----
# N.B: spatial models were tested on a subset of the data
# Computational demand of fitting spatial models to full dataset was not feasible

# Create a random subset of fparDf
df03 <- fparDf[(fparDf$YEAR == 2003),] # Get data for 2003
s1 <- df03[sample(nrow(df03), 2326),] # Take a random sample
df03 <- df03[!(df03$WDPA_PID %in% s1$WDPA_PID),] # Remove s1 from df03
fparDfSubset <- fparDf[(fparDf$WDPA_PID %in% s1$WDPA_PID),] # Get the full data set that matches the PA in s1

# Add position and group factor to fit spatial correlation parameter
fparDfSubset$pos <- numFactor(fparDfSubset$LONG,fparDfSubset$LAT)
fparDfSubset$group <- factor(1)

# Fit non-spatial model
m <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
               family = beta_family(),  
               data = fparDfSubset,
               REML = TRUE)

# Fit spatial model - gaussian
m.g <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID) + gau(0 + pos|group), 
                 family = beta_family(),  
                 data = fparDfSubset,
                 REML = TRUE)

# Fit spatial model - exponential
m.e <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID) + exp(0 + pos|group), 
               family = beta_family(),  
               data = fparDfSubset,
               REML = TRUE)

# Fit spatial model - matern
m.m <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID) + mat(0 + pos|group), 
               family = beta_family(),  
               data = fparDfSubset,
               REML = TRUE)

# Compare model fits
anova(m,m.e,m.m,m.g)
anova(m.e,m.m)
# exponential model has significantly better fit than non-spatial model 

# Check residuals
sims <- simulateResiduals(m.e)
plot(sims)
hist(resid(m.e)) 

## Fit final model ----
# N.B: The final model was fit to multiple subsets of the full dataset
# The computational demand of fitting spatial models to the full dataset was not feasible
# Each subset was randomly sampled from the full dataset
# The outputs of this model fit to multiple subdatasets were averaged to obtain the final results
finalModel <- glmmTMB(PA_FPAR_MEAN ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID) + exp(0 + pos|group), 
                family = beta_family(),  
                data = fparDfSubset,
                REML = TRUE)
                #control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

# Save model object
saveRDS(finalModel, file = "pa-fpar-final-model.rds")