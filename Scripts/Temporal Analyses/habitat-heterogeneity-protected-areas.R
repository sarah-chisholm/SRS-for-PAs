## Protected area habitat heterogeneity as a function of time
## Exploratory analysis and GLMM fitting
## August 2022

## Workspace ----
# Load libraries
library(tidyverse)
library(MuMIn)
library(DHARMa)
library(sp)
library(spaMM)

# Load data
data <- read_csv("Data/master-dataset-time(EBV).csv")
sampledPAs <- read_csv("Data/habitat-heterogeneity-corresponding-pa-IDs.csv")

## Prep dataset ----
# Create habitat heterogeneity dataset
heterogeneityDf <- data[,c(2,3,4,5,6,7,10,11,12,22,23,24)]
# Remove NAs
heterogeneityDf <- na.omit(heterogeneityDf)
rm(data)

# Set variables to factor
heterogeneityDf$BIOME25 <- factor(heterogeneityDf$BIOME25)
heterogeneityDf$WDPA_PID <- factor(heterogeneityDf$WDPA_PID)
heterogeneityDf$IUCN_GROUP <- factor(heterogeneityDf$IUCN_GROUP)
heterogeneityDf$CONTINENT <- factor(heterogeneityDf$CONTINENT)

# Identify duplicate PAs in a single year
# Duplicate coordinate points thorws an error when calculating Moran's I
df2003 <- heterogeneityDf[(heterogeneityDf$YEAR == '2003'),]
coords03 <- SpatialPoints(data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT)))
duplicates03 <- as.data.frame(zerodist(coords03,zero = 0.0))
dup.coords03 <- data.frame(LONG1 = c(coords03$LONG[430],coords03$LONG[4631],coords03$LONG[6609],coords03$LONG[138],coords03$LONG[14478],coords03$LONG[14479],coords03$LONG[14481]), 
                           LONG2 = c(coords03$LONG[1793],coords03$LONG[7246],coords03$LONG[11040],coords03$LONG[14472],coords03$LONG[14483],coords03$LONG[14484],coords03$LONG[14485]), 
                           LAT1 = c(coords03$LAT[430],coords03$LAT[4631],coords03$LAT[6609],coords03$LAT[138],coords03$LAT[14478],coords03$LAT[14479],coords03$LAT[14481]), 
                           LAT2 = c(coords03$LAT[1793],coords03$LAT[7246],coords03$LAT[11040],coords03$LAT[14472],coords03$LAT[14483],coords03$LAT[14484],coords03$LAT[14485]))

# Omit duplicates  
heterogeneityDf <- heterogeneityDf[!(heterogeneityDf$LONG %in% dup.coords03$LONG1),]
rm(coords03,duplicates03,dup.coords03,df2003)

# Check if there are the same number of PAs sampled in each year
# Assuming the PAs in first and last year are representative of PAs sampled in all other years
df03 <- heterogeneityDf[(heterogeneityDf$YEAR == 2003),]
df15 <- heterogeneityDf[(heterogeneityDf$YEAR == 2015),]
diff <- df15[!(df15$WDPA_PID %in% df03$WDPA_PID),]

# Remove PAs that were not sampled across all years
heterogeneityDf <- heterogeneityDf[!(heterogeneityDf$WDPA_PID %in% diff$WDPA_PID),]
rm(diff,df03,df15)

# Repeat to remove PAs that were sampled in 2003 but not 2015
df03 <- heterogeneityDf[(heterogeneityDf$YEAR == 2003),]
df15 <- heterogeneityDf[(heterogeneityDf$YEAR == 2015),]
diff <- df03[!(df03$WDPA_PID %in% df15$WDPA_PID),]
heterogeneityDf <- heterogeneityDf[!(heterogeneityDf$WDPA_PID %in% diff$WDPA_PID),]
rm(diff,df03,df15)

# Filter by size limit of raster resolution
heterogeneityDf <- heterogeneityDf[(heterogeneityDf$GIS_AREA >= 1),]

# Remove protected areas that were not sampled in PA ~ Buffer analyses
heterogeneityDf <- heterogeneityDf[!(heterogeneityDf$WDPA_PID %in% sampledPAs$Heterogeneity),]
rm(sampledPAs)

# Scale year (set year 1 to 0)
heterogeneityDf$YEAR_s <- heterogeneityDf$YEAR - 2003

## Exploratory analyses ---- 
# Check distribution of response variable
hist(heterogeneityDf$PA_HET)

# Check for outliers
dotchart(heterogeneityDf$PA_HET, xlab = "PA Heterogeneity", ylab = "Order of the data")
# No outliers

# Test fit of gaussian vs. poisson distribution
# N.B: glmmTMB functions would not converge, substitute with spaMM::fitme()
m <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP, 
             data = heterogeneityDf, 
             family = "gaussian",
             method = "REML")

mp <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP, 
            data = heterogeneityDf,
            family = "poisson",
            method = "REML")

mcp <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP, 
             data = heterogeneityDf,
             family = "COMPoisson",
             method = "REML")

# Compare model fits
AIC(m, mp, mcp) 
# Compois has best fit

# Check residuals
sims.p <- simulateResiduals(mp)
plot(sims.p) #Dispersion test significant. 
testDispersion(sims.p,alternative = "greater")
# ratioObsSim = 1.1935, p-value < 2.2e-16
# Data are over-dispersed

sims.cp <- simulateResiduals(mcp)
plot(sims.cp)
testDispersion(sims.cp,alternative = "greater")
# Still over-dispersed, but much better. 
# ratioObsSim = 1.0412, p-value < 2.2e-16

## Test random effects structures ----
m.w <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + (1|WDPA_PID), 
             data = HET,
             family = "COMPoisson",
             method = "REML")

m.b <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25), 
             data = HET,
             family = "COMPoisson",
             method = "REML")

m.c <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + (1|CONTINENT), 
             data = HET,
             family = "COMPoisson",
             method = "REML")

m.cb <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT), 
              data = HET,
              family = "COMPoisson",
              method = "REML")

m.bw <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|WDPA_PID), 
              data = HET,
              family = "COMPoisson",
              method = "REML")

m.cw <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + (1|CONTINENT) + (1|WDPA_PID), 
              data = HET,
              family = "COMPoisson",
              method = "REML")

m.all <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
               data = HET,
               family = "COMPoisson",
               method = "REML")

# No random intercept models converged - did not test random slopes

## Check best fit model for spatial autocorrelation ----
df2003 <- heterogeneityDf[(heterogeneityDf$YEAR == 0),]
coords <- data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT))

sims <- simulateResiduals(mcp)
groupedSims <- recalculateResiduals(sims, group = heterogeneityDf$WDPA_PID)
testSpatialAutocorrelation(groupedSims,coords$LONG,coords$LAT)
# Significant SAC detected
# observed = 1.2151e-01, expected = -8.3886e-05, sd = 9.4695e-04, p-value < 2.2e-16

## Test spatial correlation structures ----
# N.B: spatial models were tested on a subset of the data
# Computational demand of fitting spatial models to full dataset was not feasible

# Create a random subset of heterogeneityDf
df03 <- heterogeneityDf[(heterogeneityDf$YEAR == 0),] # Get data for 2003
s1 <- df03[sample(nrow(df03), 1163),] # Take a random sample
df03 <- df03[!(df03$WDPA_PID %in% s1$WDPA_PID),] # Remove s1 from df03
heterogeneityDfSubset <- heterogeneityDf[(heterogeneityDf$WDPA_PID %in% s1$WDPA_PID),] # Get the full data set that matches the PA in s1

# Add position and group factor to fit spatial correlation parameter
heterogeneityDfSubset$pos <- numFactor(heterogeneityDfSubset$LONG,heterogeneityDfSubset$LAT)
heterogeneityDfSubset$group <- factor(1)

# Fit non-spatial model
m <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP, 
           data = heterogeneityDfSubset,
           family = "COMPoisson",
           method = "REML")

# Fit spatial model - matern
m.m <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + Matern(1|LONG + LAT), 
             data = heterogeneityDfSubset,
             family = "COMPoisson",
             method = "REML")

# Compare model fits
anova(m, m.m)
# Matern model produced the best fit

# Check residuals
sims <- simulateResiduals(m.m)
plot(sims) 
hist(resid(m.m)) 

## Fit final model ----
# N.B: The final model was fit to multiple subsets of the full dataset
# The computational demand of fitting spatial models to the full dataset was not feasible
# Each subset was randomly sampled from the full dataset
# The outputs of this model fit to multiple subdatasets were averaged to obtain the final results
finalModel <- fitme(PA_HET ~ YEAR_s*IUCN_GROUP + Matern(1|LONG + LAT), 
                    data = heterogeneityDfSubset,
                    family = "COMPoisson",
                    method = "REML")

# Save model object
saveRDS(finalModel, file = "protected-areas-heterogeneity-final-model.rds")