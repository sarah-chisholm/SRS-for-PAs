## Buffer habitat heterogeneity as a function of time
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
sampledPAs <- read_csv("Data/habitat-heterogeneity-corresponding-pa-IDs.csv")

## Prep dataset ----
# Create heterogeneity dataset
heterogeneityDf <- data[,c(2,3,4,5,6,7,10,11,12,22,23,24)]
# Remove NAs
heterogeneityDf <- na.omit(heterogeneityDf)
rm(data)

# Set variables to factor
heterogeneityDf$BIOME25 <- factor(heterogeneityDf$BIOME25)
heterogeneityDf$WDPA_PID <- factor(heterogeneityDf$WDPA_PID)
heterogeneityDf$IUCN_CAT <- factor(heterogeneityDf$IUCN_CAT)
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
hist(heterogeneityDf$BUFF25_HET)

# Check for outliers
dotchart(heterogeneityDf$BUFF25_HET, xlab = "Buffer Heterogeneity", ylab = "Order of the data")
# No outliers

# Test fit of gaussian vs. poisson distribution
m <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP, 
             data = heterogeneityDf, 
             REML = TRUE)

m.p <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP, 
               data = heterogeneityDf, 
               family = "poisson", 
               REML = TRUE)

simsp <- simulateResiduals(m.p) 
testDispersion(m.p, alternative="less") 
# data:  simulationOutput
# ratioObsSim = 0.72701, p-value < 2.2e-16
# alternative hypothesis: less
# Data are under-dispersed

# Test a generalized poisson distribution to address issue of under-dispersion
m.gp <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP, 
                data = heterogeneityDf, 
                family = genpois(), 
                REML = TRUE)

# Check residuals
sims <- simulateResiduals(m.gp)
plot(sims)

# Compare model fits
anova(m,m.p,m.gp)
anova(m,m.gp)
# Generalized poisson model is significantly better than gaussian and has the lowest AIC

## Test random effects structures ----
m <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP, 
                data = heterogeneityDf, 
                family = genpois(), 
                REML = TRUE)

m.b <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

m.c <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|CONTINENT), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

m.w <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|WDPA_PID), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

m.bc <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

m.bw <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|WDPA_PID), 
             data = heterogeneityDf,  
             family = genpois(), 
             REML = TRUE)

m.cw <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|CONTINENT) + (1|WDPA_PID), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

m.all <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
                data = heterogeneityDf, 
                family = genpois(), 
                REML = TRUE)

anova(m,m.b,m.c,m.w,m.bc,m.bw,m.cw,m.all)
anova(m.cw,m.all)
# Full model has the significantly best fit. 

## Test random slopes ----
m <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
                 data = heterogeneityDf, 
                 family = genpois(), 
                 REML = TRUE)

m.b <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

m.c <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

m.bc <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
             data = heterogeneityDf, 
             family = genpois(), 
             REML = TRUE)

anova(m,m.b,m.c,m.bc)
anova(m.b,m.bc)
# Full model has a significantly better fit than others

## Check best fit model for spatial autocorrelation ----
df2003 <- heterogeneityDf[(heterogeneityDf$YEAR == 2003),]
coords <- data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT))

sims <- simulateResiduals(m.bc)
groupedSims <- recalculateResiduals(sims, group = heterogeneityDf$WDPA_PID)
testSpatialAutocorrelation(groupedSims,coords$LONG,coords$LAT)
# Significant SAC is detected in the model residuals (observed = 2.3366e-01, expected = -8.3886e-05, sd = 9.4698e-04, p-value < 2.2e-16)

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
m <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
                data = heterogeneityDfSubset, 
                family = genpois(), 
                REML = TRUE)

# Fit spatial model - gaussian
# N.B: Did not converge
m.g <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID) + gau(0 + pos|group), 
             data = heterogeneityDfSubset, 
             family = genpois(), 
             REML = TRUE)

# Fit spatial model - exponential
m.e <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID) + exp(0 + pos|group), 
             data = heterogeneityDfSubset, 
             family = genpois(), 
             REML = TRUE)

# Fit spatial model - matern
m.m <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID) + mat(0 + pos|group), 
             data = heterogeneityDfSubset, 
             family = genpois(), 
             REML = TRUE)

# Compare model fits
anova(m,m.e,m.m)
anova(m.e,m.m)
# Matern model produced the best fit. 

# Check residuals
sims <- simulateResiduals(m.m)
plot(sims) 
hist(resid(m.m)) 

## Fit final model ----
# N.B: The final model was fit to multiple subsets of the full dataset
# The computational demand of fitting spatial models to the full dataset was not feasible
# Each subset was randomly sampled from the full dataset
# The outputs of this model fit to multiple subdatasets were averaged to obtain the final results
finalModel <- glmmTMB(BUFF25_HET ~ YEAR_s*IUCN_GROUP + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID) + mat(0 + pos|group), 
                      data = heterogeneityDfSubset, 
                      family = genpois(), 
                      REML = TRUE)

# Save model object
saveRDS(finalModel, file = "buffer-heterogeneity-final-model.rds")