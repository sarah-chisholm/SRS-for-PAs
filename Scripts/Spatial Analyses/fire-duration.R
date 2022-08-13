## Mean protected area fire duration as a function of mean buffer fire duration
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

## Prep dataset ----
# Create duration dataset
durationDf <- data[,c(2,3,4,5,6,9,10,11,15,16,17)]
# Remove NAs
durationDf <- na.omit(durationDf)
rm(data)

# Set variables to factor
durationDf$BIOME25 <- factor(durationDf$BIOME25)
durationDf$WDPA_PID <- factor(durationDf$WDPA_PID)
durationDf$CONTINENT <- factor(durationDf$CONTINENT)
durationDf$IUCN_GROUP <- factor(durationDf$IUCN_GROUP)

# Remove observations for biome type 14 (3 PAs sampled) and biome 99 (5 PAs sampled)
durationDf <- durationDf[durationDf$BIOME != '14',] 
durationDf <- durationDf[durationDf$BIOME != '99',] 

# Remove PAs that didn't experience any fires across all study years
areasWithFire <- durationDf %>% 
  group_by(WDPA_PID) %>% 
  summarise(bufferSum = sum(BUFF25_DUR_MEAN),
            paSum = sum(PA_DUR_MEAN)) %>% 
  filter(bufferSum != 0 | paSum != 0) %>% 
  select(WDPA_PID)

durationDf <- durationDf %>% 
  filter(WDPA_PID %in% areasWithFire$WDPA_PID)

# Identify duplicate PAs
# Duplicate coordinate points thorws an error when calculating Moran's I
coords <- data.frame(LONG = c(durationDf$LONG), LAT = c(durationDf$LAT))
points <- SpatialPoints(coords)
duplicates <- as.data.frame(zerodist(points,zero = 0.0))
# No duplicates 
rm(coords,points,duplicates)

# Filter observations by the minimum raster resolution (0.25 decimal degrees i.e. ~780 square km)
durationDf <- durationDf[(durationDf$GIS_AREA >= 780),]

# Add a column for fire regime (1 = biome has a natural fire regime, 0 = no natural fire regime)
durationDf$FIRE <- NA
durationDf$FIRE[durationDf$BIOME == "2"|durationDf$BIOME == "3"|durationDf$BIOME == "4"|durationDf$BIOME == "5"|durationDf$BIOME == "6"|durationDf$BIOME == "7"
         |durationDf$BIOME == "8"|durationDf$BIOME == "9"|durationDf$BIOME == "12"] <- "1"

durationDf$FIRE[durationDf$BIOME == "1"|durationDf$BIOME == "10"|durationDf$BIOME == "11"|durationDf$BIOME == "13"|durationDf$BIOME == "14"|
           durationDf$BIOME == "99"] <- "0"

durationDf$FIRE <- as.factor(durationDf$FIRE)

## Exploratory analyses ----

# Check distribution of response variable
hist(durationDf$PA_DUR_MEAN) 
hist(durationDf$BUFF_DUR_MEAN)

# Check for outliers
dotchart(durationDf$PA_DUR_MEAN, xlab = "PA Fire duration", ylab = "Order of the data")
dotchart(durationDf$BUFF_DUR_MEAN, xlab = "Buffer Fire duration", ylab = "Order of the data")
# Remove outliers
durationDf <- durationDf[!(durationDf$WDPA_PID == 313739 | durationDf$WDPA_PID == 314898),]

# Test fit of gaussian distribution
m <- glmmTMB(PA_DUR_MEAN ~ BUFF_DUR_MEAN*IUCN_GROUP*FIRE,
             data = durationDf,
             REML = TRUE)

# Check residuals
sims <- simulateResiduals(m)
plot(sims)

## Test random effects structures ----
m.c <- glmmTMB(PA_DUR_MEAN ~ BUFF_DUR_MEAN*IUCN_GROUP*FIRE + (1|CONTINENT),
              data = durationDf,
              REML = TRUE)

m.b <- glmmTMB(PA_DUR_MEAN ~ BUFF_DUR_MEAN*IUCN_GROUP*FIRE + (1|BIOME),
                  data = durationDf,
                  REML = TRUE)

m.all <- glmmTMB(PA_DUR_MEAN ~ BUFF_DUR_MEAN*IUCN_GROUP*FIRE + (1|BIOME) + (1|CONTINENT),
                  data = durationDf,
                  REML = TRUE)

anova(m,m.c,m.b,m.all)
anova(m,m.c) # most simple model (m) has best fit

## Check best fit model for spatial autocorrelation ----
sims <- simulateResiduals(m)
testSpatialAutocorrelation(m, x = durationDf$LONG, y = durationDf$LAT)
# No significant SAC detected (observed = 0.0071659; p-value = 0.1806)

# Check residuals
sims <- simulateResiduals(m)
plot(sims)
hist(resid(m)) 

# Save final model
saveRDS(m.exp, file = "fire-duration-model.rds")