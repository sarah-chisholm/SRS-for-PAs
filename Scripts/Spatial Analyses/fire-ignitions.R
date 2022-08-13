## Mean protected area fire ignitions as a function of mean buffer fire ignitions
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
# Create ignitions dataset
ignitionsDf <- data[,c(2,3,4,5,6,9,10,11,18,19,20)]
# Remove NAs
ignitionsDf <- na.omit(ignitionsDf)
rm(data)

# Set variables to factor
ignitionsDf$BIOME25 <- factor(ignitionsDf$BIOME25)
ignitionsDf$WDPA_PID <- factor(ignitionsDf$WDPA_PID)
ignitionsDf$CONTINENT <- factor(ignitionsDf$CONTINENT)
ignitionsDf$IUCN_GROUP <- factor(ignitionsDf$IUCN_GROUP)

# Remove biome 99 (rock and ice) and 14
ignitionsDf <- ignitionsDf[ignitionsDf$BIOME != '99',] 
ignitionsDf <- ignitionsDf[ignitionsDf$BIOME != '14',] 

# Remove PAs that didn't experience any fires across all study years
areasWithFire <- ignitionsDf %>% 
  group_by(WDPA_PID) %>% 
  summarise(bufferSum = sum(BUFF25_DUR_MEAN),
            paSum = sum(PA_DUR_MEAN)) %>% 
  filter(bufferSum != 0 | paSum != 0) %>% 
  select(WDPA_PID)

ignitionsDf <- ignitionsDf %>% 
  filter(WDPA_PID %in% areasWithFire$WDPA_PID)

# Identify duplicate PAs
# Duplicate coordinate points thorws an error when calculating Moran's I
coords <- data.frame(LONG = c(ignitionsDf$LONG), LAT = c(ignitionsDf$LAT))
points <- SpatialPoints(coords)
duplicates <- as.data.frame(zerodist(points,zero = 0.0))
# No duplicates
rm(coords,points,duplicates)

# Filter observations by the minimum raster resolution (0.25 decimal degrees i.e. ~780 square km)
ignitionsDf <- ignitionsDf[(ignitionsDf$GIS_AREA >= 784),]

# Add a column for fire regime (1 = biome has a natural fire regime, 0 = no natural fire regime)
ignitionsDf$FIRE <- NA
ignitionsDf$FIRE[ignitionsDf$BIOME == "2"|ignitionsDf$BIOME == "3"|ignitionsDf$BIOME == "4"|ignitionsDf$BIOME == "5"|ignitionsDf$BIOME == "6"|ignitionsDf$BIOME == "7"
         |ignitionsDf$BIOME == "8"|ignitionsDf$BIOME == "9"|ignitionsDf$BIOME == "12"] <- "1"

ignitionsDf$FIRE[ignitionsDf$BIOME == "1"|ignitionsDf$BIOME == "10"|ignitionsDf$BIOME == "11"|ignitionsDf$BIOME == "13"|ignitionsDf$BIOME == "14"|
           ignitionsDf$BIOME == "99"] <- "0"

## Exploratory analyses ----

# Check distribution of response variable
hist(ignitionsDf$PA_IGN_MEAN)
hist(ignitionsDf$BUFF_IGN_MEAN)

#cleveland dot plot
dotchart(ignitionsDf$PA_IGN_MEAN, xlab = "PA Fire ignitions", ylab = "Order of the data")
dotchart(ignitionsDf$BUFF_IGN_MEAN, xlab = "Buffer Fire ignitions", ylab = "Order of the data")
# No outliers

# Test fit of gaussian distribution
m <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE,
              data = ignitionsDf,
              REML = TRUE)

# Check residuals
sims <- simulateResiduals(m)
plot(sims)

## Test random effects structures ----
m1 <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1|CONTINENT),
                   data = ignitionsDf,
                   REML = TRUE)

m2 <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1|BIOME),
                  data = ignitionsDf,
                  REML = TRUE)

m3 <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1|BIOME) + (1|CONTINENT),
               data = ignitionsDf,
               REML = TRUE)

## Test random slopes ----
m4 <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|CONTINENT) + (1 + BUFF_IGN_MEAN|BIOME),
              data = ignitionsDf,
              REML = TRUE,
              control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m5 <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1|CONTINENT) + (1 + BUFF_IGN_MEAN|BIOME),
                data = ignitionsDf,
                REML = TRUE,
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m6<- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|CONTINENT) + (1|BIOME),
                  data = ignitionsDf,
                  REML = TRUE)

m7<- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|BIOME),
                data = ignitionsDf,
                REML = TRUE,
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m8 <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|CONTINENT),
                data = ignitionsDf,
                REML = TRUE)

anova(m,m1,m2,m3,m4,m5,m6,m7,m8)
AIC(m,m1,m2,m3,m4,m5,m6,m7,m8)
# m8 has the best fot and lowest AIC

## Check best fit model for spatial autocorrelation ----

sims <- simulateResiduals(m8)
testSpatialAutocorrelation(m8, x = ignitionsDf$LONG, y = ignitionsDf$LAT)
# Significant SAC detected (observed = 0.0497434, expected = -0.0010331, sd = 0.0061169, p-value < 2.2e-16)

## Test spatial correlation structures ----

# Add position and group factor to fit spatial correlation parameter
ignitionsDf$pos <- numFactor(x = ignitionsDf$LONG, y = ignitionsDf$LAT)
ignitionsDf$group <- factor(1)

# Fit non-spatial model
m <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|CONTINENT),
             dispformula = ~0,
             data = ignitionsDf,
             REML = TRUE)

# Fit spatial model - exponential
m.exp <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|CONTINENT) + exp(0 + pos|group),
             dispformula = ~0,
             data = ignitionsDf,
             REML = TRUE)

# Fit spatial model - gaussian
m.gau <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|CONTINENT) + gau(0 + pos|group),
             dispformula = ~0,
             data = ignitionsDf,
             REML = TRUE)

# Fit spatial model - matern
m.mat <- glmmTMB(PA_IGN_MEAN ~ BUFF_IGN_MEAN*IUCN_GROUP*FIRE + (1 + BUFF_IGN_MEAN|CONTINENT) + mat(0 + pos|group),
             dispformula = ~0,
             data = ignitionsDf,
             REML = TRUE)

# Compare model fits
anova(m,m.gau,m.exp,m.mat)
# matern model has best fit

# Check residuals
sims <- simulateResiduals(m.mat)
plot(m.mat)
hist(resid(m.mat))

# Save final model
saveRDS(m.exp, file = "fire-ignitions-model.rds")