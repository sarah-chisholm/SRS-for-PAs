## Mean buffer fire duration as a function of time
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

## Prep dataset ----
# Create duration dataset
durationDf <- data[,c(2,3,4,5,6,7,10,11,12,16,17,18)] 
# Remove NAs
durationDf <- na.omit(durationDf)
rm(data)

# Set variables to factor
durationDf$BIOME25 <- factor(durationDf$BIOME25)
durationDf$WDPA_PID <- factor(durationDf$WDPA_PID)
durationDf$CONTINENT <- factor(durationDf$CONTINENT)
durationDf$IUCN_GROUP <- factor(durationDf$IUCN_GROUP)

# Remove observations for biome type 14 (3 PAs sampled) and biome 99 (5 PAs sampled)
durationDf <- durationDf[!(durationDf$BIOME25 == '14'),]
durationDf <- durationDf[!(durationDf$BIOME25 == '99'),]

# Remove PAs that didn't experience any fires across all study years
areasWithFire <- durationDf %>% 
  group_by(WDPA_PID) %>% 
  summarise(bufferSum = sum(BUFF25_DUR_MEAN),
            paSum = sum(PA_DUR_MEAN)) %>% 
  filter(bufferSum != 0 | paSum != 0) %>% 
  select(WDPA_PID)

durationDf <- durationDf %>% 
  filter(WDPA_PID %in% areasWithFire$WDPA_PID)

# Identify duplicate PAs in a single year
# Duplicate coordinate points thorws an error when calculating Moran's I
df2003 <- durationDf[(durationDf$YEAR == '2003'),]
coords03 <- SpatialPoints(data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT)))
duplicates03 <- as.data.frame(zerodist(coords03,zero = 0.0))
# No duplicates were found in the dataset

# Check if there are the same number of PAs sampled in each year
# Assuming the PAs in first and last year are representative of PAs sampled in all other years
df03 <- durationDf[(durationDf$YEAR == 2003),]
df15 <- durationDf[(durationDf$YEAR == 2015),]
rm(df03,df15,coords03,df2003,duplicates03)

# Check that all PAs were only sampled once per year
summ <- as.data.frame(table(durationDf$WDPA_PID))
summ <- summ[(summ$Freq == 13),]

durationDf <- durationDf[(durationDf$WDPA_PID %in% summ$Var1),]
rm(summ)

# Filter observations by the minimum raster resolution (0.25 decimal degrees i.e. ~780 square km)
durationDf <- durationDf[(durationDf$GIS_AREA >= 780),]

# Remove PA that was not sampled in the PA ~ BUFF dataset
durationDf <- durationDf[(durationDf$WDPA_PID != '352021'),]

# Remove two outliers that were in the PA ~ BUFF dataset
durationDf <- durationDf[!(durationDf$WDPA_PID == 313739 | durationDf$WDPA_PID == 314898),]

# Add a column for fire regime (1 = biome has a natural fire regime, 0 = no natural fire regime)
durationDf$FIRE <- NA
durationDf$FIRE[durationDf$BIOME25 == "2"|durationDf$BIOME25 == "3"|durationDf$BIOME25 == "4"|
                  durationDf$BIOME25 == "5"|durationDf$BIOME25 == "6"|durationDf$BIOME25 == "7"|
                  durationDf$BIOME25 == "8"|durationDf$BIOME25 == "9"|durationDf$BIOME25 == "12"] <- "1"

durationDf$FIRE[durationDf$BIOME25 == "1"|durationDf$BIOME25 == "10"|durationDf$BIOME25 == "11"|
                  durationDf$BIOME25 == "13"|durationDf$BIOME25 == "14"|durationDf$BIOME25 == "99"] <- "0"

durationDf$FIRE <- as.factor(durationDf$FIRE)

# Scale year (set year 1 to 0)
durationDf$YEAR_s <- durationDf$YEAR - 2003

## Exploratory analyses ---- 

# Check distribution of response variable
hist(durationDf$BUFF25_DUR_MEAN)

# Check for outliers
dotchart(durationDf$BUFF25_DUR_MEAN, xlab = "Buffer Fire Duration", ylab = "Order of the data")
# No outliers

# Test fit of gaussian vs. tweedie distribution
m <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR_s*IUCN_GROUP*FIRE, 
             data = durationDf,
             REML = TRUE)

m1 <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR_s*IUCN_GROUP*FIRE, 
              family = tweedie(),  
              data = durationDf,
              REML = TRUE)

# Check residuals
sims <- simulateResiduals(m)
plot(sims)
sims1 <- simulateResiduals(m1)
plot(sims1)

# Compare model fits
AIC(m,m1)
anova(m,m1)
# Tweedie model has a significantly better fit than the gaussian model

## Test random effects structures ----
m <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE, 
              family = tweedie(),
              data = durationDf,
              REML = TRUE)

m.bw <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|WDPA_PID), 
              family = tweedie(),
              data = durationDf,
              REML = TRUE)

m.b <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25), 
              family = tweedie(),
              data = durationDf,
              REML = TRUE)

m.cw <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|CONTINENT) + (1|WDPA_PID), 
              family = tweedie(),
              data = durationDf,
              REML = TRUE)

m.c <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|CONTINENT), 
              family = tweedie(),
              data = durationDf,
              REML = TRUE)

m.bc <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|CONTINENT), 
              family = tweedie(),
              data = durationDf,
              REML = TRUE)

m.w <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|WDPA_PID), 
                 family = tweedie(),
                 data = durationDf,
                 REML = TRUE)

m.all <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
               family = tweedie(),
               data = durationDf,
               REML = TRUE)

anova(m1,m1.b,m1.bc,m1.bw,m1.c,m1.cw,m1.w)
# m.bw has the best fit (lowest AIC), but it's not significantly better than m1  

## Test random slopes ----
m.bw <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|WDPA_PID), 
                 family = tweedie(),
                 data = durationDf,
                 REML = TRUE)

m.bw.s <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1+YEAR|BIOME25) + (1+YEAR|WDPA_PID), 
                 family = tweedie(),
                 data = durationDf,
                 REML = TRUE)

m.bw.s1 <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1+YEAR|BIOME25) + (1|WDPA_PID), 
                   family = tweedie(),
                   data = durationDf,
                   REML = TRUE)

m.bw.s2 <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1+YEAR|WDPA_PID), 
                    family = tweedie(),
                    data = durationDf,
                    REML = TRUE)

# Only m.bw.s2 finished without any warnings. 
anova(m1,m1.bw,m1.bw.s,m1.bw.s1,m1.bw.s2)
# m.bw has the lowest AIC

## Check best fit model for spatial autocorrelation ----
sims <- simulateResiduals(m.bw)
df2003 <- durationDf[(durationDf$YEAR == 2003),]
coords <- data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT))

groupedSims <- recalculateResiduals(sims, group = durationDf$WDPA_PID)
testSpatialAutocorrelation(groupedSims,coords$LONG,coords$LAT) 
# Significant SAC was detected (observed = 0.2331862; p-value < 2.2e-16). 

## Test spatial correlation structures and fit final model ----
durationDf$pos <- numFactor(durationDf$LONG,durationDf$LAT)
durationDf$group <- factor(1)

# Note:  ~0 dispersion not implemented for 'tweedie' family
m <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|WDPA_PID), 
             family = tweedie(),
             data = durationDf,
             REML = TRUE)

# Did not converge
m.gau <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|WDPA_PID) + gau(0 + pos|group), 
             family = tweedie(),
             data = durationDf,
             REML = TRUE)#,
            # control=glmmTMBControl(optimizer=optim, optArgs=list(method="CG")))

m.exp <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|WDPA_PID) + exp(0 + pos|group), 
                 family = tweedie(),
                 data = durationDf,
                 REML = TRUE)#,
                 #control=glmmTMBControl(optimizer=optim, optArgs=list(method="CG")))

# Did not converge
m.mat <- glmmTMB(BUFF25_DUR_MEAN ~ YEAR*IUCN_GROUP*FIRE + (1|BIOME25) + (1|WDPA_PID) + mat(0 + pos|group), 
                 family = tweedie(),
                 data = durationDf,
                 REML = TRUE)#,
                # control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Compare model fits
anova(m, m.exp)
# exponential model is significantly better than non-spatial model
AIC(m, m.exp)
# m.exp has the lowest AIC

# Check residuals
sims <- simulateResiduals(m.exp)
plot(sims)
hist(resid(m.exp)) 

# Save final model
saveRDS(m.exp, file = "buffer-fire-duration-model.rds")