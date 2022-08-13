## Mean buffer fire ignitions as a function of time
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
protectedAreaIds <- c('143001', '26606', '352021', '555586272', '689', '804')

# Load data
data <- read_csv("Data/master-dataset-time(EBV).csv")

## Prep dataset ----
# Create ignitions dataset
ignitionsDf <- data[,c(2,3,4,5,6,7,10,11,12,19,20,21)]  
# Remove NAs
ignitionsDf <- na.omit(ignitionsDf)
rm(data)

# Set BIOME to a factor data type.
ignitionsDf$BIOME25 <- factor(ignitionsDf$BIOME25)
ignitionsDf$WDPA_PID <- factor(ignitionsDf$WDPA_PID)
ignitionsDf$CONTINENT <- factor(ignitionsDf$CONTINENT)
ignitionsDf$IUCN_GROUP <- factor(ignitionsDf$IUCN_GROUP)

# Remove biome type 14 (only 3 PAs sampled in this biome)
# Consider removing biome 99 (only 5 PAs sampled).
ignitionsDf <- ignitionsDf[!(ignitionsDf$BIOME25 == '14'),]
ignitionsDf <- ignitionsDf[!(ignitionsDf$BIOME25 == '99'),]

# Remove PAs that didn't experience any fires across all study years
areasWithFire <- ignitionsDf %>% 
  group_by(WDPA_PID) %>% 
  summarise(bufferSum = sum(BUFF25_DUR_MEAN),
            paSum = sum(PA_DUR_MEAN)) %>% 
  filter(bufferSum != 0 | paSum != 0) %>% 
  select(WDPA_PID)

ignitionsDf <- ignitionsDf %>% 
  filter(WDPA_PID %in% areasWithFire$WDPA_PID)

# Identify duplicate PAs in a single year
# Duplicate coordinate points thorws an error when calculating Moran's I
df2003 <- ignitionsDf[(ignitionsDf$YEAR == '2003'),]
coords03 <- SpatialPoints(data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT)))
duplicates03 <- as.data.frame(zerodist(coords03,zero = 0.0))
# No duplicates

# Check if there are the same number of PAs sampled in each year
# Assuming the PAs in first and last year are representative of PAs sampled in all other years
df03 <- ignitionsDf[(ignitionsDf$YEAR == 2003),]
df15 <- ignitionsDf[(ignitionsDf$YEAR == 2015),]
rm(df03,df15,coords03,df2003,duplicates03)

# Check that all PAs were only sampled once per year
summ <- as.data.frame(table(ignitionsDf$WDPA_PID))
summ <- summ[(summ$Freq == 13),]

ignitionsDf <- ignitionsDf[(ignitionsDf$WDPA_PID %in% summ$Var1),] 
rm(summ)

# Filter observations by the minimum raster resolution (0.25 decimal degrees i.e. ~780 square km)
ignitionsDf <- ignitionsDf[(ignitionsDf$GIS_AREA >= 780),]

# Remove PAs that were not sampled in the PA ~ BUFF dataset.
ignitionsDf <- ignitionsDf[!(ignitionsDf$WDPA_PID %in% protectedAreaIds),]

# Add a column for fire regime (1 = biome has a natural fire regime, 0 = no natural fire regime)
ignitionsDf$FIRE <- NA
ignitionsDf$FIRE[ignitionsDf$BIOME25 == "2"|ignitionsDf$BIOME25 == "3"|ignitionsDf$BIOME25 == "4"|ignitionsDf$BIOME25 == "5"|ignitionsDf$BIOME25 == "6"|ignitionsDf$BIOME25 == "7"
                 |ignitionsDf$BIOME25 == "8"|ignitionsDf$BIOME25 == "9"|ignitionsDf$BIOME25 == "12"] <- "1"

ignitionsDf$FIRE[ignitionsDf$BIOME25 == "1"|ignitionsDf$BIOME25 == "10"|ignitionsDf$BIOME25 == "11"|ignitionsDf$BIOME25 == "13"|ignitionsDf$BIOME25 == "14"|
                   ignitionsDf$BIOME25 == "99"] <- "0"

ignitionsDf$FIRE <- as.factor(ignitionsDf$FIRE)

# Scale year (set year 1 to 0)
ignitionsDf$YEAR_s <- ignitionsDf$YEAR - 2003

## Exploratory analyses ---- 

# Check distribution of response variable
hist(ignitionsDf$BUFF25_IGN_MEAN) 

# Check for outliers
dotchart(ignitionsDf$BUFF25_IGN_MEAN, xlab = "Buffer Fire Ignitions", ylab = "Order of the data")
# No outliers

# Test fit of gaussian vs. tweedie distribution
m <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE, 
             data = ignitionsDf,
             REML = TRUE)

m1 <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE, 
              family = tweedie(),
              data = ignitionsDf,
              REML = TRUE)

# Compare model fits
AIC(m,m1)
anova(m,m1)
# Tweedie model has significantly better fit than the gaussian model. 

# Check residuals
sims1 <- simulateResiduals(m1)
plot(sims1)

## Test random effects structures ----
m <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE, 
                 family = tweedie(),
                 data = ignitionsDf,
                 REML = TRUE)

m.all <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
                 family = tweedie(),
                 data = ignitionsDf,
                 REML = TRUE)

m.bw <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1|WDPA_PID), 
                 family = tweedie(),
                 data = ignitionsDf,
                 REML = TRUE)

m.b <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25), 
                family = tweedie(),
                data = ignitionsDf,
                REML = TRUE)

m.cw <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|CONTINENT) + (1|WDPA_PID), 
                 family = tweedie(),
                 data = ignitionsDf,
                 REML = TRUE)

m.c <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|CONTINENT), 
                family = tweedie(),
                data = ignitionsDf,
                REML = TRUE)

m.bc <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1|CONTINENT), 
                 family = tweedie(),
                 data = ignitionsDf,
                 REML = TRUE)

m.w <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|WDPA_PID), 
                family = tweedie(),
                data = ignitionsDf,
                REML = TRUE)

anova(m,m.all,m.b,m.bc,m.bw,m.c,m.cw,m.w)
model.sel(m,m.all,m.b,m.bc,m.bw,m.c,m.cw,m.w)
# m.all has the lowest AIC and significantly better fit than others. 

## Test random slopes ----
m.s <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1+YEAR_s|WDPA_PID), 
                family = tweedie(),  
                data = ignitionsDf,
                REML = TRUE,
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m.s1 <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1+YEAR_s|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m.s2 <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1+YEAR_s|BIOME25) + (1|CONTINENT) + (1+YEAR_s|WDPA_PID), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m.s3 <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1+YEAR_s|CONTINENT) + (1+YEAR_s|WDPA_PID), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m.s4 <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1+YEAR_s|BIOME25) + (1|CONTINENT) + (1|WDPA_PID), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m.s5 <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE)
                 #control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m.s6 <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1|CONTINENT) + (1+YEAR_s|WDPA_PID), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# m.s5 was the only random slopes model to converge
anova(m.all, m.s5)
AIC(m.all, m.s5)
# m.s5 has the best fit and lowest AIC

## Check best fit model for spatial autocorrelation ----
df2003 <- ignitionsDf[(ignitionsDf$YEAR == 2003),]
coords <- data.frame(LONG = c(df2003$LONG), LAT = c(df2003$LAT))

sims <- simulateResiduals(m.s5)
groupedSims <- recalculateResiduals(sims, group = ignitionsDf$WDPA_PID)
testSpatialAutocorrelation(groupedSims,coords$LONG,coords$LAT) 
# Significant SAC was detected (observed = 0.2587097, expected = -0.0010331, sd = 0.0061177, p-value < 2.2e-16). 

## Test spatial correlation structures and fit final model ----
ignitionsDf$pos <- numFactor(ignitionsDf$LONG,ignitionsDf$LAT)
ignitionsDf$group <- factor(1)

m <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID), 
             family = tweedie(),  
             data = ignitionsDf,
             REML = TRUE)

m.exp <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1+YEAR_s|CONTINENT) + (1|WDPA_PID) +
                 exp(0 + pos|group), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE)

# Did not converge
m.gau <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID) +
                 gau(0 + pos|group), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Did not converge
m.mat <- glmmTMB(BUFF25_IGN_MEAN ~ YEAR_s*IUCN_GROUP*FIRE + (1|BIOME25) + (1|CONTINENT) + (1|WDPA_PID) +
                 mat(0 + pos|group), 
                 family = tweedie(),  
                 data = ignitionsDf,
                 REML = TRUE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Compare model fits
anova(m,m.exp)
AIC(m,m.exp)
# m.exp has a better fit than the non-spatial model. 

# Check residuals
sims <- simulateResiduals(m.mat)
plot(sims) 
hist(resid(m.mat))

# Save final model
saveRDS(m.exp, file = "pa-fire-ignitions-model.rds")