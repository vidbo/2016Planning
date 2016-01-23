# estimates linear fish densities and capture probabilities from Carmel River fall survey data

library("xlsx")
library(rstan)
library(reshape2)
library(plyr)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######################################
# Read in datasets
######################################

Carmdat <- read.xlsx2(file="/Users/david.boughton/Projects/StYnez/Electivity/Fall_Survey_Data.xlsx", sheetName="Sheet1", colClasses=c("Date", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
Carmdat$year <- as.integer(as.character(substr(Carmdat$Date, 1,4)))
Carmdat$area <- Carmdat$Station_Length * Carmdat$Station_Width * 0.092903 # conversion to square meters
Carmdat$silen <- Carmdat$Station_Length  * 0.3048 # conversion to meters
Carmdat <- Carmdat[! is.na(Carmdat$area),]  # remove sites without area measurements
Carmdat <- Carmdat[Carmdat$year > 1994,]    # remove sites for 1994 and earlier
Carmdat <- Carmdat[! (Carmdat$year==1997 & Carmdat$Site=="Redrock2"),]  # remove outlier

Carmdat$siteyear <- as.factor(paste(Carmdat$year, Carmdat$Site))
Carmdat$siyr <- as.integer(Carmdat$siteyear)
Carmdat$si <- as.integer(Carmdat$Site)
Carmdat$yr <- as.integer(Carmdat$year)

##############################################
# Carmel electrofishing data
cdatJ <- ddply(Carmdat, .(si, yr, siyr, Pass), summarize, lnlen=log(mean(length)), r=length(length), area=mean(area), silen=mean(silen))
ord <- order(cdatJ$siyr, cdatJ$Pass)
cdatJ <- cdatJ[ord,]
head(cdatJ)

cdatI <- ddply(Carmdat, .(si, yr, siyr), summarize, lnlen=log(mean(length)), r=length(length), area=mean(area), silen=mean(silen))
ord <- order(cdatI$siyr)
cdatI <- cdatI[ord,]

# study 1
dat3 <- list(
  siyr = cdatJ$siyr,
  Pass = cdatJ$Pass,
  r = cdatJ$r,
  area = cdatI$silen)     # note: habitat area is actually length here, not area
dat3$I = length(dat3$area)
dat3$J = length(dat3$siyr)

dat3

fit3 <- stan(file="depletion.stan", data=dat3,iter=5000, chains=4)

save.image("Plan01.RData")



library(rgdal)

frame8 <- readOGR(dsn=".", layer="SubSample8Frame")
frame4 <- readOGR(dsn=".", layer="SubSample4Frame")






library(rv)
var3 <- as.rv(fit3)

dens <- var3$N / dat3$silen   # density per unit length of stream channel
p <- var3$p                   # capture probability

library(foreign)
frame <- read.dbf("/Users/david.boughton/Projects/Carmel/Frame/frame.dbf")

vlen <- 1000 * frame$LengthKM[frame$Descriptio=="Valley"]   # lengths of sample units in Carmel Valley, in meters

efpm <- rvmedian(dens * p) # efpm = electrofished per meter, medians for first pass of all site-years

sample(efpm, length(vlen)) * vlen   # predicted fish sampled per site-year per pass

# tagging data from Oct 2015 MPWMD/NMFS survey
tag15 <- read.xlsx2(file="/Volumes/gis/development/Carmel/OmyData/MPWMD_Data/MPWMD_FallPopSurvey/Raw_Data/2015/MPWMD_FallPop2015_10282015.xlsx", 
                    sheetName="SH_Data", 
                    colClasses=c("Date", "character", "numeric", "numeric", "character", "numeric", "numeric", "numeric", "character", "character", "character", "character" ))

tdat <- ddply(tag15[! is.nan(tag15$Pass..),], .(Site, Date), summarize, ntag=length(PIT..))



sample(efpm, length(vlen)) * vlen   # predicted fish sampled per site-year per pass
sample(efpm, length(vlen)) * vlen   # predicted fish sampled per site-year per pass
sample(efpm, length(vlen)) * vlen   # predicted fish sampled per site-year per pass
sample(efpm, length(vlen)) * vlen   # predicted fish sampled per site-year per pass


mean(vlen) / mean(dat3$silen)  # mean number of MPWMD sample units per GRTS sample unit






