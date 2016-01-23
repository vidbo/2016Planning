

library("xlsx")


library(rstan)
library(reshape2)
library(plyr)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("stanmods.RData")
rm(Carmdat, Secodat)

library(rv)
var1 <- as.rv(fit1)
var2 <- as.rv(fit2)
var3 <- as.rv(fit3)

rm(fit1, fit2, fit3)
rm(h, rat, rat2, cdat, dat, datm, mdat, rdat)
rm(dat1, dat2, dat3)
##############################################
# fit PDFs for the pool depths at Santa Ynez River

poolsizes <- ddply(pooldat, .(year, pool.id), nrow)
poolsizes <- poolsizes[! is.na(poolsizes$pool.id),]
hist(poolsizes$V1, breaks=100)

hdr <- unique(pooldat[,c(1,3)])
depths <- list(NULL)
for(i in 1:dim(hdr)[1]) {
  idx <- hdr$pool.id[i]==pooldat$pool.id & hdr$year[i]==pooldat$year
  depths[[i]] <- sort(pooldat$depth[idx], decreasing=TRUE)
}

age0mn <- 100  # replace with rv
age0sd <- 10   # replace with rv, etc
age1mn <- 150
age1sd <- 15

age0 <- age1 <- 0 # numbers of age 0 and 1 fish in the pool
tarea <- rvnorm(1, 0.000000001, 0)  # cumulative territory area used

# prob of age 0 fish occuring at this depth = curve for age 0 fish, minus curve for age 1 fish, weighted by omega
p_age <- function(depth) {
  p <- as.rv(0*1:3)
  p[2] <- var2$omega /(1+exp(-var2$k[2]*(depth-var2$mu[2] )))   # note: omega controls relative abundance of age0 and age1
  p[1] <- (1          /(1+exp(-var2$k[1]*(depth-var2$mu[1] )))) - p[2] 
  p[3] <- 1 - p[2] - p[1]  # prob of no fish at this depth
  p
}
# retrieve prediction of n fish with body size f in area of size Tn
getn <- function(f, Tn) {
  mu1 <- exp(var3$beta*log(f) + var3$mu)   # using mu rather than alpha: assumes mean intercept among all years
  pred1 <- rvgamma(1, shape=var3$shape1*Tn, rate=var3$shape1*mu1)
  rvpois(n=1, lambda=pred1)  # put the poisson error back in
}

getT <- function(f, n, idx=NA) {
  if(is.na(idx)) {
    mu1 <- exp(var3$beta*log(f) + var3$mu)   # using mu rather than alpha: assumes mean intercept among all years
  } else {
    mu1 <- exp(var3$beta*log(f) + var3$alpha[idx])   # using alpha for a particular year
  }
  rvgamma(1, shape=var3$shape1*n, rate=var3$shape1/mu1)
}

# test the functions with some examples
# should be the same:
getn(120, 31)
getn(120, 31/5) + getn(120, 31/5) + getn(120, 31/5) + getn(120, 31/5) + getn(120, 31/5)

# should be the same
getT(120, 5)
getT(120, 1) + getT(120, 1) + getT(120, 1) + getT(120, 1) + getT(120, 1)

# should give back a dist. with mean=1
getn(120, getT(120, 1))


age0 <- age1 <- as.rv(0)
sizecrit <- 150
cellsize <- 1.0
# fix bad pool 147
depths[[147]] <- 0.001

capacities <- rvarray(0, dim=dim(hdr)[1])
for(i in 1:dim(hdr)[1]) {
  # depths
  d <- 100*depths[[i]][depths[[i]]>0]  # 100* is to convert m to cm
  # prob that each cell is used, as function of depth 
  p <- (1          /(1+exp(-var2$k[2]*(d-var2$mu[2] ))))
  # expectation of used territory area = sum of prob that each cell is used times area of cell
  ET <- sum(p*cellsize)   
  # expected number of fish in expected usable area ET
  capacities[i] <- getn(sizecrit, ET)
}
  
# individual pool capacities for each year, sorted  
cap2010 <- summary(sort(capacities[hdr$year==2010] ))
cap2011 <- summary(sort(capacities[hdr$year==2011] ))
cap <- rbind(cap2010, cap2011)
dimnames(cap)[[2]] <- paste0("cap", sub("%", "pct", dimnames(cap)[[2]]))
cap <- cbind(cap, year=hdr$year, idx=c(1:dim(cap2010)[1], 1:dim(cap2011)[1]))



h <- ggplot(cap, aes(x=idx)) + theme_bw()  + ylab("Pool Capacity") + xlab("Pools")
h <- h + geom_ribbon(aes(ymin=cap25pct, ymax=cap75pct, colour = factor(year)), fill="light grey") + geom_point(aes(y=cap50pct, colour = factor(year)))
h
  
bar50 <- ggplot(cap, aes(idx, cap50pct,
                       ymin = cap25pct, ymax=cap75pct, colour = factor(year))) + theme_bw()  + ylab("Pool Capacity") + xlab("Pools")
bar50 <- bar50 + geom_linerange() + geom_point()

bar50

bar50 + scale_y_log10()


# sum for each year
sum(capacities[hdr$year==2010])
sum(capacities[hdr$year==2011])

# predicted territory size
getT(150, 1)
getT(150, 1, idx=1)
getT(150, 1, idx=13)

getT(70, 1, idx=1)

write.xlsx(cap, file="capacities.xlsx")


cap2 <- summary(capacities)
dimnames(cap2)[[2]] <- paste0("cap", sub("%", "pct", dimnames(cap2)[[2]]))
cap2 <- cbind(hdr, cap2)
write.xlsx(cap2, file="capacities2.xlsx")
