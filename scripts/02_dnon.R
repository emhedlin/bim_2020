
# Libraries ---------------------------------------------------------------

library(R2jags)
library(tidyverse)
library(matrixStats)
library(ggplot2)
library(bayesplot)
library(tidybayes)
library(ggplot2)
library(geosphere)
library(tidyverse)
library(dplyr)
library(lattice)
library(sp)
library(nuwcru)


# Load Data ---------------------------------------------------------------


dat  <- read_csv("data/2012-2020_report.csv") %>% select(-X1)


pefa <- ifelse(dat[1:162,16:39] == "PEFA", 1, 0) # convert pefa occupied to binary
rlha <- ifelse(dat[1:162,16:39] == "RLHA", 1, 0) # convert rlha occupied to binary
rapt <- ifelse(dat[1:162,16:39] == "GYRF" |
                 dat[1:162,16:39] == "PEFA" |
                 dat[1:162,16:39] == "SNOW" |
                 dat[1:162,16:39] == "RLHA", 1, 0)


# . -----------------------------------------------------------------------
# d2non work --------------------------------------------------------------

# occupancy history - we want distance to nearest neighbour values for all SiteIDs in this file


# # check to see if there are any sites in our occupancy data that AREN'T in our locations. 
# setdiff(dat$SiteID, loc$SiteID)
# 
# # add locations to our occupancy data
# dnon <- left_join(dat, loc, by = "SiteID")
# 
# occ <- read_csv("data/DNONdat.csv")
# tmp <- filter(occ, Year %in% c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))
# 
# # build dataframe with only occupied sites for each year             
# tmp$Species <- as.factor(tmp$Species)
# levels(tmp$Species) <- c("0","1", "0", "NA", "NA", "1", "1")
# 
# # create a list of sites occupied in each year
# occ <- dplyr::select(tmp, "SiteID", "Year", "Species")
# occ <- filter(occ, Species == "1")
# 
# # add locations
# occ <- left_join(occ, loc, by = "SiteID")
# 
# # create a data frame of all sites contained within our 
# sites <- dat[,1]
# all <- left_join(sites, loc, by = "SiteID")
# all <- all[1:162,]
# md <- data.frame(SiteID = md$SiteID, 
#                    md = md$`min_dist Disturb`)
# 
# md
# md <- left_join(all, md, by = "SiteID")
# Now we have two dataframes with site locations:
# one with sites that were occupied in each year -> occ
# one with all available sites -> all

# chop our occupied sites list into years, and then calculate distances


occ <- dat %>% select(-Lat_DD, -Lon_DD) %>% pivot_longer(cols = 2:26) %>% mutate(year = str_sub(name, 2,3), survey = str_sub(name, 5,5)) %>% 
  mutate(value = ifelse(value == "not detected", 0, value)) %>%
  mutate(pefa = ifelse(value == "PEFA", 1, 0), rlha = ifelse(value == "RLHA", 1, 0)) %>% select(-name, -value) %>%
  left_join(dat %>% select(SiteID, Lat_DD, Lon_DD), by = "SiteID") 

occ <- occ %>% filter(!is.na(pefa) | !is.na(rlha)) # remove sites that were not know within the given year


year_split <- occ %>% group_by(SiteID, year) %>% summarize(occ = max(pefa > 0, na.rm = TRUE)) %>% ungroup() %>% left_join(dat %>% select(SiteID, Lat_DD, Lon_DD), by = "SiteID")
year_split <- split(year_split, year_split$year)



# Create loop to calculate dnon for each year
for (i in 1:length(year_split)){
  
  occ <-  year_split[[i]] %>% filter(occ > 0)
  all <- year_split[[i]] 
  
  
  xyAll <- SpatialPointsDataFrame(
    matrix(c(all$Lat_DD, all$Lon_DD), ncol=2), data.frame(ID=seq(1:nrow(all))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  xyOcc <- SpatialPointsDataFrame(
    matrix(c(occ$Lat_DD,occ$Lon_DD), ncol=2), data.frame(ID=seq(1:nrow(occ))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  mdist <- distm(xyAll, xyOcc)
  min.d <- apply(mdist, 1, function(x) sort(x)[2])
  year_split[[i]]$min_d <- min.d
  }

year_split
subset(year_split[[1]], year_split[[1]]$occ == 1)$Lat_DD






# . -----------------------------------------------------------------------
# Load modelling data -----------------------------------------------------
dat  <- read_csv("data/2012-2020_report.csv")
dat$SiteID <- as.character(dat$SiteID)

covs <- read_csv("data/covariates.csv")
dat <- left_join(dat, covs, by = "SiteID")
md <- md[,4]


# . -----------------------------------------------------------------------
# Occurence histories -----------------------------------------------------

# split dataframe by species and convert to presence absence
pefa <- ifelse(dat[1:171,2:25] == "PEFA", 1, 0)
rlha <- ifelse(dat[1:171,2:25] == "RLHA", 1, 0)
gyrf <- ifelse(dat[1:171,2:25] == "GYRF", 1, 0)

rapt <-   ifelse(dat[1:171,2:25] == "GYRF" |
                   dat[1:171,2:25] == "PEFA" |
                   dat[1:171,2:25] == "SNOW" |
                   dat[1:171,2:25] == "RLHA", 1, 0)

# Covariate work
year <- as.character(2012:2019)    # Year Matrix
year <- matrix(year, nrow(dat), 8, byrow=TRUE)

# . -----------------------------------------------------------------------
# ______PEFA*md______ --------------------------------------------------------

#i = site
#j = survey
#t = year
dat$md <- dat$D2D
yMat    <- as.matrix(pefa)
yArray  <- array(yMat, c(162, 3, 8))	  # This is i, j, t order
yArray  <- aperm(yArray, c(1, 3, 2))		# This is i, t, j order

# minimum distance to disturbance covariate
#dat  <- as.data.frame(md)
md <- as.vector(as.numeric(dat$md))
md[is.na(md)] <- mean(md, na.rm = TRUE)
md <-as.numeric(unlist(md))

md <- (md - mean(md)) / sd(md) # mean = 47.11172, sd = 43.70102
mdMat   <- matrix(md, nrow = length(md), ncol = 8)

# jags  --------------------------------------------------------------
dim(yArray)
dim(mdMat)

# data
covData <- list(
  y          = yArray,
  md         = mdMat,	 # covariate
  nSites     = 162, 
  nYears     = 8, 
  nOccasions = 3)

# initial values
covInits <- function() list(
  psi0 = rnorm(1), psiMd = rnorm(1), # occupancy prob
  g0   = rnorm(1), gMd   = rnorm(1), # colonization
  e0   = rnorm(1), eMd   = rnorm(1), # extincation
  p = runif(1),                      # detection
  Z = matrix(rbinom((162*8), 1, 1), 162, 8)) # starting values for Z

# parameters to be estimated
covParms <- c("psi0", "psiMd", # occupancy prob
              "g0",   "gMd",   # colonization
              "e0",   "eMd",   # extinction
              "p",             # detection
              "PAO")           # proportion occupied

# Covariate model ---------------------------------------------------------


sink("covModel.txt")
cat("
    model {

    # Priors
    psi0 ~ dnorm(0, 0.1)	# normal prior on psi intercept
    psiMd ~ dnorm(0, 0.1)  # slope of md for psi
    
    g0 ~ dnorm(0, 0.1) 		# gamma intercept
    gMd ~ dnorm(0, 0.1)
    
    e0 ~ dnorm(0, 0.1)		# epsilon intercept
    eMd ~ dnorm(0, 0.1)
    
    p ~ dunif(0, 1)
    
    
    
    # Likelihood

# state process    
    for(i in 1:nSites) {
      logit(psi[i]) <- psi0 + psiMd * md[i, 1]
      Z[i, 1] ~ dbin(psi[i], 1)
    for(t in 2:nYears) {
      logit(gamma[i, t-1]) <- g0 + gMd * md[i, t]
      logit(epsilon[i, t-1]) <- e0 + eMd * md[i, t]
      muZ[i, t] <- Z[i, t-1] * (1 - epsilon[i, t-1]) + 
      (1 - Z[i, t-1]) * gamma[i, t-1]
      Z[i, t] ~ dbin(muZ[i, t], 1)
    }

# observation process
    for(t in 1:nYears) {	
      for(j in 1:nOccasions) {
        muY[i, t, j] <- Z[i, t] * p
        y[i, t, j] ~ dbin(muY[i, t, j], 1)
        } #j
      } #t
    } #i
    
    # Derived parameters (PAO = proportion of sites occupied each year)
    
    for(t in 1:nYears) {
    PAO[t] <- sum(Z[,t]) / nSites
      }
   }
    ", fill = TRUE)
sink()


# Run covariate model -----------------------------------------------------

cov1 <- jags(covData,
             covInits,
             covParms,
             "covModel.txt", 
             n.chain=3, n.iter=8000, n.burnin=4000, n.thin=10)

cov2  <- update(cov1, 
                n.iter = 10000, 
                n.thin = 10)

# cov3  <- update(cov2, 
#                  n.iter = 10000, 
#                  n.thin = 10)


cov2
traceplot(cov2,mfrow = c(2,4), col = c(hsv(0.99,1,0.1),hsv(0.33,1,0.1),hsv(.5,1,.33)))

cov2.mcmc <- as.mcmc(cov2)
cov2McList <- as.mcmc.list(cov2.mcmc) 	# Better output format for MCMC chains
covMat <- as.matrix(cov2McList)         # convert results to a matrix
cov <- covMat[,c("eMd","gMd","psiMd")]

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(cov2.mcmc,
                facet_args = list(nrow = 5, labeller = label_parsed))
p + facet_text(size = 15)
p

# information about posterior distributions
quant <- colQuantiles(cov, probs = c(0.025, 0.975))
covPost <- data.frame(parms = colnames(cov),
                      mean = colMeans(cov),
                      lq = quant[,1],
                      uq = quant[,2])
covPost <- arrange(covPost, desc(mean))


# plot posterior distribution of md ---------------------------------------

library(ggplot2)
jpeg("figures/covMd_postdist.jpg", width = 1000, height = 600)
par(cex = 1.8)
p <- ggplot(covPost, aes(parms, mean)) + ylim(-10,10)
p <- p + geom_hline(yintercept=0) + geom_errorbar(aes(ymin=lq, ymax=uq), size=1.5, width=0,color="darkgrey") + geom_point(aes(size=2)) 

p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))
p <- p + theme_bw() + xlab("") + ylab("Posterior mean and 95% credible interval")
p <- p + theme(axis.text.x=element_text(size=rel(2)),
               axis.title.x=element_text(size=rel(2)),
               axis.text.y=element_text(size=rel(2.4)),
               axis.title.y=element_text(size=rel(2.4)),
               panel.grid.minor=element_blank(),
               panel.grid.major.x=element_blank())
p <- p + coord_flip()
p
dev.off()

# . -----------------------------------------------------------------------
# ______PEFA*nn______ --------------------------------------------------------

yMat    <- as.matrix(pefa)
yArray  <- array(yMat, c(162, 3, 8))	  # This is i, j, t order
yArray  <- aperm(yArray, c(1, 3, 2))		# This is i, t, j order

# minimum distance to disturbance covariate
all <- as.data.frame(all)
nn <- all[,4:11]

nn[,1] <- (nn[,1] - mean(nn[,1])) / sd(nn[,1]) # mean = 47.11172, sd = 43.70102
nn[,2] <- (nn[,2] - mean(nn[,2])) / sd(nn[,2]) # mean = 47.11172, sd = 43.70102
nn[,3] <- (nn[,3] - mean(nn[,3])) / sd(nn[,3]) # mean = 47.11172, sd = 43.70102
nn[,4] <- (nn[,4] - mean(nn[,4])) / sd(nn[,4]) # mean = 47.11172, sd = 43.70102
nn[,5] <- (nn[,5] - mean(nn[,5])) / sd(nn[,5]) # mean = 47.11172, sd = 43.70102
nn[,6] <- (nn[,6] - mean(nn[,6])) / sd(nn[,6]) # mean = 47.11172, sd = 43.70102
nn[,7] <- (nn[,7] - mean(nn[,7])) / sd(nn[,7]) # mean = 47.11172, sd = 43.70102
nn[,8] <- (nn[,8] - mean(nn[,8])) / sd(nn[,8]) # mean = 47.11172, sd = 43.70102
nnMat   <- as.matrix(nn, nrow = 171, ncol = 8)

# jags stuff --------------------------------------------------------------

# data
covData <- list(
  y          = yArray,
  nnMat         = nnMat,	 # covariate
  nSites     = 162, 
  nYears     = 8, 
  nOccasions = 3)

# initial values
covInits <- function() list(
  psi0 = rnorm(1), psiNN = rnorm(1), # occupancy prob
  g0   = rnorm(1), gNN   = rnorm(1), # colonization
  e0   = rnorm(1), eNN   = rnorm(1), # extincation
  p = runif(1),                      # detection
  Z = matrix(rbinom((162*8), 1, 1), 162, 8)) # starting values for Z

# parameters to be estimated
covParms <- c("psi0", "psiNN", # occupancy prob
              "g0",   "gNN",   # colonization
              "e0",   "eNN",   # extinction
              "p",             # detection
              "PAO")           # proportion occupied

# Covariate model ---------------------------------------------------------


sink("covModel.txt")
cat("
    model {
    
    # Priors
    psi0 ~ dnorm(0, 0.1)	# normal prior on psi intercept
    psiNN ~ dnorm(0, 0.1)  # slope of md for psi
    
    g0 ~ dnorm(0, 0.1) 		# gamma intercept
    gNN ~ dnorm(0, 0.1)
    
    e0 ~ dnorm(0, 0.1)		# epsilon intercept
    eNN ~ dnorm(0, 0.1)
    
    p ~ dunif(0, 1)
    
    
    
    # Likelihood
    
    # state process    
    for(i in 1:nSites) {
    logit(psi[i]) <- psi0 + psiNN * nnMat[i, 1]
    Z[i, 1] ~ dbin(psi[i], 1)
    for(t in 2:nYears) {
    logit(gamma[i, t-1]) <- g0 + gNN * nnMat[i, t]
    logit(epsilon[i, t-1]) <- e0 + eNN * nnMat[i, t]
    muZ[i, t] <- Z[i, t-1] * (1 - epsilon[i, t-1]) + 
    (1 - Z[i, t-1]) * gamma[i, t-1]
    Z[i, t] ~ dbin(muZ[i, t], 1)
    }
    
    # observation process
    for(t in 1:nYears) {	
    for(j in 1:nOccasions) {
    muY[i, t, j] <- Z[i, t] * p
    y[i, t, j] ~ dbin(muY[i, t, j], 1)
    } #j
    } #t
    } #i
    
    # Derived parameters (PAO = proportion of sites occupied each year)
    
    for(t in 1:nYears) {
    PAO[t] <- sum(Z[,t]) / nSites
    }
    }
    ", fill = TRUE)
sink()


# Run covariate model -----------------------------------------------------

cov1 <- jags(covData,
             covInits,
             covParms,
             "covModel.txt", 
             n.chain=3, n.iter=8000, n.burnin=4000, n.thin=10)

cov2  <- update(cov1, 
                n.iter = 10000, 
                n.thin = 10)

# cov3  <- update(cov2, 
#                  n.iter = 10000, 
#                  n.thin = 10)


cov2
traceplot(cov2,mfrow = c(2,4), col = c(hsv(0.99,1,0.1),hsv(0.33,1,0.1),hsv(.5,1,.33)))

cov2.mcmc <- as.mcmc(cov2)
cov2McList <- as.mcmc.list(cov2.mcmc) 	# Better output format for MCMC chains
covMat <- as.matrix(cov2McList)         # convert results to a matrix
cov <- covMat[,c("eNN","gNN","psiNN")]

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(cov2.mcmc,
                facet_args = list(nrow = 5, labeller = label_parsed))
p + facet_text(size = 15)
p
pefa
# information about posterior distributions
quant <- colQuantiles(cov, probs = c(0.025, 0.975))
covPost <- data.frame(parms = colnames(cov),
                      mean = colMeans(cov),
                      lq = quant[,1],
                      uq = quant[,2])
covPost <- arrange(covPost, desc(mean))


# plot posterior distribution of NN ---------------------------------------
graphics.off()
library(ggplot2)
jpeg("figures/covNN_postdistPefa.jpg", width = 1000, height = 600)
par(cex = 1.8)
p <- ggplot(covPost, aes(parms, mean)) + ylim(-2,2)
p <- p + geom_hline(yintercept=0) + geom_errorbar(aes(ymin=lq, ymax=uq), size=1.5, width=0,color="darkgrey") + geom_point(aes(size=2)) 

p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))
p <- p + theme_bw() + xlab("") + ylab("Posterior mean and 95% confidence interval")
p <- p + theme(axis.text.x=element_text(size=rel(2)),
               axis.title.x=element_text(size=rel(2)),
               axis.text.y=element_text(size=rel(2.4)),
               axis.title.y=element_text(size=rel(2.4)),
               panel.grid.minor=element_blank(),
               panel.grid.major.x=element_blank())
p <- p + coord_flip()
p
dev.off()

# . -----------------------------------------------------------------------
# ______RLHA*md______ --------------------------------------------------------
dim(rlha)

rlha[,3:27]
yMat    <- as.matrix(rlha[,4:27])
yArray  <- array(yMat, c(91, 3, 8))	  # This is i, j, t order
yArray  <- aperm(yArray, c(1, 3, 2))		# This is i, t, j order

# minimum distance to disturbance covariate

mean <- mean(test$md)
sd <- sd(test$md)

(x*sd) + mean

md[is.na(md)] <- mean(md, na.rm = TRUE)

md <- (md - mean(md)) / sd(md) # mean = 47.11172, sd = 43.70102
mdMat   <- matrix(md, nrow = length(md), ncol = 8)
dim(mdMat)
dim(yArray)

# delete
rlha

# jags stuff --------------------------------------------------------------

# data
covData <- list(
  y          = yArray,
  md         = mdMat,	 # covariate
  nSites     = 91, 
  nYears     = 8, 
  nOccasions = 3)

# initial values
covInits <- function() list(
  psi0 = rnorm(1), psiMd = rnorm(1), # occupancy prob
  g0   = rnorm(1), gMd   = rnorm(1), # colonization
  e0   = rnorm(1), eMd   = rnorm(1), # extincation
  p = runif(1),                      # detection
  Z = matrix(rbinom((91*8), 1, 1), 91, 8)) # starting values for Z

# parameters to be estimated
covParms <- c("psi0", "psiMd", # occupancy prob
              "g0",   "gMd",   # colonization
              "e0",   "eMd",   # extinction
              "p",             # detection
              "PAO")           # proportion occupied

# Covariate model ---------------------------------------------------------


sink("covModel.txt")
cat("
    model {
    
    # Priors
    psi0 ~ dnorm(0, 0.1)	# normal prior on psi intercept
    psiMd ~ dnorm(0, 0.1)  # slope of md for psi
    
    g0 ~ dnorm(0, 0.1) 		# gamma intercept
    gMd ~ dnorm(0, 0.1)
    
    e0 ~ dnorm(0, 0.1)		# epsilon intercept
    eMd ~ dnorm(0, 0.1)
    
    p ~ dunif(0, 1)
    
    
    
    # Likelihood
    
    # state process    
    for(i in 1:nSites) {
      logit(psi[i]) <- psi0 + psiMd * md[i, 1]
    Z[i, 1] ~ dbin(psi[i], 1)
    for(t in 2:nYears) {
      logit(gamma[i, t-1]) <- g0 + gMd * md[i, t]
      logit(epsilon[i, t-1]) <- e0 + eMd * md[i, t]
    muZ[i, t] <- Z[i, t-1] * (1 - epsilon[i, t-1]) + 
    (1 - Z[i, t-1]) * gamma[i, t-1]
    Z[i, t] ~ dbin(muZ[i, t], 1)
    }
    
    # observation process
    for(t in 1:nYears) {	
    for(j in 1:nOccasions) {
    muY[i, t, j] <- Z[i, t] * p
    y[i, t, j] ~ dbin(muY[i, t, j], 1)
    } #j
    } #t
    } #i
    
    # Derived parameters (PAO = proportion of sites occupied each year)
    
    for(t in 1:nYears) {
    PAO[t] <- sum(Z[,t]) / nSites
    }
    }
    ", fill = TRUE)
sink()


# Run covariate model -----------------------------------------------------
cov1 <- jags(covData,
             covInits,
             covParms,
             "covModel.txt", 
             n.chain=3, n.iter=8000, n.burnin=4000, n.thin=10)

cov2  <- update(cov1, 
                n.iter = 10000, 
                n.thin = 10)

# cov3  <- update(cov2, 
#                  n.iter = 10000, 
#                  n.thin = 10)


cov2
traceplot(cov2,mfrow = c(2,4), col = c(hsv(0.99,1,0.1),hsv(0.33,1,0.1),hsv(.5,1,.33)))

cov2.mcmc <- as.mcmc(cov2)
cov2McList <- as.mcmc.list(cov2.mcmc) 	# Better output format for MCMC chains
covMat <- as.matrix(cov2McList)         # convert results to a matrix
cov <- covMat[,c("eMd","gMd","psiMd")]

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(cov2.mcmc,
                facet_args = list(nrow = 5, labeller = label_parsed))
p + facet_text(size = 15)
p
cov2


mdBT <- (x*sd) + mean

range(mdBT)

md # standardized cov
mdBT # backtransformed cov

covDF <- as.data.frame(covMat)
names(covDF)

predPsi <- plogis(mean(covDF$psi0) + (mean(covDF$psiMd)*md))
plot <- cbind(mdBT, predPsi)
plot(x = plot[,1], y = plot[,2], ylim = c(0,1), las = 1)

covDF
# information about posterior distributions
quant <- colQuantiles(cov, probs = c(0.025, 0.975))
covPost <- data.frame(parms = colnames(cov),
                      mean = colMeans(cov),
                      lq = quant[,1],
                      uq = quant[,2])
covPost <- arrange(covPost, desc(mean))


paoMat <- covMat[,c("PAO[1]", "PAO[2]", "PAO[3]", "PAO[4]", "PAO[5]", "PAO[6]",
                    "PAO[7]", "PAO[8]")]



md  
# plot posterior distribution of md ---------------------------------------

library(ggplot2)
jpeg("figures/covMd_postdist_rlha.jpg", width = 1000, height = 600)
par(cex = 1.8)
p <- ggplot(covPost, aes(parms, mean)) + ylim(-2.5,2.5)
p <- p + geom_hline(yintercept=0) + geom_errorbar(aes(ymin=lq, ymax=uq), size=1.5, width=0,color="darkgrey") + geom_point(aes(size=2)) 

p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))
p <- p + theme_bw() + xlab("") + ylab("Posterior mean and 95% confidence interval")
p <- p + theme(axis.text.x=element_text(size=rel(2)),
               axis.title.x=element_text(size=rel(2)),
               axis.text.y=element_text(size=rel(2.4)),
               axis.title.y=element_text(size=rel(2.4)),
               panel.grid.minor=element_blank(),
               panel.grid.major.x=element_blank())
p <- p + coord_flip()
p
dev.off()


# Derived parameters ------------------------------------------------------

epsilonMat <- covMat[,2:3]   # epsilons (extinction probs)
gammaMat   <- covMat[,9:15]	# gammas (colonization probs)
phiMat     <- 1 - epsilonMat	# phis (survival = 1 - extinction prob)


# ** Psi for year 1 -------------------------------------------------------

psiMat           <- matrix(NA, nrow=nrow(covMat), ncol=8)	# psi (occupancy)
colnames(psiMat) <- paste("psi", 1:8, sep="")
psiMat[,1]       <- covMat[,"psi"] 						            # psi for year 1 (2012)


# ** Psi for years 2:8 ----------------------------------------------------


for(t in 2:8) {
  psiMat[,t] <- psiMat[,t-1]*phiMat[,t-1] + (1 - psiMat[,t-1]) * 
    gammaMat[,t-1]
}
psiMat
psi <- as.numeric(colMeans(psiMat))

# lambda - 0.933
lambda <- round(mean(c((psi[2] / psi[1]),  
                       (psi[3] / psi[2]),  
                       (psi[4] / psi[3]),  
                       (psi[5] / psi[4]),  
                       (psi[6] / psi[5]),  
                       (psi[7] / psi[6]),  
                       (psi[8] / psi[7]))), 2) 

# lambda SD 0.092
sd <- round(sd(c((psi[2] / psi[1]),
                 (psi[3] / psi[2]),
                 (psi[4] / psi[3]),
                 (psi[5] / psi[4]),
                 (psi[6] / psi[5]),
                 (psi[7] / psi[6]),
                 (psi[8] / psi[7]))), 2)


# Plot Psi ----------------------------------------------------------------
psiCI <- apply(psiMat, 2, quantile, probs=c(0.025, 0.975)) # calculate confidence

jpeg("figures/rlha_Psi.jpg", width = 1000, height = 800)
par(cex = 1.8)
plot(2012:2019, colMeans(paoMat), ylim=c(0,1), las=1,pch="",cex = 1.8,
     xlab="Year", ylab=expression(paste('Probability of Occupancy ( ',psi,' )')), font.lab=2, bty="l",xaxt="n",yaxt="n", cex.lab = 1.2)
lines(2012:2019,obsRlha$prop,lwd=3,lty=2, col = "darkgrey")
lines(2012:2019,colMeans(psiMat),lwd=3,lty=1)
segments(2012:2019, psiCI[1,], 2012:2019, psiCI[2,], lwd = 3, col=rgb(0,0,0,0.4))
points(2012:2019, colMeans(psiMat), ylim=c(0,1), las=1,pch=16,cex = 1.8)
# plot turnover in the background
#points(2013:2019,colMeans(turnoverMat),col = "grey")
#lines(2013:2019,colMeans(turnoverMat),lwd=2,lty=3, col = "grey")
#segments(2013:2019, eqCI[1,], 2013:2019, eqCI[2,], lwd = 2, col=rgb(0,0,0,0.15))

axis(1, at=seq(2012,2019,1),labels=seq(2012,2019,1),tck=0.02,cex.axis=1)
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
axis(3, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
legend("topright", inset = c(0,.15), title = expression(paste( '   ',lambda,'')), legend = bquote(.(lambda)~" +/- "~.(sd)), bty = "n")
legend("topright", inset = c(0,0), legend = "RLHA", bty = "n", cex = 2)
legend("topleft", inset = c(0,0), legend = c("Estimated", "Observed"), lty = c(1, 2), col = c("black", "darkgrey"), lwd = c(3,3), bty = "n", cex = 1)
dev.off()

# ** Growth rate ----------------------------------------------------------

growthMat <- matrix(NA, nrow=nrow(nullMat), ncol=7)
colnames(growthMat) <- paste("growth", 1:7, sep="")
head(growthMat)

for(t in 2:8) {
  growthMat[,t-1] <- psiMat[,t] / psiMat[,t-1]
}



# ** Turnover -------------------------------------------------------------

turnoverMat <- matrix(NA, nrow=nrow(nullMat), ncol=7)
colnames(turnoverMat) <- paste("turnover", 1:7, sep="")


for(t in 2:8) {
  turnoverMat[,t-1] <- (gammaMat[,t-1] * (1 - psiMat[,t-1])) / psiMat[,t]
}

eqCI <- apply(turnoverMat, 2, quantile, probs=c(0.025, 0.975)) # calculate confidence


# Plot equilibrium probability from 2013 to 2019
plot(2013:2019, colMeans(turnoverMat), ylim=c(0,1), las=1,pch=16,cex = 1.8,font.lab=2,
     xlab="Year", ylab= "Turnover", bty="l",xaxt="n",yaxt="n", cex.lab = 1.2)

#lines(2013:2019, gamma,lwd=2,lty=1, col = "grey")
#lines(2013:2019,epsilon,lwd=2,lty=1, col = "grey")
lines(2013:2019,colMeans(turnoverMat),lwd=2,lty=2)
segments(2013:2019, eqCI[1,], 2013:2019, eqCI[2,], lwd = 2, col=rgb(0,0,0,0.4))
axis(1, at=seq(2013,2019,1),labels=seq(2013,2019,1),tck=0.02,cex.axis=1)
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
legend("bottomright", inset = c(0,.13), legend = "Turnover", bty = "n")
legend("bottomright", inset = c(0,0), legend = "RLHA", bty = "n", cex = 2)

# ** Equilibrium ----------------------------------------------------------

equilibMat <- matrix(NA, nrow=nrow(nullMat), ncol=7)
colnames(equilibMat) <- paste("equilib", 1:7, sep="")

for(t in 1:7) {
  equilibMat[,t] <- gammaMat[,t] / (gammaMat[,t] + (1 - phiMat[,t])) 
}

eqCI <- apply(equilibMat, 2, quantile, probs=c(0.025, 0.975)) # calculate confidence


# Plot equilibrium probability from 2013 to 2019
plot(2013:2019, colMeans(equilibMat), ylim=c(0,1), las=1,pch=16,cex = 1.8,font.lab=2,
     xlab="Year", ylab= "Equilibrium", bty="l",xaxt="n",yaxt="n", cex.lab = 1.2)

#lines(2013:2019, gamma,lwd=2,lty=1, col = "grey")
#lines(2013:2019,epsilon,lwd=2,lty=1, col = "grey")
lines(2013:2019,colMeans(equilibMat),lwd=2,lty=2)
segments(2013:2019, eqCI[1,], 2013:2019, eqCI[2,], lwd = 2, col=rgb(0,0,0,0.4))
axis(1, at=seq(2013,2019,1),labels=seq(2013,2019,1),tck=0.02,cex.axis=1)
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
legend("topright", inset = c(0,.15), legend = "equilibrium probability", bty = "n")
legend("topright", inset = c(0,0), legend = "RLHA", bty = "n", cex = 2)

# . -----------------------------------------------------------------------
# ______RLHA*nn______ --------------------------------------------------------

yMat    <- as.matrix(rlha)
yArray  <- array(yMat, c(162, 3, 8))	  # This is i, j, t order
yArray  <- aperm(yArray, c(1, 3, 2))		# This is i, t, j order

# minimum distance to disturbance covariate
all <- as.data.frame(all)
nn <- all[,4:11]

nn[,1] <- (nn[,1] - mean(nn[,1])) / sd(nn[,1]) # mean = 47.11172, sd = 43.70102
nn[,2] <- (nn[,2] - mean(nn[,2])) / sd(nn[,2]) # mean = 47.11172, sd = 43.70102
nn[,3] <- (nn[,3] - mean(nn[,3])) / sd(nn[,3]) # mean = 47.11172, sd = 43.70102
nn[,4] <- (nn[,4] - mean(nn[,4])) / sd(nn[,4]) # mean = 47.11172, sd = 43.70102
nn[,5] <- (nn[,5] - mean(nn[,5])) / sd(nn[,5]) # mean = 47.11172, sd = 43.70102
nn[,6] <- (nn[,6] - mean(nn[,6])) / sd(nn[,6]) # mean = 47.11172, sd = 43.70102
nn[,7] <- (nn[,7] - mean(nn[,7])) / sd(nn[,7]) # mean = 47.11172, sd = 43.70102
nn[,8] <- (nn[,8] - mean(nn[,8])) / sd(nn[,8]) # mean = 47.11172, sd = 43.70102
nnMat   <- as.matrix(nn, nrow = 171, ncol = 8)

# jags stuff --------------------------------------------------------------

# data
covData <- list(
  y          = yArray,
  nnMat         = nnMat,	 # covariate
  nSites     = 162, 
  nYears     = 8, 
  nOccasions = 3)

# initial values
covInits <- function() list(
  psi0 = rnorm(1), psiNN = rnorm(1), # occupancy prob
  g0   = rnorm(1), gNN   = rnorm(1), # colonization
  e0   = rnorm(1), eNN   = rnorm(1), # extincation
  p = runif(1),                      # detection
  Z = matrix(rbinom((162*8), 1, 1), 162, 8)) # starting values for Z

# parameters to be estimated
covParms <- c("psi0", "psiNN", # occupancy prob
              "g0",   "gNN",   # colonization
              "e0",   "eNN",   # extinction
              "p",             # detection
              "PAO")           # proportion occupied

# Covariate model ---------------------------------------------------------


sink("covModel.txt")
cat("
    model {
    
    # Priors
    psi0 ~ dnorm(0, 0.1)	# normal prior on psi intercept
    psiNN ~ dnorm(0, 0.1)  # slope of md for psi
    
    g0 ~ dnorm(0, 0.1) 		# gamma intercept
    gNN ~ dnorm(0, 0.1)
    
    e0 ~ dnorm(0, 0.1)		# epsilon intercept
    eNN ~ dnorm(0, 0.1)
    
    p ~ dunif(0, 1)
    
    
    
    # Likelihood
    
    # state process    
    for(i in 1:nSites) {
    logit(psi[i]) <- psi0 + psiNN * nnMat[i, 1]
    Z[i, 1] ~ dbin(psi[i], 1)
    for(t in 2:nYears) {
    logit(gamma[i, t-1]) <- g0 + gNN * nnMat[i, t]
    logit(epsilon[i, t-1]) <- e0 + eNN * nnMat[i, t]
    muZ[i, t] <- Z[i, t-1] * (1 - epsilon[i, t-1]) + 
    (1 - Z[i, t-1]) * gamma[i, t-1]
    Z[i, t] ~ dbin(muZ[i, t], 1)
    }
    
    # observation process
    for(t in 1:nYears) {	
    for(j in 1:nOccasions) {
    muY[i, t, j] <- Z[i, t] * p
    y[i, t, j] ~ dbin(muY[i, t, j], 1)
    } #j
    } #t
    } #i
    
    # Derived parameters (PAO = proportion of sites occupied each year)
    
    for(t in 1:nYears) {
    PAO[t] <- sum(Z[,t]) / nSites
    }
    }
    ", fill = TRUE)
sink()


# Run covariate model -----------------------------------------------------

cov1 <- jags(covData,
             covInits,
             covParms,
             "covModel.txt", 
             n.chain=3, n.iter=8000, n.burnin=4000, n.thin=10)

cov2  <- update(cov1, 
                n.iter = 10000, 
                n.thin = 10)

#cov3  <- update(cov2, 
#                 n.iter = 10000, 
#                 n.thin = 10)
#

cov2

cov2.mcmc <- as.mcmc(cov2)
cov2McList <- as.mcmc.list(cov2.mcmc) 	# Better output format for MCMC chains
covMat <- as.matrix(cov2McList)         # convert results to a matrix
cov <- covMat[,c("eNN","gNN","psiNN")]

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(cov2.mcmc,
                facet_args = list(nrow = 5, labeller = label_parsed))
p + facet_text(size = 15)
p
rlha

# information about posterior distributions
quant <- colQuantiles(cov, probs = c(0.025, 0.975))
covPost <- data.frame(parms = colnames(cov),
                      mean = colMeans(cov),
                      lq = quant[,1],
                      uq = quant[,2])
covPost <- arrange(covPost, desc(mean))


# plot posterior distribution of NN ---------------------------------------

library(ggplot2)
jpeg("figures/covNN_postdistPefa.jpg", width = 1000, height = 600)
par(cex = 1.8)
p <- ggplot(covPost, aes(parms, mean)) + ylim(-2,2)
p <- p + geom_hline(yintercept=0) + geom_errorbar(aes(ymin=lq, ymax=uq), size=1.5, width=0,color="darkgrey") + geom_point(aes(size=2)) 

p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))
p <- p + theme_bw() + xlab("") + ylab("Posterior mean and 95% confidence interval")
p <- p + theme(axis.text.x=element_text(size=rel(2)),
               axis.title.x=element_text(size=rel(2)),
               axis.text.y=element_text(size=rel(2.4)),
               axis.title.y=element_text(size=rel(2.4)),
               panel.grid.minor=element_blank(),
               panel.grid.major.x=element_blank())
p <- p + coord_flip()
p
dev.off()
# . -----------------------------------------------------------------------
# ______Guild*md______ --------------------------------------------------------

yMat    <- as.matrix(rapt)
yArray  <- array(yMat, c(162, 3, 8))	  # This is i, j, t order
yArray  <- aperm(yArray, c(1, 3, 2))		# This is i, t, j order

# minimum distance to disturbance covariate
dat  <- as.data.frame(dat)
md   <- dat$D2D
md[is.na(md)] <- mean(md, na.rm = TRUE)

md <- (md - mean(md)) / sd(md) # mean = 47.11172, sd = 43.70102
mdMat   <- matrix(md, nrow = length(md), ncol = 8)

# jags stuff --------------------------------------------------------------

# data
covData <- list(
  y          = yArray,
  md         = mdMat,	 # covariate
  nSites     = 162, 
  nYears     = 8, 
  nOccasions = 3)

# initial values
covInits <- function() list(
  psi0 = rnorm(1), psiMd = rnorm(1), # occupancy prob
  g0   = rnorm(1), gMd   = rnorm(1), # colonization
  e0   = rnorm(1), eMd   = rnorm(1), # extincation
  p = runif(1),                      # detection
  Z = matrix(rbinom((162*8), 1, 1), 162, 8)) # starting values for Z

# parameters to be estimated
covParms <- c("psi0", "psiMd", # occupancy prob
              "g0",   "gMd",   # colonization
              "e0",   "eMd",   # extinction
              "p",             # detection
              "PAO")           # proportion occupied

# Covariate model ---------------------------------------------------------


sink("covModel.txt")
cat("
    model {
    
    # Priors
    psi0 ~ dnorm(0, 0.1)	# normal prior on psi intercept
    psiMd ~ dnorm(0, 0.1)  # slope of md for psi
    
    g0 ~ dnorm(0, 0.1) 		# gamma intercept
    gMd ~ dnorm(0, 0.1)
    
    e0 ~ dnorm(0, 0.1)		# epsilon intercept
    eMd ~ dnorm(0, 0.1)
    
    p ~ dunif(0, 1)
    
    
    
    # Likelihood
    
    # state process    
    for(i in 1:nSites) {
      logit(psi[i]) <- psi0 + psiMd * md[i, 1]
    Z[i, 1] ~ dbin(psi[i], 1)
    for(t in 2:nYears) {
      logit(gamma[i, t-1]) <- g0 + gMd * md[i, t]
      logit(epsilon[i, t-1]) <- e0 + eMd * md[i, t]
    muZ[i, t] <- Z[i, t-1] * (1 - epsilon[i, t-1]) + 
    (1 - Z[i, t-1]) * gamma[i, t-1]
    Z[i, t] ~ dbin(muZ[i, t], 1)
    }
    
    # observation process
    for(t in 1:nYears) {	
    for(j in 1:nOccasions) {
    muY[i, t, j] <- Z[i, t] * p
    y[i, t, j] ~ dbin(muY[i, t, j], 1)
    } #j
    } #t
    } #i
    
    # Derived parameters (PAO = proportion of sites occupied each year)
    
    for(t in 1:nYears) {
    PAO[t] <- sum(Z[,t]) / nSites
    }
    }
    ", fill = TRUE)

sink()


# Run covariate model -----------------------------------------------------

cov1 <- jags(covData,
             covInits,
             covParms,
             "covModel.txt", 
             n.chain=3, n.iter=8000, n.burnin=4000, n.thin=10)

cov2  <- update(cov1, 
                n.iter = 10000, 
                n.thin = 10)

# cov3  <- update(cov2, 
#                  n.iter = 10000, 
#                  n.thin = 10)




traceplot(cov2,mfrow = c(2,4), col = c(hsv(0.99,1,0.1),hsv(0.33,1,0.1),hsv(.5,1,.33)))

cov2.mcmc <- as.mcmc(cov2)
cov2McList <- as.mcmc.list(cov2.mcmc) 	# Better output format for MCMC chains
covMat <- as.matrix(cov2McList)         # convert results to a matrix
cov <- covMat[,c("eMd","gMd","psiMd")]

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(cov2.mcmc,
                facet_args = list(nrow = 5, labeller = label_parsed))
p + facet_text(size = 15)
p
cov2

# information about posterior distributions
quant <- colQuantiles(cov, probs = c(0.025, 0.975))
covPost <- data.frame(parms = colnames(cov),
                      mean = colMeans(cov),
                      lq = quant[,1],
                      uq = quant[,2])
covPost <- arrange(covPost, desc(mean))


# plot posterior distribution of md ---------------------------------------

library(ggplot2)
jpeg("figures/covMd_postdist_rlha.jpg", width = 1000, height = 600)
par(cex = 1.8)
p <- ggplot(covPost, aes(parms, mean)) + ylim(-2.5,2.5)
p <- p + geom_hline(yintercept=0) + geom_errorbar(aes(ymin=lq, ymax=uq), size=1.5, width=0,color="darkgrey") + geom_point(aes(size=2)) 

p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))
p <- p + theme_bw() + xlab("") + ylab("Posterior mean and 95% confidence interval")
p <- p + theme(axis.text.x=element_text(size=rel(2)),
               axis.title.x=element_text(size=rel(2)),
               axis.text.y=element_text(size=rel(2.4)),
               axis.title.y=element_text(size=rel(2.4)),
               panel.grid.minor=element_blank(),
               panel.grid.major.x=element_blank())
p <- p + coord_flip()
p
dev.off()


# Derived parameters ------------------------------------------------------

epsilonMat <- covMat[,2:3]   # epsilons (extinction probs)
gammaMat   <- nullMat[,9:15]	# gammas (colonization probs)
phiMat     <- 1 - epsilonMat	# phis (survival = 1 - extinction prob)


# ** Psi for year 1 -------------------------------------------------------

psiMat           <- matrix(NA, nrow=nrow(nullMat), ncol=8)	# psi (occupancy)
colnames(psiMat) <- paste("psi", 1:8, sep="")
psiMat[,1]       <- nullMat[,"psi"] 						            # psi for year 1 (2012)


# ** Psi for years 2:8 ----------------------------------------------------


for(t in 2:8) {
  psiMat[,t] <- psiMat[,t-1]*phiMat[,t-1] + (1 - psiMat[,t-1]) * 
    gammaMat[,t-1]
}
psiMat
psi <- as.numeric(colMeans(psiMat))

# lambda - 0.933
lambda <- round(mean(c((psi[2] / psi[1]),  
                       (psi[3] / psi[2]),  
                       (psi[4] / psi[3]),  
                       (psi[5] / psi[4]),  
                       (psi[6] / psi[5]),  
                       (psi[7] / psi[6]),  
                       (psi[8] / psi[7]))), 2) 

# lambda SD 0.092
sd <- round(sd(c((psi[2] / psi[1]),
                 (psi[3] / psi[2]),
                 (psi[4] / psi[3]),
                 (psi[5] / psi[4]),
                 (psi[6] / psi[5]),
                 (psi[7] / psi[6]),
                 (psi[8] / psi[7]))), 2)


# Plot Psi ----------------------------------------------------------------
psiCI <- apply(psiMat, 2, quantile, probs=c(0.025, 0.975)) # calculate confidence

jpeg("figures/rlha_Psi.jpg", width = 1000, height = 800)
par(cex = 1.8)
plot(2012:2019, colMeans(paoMat), ylim=c(0,1), las=1,pch="",cex = 1.8,
     xlab="Year", ylab=expression(paste('Probability of Occupancy ( ',psi,' )')), font.lab=2, bty="l",xaxt="n",yaxt="n", cex.lab = 1.2)
lines(2012:2019,obsRlha$prop,lwd=3,lty=2, col = "darkgrey")
lines(2012:2019,colMeans(psiMat),lwd=3,lty=1)
segments(2012:2019, psiCI[1,], 2012:2019, psiCI[2,], lwd = 3, col=rgb(0,0,0,0.4))
points(2012:2019, colMeans(psiMat), ylim=c(0,1), las=1,pch=16,cex = 1.8)
# plot turnover in the background
#points(2013:2019,colMeans(turnoverMat),col = "grey")
#lines(2013:2019,colMeans(turnoverMat),lwd=2,lty=3, col = "grey")
#segments(2013:2019, eqCI[1,], 2013:2019, eqCI[2,], lwd = 2, col=rgb(0,0,0,0.15))

axis(1, at=seq(2012,2019,1),labels=seq(2012,2019,1),tck=0.02,cex.axis=1)
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
axis(3, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
legend("topright", inset = c(0,.15), title = expression(paste( '   ',lambda,'')), legend = bquote(.(lambda)~" +/- "~.(sd)), bty = "n")
legend("topright", inset = c(0,0), legend = "RLHA", bty = "n", cex = 2)
legend("topleft", inset = c(0,0), legend = c("Estimated", "Observed"), lty = c(1, 2), col = c("black", "darkgrey"), lwd = c(3,3), bty = "n", cex = 1)
dev.off()

# ** Growth rate ----------------------------------------------------------

growthMat <- matrix(NA, nrow=nrow(nullMat), ncol=7)
colnames(growthMat) <- paste("growth", 1:7, sep="")
head(growthMat)

for(t in 2:8) {
  growthMat[,t-1] <- psiMat[,t] / psiMat[,t-1]
}



# ** Turnover -------------------------------------------------------------

turnoverMat <- matrix(NA, nrow=nrow(nullMat), ncol=7)
colnames(turnoverMat) <- paste("turnover", 1:7, sep="")


for(t in 2:8) {
  turnoverMat[,t-1] <- (gammaMat[,t-1] * (1 - psiMat[,t-1])) / psiMat[,t]
}

eqCI <- apply(turnoverMat, 2, quantile, probs=c(0.025, 0.975)) # calculate confidence


# Plot equilibrium probability from 2013 to 2019
plot(2013:2019, colMeans(turnoverMat), ylim=c(0,1), las=1,pch=16,cex = 1.8,font.lab=2,
     xlab="Year", ylab= "Turnover", bty="l",xaxt="n",yaxt="n", cex.lab = 1.2)

#lines(2013:2019, gamma,lwd=2,lty=1, col = "grey")
#lines(2013:2019,epsilon,lwd=2,lty=1, col = "grey")
lines(2013:2019,colMeans(turnoverMat),lwd=2,lty=2)
segments(2013:2019, eqCI[1,], 2013:2019, eqCI[2,], lwd = 2, col=rgb(0,0,0,0.4))
axis(1, at=seq(2013,2019,1),labels=seq(2013,2019,1),tck=0.02,cex.axis=1)
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
legend("bottomright", inset = c(0,.13), legend = "Turnover", bty = "n")
legend("bottomright", inset = c(0,0), legend = "RLHA", bty = "n", cex = 2)

# ** Equilibrium ----------------------------------------------------------

equilibMat <- matrix(NA, nrow=nrow(nullMat), ncol=7)
colnames(equilibMat) <- paste("equilib", 1:7, sep="")

for(t in 1:7) {
  equilibMat[,t] <- gammaMat[,t] / (gammaMat[,t] + (1 - phiMat[,t])) 
}

eqCI <- apply(equilibMat, 2, quantile, probs=c(0.025, 0.975)) # calculate confidence


# Plot equilibrium probability from 2013 to 2019
plot(2013:2019, colMeans(equilibMat), ylim=c(0,1), las=1,pch=16,cex = 1.8,font.lab=2,
     xlab="Year", ylab= "Equilibrium", bty="l",xaxt="n",yaxt="n", cex.lab = 1.2)

#lines(2013:2019, gamma,lwd=2,lty=1, col = "grey")
#lines(2013:2019,epsilon,lwd=2,lty=1, col = "grey")
lines(2013:2019,colMeans(equilibMat),lwd=2,lty=2)
segments(2013:2019, eqCI[1,], 2013:2019, eqCI[2,], lwd = 2, col=rgb(0,0,0,0.4))
axis(1, at=seq(2013,2019,1),labels=seq(2013,2019,1),tck=0.02,cex.axis=1)
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1),las=1,tck=0.02,cex.axis=1, cex=1.4 )
legend("topright", inset = c(0,.15), legend = "equilibrium probability", bty = "n")
legend("topright", inset = c(0,0), legend = "RLHA", bty = "n", cex = 2)

# . -----------------------------------------------------------------------
# ______Guild*nn______ --------------------------------------------------------

yMat    <- as.matrix(rapt)
yArray  <- array(yMat, c(424, 3, 8))	  # This is i, j, t order
yArray  <- aperm(yArray, c(1, 3, 2))		# This is i, t, j order

# minimum distance to disturbance covariate
all <- as.data.frame(all)
nn <- all[,4:11]

nn[,1] <- (nn[,1] - mean(nn[,1])) / sd(nn[,1]) # mean = 47.11172, sd = 43.70102
nn[,2] <- (nn[,2] - mean(nn[,2])) / sd(nn[,2]) # mean = 47.11172, sd = 43.70102
nn[,3] <- (nn[,3] - mean(nn[,3])) / sd(nn[,3]) # mean = 47.11172, sd = 43.70102
nn[,4] <- (nn[,4] - mean(nn[,4])) / sd(nn[,4]) # mean = 47.11172, sd = 43.70102
nn[,5] <- (nn[,5] - mean(nn[,5])) / sd(nn[,5]) # mean = 47.11172, sd = 43.70102
nn[,6] <- (nn[,6] - mean(nn[,6])) / sd(nn[,6]) # mean = 47.11172, sd = 43.70102
nn[,7] <- (nn[,7] - mean(nn[,7])) / sd(nn[,7]) # mean = 47.11172, sd = 43.70102
nn[,8] <- (nn[,8] - mean(nn[,8])) / sd(nn[,8]) # mean = 47.11172, sd = 43.70102
nnMat   <- as.matrix(nn, nrow = 171, ncol = 8)

# jags stuff --------------------------------------------------------------

# data
covData <- list(
  y          = yArray,
  nnMat         = nnMat,	 # covariate
  nSites     = 162, 
  nYears     = 8, 
  nOccasions = 3)

# initial values
covInits <- function() list(
  psi0 = rnorm(1), psiNN = rnorm(1), # occupancy prob
  g0   = rnorm(1), gNN   = rnorm(1), # colonization
  e0   = rnorm(1), eNN   = rnorm(1), # extincation
  p = runif(1),                      # detection
  Z = matrix(rbinom((162*8), 1, 1), 162, 8)) # starting values for Z

# parameters to be estimated
covParms <- c("psi0", "psiNN", # occupancy prob
              "g0",   "gNN",   # colonization
              "e0",   "eNN",   # extinction
              "p",             # detection
              "PAO")           # proportion occupied

# Covariate model ---------------------------------------------------------


sink("covModel.txt")
cat("
    model {
    
    # Priors
    psi0 ~ dnorm(0, 0.1)	# normal prior on psi intercept
    psiNN ~ dnorm(0, 0.1)  # slope of md for psi
    
    g0 ~ dnorm(0, 0.1) 		# gamma intercept
    gNN ~ dnorm(0, 0.1)
    
    e0 ~ dnorm(0, 0.1)		# epsilon intercept
    eNN ~ dnorm(0, 0.1)
    
    p ~ dunif(0, 1)
    
    
    
    # Likelihood
    
    # state process    
    for(i in 1:nSites) {
    logit(psi[i]) <- psi0 + psiNN * nnMat[i, 1]
    Z[i, 1] ~ dbin(psi[i], 1)
    for(t in 2:nYears) {
    logit(gamma[i, t-1]) <- g0 + gNN * nnMat[i, t]
    logit(epsilon[i, t-1]) <- e0 + eNN * nnMat[i, t]
    muZ[i, t] <- Z[i, t-1] * (1 - epsilon[i, t-1]) + 
    (1 - Z[i, t-1]) * gamma[i, t-1]
    Z[i, t] ~ dbin(muZ[i, t], 1)
    }
    
    # observation process
    for(t in 1:nYears) {	
    for(j in 1:nOccasions) {
    muY[i, t, j] <- Z[i, t] * p
    y[i, t, j] ~ dbin(muY[i, t, j], 1)
    } #j
    } #t
    } #i
    
    # Derived parameters (PAO = proportion of sites occupied each year)
    
    for(t in 1:nYears) {
    PAO[t] <- sum(Z[,t]) / nSites
    }
    }
    ", fill = TRUE)
sink()


# Run covariate model -----------------------------------------------------

cov1 <- jags(covData,
             covInits,
             covParms,
             "covModel.txt", 
             n.chain=3, n.iter=8000, n.burnin=4000, n.thin=10)

cov2  <- update(cov1, 
                n.iter = 10000, 
                n.thin = 10)

cov3  <- update(cov2, 
                n.iter = 10000, 
                n.thin = 10)


cov2

cov2.mcmc <- as.mcmc(cov2)
cov2McList <- as.mcmc.list(cov2.mcmc) 	# Better output format for MCMC chains
covMat <- as.matrix(cov2McList)         # convert results to a matrix
cov <- covMat[,c("eNN","gNN","psiNN")]

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(cov2.mcmc,
                facet_args = list(nrow = 5, labeller = label_parsed))
p + facet_text(size = 15)
p
rapt

cov3
# information about posterior distributions
quant <- colQuantiles(cov, probs = c(0.025, 0.975))
covPost <- data.frame(parms = colnames(cov),
                      mean = colMeans(cov),
                      lq = quant[,1],
                      uq = quant[,2])
covPost <- arrange(covPost, desc(mean))


# plot posterior distribution of NN ---------------------------------------

library(ggplot2)
jpeg("figures/covNN_postdistPefa.jpg", width = 1000, height = 600)
par(cex = 1.8)
p <- ggplot(covPost, aes(parms, mean)) + ylim(-2,2)
p <- p + geom_hline(yintercept=0) + geom_errorbar(aes(ymin=lq, ymax=uq), size=1.5, width=0,color="darkgrey") + geom_point(aes(size=2)) 

p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))
p <- p + theme_bw() + xlab("") + ylab("Posterior mean and 95% confidence interval")
p <- p + theme(axis.text.x=element_text(size=rel(2)),
               axis.title.x=element_text(size=rel(2)),
               axis.text.y=element_text(size=rel(2.4)),
               axis.title.y=element_text(size=rel(2.4)),
               panel.grid.minor=element_blank(),
               panel.grid.major.x=element_blank())
p <- p + coord_flip()
p
dev.off()

