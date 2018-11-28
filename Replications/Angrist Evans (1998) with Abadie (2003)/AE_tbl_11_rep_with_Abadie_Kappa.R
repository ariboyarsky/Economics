# Ariel Boyarsky
# aboyarsky@uchicago.edu
#
# Code to replicate Angrist Evans (1998). The code also reestimates
# the coefficents using Abadie's (2003) kappa method.
#

# clean up, load libraries, set dir

library(ggplot2)
library(data.table)
library(knitr)
library(Hmisc) # cut2
library(Matrix)
library(stargazer)

# clean
rm(list = ls())

# make current dir wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load functions file
source("AE_rep_functions.R")

# load data
dta <- read.csv("Data/1980-AllWomen.csv")

# Table 11
res <- data.frame()

# Set up vars
cntrls <- c("agem1", "agefstm","AGEQK","AGEQ2ND","boy1st", "boy2nd", 
            "blackm", "hispm", "othracem")
C <- dta[,cntrls]
X <- dta[,c("morekids")]

################################################################
################################################################
# all women

# iterate over these vars
deps <- c("workedm", "weeksm1", "hourswm","incomem", "famincl", "nonmomil")

Z <- dta[,c("samesex")]
allw_SS <- c()
for (yvar in deps){
  Y <- dta[,c(yvar)]
  r <- ols.iv(Y,Z,X,C)
  allw_SS <- c(allw_SS, r$beta.iv[c("Z"),], r$SE[2])
}


# twins instrument
Z <- dta[,c("multi2nd")]
allw_twins <- c()
for (yvar in deps){
  Y <- dta[,c(yvar)]
  r <- ols.iv(Y,Z,X,C)
  allw_twins <- c(allw_twins, r$beta.iv[c("Z"),], r$SE[2])
}


################################################################
################################################################

# Now limit to married women
married <- subset(dta, msample==1)
# Set up vars
cntrls <- c("agem1", "agefstm","AAGE","ASEX2ND","boy1st", "boy2nd", 
            "blackm", "hispm", "othracem")
C <- married[,cntrls]
X <- married[,c("morekids")]

deps <- c("workedm", "weeksm1", "hourswm","incomem", "famincl", "nonmomil")

Z <- married[,c("samesex")]
mar_SS <- c()
for (yvar in deps){
  Y <- married[,c(yvar)]
  r <- ols.iv(Y,Z,X,C)
  mar_SS <- c(mar_SS, r$beta.iv[c("Z"),], r$SE[2])
}


# twins instrument
Z <- married[,c("multi2nd")]
mar_twins <- c()
for (yvar in deps){
  Y <- married[,c(yvar)]
  r <- ols.iv(Y,Z,X,C)
  mar_twins <- c(mar_twins, r$beta.iv[c("Z"),], r$SE[2])
}

################################################################
################################################################
# Now look at husband outcomes

# all husband dependent vars
deps <- c("workedd", "weeksd1", "hourswd","incomed")

Z <- married[,c("samesex")]
husband_SS <- c()
for (yvar in deps){
  Y <- married[,c(yvar)]
  r <- ols.iv(Y,Z,X,C)
  husband_SS <- c(husband_SS, r$beta.iv[c("Z"),], r$SE[2])
}


# twins instrument
Z <- married[,c("multi2nd")]
husband_twins <- c()
for (yvar in deps){
  Y <- married[,c(yvar)]
  r <- ols.iv(Y,Z,X,C)
  husband_twins <- c(husband_twins, r$beta.iv[c("Z"),], r$SE[2])
}


################################################################
################################################################
################################################################

# Abadie (2003) kappa estimation
# Bootsrap Parameters
B <- 20

# Panel: All Women
# Z = samesex, X = regular covariates
# Set up vars
cntrls <- c("agem1", "agefstm","AGEQK","AGEQ2ND","boy1st", "boy2nd", 
            "blackm", "hispm", "othracem")
C <- dta[,cntrls]
D <- dta[,c("morekids")]

# estimate P[Z=1|X] with samesex instrument
Z <- dta[,c("samesex")]

logit <- log.logit(C, Z)
P <- logit$Pred
dta$P <- P
# calculate abadie's kappa
dta$kappa <- 1 - (dta$morekids*(1-dta$samesex))/(1-dta$P) - 
  ((1-dta$morekids)*dta$samesex)/(dta$P)

X <- as.matrix(cbind(D, C))


# all women dependent vars
deps <- c("workedm", "weeksm1", "hourswm","incomem", "famincl", "nonmomil")

# loop through and fill up abadie results for all women panel
# Z = samesex
abadie_allw_ss <- c()
for (yvar in deps){
  bootstraps <- sapply(X   = seq(1, B),
                       FUN = boot.resample,
                       dta = dta,
                       Y  = yvar,
                       Controls  = cntrls,
                       D  = "morekids",
                       Z  = "samesex")
  
  Y <- dta[,c(yvar)]
  SE <- sqrt(var(bootstraps[2,]))
  abadie_allw_ss <- c(abadie_allw_ss, 
                      linreg(Y,X, weights = dta[,c("kappa")])$beta[2], 
                      SE)
}

# Loop through and fill up abadie results for all women panel
# Z = samesex
Z <- dta[,c("multi2nd")]

# P[Z=1|X]
logit <- log.logit(C, Z)
P <- logit$Pred
dta$P <- P

# calculate abadie's kappa
dta$kappa <- 1 - (dta$morekids*(1-dta$multi2nd))/(1-dta$P) - 
  ((1-dta$morekids)*dta$multi2nd)/(dta$P)

# X mat for beta
X <- as.matrix(cbind(D, C))

# fill up all women col with Z = twins 
abadie_allw_twins <- c()
for (yvar in deps){
  bootstraps <- sapply(X   = seq(1, B),
                       FUN = boot.resample,
                       dta = dta,
                       Y  = yvar,
                       Controls  = cntrls,
                       D  = "morekids",
                       Z  = "multi2nd")
  
  Y <- dta[,c(yvar)]
  SE <- sqrt(var(bootstraps[2,]))
  abadie_allw_twins <- c(abadie_allw_twins, 
                      linreg(Y,X, weights = dta[,c("kappa")])$beta[2], 
                      SE)
}
#############################################################################
#############################################################################

# Now we do married women panel
married <- subset(dta, msample==1)

C <- married[,cntrls]
D <- married[,c("morekids")]

# estimate P[Z=1|X] with samesex instrument
Z <- married[,c("samesex")]

logit <- log.logit(C, Z)
P <- logit$Pred
married$P <- P
# calculate abadie's kappa
married$kappa <- 1 - (married$morekids*(1-married$samesex))/(1-married$P) - 
  ((1-married$morekids)*married$samesex)/(married$P)

X <- as.matrix(cbind(D, C))


# loop through and fill up abadie results for all women panel
# Z = samesex
abadie_mar_ss <- c()
for (yvar in deps){
  bootstraps <- sapply(X   = seq(1, B),
                       FUN = boot.resample,
                       dta = married,
                       Y  = yvar,
                       Controls  = cntrls,
                       D  = "morekids",
                       Z  = "samesex")
  
  Y <- married[,c(yvar)]
  SE <- sqrt(var(bootstraps[2,]))
  abadie_mar_ss <- c(abadie_mar_ss, 
                      linreg(Y,X, weights = married[,c("kappa")])$beta[2], 
                      SE)
}

# loop through and fill up abadie results for all women panel
# Z = twnins
Z <- married[,c("multi2nd")]

# P[Z=1|X]
logit <- log.logit(C, Z)
P <- logit$Pred
married$P <- P

# calculate abadie's kappa
married$kappa <- 1 - (married$morekids*(1-married$multi2nd))/(1-married$P) - 
  ((1-married$morekids)*married$multi2nd)/(married$P)

# X mat for beta
X <- as.matrix(cbind(D, C))

# fill up all women col with Z = twins 
abadie_mar_twins <- c()
for (yvar in deps){
  bootstraps <- sapply(X   = seq(1, B),
                       FUN = boot.resample,
                       dta = married,
                       Y  = yvar,
                       Controls  = cntrls,
                       D  = "morekids",
                       Z  = "multi2nd")
  
  Y <- married[,c(yvar)]
  SE <- sqrt(var(bootstraps[2,]))
  abadie_mar_twins <- c(abadie_mar_twins, 
                         linreg(Y,X, weights = married[,c("kappa")])$beta[2], 
                         SE)
}

#############################################################################
#############################################################################
# Now we do husbands  panel so we need only switch dependent vars

deps <- c("workedd", "weeksd1", "hourswd","incomed")

# estimate P[Z=1|X] with samesex instrument
Z <- married[,c("samesex")]

logit <- log.logit(C, Z)
P <- logit$Pred
married$P <- P
# calculate abadie's kappa
married$kappa <- 1 - (married$morekids*(1-married$samesex))/(1-married$P) - 
  ((1-married$morekids)*married$samesex)/(married$P)

X <- as.matrix(cbind(D, C))


# loop through and fill up abadie results for hus panel
# Z = samesex
abadie_husband_ss <- c()
for (yvar in deps){
  bootstraps <- sapply(X   = seq(1, B),
                       FUN = boot.resample,
                       dta = married,
                       Y  = yvar,
                       Controls  = cntrls,
                       D  = "morekids",
                       Z  = "samesex")
  
  Y <- married[,c(yvar)]
  SE <- sqrt(var(bootstraps[2,]))
  abadie_husband_ss <- c(abadie_husband_ss, 
                      linreg(Y,X, weights = married[,c("kappa")])$beta[2], 
                      SE)
}

# loop through and fill up abadie results for all women panel
# Z = twnins
Z <- married[,c("multi2nd")]

# P[Z=1|X]
logit <- log.logit(C, Z)
P <- logit$Pred
married$P <- P

# calculate abadie's kappa
married$kappa <- 1 - (married$morekids*(1-married$multi2nd))/(1-married$P) - 
  ((1-married$morekids)*married$multi2nd)/(married$P)

# X mat for beta
X <- as.matrix(cbind(D, C))

# fill up hus col with Z = twins 
abadie_husband_twins <- c()
for (yvar in deps){
  bootstraps <- sapply(X   = seq(1, B),
                       FUN = boot.resample,
                       dta = married,
                       Y  = yvar,
                       Controls  = cntrls,
                       D  = "morekids",
                       Z  = "multi2nd")
  
  Y <- married[,c(yvar)]
  SE <- sqrt(var(bootstraps[2,]))
  abadie_husband_twins <- c(abadie_husband_twins, 
                         linreg(Y,X, weights = married[,c("kappa")])$beta[2], 
                         SE)
}

# pad husband with NAs
husband_SS <- c(husband_SS, NA, NA, NA, NA)
husband_twins <- c(husband_twins, NA, NA, NA, NA)

tbl11rep <- data.frame("All Women (Z = Same Sex)" = allw_SS, 
                       "All Women (Z = Twins)" = allw_twins,
                       "Married Women (Z = Same Sex)" = mar_SS,
                       "Married Women (Z = Twins)" = mar_twins,
                       "Husband Outcomes (Z = Same Sex)" = husband_SS,
                       "Husband Outcomes  (Z = Twins)" = husband_twins)

stargazer(tbl11rep, summary=FALSE)

# pad husband with NAs
abadie_husband_ss <- c(abadie_husband_ss, NA, NA, NA, NA)
abadie_husband_twins <- c(abadie_husband_twins, NA, NA, NA, NA)

abadie <- data.frame("All Women (Z = Same Sex)" = abadie_allw_ss, 
                       "All Women (Z = Twins)" = abadie_allw_twins,
                       "Married Women (Z = Same Sex)" = abadie_mar_ss,
                       "Married Women (Z = Twins)" = abadie_mar_twins,
                       "Husband Outcomes (Z = Same Sex)" = abadie_husband_ss,
                       "Husband Outcomes  (Z = Twins)" = abadie_husband_twins)

stargazer(abadie, summary=FALSE)
