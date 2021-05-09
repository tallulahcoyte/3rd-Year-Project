# loading in packages

install.packages("piecewiseSEM")
library(piecewiseSEM)
install.packages("nlme")
library(nlme)
install.packages("lme4")
library(lme4)
install.packages("ggplot2")
library(ggplot2)
install.packages("lmtest")
library(lmtest)
install.packages("glmmADMB", repos="http://www.math.mcmaster.ca/bolker/R",type="source")
library(glmmADMB)
install.packages("devtools")
library(devtools)
install.packages("usethis")
library(usethis)
devtools::install_github("hohenstein/remef")
library(remef)
install.packages("ggeffects")
library(ggeffects)
install.packages("effects")
library(effects)

#  loading in data 

shipley <- read.csv("shipley.csv")
View (shipley)

hist(shipley$DD, xlab = "Degree day of bud burst", ylab = "Frequency", main = "Histogram of degree day of bud burst")
hist(shipley$Date, xlab = "Julian date of bud burst", ylab = "Frequency", main = "Histogram of Julian date of bud burst")
hist(shipley$Growth, xlab = "Increase in diameter growth of the tree", ylab = "Frequency", main = "Histogram of increase in diameter growth of tree")
hist(shipley$Survival, xlab = "Probability of survival until the next year", ylab = "Frequency", main = "Histogram of probability of survival until the next year")

# Fitting independence claims by 'hand'

#1 (lat,Date) | DD
fit1<-lme(Date~DD+lat,data=shipley,random=~1|site/tree,na.action=na.omit) 
summary(fit1)
#2 (lat, Growth) | DD
fit2<-lme(Growth~DD+lat,data=shipley,random=~1|site/tree,na.action=na.omit) 
summary(fit2)
#3 (lat, Survival) | Growth
fit3<-lme(Survival~Growth+lat,data=shipley,random=~1|site/tree,na.action=na.omit)
summary(fit3)
#4 (DD, Growth) | (lat, Date)
fit4<-lme(Growth~DD+lat+Date,data=shipley,random=~1|site/tree, na.action=na.omit)
summary(fit4)
#5 (DD, Survival) | (lat, Growth)
fit5<-lme(Survival~DD+lat+Growth,data=shipley,random=~1|site/tree,na.action=na.omit)
summary(fit5)
#6 (Date, Survival) | (DD, Growth)
fit6<-lme(Survival~DD+Date+Growth,data=shipley,random=~1|site/tree,na.action=na.omit)
summary(fit6)

# Fitting independence claims with piecewiseSEM package

shipley.list <- list(
  
  lme(DD ~ lat, random = ~ 1 | site / tree, na.action = na.omit, 
      data = shipley),
  
  lme(Date ~ DD, random = ~ 1 | site / tree, na.action = na.omit, 
      data = shipley),
  
  lme(Growth ~ Date, random = ~ 1 | site / tree, na.action = na.omit, 
      data = shipley),
  
  glmer(Live ~ Growth + (1 | site) + (1 | tree), 
        family = binomial(link = "logit"), data = shipley) 
  
)
shipley.psem <- as.psem(shipley.list)
new.summary <- summary(shipley.psem, .progressBar = F)
new.summary

# PIECEWISE SEM FIRST MODEL

caseysemdata <- read.csv("caseysemdata.csv")
caseysemdata$numReef <- factor(caseysemdata$Reef, levels = c("HP3", "Bell", "HP2", "21-466", "21-544", "21-507", "Recreation", "21-500", "Frigate", "Carter", "Hilder", "Yonge", "Day", "Jewell", "Hicks"), labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
caseysemdata$numReef <- as.numeric(caseysemdata$numReef)
caseysemdata <- na.omit(caseysemdata)
caseysemdata$All_Algae <- caseysemdata$CCA + caseysemdata$Macro + caseysemdata$Turf
caseysemdata$Pred <- caseysemdata$Tpred + caseysemdata$Ntpred
View(caseysemdata)

lmApex <- lm(Apex ~ Zone, data=caseysemdata)
summary(lmApex)
caseysemdata$apexresiduals <- lmApex$residuals
var.func.apex <- lm(apexresiduals^2 ~ Zone, data=caseysemdata)
summary(var.func.apex)

lmmeso <- lm(Tpred ~ Zone, data=caseysemdata)
caseysemdata$mesoresiduals <- lmmeso$residuals
var.func.meso <- lm(mesoresiduals^2 ~ Zone, data=caseysemdata)
summary(var.func.meso)

lmMHerb <- lm(Mherb ~ Apex + Tpred, data=caseysemdata)
summary(lmMHerb)
caseysemdata$mherbresiduals <- lmMHerb$residuals
var.func.mherb <- lm(mherbresiduals^2 ~ Apex + Tpred, data=caseysemdata)
summary(var.func.mherb)

lmTgraz <- lm(Tgraz ~ Apex + Tpred + Mherb, data=caseysemdata)
caseysemdata$tgrazresiduals <- lmTgraz$residuals
var.func.tgraz <- lm(tgrazresiduals^2 ~ Apex + Tpred + Mherb, data=caseysemdata)
summary(var.func.tgraz)

lmAlgae <- lm( All_Algae ~ Mherb + Tgraz, data = caseysemdata)
caseysemdata$algaeresidual <- lmAlgae$residuals
var.func.algae <- lm(algaeresidual^2 ~ Mherb + Tgraz, data = caseysemdata)
summary(var.func.algae)

lmCoral <- lm( Coral ~ Mherb + Tgraz + All_Algae, data = caseysemdata)
caseysemdata$coralresidual <- lmCoral$residuals
var.func.coral <- lm(coralresidual^2 ~ Mherb + Tgraz + All_Algae, data = caseysemdata)
summary(var.func.coral)

lmPred <- lm(Pred ~ Apex, data=caseysemdata)
caseysemdata$predresidual <- lmPred$residuals
var.func.pred <- lm(predresidual^2 ~ Apex, data=caseysemdata)
summary(var.func.pred)

lmNtpred <- lm(Ntpred ~ Zone + Tpred, data=caseysemdata)
caseysemdata$ntpredresiduals <- lmNtpred$residuals
var.func.ntpred <- lm(ntpredresiduals^2 ~ Zone + Tpred, data=caseysemdata)
summary(var.func.ntpred)

bptest(lmMHerb)
bptest(lmApex)
bptest(lmTgraz)
bptest(lmmeso)
bptest(lmPred)
bptest(lmNtpred)

qchisq(.95, df = 1)
qchisq(.95, df = 2)
qchisq(.95, df = 3)

ggplot(data = caseysemdata, aes(y = mherbresiduals, x = Apex + Tpred), xlab = "Predictor variables",  ylab = "Residual for Mobile Herbivores") + geom_point(col = 'blue') + geom_abline(slope = 0)  
ggplot(data = caseysemdata, aes(y = apexresiduals, x =  Zone)) + geom_point(col = 'red') + geom_abline(slope = 0)
ggplot(data = caseysemdata, aes(y = tgrazresiduals, x =  Apex + Tpred + Mherb)) + geom_point(col = 'green') + geom_abline(slope = 0)
ggplot(data = caseysemdata, aes(y = mesoresiduals, x = Zone)) + geom_point(col = 'orange') + geom_abline(slope = 0)
ggplot(data = caseysemdata, aes(y = algaeresidual, x =  Mherb + Tgraz)) + geom_point(col = 'pink') + geom_abline(slope = 0)
ggplot(data = caseysemdata, aes(y = coralresidual, x = Mherb + Tgraz + All_Algae)) + geom_point(col = 'blueviolet') + geom_abline(slope = 0)

# TRANSFORMING VARIABLES

hist(caseysemdata$Pred, xlab = "All mesopredator biomass", ylab = "Frequency", main = "Histogram of all mesopredator biomass")
hist(caseysemdata$Ntpred, xlab = "Non-targeted mesopredator biomass", ylab = "Frequency", main = "Histogram of Non-targeted mesopredator biomass")

caseysemdata$CoralSR <- sqrt(caseysemdata$Coral)

caseysemdata$Pred <- caseysemdata$Tpred + caseysemdata$Ntpred


# ALL ZONES, TARGETED MESO


# (X1, X4) | (X2, X3)

lme(MherbSR ~ Zone + ApexSR + TpredSR, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X1, X5) | { X2, X3 }

lme(TgrazSR ~ Zone + ApexSR + TpredSR, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X1, X6) | { X4, X5 }

glmer( All_Algae ~ Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X1, X7) | { X4, X5, X6 }

lme(CoralSR ~ Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X2, X6) | { X1, X4, X5 }

glmer( All_Algae ~ ApexSR + Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X2, X7) | { X1, X4, X5, X6 }

lme(CoralSR ~ ApexSR + Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X3, X6) | { X1, X2, X4, X5 }

glmer( All_Algae ~ TpredSR + ApexSR + Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X3, X7) | { X1, X2, X4, X5, X6 }

lme(CoralSR ~ TpredSR + ApexSR + Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)


# All meso 


# (X1, X4) | (X2, X9)

lme(MherbSR ~ Zone + ApexSR + Pred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X1, X5) | { X2, X9 }

lme(TgrazSR ~ Zone + ApexSR + Pred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X1, X6) | { X4, X5 }

glmer( All_Algae ~ Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X1, X7) | { X4, X5, X6 }

lme(CoralSR ~ Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X2, X6) | { X1, X4, X5 }

glmer( All_Algae ~ ApexSR + Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X2, X7) | { X1, X4, X5, X6 }

lme(CoralSR ~ ApexSR + Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X9, X6) | { X1, X2, X4, X5 }

glmer( All_Algae ~ Pred + ApexSR + Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X9, X7) | { X1, X2, X4, X5, X6 }

lme(CoralSR ~ Pred + ApexSR + Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)


# Targeted and non targeted mesopred 


# (X1, X4) | (X2, X3, X8)

lme(MherbSR ~ Zone + ApexSR + Pred + Ntpred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X1, X5) | { X2, X3, X8 }

lme(TgrazSR ~ Zone + ApexSR + Pred + Ntpred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X1, X6) | { X4, X5 }

glmer( All_Algae ~ Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X1, X7) | { X4, X5, X6 }

lme(CoralSR ~ Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X2, X6) | { X1, X4, X5 }

glmer( All_Algae ~ ApexSR + Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X2, X7) | { X1, X4, X5, X6 }

lme(CoralSR ~ ApexSR + Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X3, X6) | { X1, X2, X4, X5 }

glmer( All_Algae ~ Pred + ApexSR + Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X8, X6) | { X1, X2, X4, X5 }

glmer( All_Algae ~ Ntpred + ApexSR + Zone + MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata)

# (X3, X7) | { X1, X2, X4, X5, X6 }

lme(CoralSR ~ Pred + ApexSR + Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# (X8, X7) | { X1, X2, X4, X5, X6 }

lme(CoralSR ~ Ntpred + ApexSR + Zone + MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata)

# first model
one.trophic.fished.list <- list(
  glm( ApexSR ~ Zone, family=poisson(link="log"), na.action=na.omit, data = caseysemdata),
  lme( TpredSR ~ Zone + ApexSR, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  lme( MherbSR ~ ApexSR + TpredSR, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  lme( TgrazSR ~ ApexSR + TpredSR, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  glmer( All_Algae ~ MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata),
  lme( CoralSR ~ MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  ApexSR %~~% TpredSR,
  MherbSR %~~% TgrazSR
)
one.trophic.fished.psem <- as.psem(one.trophic.fished.list)
one.summary <- summary(one.trophic.fished.psem, .progressBar = F)
one.summary

# Fisher's C = 21.456 with P-value = 0.162 and on 16 degrees of freedom

# second model
two.trophic.fished.list <- list(
  glm( ApexSR ~ Zone, family=poisson(link="log"), na.action=na.omit, data = caseysemdata),
  glmer( Pred ~ Zone + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata),
  lme( MherbSR ~ ApexSR + Pred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  lme( TgrazSR ~ ApexSR + Pred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  glmer( All_Algae ~ MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata),
  lme( CoralSR ~ MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  ApexSR %~~% Pred,
  MherbSR %~~% TgrazSR
)
two.trophic.fished.psem <- as.psem(two.trophic.fished.list)
two.summary <- summary(two.trophic.fished.psem, .progressBar = F)
two.summary
# All_Algae ~ ApexSR + ...      coef 226    -2.7304  0.0063 **
# Fisher's C = 25.356 with P-value = 0.064 and on 16 degrees of freedom

# third model 
three.trophic.fished.list <- list(
  glm( ApexSR ~ Zone, family=poisson(link="log"), na.action=na.omit, data = caseysemdata),
  lme( TpredSR ~ Zone + ApexSR, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  lme( Ntpred ~ Zone + TpredSR, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  lme( MherbSR ~ ApexSR + TpredSR + Ntpred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  lme( TgrazSR ~ ApexSR + TpredSR + Ntpred, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  glmer( All_Algae ~ MherbSR + TgrazSR + (1 | numReef / Site), family = "poisson", na.action = na.omit, data = caseysemdata),
  lme( CoralSR ~ MherbSR + TgrazSR + All_Algae, random = ~ 1 | Reef / Site, na.action = na.omit, data = caseysemdata),
  ApexSR %~~% TpredSR,
  ApexSR %~~% Ntpred,
  MherbSR %~~% TgrazSR
)
three.trophic.fished.psem <- as.psem(three.trophic.fished.list)
three.summary <- summary(three.trophic.fished.psem, .progressBar = F)
three.summary
trophic$numReef <- as.numeric(trophic$numReef)
# All_Algae ~ ApexSR + ...      coef 226    -2.7304  0.0063 **
# Fisher's C = 25.368 with P-value = 0.188 and on 20 degrees of freedom


