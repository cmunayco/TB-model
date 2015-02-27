library(sensitivity)
library(rgl)
library(Rglpk)
library(pse)
library(FME)

pars <- list(bpi = 4400, beta = 0.00005,
 mu = 0.0222, v = 0.00256, p = 0.05, f=0.70, 
 q=0.85, mut=0.139, w=0.005, c=0.058)

## solving model
solveTB <- function(pars, times=seq(0,100,by=1)) {
 derivs <- function(t, state, pars) { # returns rate of change
 with(as.list(c(state, pars)), {
 dS <- bpi - beta*S*Ti - mu*S
 dL <- (1 - p)*beta*S*Ti - (v + mu)*L
 dTi <- p*f*beta*S*Ti + q*v*L + w*R - (mu + mut + c)*Ti
 dTn <- p*(1 - f)*beta*S*Ti + (1 - q)*v*L + w*R - (mu + mut + c)*Tn
 dR <- c*Ti + c*Tn - (2*w + mu)*R 
 return(list(c(dS, dL, dTi, dTn, dR), TOC = S + L + Ti + Tn + R))
 })
 }
 state <- c(S = 199999, L = 0, Ti=1, Tn=0, R=0)
 ## ode solves the model by integration...
 return(ode(y = state, times = times, func = derivs, parms = pars))
 }

out<-solveTB(pars) ## output

## Figure 1 All components
quartz(width=10, height=6, pointsize=10)
matplot(out[,1], out[,-1], type = "l", lty = 1:6, lwd = c(2, 2, 2,2,2,1),
col = c("red", "blue", "green", "purple","yellow", "black"), xlab = "time, year", ylab = "Number of Individuals")
legend("topright", c("Susceptible", "Latently infected","Infectious Tuberculosis", 
										 "Non-Infectious Tuberculosis","Recovered", "TOC"), 
lty = 1:6, lwd = c(2, 2, 2, 2, 2, 1), col = c("red", "blue", "green", "purple","yellow", "black"))


## Figure 2. Infectious Tuberculosis and Non-Infectious Tuberculosis cases

quartz(width=10, height=6, pointsize=10)
matplot(out[,1], out[,c(4,5)], type = "l", lty = 1:2, lwd = c(2,1),
col = c("red", "blue"), xlab = "time, year", ylab = "Number of Individuals")
legend("topright", c("Infectious Tuberculosis", 
										 "Non-Infectious Tuberculosis"), 
lty = 1:6, lwd = c(2, 1), col = c("red", "blue"))


quartz(width=10, height=6, pointsize=10)
matplot(out[,1], ((out[,3])/out[,2])*10000, type = "l", lty = 1, lwd = 1,
col = "blue", xlab = "time, year", ylab = "Incidence of Infection per 10,000 pop")
legend("topright", "L/S", 
lty = 1, lwd = 1, col="blue")

quartz(width=10, height=6, pointsize=10)
matplot(out[,1], ((out[,3]+out[,4]+out[,5]+out[,6])/out[,7])*100, type = "l", lty = 1, lwd = 1,
col = "blue", xlab = "time, year", ylab = "Prevalence of Infection (Percent)")
legend("topright", "(L+Ti+Tn+R)/N", 
lty = 1, lwd = 1, col="blue")


quartz(width=10, height=6, pointsize=10)
matplot(out[,1], ((out[,4] + out[,5] + out[,6])/out[,7])*100, type = "l", lty = 1, lwd = 1,
col = "red", xlab = "time, year", ylab = "Prevalence of Disease")
legend("topright", "(Ti + Tn + R)/S", 
lty = 1, lwd = 1, col="red")

quartz(width=10, height=6, pointsize=10)
matplot(out[,1], ((out[,4] + out[,5])/out[,2])*10000, type = "l", lty = 1:2, lwd = c(2, 1),
col = "red", xlab = "time, year", ylab = "Incidence of Disease per 10,000 pop")
legend("topright", "(Ti + Tn)/S", 
		lty = 1, lwd = 1, col = "red")


## GLobal sensitivity

parRanges <- data.frame(min = c(3000, 0.00001, 0.01,0.00256, 0.0, 0.5,0.5, 0.058, 0.0, 0.021), 
												max = c(6000, 0.00009, 0.04,0.00527, 0.30, 0.85,0.85, 0.461, 0.01, 0.086))
rownames(parRanges) <- c("bpi", "beta", "mu", "v", "p", "f", "q", "mut", "w", "c")
parRanges

## Local sensitivity for "beta"

tout <- 0:50
print(system.time(
beta <- sensRange(func = solveTB, parms = pars, dist = "grid",
sensvar = c("L","Ti","Tn"), parRange = parRanges[2,], num = 1000)
))
head(summary(beta))

summ.beta <- summary(beta)
quartz(width=10, height=6, pointsize=10)
par(mfrow=c(2, 3))
plot(summ.beta, xlab = "time, year", ylab = "Number of individuals",
 legpos = "topright", mfrow = NULL)
plot(summ.beta, xlab = "time, year", ylab = "Number of individuals", mfrow = NULL,
 quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, "Sensitivity to beta (Transmission coefficient)", cex = 1.25)
par(mfrow = c(1, 1))

## Local sensitivity for "p"

tout <- 0:50
print(system.time(
p <- sensRange(func = solveTB, parms = pars, dist = "grid",
sensvar = c("L","Ti","Tn"), parRange = parRanges[5,], num = 1000)
))
head(summary(p))

summ.p <- summary(p)
quartz(width=10, height=6, pointsize=10)
par(mfrow=c(2, 3))
plot(summ.p, xlab = "time, year", ylab = "Number of individuals",
 legpos = "topright", mfrow = NULL)
plot(summ.p, xlab = "time, year", ylab = "Number of individuals", mfrow = NULL,
 quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, 
"Sensitivity to p (Proportion of new infections that become TB disease within 1 year)", cex = 1.25)
par(mfrow = c(1, 1))

## Local sensitivity for "f"

tout <- 0:50
print(system.time(
f <- sensRange(func = solveTB, parms = pars, dist = "grid",
sensvar = c("L","Ti","Tn"), parRange = parRanges[6,], num = 1000)
))
head(summary(f))

summ.f <- summary(f)
quartz(width=10, height=6, pointsize=10)
par(mfrow=c(2, 3))
plot(summ.f, xlab = "time, year", ylab = "Number of individuals",
 legpos = "topright", mfrow = NULL)
plot(summ.f, xlab = "time, year", ylab = "Number of individuals", mfrow = NULL,
 quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, 
"Sensitivity to f (Probability of developing infectious TB, given that one develops fast TB)", cex = 1.25)
par(mfrow = c(1, 1))

## Local sensitivity for "q"

tout <- 0:50
print(system.time(
q <- sensRange(func = solveTB, parms = pars, dist = "grid",
sensvar = c("L","Ti","Tn"), parRange = parRanges[7,], num = 1000)
))
head(summary(q))

summ.q <- summary(q)
quartz(width=10, height=6, pointsize=10)
par(mfrow=c(2, 3))
plot(summ.q, xlab = "time, year", ylab = "Number of individuals",
 legpos = "topright", mfrow = NULL)
plot(summ.q, xlab = "time, year", ylab = "Number of individuals", mfrow = NULL,
 quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, 
"Sensitivity to q (Probability of developing infectious TB, given that one develops slow TB)", cex = 1.25)
par(mfrow = c(1, 1))


## Local sensitivity for "mut"
tout <- 0:50
print(system.time(
mut <- sensRange(func = solveTB, parms = pars, dist = "grid",
sensvar = c("L","Ti","Tn"), parRange = parRanges[8,], num = 1000)
))
head(summary(mut))

summ.mut <- summary(mut)
quartz(width=10, height=6, pointsize=10)
par(mfrow=c(2, 3))
plot(summ.mut, xlab = "time, year", ylab = "Number of individuals",
 legpos = "topright", mfrow = NULL)
plot(summ.mut, xlab = "time, year", ylab = "Number of individuals", mfrow = NULL,
 quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, 
"Sensitivity to mut (Mortality rate due to TB (per capita))", cex = 1.25)
par(mfrow = c(1, 1))

## Local sensitivity for "w"
tout <- 0:50
print(system.time(
w <- sensRange(func = solveTB, parms = pars, dist = "grid",
sensvar = c("L","Ti","Tn"), parRange = parRanges[9,], num = 1000)
))
head(summary(w))

summ.w <- summary(w)
quartz(width=10, height=6, pointsize=10)
par(mfrow=c(2, 3))
plot(summ.w, xlab = "time, year", ylab = "Number of individuals",
 legpos = "topright", mfrow = NULL)
plot(summ.w, xlab = "time, year", ylab = "Number of individuals", mfrow = NULL,
 quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, 
"Sensitivity to w (Relapse rate to active TB (per capita))", cex = 1.25)
par(mfrow = c(1, 1))


## Local sensitivity for "c"
tout <- 0:50
print(system.time(
c <- sensRange(func = solveTB, parms = pars, dist = "grid",
sensvar = c("L","Ti","Tn"), parRange = parRanges[10,], num = 1000)
))
head(summary(c))

summ.c <- summary(c)
quartz(width=10, height=6, pointsize=10)
par(mfrow=c(2, 3))
plot(summ.c, xlab = "time, year", ylab = "Number of individuals",
 legpos = "topright", mfrow = NULL)
plot(summ.c, xlab = "time, year", ylab = "Number of individuals", mfrow = NULL,
 quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, 
"Sensitivity to c (Natural cure rate (per capita))", cex = 1.25)
par(mfrow = c(1, 1))

## GLobal sensitivity analysis to Ti

Sens2 <- summary(sensRange(func = solveTB, parms = pars,
 dist = "latin", sensvar = "Ti", parRange = parRanges, num = 1000))
quartz(width=10, height=6, pointsize=10)
plot(Sens2, main = "Sensitivity bpi, beta, mu, v, p, f, q, mut, w, c", xlab = "time, year",
ylab = "Number of individuals")

## GLobal sensitivity analysis to L

Sens2 <- summary(sensRange(func = solveTB, parms = pars,
 dist = "latin", sensvar = "L", parRange = parRanges, num = 1000))
quartz(width=10, height=6, pointsize=10)
plot(Sens2, main = "Sensitivity bpi, beta, mu, v, p, f, q, mut, w, c", xlab = "time, year",
ylab = "Number of individuals")


## Local sensitivity

##Sensitivity functions are generated with sensFun, and estimate the effect of a 
## selection of parameters

## Sensitivity to Ti
SnsTB<- sensFun(func = solveTB, parms = pars,
 sensvar = "Ti", varscale = 1)
head(SnsTB)

quartz(width=10, height=6, pointsize=10)
plot(SnsTB, main="Sensitivity of Ti (Infectious tuberculosis)")

## Sensitivity to L
SnsTB<- sensFun(func = solveTB, parms = pars,
 sensvar = "L", varscale = 1)
head(SnsTB)

quartz(width=10, height=6, pointsize=10)
plot(SnsTB, main="Sensitivity of L (Latently infected)")


## univariate sensitivity
summary(SnsTB)
summary(sensFun(solveTB, pars, varscale = 1), var = TRUE)

## Bivariate sensitivity

cor(SnsTB[ ,-(1:2)])

quartz(width=10, height=6, pointsize=10)
pairs(SnsTB)


## Monte Carlo runs
SF <- function (pars) {
 out <- solveTB(pars)
 return(out[nrow(out), 3:4])
 }
CRL <- modCRL(func = SF, parms = pars, parRange = parRanges[2,])

quartz(width=10, height=6, pointsize=10)
plot(CRL)

#CRL2 <- modCRL(func = SF, parms = pars, parMean = c(beta = 0.5, v = 0.6),
# parCovar = matrix(nr = 2, data = c(0.02, 0.02, 0.02, 0.05)),
# dist = "norm", sensvar = "Ti", num = 150)

#quartz(width=10, height=6, pointsize=10)
#pairs(CRL2)

## Multivariate sensitivity analysis

Coll <- collin(SnsTB)
Coll
Coll [Coll[,"collinearity"] < 20 & Coll[ ,"N"] == 3, ]
collin(SnsTB, parset = 1:2)

### Fitting the model to data

## Model fitting

pars <- list(bpi = 78000, beta = 0.00005,
 mu = 1/40, v = 0.00050, p = 0.05, f=0.70, 
 q=0.85, mut=0.136, w=0.005, c=0.058)

## solving model
solveTB <- function(pars, times=seq(0,20,by=1)) {
 derivs <- function(t, state, pars) { # returns rate of change
 with(as.list(c(state, pars)), {
 dS <- bpi - beta*S*Ti - mu*S
 dL <- (1 - p)*beta*S*Ti - (v + mu)*L
 dTi <- p*f*beta*S*Ti + q*v*L + w*R - (mu + mut + c)*Ti
 dTn <- p*(1 - f)*beta*S*Ti + (1 - q)*v*L + w*R - (mu + mut + c)*Tn
 dR <- c*Ti + c*Tn - (2*w + mu)*R 
 return(list(c(dS, dL, dTi, dTn, dR), TOC = S + L + Ti + Tn + R))
 })
 }
 state <- c(S = 25000000, L = 5000000, Ti=35687, Tn=10000, R=0)
 ## ode solves the model by integration...
 return(ode(y = state, times = times, func = derivs, parms = pars))
 }

out<-solveTB(pars) ## output

setwd("/Users/cesarmunayco/Documents/Tuberculosis/Analisis Epidemiologico de la TB/Spatial and temporal heterogeneity")
getwd()
Data<-read.csv("tb_peru.csv")
Data<-data.frame(time=1:22,R=Data[50:71,3])
Data<-as.matrix(Data)

#Data2newbk<-data.frame(time=1:24,Ti=Data[48:71,4])
#nrow(Data1new)
#nrow(Data2newbk)

## Real data
quartz(width=10, height=6, pointsize=10)
par(bg = "white") ## add background to a plot
plot(Data,pch = 18, cex = 2, xlab="time, year", ylab="Number of TB cases")

## Model fitting

Objective <- function(x, parset = names(x)) {
 pars[parset] <- x
 tout <- seq(0, 22, by = 1)
 ## output times
 out <- solveTB(pars, tout)
 ## Model cost
 return(modCost(obs = Data, model = out))
 }

## Collinearity analysis
Coll <- collin(sF <- sensFun(func = Objective, parms = pars, varscale = 1))
Coll

quartz(width=10, height=6, pointsize=10)
plot(Coll, log = "y")
abline(h = 20, col = "red")

collin(sF,parset=1:3)


## We now use function modFit to locate the minimum

## First fit beta and gamma
print(system.time(Fit <- modFit(p = c(bpi = 78000, beta = 0.000005,
 mu = 1/40, v = 0.00050, p = 0.05, f=0.70, 
 q=0.85, mut=0.136, w=0.005, c=0.058),
 f = Objective, lower = c(50000,0.000001,0.001,0.00010,0.001,0.60,0.50,0.100,0.001,0.001))))
summary(Fit)


init <- solveTB(pars)
plot(init)
pars[c("bpi","beta","mu","v","p","f","q","mut","w","c")] <- Fit$par
out <- solveTB(pars)
Cost <- modCost(obs = Data, model = out)
Cost

quartz(width=10, height=6, pointsize=10)
plot(out, init, xlab = "time, year", ylab = "Number of TB cases", lwd = 2,
 obs = Data, obspar = list(cex = 2, pch = 18))
legend ("bottomright", lwd = 2, col = 1:2, lty = 1:2, c("fitted", "original"))

quartz(width=10, height=6, pointsize=10)
plot(Cost, xlab = "time", ylab = "", main = "residuals")

## Markov chain Monte Carlo

SF<-summary(Fit)
SF
SF[]

Var0 <- SF$modVariance
covIni <- SF$cov.scaled *2.4^2/2

MCMC <- modMCMC(p = coef(Fit), f = Objective, jump = covIni,
 var0 = Var0, wvar0 = 1)

quartz(width=10, height=6, pointsize=10)
plot(MCMC, Full = TRUE) ## The plot method shows the trace of the parameters and, in Full is TRUE, also the model function.

quartz(width=10, height=6, pointsize=10)
pairs(MCMC)

MC <- as.mcmc(MCMC$pars)

quartz(width=10, height=6, pointsize=10)
cumuplot(MC)

cov(MCMC$pars) ## compare the covariances based on generated parameters with the ones from the fit

covIni

## Distributions
quartz(width=10, height=6, pointsize=10)
par(mfrow = c(2, 2))
Minmax <- data.frame(min = c(1, 2), max = c(2, 3))
rownames(Minmax) <- c("par1", "par2")
Mean <- c(par1 = 1.5, par2 = 2.5)
Covar <- matrix(nr = 2, data = c(2, 2, 2, 3))
plot(Unif(Minmax, 100), main = "Unif", xlim = c(1, 2), ylim = c(2, 3))
plot(Grid(Minmax, 100), main = "Grid", xlim = c(1, 2), ylim = c(2, 3))
plot(Latinhyper(Minmax, 5), main = "Latin hypercube", xlim = c(1, 2),
 ylim = c(2, 3))
grid()
plot(Norm(parMean = Mean, parCovar = Covar, num = 1000),
 main = "multi normal")