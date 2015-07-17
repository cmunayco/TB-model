library(sensitivity)
library(rgl)
library(Rglpk)
library(pse)
library(FME)
require(ggplot2)
require(reshape)


## Real data
Data<-read.csv("tb_peru.csv")
Data1<-data.frame(time=1:23,R=Data[48:70,6])
Data1<-as.matrix(Data1)



quartz(width=10, height=6, pointsize=10)
par(mar=c(5,4,4,2)+0.1)
#par(bg = "white") ## add background to a plot
plot(Data1[,1],Data1[,2],pch = 18, cex = 2, xlab="time, year", ylab="Number of TB cases")



##### The TB model
sibr.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  L <- x[2]
  Ti <- x[3]
  Tn <- x[4]
  R <- x[5]
  ## now extract the parameters
  br <- params["br"]
  beta <- params["beta"]
  mug <- params["mug"]
  v <- params["v"]
  pp <- params["pp"]
  f <- params["f"]
  q <- params["q"]
  mut <- params["mut"]
  w <- params["w"]
  c <- params["c"]
  N <- 9931000;
  
  ## now code the model equations
  dS.dt <- br - beta*S*Ti - mug*S
  dL.dt <- (1 - pp)*beta*S*Ti - (v + mug)*L
  dTi.dt <- pp*f*beta*S*Ti + q*v*L + w*R - (mug + mut + c)*Ti
  dTn.dt <- pp*(1 - f)*beta*S*Ti + (1 - q)*v*L + w*R - (mug + mut + c)*Tn
  dR.dt <- c*Ti + c*Tn - (2*w + mug)*R 
  ## combine results into a single vector
  dxdt <- c(dS.dt,dL.dt,dTi.dt,dTn.dt,dR.dt)
  ## return result as a list!
  list(dxdt)
}

params <- c(br = 278068, beta = 0.000001,
            mug = 1/40, v = 0.00256, pp = 0.05, f=0.50, 
            q=0.5, mut=0.136, w=0.005, c=0.058)

times <- seq(from=0,to=22,by=1/4) ## returns a sequence
xstart <- c(S =6907085 , L = 3000000, Ti=21525, Tn=2390, R=24000) ## initial conditions



out <- ode(
  func=sibr.model,
  y=xstart,
  times=times,
  parms=params
)
class(out)

head(out)

out <- as.data.frame(
  ode(
    func=sibr.model,
    y=xstart,
    times=times,
    parms=params
  )
)


quartz(width=10, height=6, pointsize=10)
out <- subset(out,select=c(S,L,Ti,Tn,R,time))
ggplot(data=melt(out,id.var="time"),
       mapping=aes(x=time,y=value,group=variable,color=variable))+
  geom_line()


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
matplot(out[,1], (out[,c(4)]+ out[,c(5)]), type = "l", lty = 1:2, lwd = c(2,1),
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


## Real data
quartz(width=10, height=6, pointsize=10)
par(mar=c(5,4,4,2)+0.1)
#par(bg = "white") ## add background to a plot
plot(1:23,Data1[,2],pch = 18, cex = 2, xlab="time, year", ylab="Number of TB cases",
     ylim=c(0,250000))
par(new=T)
plot(1:23, (out[c(3:25),c(4)]+ out[c(3:25),c(5)]),xlab="",ylab="", col="red",
     ylim=c(0,250000))


##### Fitting the model

times <- seq(from=0,to=23,by=1)
params <- c(br = 278068, beta = 0.000001,
            mug = 1/40, v = 0.00256, pp = 0.05, f=0.50, 
            q=0.5, mut=0.136, w=0.005, c=0.058,N=9931000,p=0.3,k=1000)
out <- as.data.frame(
  ode(
    func=sibr.model,
    y=xstart,
    times=times,
    parms=params
  )
)
within(
  out,
  C <- rnbinom(n=length(R),mu=params["p"]*R,size=params["k"])
) -> out
quartz(width=10, height=6, pointsize=10)
ggplot(data=out,mapping=aes(x=time,y=C))+geom_point()+geom_line()



sibr.nll <- function (br, beta, mug, v, pp, f, q, mut, w, c, I.0, p, k) {
  times <- c(Data1[,1]-1,Data1[,1])
  ode.params <- c(br=br,beta=beta, mug=mug, v=v,p=p,f=f,q=q,mut=mut,w=w,c=c)
  xstart <- c(S=6907085-I.0,L = 3000000, I=I.0,Tn=2390,R=24000)
  out <- ode(
    func=sibr.model,
    y=xstart,
    times=times,
    parms=ode.params
  )
  ## 'out' is a matrix
  ll <- dnbinom(x=Data1[,2],size=params["k"],mu=p*out[-1,"R"],log=TRUE)
  -sum(ll)
}


nll <- function (par) {
  sibr.nll(beta=par[1],br=278068,mug = 1/40, v = 0.00256, pp = 0.05, f=0.50, 
           q=0.5, mut=0.136, w=0.005, c=0.058,
           I.0=21525,p=0.5,k=100)
}
betacurve <- data.frame(beta=seq(0.0000001,0.00001,length=100))
within(betacurve,nll <- sapply(beta,nll)) -> betacurve
ggplot(data=betacurve,mapping=aes(x=beta,y=nll))+geom_line()


fit <- optim(fn=nll,par=0.0000001,method="Brent",lower=0.00000001,upper=0.000001)
fit

write.csv(Data1,"Data1.csv")
