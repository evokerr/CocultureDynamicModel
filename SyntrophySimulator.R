####################################################################
# These are the simultaions for a DM mutualism model 
####################################################################

install.packages("deSolve")
library(deSolve)

#---------------------
# PARAMETERS:
#---------------------

# Equation parameters:
mu.D.1 <- .06
mu.D.2 <- .07
mu.M.1 <- .05
mu.M.2 <- .05
k.D.1 <- 1.5
k.D.2 <- 1.5
k.M.1 <- .00005
k.M.2 <- .00005
delta.D.1 <- 0
delta.D.2 <- 0
delta.M.1 <- 0
delta.M.2 <- 0
y.D.1 <- 2.5*(10^6)
y.D.2 <- 2.5*(10^5)
y.M.1 <- 2.62*(10^6)
y.M.2 <- 2.62*(10^6)
alpha.D.1 <- 1.906
alpha.D.2 <- 1.906
P.star.1 <- .75
P.star.2 <- .75

# Parameter vector:
parameters <- c(mu.D.1, mu.D.2, mu.M.1, mu.M.2,
                k.D.1, k.D.2, k.M.1, k.M.2,
                delta.D.1, delta.D.2, delta.M.1, delta.M.2,
                y.D.1, y.D.2, y.M.1, y.M.2, alpha.D.1, alpha.D.2,
                P.star.1, P.star.2)

# Initial densities:
D.1_0 <- 0  
D.2_0 <- 10^6
M.1_0 <- .9*10^6 
M.2_0 <- 0
S_0 <- 40
P_0 <- 0

#----------------------------
#FUNCTIONS & STATE VARIABLES:
#----------------------------

#This function returns the growth rate given maximal growth rate,
#half saturation, and the amount of resource according a basic
#Monod function
gamma <- function(mu, k, R) {
  return((mu*R)/(k+R))
}

#This vector defines the state variables 
state <- c(D1 = D.1_0,
           D2 = D.2_0,
           M1 = M.1_0,
           M2 = M.2_0,
           S = S_0,
           P = P_0)

#This function lays out the differential equations that
#define the DM model.
DM_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dD1 <- (1-(P/P.star.1))*gamma(mu.D.1,k.D.1,S)*D1 - delta.D.1*D1
    dD2 <- (1-(P/P.star.2))*gamma(mu.D.2,k.D.2,S)*D2 - delta.D.2*D2
    dM1 <- gamma(mu.M.1,k.M.1,P)*M1 - delta.M.1*M1
    dM2 <- gamma(mu.M.2,k.M.2,P)*M2 - delta.M.2*M2
    dS <- (-(1-(P/P.star.1))*gamma(mu.D.1,k.D.1,S)*D1*(1/y.D.1) - (1-(P/P.star.2))*gamma(mu.D.2,k.D.2,S)*D2*(1/y.D.2))
    dP <- (1-(P/P.star.1))*gamma(mu.D.1,k.D.1,S)*D1*(1/y.D.1)*alpha.D.1 + (1-(P/P.star.2))*gamma(mu.D.2,k.D.2,S)*D2*(1/y.D.2)*alpha.D.2 - gamma(mu.M.1,k.M.1,P)*M1*(1/y.M.1) - gamma(mu.M.2,k.M.2,P)*M2*(1/y.M.2)
    list(c(dD1, dD2, dM1, dM2, dS, dP))
  }) 
}


#Using the deSolve function "ode", this function numerically simulates
#one of the parametric cases for a certain amount of time. It returns
#a data frame with all the time courses for each host type.
simulate <- function(total.time, time.increment) {
  times <- seq(0, total.time, by = time.increment)
  out <- ode(y = state, times = times, func = DM_model, parms = parameters)
  colnames(out) <- c("time","Dv_1","Dv_2","Mm_1","Mm_2","S","P")
  mod.df<-as.data.frame(out)
  return(mod.df)
}

#NEW PARAMETERS
N.demes <- 200
Deme.vol <- 1
Pool.dil <- 1/50
time.per.round <- 24
#N.rounds<-50

Cumulative.N.strains <- c(D.1_0*N.demes,
                          D.2_0*N.demes,
                          M.1_0*N.demes,
                          M.2_0*N.demes)
Cumulative.vol <- N.demes * Deme.vol
simulate.metapopulation <- function(N.rounds) {
  D.1 <- rep(0,N.rounds)
  D.2 <- rep(0,N.rounds)
  M.1 <- rep(0,N.rounds)
  M.2 <- rep(0,N.rounds)
  traj.df <- data.frame(D.1,D.2,M.1,M.2)
  for(r in 1:N.rounds) { 
    cell.count <- sum(Cumulative.N.strains)
    mean.cell.number <- (cell.count/Cumulative.vol)*Deme.vol
    init.cells <- rpois(N.demes, mean.cell.number)
    prob.strains <- Cumulative.N.strains/cell.count
    Cumulative.N.strains <- c(0,0,0,0)
    for(deme in 1:N.demes) { #deme<-1
      init.strains <- rmultinom(1,init.cells[deme],prob.strains)
      state <- c(D1 = init.strains[1,1]/Deme.vol,
                 D2 = init.strains[2,1]/Deme.vol,
                 M1 = init.strains[3,1]/Deme.vol,
                 M2 = init.strains[4,1]/Deme.vol,
                 S = S_0,
                 P = P_0)
      print(init.strains[,1])
      times <- seq(0, time.per.round, by = time.per.round)
      out <- ode(y = state, times = times, func = DM_model, parms = parameters)
      colnames(out) <- c("time","Dv_1","Dv_2","Mm_1","Mm_2","S","P")
      mod.df<-as.data.frame(out)
      Cumulative.N.strains <- Cumulative.N.strains + Deme.vol*(c(max(0.0,mod.df[mod.df$time==time.per.round,]$Dv_1),
                                                                 max(0.0,mod.df[mod.df$time==time.per.round,]$Dv_2),
                                                                 max(0.0,mod.df[mod.df$time==time.per.round,]$Mm_1),
                                                                 max(0.0,mod.df[mod.df$time==time.per.round,]$Mm_2)))
    }
    traj.df[r,1] <- Cumulative.N.strains[1]/Cumulative.vol
    traj.df[r,2] <- Cumulative.N.strains[2]/Cumulative.vol
    traj.df[r,3] <- Cumulative.N.strains[3]/Cumulative.vol
    traj.df[r,4] <- Cumulative.N.strains[4]/Cumulative.vol
    Cumulative.N.strains <- Cumulative.N.strains*Pool.dil
  }
  return(traj.df)
}

#This function plots (or adds a line to an existing plot, if first==FALSE)
#the trajectory defined by vectors x and y. A color vector (rgbv) as well
#as labels and limits for the y axis are also arguments.
add.curve <- function(x, y, is.first, rgbv, ylabel, ylimit, is.shaded) {
  if(is.first) {
    par(mar = c(4, 4, .5, .5),mfrow=c(1,1))
    plot(x, y, type="l", lwd=2, xlab = "time", ylab = ylabel, 
         col=rgb(rgbv[1],rgbv[2],rgbv[3]),ylim=c(ylimit[1],ylimit[2]))
    if(is.shaded) {
      polygon(c(x,rev(x)),c(y,rep(0,length(x))),
              col=rgb(rgbv[1],rgbv[2],rgbv[3],0.5),border=NA)
    }
    lines(x, y, lwd=2, col=rgb(rgbv[1],rgbv[2],rgbv[3]))
  } else {
    if(is.shaded) {
      polygon(c(x,rev(x)),c(y,rep(0,length(x))),
            col=rgb(rgbv[1],rgbv[2],rgbv[3],0.5),border=NA)
    }
    lines(x, y, lwd=2, col=rgb(rgbv[1],rgbv[2],rgbv[3]))
  }
}


#-------------------------------
#RUNNING & PLOTTING SIMULATIONS:
#-------------------------------

#Plotting parameters:
plot.upper.limit <- 9

mod[5600,]


#Draw dynamics
mod<-simulate(120,1/60)
add.curve(mod$time, log10(mod$Dv_1+1), TRUE, c(1,0,0), "log(density)", c(0,plot.upper.limit), FALSE)
add.curve(mod$time, log10(mod$Dv_2+1), FALSE, c(.75,0,0), "density", c(0,plot.upper.limit), FALSE)
add.curve(mod$time, log10(mod$Mm_1+1), FALSE, c(0,0,1), "density", c(0,plot.upper.limit), FALSE)
add.curve(mod$time, log10(mod$P+1), FALSE, c(0,1,0), "density", c(0,plot.upper.limit), FALSE)

#add.curve(mod$time, mod$Dv_2, FALSE, c(.75,0,0), "density", c(0,plot.upper.limit), FALSE)
#add.curve(mod$time, mod$Mm_1, FALSE, c(0,0,1), "density", c(0,plot.upper.limit), FALSE)
#add.curve(mod$time, mod$P, FALSE, c(0,1,0), "density", c(0,plot.upper.limit), FALSE)

#add.curve(mod$time, mod$S, FALSE, c(.5,.5,.5), "density", c(0,plot.upper.limit), FALSE)
#add.curve(mod$time, mod$P, FALSE, c(0,1,0), "density", c(0,plot.upper.limit), FALSE)
#lines(c(0,24),c(700,700),lty=2,col="gray")
#add.curve(mod$time, mod$Mm_2, FALSE, c(0,0,.75), "density", c(0,plot.upper.limit), FALSE)

#met.df<-simulate.metapopulation(30)
#add.curve(1:30, met.df$D.1, TRUE, c(1,0,0), "density", c(0,plot.upper.limit), FALSE)
#add.curve(1:30, met.df$D.2, FALSE, c(.75,0,0), "density", c(0,plot.upper.limit), FALSE)
#add.curve(1:30, met.df$M.1, FALSE, c(0,0,1), "density", c(0,plot.upper.limit), FALSE)

#met2.df<-simulate.metapopulation(30)
#add.curve(1:30, met2.df$D.1, TRUE, c(1,0,0), "density", c(0,plot.upper.limit), FALSE)
#add.curve(1:30, met2.df$D.2, FALSE, c(.75,0,0), "density", c(0,plot.upper.limit), FALSE)
#add.curve(1:30, met2.df$M.1, FALSE, c(0,0,1), "density", c(0,plot.upper.limit), FALSE)

