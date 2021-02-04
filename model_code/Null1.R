library(tgp)
# sh1.daily: mean daily specific humidity
# mean.temp: mean daily temperature
# parms(colnames):'R0','D','L','S0','I0','p'
# discrete=T

# calculate R0
num_ens=dim(parms)[1]
R0 <- matrix(parms[, 1], length(sh1.daily), num_ens, byrow = TRUE)

# Get other params:
D <- parms[, 'D']; S0 <- parms[, 'S0']; I0 <- parms[, 'I0']; expoI <- parms[, 'p']; L <- parms[,"L"]

tm_strt <- 1; tm_end <- length(sh1.daily) 
tm_step <- 1; tm.range <- tm_strt:tm_end

# Calculate betas:
beta <- sapply(1:num_ens, function(ix) {
  R0[, ix] / D[ix]
})
rm(R0)

# SIRS model
SIRS_Null1 <- function(tm_strt, tm_end, tm_step, S0, I0, N, D, L, beta, expoI, realdata=FALSE){
  # function to integrate to the next time step
  # use SIR model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         D: infection period, day; L: immune period, day; 
  #         alpha: rate from exposed to infectious; beta: transmission matrix
  # output: S, I for all time steps
  cnt=1;
  # beta stores only data during the time used for the truth
  tm_strt=tm_strt-tm.range[1]+1; # adjust the index to match beta
  tm_end=tm_end-tm.range[1]+1;
  tm_vec=seq(tm_strt,tm_end,by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  Np=length(S0); # number of particles
  S=I=newI=matrix(0,Np,tm_sz)
  S[,1]=S0; I[,1]=I0; 
  newI[,1]=0;
  if(! exists("discrete")) discrete=FALSE; 
  
  print(discrete)
  
  if (discrete){
    S[,1]=round(S0,0); I[,1]=round(I0,0); 
    for (t in tm_vec){
      cnt=cnt+1;
      
      # toggle pandemic start:
      if (t == 4240) { 
        S[, cnt - 1] <- round(as.vector(lhs(num_ens, pdmSinit*N)), 0) 
      }
      #########################
      
      Eimmloss=tm_step*(1/L*(N-S[,cnt-1]-I[,cnt-1]))
      Einf=tm_step*(beta[t,]*pmin(I[,cnt-1], I[,cnt-1]^expoI)*S[,cnt-1]/N)
      Erecov=tm_step*(1/D*I[,cnt-1])
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci
      Ts1=S[,cnt-1]+round(sk1/2,0)
      Ti1=I[,cnt-1]+round(ik1/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]*pmin(Ti1, Ti1^expoI)*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0  
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[,cnt-1]+round(sk2/2,0)
      Ti2=I[,cnt-1]+round(ik2/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]*pmin(Ti2, Ti2^expoI)*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0   
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0                   
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[,cnt-1]+round(sk3,0)
      Ti3=I[,cnt-1]+round(ik3,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]*pmin(Ti3, Ti3^expoI)*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0   
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      seed=rpois(Np,.1)
      mu = birth.rate.HK 
      S[,cnt]=S[,cnt-1]+round(sk1/6+sk2/3+sk3/3+sk4/6 + mu * N - mu * S[, cnt-1],0)-seed
      I[,cnt]=I[,cnt-1]+round(ik1/6+ik2/3+ik3/3+ik4/6 - mu * I[, cnt-1],0)+seed
      newI[,cnt]=round(newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0);
    }
  } else {
    # run continuously
    for (t in tm_vec){
      cnt=cnt+1;
      
      # toggle pandemic start: 
      if (t == 4240) {
        S[, cnt - 1] <- as.vector(lhs(num_ens, c(0.6 * N, 0.8 * N)))
      }
      #########################
      
      Eimmloss=tm_step*(1/L*(N-S[,cnt-1]-I[,cnt-1]))
      Einf=tm_step*(beta[t,]*I[,cnt-1]*S[,cnt-1]/N)
      Erecov=tm_step*(1/D*I[,cnt-1])
      Eimmloss[Eimmloss<0]=0   
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci
      Ts1=S[,cnt-1]+sk1/2
      Ti1=I[,cnt-1]+ik1/2
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]*Ti1*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0  
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[,cnt-1]+sk2/2
      Ti2=I[,cnt-1]+ik2/2
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]*Ti2*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0   
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0                   
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[,cnt-1]+sk3
      Ti3=I[,cnt-1]+ik3
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]*Ti3*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0  
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      seed = 1 
      
      mu = birth.rate.HK 
      
      S[,cnt]=S[,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed + mu * N - mu * S[,cnt-1]
      I[,cnt]=I[,cnt-1]+ik1/6+ik2/3+ik3/3+ik4/6+seed - mu * I[,cnt-1] # natural mortality
      newI[,cnt]=newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed;
    }
  }
  
  S=t(S); I=t(I); newI=t(newI);
  if (realdata==FALSE){
    rec=list(S=S,I=I); 
  } else {
    rec=list(S=S,I=I,newI=newI); 
  }
  rec;
}
sim_Null1 <- SIRS_Null1(tm_strt, tm_end, tm_step, S0, I0, N, D, L, beta, expoI, realdata=T)
