library(tgp)
# sh1.daily: mean daily specific humidity
# mean.temp: mean daily temperature
# parms(colnames):'qmin','qmax','qmid','R0max','R0diff','Tc','Tdiff','Texp','D','Lshort','L','rho','S0','I0','p'
# discrete=T

# calculate R0
calc_R0_AH_T_Vary <- function(in.parms, num.ens, sh1, mean.temp) {
  
  # Create matrices for storing/returning results:
  res.temp = res.temp.red = matrix(0, length(sh1), num.ens)
  
  # Loop through all ensemble members:
  for (ix in 1:num.ens) {
    
    # Assign parameters:
    q.mn <- in.parms[ix, 'qmin']; q.mx <- in.parms[ix, 'qmax']; q.md <- in.parms[ix, 'qmid']
    R0.max <- in.parms[ix, 'R0max']; R0.diff <- in.parms[ix, 'R0diff']
    Tc <- in.parms[ix, 'Tc']; Tdiff <- in.parms[ix, "Tdiff"]; t.exp <- in.parms[ix, 'Texp']
    
    if(F){
      if (dim(in.parms)[2] == 9) {
        q.mn.cut <- in.parms[ix, 9]
      } else {
        q.mn.cut <- q.mn
      }
    }
    
    q.mn.cut <- q.mn
    
    # Calculate and correct R0.min
    R0.min <- R0.max - R0.diff
    if (R0.min < 0) {
      R0.min <- 0.1
    }
    Tmin <- Tc - Tdiff
    
    # Calculate parabola params:
    if(F){
      b <- ((R0.max - R0.min) * (q.mx + q.mn)) / ((q.mx - q.md) * (q.mn - q.md))
      a <- (-1 * b) / (q.mx + q.mn)
      c <- R0.min - a * q.md ** 2 - b * q.md
    }
    
    # given the symmetry:
    q.mx.left = 2 * q.md - q.mn; 
    b.left <- ((R0.max - R0.min) * (q.mx.left + q.mn)) / ((q.mx.left - q.md) * (q.mn - q.md))
    a.left <- (-1 * b.left) / (q.mx.left + q.mn)
    c.left <- R0.min - a.left * q.md ** 2 - b.left * q.md
    
    q.mn.right = 2 * q.md - q.mx
    b.right <- ((R0.max - R0.min) * (q.mx + q.mn.right)) / ((q.mx - q.md) * (q.mn.right - q.md))
    a.right <- (-1 * b.right) / (q.mx + q.mn.right)
    c.right <- R0.min - a.right * q.md ** 2 - b.right * q.md
    
    fit1 = fit2 =numeric(length(sh1))
    
    # split the data into two sets (those >=q.md, and those <q.md)
    idx.left = which(sh1 < q.md); idx.right = which(sh1 >= q.md)
    
    # Full model:
    q1.left <- sh1[idx.left]; q1.right = sh1[idx.right]
    t1.left <- mean.temp[idx.left]; t1.right = mean.temp[idx.right]
    fit1[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
    fit1[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
    
    # Reduced model:
    q1 <- sh1
    q1[q1 < q.mn.cut] <- q.mn.cut; q1[q1 > q.mx] <- q.mx
    t1 <- mean.temp; t1[t1 < Tmin] <- Tmin
    q1.left <- q1[idx.left]; q1.right = q1[idx.right]
    t1.left <- t1[idx.left]; t1.right = t1[idx.right]
    fit2[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
    fit2[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
    
    
    # Store results in matrices:
    res.temp[, ix] <- fit1; res.temp.red[, ix] <- fit2
  }
  
  # Return results:
  return(list(res.temp, res.temp.red))
}
num_ens=dim(parms)[1]
R0 <- calc_R0_AH_T_Vary(parms, num_ens, sh1.daily, mean.temp.daily)[[2]]

# Get other params:
D <- parms[, 'D']; S0 <- parms[, 'S0']; I0 <- parms[, 'I0']; expoI <- parms[, 'p']; 
L <- parms[,"L"]; Lshort <- parms[,'Lshort']; Llong <- parms[,'L']; percLshort <- parms[,'rho']

tm_strt <- 1; tm_end <- length(sh1.daily) 
tm_step <- 1; tm.range <- tm_strt:tm_end

# Calculate betas:
beta <- sapply(1:num_ens, function(ix) {
  R0[, ix] / D[ix]
})
rm(R0)

# SIRS model
SIRS_AH_T_Vary <- function(tm_strt, tm_end, tm_step, S0, I0, N, D, Lshort,Llong,percLshort, beta,birth.rate = birth.rate.HK, expoI=1, expoS=1, newIpad=NULL, realdata=FALSE){
  # run stochastically, using random Poisson draws
  # immunity lose modeled as a step function, rather than gradully
  # function to integrate to the next time step
  # use SIR model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         D: infection period, day; L: immune period, day; 
  #         alpha: rate from exposed to infectious; beta: transmission matrix
  # output: S, I for all time steps
  cnt=1;
  # beta stores only data during the time used for the truth
  tm_vec=seq(tm_strt,tm_end,by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  Np=length(S0); # number of particles
  S=I=newI=newi=matrix(0,Np,tm_sz)
  S[,1]=S0; I[,1]=I0; 
  newI[,1]=0; newi[,1]=0
  if(! is.null(newIpad) & is.null(dim(newIpad))) newIpad=matrix(newIpad,Np,length(newIpad),byrow=T) # if pad newI at the very begining
  # run continuously
  for (t in tm_vec){
    cnt=cnt+1;
    
    # toggle pandemic start:
    if (t == 4240) { 
      S[, cnt - 1] <- round(as.vector(lhs(num_ens, pdmSinit*N)), 0) 
    }
    #########################
    
    Eimmloss=tm_step*(1/Llong*(N-S[,cnt-1]-I[,cnt-1]) ) * (1-percLshort) # long term immunity * (1-percLshort)
    Einf=tm_step*(beta[t,]* pmin(I[,cnt-1], I[,cnt-1]^expoI) * S[,cnt-1]^expoS /N) # raise I to expoI - nonlinear
    Erecov=tm_step*(1/D*I[,cnt-1])
    Eimmloss[Eimmloss<0]=0   
    Einf[Einf<0 | is.na(Einf)]=0
    Erecov[Erecov<0]=0
    smcl=rpois(Np,Eimmloss)
    smci=rpois(Np,Einf)
    smcr=rpois(Np,Erecov)
    
    sk1=smcl-smci
    ik1=smci-smcr
    ik1a=smci;
    Ts1=S[,cnt-1]+sk1/2
    Ti1=I[,cnt-1]+ik1/2
    
    Eimmloss=tm_step*(1/Llong*(N-Ts1-Ti1)) * (1-percLshort)
    Einf=tm_step*(beta[t,]*pmin(Ti1, Ti1^expoI)*Ts1^expoS/N)
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
    Ts2=S[,cnt-1]+sk2/2
    Ti2=I[,cnt-1]+ik2/2
    
    Eimmloss=tm_step*(1/Llong*(N-Ts2-Ti2)) * (1-percLshort)
    Einf=tm_step*(beta[t,]*pmin(Ti2,Ti2^expoI)*Ts2^expoS/N)
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
    Ts3=S[,cnt-1]+sk3
    Ti3=I[,cnt-1]+ik3
    
    Eimmloss=tm_step*(1/Llong*(N-Ts3-Ti3)) * (1-percLshort)
    Einf=tm_step*(beta[t,]*pmin(Ti3,Ti3^expoI)*Ts3^expoS/N)
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
    
    seed=rpois(Np,.1); 
    S[,cnt]=round(S[,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed + birth.rate*(N-S[,cnt-1]),0)
    I[,cnt]=round(I[,cnt-1]+ik1/6+ik2/3+ik3/3+ik4/6+seed - birth.rate*I[,cnt-1],0) # death
    # add immune loss
    Tidx=pmax(0,cnt-round(Lshort,0))
    J=1:Np
    Jnon0=J[Tidx!=0]
    if(any(Jnon0)){
      for(jj in Jnon0){
        S[jj,cnt]=S[jj,cnt] + newi[jj,Tidx[jj]] * percLshort[jj]
      }
    }
    if(!is.null(newIpad)){
      Tidx=cnt-round(Lshort,0)
      J=1:Np
      Jneg=J[Tidx<=0]
      if(any(Jneg)){
        for(jj in Jneg){
          S[jj,cnt]=S[jj,cnt] + newIpad[jj,Tidx[jj]+365] * percLshort[jj]
        }
      }
    }
    S[,cnt]=pmin(S[,cnt],N)
    newI[,cnt]=round(newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0);
    newi[,cnt]=round(ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0); # new cases on that day
  }
  
  S=t(S); I=t(I); newI=t(newI);
  if (realdata==FALSE){
    rec=list(S=S,I=I); 
  } else {
    rec=list(S=S,I=I,newI=newI); 
  }
  rec;
}
sim_AH_T_Vary <- SIRS_AH_T_Vary(tm_strt, tm_end, tm_step, S0, I0, N, D, Lshort,Llong,percLshort, beta,birth.rate = birth.rate.HK, expoI=1, expoS=1,realdata=T)
