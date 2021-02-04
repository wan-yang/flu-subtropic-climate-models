# sh1.daily: mean daily specific humidity
library(tgp)
# mean.temp: mean daily temperature
# parms(colnames):'qmin','qmax','qmid','R0max','R0diff','Tc','Tdiff','Texp','D','L','S0','I0','p'
# discrete=T

# calculate R0
calc_R0_AH_T_Strain <- function(in.parms, num.ens, sh1, mean.temp) {
  
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
R0 <- calc_R0_AH_T_Strain(parms, num_ens, sh1.daily, mean.temp.daily)[[2]]

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
# Get seeding times:
load(paste0(dir_code, 'seeding_lists.RData'))
s.list <- seeding.list[[3]]
which.from.slist <- ceiling(tm.range / 7); which.from.slist[7505] <- which.from.slist[7504]
s.list <- s.list[which.from.slist]; rm(which.from.slist)
# the last day is missing b/c no data, set it to the same as previous day
if(is.na(s.list[tm_end])) s.list[tm_end] = s.list[tm_end-1]
SIRS_AH_T_Strain <- function(tm_strt, tm_end, tm_step, S0, I0, N, D, L, beta, expoI, realdata=FALSE){
  cnt=1;
  
  tm_strt=tm_strt-tm.range[1]+1; # adjust the index to match beta
  tm_end=tm_end-tm.range[1]+1;
  tm_vec=seq(tm_strt,tm_end,by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  Np=length(S0); # number of particles
  S=I=newI=array(0,c(Np,tm_sz,3)) # two strains, plus total
  S[, 1, 1:2]=S0; I[, 1, 1:2]=I0;
  newI[, 1, ]=0;
  if(! exists("discrete")) discrete=FALSE; 
  
  if (discrete){
    S[, 1, 1:2]=round(S0,0); I[, 1, 1:2]=round(I0, 0)
    for (t in tm_vec) {
      cnt=cnt+1;
      
      ### START OF PANDEMIC ###
      # only applies to "strain 2" - H1N1/B
      if (t == 4240) { 
        S[, cnt - 1, 2] <- round(as.vector(lhs(num_ens, pdmSinit*N)), 0) 
      }
      #########################
      
      for (j in 1:2) { # for each "strain"
        
        # Step 1:
        Eimmloss <- tm_step * (1 / L) * (N - S[, cnt - 1, j] - I[, cnt - 1, j])
        Einf <- tm_step * (beta[t, ] * pmin(I[, cnt-1, j], I[, cnt-1, j] ^ expoI) * S[, cnt-1, j] / N)
        Erecov <- tm_step * (1 / D) * I[, cnt - 1, j]
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf<0 | is.na(Einf)]=0
        Erecov[Erecov < 0] <- 0
        
        smcl <- rpois(Np, Eimmloss)
        smci <- rpois(Np, Einf)
        smcr <- rpois(Np, Erecov)
        
        sk1 <- smcl - smci
        ik1 <- smci - smcr
        ik1a <- smci
        
        Ts1 <- S[, cnt - 1, j] + round(sk1 / 2, 0)
        Ti1 <- I[, cnt - 1, j] + round(ik1 / 2, 0)
        
        # Step 2:
        Eimmloss <- tm_step * (1 / L) * (N - Ts1 - Ti1)
        Einf <- tm_step * (beta[t, ] * pmin(Ti1, Ti1 ^ expoI) * Ts1 / N)
        Erecov <- tm_step * (1 / D) * Ti1
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf<0 | is.na(Einf)]=0
        Erecov[Erecov < 0] <- 0
        
        smcl <- rpois(Np, Eimmloss)
        smci <- rpois(Np, Einf)
        smcr <- rpois(Np, Erecov)
        
        sk2 <- smcl - smci
        ik2 <- smci - smcr
        ik2a <- smci
        
        Ts2 <- S[, cnt - 1, j] + round(sk2 / 2, 0)
        Ti2 <- I[, cnt - 1, j] + round(ik2 / 2, 0)
        
        # Step 3:
        Eimmloss <- tm_step * (1 / L) * (N - Ts2 - Ti2)
        Einf <- tm_step * (beta[t, ] * pmin(Ti2, Ti2 ^ expoI) * Ts2 / N)
        Erecov <- tm_step * (1 / D) * Ti2
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf<0 | is.na(Einf)]=0
        Erecov[Erecov < 0] <- 0
        
        smcl <- rpois(Np, Eimmloss)
        smci <- rpois(Np, Einf)
        smcr <- rpois(Np, Erecov)
        
        sk3 <- smcl - smci
        ik3 <- smci - smcr
        ik3a <- smci
        
        Ts3 <- S[, cnt - 1, j] + round(sk3 / 2, 0)
        Ti3 <- I[, cnt - 1, j] + round(ik3 / 2, 0)
        
        # Step 4:
        Eimmloss <- tm_step * (1 / L) * (N - Ts3 - Ti3)
        Einf <- tm_step * (beta[t, ] * pmin(Ti3, Ti3 ^ expoI) * Ts3 / N)
        Erecov <- tm_step * (1 / D) * Ti3
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf<0 | is.na(Einf)]=0
        Erecov[Erecov < 0] <- 0
        
        smcl <- rpois(Np, Eimmloss)
        smci <- rpois(Np, Einf)
        smcr <- rpois(Np, Erecov)
        
        sk4 <- smcl - smci
        ik4 <- smci - smcr
        ik4a <- smci
        
        mu = birth.rate.HK 
        
        S[, cnt, j] <- S[, cnt - 1, j] + round(sk1 / 6 + sk2 / 3 + sk3 / 3 + sk4 / 6 + mu * N - mu * S[, cnt - 1, j], 0)
        I[, cnt, j] <- I[, cnt - 1, j] + round(ik1 / 6 + ik2 / 3 + ik3 / 3 + ik4 / 6 - mu * I[, cnt - 1, j], 0)
        newI[, cnt, j] <- newI[, cnt - 1, j] + round(ik1a / 6 + ik2a / 3 + ik3a / 3 + ik4a / 6, 0)
      }
      
      # Random seeding:
      seed <- rpois(Np, 0.1)
      if (s.list[t] == 1) { # seed first "strain"
        S[, cnt, 1] <- S[, cnt, 1] - seed
        I[, cnt, 1] <- I[, cnt, 1] + seed
      } else if (s.list[t] == 2) { # seed second "strain"
        S[, cnt, 2] <- S[, cnt, 2] - seed
        I[, cnt, 2] <- I[, cnt, 2] + seed
      } else if (s.list[t] == 0) { # randomly choose which strain to seed - equally present
        j <- sample(c(1, 2), 1)
        S[, cnt, j] <- S[, cnt, j] - seed
        I[, cnt, j] <- I[, cnt, j] + seed
      }
    }
  } else {
    # run continuously
    for (t in tm_vec) { 
      cnt=cnt+1;
      
      ### START OF PANDEMIC ###
      # only applies to "strain 2" - H1N1/B
      if (t == 4240) { 
        S[, cnt - 1, 2] <- round(as.vector(lhs(num_ens, pdmSinit*N)), 0) 
      }
      #########################
      
      for (j in 1:2) { # for each "strain"
        
        # Step 1:
        Eimmloss <- tm_step * (1 / L) * (N - S[, cnt - 1, j] - I[, cnt - 1, j])
        Einf <- tm_step * (beta[t, ] * I[, cnt - 1, j] * S[, cnt - 1, j] / N)
        Erecov <- tm_step * (1 / D) * I[, cnt - 1, j]
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf < 0] <- 0
        Erecov[Erecov < 0] <- 0
        
        smcl <- Eimmloss
        smci <- Einf
        smcr <- Erecov
        
        sk1 <- smcl - smci
        ik1 <- smci - smcr
        ik1a <- smci
        
        Ts1 <- S[, cnt - 1, j] + sk1 / 2
        Ti1 <- I[, cnt - 1, j] + ik1 / 2
        
        # Step 2:
        Eimmloss <- tm_step * (1 / L) * (N - Ts1 - Ti1)
        Einf <- tm_step * (beta[t, ] * Ti1 * Ts1 / N)
        Erecov <- tm_step * (1 / D) * Ti1
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf < 0] <- 0
        Erecov[Erecov < 0] <- 0
        
        smcl <- Eimmloss
        smci <- Einf
        smcr <- Erecov
        
        sk2 <- smcl - smci
        ik2 <- smci - smcr
        ik2a <- smci
        
        Ts2 <- S[, cnt - 1, j] + sk2 / 2
        Ti2 <- I[, cnt - 1, j] + ik2 / 2
        
        # Step 3:
        Eimmloss <- tm_step * (1 / L) * (N - Ts2 - Ti2)
        Einf <- tm_step * (beta[t, ] * Ti2 * Ts2 / N)
        Erecov <- tm_step * (1 / D) * Ti2
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf < 0] <- 0
        Erecov[Erecov < 0] <- 0
        
        smcl <- Eimmloss
        smci <- Einf
        smcr <- Erecov
        
        sk3 <- smcl - smci
        ik3 <- smci - smcr
        ik3a <- smci
        
        Ts3 <- S[, cnt - 1, j] + sk3 / 2
        Ti3 <- I[, cnt - 1, j] + ik3 / 2
        
        # Step 4:
        Eimmloss <- tm_step * (1 / L) * (N - Ts3 - Ti3)
        Einf <- tm_step * (beta[t, ] * Ti3 * Ts3 / N)
        Erecov <- tm_step * (1 / D) * Ti3
        
        Eimmloss[Eimmloss < 0] <- 0
        Einf[Einf < 0] <- 0
        Erecov[Erecov < 0] <- 0
        
        smcl <- Eimmloss
        smci <- Einf
        smcr <- Erecov
        
        sk4 <- smcl - smci
        ik4 <- smci - smcr
        ik4a <- smci
        
        mu = birth.rate.HK 
        
        S[, cnt, j] <- S[, cnt - 1, j] + sk1 / 6 + sk2 / 3 + sk3 / 3 + sk4 / 6 + mu * N - mu * S[, cnt - 1, j]
        I[, cnt, j] <- I[, cnt - 1, j] + ik1 / 6 + ik2 / 3 + ik3 / 3 + ik4 / 6 - mu * I[, cnt - 1, j]
        newI[, cnt, j] <- newI[, cnt - 1, j] + ik1a / 6 + ik2a / 3 + ik3a / 3 + ik4a / 6
        
        I[, cnt, j][I[, cnt, j] < 1] <- 0 
        newI[, cnt, j][I[, cnt, j] == 0] <- newI[, cnt - 1, j][I[, cnt, j] == 0]
      }
      
      # Random seeding:
      seed <- rpois(Np, 0.1)
      if (s.list[t] == 1) { # seed first "strain"
        S[, cnt, 1] <- S[, cnt, 1] - seed
        I[, cnt, 1] <- I[, cnt, 1] + seed
      } else if (s.list[t] == 2) { # seed second "strain"
        S[, cnt, 2] <- S[, cnt, 2] - seed
        I[, cnt, 2] <- I[, cnt, 2] + seed
      } else if (s.list[t] == 0) { # randomly choose which strain to seed - equally present
        j <- sample(c(1, 2), 1)
        S[, cnt, j] <- S[, cnt, j] - seed
        I[, cnt, j] <- I[, cnt, j] + seed
      }
      
    }
  }
  
  # Calculate total influenza:
  S[,, 3] <- S[,, 1] + S[,, 2]
  I[,, 3] <- I[,, 1] + I[,, 2]
  newI[,, 3] <- newI[,, 1] + newI[,, 2]
  
  S1 <- S[,, 1]; I1 <- I[,, 1]; newI1 <- newI[,, 1]
  S2 <- S[,, 2]; I2 <- I[,, 2]; newI2 <- newI[,, 2]
  S <- S[,, 3]; I <- I[,, 3]; newI <- newI[,, 3]
  
  S1=t(S1); I1=t(I1); newI1=t(newI1);
  S2=t(S2); I2=t(I2); newI2=t(newI2);
  S=t(S); I=t(I); newI=t(newI);
  if (realdata==FALSE){
    rec=list(S=S,I=I); 
  } else {
    rec=list(S=S,I=I,newI=newI);
  }
  rec;
}
sim_AH_T_Strain <- SIRS_AH_T_Strain(tm_strt, tm_end, tm_step, S0, I0, N, D, L, beta.red, expoI, realdata = T)
