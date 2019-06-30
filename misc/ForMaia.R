R0 <- 100
S <- exp(-0.2)
Amax <- 7
Neqn <- rep(0,Amax)
Matt <- c(0,0.2,0.5,0.8,1,1,1)  ##maturity at time
# Dummy - I am lazy (this is expected weight)
Wght <- c(0,1,1,2,2,2,2)  

Neqn[1] <- R0
for (Iage in 2:Amax) Neqn[Iage] <- Neqn[Iage-1]*S ## survivorship
Neqn[Amax] <- Neqn[Amax]/(1-S) ## force plus group
print(Neqn)

SSB0 <- sum(Neqn*Matt*Wght) ## this is what we did too
# SSB
print(sum(Neqn))
# Nmature (for checking)
print(sum(Neqn*Matt))
NCnt <-20*sum(Neqn) ## never used again
NIBM <- matrix(NA,nrow=20*sum(Neqn),ncol=4) ## empty matrix

Ipnt <- 0
for (Iage in 1:Amax)
 {
  # Integer numbers
  Nage <- round(Neqn[Iage]) 
  for (Inum in 1:Nage)
   {
    Ipnt <- Ipnt + 1
    NIBM[Ipnt,1] <- Iage-1 ## age
    # Length (initial) - don't worry this will burn out quickly
    NIBM[Ipnt,2] <- 100*(1.0-exp(-0.2*(Iage-1)))*exp(rnorm(1,0,0.2)-0.2^2/2.0) ## bias corr length
    # Is this animal mature
    NIBM[Ipnt,3] <- runif(1,0,1) < Matt[Iage]   ## logistic ogive                       
    # Weight
    NIBM[Ipnt,4] <- 0.0001*(NIBM[Ipnt,2])^3
   }  
 }  
print(head(NIBM))
# How many animals
print(Ipnt)
# How many mature animals
print(sum(NIBM[1:Ipnt,3]))

# COmpute the Ps
P <- rep(1,Amax); P[1] <- 0
for (Iage in 2:Amax) 
 if (Matt[Iage-1] < 1)  
  P[Iage] <- (Matt[Iage]-Matt[Iage-1])/(1-Matt[Iage-1])
print(P)

# Now loop forward
Nyear <- 123
for (Iyear in 1:Nyear)
 {
  Alive <- !is.na(NIBM[,1]) & NIBM[,1] > -1
 cat("Alive",sum(Alive),"\n")
  
  # Compute SSB
  SSB <- sum(NIBM[Alive,3]*NIBM[Alive,4])
  # COmpute recruits here
  for (II in 1:Ipnt)
   if (Alive[II])  
    {
     MyAge <- NIBM[II,1]
     # Died this year
     if (ruinf(1,0,1) > S) NIBM[II,1] <- -1
     # Now grow the animal and update its weight
     
     # Mature this year
     if (NIBM[II,3] == 0)
       if (runif(1,0,1)<P[MyAge]) NIBM[II,3] <- 1
     
    } 
  
  
  # Add recruits starting from slots Ipnt + 1
  
  
  
  AAA
  # TO speed up remove all Age -1 animals at each time-step
  
  
 }  


   