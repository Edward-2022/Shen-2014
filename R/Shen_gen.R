
# Shen's data generation method from his 2014 paper.

Shen_gen <- function(n= 100, mu=0.5, d0=1, mc=6 )
{
  # Genereate the times for event of interest from random exponenital dist.
  Ti = rexp(n, rate = 1)
  
  # Generate the time for entering the study from random exp dist.
  V = rexp(n, rate = 1/mu)
  
  # Only keep samples if Vi* <=  Ti* and generate values until we have n values where Vi <= Ti
  indexes = which(V <= Ti)
  Ti = Ti[indexes]
  V  = V[indexes]
  while( length(Ti) < n )
  {
    new_Ti  = rexp(n, rate = 1)
    new_V   = rexp(n, rate = mu)
    indexes = which(new_V <= new_Ti)
    new_Ti  = new_Ti[indexes]
    new_V   = new_V[indexes]
    Ti      = append(Ti, new_Ti)
    V       = append(V, new_V)
  }
  Ti = Ti[1:n]
  V  = V[1:n]
  
  
  # generate censoring
  C = V + d0
  
  # Generate a random number of visits for each patient and group from binomial dist.
  N = 2 + rbinom(n, mc, 0.5)
  
  # Generate the X and Y vectors for each patient i in the first group
  X <- list()
  Y <- list()
  for(i in 1:n)
  {
    k = N[i]
    temp = runif(k,0,1)
    X <- append(X, list(temp))
    Y <- append(Y, list( cumsum(temp) + V[i] ))
  }
  
  
  L <- rep(-1,n)
  R <- rep(-1,n)
  
  # iterate through each patient i
  for(i in 1:n)
  {
    # iterate through each interval for patient i
    # The loop ends at Ni - 1, since Y_(j+1) is used 
    intervalCensored = FALSE
    for(j in 1:(N[i]-1))
    {
      # If Ti is in one of the Y intervals ...
      if( (Y[[i]][[j]] <= Ti[i]) && (Ti[i] <= Y[[i]][[j+1]])  )
      {
        # ... and the censoring time comes after the interval, then interval censored
        if(Y[[i]][[j+1]] < C[i] || Ti[i] < C[i])
        {
          L[i] = Y[[i]][[j]]
          R[i] = Y[[i]][[j+1]]
          intervalCensored = TRUE
          break
        }
      }
    }
    
    # if C < T   OR  if C and T are in same interval, then right censored
    if( ! intervalCensored )
    {
      # Right Censoring
      if( Ti[i] > C[i] || ( Y[[ i ]][[ N[i] ]] < Ti[i]  && Ti[i] < C[i] )   )
      {
        L[i] = C[i]
        R[i] = Inf
      }
      # Left censoring
      else if(Ti[i] < Y[[i]][[1]])
      {
        L[i] = V[i]
        R[i] = Y[[i]][[1]]
      }      
    }
  }
  
  A <- as.matrix( cbind(C, L, R, Ti) )
  A
}


