Shen_gen <- function(n= 100, bt=500, iterations = 5000, mu=0.5, d0=1, mc=6 )
{
  # Genereate the times for event of interest from random exponenital dist.
  Ti = rexp(n, rate = 1)

  # Generate the time for entering the study from random exp dist.
  V = rexp(n, rate = mu)

  # Only keep samples if Vi* <=  Ti* and generate values until we have n values where Vi <= Ti
  indexes = which(V <= Ti)
  Ti = Ti[indexes]
  V = V[indexes]
  while( length(Ti) < n )
  {
    new_Ti = rexp(n, rate = 1)
    new_V = rexp(n, rate = mu)
    indexes = which(new_V <= new_Ti)
    new_Ti = new_Ti[indexes]
    new_V = new_V[indexes]
    Ti = append(Ti, new_Ti)
    V = append(V, new_V)
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
  
  for(i in 1:n)
  {
    sameIntervalFlag = FALSE
    ### interval censoring
    for(j in 1:(N[i]-1))
    {
      # If Ti is on one of the Y intervals ...
      if( (Y[[i]][[j]] <= Ti[i]) && (Ti[i] <= Y[[i]][[j+1]])  )
      {
        # ... and the censoring time is in the same interval, then right censored
        if( (Y[[i]][[j]] <= C[i]) && (C[i] <= Y[[i]][[j+1]])  )
        {
          sameIntervalFlag = TRUE
        }
        # ... and the censoring time comes after the interval, then interval censored
        else if(Y[[i]][[j+1]] < C[i])
        {
          L[i] = Y[[i]][[j]]
          R[i] = Y[[i]][[j+1]]
        }
      }
    }
    
    # if C < T   OR  if C and T are in same interval, then right censored
    ### Right Censoring
    if(Ti[i] > C[i] || sameIntervalFlag)
    {
      L[i] = C[i]
      R[i] = 10000
    }
    
    ### Left censoring
    if(Ti[i] < Y[[i]][[1]])
    {
      L[i] = V[i]
      R[i] = Y[[i]][[1]]
    }
    
  }
  
  A <- cbind(X,Y,V,Ti,L,R,C)
  A
}



Group1 = Shen_gen(mu=0.5, d0=1,mc = 6)
print(Group1)

Group2 = Shen_gen(mu=2, d0=2,mc = 10)
print(Group2)



