# Title: inbreed function
#
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150401
# Last Updated: 20150607

### Intention: Uses model from Crow and Kimura (1970; pg 269) to determine theoretical loss of heterozygosity (relative to that of initial population size) through generations. Takes in a matrix of population sizes and a matrix of migrant numbers. Calculates the inbreeding coefficient (F) through time representing the probability that two alleles are identical by descent (common ancestry).

# Note that this calculation generates the Ft value by NOT assuming self fertilization is possible. Thus, rather than the probability that both alleles came from a non migrant being (1-m.rate)*(1-m.rate) where m.rate is the migration RATE (m/pop.size), the probability is is instead (1-(m/pop.size))*(1-((m-1)/pop.size)). We subtract 1 from the number of migrants for the second allele because it must come from a different individual.
# This affects the calculation more strongly for small populations, which is relevant here. For example, if pop.size is 4 and 2 individuals are migrants, the migration rate is 0.5. Using the first formulation, the probability that both alleles come from a non migrant is (1-m.rate)^2 = (1-0.5)^2 = 0.25. Using the second formulation, the probability is (1-(2/4))*(1-((2-1)/4)) = 0.125

# pop.size should be a vector of population sizes that includes the migrants. For this experiment, that means the Nt dataframe.
# m should be a vector of migrant numbers
# Optional input F0 describes the inbreeding coefficient of the large external population. Perhaps, through the course of dozens of generations in the lab, this large population is somewhat inbred as well.

# inbreed <- function(N, m, F0=0)
# {
#   N <- as.numeric(N)
#   m <- as.numeric(m)
#   
#   F <- numeric(length(N)+1)
#   F[1:2] <- F0
#   
#   for (ft in 3:length(F))
#   {
#     t <- ft-1 # Time from perspective of population is offset by 1 generation
#     F[ft] <- F[ft-1] + ( (1-2*F[ft-1]+F[ft-2]) / (2*N[t]) )
#     migrant.correction <- ((N[t]-m[t])/N[t])*((N[t]-m[t]-1)/(N[t]-1))
#     F[ft] <- F[ft] * migrant.correction
#     if ((N[t] == 1) & (m[t] > 0)) F[ft] <- 1
#    }
#   
#   return(F[-1])
# }

# Update to fix bug where N[t]==1 would return an NaN due to dividing by 0 and then propagate through the rest of the time series
inbreed <- function(N, m, F0=0)
{
  N <- as.numeric(N)
  m <- as.numeric(m)
  
  F <- numeric(length(N)+1)
  F[1:2] <- F0
  
  for (ft in 3:length(F))
  {
    t <- ft-1 # Time from perspective of population is offset by 1 generation
    F[ft] <- F[ft-1] + ( (1-2*F[ft-1]+F[ft-2]) / (2*N[t]) )
    migrant.correction <- ((N[t]-m[t])/N[t])*((N[t]-m[t]-1)/(N[t]-1))
    # If only one individual, the migrant correction breaks down because it tries to divide by 0. If the 1 individual is a migrant, the probability of 2 of its alleles being identical by descent is 0, so it gets a 1. All cases should get a 0 here, because there was never a time when only 1 migrant was added.
    if (N[t] == 1) {migrant.correction <- ifelse(m[t] == 1, yes=1, no=0)}
    F[ft] <- F[ft] * migrant.correction
  }
  
  return(F[-1])
}