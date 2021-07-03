choice_probability <- function(choiceset, parameters) {
# A function to calculate the probability of taking a risky option
# using a prospect theory model.
# Assumes parameters are [rho, lambda, mu] as used in S-H 2009, 2013, 2015, etc.
# Assumes choiceset has columns riskygain, riskyloss, and certainalternative.
# Creates binary choices, with True/1 = risky, False/0 = safe.
#
# Peter Sokol-Hessner
# July 2021

  # extract  parameters
  rho = as.double(parameters[1]); # risk attitudes
  lambda = as.double(parameters[2]); # loss aversion
  mu = as.double(parameters[3]); # choice consistency
  
  # calculate utility of the two options
  utility_risky_option = 0.5 * choiceset$riskygain^rho + 
                        -0.5 * lambda * abs(choiceset$riskyloss)^rho;
  utility_safe_option = choiceset$certainalternative^rho;
  
  # normalize values using this term
  div <- max(choiceset)^rho; # decorrelates rho & mu
  
  # calculate the probability of selecting the risky option
  p = 1/(1+exp(-mu/div*(utility_risky_option - utility_safe_option)));
  
  return(p)
  }
