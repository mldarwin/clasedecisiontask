binary_choice_from_probability <- function(probabilities) {
# A simple function to stochastically convert choice probabilities to logicals.
# Assumes probabilities are given in a 1-dimensional vector.
#
# Peter Sokol-Hessner
# July 2021
  
  binary_choices = probabilities > runif(length(probabilities));
  # If 'p' value is low, e.g. 0.2, then it will only be greater 
  # than a random number uniformly-distributed between 0 and 1
  # 20% of the time, that is, with p = 0.2
  
  return(binary_choices)
}
