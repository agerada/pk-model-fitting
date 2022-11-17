# A quick detour into the curse of dimensionality
# I am just going to revisit the previous grid approximation method, but this 
# time also attempt to estimate the intercept 

## pre-analytics ##
# Import data
d <- datasets::Theoph
d <- d[d$Subject==1,] # keep first Subject only
d <- d[4:nrow(d),] # keep only the linear drug concentrations
model <- lm(conc ~ Time, data = d)
plot(conc ~ Time, data = d)
abline(model)

X <- t(t(d$Time))
Y <- t(t(d$conc))
## end pre-analytics ##

# This time, we do not scale the variables (since we want to estimate intercept)

# The core algorithm will be a nested loop, iterating through every permutation 
# of possible parameters for both m and c
# Since we now are estimating two parameters instead of 1, the intuition is that
# it will take twice as long to estimate. Let's see if the intuition holds up... 
# Therefore we need to start with sequences for m and c, to represent the 
# parameter range/sweep. 
param_m <- seq(from = -10, to = 10, length.out = 100) # spaced out to give a sequence of 100
param_c <- seq(from = 5, to = 10, length.out = 100)

# This is the key bit. This is performed first in a functional style. I am using 
# the purrr package because it offers more flexibility compared to lapply, but this
# could similarly be done using lapply. 
# This applies the inner function onto every combination of m and c within the
# sequence. 
sweep <- purrr::map_dfr(param_m, 
                        function(m) purrr::map(param_c, 
                                               function(c) 
                                                 {
                                                   prediction = linear_pred(X, m, c)
                                                   return
                                                   (
                                                     list(m_par = m, 
                                                          c_par = c, 
                                                          loss = sum_sq(linear_pred(X, m, c), Y))
                                                   )
                                                 }
                                               )
                        )

# Here is the same process done using procedural programming, which may be 
# more familiar. The results are the same. Computer programmers argue about which
# style is better or easier to read... 
m_par <- NULL
c_par <- NULL
loss <- NULL
for(m in param_m){
  for(c in param_c){
    prediction <- linear_pred(X, m, c)
    m_par <- c(m_par, m)
    c_par <- c(c_par, c)
    loss <- c(loss,sum_sq(linear_pred(X, m, c), Y))
  }
}
list(m_par = m_par, c_par=c_par, loss=loss)

# Now we plot the results (we are interested in the slope and intercept)
plot(conc ~ Time, data = d)
# We just need to find out which parameters led to the lowest loss function
# The which.min function gives us the location of this
minimum_loss_sweep <- sweep[which.min(sweep$loss),]

abline(minimum_loss_sweep$c_par, minimum_loss_sweep$m_par)

# Therefore we do manage to get the result, but even though we added 100 estimations of 
# one parameter, we end up doing 10,000 calculations in total (100 * 100)! 
# this is already becoming unmanageable, and our intuition was incorrect. 
# This is the curse of dimensionality. 

# we can also see that the problem is already less identifiable. There are many
# values of c that give a low loss function and and we are inefficiently
# trying to find the lowest value in a parabola shaped like a bathtub.  
library(tidyverse)
sweep %>% 
  distinct(m_par, c_par, .keep_all = T) %>% 
  ggplot(aes(x = m_par, y = c_par, z = loss)) + geom_contour_filled()
