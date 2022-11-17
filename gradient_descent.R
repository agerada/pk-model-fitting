## pre-analytics ##
## Just re-running these to ensure we have the correct data
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

# We can recall that the linear function is y = mx + c
# and loss function is sum(y_obs - y_pred)**2
# or sum(y_obs - (mx + c))**2
# We want to calculate the gradient of the loss function (direction of steepest
# ascent). 
# We need the partial derivatives with respect to m and c: 
# L(m,c) = sum(y_obs - (mx+c))**2
# dLdm = 2*(sum(y_obs - (m*x+c))) * -x = -2*x*(y-(w*x+c))
# dLdc = -2*(sum(y_obs - (m*x+c)))

gradient_descent_linear <- function(x, y, learning_rate, epochs, 
                             linear_pred_f, 
                             starting_m = rnorm(1, 0), 
                             starting_c = rnorm(1,0)) {
  # Gradient descent function for a linear model
  # Parameters: 
  # x = data in (time points), y = observed values for y
  # learning_rate = this is a hyperparameter that indicates the size of our steps
  # if learning rate is high, we progress quickly, but risk overshooting
  # if it is too low, algorithm will not find minimum in a reasonable time
  # epochs = number of iterations to attempt
  # starting_m = starting m for algorithm, default is Normal(mean=0,sd=1)
  # starting_c = starting c for algorithm, default is Normal(mean=0,sd=1)
  # Function returns a list with the m and c parameters
  m <- starting_m
  c <- starting_c
  n <- length(x)
  dldm <- 0
  dldc <- 0
  i <- 1
  while (i < epochs) {
    y_pred <- linear_pred_f(x,m,c)
    dldm <- (-2 / n) * sum(x * (y - y_pred))
    dldc <- (-2 / n) * sum(y - y_pred)
    # Note that below we take the negative of learning_rate * derivative, since
    # we want to go down the slope
    m <- m - learning_rate * dldm
    c <- c - learning_rate * dldc
    i <- i+1
  }
  return(list('m' = m, 'c' = c))
}

# Let's try this out
gradient_descent_linear(X, Y, learning_rate = 0.001, epochs = 1000, 
                 linear_pred_f = linear_pred, 
                 starting_m = 1, 
                 starting_c = 10)
# Great, we get very similar values to the lm model
# What happens if we choose a bigger learning rate? 
gradient_descent_linear(X, Y, learning_rate = 0.01, epochs = 1000, 
                 linear_pred_f = linear_pred, 
                 starting_m = 1, 
                 starting_c = 10)
# We get ridiculous values. The algorithm was swinging everywhere! 

# We can now visualise the gradient descent 
# We need to wrap the gradient descent function around another while loop with 
# another epochs value. Then we run the gradient descent twice per iteration. 
epochs_extern <- 50
j <- 1
params <- list('m' = -0.2, 'c' =10) # starting parameters

# the following lines of code are not new, they just plot the loss function
# parabola
param_sweep <- seq(from = -0.5, to = -0.1, length.out = 100)
sweep <- lapply(param_sweep, function(param) linear_pred(X, param, params[['c']]))
sweep_loss <- sapply(sweep, function(s) sum_sq(s, Y))
plot(sweep_loss~param_sweep, pch=1, cex =0.2)

# I am going to add another feature to the algorithm, convergence. If the 
# algorithm finds that we are not making significant changes to parameters with 
# further iterations, then we assume we have converged, and terminate. 
converge_threshold <- 0.0001
while (j < epochs_extern) {
  # with every iteration, we calculate the current loss function with the current
  # parameters, and plot this as a big dot
  loss <- sum_sq(linear_pred(X,params[['m']],params[['c']]), Y)
  points(params[['m']], loss, cex = 2)
  
  # now we use gradient descent to take a step downwards
  params_new <- gradient_descent_linear(X, Y, learning_rate = 0.001, epochs = 2, 
                   linear_pred_f = linear_pred, 
                   starting_m = params[['m']], 
                   starting_c = params[['c']])
  
  # store the difference between new and old parameters
  m_difference <- abs(params$m - params_new$m)
  c_difference <- abs(params$c - params_new$c)
  # if the difference is insignificant, we have converged
  if (m_difference  < converge_threshold & c_difference < converge_threshold)
  {
    print(paste0("Converged at ", j, " iterations"))
    break # stop iterating
  }
  # otherwise, set the working parameter list to the new estimation and keep going
  params <- params_new
  j <- j+1
  Sys.sleep(0.25) # we must pause of 0.25 of a second, otherwise too quick to visualise
}

# Looks good! Try changing parameters, learning rate, convergence, epochs, etc
# and see what effect it has.. 