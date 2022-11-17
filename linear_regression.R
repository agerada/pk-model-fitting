# Import data
d <- datasets::Theoph
d <- d[d$Subject==1,] # keep first Subject only
d <- d[4:nrow(d),] # keep only the linear drug concentrations

# Linear model using built in R function
model <- lm(conc ~ Time, data = d)
plot(conc ~ Time, data = d)
abline(model)

# Scale variables
# Intercept is now = 0
d$Time_scaled <- scale(d$Time)
d$conc_scaled <- scale(d$conc)

# Store the independant and dependent variables in matrix notation
X <- d$Time_scaled
Y <- d$conc_scaled

# Closed form mathematical solution
closed <- solve(t(X) %*% X) %*% t(X) %*% Y
 
# Sum of squares loss function
sum_sq <- function(y_pred, y_measured) {
  sum((y_pred - y_measured)**2)
}

# Linear prediction function (y = mx  + c)
linear_pred <- function(x, m, c) {
  (x %*% m) + c
}

# Estimate m through brute force
i <- 1
n <- 1000
random_m <- runif(1, min = -2, max = 2)
loss_function <- sum_sq(linear_pred(X, random_m, 0), Y)
while (i < n) {
  temp_b <- runif(1, min = -5, max = 5)
  loss_function_temp <- sum_sq(linear_pred(X, temp_b, 0), Y)
  if (loss_function_temp < loss_function) {
    random_m <- temp_b
    loss_function <- loss_function_temp
  }
  i <- i + 1
}
# We end up with random_m that is optimised for the minimum loss function 
# within the iteration threshold 

# Let us visualise the the loss function against a broad range of 
# parameter m values (-10 to 10). This is called a parameter sweep or
# grid approximation. 
param_sweep <- seq(from = -10, to = 10, length.out = 100)
sweep <- lapply(param_sweep, function(param) linear_pred(X, param, 0))
sweep_loss <- sapply(sweep, function(s) sum_sq(s, Y))
plot(sweep_loss~param_sweep)

# We have a clear parabola that we are trying to minimise. Lets zoom 
# in on the area between -2 and 0. 
param_sweep_fine <- seq(from = -2, to = 0, length.out = 100)
sweep_fine <- lapply(param_sweep_fine, function(param) linear_pred(X, param, 0))
sweep_fine_loss <- sapply(sweep_fine, function(s) sum_sq(s, Y)) 
plot(sweep_fine_loss~param_sweep_fine)

# We can minimise the loss function using this method as well: 
param_grid_approx <- seq(from = -2, to = 0, length.out = 100)
grid_approx_sweep <- lapply(param_grid_approx, 
                            function(param) linear_pred(X, param, 0))
grid_approx <- sapply(grid_approx_sweep, function(s) sum_sq(s, Y))
# Using this method, this is the parameter that gives the lowest loss fn: 
param_grid_approx[which.min(grid_approx)]

# Let us now log the brute force estimations, and see when we reach 
# optimisation
# I have wrapped the simulation in a function that returns a dataframe of
# random optimisations, with the running minimum loss function (and respective
# parameter)
# Parameters: 
# loss_f = loss function that will be optimised (by default minimised)
# linear_pred_f = linear prediction model function
# X = independent variable used for linear model
# Y = dependent variable (observed)
# intercept = intercept used for linear model, default = 0
# iters = number of iterations to run optimisation, default = 1000
# lower_bound = lower bound for parameter space
# higher_bound = higher bound for parameter space
stochastic_optimisation <- function(loss_f, linear_pred_f, 
                                    X, Y, intercept = 0, iters = 1000, 
                                    lower_bound, higher_bound) {

  output_i <- seq(1,iters)
  i <- 1
  random_m <- runif(1, min = lower_bound, max = higher_bound)
  loss <- loss_f(linear_pred_f(X, random_m, intercept), Y)
  current_min_m <- random_m
  current_min_loss <- loss
  output_m_log <- current_min_m
  output_min_loss_log <- current_min_loss
  while (i < iters) {
    temp_m <- runif(1, min = -2, max = 2)
    loss_temp <- loss_f(linear_pred_f(X, temp_m, intercept), Y)
    if (loss_temp < current_min_loss) {
      current_min_m <- temp_m
      current_min_loss <- loss_temp
    } 
    output_m_log <- c(output_m_log, current_min_m)
    output_min_loss_log <- c(output_min_loss_log, current_min_loss)
    i <- i + 1
  }
  data.frame(iter = output_i, loss = output_min_loss_log, parameter = output_m_log)
}

# Let's try this out: 
n_runs <- 100
X <- t(t(d$Time)) # Running this on unscaled variables
Y <- t(t(d$conc))
# Let's try one attempt first, with default 1000 iterations
run1 <- stochastic_optimisation(sum_sq, linear_pred, X = X, Y = Y, 
                                intercept = 9.9670, # since unscaled, provided intercept from lm model
                                lower_bound = -2, 
                                higher_bound = 2)

# Seems to work, let's now repeat the experiment 100 (n_runs) times
run_many <- lapply(seq(n_runs), function(x) 
  stochastic_optimisation(sum_sq, linear_pred, X = X, Y = Y, 
                          lower_bound = -2, 
                          higher_bound = -2))

# And plot these
plot(NULL,xlim=c(0,250), ylim=c(0,2000), type='n', xlab='Iteration', ylab='Minimum loss function')
i <- 1
for(i in seq_along(run_many)) {lines(run_many[[i]]$loss ~ run_many[[i]]$iter, type='l')}
# As we go left to right, we can see the experiments minimising the loss function. 
# We see that using this method, we generally have optimised the loss function
# in about 75 iterations. 