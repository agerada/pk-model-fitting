# We will now try to expand our model to include the earlier doses and 
# concentrations. These are: 
## pre-analytics ##
# Import data

d <- datasets::Theoph
d <- d[d$Subject==1,] # keep first Subject only
plot(conc ~ Time, data = d)

## end pre-analytics ##

X <- t(t(d$Time))
Y <- t(t(d$conc))
Y2 <- rexp(length(X))

# A one compartment model with first order absorption may be an appropriate
# model for this problem. 
# Firstly let's describe the model using R functions: 
dXgidt <- function(ka, Xgi) {
  - (ka * Xgi)
}
dXcdt <- function(ka, kel, Xc, Xgi){
  (ka * Xgi) - (kel * Xc)
}

# Using laplace transform, we can convert these equations to compute the drug amounts 
# Since we only have central concentrations to fit, we only need the laplace
# transform of dXcdt, which is: 
c_laplace <- function(t, ka, kel, V, dose){
  left_part <- (dose * ka) / (V * (ka - kel)) 
  right_part <- (exp( - kel * t) - exp( - ka * t))
  return(left_part * right_part)
}

# We now have a model that can give us calculations for the drug concentration 
# at input time points
# The next step is to attempt gradient descent again, to find the optimal paramaters. 
# In the previous lecture, we noted that we need to calculate the gradient for this, using 
# the partial derivatives of each parameter. 
# This is the first hurdle. The partial derivatives of models of such complexity 
# are very difficult to compute by hand. We need to resort to a method that allows us
# to numerically estimate these derivatives. For this we require the Finite 
# Differences method. We can estimate the partial derivative with respect to a 
# parameter x by evaluating the function at x + h (sometimes notated as epsilon)
# For example, for f(x), we can estimate f'(x): 
# f'(x) = f(x + h) - f(x) / h
# Remember that the function that we are trying to optimise is the loss function, 
# NOT the model. Otherwise the gradient descent method is not really changed. 
# I have updated the wrapper below, since we now need to introduce another 
# variable/hyperparameter - epsilon
gradient_descent_analytical <- function(t, y, learning_rate, epochs, epsilon,  
         analytical_function, 
         starting_params, 
         dose, 
         converge_threshold = 0.001
         ) {
  i <- 1
  params <- append(starting_params, list('dose' = dose, 't' = t))
  y_pred <- do.call(analytical_function, params)
  old_loss <- sum((y - y_pred)**2)
  
  while (i < epochs) {
    for(j in seq_along(starting_params)){
      y_pred <- with(
        data.frame(
          append(starting_params, list('dose' = dose))
          ), 
        analytical_function())
      old_param <- starting_params[[j]]
      starting_params[[j]] <- starting_params[[j]] + epsilon
      y_pred_plus_e <- with(
        data.frame(
          append(starting_params, list('dose' = dose))
          ), 
        analytical_function())
      loss <- sum((y - y_pred)**2)
      loss_plus_e <- sum((y - y_pred_plus_e)**2)
      starting_params[[j]] <- old_param - learning_rate * ((loss_plus_e - loss) / epsilon)
    }
    if(abs(loss - old_loss) < converge_threshold){
      print(paste("Converged at", i, "iterations"))
      break; 
    } else{ old_loss <- loss}
    i <- i+1
  }
  return(starting_params)
}

# Let us try out this gradient descent function
# One of the inputs is the dose, which is present in the dataset, however 
# please note that that is mg/kg, therefore we need to multiple by the weight!
set.seed(25) # this is just so that you get the same results on your PC
dose_pt_1 <- 4.02 * 79.6
estimated_params <- gradient_descent_analytical(X, Y, learning_rate = 0.0001, epochs = 10000, 
                            epsilon = 0.001, 
                            analytical_function = c_laplace, 
                            starting_params = list(
                              'ka' = runif(1,0.1, 2), 'kel' = runif(1, 0.1,2), 'V' = runif(1, 0.001,10)
                            ), 
                            dose = dose_pt_1)
# We now need to use these "optimised" parameters to draw a curve and plot this
# with the observed concentrations, in order to assess the model fit.
# We do this by using the Laplace function for the model, using these parameters
# as inputs. The with() function allows us to do this easily. This is an advanced 
# R function, and is not necessary (we could manually call c_laplace and pass
# references to the new parameters). We can then plot this curve. 
predict_y <- with(
  data.frame(
    (estimated_params)
  ), sapply(predict_x, function(x) c_laplace(x, ka, kel, V, dose)
))
plot(conc ~ Time, data = d)
lines(predict_x, predict_y)
# We get a pretty poor fit. Almost certainly a local minium. Our volume of 
# distribution is also probably too low. As in the previous session, we should 
# try multiple starting values for the parameters, and we can also widen the 
# random starting values. 
runs <- 10
i <- 1
ka <- runif(runs,0.01, 2)
kel <- runif(runs, 0.01, 2)
V <- runif(runs, 0.01,30)
# the following will just be used as vectors to store the parameter estimates
ka_est <- NULL
kel_est <- NULL
V_est <- NULL
loss_log <- NULL

while(i <= runs) {
  print(paste("Run", i, "of", runs))
  ka <- runif(runs,0.01, 2)
  kel <- runif(runs, 0.01, 2)
  V <- runif(runs, 0.1,30)
  starting_params = list(
    'ka' = ka[[i]], 'kel' = kel[[i]], 'V' = V[[i]]
  )
  starting_params <- gradient_descent_analytical(X, Y, learning_rate = 0.0001, epochs = 10000, 
                                                  epsilon = 0.001, 
                                                  analytical_function = c_laplace, 
                                                  starting_params, 
                                                  dose = 4.02 * 79.6)
  
    predict_y <- with(
      data.frame(
        (starting_params)
      ), sapply(predict_x, function(x) c_laplace(x, ka, kel, V, dose)
      ))
    ka_est <- c(ka_est, starting_params[['ka']])
    kel_est <- c(kel_est, starting_params[['kel']])
    V_est <- c(V_est, starting_params[['V']])
    loss_log <- c(loss_log, sum((c_laplace(X, starting_params[['ka']], 
                                           starting_params[['kel']], 
                                           starting_params[['V']], 
                                           dose))
                                - Y)**2)
    # if you would like to see the plots as they are estimated, just 
    # uncomment the following two lines
    #plot(conc ~ Time, data = d)
    #lines(predict_x, predict_y)
    i <- i+1
}

# The next few lines just pick the parameters that produce the minimum loss
# function and plot out the model
results <- data.frame(ka=ka_est, kel=kel_est, V=V_est, loss_log)
predict_y <- with(
  results[which.min(results$loss_log),], 
  sapply(predict_x, function(x) c_laplace(x, ka, kel, V, dose)
  ))
plot(conc ~ Time, data = d)
lines(predict_x, predict_y)
# You should be able to see that we get a decent fitting model, visually
# Just before concluding, let us try to visualise the loss function over 
# a parameter sweep: 
param_ka <- seq(from=0.01, to=2, length.out=100)
param_kel <- seq(from=0.01, to=2, length.out=100)
param_v <- seq(from=1, to = 50, length.out=100)
param_grid <- expand.grid(ka=param_ka,kel=param_kel,v=param_v)
param_grid$dose <- dose_pt_1

loss_log <- vector(mode='numeric', length=nrow(param_grid))
for(i in 1:nrow(param_grid)) {
  loss_log[[i]] <- sum( (c_laplace(X, param_grid[i,'ka'], 
                                  param_grid[i,'kel'], 
                                  param_grid[i,'v'], 
                                  dose) - Y) **2)
}

sample_vector <- sort(sample.int(length(loss_log), 10000))
plot(log(loss_log[sample_vector]) ~ param_grid[sample_vector,'v'])
# We do see that we seem to be minimising the loss function around values 
# approx 22-30L for Volume. 
# We can use 3d plots to visualise the impact of 2 parameters on the loss 
# function: 
param_ka <- seq(from=0.01, to=2, length.out=100)
param_kel <- seq(from=0.01, to=2, length.out=100)
param_v <- results[which.min(results$loss_log),"V"]
param_grid <- expand.grid(ka=param_ka,kel=param_kel,v=param_v)
param_grid$dose <- dose

loss_log <- vector(mode='numeric', length=nrow(param_grid))
loss_matrix <- matrix(nrow=length(param_ka), ncol=length(param_kel))
for(i in seq_along(param_ka)) {
  for(j in seq_along(param_kel)){
    loss_matrix[i,j] <- sum( (c_laplace(X, param_ka[i], 
                                     param_kel[j], 
                                     param_v, 
                                     dose) - Y) **2)
  }
  
}

library(plotly)
plot_ly(data.frame(param_ka, param_kel), z = loss_matrix) %>% add_surface()
