# Start by loading up the data (ref Rowland and Tozer 5th ed, page 69)
dose <- 100
data <- data.frame(time = c(0.25, 0.5, 0.75, 1, 1.5, 2, 4, 6, 8, 12, 24, 30, 36, 48, 60, 72), 
                   plasma = c(3.0, 2.8, 2.4, 2.2, 2.0, 1.8, 1.4, 1.2, 1.1, 0.9, 0.55, 0.45, 0.36, 
                              0.23, 0.15, 0.1), 
                   urine_conc = c(rep(NA, 9), 13.3, 7.05, NA, 3.35, 3.32, NA, 1.91), 
                   urine_amount_interval = c(rep(NA, 9), 12.06, 6.70, NA, 4.13, 2.6, NA, 2.73))
# Record the cumulative total of urine amounts (this is a calculated field
# in the table)
data$urine_cumulative_amount <- cumsum(tidyr::replace_na(data$urine_amount_interval, 0) )

# We will use a one compartment model with two first order elimination pathways
p_laplace <- function(t, ke, km, v, dose) {
  # analytical solution for plasma concentration
  p_amount <- dose * exp( - ( (ke + km) * t) )
  p_amount / v
}

u_laplace <- function(t, ke, km, dose){
  # analytical solution for amount in urine "compartment"
  u_amount <- ( (ke * dose) / (ke+km) ) * (1 - exp( - ( (ke+km) * t ) ) )
  u_amount
}

u_integral <- function(from, to, ke, km, v, dose) {
  # integral to calculate urine amount released in a period of time
  to_calc <- ( - ( ke * v * dose) / (ke + km) ) * exp( - ( (ke + km) * to ) )
  from_calc <- ( - ( ke * v * dose) / (ke + km) ) * exp( - ( (ke + km) * from ) )
  to_calc - from_calc
}

# for this session, we will use internal R optimisation functions such as nlm()
# but in principle these are iterative gradient descent methods, so there is 
# no new concept there. This is just easier to use for now since we do not
# have to implement a more generalized gradient descent function 

# before progressing to the full model, let's just start with a simple model
# and just fit the plasma compartment data. 

loss_one_comp <- function(x) {
  # x = vector of parameters to optimise
  # x[1] = ke
  # x[2] = km
  # x[3] = v
  ke <- x[1]
  km <- x[2]
  v <- x[3]
  p_pred <- p_laplace(data$time, ke, km, v, dose)
  sum((p_pred - data$plasma)^2)
}

model1 <- nlm(loss_one_comp, p=c(0.2, 0.2, 10))
model1_pred <- p_laplace(data$time, model1$estimate[1], model1$estimate[2], model1$estimate[3], dose)
# we can access the estimated parameters using model1$estimate
# clearly the parameters are way off. nlm() function can be highly sensitive to
# starting values, so let's try with optim() instead. This uses algorithms that
# are less susceptible to bad starting values
model2 <- optim(c(1,1,1), loss_one_comp)
plot(data$plasma ~ data$time)
model2_pred <- p_laplace(data$time, model2$par[1], model2$par[2], model2$par[3], dose)
lines(model2_pred ~ data$time)

# ok we are still getting weird parameter values, probably because we have too many 
# elimination parameters on this model. Simple tweak should sort this out, 
# just using one elimination constant. This is a typical problem of non-
# identifiability. We need to first redefine a Laplace model with one elimination 
# constant, and then tweak the loss function slightly. 
one_comp_laplace <- function(t, ke, v, dose) {
  p_amount <- dose * exp( - ( (ke) * t) )
  p_amount / v
}
one_comp_loss <- function(x) { 
  # x = vector of parameters to optimise
  # x[1] = ke
  # x[2] = v
  ke <- x[1]
  v <- x[2]
  p_pred <- one_comp_laplace(data$time, ke, v, dose)
  sum((p_pred - data$plasma)^2)
}
model3 <- optim(c(1,1), one_comp_loss)
plot(data$plasma ~ data$time)
model3_pred <- one_comp_laplace(data$time, model3$par[1], model3$par[2], dose)
lines(model3_pred ~ data$time)
# now we have more reasonable parameters of ke = 0.107 and v = 38.4

# Let's proceed to the model that incorporates both urine and plasma
# The model functions are already defined, but we just need a loss function
# that can optimise both models to both sets of observations concurrently. 
# For the plasma concentrations, calculating the residuals is straight-forwarded, 
# and we have already seen this. 
# For urine this is more difficult, since we only have large volumes of urine 
# collected at particular intervals. 
# I decided to make the urine predictions by calculating the cumulative total
# of urine amounts, using the Laplace transform of urine compartment
# This is done using the cumsum() function 
loss <- function(x) { 
  # loss function for optimisation
  # x = vector of parameters to optimise
  # x[1] = ke
  # x[2] = km
  # x[3] = v
  ke <- x[1]
  km <- x[2]
  v <- x[3]
  p_pred <- p_laplace(data$time, ke, km, v, dose)
  u_pred <- u_laplace(data$time, ke, km, dose)
  u_pred <- cumsum(u_pred)
  p_residual <- p_pred - data$plasma
  u_residual <- ifelse(!is.na(data$urine_amount_interval), (u_pred - data$urine_cumulative_amount), 0)
  sum((p_residual - u_residual)^2)
}

# Let's try some plots. Remember that we have "missing" observations for the first 10
# data points when no urine was yet collected. We should skip these. 
model4 <- optim(c(0.1,0.1,10), loss)
plot(data$plasma ~ data$time, main="Plasma concentration (dots = measured, line = predicted", xlab = "Time", 
     ylab = "Concentration")
model4_pred_plasma <- p_laplace(data$time, model4$par[1], model4$par[2], model4$par[3], dose)
lines(model4_pred_plasma ~ data$time)

# OK - this looks weird. Let's look at urine: 
plot(data$urine_cumulative_amount[10:nrow(data)] ~ data$time[10:nrow(data)], 
     main="Urine cumulative amounts (predicted and observed", 
     xlab="Time", 
     ylab="Concentration")
lines(model4_pred_urine_cum ~ data$time)
# These look better, so let's look at the residuals for urine and plasma: 
model4_pred_urine <- u_laplace(data$time, model4$par[1], model4$par[2], dose)
model4_pred_urine_cum <- cumsum(model4_pred_urine)
plot(data$urine_cumulative_amount[10:nrow(data)] ~ model4_pred_urine_cum[10:nrow(data)], 
     main="Urine residuals", xlab="Predicted", ylab="Observed")
plot(data$plasma ~ model4_pred_plasma, main="Plasma residuals", 
     xlab="Predicted", ylab="Observed")

# It seems like we are getting reasonable results for urine, but the plasma 
# predictions are not good. I think this is because within the loss function, 
# we are comparing plasma concentrations and urine cumulative amounts. The 
# magnitude of the latter is much more than the former, so the loss function
# is heavily weighted towards giving us a good fit for the urine data, at the 
# expense of the plasma data. 
# There are probably different solutions to this problem, but using a mean 
# squared error for the loss function does equalise things a bit. The code
# below is the same as above, apart from the last line in the loss function. 
loss_mse <- function(x) { 
  # loss function for optimisation
  # x = vector of parameters to optimise
  # x[1] = ke
  # x[2] = km
  # x[3] = v
  ke <- x[1]
  km <- x[2]
  v <- x[3]
  p_pred <- p_laplace(data$time, ke, km, v, dose)
  u_pred <- u_laplace(data$time, ke, km, dose)
  u_pred <- cumsum(u_pred)
  p_residual <- p_pred - data$plasma
  u_residual <- ifelse(!is.na(data$urine_amount_interval), (u_pred - data$urine_cumulative_amount), 0)
  mean(sqrt(sum((p_residual - u_residual)^2)))
}

model5 <- optim(c(0.1,0.1,10), loss_mse)
plot(data$plasma ~ data$time, main="Plasma concentration (dots = measured, line = predicted", xlab = "Time", 
     ylab = "Concentration")
model5_pred_plasma <- p_laplace(data$time, model5$par[1], model5$par[2], model5$par[3], dose)
lines(model5_pred_plasma ~ data$time)

plot(data$urine_cumulative_amount[10:nrow(data)] ~ data$time[10:nrow(data)], 
     main="Urine cumulative amounts (predicted and observed", 
     xlab="Time", 
     ylab="Concentration")
lines(model5_pred_urine_cum ~ data$time)

model5_pred_urine <- u_laplace(data$time, model5$par[1], model5$par[2], dose)
model5_pred_urine_cum <- cumsum(model5_pred_urine)
plot(data$urine_cumulative_amount[10:nrow(data)] ~ model5_pred_urine_cum[10:nrow(data)], 
     main="Urine residuals", xlab="Predicted", ylab="Observed")
plot(data$plasma ~ model5_pred_plasma, main="Plasma residuals", 
     xlab="Predicted", ylab="Observed")
# Looks better, but still not quite there. Maybe adding a tissue compartment
# might help? That is one for next time. 