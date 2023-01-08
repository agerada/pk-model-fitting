# Desing and fitting of PK models

## Summary of files

* linear_regression.R = R script that implements a few important functions: 
    - least squares loss function
    - analytical solution for linear regression slope
    - grid and random fitting function

* curse_of_dimensionality.R = this script explores the impact of adding one 
parameter to our grid approximation

* gradient_descent.R = this script implements a gradient descent algorithm that 
looks for the best fitting parameters for a linear model. We can also visualise
the fit as it is optimised. The main function within this script is: 
    - gradient_desent_linear
    
* compartments.nlogo = this is a visual representation of a one compartment 
model with first order elimination, using the system dynamics module in  [NetLogo](https://ccl.northwestern.edu/netlogo/)

* compartments2.nlogo = same as above, but added a gut compartment for first-
order absorption

* compartments1.R = this script implements the first ODE compartment model on the 
Theopylline dataset, and introduces two new important concepts: 
    - analytical solutions using Laplace transform
    - numerical estimation of loss function partial derivatives using the 
    Finite Differences method
    
* urine.R = implements a one compartment model with two first order elimination
pathways for a drug with urine measurements and presumed extra-renal metabolism 
