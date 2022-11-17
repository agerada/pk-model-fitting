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