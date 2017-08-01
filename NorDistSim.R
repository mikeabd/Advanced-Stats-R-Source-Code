#Installing required packages
install.packages(c("dplyr", "ggplot2", "forecast", "lubridate", "ROCR"
                   , "rjson", "timeSeries", "xts", "highfrequency"
                   , "manipulate", "tseries", "bsts", "boot", "mixtools"))

# Loading Required Packages
library(dplyr)
library(ggplot2)
library(forecast)
library(lubridate)
library(ROCR)
library(rjson)
library(timeSeries)
library(xts)
library(highfrequency)
library(manipulate)
library(tseries)
library(bsts)
library(boot)
library(mixtools)

# Setting work directory
setwd("C:/Users/hosse/Desktop/ML/NormalDistSimulation")
paste("Work Directory is Set Up to", getwd())

# Set seed value
set.seed(20)

#######################################################
####### Bivariate Normal Distribution Simulation ######
#######################################################

# Number of random samples
N <- 200

# Normal Distribution Parameter Setting
rho <- -0.6
mu1 <- 1; s1 <- 2
mu2 <- 1; s2 <- 8

Mu <- c(mu1, mu2) 
Sigmas <- matrix(c(s1^2, s1*s2*rho, s2*s1*rho, s2^2), 2)

# Function to draw ellipse for bivariate normal data
ellipse_bvn <- function(bvn, alpha, col){
  Xbar <- apply(bvn,2,mean)
  S <- cov(bvn)
  ellipse(Xbar, S, alpha = alpha, col=col, npoints = 250, newplot = FALSE, draw = TRUE)
}

# Sampling data for Simulation
gibbs <- function (n, mu1, s1, mu2, s2, rho) 
{
  mat <- matrix(ncol = 2, nrow = n)
  x <- 0
  y <- 0
  mat[1, ] <- c(x, y)
  for (i in 2:n) {
    x <- rnorm(1, mu1+(s1/s2)*rho*(y - mu2), sqrt((1 - rho^2)*s1^2))
    y <- rnorm(1, mu2+(s2/s1)*rho*(x - mu1), sqrt((1 - rho^2)*s2^2))
    mat[i, ] <- c(x, y)
  }
  mat
}

bvn <- gibbs(N,mu1,s1,mu2,s2,rho)
colnames(bvn) <- c("X1","X2")

# Plotting Results
plot(bvn,xlab="X1", ylab="X2", main = " Bivariate Normal Distribution Simulation Plot "
     , col="black")
ellipse_bvn(bvn, 0.5, col = "red")
ellipse_bvn(bvn, 0.05, col = "blue")


###########################################################
####### Multi-variate Normal Distribution Simulation ######
###########################################################

install.packages("mvtnorm")
library(mvtnorm)

# Function to claculate mean and cov of a given data
sampleSize <- 1000
nbrVars <- 3
data_mat <- function (nbrVars, sampleSize, LB, UB)
{
  data <- matrix(ncol = nbrVars, nrow = sampleSize)
  mean_mat <- matrix(ncol = nbrVars, nrow = 1)
  covar_mat <- matrix(ncol = nbrVars, nrow = nbrVars)
  sigma_mat <- matrix(ncol = nbrVars, nrow = 1)
  
  for(i in 1:nbrVars){
    data[, i] <- sample(LB:UB, sampleSize, replace = TRUE)
    mean_mat[, i] <- mean(data[, i])
    sigma_mat[, i] <- var(data[, i])
  }
  covar_mat <- cov(data, method = "pearson")
  output <- list("Data" = data, "Mean Matrix" = mean_mat
                 , "Covariance Matrix" = covar_mat, "Variance" = sigma_mat)
  return(output)
}

# Varibale to hold the simulation results
mvn <- rmvnorm(n = 1000, unlist(as.list(data_mat(nbrVars, sampleSize, 1, 10)[2])))

for(i in 1:(ncol(mvn)-1)){
  for(j in 1:ncol(mvn)){
    if(j > i){
      plot_main <- paste("Var", i, " and Var ", j, " Scatter Plot", sep="")
      plot(x = mvn[, i], y = mvn[, j], xlab= paste("X", i), ylab= paste("X", j)
           , main = plot_main , col= i*j)
    }
  }
}



