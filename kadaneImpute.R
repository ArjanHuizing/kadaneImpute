# Description
  # Based on the 2001 Moriarity & Scheuren paper on data fusion. 
  # This is an imputation method for the package `mice` that allows a user to set 
  # a correlation between two or more outcomes which have no overlap.

# necessary libraries
library(mice)
library(Rcpp)
library(matrixcalc)

# impute.kadane function
mice.impute.kadane <- function(data, format = "imputes", kadane.corr = 0, ...){
  kad <- kadaneImpute(data, kadane.corr = kadane.corr, ...)
  return(mice:::single2imputes(kad, is.na(data)))
}

# covMiss function
sourceCpp(file = "covMiss.cpp")

# Kadane imputation function
kadaneImpute <- function(data, kadane.corr, kadane.match = TRUE, donors = 5L, ...){
  # find variables that need to be imputed
  miss <- which(colMeans(is.na(data)) > 0)
  np <- length(miss)
  if(np < 2){
    stop("Less than 2 potential outcomes used.")
  }
  if(length(kadane.corr) != choose(length(miss), 2) & length(kadane.corr) != 1){
    stop("Number of correlations `kadane.corr` needs to be equal to 1 or the number of outcomes 
         with missing values specified in the impute.kadane block")
  }
  
  # fill in missing covars with used-specified correlation
  covars <- covMiss(as.matrix(data))
  covarmiss <- combn(diag(covars)[miss], 2, FUN = function(x){(x[1]*x[2])^(1/2)*kadane.corr})
  fill <- combn(miss, 2)
  for(i in 1:ncol(fill)){
    covars[fill[1, i], fill[2, i]] <- covarmiss[i]
    covars[fill[2, i], fill[1, i]] <- covarmiss[i]
  }
  # check positive definite
  if(!matrixcalc::is.positive.definite(covars)){
    warning("Covariance matrix not positive-definite.")
  }
  
  # regression step (moriarity & scheuren 2001, pg 411)
  predicted <- matrix(NA, nrow = nrow(data), ncol = length(miss))
  for(i in 1:length(miss)){
    for(j in 1:nrow(data)){
      miss.j <- which(colMeans(is.na(data[j, , drop = FALSE])) == 1)
      # all outcomes i predicted for participant j, using participant j's
      # observed outcome + covariates.
      predicted[j, i] <- mean(data[, miss[i]], na.rm = TRUE) +
        (covars[miss[i], -miss.j]) %*% solve(covars[-miss.j, -miss.j]) %*% 
      t(data[j, -miss.j] - colMeans(data[, -miss.j, drop = FALSE], na.rm = TRUE))
    }
  }
  
  # Residual variance (moriarity & scheuren 2001, pg 413)
  resids <- numeric(0)
  for(i in 1:length(miss)){
    resids[i] <- covars[miss[i], miss[i]] - 
      (covars[miss[i], -miss[i]]) %*% solve(covars[-miss[i], -miss[i]]) %*%
      t(t(covars[miss[i], -miss[i]]))
  }
  if(min(resids) < 0){
    warning("negative variance detected in: ", 
            paste(names(miss[which(resids < 0)]), collapse = ", "))
  }
  
  # Add residual variance to predicted values, replace predicted with observed.
  imputed <- data
  for(i in 1:length(miss)){
    imputed[is.na(imputed[, miss[i]]), miss[i]] <- predicted[is.na(imputed[, miss[i]]), i] + 
      c(rnorm(n = nrow(predicted[is.na(imputed[, miss[i]]),]), 0, resids[i]))
  }
  
  # Match with predictive mean matching
  if(kadane.match){
    sets <- combn(miss, 2) 
   for(i in 1:ncol(sets)){
     #mahal <- mahalanobis(imputed[sets[, i]], 
      #                    center = colMeans(imputed[sets[, i]]), 
       #                   cov = covars[sets[, i],sets[, i]])

     set1 <- which(!is.na(data[, sets[1, i]]))
     set2 <- which(!is.na(data[, sets[2, i]]))
     match1 <- mice:::matcher(obs = data[set1, sets[1, i]], 
                              mis = imputed[set2, sets[1, i]], 
                              k = donors)
     match2 <- mice:::matcher(obs = data[set2, sets[2, i]], 
                              mis = imputed[set1, sets[2, i]], 
                              k = donors)
     data[set2, sets[1, i]] <- data[set1[match1], sets[1, i]]
     data[set1, sets[2, i]] <- data[set2[match2], sets[2, i]]
    return(data)
   } # end sets of outcomes loop
  } else{
    return(imputed)
  }# end matching conditional
} # end function
