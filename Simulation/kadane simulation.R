# Testing the method
library(MASS)

# simulation parameters
res <- data.frame(set.r = NA, match = NA, r = NA, bias = NA, ci = NA, width = NA)
nsim <- 5000
match <- c(FALSE, TRUE)
corr <- c(seq(0, 0.9, by = 0.1), 0.99)
set.seed(134589)

pb <- txtProgressBar(min = 0, max = length(match)*length(corr)*nsim, style = 3)
  for(j in 1:length(match)){
    for(r in 1:length(corr)){
      simres <- data.frame(r = NA, bias = NA, ci = NA, width = NA)
      for(i in 1:nsim){
        setTxtProgressBar(pb, ((j-1)*(length(corr)*nsim)) + ((r-1)*nsim) + i)
        data <- mvrnorm(n = 100, mu = c(0, 0.5, 1), Sigma = matrix(c(1, corr[r], 0.5, 
                                                                     corr[r], 1, 0.5,
                                                                     0.5, 0.5, 1), nrow = 3))
        obs <- as.data.frame(data)
        obs[1:50, 1] <- NA
        obs[51:100, 2] <- NA
        
        # impute
        kadaneimp <- mice(obs, method = c("kadane", ""), kadane.match = match[j], 
                          kadane.corr = corr[r], blocks = list(c("V1", "V2"), c("V3")), 
                          maxit = 1, m = 1, printFlag = FALSE)
        imp <- complete(kadaneimp, action = "long")
        # evaluate - bias, ci, width, realised correlation
        simres[i, "r"] <- cor(imp[3:5])[1,2]
        biases <- c(imp[1:50, "V1"] - data[1:50, 1], imp[51:100, "V2"] - data[51:100, 2])
        simres[i, "bias"] <- mean(biases)
        ci <- quantile(biases, probs = c(0.025, 0.975), na.rm = TRUE)
        simres[i, "ci"] <- ifelse(ci[1] < 0 & ci[2] > 0, 1, 0)
        simres[i, "width"] <- abs(ci[1] - ci[2])
        
      }
      store <- ifelse(j == 1, r, r + length(corr))
      res[store, "set.r"] <- corr[r]
      res[store, "match"] <- match[j]
      res[store, "r"] <- mean(simres[, "r"], na.rm = T)
      res[store, "bias"] <- mean(simres[, "bias"], na.rm = T)
      res[store, "ci"] <- mean(simres[, "ci"], na.rm = T)
      res[store, "width"] <- mean(simres[, "width"], na.rm = T)
  }
}
close(pb)

# Results
res

# Plot it
library(ggplot2)
library(cowplot)

plotIt <- ggplot(res, aes(x = set.r, colour = match)) + theme_minimal()
plot_grid(plotIt + geom_line(aes(y = r)) + labs(y = "imputed correlation"),
          plotIt + geom_line(aes(y = bias)),
          plotIt + geom_line(aes(y = ci)) + lims(y = c(0.9, 1)) + labs(y = "coverage rate") +
            geom_hline(yintercept = 0.95, linetype = 2),
          plotIt + geom_line(aes(y = width)) + labs(y = "average width"),
          nrow = 2)
