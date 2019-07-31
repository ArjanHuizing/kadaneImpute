# kadaneImpute
Based on the 2001 paper 'Statistical matching: A paradigm for assesing the uncertainty in the procedure' by Moriarity and Scheuren. kadaneImpute imputes missing values through linear regression, with the added ability to specify a correlation between outcomes. This can be useful when outcomes are never jointly observed, for example with potential outcomes.

The procedure can be run inside the R package `mice`, and should be specified for a block of atleast 2 continous variables with missings.

# To-do
- Use maximum likelihood to estimate the covariance matrix to prevent serious errors when estimating more complex data.
- Use robust residual variance estimation
- More options for correlation specification
