#' Fit spatio-temporal Model
#' 
#' \code{idpSTM} fits an AR(1) temporal model on grouped data with information of space and time.
#' The input data is tilled as an n by n lattice and turned into a count data (per tile) 
#' before fitting a generallised linear model (GLM) with a logarithmic link, 
#' see details in Qiao et al. (2018).
#' The conditional distribution can be chosen to be either Poisson or negative binomial.
#' 
#' @references {
#' P. Qiao, C. Mølck, D. Ferrari, and F. Hollande. (2018) 
#' A spatio-temporal model and inference tools for 
#' longitudinal count data on multicolor cell growth. 
#' The International Journal of Biostatistics.
#' }
#' 
#' @param data A matrix with four cloumns (in order): time, x, y, group, where x and y are coordinates of location.
#' See the cell growth data as an example.
#' @param n An integer number indicating size of the grids. The image is tiled into an n*n grid. (Suggest n >= 6, otherwise neighbourhoods mostly overlapping with each other will lead to highly correlated covariates)
#' @param marginal Conditional distribution: "pois" for Poisson, "nbinom" for negative binomial. 
#' @param maxit Maximum number of iterations of maximum likelihood estimation (default 50). 
#' @param fit Logical, if TRUE, return fitted values. (only available for Poisson for now)
#' 
#' @return A list with components
#' \itemize{
#' \item \code{likelihood:}   Maximized log-likelihood. 
#' \item \code{coefficients:} A list of parameter estimates: 
#'                     \itemize{
#'                       \item \code{intercepts:}  A separate intercept for each group,
#'                        \item \code{main_effects:}  Regression coefficients orgainized in a K by K matrix, where K is the number of groups. 
#'                        The entry in the ith row, jth column is interpreted as the temporal impact of the precence of the ith group in the neighbourhood
#'                         on the growth of the jth group. 
#'                        \item \code{dispersion:} Overdispersion parameter if marginal == "nbinom". }
#' \item \code{standard_error:}  A list of estimated standard errors, in the same format as \code{coefficients}. 
#' }
#' 
#' @examples
#' n <- 25
#' data("cell_growth_data")
#' est_pois <- idpSTM(cell_growth_data, n, "pois")
#' print(est_pois$coefficients$main_effects)
#' print(est_pois$standard_error$main_effects)
#' 
#' \dontrun{
#' ## It's normal to have a few warnings related to sigular matrix inversion
#' ## when marginal == "nbinom". 
#' est_nbinom <- idpSTM(cell_growth_data, n, "nbinom")
#' print(est_nbinom$coefficients$main_effects)
#' print(est_nbinom$standard_error$main_effects)
#' }


idpSTM <- function(data, n, marginal, maxit = 50, fit = FALSE){
  
  K <- max(data[, 4])
  marginal <- which(c("pois", "nbinom") == marginal)
  
  dat <- tilling(data, n)
  res <- idptSTM_cpp(dat, n, marginal, maxit, fit)
  
  if(marginal == 1)
    ret_est <- list("intercept" = res$intercept, "main_effects" = res$main_effects)
  
  if(marginal == 2)
    ret_est <- list("intercept" = res$intercept, "main_effects" = res$main_effects, 
                    "dispersion" = res$dispersion)
  
  se <- matrix(sqrt(res$se), K + 1, K)
  if(marginal == 1) ret_se <- list("intercept" = se[1, ], "main_effects" = se[-1, ])
  
  if(marginal == 2)
    ret_se <- list("intercept" = se[1, ], "main_effects" = se[-1, ], 
                   "dispersion" = sqrt(res$se_dispersion))
    
  
  
  if(marginal == 1 && fit){
    return(list( "coefficients" = ret_est,
                 "likelihood" = res$likelihood, 
                 "standard_error" = ret_se, 
                 "fitted" = res$fit, 
                 "observed" = res$obs))
  }else{
    return(list( "coefficients" = ret_est,
                 "likelihood" = res$likelihood, 
                 "standard_error" = ret_se))
  }
  
}


#' Model selsction of spatio-temporal model
#'
#' \code{idpSTMSelect} selects the best model on \code{idpSTM} using \code{logGLMselect} 
#'
#' @param data A matrix with four cloumns: time, x, y, group, where x and y are coordinates of location.
#' @param n An integer number indicating size of the grids.
#'  The image is tiled into an n*n grid. (Suggest n >= 6, otherwise neighbourhoods mostly overlapping with each other will lead to highly correlated covariates)
#' @param marginal Conditional distribution: "pois" for Poisson, "nbinom" for negative binomial. 
#' @param maxit Maximum number of iterations of maximum likelihood estimation (default is 50). 
#' @param ModelCnt The number of models to be generated via Gibbs samplling (default is 100). 
#' @param Message Logical, if TRUE, prints the top 5 most frequently generated models. 
#'
#' @return A list with components
#' \itemize{
#' \item \code{likelihood:}   Maximized log-likelihood. 
#' \item \code{coefficients:} A list of parameter estimates: 
#'                     \itemize{
#'                       \item \code{intercepts:}  A separate intercept for each group,
#'                        \item \code{main_effects:}  Regression parameters indicating temporal impacts on growth between groups,
#'                        \item \code{dispersion:}  Overdispersion parameters if marginal = "nbinom".
#'                        }
#' \item \code{selected_model:}  A binary vector with 1 indicating selected variables and 0 otherwise.
#' \item \code{standard_error:}  A list of estimated standard errors, in the same format as \code{coefficients}. 
#'                               Standard errors are reported only for selected parameters, and are 0's otherwise. 
#' }
#' 
#' @examples
#' n <- 25
#' data("cell_growth_data") 
#' select_est <- idpSTMSelect(cell_growth_data, n, marginal = "pois", Message = FALSE)
#' print(select_est$selected_model)
#' 
#' ## It's normal to have  warnings related to sigular matrix inversion in solve()
#' ## when marginal == "nbinom". 
#' \dontrun{
#' select_est <- idpSTMSelect(cell_growth_data, n, marginal = "nbinom", Message = FALSE)
#' print(select_est$selected_model)
#' }



idpSTMSelect <- function(data, n, marginal, ModelCnt = 100, maxit = 50, Message = F){
  
  dat <- tilling(data, n)
  
  marginal <- which(c("pois", "nbinom") == marginal)
  res <- idpSTModelSelection_cpp(dat, n, marginal, maxit, ModelCnt, Message)
  
  if(marginal == 1){
    ret_est <- list("intercept" = res$parameters[1, ], "main_effects" = res$parameters[-1, ])
    ret_se <- list("intercept" = sqrt(res$se[1, ]), "main_effects" = sqrt(res$se[-1, ]))
    selected <- res$v
  }else{
    ret_est <- list("intercept" = res$parameters[1, ], "main_effects" = res$parameters[-c(1, nrow(res$parameters)), ], 
                    "dispersion" = res$parameters[nrow(res$parameters), ])
    ret_se <- list("intercept" = sqrt(res$se[1, ]), "main_effects" = sqrt(res$se[-c(1, nrow(res$se)), ]), 
                   "dispersion" = res$se[nrow(res$se), ])
    selected <- res$v[-nrow(res$v), ]
  }
  
  
  
  return(list("coefficients" = ret_est, "likelihood" = res$likelihood, "standard_error" = ret_se, 
              "selected_model" = selected))
}


#' Model selsction of log-GLM
#'
#' \code{logGLMselect} finds the best model among all possible models. 
#' Models are fitted with the generallised linear model (GLM) with logarithmic link, 
#' response distribution can be specified as Poisson or Negative binomial. 
#' Models are ranked with the Bayesian Information Criterion (BIC). 
#' The best models are found using a Gibbs samplling method
#' introduced by Qian and Field (2002), which allows very large candidate sets to be adressed.
#'
#' @param y A vector of count data response.  
#' @param x A matrix of covariates. 
#' @param marginal Conditional distribution: "pois" for Poisson, "nbinom" for negative binomial. 
#' @param maxit Maximum number of iterations of maximum likelihood estimation (default 50). 
#' @param skip A vector of indices of predictors to be forced in the selected model.
#' @param ModelCnt The number of models to be generated via Gibbs samplling (default 100).  
#' @param Message Logical, if TRUE, prints the top 5 most frequently generated models. (Suggest on)
#'
#' @return A list with components
#' \itemize{
#' \item \code{likelihood:}   Maximized log-likelihood. 
#' \item \code{coefficients:} A list of regression parameter estimates.
#' \itemize{
#'                       \item \code{intercepts:} ,
#'                        \item \code{main_effects:}  Regression coefficients,
#'                        \item \code{dispersion:}  Overdispersion parameter if marginal = "nbinom".
#'                        } 
#' \item \code{selected_model:}  A binary vector with 1 indicating selected variables and 0 otherwise.
#' \item \code{standard_error:}  A list of estimated standard errors in the same formatt as coefficients. 
#' }
#' 
#' @references {
#' G. Qian and C. Field. Using mcmc for logistic regression model selection involving large
#' number of candidate models. In Monte Carlo and Quasi-Monte Carlo Methods 2000, pages 460–474. Springer, 2002.
#' }
#'
#' @examples
#' ## Poisson 
#' nn <- 1000
#' x <- matrix(rnorm(5*nn), nn, 5)
#' intercept <- 1
#' true_beta <- c(1, 0, 0, 2, 3)
#' mu <- exp(x %*% true_beta + intercept)
#' y <- sapply(mu, rpois, n = 1)
#' 
#' res <- suppressWarnings(logGLMselect(y, x, "pois", maxit = 200, ModelCnt = 100, Message = TRUE))
#' cat("true model:", c(1, as.numeric(as.logical(true_beta))), "\n")

## Negative Binomial 
#' dispersion <- 1.5
#' y <- rep(0, nn)
#' for(i in 1:nn){
#'   y[i] <- rnbinom(1, size = dispersion, mu = mu[i])
#' }
#' res <- suppressWarnings(logGLMselect(y, x, "nbinom", maxit = 200, ModelCnt = 100, Message = TRUE))
#' cat("true model:", c(1, as.numeric(as.logical(true_beta)), 1), "\n")
#' cat("dispersion parameter: ", res$coefficients$dispersion, "\n")

logGLMselect <- function(y, x, marginal, 
                         maxit = 50, skip = NULL, ModelCnt = 100, Message = T){
  
  marginal <- which(c("pois", "nbinom") == marginal)
  
  if(nrow(x) != length(y)){ stop(" x must have the same number of rows as y")}
  if(length(unique(x[, 1])) > 1) x <- cbind(rep(1, nrow(x)), x)
  
  if(!length(skip)){ 
    skip <- 0 ## if not specified, skip intercept
    if(marginal == 2){
      skip <- c(skip, ncol(x)) ## skip also dispersion parameter
    } 
      }else{ 
    if(range(skip)[1] < 1 || range(skip)[2] > ncol(x)) stop("Indices in skip are out of range")
    skip = skip - 1} 
  
  res <- logGLMselect_cpp(y, x, marginal, maxit, skip, ModelCnt, Message)
  
  res_se <- rep(0, length(res$v))
  res_se[which(res$v == 1)] <- sqrt(res$se)
  if(marginal == 1){
    est <- list("intercept" = res$parameters[1], "coefficients" = res$parameters[-1])
    ret_se <- list("intercept" = res_se[1], "coefficients" = res_se[-1])
  }else{
    est <- list("intercept" = res$parameters[1], "coefficients" = res$parameters[-c(1, length(res$pparameters))], 
                "dispersion" = res$parameters[length(res$parameters)])
    ret_se <- list("intercept" = res_se[1], "coefficients" = res_se[-c(1, length(res_se))], 
                   "dispersion" = res_se[length(res_se)])
  }
  
  ret <- list("likelihood" = res$likelihood, "coefficients" = est, 
              "selected_model" = res$v, "standard_error" = ret_se)
 
  return(ret)
}




