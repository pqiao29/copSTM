#' Fit spatio-temporal Model
#' 
#' \code{idpSTM} takes a grouped data with information of spatial location and time point of each observation (typically a cell growth data), 
#' and fit a model of spatial neighbourhood impact between groups over time.
#' 
#' @param data A matrix with four cloumns: time, x, y, group, where x and y are coordinates of location.
#' @param n An integer number indicating size of the grids. The image is tiled into an n*n grid. (Suggest n >= 6, otherwise neighbourhoods mostly overlapping with each other will lead to highly correlated covariates)
#' @param maxit Maximum number of iterations of maximum likelihood estimation. 
#' @param fit Logical, if TRUE, return fitted values.
#' 
#' @return A list with components
#' \itemize{
#' \item \code{likelihood:}   Maximized log-likelihood. 
#' \item \code{coefficients:} A list of parameter estimates: 
#'                     \itemize{
#'                       \item \code{intercepts:}  A separate intercept for each group,
#'                        \item \code{main_effects:}  Regression parameters indicating temporal impacts on growth between groups.}
#' \item \code{standard_error:}  A list of estimated standard errors, in the same format as \code{coefficients}. 
#' }
#' 
#' @examples
#' n <- 25
#' data("cell_growth_data")
#' est <- idpSTM(cell_growth_data, n, maxit = 30)
#' print("======= Estimates =======")
#' print(est$coefficients)
#' print("======= Standard errors =======")
#' print(est$standard_error)

idpSTM <- function(data, n, maxit, fit = FALSE){
  
  K <- max(data[, 4])
  
  dat <- tilling(data, n)
  res <- idptSTM_cpp(dat, n, maxit, fit)
  
  ret_est <- list("intercept" = res$intercept, "main_effects" = res$main_effects)
  
  se <- matrix(sqrt(res$se), K + 1, K)
  ret_se <- list("intercept" = se[1, ], "main_effects" = se[-1, ])
  
  if(fit){
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
#' \code{idpSTMSelect} performs variable selection on model \code{idpSTM} using BIC via Gibbs samplling. 
#'
#' @param data A matrix with four cloumns: time, x, y, group, where x and y are coordinates of location.
#' @param n An integer number indicating size of the grids. The image is tiled into an n*n grid. (Suggest n >= 6, otherwise neighbourhoods mostly overlapping with each other will lead to highly correlated covariates)
#' @param maxit Maximum number of iterations of maximum likelihood estimation. 
#' @param ModelCnt The number of models to be generated via Gibbs samplling. 
#' @param Message Logical, if TRUE, prints the top 5 most frequently generated models. 
#'
#' @return A list with components
#' \itemize{
#' \item \code{likelihood:}   Maximized log-likelihood. 
#' \item \code{coefficients:} A list of parameter estimates: 
#'                     \itemize{
#'                       \item \code{intercepts:}  A separate intercept for each group,
#'                        \item \code{main_effects:}  Regression parameters indicating temporal impacts on growth between groups,
#'                        }
#' \item \code{selected_model:}  A binary vector with 1 indicating selected variables and 0 otherwise.
#' \item \code{standard_error:}  A list of estimated standard errors, in the same format as \code{coefficients}. 
#' }
#' 
#' @examples
#' n <- 25
#' data("cell_growth_data") 
#' select_est <- idpSTMSelect(cell_growth_data, n, maxit = 30, ModelCnt = 20, Message = FALSE)
#' print("======= Selected model =======")
#' print(select_est$selected_model)
#' print("======= Estimates =======")
#' print(select_est$coefficients)
#' print("======= Standard errors =======")
#' print(select_est$standard_error)

idpSTMSelect <- function(data, n, maxit, ModelCnt, Message = F){
  
  dat <- tilling(data, n)
  res <- idpSTModelSelection_cpp(dat, n, maxit, ModelCnt, Message)
  
  ret_est <- list("intercept" = res$parameters[1, ], "main_effects" = res$parameters[-1, ])
  ret_se <- list("intercept" = sqrt(res$se[1, ]), "main_effects" = sqrt(res$se[-1, ]))
  
  return(list("coefficients" = ret_est, "likelihood" = res$likelihood, "standard_error" = ret_se, 
              "selected_model" = res$v))
}


#' Model selsction of log-GLM
#'
#' \code{logGLMselect} performs variable selection on log-GLM using BIC via Gibbs samplling. 
#'
#' @param y A vector of count data response.  
#' @param x A matrix of covariates. 
#' @param maxit Maximum number of iterations of maximum likelihood estimation. 
#' @param skip A vector of indices of predictors to be forced in the selected model.
#' @param ModelCnt The number of models to be generated via Gibbs samplling.  
#' @param Message Logical, if TRUE, prints the top 5 most frequently generated models. (Suggest on)
#'
#' @return A list with components
#' \itemize{
#' \item \code{likelihood:}   Maximized log-likelihood. 
#' \item \code{coefficients:} A vector of parameter estimates (the first one being intercept).
#' \item \code{selected_model:}  A binary vector with 1 indicating selected variables and 0 otherwise.
#' \item \code{standard_error:}  A vector of estimated standard errors. 
#' }
#'
#' @examples
#' set.seed(777)  
#' p <- 10      # dimension
#' nn <- 3000   # sample size
#' 
#' x <- matrix(0, nn, p)
#' for(c in 1:p){
#'   x[, c] <- rnorm(nn, (p/5))
#' }
#' intercept <- 5
#' true_beta <- sample(c(-0.2, 0, 0.1), p, replace = TRUE)
#' y <- x %*% true_beta + intercept 
#' y <- sapply(y, exp)
#' y <- sapply(y, rpois, n = 1)
#' 
#' obj <- logGLMselect(y, x, skip = NULL, ModelCnt = 500, Message = TRUE)
#' cat("true model:", c(1, as.numeric(as.logical(true_beta))), "\n")

logGLMselect <- function(y, x, maxit = 50, skip = NULL, ModelCnt = 20, Message = T){
  
  if(nrow(x) != length(y)){ stop(" x must have the same number of rows as y")}
  if(!length(skip)){ skip = 0 }else{ 
    if(range(skip)[1] < 1 || range(skip)[2] > ncol(x)) stop("Indices in skip are out of range")
    skip = skip - 1} 
  if(length(unique(x[, 1])) > 1) x <- cbind(rep(1, nrow(x)), x)
  
  res <- logGLMselect_cpp(y, x, maxit, skip, ModelCnt, Message)
  
  ret_se <- rep(0, ncol(x))
  ret_se[which(res$v == 1)] <- sqrt(res$se)
 
  return(list("likelihood" = res$likelihood, "coefficients" = res$parameters, 
               "selected_model" = res$v, "standard_erro" = ret_se))
}




