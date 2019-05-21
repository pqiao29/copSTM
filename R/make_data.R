
tilling <- function(data, n){
  cordi_x <- data[, 2]
  cordi_y <- data[, 3]
  
  Mx <- max(cordi_x) + 1.e-8
  mx <- min(cordi_x) - 1.e-8
  My <- max(cordi_y) + 1.e-8
  my <- min(cordi_y) - 1.e-8
  l.x <- (Mx - mx)/n
  l.y <- (My - my)/n
  
  loc_x <- findInterval(cordi_x, seq(mx, Mx, l.x))
  loc_y <- findInterval(cordi_y, seq(my, My, l.y))
  tmp_tile  <- (loc_y - 1) * n + loc_x
  
  cbind(data[, c(1, 4)], tmp_tile)
}


#' Data organization
#'
#' \code{make_data} takes a matrix of grouped spatio-temporal data and returns a reorganized data
#' in the form of response and covariates, 
#' ready for input of \code{copSTM} or \code{copSTModelSelect}.
#'
#' @param data Same as the first input in \code{\link{idpSTM}} and \code{\link{idpSTMSelect}}. 
#' A matrix with four cloumns (in order): time, x, y, group, where x and y are coordinates of location.
#' See the cell growth data as an example.
#' @param n An integer number indicating size of the grids. The image is tiled into an n*n grid. (Suggest n >= 6, otherwise neighbourhoods mostly overlapping with each other will lead to highly correlated covariates)
#'
#' @return A list with components: response, covariates, T (the number of time points) and K (the number of groups).
#'
#' @seealso \code{\link{copSTM}}, \code{\link{copSTMSelect}}

make_data <- function(data, n){
  dat = tilling(data, n)
  return(data_cor(dat, n))
}



#' Simulate grouped spatio-temporal count data
#'
#' Simulate correlated count data for testing \code{copSTM} and \code{copSTMSelect}.
#' Generated data could be fit into any GLM-based models with logarithm link for regression coefficients. 
#' Fit into \code{copSTM} or \code{copSTMSelect} to get estimates of correlation parameters as well.
#' Data can be simulated under two settings (see argument \code{temporal}): 
#' 1. The temporal setting: Covariates of \eqn{y_t} is generated based on \eqn{y_{t - 1}},
#' where \eqn{t} denotes time point. 
#' Regression coefficients have the same interpretation as main effects in \code{\link{idpSTM}}.
#' 2. The non-temporal setting: \eqn{y_t} is considered independent for different t,  
#' and covariates are randomly generated from standard Gaussian distribution. 
#' For both settings, \eqn{y_t} is generated as a \eqn{K*n^2}-dimensional vector 
#' through a Gaussian copula 
#' with correlation structure specified in Arguments \code{cor_type} and \code{rho}.
#'
#' @param n An integer spatial parameter. Data is simulated as n*n grids. 
#' @param K The number of groups. 
#' @param t_size The number time points to be generated.
#' @param beta  True values of regression coefficients.
#' If temporal == TRUE, beta needs to have length K*(K + 1) 
#' (i.e. K intercepts for K groups and K*K impact parameters, see main_effects in values of \code{\link{idpSTM}} ).
#' If temporal == FALSE, beta can have any length greater than 1. 
#' @param marginal Distribution of response variable: "pois" for Poisson, "nbinom" for negative binomial. 
#' @param cor_type Correlation type. 
#'  \itemize{
#'  \item \code{sp}: Spatial only, consider only correlation of the same group from neighbouring tiles,
#'  \item \code{mv}: Multivariate only, consider only correlation between different groups in the same tile,
#'  \item \code{both}: Both "sp" and "mv", as well as between-group correlation from neighbouring tiles.
#'  }
#' @param rho  True values of correlations. 
#' If cor_type = "sp", needs to have length K, 
#' If cor_type = "mv", needs to have length K*(K + 1)/2, 
#' If cor_type = "both", needs to have length K + K*(K + 1), in the order of:  
#' ("sp", "mv", between-group correlation from neighbouring tiles).
#' @param temporal Logical, if TRUE, generates response in temporal setting.
#' @param dispersion Overdispersion parameter, useful only when marginal == "nbinom". 
#' Required to be strictly positive. 
#' @param y_ini An integer initial value for the first time point.
#'
#' @return A list of response and covariates separately. 
#'          
#' @seealso \code{\link{copSTM}}, \code{\link{copSTMSelect}}

sim_data <- function(n, K, t_size, beta, marginal, 
                     cor_type, rho, temporal, 
                     dispersion = 1, y_ini = 10){

  if(cor_type == "sp") ct <- 1
  if(cor_type == "mv") ct <- 2 
  if(cor_type == "both") ct <- 3 
  
  marginal <- which(c("pois", "nbinom") == marginal)
  
  return(sim_data_cpp(n, K, temporal, t_size, beta, rho, ct, y_ini, marginal, dispersion))
}
