
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
#' \code{make_data} takes a matrix of spatio-temporal data and returns a reorganized data ready for input of \code{copSTM} or \code{copSTModelSelect}.
#'
#' @param data A matrix with four cloumns: time, x, y, group, where x and y are coordinates of location. 
#' @param n An integer number indicating size of the grids. The image is tiled into an n*n grid. (Suggest n >= 6, otherwise neighbourhoods mostly overlapping with each other will lead to highly correlated covariates)
#'
#' @return A list with components: response, covariates, T (the number of time points) and K (the number of groups).
#'
#' @seealso \code{\link{copSTM}}, \code{\link{copSTMSelect}}

make_data <- function(data, n){
  dat = tilling(data, n)
  return(data_cor(dat, n))
}



#' Simulate spatio-temporal data
#'
#' Simulate data for testing \code{copSTM} and \code{copSTMSelect}.
#'
#' @param n An integer spatial parameter. Data is simulated as n*n grids. 
#' @param K The number of groups. 
#' @param t_size The number time points to be generated.
#' @param beta  True values of regression parameters. 
#' @param rho_v  True values of correlations. 
#' @param cor_type A character string, one of "sp", "mv" and "both". 
#' @param temporal Logical, if TRUE, generates response in temporal setting.
#'  Otherwise, each d-dimnesional $y_i (i = 1, \dots, t_size)$ are generated independently.  
#' @param y_ini An integer initial value for the first time point.
#'  Specifically, if "sp" (short for spatial), only spatial correlations of the same group from neighbouring tiles are captured; 
#'  if "mv" (short for multivariate), only the correlations between groups in the same tile are captured; 
#'  if "both", the model computes "sp" and "mv", as well as between-group correlation from neighbouring tiles.
#'
#' @return Data ready to input to \code{copSTM} or \code{copSTMSelect}.
#' 
#' @seealso \code{\link{copSTM}}, \code{\link{copSTMSelect}}

sim_data <- function(n, K, t_size, beta, rho_v, 
                     cor_type = 3, temporal = TRUE, y_ini = 10){
  if(cor_type == "sp") ct <- 1
  if(cor_type == "mv") ct <- 2 
  if(cor_type == "both") ct <- 3 
  return(sim_data_cpp(n, K, temporal, t_size, beta, rho_v, ct, y_ini))
}