
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
#' @seealso \code{\link{copSTM_fit}}, \code{\link{copSTModelSelect}}

make_data <- function(data, n){
  dat = tilling(data, n)
  return(data_cor(dat, n))
}



#' Simulate spatio-temporal data
#'
#' Simulate data for testing \code{copSTM} and \code{copSTModelSelect}.
#'
#' @param y_ini An integer initial value for the first time point.
#' @param n An integer spatial parameter. Data is simulated as n*n grids. 
#' @param K The number of groups. 
#' @param t_size The number time points to be generated.
#' @param beta  True values of regression parameters. 
#' @param rho_v  True values of correlations. 
#'
#' @return Data ready to input to \code{copSTM} or \code{copSTModelSelect}.
#' 
#' @seealso \code{\link{copSTM_fit}}, \code{\link{copSTModelSelect}}

sim_data <- function(y_ini, n, K, t_size, beta, rho_v){
  return(sim_data_cpp(y_ini, n, K, t_size, beta, rho_v))
}