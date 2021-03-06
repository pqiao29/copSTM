# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

bd_unary <- function(xx_row, y_i, upper, beta, marginal, dispersion) {
    .Call(`_copSTM_bd_unary`, xx_row, y_i, upper, beta, marginal, dispersion)
}

bd_deriv_unary_arma <- function(xx_row, y_i, upper, beta, marginal, dispersion) {
    .Call(`_copSTM_bd_deriv_unary_arma`, xx_row, y_i, upper, beta, marginal, dispersion)
}

copSTModelSelect_cpp <- function(x, y, cor_type, K, n, marginal, temporal, ModelCnt, B, maxit1, maxit2, add_penalty = 0, Message_prog = TRUE, Message_res = TRUE, eps = 0.1) {
    .Call(`_copSTM_copSTModelSelect_cpp`, x, y, cor_type, K, n, marginal, temporal, ModelCnt, B, maxit1, maxit2, add_penalty, Message_prog, Message_res, eps)
}

copSTM_cpp <- function(x, y, marginal, temporal, cor_type, K, n, maxit, eps, std_err, B = 0L, Message_prog = FALSE) {
    .Call(`_copSTM_copSTM_cpp`, x, y, marginal, temporal, cor_type, K, n, maxit, eps, std_err, B, Message_prog)
}

logGLMselect_cpp <- function(y, x, marginal, maxit, skip, ModelCnt, Message) {
    .Call(`_copSTM_logGLMselect_cpp`, y, x, marginal, maxit, skip, ModelCnt, Message)
}

idpSTModelSelection_cpp <- function(dat, n, marginal, maxit, ModelCnt, Message) {
    .Call(`_copSTM_idpSTModelSelection_cpp`, dat, n, marginal, maxit, ModelCnt, Message)
}

idptSTM_cpp <- function(dat, n, marginal, maxit, fit_plot) {
    .Call(`_copSTM_idptSTM_cpp`, dat, n, marginal, maxit, fit_plot)
}

make_cor_label <- function(K, n, p_rho, cor_type) {
    .Call(`_copSTM_make_cor_label`, K, n, p_rho, cor_type)
}

data_cor <- function(dat, n) {
    .Call(`_copSTM_data_cor`, dat, n)
}

sim_data_cpp <- function(n, K, temporal, t_size, beta, rho_v, cor_type, y_ini, marginal, dispersion) {
    .Call(`_copSTM_sim_data_cpp`, n, K, temporal, t_size, beta, rho_v, cor_type, y_ini, marginal, dispersion)
}

