#' @title MLE of KK-NNT in Exponential distribution
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the MLE of the Kraemer & Kupfer's type NNT for Exponential distribution.
#' @param treat a numeric vector of the treatment arm results
#' @param control a numeric vector of the control arm results
#' @param yt_bar mean value of the treatment arm vector
#' @param yc_bar mean value of the control arm vector
#' @param p_t.boot BS estimator of the sample proportion of success in the treatment arm
#' @param p_c.boot BS estimator of the sample proportion of success in the control arm
#' @param n_c number of observations in the control arm
#' @param n_t number of observations in the treatment arm
#' @keywords internal
#' @return MLE of the KK-NNT for the Exponential distribution

##############################
### NNT KK EXPON INCREASE ###
##############################

nntkk_exp_inc = function( treat, control, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) {

  # point est.

  nnt.kk            = ifelse( ( 1 / yc_bar + 1 / yt_bar ) / ( ( 1 / yc_bar - 1 / yt_bar ) ) > 0,
                              ( 1 / yc_bar + 1 / yt_bar ) / ( ( 1 / yc_bar - 1 / yt_bar ) ),
                              Inf )

  nnt.kk.bs         =   ifelse( ( 1 / p_c.boot$t[, 3] + 1 / p_t.boot$t[, 3] ) /  ( 1 / p_c.boot$t[, 3] - 1 / p_t.boot$t[, 3] ) > 0,
                           ( 1 / p_c.boot$t[, 3] + 1 / p_t.boot$t[, 3] ) /  ( 1 / p_c.boot$t[, 3] - 1 / p_t.boot$t[, 3] ),
                           Inf )

  # DELTA's CI

  grad.exp          =  2  / ( 1 / yc_bar - 1 / yt_bar ) ^ 2 * c( 1 / yc_bar, - 1 / yt_bar )

  # inverse Fisher matrix

  inv_fisher.exp    = diag( c( ( yt_bar ) ^ ( - 2 ),
                               ( yc_bar ) ^ ( - 2 )  ) , 2 )

  # variance of nnt.kk parametric

  var_nnt.kk       = t( grad.exp ) %*% inv_fisher.exp %*% grad.exp * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.kk par delta CI

  ci.d.nnt_kk      = c( max( nnt.kk - 1.96 * sqrt( var_nnt.kk ), 1),
                        nnt.kk + 1.96 * sqrt( var_nnt.kk ) )
  # BS CI

  ci.bs          =  c( max( quantile(nnt.kk.bs, .025), 1), quantile(nnt.kk.bs, .975) )

  output           = cbind( nnt.kk, t( ci.d.nnt_kk ), t( ci.bs ) )

  colnames(output) = c( "KK-NNT MLE",
                        "CI DL L",     "CI DL U",
                        "CI NBS L",    "CI NBS U" )

  return( output )
}

##############################
### NNT KK EXPON DECREASE ###
##############################

nntkk_exp_dec = function( treat, control, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) {

  # point est.

  nnt.kk            = ifelse( ( 1 / yc_bar + 1 / yt_bar ) / ( ( - 1 / yc_bar + 1 / yt_bar ) ) > 0,
                              ( 1 / yc_bar + 1 / yt_bar ) / ( ( - 1 / yc_bar + 1 / yt_bar ) ),
                              Inf )

  nnt.kk.bs         =   ifelse( ( 1 / p_c.boot$t[, 3] + 1 / p_t.boot$t[, 3] ) /  ( - 1 / p_c.boot$t[, 3] + 1 / p_t.boot$t[, 3] ) > 0,
                                ( 1 / p_c.boot$t[, 3] + 1 / p_t.boot$t[, 3] ) /  ( - 1 / p_c.boot$t[, 3] + 1 / p_t.boot$t[, 3] ),
                                Inf )

  # DELTA's CI

  grad.exp          =  2  / ( 1 / yc_bar - 1 / yt_bar ) ^ 2 * c(- 1 / yc_bar, - 1 / yt_bar )

  # inverse Fisher matrix

  inv_fisher.exp    = diag( c( ( yt_bar ) ^ ( - 2 ),
                               ( yc_bar ) ^ ( - 2 )  ) , 2 )

  # variance of nnt.kk parametric

  var_nnt.kk       = t( grad.exp ) %*% inv_fisher.exp %*% grad.exp * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.kk par delta CI

  ci.d.nnt_kk      = c( max( nnt.kk - 1.96 * sqrt( var_nnt.kk ), 1),
                        nnt.kk + 1.96 * sqrt( var_nnt.kk ) )
  # BS CI

  ci.bs          =  c( max( quantile(nnt.kk.bs, .025), 1), quantile(nnt.kk.bs, .975) )

  output           = cbind( nnt.kk, t( ci.d.nnt_kk ), t( ci.bs ) )

  colnames(output) = c( "KK-NNT MLE",
                        "CI DL L", "CI DL U",
                        "CI NBS L", "CI NBS U" )

  return( output )
}
