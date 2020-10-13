#' @title MLE of Laupacis' type NNT in Exponential dist.
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the MLE of the Laupacis' type NNT in Exponential distribution.
#' @param treat a numeric vector of the treatment arm results
#' @param control a numeric vector of the control arm results
#' @param tau a scalar that indicates the MCID
#' @param yt_bar mean value of the treatment arm vector
#' @param yc_bar mean value of the control arm vector
#' @param p_c sample proportion of success in the control arm
#' @param p_t sample proportion of success in the treatment arm
#' @param p_t.boot BS estimator of the sample proportion of success in the treatment arm
#' @param p_c.boot BS estimator of the sample proportion of success in the control arm
#' @param n_c number of observations in the control arm
#' @param n_t number of observations in the treatment arm
#' @return MLE of the Laupacis' type NNT for Exponential distribution
#' @keywords internal

##############################
### NNT MLE EXPON INCREASE ###
##############################

mle_exp_inc = function( treat, control, tau, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) {

  nnt.v      = ifelse( ( exp( - yt_bar ^ (-1) * tau ) - exp( - yc_bar ^ (-1) * tau ) ) ^ (-1) > 1,
                       ( exp( - yt_bar ^ (-1) * tau ) - exp( - yc_bar ^ (-1) * tau ) ) ^ (-1),
                       Inf )

  nnt.v.bs   = ifelse( ( exp( - p_t.boot$t[ ,3] ^ (-1) * tau ) - exp( - p_c.boot$t[ ,3] ^ (-1) * tau ) ) ^ (-1) > 0,
                       ( exp( - p_t.boot$t[ ,3] ^ (-1) * tau ) - exp( - p_c.boot$t[ ,3] ^ (-1) * tau ) ) ^ (-1),
                       Inf )

  # DELTA's CI
  # gradient of nnt.v

  grad.nnt       = c(   nnt.v ^ 2 * exp( - ( yt_bar ) ^ (-1) * tau ) * tau,
                      - nnt.v ^ 2 * exp( - ( yc_bar ) ^ (-1) * tau ) * tau )

  # inverse Fisher inf. matrix of lambda_t, lambda_c

  inv_fisher     =   diag( c( ( yt_bar ) ^ (-2),
                              ( yc_bar ) ^ (-2)  ) , 2 )

  # variance of nnt.v
  var_nnt.v      = t( grad.nnt ) %*% inv_fisher %*% grad.nnt * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.v delta CI
  ci.d.mle       = c( max( nnt.v - 1.96 * sqrt( var_nnt.v ), 1),
                      nnt.v + 1.96 * sqrt( var_nnt.v )   )

  # BS CI

  ci.bs          =  c( max( quantile(nnt.v.bs, .025), 1), quantile(nnt.v.bs, .975) )

  output           = cbind( nnt.v, t( ci.d.mle ), t( ci.bs ) )

  colnames(output) = c( "NNT MLE",
                        "CI DELTA L", "CI DELTA R",
                        "CI BS L",    "CI BS R" )

  return( output )
}

##############################
### NNT MLE EXPON DECREASE ###
##############################

mle_exp_dec = function( treat, control, tau, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) {

  nnt.v      = ifelse( ( - exp( - yt_bar ^ (-1) * tau ) + exp( - yc_bar ^ (-1) * tau ) ) ^ (-1) > 1,
                       ( - exp( - yt_bar ^ (-1) * tau ) + exp( - yc_bar ^ (-1) * tau ) ) ^ (-1),
                       Inf )

  nnt.v.bs   = ifelse( ( - exp( - p_t.boot$t[ ,3] ^ (-1) * tau ) + exp( - p_c.boot$t[ ,3] ^ (-1) * tau ) ) ^ (-1) > 0,
                       ( - exp( - p_t.boot$t[ ,3] ^ (-1) * tau ) + exp( - p_c.boot$t[ ,3] ^ (-1) * tau ) ) ^ (-1),
                       Inf )

  # DELTA's CI
  # gradient of nnt.v

  grad.nnt       = c(   nnt.v ^ 2 * exp( - ( yt_bar ) ^ (-1) * tau ) * tau,
                      - nnt.v ^ 2 * exp( - ( yc_bar ) ^ (-1) * tau ) * tau )

  # inverse Fisher inf. matrix of lambda_t, lambda_c

  inv_fisher     =   diag( c( ( yt_bar ) ^ (-2),
                              ( yc_bar ) ^ (-2)  ) , 2 )

  # variance of nnt.v
  var_nnt.v      = t( grad.nnt ) %*% inv_fisher %*% grad.nnt * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.v delta CI
  ci.d.mle       = c( max( nnt.v - 1.96 * sqrt( var_nnt.v ), 1),
                      nnt.v + 1.96 * sqrt( var_nnt.v )   )

  # BS CI

  ci.bs          =  c( max( quantile(nnt.v.bs, .025), 1), quantile(nnt.v.bs, .975) )

  output           = cbind( nnt.v, t( ci.d.mle ), t( ci.bs ) )

  colnames(output) = c( "NNT MLE",
                        "CI DELTA L", "CI DELTA R",
                        "CI BS L",    "CI BS R" )

  return( output )
}
