#' @title MLE of NNT in unequal sd Normal distribution
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the MLE of the Laupacis' type NNT for Normal distribution with unequal sd.
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
#' @param s_t sd of the variance MLE in the treatment arm vector
#' @param s_c sd of the variance MLE in the control arm vector
#' @param s_t.bs BS estimator of the sd of the variance MLE in the treatment arm vector
#' @param s_c.bs BS estimator of the sd of the variance MLE in the control arm vector
#' @keywords internal
#' @return MLE of the Laupacis' type NNT for Normal distribution with unequal sd

##############################
### MLE NORM UNEQ INCREASE ###
##############################

# point est.

mle_norm_uneq_inc = function( treat, control, tau, yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, s_t.bs, s_c.bs, n_c, n_t ) {

  nnt.v      = ifelse( (    pnorm( ( yt_bar - tau ) / s_t )
                          - pnorm( ( yc_bar - tau ) / s_c )  ) ^ ( - 1 )  > 0,
                       (    pnorm( ( yt_bar - tau ) / s_t )
                          - pnorm( ( yc_bar - tau ) / s_c )  ) ^ ( - 1 ),
                       Inf )


  nnt.v.bs   = ifelse( (    pnorm( ( p_t.boot$t[ ,3] - tau ) / s_t.bs )
                          - pnorm( ( p_c.boot$t[ ,3] - tau ) / s_c.bs )  ) ^ ( - 1 ) > 0,
                       (    pnorm( ( p_t.boot$t[ ,3] - tau ) / s_t.bs )
                          - pnorm( ( p_c.boot$t[ ,3] - tau ) / s_c.bs )  ) ^ ( - 1 ),
                       Inf )

  # DELTA's CI
  # gradient of nnt.v

  grad.nnt       = nnt.v ^ 2 * c( dnorm( ( yt_bar - tau ) / s_t ) * 1 / s_t,
                                        - dnorm( ( yc_bar - tau ) / s_c ) * 1 / s_c,
                                          dnorm( ( yt_bar - tau ) / s_t ) * ( tau - yt_bar ) / s_t ^ 2,
                                          dnorm( ( yc_bar - tau ) / s_t ) * ( yc_bar - tau ) / s_c ^ 2 )

  # inverse Fisher inf. matrix of mu_t, mu_c and sigma_t, sigma_c

  inv_fisher     = diag( c( s_t ^ 2,
                            s_c ^ 2,
                            s_t ^ 2 / 2,
                            s_c ^ 2 / 2  ) , 4 )

  # variance of nnt.v
  var_nnt.v      = t( grad.nnt ) %*% inv_fisher %*% grad.nnt * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.v delta CI
  ci.d.mle   = c( max( nnt.v - 1.96 * sqrt( var_nnt.v ), 1),
                       nnt.v + 1.96 * sqrt( var_nnt.v )   )

  # BS CI

  ci.bs      =  c( max( quantile(nnt.v.bs, .025), 1), quantile(nnt.v.bs, .975) )

  output           = cbind( nnt.v, t( ci.d.mle ), t( ci.bs ) )

  colnames(output) = c( "NNT MLE",
                        "CI DELTA L", "CI DELTA R",
                        "CI BS L",    "CI BS R" )
  return( output )
}

##############################
### MLE NORM UNEQ DECREASE ###
##############################

mle_norm_uneq_dec = function( treat, control, tau, yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, s_t.bs, s_c.bs, n_c, n_t ) {

  nnt.v      = ifelse( (    pnorm( ( - yt_bar + tau ) / s_t )
                          - pnorm( ( - yc_bar + tau ) / s_c )  ) ^ ( - 1 )  > 0,
                       (    pnorm( ( - yt_bar + tau ) / s_t )
                          - pnorm( ( - yc_bar + tau ) / s_c )  ) ^ ( - 1 ),
                       Inf )

  nnt.v.bs   = ifelse( (    pnorm( ( - p_t.boot$t[ ,3] + tau ) / s_t.bs )
                          - pnorm( ( - p_c.boot$t[ ,3] + tau ) / s_c.bs )  ) ^ ( - 1 ) > 0,
                       (    pnorm( ( - p_t.boot$t[ ,3] + tau ) / s_t.bs )
                          - pnorm( ( - p_c.boot$t[ ,3] + tau ) / s_c.bs )  ) ^ ( - 1 ),
                       Inf )

  # DELTA's CI
  # gradient of nnt.v

  grad.nnt       = nnt.v ^ 2 * c( dnorm( ( yt_bar - tau ) / s_t ) * 1 / s_t,
                                         - dnorm( ( yc_bar - tau ) / s_c ) * 1 / s_c,
                                         dnorm( ( yt_bar - tau ) / s_t ) * ( tau - yt_bar ) / s_t ^ 2,
                                         dnorm( ( yc_bar - tau ) / s_t ) * ( yc_bar - tau ) / s_c ^ 2 )

  # inverse Fisher inf. matrix of mu_t, mu_c and sigma_t, sigma_c

  inv_fisher    = diag( c( s_t ^ 2,
                            s_c ^ 2,
                            s_t ^ 2 / 2,
                            s_c ^ 2 / 2  ) , 4 )

  # variance of nnt.v
  var_nnt.v  = t( grad.nnt ) %*% inv_fisher %*% grad.nnt * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.v delta CI
  ci.d.mle   = c( max( nnt.v - 1.96 * sqrt( var_nnt.v ), 1),
                  nnt.v + 1.96 * sqrt( var_nnt.v )   )

  # BS CI

  ci.bs      =  c( max( quantile(nnt.v.bs, .025), 1), quantile(nnt.v.bs, .975) )

  output     = cbind( nnt.v, t( ci.d.mle ), t( ci.bs ) )

  colnames(output) = c( "NNT MLE",
                        "CI DELTA L", "CI DELTA R",
                        "CI BS L",    "CI BS R" )
  return( output )
}
