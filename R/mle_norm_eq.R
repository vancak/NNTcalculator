#' @title MLE of NNT in equal sd Normal distribution
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the MLE of the Laupacis' type NNT for Normal distribution with equal variances.
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
#' @param s_ml sd of the the pooled variance MLE
#' @param s_ml.bs BS estimator of the sd of the the pooled variance MLE
#' @return MLE of the Laupacis' type NNT for Normal distribution with equal sd
#' @keywords internal
################################
### MLE NORM EQ VAR INCREASE ###
################################

mle_norm_eq_inc = function( treat, control, tau, yt_bar, yc_bar, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_c, n_t ) {

nnt.v      = ifelse( (  pnorm( ( yt_bar - tau ) / s_ml )
                      - pnorm( ( yc_bar - tau ) / s_ml )  ) ^ ( - 1 )  > 0,
                     (  pnorm( ( yt_bar - tau ) / s_ml )
                      - pnorm( ( yc_bar - tau ) / s_ml )  ) ^ ( - 1 ),
                     Inf )


nnt.v.bs   = ifelse( (  pnorm( ( p_t.boot$t[ ,3] - tau ) / s_ml.bs )
                      - pnorm( ( p_c.boot$t[ ,3] - tau ) / s_ml.bs )  ) ^ ( - 1 ) > 0,
                     (  pnorm( ( p_t.boot$t[ ,3] - tau ) / s_ml.bs )
                      - pnorm( ( p_c.boot$t[ ,3] - tau ) / s_ml.bs )  ) ^ ( - 1 ),
                     Inf )

# DELTA's CI
# gradient of nnt.v
grad.nnt   = c( - nnt.v ^ 2 / s_ml * dnorm( ( yt_bar - tau)/s_ml ),
                  nnt.v ^ 2 / s_ml * dnorm( ( yc_bar - tau)/s_ml ),
                  nnt.v ^ 2 / ( s_ml ^ 2 ) * ( dnorm( (yc_bar - tau)/s_ml ) * (yc_bar - tau)
                                             - dnorm( (yt_bar - tau)/s_ml ) * (yt_bar - tau) )  )

# inverse Fisher inf. matrix of mu_t, mu_c and sigma
inv_fisher = diag( c(     s_ml ^ 2,
                          s_ml ^ 2,
                      2 * s_ml ^ 2 / 2 ), 3 )

# variance of nnt.v
var_nnt.v  = t( grad.nnt ) %*% inv_fisher %*% grad.nnt * ( n_c + n_t ) / ( 2 * n_t * n_c )

# nnt.v delta CI
ci.d.mle   = c( max( nnt.v - 1.96 * sqrt( var_nnt.v ), 1),
                     nnt.v + 1.96 * sqrt( var_nnt.v )   )

# BS CI
ci.bs      =  c( max( quantile(nnt.v.bs, .025), 1), quantile(nnt.v.bs, .975) )

output           = cbind( nnt.v, t( ci.d.mle ), t( ci.bs ) )

colnames(output) = c( "NNT MLE",
                      "CI DL L", "CI DL U",
                      "CI NBS L", "CI NBS U" )

return( output )
}

################################
### MLE NORM EQ VAR DECREASE ###
################################

mle_norm_eq_dec = function( treat, control, tau, yt_bar, yc_bar, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_c, n_t ) {

  nnt.v      = ifelse( (    pnorm( ( - yt_bar + tau ) / s_ml )
                          - pnorm( ( - yc_bar + tau ) / s_ml )  ) ^ ( - 1 )  > 0,
                       (    pnorm( ( - yt_bar + tau ) / s_ml )
                          - pnorm( ( - yc_bar + tau ) / s_ml )  ) ^ ( - 1 ),
                       Inf )


  nnt.v.bs   = ifelse( (    pnorm( ( - p_t.boot$t[ ,3] + tau ) / s_ml.bs )
                          - pnorm( ( - p_c.boot$t[ ,3] + tau ) / s_ml.bs )  ) ^ ( - 1 ) > 0,
                       (    pnorm( ( - p_t.boot$t[ ,3] + tau ) / s_ml.bs )
                          - pnorm( ( - p_c.boot$t[ ,3] + tau ) / s_ml.bs )  ) ^ ( - 1 ),
                       Inf )

  # DELTA's CI
  # gradient of nnt.v
  grad.nnt   = c( - nnt.v ^ 2 / s_ml * dnorm( ( yt_bar - tau)/s_ml ),
                    nnt.v ^ 2 / s_ml * dnorm( ( yc_bar - tau)/s_ml ),
                    nnt.v ^ 2 / ( s_ml ^ 2 ) * ( dnorm( (yc_bar - tau)/s_ml ) * (yc_bar - tau)
                                               - dnorm( (yt_bar - tau)/s_ml ) * (yt_bar - tau) )  )

  # inverse Fisher inf. matrix of mu_t, mu_c and sigma
  inv_fisher = diag( c(     s_ml ^ 2,
                            s_ml ^ 2,
                            2 * s_ml ^ 2 / 2 ), 3 )

  # variance of nnt.v
  var_nnt.v  = t( grad.nnt ) %*% inv_fisher %*% grad.nnt * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.v delta CI
  ci.d.mle   = c( max( nnt.v - 1.96 * sqrt( var_nnt.v ), 1),
                  nnt.v + 1.96 * sqrt( var_nnt.v )   )

  # BS CI
  ci.bs      =  c( max( quantile(nnt.v.bs, .025), 1), quantile(nnt.v.bs, .975) )

  output           = cbind( nnt.v, t( ci.d.mle ), t( ci.bs ) )

  colnames(output) = c( "NNT MLE",
                        "CI DL L",    "CI DL  U",
                        "CI NBS L",   "CI NBS U" )

  return( output )
}
