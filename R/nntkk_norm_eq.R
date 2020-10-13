#' @title MLE of the KK-NNT in Normal distribution with equal variances
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the MLE of the Kraemer & Kupfer's type NNT for Normal distribution with equal sd.
#' @param treat a numeric vector of the treatment arm results
#' @param control a numeric vector of the control arm results
#' @param d Cohen's delta MLE
#' @param s_ml pooled sd MLE
#' @param s_ml.bs BS of the pooled sd MLE
#' @param p_t.boot BS estimator of the sample proportion of success in the treatment arm
#' @param p_c.boot BS estimator of the sample proportion of success in the control arm
#' @param n_c number of observations in the control arm
#' @param n_t number of observations in the treatment arm
#' @keywords internal
#' @return MLE of the KK-NNT for the Normal distribution with equal sd

#######################################
### KK-NNT MLE NORM EQ VAR INCREASE ###
#######################################

nntkk_norm_eq_inc = function( treat, control, d, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_c, n_t ) {

  nnt.kk       = ifelse( ( 2 * pnorm( d / sqrt ( 2 ) ) - 1 ) ^ ( - 1 )  > 0,
                         ( 2 * pnorm( d / sqrt ( 2 ) ) - 1 ) ^ ( - 1 ),
                        Inf )

  d.bs         = ( p_t.boot$t[ ,3] - p_c.boot$t[ ,3] ) / s_ml.bs

  nnt.kk.bs    =   ifelse( ( 2 * pnorm( d.bs / sqrt ( 2 ) ) - 1 ) ^ ( - 1 )  > 0,
                           ( 2 * pnorm( d.bs / sqrt ( 2 ) ) - 1 ) ^ ( - 1 ),
                           Inf )

  # DELTA's CI

  phi_deriv   = - ( nnt.kk ) ^ 2 * exp( - d ^ 2 / 4 ) / sqrt( pi )

  # variance of nnt.kk

  var_nnt.kk   = ( n_t + n_c ) / ( 2 * n_t * n_c )  * ( phi_deriv ) ^ 2

  # nnt.kk BS CI

  ci.bs.nnt_k  = c( max( quantile(nnt.kk.bs, .025), 1), quantile(nnt.kk.bs, .975) )

  # nnt.kk function CI
  ci.fun.nnt_k = c( max( ( 2 * pnorm( ( d + 1.96 * sqrt( ( n_t + n_c ) / ( 2 * n_t * n_c ) ) * sqrt( ( 4 + d ^ 2 ) / 2 ) ) / sqrt( 2 ) )  - 1 ) ^ ( - 1 ), 1),
                         ( 2 * pnorm( ( d - 1.96 * sqrt( ( n_t + n_c ) / ( 2 * n_t * n_c ) ) * sqrt( ( 4 + d ^ 2 ) / 2 ) ) / sqrt( 2 ) )  - 1 ) ^ ( - 1 ) )

  # nnt.kk delta CI
  ci.d.nnt_k   = c( max( nnt.kk - 1.96 * sqrt( var_nnt.kk ) * sqrt( ( 4 + d ^ 2 ) / 2 ), 1),
                    nnt.kk + 1.96 * sqrt( var_nnt.kk ) * sqrt( ( 4 + d ^ 2 ) / 2 ))

  output           = cbind( nnt.kk, t( ci.fun.nnt_k ), t( ci.d.nnt_k ), t( ci.bs.nnt_k ) )

  colnames(output) = c( "NNT KK",
                        "CI COHEN L", "CI COHEN R",
                        "CI DELTA L", "CI DELTA R",
                        "CI BS L",    "CI BS R" )
  return( output )
}

#######################################
### KK-NNT MLE NORM EQ VAR DECREASE ###
#######################################

nntkk_norm_eq_dec = function( treat, control, d, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_c, n_t ) {

  nnt.kk       = ifelse( ( 1 - 2 * pnorm( d / sqrt ( 2 ) ) ) ^ ( - 1 )  > 0,
                         ( 1 - 2 * pnorm( d / sqrt ( 2 ) ) ) ^ ( - 1 ),
                         Inf )

  d.bs         = ( p_t.boot$t[ ,3] - p_c.boot$t[ ,3] ) / s_ml.bs

  nnt.kk.bs    =   ifelse( ( 1 - 2 * pnorm( d.bs / sqrt ( 2 ) ) ) ^ ( - 1 )  > 0,
                           ( 1 - 2 * pnorm( d.bs / sqrt ( 2 ) ) ) ^ ( - 1 ),
                           Inf )

  # DELTA's CI

  phi_deriv   = - ( nnt.kk ) ^ 2 * exp( - d ^ 2 / 4 ) / sqrt( pi )

  # variance of nnt.kk

  var_nnt.kk   = ( n_t + n_c ) / ( 2 * n_t * n_c )  * ( phi_deriv ) ^ 2

  # nnt.kk BS CI

  ci.bs.nnt_k  = c( max( quantile(nnt.kk.bs, .025), 1), quantile(nnt.kk.bs, .975) )

  # nnt.kk function CI
  ci.fun.nnt_k = c( max( - ( 2 * pnorm( ( d - 1.96 * sqrt( ( n_t + n_c ) / ( 2 * n_t * n_c ) ) * sqrt( ( 4 + d ^ 2 ) / 2 ) ) / sqrt( 2 ) )  - 1 ) ^ ( - 1 ),  1),
                         - ( 2 * pnorm( ( d + 1.96 * sqrt( ( n_t + n_c ) / ( 2 * n_t * n_c ) ) * sqrt( ( 4 + d ^ 2 ) / 2 ) ) / sqrt( 2 ) )  - 1 ) ^ ( - 1 )  )

  # nnt.kk delta CI
  ci.d.nnt_k   = c( max( nnt.kk - 1.96 * sqrt( var_nnt.kk ) * sqrt( ( 4 + d ^ 2 ) / 2 ), 1),
                         nnt.kk + 1.96 * sqrt( var_nnt.kk ) * sqrt( ( 4 + d ^ 2 ) / 2 ) )

  output           = cbind( nnt.kk, t( ci.fun.nnt_k ), t( ci.d.nnt_k ), t( ci.bs.nnt_k ) )

  colnames(output) = c( "NNT KK",
                        "CI COHEN L", "CI COHEN R",
                        "CI DELTA L", "CI DELTA R",
                        "CI BS L",    "CI BS R" )

  return( output )
}
