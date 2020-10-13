#' @title Non-paramatric Laupacis' NNT
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the non-parametric MLE.
#' @param treat a numeric vector of the treatment arm results
#' @param control a numeric vector of the control arm results
#' @param p_c sample proportion of success in the control arm
#' @param p_t sample proportion of success in the treatment arm
#' @param p_t.boot BS estimator of the sample proportion of success in the treatment arm
#' @param p_c.boot BS estimator of the sample proportion of success in the control arm
#' @param n_c number of observations in the control arm
#' @param n_t number of observations in the treatment arm
#' @return Laupacis' type NNT non-parametric MLE
#' @keywords internal
####################
### LAUPACIS NNT ###
####################

laupacis <- function( treat, control, p_t, p_c, n_t, n_c, p_t.boot, p_c.boot ) {

  # point est.
  nntl        = ifelse( 1 / ( p_t - p_c ) > 0,
                        1 / ( p_t - p_c ),
                        Inf )

  nntl.bs     = ifelse( 1 / ( p_t.boot$t[ ,1] - p_c.boot$t[ ,1] ) > 0,
                        1 / ( p_t.boot$t[ ,1] - p_c.boot$t[ ,1] ),
                        Inf )

  # Wald's CI
  sd.wald       = sqrt( p_t * ( 1 - p_t ) / n_t + p_c * ( 1 - p_c ) / n_c )

  ci_w          =    c( max( 1 / (  p_t - p_c + qnorm(.975) * sd.wald ), 1),
                        ifelse( 1 / (  p_t - p_c - qnorm(.975) * sd.wald ) > 0,
                                1 / (  p_t - p_c - qnorm(.975) * sd.wald ),
                                Inf ) )

  # DELTA's CI
  sd.delta      = ( 1 / (p_t - p_c) ^ 2 ) * sqrt(   p_t * (1 - p_t) / n_t
                                                    + p_c * (1 - p_c) / n_c )
  ci_d          = c( max( nntl - qnorm(.975) * sd.delta, 1),
                     nntl + qnorm(.975) * sd.delta )


  ci_bs         = c( max(quantile(nntl.bs, .025), 1),  quantile(nntl.bs, .975) )

  output           = cbind(nntl, t( ci_w ), t( ci_d ) , t( ci_bs ) )

   colnames(output) = c( "NNTL",
                        "CI WALD L",  "CI WALD R",
                        "CI DELTA L", "CI DELTA R",
                        "CI BS L",    "CI BS R" )
  return( output )
 }

