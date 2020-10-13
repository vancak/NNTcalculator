#' @title Laupacis' NNT FL type estimator
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results and returns the Furukawa & Leucht's NNT.
#' @param treat a numeric vector of the treatment arm results
#' @param control a numeric vector of the control arm results
#' @param tau a scalar that indicates the MCID
#' @param d Cohen's delta MLE
#' @param d.bs Cohen's delta BS of the MLE estimator
#' @param p_c sample proportion of success in the control arm
#' @param p_t sample proportion of success in the treatment arm
#' @param p_t.boot BS estimator of the sample proportion of success in the treatment arm
#' @param p_c.boot BS estimator of the sample proportion of success in the control arm
#' @param n_c number of observations in the control arm
#' @param n_t number of observations in the treatment arm
#' @keywords internal
#' @return Laupacis' NNT using Furukawa & Leucht's estimator


################################
### FURUKAWA EQ VAR INCREASE ###
################################

furukawa_inc = function( treat, control, tau, d, d.bs, p_c, p_t.boot, p_c.boot, n_c, n_t ) {

  # point est.

  nnt.f            = ifelse( 1 / ( pnorm(  d   + qnorm(p_c) )  - p_c ) > 0,
                                       1 / ( pnorm(  d   + qnorm(p_c) )  - p_c ),
                                       Inf )

  nnt.f.bs         = ifelse( 1 / ( pnorm(  d.bs + qnorm(p_c.boot$t[ ,1]) )  - p_c.boot$t[ ,1] ) > 0,
                                       1 / ( pnorm(  d.bs + qnorm(p_c.boot$t[ ,1]) )  - p_c.boot$t[ ,1] ),
                                       Inf )

  # BS CI

  ci.bs            =  c( max( quantile(nnt.f.bs, .025), 1), quantile(nnt.f.bs, .975) )

  output           =  cbind( nnt.f, t( ci.bs ) )

  colnames(output) = c( "NNT FL",
                        "CI BS L",    "CI BS R" )

  return( output )
}

################################
### FURUKAWA EQ VAR DECREASE ###
################################


furukawa_dec = function( treat, control, tau, d, d.bs, p_c, p_t.boot, p_c.boot, n_c, n_t ) {

  # point est.

  nnt.f            = ifelse( 1 / ( pnorm( - d   + qnorm(p_c) )  - p_c ) > 0,
                                       1 / ( pnorm( - d   + qnorm(p_c) )  - p_c ),
                                       Inf )

  nnt.f.bs         = ifelse( 1 / ( pnorm( - d.bs + qnorm(p_c.boot$t[ ,1]) )  - p_c.boot$t[ ,1] ) > 0,
                                       1 / ( pnorm( - d.bs + qnorm(p_c.boot$t[ ,1]) )  - p_c.boot$t[ ,1] ),
                                       Inf )

  # BS CI

  ci.bs            =  c( max( quantile(nnt.f.bs, .025), 1), quantile(nnt.f.bs, .975) )

  output           =  cbind( nnt.f, t( ci.bs ) )

  colnames(output) = c( "NNT FL",
                        "CI BS L",    "CI BS R" )

  return( output )
}
