#' @title Non-Parametric MLE of KK-NNT in Exponential distribution
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the non-parametric MLE of the Kraemer & Kupfer's NNT.
#' @param treat a numeric vector of the treatment arm results
#' @param control a numeric vector of the control arm results
#' @param p_c.boot BS estimator of the sample proportion of success in the control arm
#' @param p_tk sample estimator of P(Y_t > Y_c)
#' @keywords internal
#' @return Non-Paramaetric MLE of the KK-NNT

########################
### NON-PARAM KK-NNT ###
########################

nntkk_nonparam<- function( treat, control, p_c.boot, p_tk ) {

  nnt.kkn         = ifelse( ( 2 * p_tk - 1 ) ^ ( - 1 ) > 0,
                            ( 2 * p_tk - 1 ) ^ ( - 1 ),
                            Inf )

  # point est.
  nnt.kkn.bs   = ifelse( ( 2 * p_c.boot$t[ ,4] - 1 ) ^  ( - 1 ) > 0,
                         ( 2 * p_c.boot$t[ ,4] - 1 ) ^  ( - 1 ),
                         Inf )

  ci_bs         = c( max(quantile(nnt.kkn.bs, .025), 1),  quantile(nnt.kkn.bs, .975) )

  output           = cbind( nnt.kkn ,t( ci_bs ) )

  colnames(output) = c( "KK-NNT",
                        "CI BS L",  "CI BS R" )
  return( output )
}

