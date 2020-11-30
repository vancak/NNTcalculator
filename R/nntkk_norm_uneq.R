#' @title MLE of the KK-NNT in Normal distribution with unequal sd
#'
#' @description Internal function. Not for users. Takes two numeric vectors of control and treatment results
#'  and returns the MLE of the Kraemer & Kupfer's type NNT for Normal distribution with unequal sd.
#' @param s_t sd of the treatment arm vector
#' @param s_c sd of the control arm vector
#' @param p_t.boot BS estimator of the sample proportion of success in the treatment arm
#' @param p_c.boot BS estimator of the sample proportion of success in the control arm
#' @param n_c number of observations in the control arm
#' @param n_t number of observations in the treatment arm
#' @param yt_bar mean value of the treatment arm vector
#' @param yc_bar mean value of the control arm vector
#' @keywords internal
#' @return MLE of the KK-NNT for the Normal distribution with unequal sd

#########################################
### KK-NNT MLE NORM UNEQ VAR INCREASE ###
#########################################

nntkk_norm_uneq_inc = function( yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, n_c, n_t ) {

  # point est.

  nnt.kk       = ifelse( ( 2 * pnorm( ( yt_bar  - yc_bar ) / sqrt( s_t ^ 2 + s_c ^ 2 ) ) - 1 ) ^ ( - 1 ) > 0,
                         ( 2 * pnorm( ( yt_bar  - yc_bar ) / sqrt( s_t ^ 2 + s_c ^ 2 ) ) - 1 ) ^ ( - 1 ),
                         Inf )

  nnt.kk.bs    =  ifelse( ( 2 * pnorm( ( p_t.boot$t[ ,3]  - p_c.boot$t[ ,3] ) / sqrt ( p_c.boot$t[ ,2] + p_t.boot$t[ ,2]  ) ) - 1 ) ^ ( - 1 )  > 0,
                          ( 2 * pnorm( ( p_t.boot$t[ ,3]  - p_c.boot$t[ ,3] ) / sqrt ( p_c.boot$t[ ,2] + p_t.boot$t[ ,2]  ) ) - 1 ) ^ ( - 1 ),
                          Inf )

  # DELTA's CI

  grad.exp          = sqrt( 2 ) * nnt.kk ^ 2 * dnorm( ( yt_bar  - yc_bar ) / sqrt( s_t ^ 2 + s_c ^ 2 ) ) / sqrt( s_t ^ 2 + s_c ^ 2 ) * c( - 1,
                                                                                                                                  1,
                                                                                                                                  (yt_bar - yc_bar) * s_t / ( s_t ^ 2 + s_c ^ 2 ),
                                                                                                                                  (yt_bar - yc_bar) * s_c / ( s_t ^ 2 + s_c ^ 2 ) )

  # inverse Fisher matrix

  inv_fisher.exp   = diag( c( s_t ^ 2, s_c ^ 2, s_t ^ 2 / 2, s_c ^ 2 / 2  ) , 4 )

  # variance of nnt.kk parametric

  var_nnt.kk       = t( grad.exp ) %*% inv_fisher.exp %*% grad.exp * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.kk par delta CI

  ci.d.nnt_kk      = c( max( nnt.kk - 1.96 * sqrt( var_nnt.kk ), 1),
                        nnt.kk + 1.96 * sqrt( var_nnt.kk ) )

  # nnt.kk BS CI

  ci.bs.nnt_k  = c( max( quantile(nnt.kk.bs, .025), 1), quantile(nnt.kk.bs, .975) )

  output           = cbind( nnt.kk, t( ci.d.nnt_kk ), t( ci.bs.nnt_k ) )

  colnames(output) = c( "KK-NNT",
                        "CI DL L", "CI DL U",
                        "CI NBS L","CI NBS U" )

  return( output )
}

#########################################
### KK-NNT MLE NORM UNEQ VAR DECREASE ###
#########################################

nntkk_norm_uneq_dec = function( yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, n_c, n_t ) {

  # point est.

  nnt.kk       = ifelse( ( 2 * pnorm( - ( yt_bar  - yc_bar ) / sqrt( s_t ^ 2 + s_c ^ 2 ) ) - 1 ) ^ ( - 1 ) > 0,
                         ( 2 * pnorm( - ( yt_bar  - yc_bar ) / sqrt( s_t ^ 2 + s_c ^ 2 ) ) - 1 ) ^ ( - 1 ),
                         Inf )

  nnt.kk.bs    =  ifelse( ( 2 * pnorm( - ( p_t.boot$t[ ,3]  - p_c.boot$t[ ,3] ) / sqrt ( p_c.boot$t[ ,2] + p_t.boot$t[ ,2]  ) ) - 1 ) ^ ( - 1 )  > 0,
                          ( 2 * pnorm( - ( p_t.boot$t[ ,3]  - p_c.boot$t[ ,3] ) / sqrt ( p_c.boot$t[ ,2] + p_t.boot$t[ ,2]  ) ) - 1 ) ^ ( - 1 ),
                          Inf )

  # DELTA's CI

  grad.exp          = sqrt( 2 ) * nnt.kk ^ 2 * dnorm( ( yt_bar  - yc_bar ) / sqrt( s_t ^ 2 + s_c ^ 2 ) ) / sqrt( s_t ^ 2 + s_c ^ 2 ) * c( - 1,
                                                                                                                                  1,
                                                                                                                                  (yt_bar - yc_bar) * s_t / ( s_t ^ 2 + s_c ^ 2 ),
                                                                                                                                  (yt_bar - yc_bar) * s_c / ( s_t ^ 2 + s_c ^ 2 ) )

  # inverse Fisher matrix

  inv_fisher.exp   = diag( c( s_t ^ 2, s_c ^ 2, s_t ^ 2 / 2, s_c ^ 2 / 2  ) , 4 )

  # variance of nnt.kk parametric

  var_nnt.kk       = t( grad.exp ) %*% inv_fisher.exp %*% grad.exp * ( n_c + n_t ) / ( 2 * n_t * n_c )

  # nnt.kk par delta CI

  ci.d.nnt_kk      = c( max( nnt.kk - 1.96 * sqrt( var_nnt.kk ), 1),
                        nnt.kk + 1.96 * sqrt( var_nnt.kk ) )

  # nnt.kk BS CI

  ci.bs.nnt_k  = c( max( quantile(nnt.kk.bs, .025), 1), quantile(nnt.kk.bs, .975) )

  output           = cbind( nnt.kk, t( ci.d.nnt_kk ), t( ci.bs.nnt_k ) )

  colnames(output) = c( "KK-NNT",
                        "CI DL L", "CI DL U",
                        "CI NBS L", "CI NBS U" )

  return( output )
}
