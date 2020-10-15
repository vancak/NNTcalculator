#' @title Kraemer & Kupfer's NNT calculator
#'
#' @description Calculates the parametric and the non-parametric KK-NNT.
#' Takes two numeric vectors that are the outcomes of the treatment and control arms,
#' and returns the estimated KK-NNT using the specified method.
#'
#' @param type           specification of the estimation method; 'mle' (for the Maximum Likelihood estimator)
#'  and 'non-param' (for the non-parametric MLE)
#' @param treat          vector of response variable of the treatment group
#' @param control        vector of response variable of the control group
#' @param decrease       logical TRUE or FALSE. Indicates whether the MCID change is decrease in the response variable
#' @param dist           distribution type (only for the 'mle' type estimator);
#'  "norm" (Normal), "expon" (Exponential). The default value is 'none'
#' @param equal.var      logical TRUE or FALSE; Indicates whether the variances are equal - for normal distribution only. The default value is TRUE
#' @import stats
#' @import boot
#' @return The estimated Kraemer & Kupfer's NNT and its confidence intervals using the specified estimation method.
#' @references Kraemer, H. C., & Kupfer, D. J. (2006). Size of treatment effects and their importance to clinical research and practice. Biological psychiatry, 59(11), 990-996.
#' @references Vancak, V., Goldberg, Y., & Levine, S. Z. (2020). Systematic analysis of the number needed to treat. Statistical Methods in Medical Research, 0962280219890635.
#' @examples
#' nnt_kk( type     = "non-param",
#'        treat     = rnorm(100, 100, 10),
#'        control   = rnorm(100, 110, 20),
#'        dist      = "none",
#'        decrease  = TRUE )
#' @export
nnt_kk <- function( type,           # estimator type; parametric or non-parametric
                    treat,          # vector of response variable in the treatment group
                    control,        # vector of response variable in the control group
                    decrease,       # whether the MCID change is decrease in the response variable (TRUE\FALSE)
                    dist = "none",  # distribution type (if parametric)
                    equal.var  ) {  # whether the variances are equal - only for normal dist.

  if( !(type %in% c("mle", "non-param") ) ) warning("type can be 'non-param' or 'mle' only")
  if( !is.numeric(treat) )                       warning("treatment arm vector (treat) must be a numeric vector")
  if( !is.numeric(control) )                     warning("control arm vector (treat) must be a numeric vector")
  if( !(decrease %in% c(T, F, TRUE, FALSE) ) )   warning("decrease must be TRUE or FALSE")
  if( type == "mle" & !(dist %in% c("normal", "expon") ) ) warning("distribution (dist) can be 'expon' (Exponential) or 'normal' only")
#  if( type == "mle" & dist == "normal" & !(equal.var %in% c(T, F, TRUE, FALSE) ) )  warning("equal.var must be TRUE or FALSE")
  if( type == "non-param" ) message("if 'dist' and/or 'equal.var' arguments were set they are ignored")

  # remove NAs
  treat   <- treat[   !is.na( treat )  ]
  control <- control[ !is.na( control )]

  # initial values
  n_c     = length( control )
  n_t     = length( treat )

  yc_bar  = mean( control )
  yt_bar  = mean( treat )

  s_ml    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * var( control ) + (n_t - 1) * var( treat ) ) ) ^ ( 1 / 2 )

  d       = ( yt_bar - yc_bar ) / s_ml

  s_t     = ( ( n_t - 1 ) * var( treat ) /  n_t ) ^ ( 1 / 2 )

  s_c     = ( ( n_c - 1 ) * var( control ) / n_c ) ^ ( 1 / 2 )

  ### BOOTSTRAP FUNCTIONS (INCREASE) ###

  p_c1   = function(data, indices) {
    p_tc    = sum(  outer(treat, data[indices], ">") ) / ( length( data[indices] ) * n_t )
    var.pc  = var(  data[indices] )
    mean.pc = mean( data[indices] )
    return( c(NA, var.pc, mean.pc, p_tc) )
  }

  p_t1 = function(data, indices) {
    var.pt  = var( data[indices] )
    mean.pt = mean( data[indices] )
    return( c(NA, var.pt, mean.pt, NA) )
  }

  ############################
  ### NON-PARAM - ANY DIST ###
  ############################

  ### INCREASE ###

  if( type == "non-param" & decrease == F & ( dist == "none" | missing(dist) ) & missing(equal.var) ) {

    # prob of interest
    p_tk     = sum( outer( treat, control, FUN = ">" ) ) / ( n_t * n_c )

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)

    return( nntkk_nonparam(  treat, control, p_c.boot, p_tk ) )

  }

  ### DECREASE ###

  if( type == "non-param" & decrease == T & ( dist == "none" | missing(dist) ) & missing(equal.var) ) {

    # prob of interest
    p_tk   = sum( outer( control, treat, FUN = ">" ) ) / ( n_t * n_c )

    ### BOOTSTRAP FUNCTIONS ###
    p_c1   = function(data, indices) {
      p_tc    = sum(  outer(control, data[indices], ">") ) / ( length(data[indices] ) * n_t )
      var.pc  = var(  data[indices] )
      mean.pc = mean( data[indices] )
      return( c(NA, var.pc, mean.pc, p_tc) )
    }

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = treat, statistic = p_c1, R = 1000)

    ### BOOTSTRAP ESTIMATORS ###

    return( nntkk_nonparam( treat, control, p_c.boot, p_tk ) )
  }

  #################
  ### MLE EXPON ###
  #################

  ### INCREASE ###

  if( type == "mle" & decrease == F & dist == "expon" ) {

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    return( nntkk_exp_inc(  treat, control, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) )
  }

  ### DECREASE ###

  if( type == "mle" & decrease == T & dist == "expon" ) {

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    return( nntkk_exp_dec(  treat, control, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) )
  }

  #########################
  ### NORM EQ VAR - MLE ###
  #########################

  ### INCREASE ###

  if( type == "mle" & decrease == F & dist == "normal" & equal.var == T ) {

    ### BOOTSTRAP ESTIMATORS ###

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^ ( 1 / 2 )

    return( nntkk_norm_eq_inc(treat, control, d, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_c, n_t) )
  }

  ### DECREASE ###

  if( type == "mle" & decrease == T & dist == "normal" & equal.var == T ) {

    ### BOOTSTRAP ESTIMATORS ###

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^ ( 1 / 2 )

    return( nntkk_norm_eq_dec( treat, control, d, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_c, n_t ) )
  }

  #########################
  ### MLE NORM UNEQ VAR ###
  #########################

  ### INCREASE ###

  if( type == "mle" & decrease == F & dist == "normal" & equal.var == F ) {

    ### BOOTSTRAP ESTIMATORS ###

    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    return( nntkk_norm_uneq_inc( yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, n_c, n_t ) )
  }

  ### DECREASE ###

  if( type == "mle" & decrease == T & dist == "normal" & equal.var == F ) {

    ### BOOTSTRAP FUNCTIONS ###

    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    return( nntkk_norm_uneq_dec( yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, n_c, n_t )  )
  }
}

# set.seed(1)
#
#  nnt_kk(  type     = "non-param",
#          treat     = rnorm(100, 100, 10),
#          control   = rnorm(100, 110, 20),
#          dist      = "expon",
#          equal.var = T,
#          decrease  = T )
