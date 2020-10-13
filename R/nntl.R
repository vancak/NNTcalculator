#' @title Laupacis' unadjusted NNT calculator
#'
#' @description Calculates Laupacis' type NNT.
#' Takes two numeric vectors that are the outcomes of the treatment and control arms,
#' and returns the estimated NNT using the specified estimation method.
#' @param type           specification of the estimation method; 'mle' (for the Maximum Likelihood estimator),
#'  'fl' (for Furukawa & Leucht's estimator),
#'  'laupacis' (for the non-parametric MLE estimator)
#' @param treat          vector of response variable of the treatment group
#' @param control        vector of response variable of the control group
#' @param cutoff         a scalar that is the MCID
#' @param decrease       TRUE or FALSE. Indicates whether the MCID change is decrease in the response variable
#' @param dist           distribution type (if specified); "normal" (Normal), "expon" (Exponential).
#' The default value is 'none'.
#' @param equal.var      TRUE or FALSE; Indicates whether the variances are equal - for normal distribution only.
#' The default value is TRUE.
#' @return The estimated Laupacis' NNT and its confidence intervals using the specified estimation method.
#' @references Laupacis, A., Sackett, D. L., & Roberts, R. S. (1988). An assessment of clinically useful measures of the consequences of treatment. New England journal of medicine, 318(26), 1728-1733.
#' @references Vancak, V., Goldberg, Y., & Levine, S. Z. (2020). Systematic analysis of the number needed to treat. Statistical Methods in Medical Research, 0962280219890635.
#' @references Furukawa, T. A., & Leucht, S. (2011). How to obtain NNT from Cohen's d: comparison of two methods. PloS one, 6(4), e19070.
#' @examples
#'  nnt_l( type      = "mle",
#'         treat     = rnorm(1000, 110, 10),
#'         control   = rnorm(1000, 100, 10),
#'         cutoff    = 100,
#'         equal.var = TRUE,
#'         dist      = "normal",
#'         decrease  = FALSE )
#'
#' @export
#'
nnt_l <- function( type,           # estimator type, e.g., laupacis, mle, furukawa
                   treat,          # vector of response variable in the treament grp
                   control,        # vector of response variable in the control grp
                   cutoff,         # the MCID
                   decrease,       # whether the MCID change is decrease in the response variable (TRUE\FALSE)
                   dist = "none",  # distribution type (if specified). The default value is 'none'.
                   equal.var = TRUE ) {   # whether the variances are equal - applies only for normal dist. The default value is TRUE.

  if( !(type %in% c("laupacis", "mle", "fl") ) ) warning("type can be 'laupacis', 'mle' or 'fl' only")
  if( !is.numeric(treat) )                       warning("treatment arm vector (treat) must be a numeric vector")
  if( !is.numeric(control) )                     warning("control arm vector (treat) must be a numeric vector")
  if( !is.numeric(cutoff) | length(cutoff) > 1 ) warning("cutoff must be a scalar")
  if( !(decrease %in% c(T, F, TRUE, FALSE) ) )   warning("decrease must be TRUE or FALSE")
  if( type == "mle" & !(dist %in% c("normal", "expon") ) ) warning("distribution (dist) can be 'expon' (Exponential) or 'normal' only")
#  if( type == "mle" & dist == "normal" & !(equal.var %in% c(T, F, TRUE, FALSE) ) )  warning("equal.var must be TRUE or FALSE")
#  if( type == "laupacis" )                       message("if 'dist' and/or 'equal.var' arguments were set they are ignored")
  if( type == "fl" & dist != "normal" )          warning("for type = fl distribution (dist) must be normal")
  if( type == "fl" & equal.var != T )            warning("for fl type estimator equal.var must be TRUE")

  # remove NAs

  treat   <- treat[   !is.na( treat )  ]
  control <- control[ !is.na( control )]

      # initial values
  n_c     = length( control )
  n_t     = length( treat )

  yc_bar  = mean( control )
  yt_bar  = mean( treat )

  tau     = cutoff

  s_ml    = ( 1 / ( n_c + n_t)  * ( (n_c - 1) * var( control ) + (n_t - 1) * var( treat ) ) ) ^ ( 1/2 )

  d       = ( mean(treat) - mean(control) ) / s_ml

  s_t     = ( ( n_t - 1 ) * var( treat ) /  n_t ) ^ ( 1 / 2 )

  s_c     = ( ( n_c - 1 ) * var( control ) / n_c ) ^ ( 1 / 2 )

  ### BOOTSTRAP FUNCTIONS (INCREASE) ###

  p_c1   = function(data, indices) {
    p.c     = ifelse( mean( data[indices] > cutoff, na.rm = T ) < 0.001, # prob of success in the control grp
                      0.001,
                      mean( data[indices] > cutoff, na.rm = T ))
    var.pc  = var( data[indices] )                                       # variance of sampled control values
    mean.pc = mean( data[indices] )                                      # mean of sampled control values
    return(c(p.c, var.pc, mean.pc))
  }

  p_t1 = function(data, indices) {
    p.t     = ifelse( mean( data[indices] > cutoff, na.rm = T ) > 0.999, # prob of success in the treatment grp
                      0.999,
                      mean( data[indices] > cutoff, na.rm = T  ) )
    var.pt  = var( data[indices] )                                       # variance of sampled treatment values
    mean.pt = mean( data[indices] )                                      # mean of sampled treatment values
    return(c(p.t, var.pt, mean.pt))
  }

   ###########################
   ### LAUPACIS - ANY DIST ###
   ###########################

  ### INCREASE ###

  if( type == "laupacis" & decrease == F & ( dist == "none" || missing(dist) == TRUE ) ) {

     # probability treatment
   p_t     =  ifelse( mean( treat   > cutoff, na.rm = T )  > 0.999,
                      0.999,
                      mean( treat   > cutoff, na.rm = T ) )

     # probability control
   p_c     =  ifelse( mean( control > cutoff, na.rm = T ) < 0.001,
                     0.001,
                     mean( control > cutoff, na.rm = T) )


  # 95% quantile BS confidence interval
  p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
  p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

  ### BOOTSTRAP ESTIMATORS ###

  s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^ ( 1 / 2 )

  s_t.bs     = ( ( n_t - 1 )   * p_t.boot$t[, 2] / n_t ) ^ ( 1 / 2 )

  s_c.bs     = ( ( n_c - 1 )   * p_c.boot$t[, 2] / n_c ) ^ ( 1 / 2 )

  return( laupacis( treat, control, p_t, p_c, n_t, n_c, p_t.boot, p_c.boot ) )

    }

  ### DECREASE ###

   if( type == "laupacis" & decrease == T & ( dist == "none" || missing(dist) == TRUE ) ) {

     # prob treatment
     p_t     =  ifelse( mean( treat   < cutoff, na.rm = T ) > 0.999,
                       0.999,
                       mean( treat   < cutoff, na.rm = T )  )

     # prob control
     p_c     =  ifelse( mean( control < cutoff, na.rm = T ) < 0.001,
                       0.001,
                       mean( control < cutoff, na.rm = T)  )

    ### BOOTSTRAP FUNCTIONS ###

    p_c1   = function(data, indices) {
      p.c     = ifelse( mean(ifelse( data[indices] < cutoff, 1, 0)) < 0.001,
                        0.001,
                        mean(ifelse( data[indices] < cutoff, 1, 0)) )
      var.pc  = var( data[indices] )
      mean.pc = mean( data[indices] )
      return(c(p.c, var.pc, mean.pc))
    }

    p_t1 = function(data, indices) {
      p.t     = ifelse( mean( ifelse( data[indices] < cutoff, 1, 0) ) > 0.999,
                        0.999,
                        mean( ifelse( data[indices] < cutoff, 1, 0) ) )
      var.pt  = var( data[indices] )
      mean.pt = mean( data[indices] )
      return(c(p.t, var.pt, mean.pt))
    }

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    ### BOOTSTRAP ESTIMATORS ###

    s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^ ( 1 / 2 )

    s_t.bs     = ( ( n_t - 1 )   * p_t.boot$t[, 2] / n_t ) ^ ( 1 / 2 )

    s_c.bs     = ( ( n_c - 1 )   * p_c.boot$t[, 2] / n_c ) ^ ( 1 / 2 )

    return( laupacis( treat, control, p_t, p_c, n_t, n_c, p_t.boot, p_c.boot ) )
    }

  #################
  ### MLE EXPON ###
  #################

  ### INCREASE ###

  if( type == "mle" & decrease == F & dist == "expon" ) {

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    return( mle_exp_inc( treat, control, tau, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) )
  }

  ### DECREASE ###

  if( type == "mle" & decrease == T & dist == "expon" ) {

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    return( mle_exp_dec( treat, control, tau, yt_bar, yc_bar, p_t.boot, p_c.boot, n_c, n_t ) )
  }

  ##############################
  ### NORM EQ VAR - FURUKAWA ###
  ##############################

  ### INCREASE ###

  if( type == "fl" & decrease == F & dist == "normal" ) {

    # probability control
    p_c     =  ifelse( mean( control > cutoff, na.rm = T ) < 0.001,
                       0.001,
                       mean( control > cutoff, na.rm = T) )

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    ### BOOTSTRAP ESTIMATORS ###

    s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^ ( 1 / 2 )

    d.bs       = ( p_t.boot$t[ ,3] - p_c.boot$t[ ,3] ) / s_ml.bs

    return( furukawa_inc( treat, control, tau, d, d.bs, p_c, p_t.boot,  p_c.boot, n_c, n_t ) )
  }

  ### DECREASE ###

  if( type == "fl" & decrease == T & dist == "normal" ) {

    p_c     =  ifelse( mean( control < cutoff, na.rm = T ) < 0.001,
                       0.001,
                       mean( control < cutoff, na.rm = T)  )

    ### BOOTSTRAP FUNCTIONS ###

    p_c1   = function(data, indices) {
      p.c     = ifelse( mean(ifelse( data[indices] < cutoff, 1, 0)) < 0.001,
                        0.001,
                        mean(ifelse( data[indices] < cutoff, 1, 0)) )
      var.pc  = var( data[indices] )
      mean.pc = mean( data[indices] )
      return(c(p.c, var.pc, mean.pc))
    }

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^ ( 1 / 2 )

    d.bs       = ( p_t.boot$t[ ,3] - p_c.boot$t[ ,3] ) / s_ml.bs

    return( furukawa_dec( treat, control, tau, d, d.bs, p_c, p_t.boot,  p_c.boot, n_c, n_t ) )
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

    s_t.bs     = ( ( n_t - 1 )   * p_t.boot$t[, 2] / n_t ) ^ ( 1 / 2 )

    s_c.bs     = ( ( n_c - 1 )   * p_c.boot$t[, 2] / n_c ) ^ ( 1 / 2 )

     return( mle_norm_eq_inc( treat, control, tau, yt_bar, yc_bar, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_t, n_c ) )
    }

  ### DECREASE ###

   if( type == "mle" & decrease == T & dist == "normal" & equal.var == T ) {

     ### BOOTSTRAP FUNCTIONS ###

    # 95% quantile BS confidence interval
    p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
    p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

    ### BOOTSTRAP ESTIMATORS ###

    s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^ ( 1 / 2 )

    s_t.bs     = ( ( n_t - 1 )   * p_t.boot$t[, 2] / n_t ) ^ ( 1 / 2 )

    s_c.bs     = ( ( n_c - 1 )   * p_c.boot$t[, 2] / n_c ) ^ ( 1 / 2 )

    return( mle_norm_eq_dec( treat, control, tau, yt_bar, yc_bar, s_ml, p_t.boot, p_c.boot, s_ml.bs, n_t, n_c ) )
  }

    #########################
    ### MLE NORM UNEQ VAR ###
    #########################

    ### INCREASE ###

    if( type == "mle" & decrease == F & dist == "normal" & equal.var == F ) {

      # 95% quantile BS confidence interval
      p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
      p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

      ### BOOTSTRAP ESTIMATORS ###

      s_t.bs     = ( ( n_t - 1 )   * p_t.boot$t[, 2] / n_t ) ^ ( 1 / 2 )

      s_c.bs     = ( ( n_c - 1 )   * p_c.boot$t[, 2] / n_c ) ^ ( 1 / 2 )

      return( mle_norm_uneq_inc( treat, control, tau, yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, s_t.bs, s_c.bs, n_t, n_c ) )
    }

   ### DECREASE ###

    if( type == "mle" & decrease == T & dist == "normal" & equal.var == F ) {

      ### BOOTSTRAP FUNCTIONS ###

      # 95% quantile BS confidence interval
      p_c.boot = boot(data = control, statistic = p_c1, R = 1000)
      p_t.boot = boot(data = treat,   statistic = p_t1, R = 1000)

      s_ml.bs    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot$t[ ,2] + (n_t - 1) * p_t.boot$t[ ,2] ) ) ^( 1 / 2 )

      s_t.bs     = ( ( n_t - 1 )   * p_t.boot$t[, 2] / n_t ) ^ ( 1 / 2 )

      s_c.bs     = ( ( n_c - 1 )   * p_c.boot$t[, 2] / n_c ) ^ ( 1 / 2 )

      return( mle_norm_uneq_dec(  treat, control, tau, yt_bar, yc_bar, s_t, s_c, p_t.boot, p_c.boot, s_t.bs, s_c.bs, n_t, n_c ) )
    }
 }

