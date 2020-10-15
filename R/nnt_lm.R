#' @title Laupacis' adjusted and marginal NNT in linear regression model
#'
#' @description Internal function. Not for users. Calculates Laupacis' adjusted and marginal NNT
#' in linear regression model.
#' Takes a data-set suitable for a simple linear regression with a dummy variable and an
#' interaction term, and returns the estimated adjusted NNT(x) and the marginal NNT.
#' @param response      vector of the response variable (i.e., the dependent variable).
#' @param x             vector of the explanatory variable.
#' @param cutoff        the MCID threshold.
#' @param decrease      TRUE or FALSE. Indicates whether the MCID change is decrease in the response variable
#' @param group         allocated arm variable where 1 corresponds to the treatment arm, and 0 to the control arm.
#' @param adj           value that the NNT need to be adjusted for. The default value is the mean of x.
#' @param data          data frame that contains the required variables for the computations.
#' @keywords internal
#' @return The estimated marginal and adjusted NNT with their corresponding 95 percent confidence intervals.
#'
#' data(linreg_data)
#'
# nnt_lm( response = linreg_data$y,
#         x        = linreg_data$x_var,
#         cutoff   = 3,
#         group    = linreg_data$gr,
#         decrease = F,
#         adj      = 2.6,
#         data     = linreg_data )
#
# inv_dat = data.frame( y = linreg_data$y, x_var = linreg_data$x_var, gr = 1 - linreg_data$gr )
#
# nnt_lm( response = inv_dat$y,
#         x        = inv_dat$x_var,
#         cutoff   = 3,
#         group    = inv_dat$gr,
#         decrease = T,
#         adj      = 2.6,
#         data     = inv_dat )
nnt_lm <- function( response,       # vector of the response variable
                    x,              # vector of the explanatory variable
                    cutoff,         # the MCID
                    group,          # vector of 0/1 indicators for control and treatment arms, respectively
                    decrease,       # whether the MCID change is decrease in the response variable (TRUE\FALSE)
                    adj,            # specific value that the NNT need to be adjusted for. The default value is mean of x
                    data ) {        # the analyzed data-set
# library(boot)

  ### initial values ###
  dat1    = data

  dat1    = data.frame( y = response, x = x, gr = group )

  attach(dat1, warn.conflicts = F)

  gr      = dat1$gr

  fun_g   = function(h){ ifelse( h > 0, 1/h, Inf ) }

  tau     = cutoff

  adj     = ifelse( !missing(adj), adj, round( mean(dat1$x, na.rm = T), 2) )

  ### control arm ###
  sum = summary( lm( y ~  gr + x + x * gr, data = dat1 ) )
  sum

  sig_m   = sum$sigma


  ### INCREASE MARGINAL ###

  if( decrease == F ) {

    #####################
    ### ADJUSTED NNTs ###
    #####################

    ### ESTIMATION ###
    ps_x  = function( h ){  (1 - pnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m ) ) -
                            (1 - pnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m ) )   }

    NNT_X = fun_g( ps_x(adj) )

    #####################
    ### MARGINAL NNTs ###
    #####################

    av_ps     = mean( ps_x(x) )

    NNT_MLE   = fun_g( av_ps )

    #########################
    ### ADJUSTED NNT CIs  ###
    #########################

    ### PARTIAL DERIVATIVES ###

    ### ps(x) deriv ###
    deriv_psb0x     = function(h){ dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( 1 / sig_m ) -
                                   dnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m )   * ( 1 / sig_m ) }

    deriv_psbgrx    = function(h){ dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( 1 / sig_m ) }

    deriv_psbvarx   = function(h){ dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( h / sig_m ) -
                                   dnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m )   * ( h / sig_m ) }

    deriv_psbgrvarx = function(h){ dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( h / sig_m ) }

    deriv_pssigx    = function(h){ dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( - ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m ^ 2 ) -
                                   dnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m ) * ( - ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m ^ 2 ) }

    grad_psx        = function(h) { c( deriv_psb0x(h), deriv_psbgrx(h), deriv_psbvarx(h), deriv_psbgrvarx(h), deriv_pssigx(h) ) }

    cov_b           = rbind( cbind( vcov( sum ), rep(0, 4) ),
                             c( rep(0, 4), 2 * sig_m ^ 4 / ( nrow(dat1) - 4 ) ) )

    citr_l  = numeric()
    citr_r  = numeric()

    cidl_l  = numeric()
    cidl_r  = numeric()

    cinbs_l = numeric()
    cinbs_r = numeric()

    cipbs_l = numeric()
    cipbs_r = numeric()

    ### TR & DL CIs ###
    citr_l   =  max( fun_g( ps_x(adj)  + 1.96  * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) ) ), 1)
    citr_r   =  fun_g( ps_x(adj)       - 1.96  * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) ) )

    cidl_l   =  max( NNT_X - 1.96  *  NNT_X ^ 2 * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) ), 1 )
    cidl_r   =  NNT_X      + 1.96  *  NNT_X ^ 2 * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) )

    ### NBS CIs ###

    # function to obtain NNT(K) from the data
    nnt_1  = function(data, indices) {

      d        = data[indices,]

      fit      = aov(y ~ gr + x + x * gr , data = d)

      sumb     = summary.lm(fit)

      ps_bx    = function( h ) { (1 - pnorm( ( tau - sumb$coef[1] - sumb$coef[2] - ( sumb$coef[3] + sumb$coef[4] ) * h ) / sumb$sigma ) )  -
                                (1 - pnorm( ( tau - sumb$coef[1] - sumb$coef[3] * h ) / sumb$sigma ) ) }

      NNT_XB   = fun_g( ps_bx( adj ) )

      av_psb   = mean( ps_bx( x ) )

      NNT_UNB  = fun_g( av_psb )

      return(c( NNT_XB, NNT_UNB ))

        }

    # bootstrapping with 1000 replications
    results <- boot(data      = dat1,
                    statistic = nnt_1,
                    R         = 1000)


    # 95% NBS confidence intervals

      cinbs_l = quantile(results$t[,1], probs = c(0.025))
      cinbs_r = quantile(results$t[,1], probs = c(0.975))

    ### Parametric Bootstrap 95% CI ###
 #   library(mvtnorm)

    pbs = rmvnorm(1000,
                  mean  = c( sum$coef[1:4], sig_m ^ 2 ),
                  sigma =  cov_b )

    pbs_bx   = function(h) { (1 - pnorm( ( tau - pbs[,1] - pbs[,2] - ( pbs[,3] + pbs[,4] ) * h ) / sqrt(pbs[,5]) ) )  -
                             (1 - pnorm( ( tau - pbs[,1] - pbs[,3] * h ) / sqrt(pbs[,5]) ) ) }

    NNT_PBS  = fun_g( pbs_bx( adj ) )

    cipbs_l  = max( NNT_X - 1.96 * sd(NNT_PBS, na.rm = T), 1 )
    cipbs_r  = NNT_X      + 1.96 * sd(NNT_PBS, na.rm = T)

    av_ppbs = numeric()

    for(k in 1:1000){
      av_ppbs[k] = mean(  (1 - pnorm( ( tau - pbs[k,1] - pbs[k,2] - ( pbs[k,3] + pbs[k,4] ) * x ) / sqrt(pbs[k,5]) ) )  -
                          (1 - pnorm( ( tau - pbs[k,1] - pbs[k,3] * x ) / sqrt(pbs[k,5]) ) ) )
    }

    NNT_UNPBS = fun_g( av_ppbs )


    ########################
    ### MARGINAL NNT CIs ###
    ########################

    grad_unps  =  c( mean( deriv_psb0x(x)     ),
                     mean( deriv_psbgrx(x)    ),
                     mean( deriv_psbvarx(x)   ),
                     mean( deriv_psbgrvarx(x) ),
                     mean( deriv_pssigx(x)    )  )

    grad_ml = - NNT_MLE ^ 2  * grad_unps

    citr_mlel   =  max( fun_g( av_ps + 1.96  * sqrt( grad_unps %*% cov_b %*% grad_unps ) ), 1 )
    citr_mler   =  fun_g( av_ps - 1.96  * sqrt(  grad_unps %*% cov_b %*% grad_unps ) )

    cidl_mlel   =  max( NNT_MLE - 1.96  * sqrt( grad_ml %*% cov_b %*% grad_ml ), 1 )
    cidl_mler   =  NNT_MLE      + 1.96  * sqrt( grad_ml %*% cov_b %*% grad_ml )

    cinbs_mlel  = quantile(results$t[,2],  probs = c(0.025))
    cinbs_mler  = quantile(results$t[,2], probs = c(0.975))

    cipbs_mlel   =  max( NNT_MLE - 1.96  * sd( NNT_UNPBS, na.rm = T ) , 1 )
    cipbs_mler   =  NNT_MLE      + 1.96  * sd( NNT_UNPBS, na.rm = T )

    ### MARGINAL NNT L CIs ###

    nntl =  nnt_l( type      = "laupacis",
                   treat     = dat1[ gr == 1, 1 ],
                   control   = dat1[ gr == 0, 1 ],
                   cutoff    = tau,
                   dist      = "none",
                   equal.var = T,
                   decrease  = decrease )


    ### FINAL CIs DATA FRAME ###

    out_list = list()

    out_list[[1]] = c( nntl, c(NA, NA) )

    out_list[[2]] = t(c(NNT_MLE,
                        citr_mlel,  citr_mler,
                        cidl_mlel,  cidl_mler,
                        cinbs_mlel, cinbs_mler,
                        cipbs_mlel, cipbs_mler))


      out_list[[3]] =  t(c(NNT_X,
                             citr_l,  citr_r,
                             cidl_l,  cidl_r,
                             cinbs_l, cinbs_r,
                             cipbs_l, cipbs_r))


    out_list1 = matrix( unlist( out_list ), ncol = 9, byrow = T )

        colnames(out_list1) = c("NNT",
                                "CI TR L",  "CI TR R",
                                "CI DL L",  "CI DL R",
                                "CI NBS L", "CI NBS R",
                                "CI PBS L", "CI PBS R"  )

    rownames(out_list1) = c("NNT L", "NNT MLE", paste(c("NNT", "(", adj, ")" ), sep = "", collapse = "") )

    return( out_list1 )

  }

############################################################################################################
##########################
### SUCCESS = DECREASE ###
##########################

  if( decrease == T ) {

    #####################
    ### ADJUSTED NNTs ###
    #####################

    ### ESTIMATION ###
    ps_x  = function( h ){  - (1 - pnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m ) ) +
        (1 - pnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m ) )   }

    NNT_X = fun_g( ps_x(adj) )

    #####################
    ### MARGINAL NNTs ###
    #####################

    av_ps     = mean( ps_x(x) )

    NNT_MLE   = fun_g( av_ps )

    #########################
    ### ADJUSTED NNT CIs  ###
    #########################

    ### PARTIAL DERIVATIVES ###

    ### ps(x) deriv ###
    deriv_psb0x     = function(h){ - dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( 1 / sig_m ) +
        dnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m )   * ( 1 / sig_m ) }

    deriv_psbgrx    = function(h){ dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( 1 / sig_m ) }

    deriv_psbvarx   = function(h){ - dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( h / sig_m ) +
        dnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m )   * ( h / sig_m ) }

    deriv_psbgrvarx = function(h){ dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( h / sig_m ) }

    deriv_pssigx    = function(h){ - dnorm( ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m )  * ( - ( tau - sum$coef[1] - sum$coef[2] - ( sum$coef[3] + sum$coef[4] ) * h ) / sig_m ^ 2 ) +
        dnorm( ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m ) * ( - ( tau - sum$coef[1] - sum$coef[3] * h ) / sig_m ^ 2 ) }

    grad_psx        = function(h) { c( deriv_psb0x(h), deriv_psbgrx(h), deriv_psbvarx(h), deriv_psbgrvarx(h), deriv_pssigx(h) ) }

    cov_b           = rbind( cbind( vcov( sum ), rep(0, 4) ),
                             c( rep(0, 4), 2 * sig_m ^ 4 / ( nrow(dat1) - 4 ) ) )

    citr_l  = numeric()
    citr_r  = numeric()

    cidl_l  = numeric()
    cidl_r  = numeric()

    cinbs_l = numeric()
    cinbs_r = numeric()

    cipbs_l = numeric()
    cipbs_r = numeric()

    ### TR & DL CIs ###
    citr_l   =  max( fun_g( ps_x(adj)  + 1.96  * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) ) ), 1)
    citr_r   =  fun_g( ps_x(adj)       - 1.96  * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) ) )

    cidl_l   =  max( NNT_X - 1.96  *  NNT_X ^ 2 * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) ), 1 )
    cidl_r   =  NNT_X      + 1.96  *  NNT_X ^ 2 * sqrt( grad_psx(adj) %*% cov_b %*% grad_psx(adj) )

    ### NBS CIs ###

    # function to obtain NNT(K) from the data
    nnt_1  = function(data, indices) {

      d        = data[indices,]

      fit      = aov(y ~ gr + x + x * gr , data = d)

      sumb     = summary.lm(fit)

      ps_bx    = function( h ) { - (1 - pnorm( ( tau - sumb$coef[1] - sumb$coef[2] - ( sumb$coef[3] + sumb$coef[4] ) * h ) / sumb$sigma ) )  +
          (1 - pnorm( ( tau - sumb$coef[1] - sumb$coef[3] * h ) / sumb$sigma ) ) }

      NNT_XB   = fun_g( ps_bx( adj ) )

      av_psb   = mean( ps_bx( x ) )

      NNT_UNB  = fun_g( av_psb )

      return(c( NNT_XB, NNT_UNB ))

    }

    # bootstrapping with 1000 replications
    results <- boot(data      = dat1,
                    statistic = nnt_1,
                    R         = 1000)


    # 95% NBS confidence intervals

    cinbs_l = quantile(results$t[,1], probs = c(0.025))
    cinbs_r = quantile(results$t[,1], probs = c(0.975))

    ### Parametric Bootstrap 95% CI ###
  #  library(mvtnorm)

    pbs = rmvnorm(1000,
                  mean  = c( sum$coef[1:4], sig_m ^ 2 ),
                  sigma =  cov_b )

    pbs_bx   = function(h) { - (1 - pnorm( ( tau - pbs[,1] - pbs[,2] - ( pbs[,3] + pbs[,4] ) * h ) / sqrt(pbs[,5]) ) )  +
        (1 - pnorm( ( tau - pbs[,1] - pbs[,3] * h ) / sqrt(pbs[,5]) ) ) }

    NNT_PBS  = fun_g( pbs_bx( adj ) )

    cipbs_l  = max( NNT_X - 1.96 * sd(NNT_PBS, na.rm = T), 1 )
    cipbs_r  = NNT_X      + 1.96 * sd(NNT_PBS, na.rm = T)

    av_ppbs = numeric()

    for(k in 1:1000){
      av_ppbs[k] = mean(  - (1 - pnorm( ( tau - pbs[k,1] - pbs[k,2] - ( pbs[k,3] + pbs[k,4] ) * x ) / sqrt(pbs[k,5]) ) )  +
                            (1 - pnorm( ( tau - pbs[k,1] - pbs[k,3] * x ) / sqrt(pbs[k,5]) ) ) )
    }

    NNT_UNPBS = fun_g( av_ppbs )


    ########################
    ### MARGINAL NNT CIs ###
    ########################

    grad_unps  =  c( mean( deriv_psb0x(x)     ),
                     mean( deriv_psbgrx(x)    ),
                     mean( deriv_psbvarx(x)   ),
                     mean( deriv_psbgrvarx(x) ),
                     mean( deriv_pssigx(x)    )  )

    grad_ml = - NNT_MLE ^ 2  * grad_unps

    citr_mlel   =  max( fun_g( av_ps + 1.96  * sqrt( grad_unps %*% cov_b %*% grad_unps ) ), 1 )
    citr_mler   =  fun_g( av_ps - 1.96  * sqrt(  grad_unps %*% cov_b %*% grad_unps ) )

    cidl_mlel   =  max( NNT_MLE - 1.96  * sqrt( grad_ml %*% cov_b %*% grad_ml ), 1 )
    cidl_mler   =  NNT_MLE      + 1.96  * sqrt( grad_ml %*% cov_b %*% grad_ml )

    cinbs_mlel  = quantile(results$t[,2],  probs = c(0.025))
    cinbs_mler  = quantile(results$t[,2], probs = c(0.975))

    cipbs_mlel   =  max( NNT_MLE - 1.96  * sd( NNT_UNPBS, na.rm = T ) , 1 )
    cipbs_mler   =  NNT_MLE      + 1.96  * sd( NNT_UNPBS, na.rm = T )

    ### MARGINAL NNT L CIs ###

    nntl =  nnt_l( type      = "laupacis",
                   treat     = dat1[ gr == 1, 1 ],
                   control   = dat1[ gr == 0, 1 ],
                   cutoff    = tau,
                   dist      = "none",
                   equal.var = T,
                   decrease  = decrease )


    ### FINAL CIs DATA FRAME ###

    out_list = list()

    out_list[[1]] = c( nntl, c(NA, NA) )

    out_list[[2]] = t(c(NNT_MLE,
                        citr_mlel,  citr_mler,
                        cidl_mlel,  cidl_mler,
                        cinbs_mlel, cinbs_mler,
                        cipbs_mlel, cipbs_mler))


    out_list[[3]] =  t(c(NNT_X,
                         citr_l,  citr_r,
                         cidl_l,  cidl_r,
                         cinbs_l, cinbs_r,
                         cipbs_l, cipbs_r))


    out_list1 = matrix( unlist( out_list ), ncol = 9, byrow = T )

    colnames(out_list1) = c("NNT",
                            "CI TR L",  "CI TR R",
                            "CI DL L",  "CI DL R",
                            "CI NBS L", "CI NBS R",
                            "CI PBS L", "CI PBS R"  )

    rownames(out_list1) = c("NNT L", "NNT MLE", paste(c("NNT", "(", adj, ")" ), sep = "", collapse = "") )

    return( out_list1 )

  }}

