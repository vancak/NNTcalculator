#' @title Laupacis' adjusted and marginal NNT in logistic regression model
#'
#' @description Internal function. Not for users. Calculates Laupacis' adjusted and marginal NNT
#' in logistic regression model.
#' Takes a data-set suitable for a logistic regression with a dummy variable and an
#' interaction term, and returns the estimated adjusted NNT(x) and the marginal NNT.
#' @param response      vector of the response variable (i.e., the dependent variable).
#' @param x             vector of the explanatory variable.
#' @param group         allocated arm variable where 1 corresponds to the treatment arm, and 0 to the control arm.
#' @param adj           value that the NNT need to be adjusted for. The default value is the mean of x.
#' @param data          data frame that contains the required variables for the computations.
#'
#' @return The estimated marginal and adjusted NNT with their corresponding 95 percent confidence intervals.
#' @keywords internal
#'
#' load("logreg_data.RData" )
#'
#' nnt_logreg( response = logreg_data$y,
#'             x        = logreg_data$x_var,
#'             group    = logreg_data$gr,
#'             adj      = 2,
#'             data     = logreg_data )
#'
nnt_logreg <- function( response,       # vector of the response variable
                        x,              # vector of the explanatory variable
                        group,          # vector of 0/1 indicators for control and treatment arms, respectively
                        adj,            # specific value that the NNT need to be adjusted for. The default value is mean of x
                        data ) {        # the analyzed data-set
# library(boot)

  ### initial values ###
  dat1    = data

  dat1    = data.frame( y = response, x = x, gr = group )

  attach(dat1, warn.conflicts = F)

  gr      = dat1$gr

  adj     = ifelse( !missing(adj), adj, round( mean(dat1$x, na.rm = T), 2) )

  fun_g   = function(h){ ifelse( h > 0, 1/h, Inf ) }

  sum = summary( glm( y ~ gr + x + x * gr, data = dat1, family = "binomial" ) )
  sum

    #####################
    ### ADJUSTED NNTs ###
    #####################

    ### ESTIMATION ###
  p_t_mle = function(h){ 1 / ( 1 + exp( - ( coef(sum)[1, 1] + coef(sum)["gr", 1] + (coef(sum)["gr:x", 1] + coef(sum)["x", 1]) * h ) ) ) }

  p_c_mle = function(h){ 1 / ( 1 + exp( - ( coef(sum)[1, 1] + coef(sum)["x", 1] * h ) ) ) }

  ps_x    = function(h){  p_t_mle(h) - p_c_mle(h)  }

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
    deriv_psb0x      = function(h){   p_c_mle(h) * ( 1 - p_c_mle(h) ) - p_t_mle(h) * ( 1 - p_t_mle(h) ) }
    deriv_psbgrx     = function(h){ - p_t_mle(h) * ( 1 - p_t_mle(h) ) }
    deriv_psbvarx    = function(h){ ( p_c_mle(h) * ( 1 - p_c_mle(h) ) - p_t_mle(h) * ( 1 - p_t_mle(h) ) ) * h }
    deriv_psbgrvarx  = function(h){ - p_t_mle(h) * ( 1 - p_t_mle(h) * h ) }

    grad_psx        = function(h) { c( deriv_psb0x(h), deriv_psbgrx(h), deriv_psbvarx(h), deriv_psbgrvarx(h) ) }

    cov_b           =  vcov( sum )

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

      fit      = glm( y ~ gr + x + x * gr,
                      data   = d,
                      family = "binomial" )

      sumb     = summary.lm(fit)

      pbt_mle  = function(h){ 1 / ( 1 + exp( - ( coef(sumb)[1, 1] + coef(sumb)["gr", 1] + (coef(sumb)["gr:x", 1] + coef(sumb)["x", 1]) * h ) ) ) }

      pbc_mle  = function(h){ 1 / ( 1 + exp( - ( coef(sumb)[1, 1] + coef(sumb)["x", 1] * h ) ) ) }

      ps_bx     = function(h){  pbt_mle(h) - pbc_mle(h)  }

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
                  mean  = c( sum$coef[1:4] ),
                  sigma =  cov_b )

    pbt_mle  = function(h){ 1 / ( 1 + exp( - ( pbs[,1] + pbs[,2] + (pbs[,3] + pbs[,4]) * h ) ) ) }

    pbc_mle  = function(h){ 1 / ( 1 + exp( - ( pbs[,1] + pbs[,3] * h ) ) ) }

    pbs_bx   = function(h){  pbt_mle(h) - pbc_mle(h)  }

    NNT_PBS  = fun_g( pbs_bx( adj ) )

    cipbs_l  = max( NNT_X - 1.96 * sd(NNT_PBS, na.rm = T), 1 )
    cipbs_r  = NNT_X      + 1.96 * sd(NNT_PBS, na.rm = T)

    av_ppbs  = numeric()

    for(k in 1:1000){
      av_ppbs[k] = mean( 1 / ( 1 + exp( - ( pbs[k,1] + pbs[k,2] + (pbs[k,3] + pbs[k,4]) * x ) ) ) -
                         1 / ( 1 + exp( - ( pbs[k,1] + pbs[k,3] * x ) ) ) )
    }

    NNT_UNPBS = fun_g( av_ppbs )


    ########################
    ### MARGINAL NNT CIs ###
    ########################

    grad_unps   =  c( mean( deriv_psb0x(x)     ),
                      mean( deriv_psbgrx(x)    ),
                      mean( deriv_psbvarx(x)   ),
                      mean( deriv_psbgrvarx(x) ) )

    grad_ml     = - NNT_MLE ^ 2  * grad_unps

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
                   cutoff    = 0.5,
                   dist      = "none",
                   equal.var = T,
                   decrease  = FALSE )


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
                            "CI TR L",  "CI TR U",
                            "CI DL L",  "CI DL U",
                            "CI NBS L", "CI NBS U",
                            "CI PBS L", "CI PBS U"  )

    rownames(out_list1) = c("NNT L", "NNT MLE", paste(c("NNT", "(", adj, ")" ), sep = "", collapse = "") )

    return( out_list1 )

  }

