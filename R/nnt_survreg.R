#' @title NNT-KM and NNT-COX calculator
#'
#' @description Calculates Laupacis' type adjusted and marginal NNT(y|x) for survival data.
#' Takes a data-set suitable for a survival analysis, a fixed time point y,
#' and an explanatory variables x,  and returns
#' the estimated adjusted NNT(y|x), and the estimated unadjusted and marginal NNT(y).
#' In the Cox-model-based estimators, the baseline hazard is estimated by the Breslow's nonparametric MLE.
#' In the unadjusted estimator of NNT(y), the survival probabilities are estimated by the Kaplan-Meier nonparametric MLE.
#' @param response      vector of the response variable; times of the events/censoring
#' @param status        column that contains 0/1 indicator, where 1 is failure time, and 0 is censoring time.
#' @param x             vector of the explanatory variable.
#' @param group         allocated arm variable where 1 corresponds to the treatment arm, and 0 to the control arm.
#' @param adj           the x value that the NNT(y|x) need to be adjusted for. The default value is the mean of x.
#' @param time.point    the fixed time point y for NNT(y|x), and NNT(y)
#' @param data          analyzed data frame that contains the required variables for the computations.
#'
#' @return The estimated unadjusted, marginal and adjusted NNTs with their corresponding 95 percent confidence intervals.
#' @references Therneau, T. M., & Lumley, T. (2014). Package survival. Survival analysis Published on CRAN, 2, 3.
#' @references Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. Journal of the American statistical association, 53(282), 457-481.
#' @references Cox, D. R., & Oakes, D. (1984). Analysis of survival data (Vol. 21). CRC Press.
#' @references Cox, D. R. (1972). Regression models and life-tables. Journal of the Royal Statistical Society: Series B (Methodological), 34(2), 187-202.
#' @examples
#' data(survreg_data)
#'
#' nnt_survreg( response = survreg_data$stop,
#'              status   = survreg_data$status,
#'                    x  = survreg_data$x.1,
#'                group  = survreg_data$x,
#'                   adj = -1,
#'            time.point = 0.5,
#'                  data = survreg_data )
#' @export
nnt_survreg <- function( response,       # vector of the response variable
                         status,         # column that contains 0/1 indicator, where 1 is failure time, and 0 is censoring time.
                         x,              # vector of the explanatory variable
                         group,          # vector of 0/1 indicators for control and treatment arms, respectively
                         adj,            # specific value that the NNT need to be adjusted for. # The default value is the mean of x
                         time.point,     # the fixed time point y for adjusted and marginal NNTs
                         data ) {        # the analyzed data-set

  # require(survival)
  # require(boot)
  # require(MASS)

  ### initial values ###
  dat1    = data

  dat1    = data.frame( time = response, status = status, x = x, gr = group )

  t       = time.point

  adj     = ifelse( !missing(adj), adj, round( mean(dat1$x, na.rm = T), 2) )

  attach(dat1, warn.conflicts = F)

  fun_g   = function(h){ ifelse( h > 0, 1/h, Inf ) }

  coxph_both = coxph( Surv(time, status) ~  gr + x,
                      method = "breslow",
                      data   = dat1 )

  # baseline hazard

  surv0 = round( basehaz( coxph_both, centered = F ), 10 )

  # times of the baseline hazard

  dat_both = as.data.frame( cbind( surv0$time, exp(- surv0$hazard) ) )

  names( dat_both ) = c( "time", "surv_base" )

  diff_abs = abs( dat_both$time -  t  )

  #####################
  ### MARGINAL NNTs ###
  #####################

  ### NONPARAMETRIC NNT-KM ###

  ps_km    = function(t){

    diff_abst = abs( t - km_t$time )

    diff_absc = abs( t - km_c$time )

    p_t_km = km_t$surv[ which(diff_abst == min( diff_abst ) ) ]

    p_c_km = km_c$surv[ which(diff_absc == min( diff_absc ) ) ]

    return( p_t_km - p_c_km )
  }


  km_t = survfit( Surv(dat1[dat1$gr == 1, "time"], dat1[dat1$gr == 1, "status"] ) ~ 1)

  km_c = survfit( Surv(dat1[dat1$gr == 0, "time"], dat1[dat1$gr == 0, "status"] ) ~ 1 )

  NNT_KM = fun_g( ps_km(t) )

  ### PARAMETRIC NNT-COX ###

  tab <- data.frame(table(dat1[dat1$status == 1, "time"]))
  y   <- as.numeric(levels(tab[, 1]))[tab[, 1]]     #ordered distinct event times
  d   <- tab[, 2]

  fit = coxph_both

  cumhaz0_t =  function(t){

    h0 <- rep(NA, length(y))

    for(l in 1:length(y))
    {
      h0[l] <- d[l] / sum(exp(dat1[dat1$time >= y[l], "gr"] * coef(fit)[1] +
                              dat1[dat1$time >= y[l], "x"]  * coef(fit)[2]  ) )
    }

    diff_abs1 = function(t){  abs(t - y)  }

    return( cumsum(h0)[ which( diff_abs1(t) == min( diff_abs1(t) ) )] )
  }

ps_x = function( t, x ) {

      survt = exp( - cumhaz0_t( t ) ) ^ ( exp(  coef(fit)[1] + coef(fit)[2] * x ) )
      survc = exp( - cumhaz0_t( t ) ) ^ ( exp(  coef(fit)[2]  * x ) )

    return( survt - survc )
  }

  NNT_X   = fun_g( ps_x( t, adj ) )

  av_psb  = mean( ps_x(t, dat1$x) )

  NNT_COX = fun_g( av_psb )

  ############################
  ### CONFIDENCE INTERVALS ###
  ############################

  citr_kml    = NA
  citr_kmr    = NA

  cidl_kml    = NA
  cidl_kmr    = NA

  cinbs_kml   = numeric()
  cinbs_kmr   = numeric()

  cipbs_kml   = NA
  cipbs_kmr   = NA

  citr_coxl   =  NA
  citr_coxr   =  NA

  cidl_coxl   =  NA
  cidl_coxr   =  NA

  cinbs_coxl  =  numeric()
  cinbs_coxr  =  numeric()

  cipbs_coxl  =  NA
  cipbs_coxr  =  NA

  citr_l  = NA
  citr_r  = NA

  cidl_l  = NA
  cidl_r  = NA

  cinbs_l = numeric()
  cinbs_r = numeric()

  cipbs_l = NA
  cipbs_r = NA

  ### PARTIAL DERIVATIVES ###

  cumhaz0_tderiv = function(t, x){

    diff_abs1 = function(t){ abs(t - y) }

    cumhazderiv1 <- rep(NA, length(y))
    cumhazderiv2 <- rep(NA, length(y))

    for(l in 1:length(y))
    {
      cumhazderiv1[l] <- ( d[l] / sum( exp(dat1[ dat1$time >= y[l], "gr"] * coef(fit)[1] +
                                           dat1[ dat1$time >= y[l], "x"]  * coef(fit)[2]  ) ) ) ^ 2 *
                                  sum( exp(dat1[ dat1$time >= y[l], "gr"] * coef(fit)[1] +
                                           dat1[ dat1$time >= y[l], "x"]  * coef(fit)[2]  ) * dat1[ dat1$time >= y[l], "gr"] )

      cumhazderiv2[l] <- ( d[l] / sum( exp(dat1[ dat1$time >= y[l], "gr"] * coef(fit)[1] +
                                           dat1[ dat1$time >= y[l], "x"]  * coef(fit)[2]  ) ) ) ^ 2 *
                                  sum( exp(dat1[ dat1$time >= y[l], "gr"] * coef(fit)[1] +
                                           dat1[ dat1$time >= y[l], "x"]  * coef(fit)[2]  ) * dat1[ dat1$time >= y[l], "x" ] )
    }

    return( cbind( cumhazderiv1, cumhazderiv2 ) )
  }


  b1_deriv = function( t, x ){

    h0 <- rep(NA, length(y))

    for(l in 1:length(y))
     {
       h0[l] <- d[l] / sum(exp(dat1[ dat1$time >= y[l], "gr"] * coef(fit)[1] +
                               dat1[ dat1$time >= y[l], "x"]  * coef(fit)[2]  ) )
     }

    diff_abs1 = function(t){  abs(t - y)  }

    survt = exp( - cumhaz0_t( t ) ) ^ ( exp( coef(fit)[1] + coef(fit)[2] * x ) )
    survc = exp( - cumhaz0_t( t ) ) ^ ( exp( coef(fit)[2] * x ) )

    (- 1) * cumsum(h0)[which( diff_abs1(t) == min( diff_abs1(t) ) )] *
      ( survt * exp( coef(fit)[1] + coef(fit)[2] * x ) - survc * exp( coef(fit)[2] * x ) ) * x +
      ( survt * exp( coef(fit)[1] + coef(fit)[2] * x ) - survc * exp( coef(fit)[2] * x ) ) *
      cumsum(cumhaz0_tderiv( t, x )[ ,"cumhazderiv1"])[ which( diff_abs1(t) == min( diff_abs1(t) ) ) ]
  }


  b2_deriv = function( t, x ){

    h0 <- rep(NA, length(y))

    for(l in 1:length(y))
    {
      h0[l] <- d[l] / sum(exp(dat1[ dat1$time >= y[l], "gr"] * coef(fit)[1] +
                              dat1[ dat1$time >= y[l], "x"]  * coef(fit)[2]  ) )
    }

    diff_abs1 = function(t){  abs(t - y)  }

    survt = exp( - cumhaz0_t( t ) ) ^ ( exp(  coef(fit)[1] + coef(fit)[2] * x ) )
    survc = exp( - cumhaz0_t( t ) ) ^ ( exp(  coef(fit)[2]  * x ) )

    (- 1) * cumsum(h0)[which( diff_abs1(t) == min( diff_abs1(t) ) )] *
      ( survt * exp( coef(fit)[1] + coef(fit)[2] * x )  ) * x +
      ( survt * exp( coef(fit)[1] + coef(fit)[2] * x )  ) *
      cumsum(cumhaz0_tderiv( t, x )[ ,"cumhazderiv2"])[ which( diff_abs1(t) == min( diff_abs1(t) ) )]
  }


  ### gradient of p_s ###

  grad_psx  = cbind( b1_deriv( t, adj ), b2_deriv( t, adj ) )

  grad_unps = cbind( mean( b1_deriv( t, x = dat1$x ) ), mean( b2_deriv( t, x = dat1$x ) ) )

  var_psx   = grad_psx %*% vcov(fit) %*% t( grad_psx )

  var_unps  = grad_unps %*% vcov(fit) %*% t( grad_unps )

########### TRANS CIs ###########
  citr_l   =  max( fun_g( ps_x( t, adj ) + 1.96 * sqrt(var_psx) ), 1 )
  citr_r   =  fun_g( ps_x( t, adj ) - 1.96 * sqrt(var_psx) )

  cidl_l   =  max( NNT_X - 1.96 * NNT_X ^ 2 * sqrt(var_psx), 1 )
  cidl_r   =  NNT_X + 1.96 * NNT_X ^ 2 * sqrt(var_psx)

  citr_coxl   =  max( fun_g( mean( ps_x(t, dat1$x) ) + 1.96 * sqrt(var_unps) ), 1 )
  citr_coxr   =  fun_g( mean( ps_x(t, dat1$x) )      - 1.96 * sqrt(var_unps) )

  cidl_coxl   =  max( NNT_COX - 1.96 * NNT_COX ^ 2 * sqrt(var_unps), 1 )
  cidl_coxr   =  NNT_COX +      1.96 * NNT_COX ^ 2 * sqrt(var_unps)

######### NBS CIs #########

  nnt_1  = function(data, indices) {

    d        = data[indices,]

    fit      = coxph( Surv(time, status) ~  gr + x,
                      method = "breslow",
                      data   = d )

    # baseline hazard

    surv0_nbs    = round( basehaz( fit, centered = F ), 10 )

    # times of the baseline hazard

    dat_both_nbs = as.data.frame( cbind( time = surv0_nbs$time, surv_base = exp(- surv0_nbs$hazard) ) )

    diff_abs_nbs = abs( dat_both_nbs$time -  t  )

    pbs_x    = function(t, x){

      pbt_mle = dat_both_nbs$surv_base[ which( diff_abs_nbs == min( diff_abs_nbs ) ) ] ^ ( exp(  coef(fit)["gr"] + coef(fit)["x"] * x ) )

      pbc_mle = dat_both_nbs$surv_base[ which( diff_abs_nbs == min( diff_abs_nbs ) ) ] ^ ( exp( ( coef(fit)["x"] ) * x ) )

      return( pbt_mle - pbc_mle )
    }

    ps_km    = function(t){

    km_t   = survfit( Surv( d[d$gr == 1, "time"], d[d$gr == 1, "status"] ) ~ 1 )

    km_c   = survfit( Surv( d[d$gr == 0, "time"], d[d$gr == 0, "status"] ) ~ 1 )

    diff_abst = abs(t - km_t$time)
    diff_absc = abs(t - km_c$time)

    pbt_km = km_t$surv[ which( diff_abst == min( diff_abst ) ) ]

    pbc_km = km_c$surv[ which( diff_absc == min( diff_absc ) ) ]

    return( pbt_km - pbc_km )
    }

    av_psb  = mean( pbs_x(t, d$x) )

    NNT_KM_NBS   = fun_g( ps_km(t) )
    NNT_COX_NBS  = fun_g( av_psb )
    NNT_X_NBS    = fun_g( pbs_x(t, adj) )

    return(c( NNT_KM_NBS, NNT_COX_NBS, NNT_X_NBS ))

  }

  # bootstrapping with 1000 replications
  results <- boot(data      = dat1,
                  statistic = nnt_1,
                  R         = 1000)

  # 95% NBS confidence intervals
  cinbs_kml  = quantile(results$t[,1], probs = c(0.025))
  cinbs_kmr  = quantile(results$t[,1], probs = c(0.975))

  cinbs_coxl  = quantile(results$t[,2], probs = c(0.025))
  cinbs_coxr  = quantile(results$t[,2], probs = c(0.975))

  cinbs_l = quantile(results$t[,3], probs = c(0.025))
  cinbs_r = quantile(results$t[,3], probs = c(0.975))

########## PBS CIs ##########

  B = 30

  cum_haz0_t =  function(t, k){

    tab <- data.frame(table(dat1[dat1$status == 1, "time"]))
    y   <- as.numeric(levels(tab[, 1]))[tab[, 1]] # ordered distinct event times
    d   <- tab[, 2]

    h0 <- rep(NA, length(y))

    for(l in 1:length(y))
    {
      h0[l] <- d[l] / sum(exp(dat1[dat1$time >= y[l], "gr"] * beth[k, 1] +
                              dat1[dat1$time >= y[l], "x"]  * beth[k, 2]  ) )
    }

    diff_abs1 = function(t){  abs(t -  y) }

    return( cumsum(h0)[which( diff_abs1(t) == min( diff_abs1(t) ) )] )
  }

  ###  PBS CI LOOP ###

  nnt_pbsx   = numeric()

  nnt_pbscox = numeric()

    beth = mvrnorm(n     =  B,
                   mu    =  coef(coxph_both),
                   Sigma =  coxph_both$var,
                   empirical = F)

    ### nnt_bs matrix ###

    for( k in 1:B ) {

      ### PROGRESS BAR ###

      cat(paste0(round(k / B * 100), '% completed'))
      Sys.sleep(2)
      if (k == B) cat(" ", sep="\n")
      else cat('\014')

      survt_bsx = exp( - cum_haz0_t( t, k ) ) ^ ( exp( beth[k, 1] + beth[k, 2] * adj  ) )
      survc_bsx = exp( - cum_haz0_t( t, k ) ) ^ ( exp( beth[k, 2] * adj ) )

      survt_bsm = exp( - cum_haz0_t( t, k ) ) ^ ( exp( beth[k, 1] + beth[k, 2] * x ) )
      survc_bsm = exp( - cum_haz0_t( t, k ) ) ^ ( exp( beth[k, 2] * x ) )

      nnt_pbsx[k]     = fun_g( survt_bsx - survc_bsx )
      nnt_pbscox[k]   = fun_g( mean( survt_bsm - survc_bsm ) )
    }

    cipbs_coxl  =  max( NNT_COX - 1.96 * sd( nnt_pbscox, na.rm = T ), 1 )
    cipbs_coxr  =  NNT_COX + 1.96 * sd( nnt_pbscox, na.rm = T )

    cipbs_l = max( NNT_X - 1.96 * sd( nnt_pbsx, na.rm = T ), 1 )
    cipbs_r = NNT_X + 1.96 * sd( nnt_pbsx, na.rm = T )

##########################
### NNT-KM TR & DL CIs ###
##########################

    diff_abst = abs( t - km_t$time )
    diff_absc = abs( t - km_c$time )

    var_pskm = km_c$std.err[which(diff_absc == min( diff_absc ) ) ] ^ 2 + km_t$std.err[which(diff_abst == min( diff_abst ) ) ] ^ 2

    citr_kml    = fun_g( ps_km(t)      + 1.96 * sqrt( var_pskm ) )
    citr_kmr    = max( fun_g( ps_km(t) - 1.96 * sqrt( var_pskm ) ), 1 )

    cidl_kml =  max( NNT_KM - 1.96 * NNT_KM ^ 2 * sqrt( var_pskm ), 1 )
    cidl_kmr =  NNT_KM + 1.96 * NNT_KM ^ 2 * sqrt( var_pskm )

  ### FINAL CIs DATA FRAME ###

  out_list = list()

  out_list[[1]] = t(c(NNT_KM,
                      citr_kml,   citr_kmr,
                      cidl_kml,   cidl_kmr,
                      cinbs_kml,  cinbs_kmr,
                      cipbs_kml,  cipbs_kmr))

  out_list[[2]] = t(c(NNT_COX,
                      citr_coxl,  citr_coxr,
                      cidl_coxl,  cidl_coxr,
                      cinbs_coxl, cinbs_coxr,
                      cipbs_coxl, cipbs_coxr))


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

  rownames(out_list1) = c("NNT KM", "NNT COX", paste(c("NNT", "(", t, "|", adj, ")" ), sep = "", collapse = "") )


  return( out_list1 )

}

