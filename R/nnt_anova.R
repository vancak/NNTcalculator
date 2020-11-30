#' @title Laupacis' adjusted and marginal NNT in one-way ANOVA
#'
#' @description Internal function. Not for users. Calculates Laupacis' adjusted and marginal NNT in one-way ANOVA.
#' Takes a data-set suitable for a one-way ANOVA analysis and returns
#' the estimated adjusted NNT(k) and the marginal NNT.
#' @param response      vector of the response variable (i.e., the dependent variable).
#' @param x             vector of the explanatory variable, i.e., the factor of different treatments.
#' @param cutoff        the MCID threshold.
#' @param decrease      TRUE or FALSE. Indicates whether the MCID change is decrease in the response variable.
#' @param base          control group of the x variable in the one-way ANOVA model.
#' @param data          data frame that contains the required variables for the computations.
#'
#' @return The estimated marginal and adjusted NNT with their corresponding 95 percent confidence intervals.
#'
#' @keywords internal
#'
#' data(anova_data)
#'
#' ### SUCCESS = INCREASE
# nnt_anova(  response  = anova_data$y,
#             x         = anova_data$gr,
#             cutoff    = 0,
#             base      = 1,
#             decrease  = FALSE,
#             data      = anova_data)
#
# ### SUCCESS = DECREASE
# nnt_anova(     response  = anova_data$y,
#                x         = anova_data$gr,
#                cutoff    = 2,
#                base      = 4,
#                decrease  = TRUE,
#                data      = anova_data)
nnt_anova <- function( response,       # vector of the response variable
                       x,              # vector of the explanatory variable
                       cutoff,         # the MCID
                       base,           # control group of the x variable in the one-way ANOVA model.
                       decrease,       # whether the MCID change is decrease in the response variable (TRUE\FALSE)
                       data ) {        # the analyzed data-set

#  library(boot)

    ### initial values ###
    dat1    = data

    dat1    = data.frame( y = response, x = x )

    attach(dat1, warn.conflicts = F)

    fun_g   = function(x){ ifelse( x > 0, 1/x, Inf ) }

    tau     = cutoff

    ### control arm ###

    m1      = aov( y ~ relevel( factor( x ), ref = base ), data  = dat1 )

    sum     = summary.lm(m1)

    sig_m   = sum$sigma

    K       = 2:length(unique(x))

    ### INCREASE MARGINAL ###

    if( decrease == F ) {

      #####################
      ### ADJUSTED NNTs ###
      #####################

      ### ESTIMATION ###
      ps    = pnorm( ( tau - sum$coef[1] ) / sig_m ) - pnorm( ( tau - ( sum$coef[1] + sum$coef[K] ) ) / sig_m )

      NNT_K = fun_g( ps )

      #####################
      ### MARGINAL NNTs ###
      #####################

      ### ANOVA NNT ###
      av_ps    = 1 / ( length(unique(x)) - 1 ) * ( sum(ps) )

      NNT_MLE  = fun_g(av_ps)

      #########################
      ### ADJUSTED NNT CIs  ###
      #########################

      ### PARTIAL DERIVATIVES ###

      deriv_psb0  = dnorm( (tau - sum$coef[1]) / sig_m )  * ( - 1 / sig_m ) -
        dnorm( (tau - sum$coef[1] - sum$coef[K]) / sig_m )  * ( - 1 / sig_m )

      deriv_psbk  =  - dnorm( (tau - sum$coef[1] - sum$coef[K] )  / sig_m  ) * ( - 1 / sig_m )

      deriv_pssig =  dnorm( (tau - sum$coef[1]) / sig_m  ) * ( - ( tau - sum$coef[1]) / sig_m ^ 2 ) -
        dnorm( (tau - sum$coef[1] - sum$coef[K] )  / sig_m  )  * ( - (tau - sum$coef[1] - sum$coef[K] ) / sig_m ^ 2 )

      cov_b   = list()
      grad    = list()

      citr_l  = numeric()
      citr_r  = numeric()

      cidl_l  = numeric()
      cidl_r  = numeric()

      cinbs_l = numeric()
      cinbs_r = numeric()

      cipbs_l = numeric()
      cipbs_r = numeric()

      ### TR & DL CIs ###
      for( j in 1:length(K) ){

        grad[[j]]   = c( deriv_psb0[j], deriv_psbk[j], deriv_pssig[j] )

        cov_b[[j]]  =  rbind( c( vcov( m1 )[1,     c(1, (j+1))], 0 ),
                              c( vcov( m1 )[(j+1), c(1, (j+1))], 0 ),
                              c( 0, 0, 2 * sig_m ^ 4 / ( nrow(dat1) - ( length(K) + 1 ) ) ) )

        citr_l[j]   =  max( fun_g( ps[j]  + 1.96  * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] ) ), 1)
        citr_r[j]   =  fun_g( ps[j] - 1.96  * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] ) )

        cidl_l[j]   =  max( NNT_K[j] - 1.96  *  NNT_K[j] ^ 2 * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] ), 1 )
        cidl_r[j]   =  NNT_K[j]      + 1.96  *  NNT_K[j] ^ 2 * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] )

      }

      ### NBS CIs ###

      # function to obtain NNT(K) from the data
      nnt_1  = function(data, indices) {

        d       = data[indices,]

        fit     = aov(y ~ relevel( factor( x ), ref = base ), data = d)

        sumb    = summary.lm(fit)

        psb     = pnorm( ( tau - sumb$coef[1] ) / sumb$sigma ) - pnorm( ( tau - ( sumb$coef[1] + sumb$coef[K] ) ) / sumb$sigma )

        NNT_KB  = fun_g( psb )

        av_psb  = 1 / length(K) * sum( psb )

        NNT_UNB = fun_g(av_psb)

        return(c( NNT_KB, NNT_UNB ))
      }

      # bootstrapping with 1000 replications
      results <- boot(data      = dat1,
                      statistic = nnt_1,
                      R         = 1000)


      # 95% NBS confidence intervals

      for( j in 1:length(K) ) {

        cinbs_l[j] = quantile(results$t[,j], probs = c(0.025))
        cinbs_r[j] = quantile(results$t[,j], probs = c(0.975))

      }

      ### Parametric Bootstrap 95% CI ###
   #   library(mvtnorm)

      cov_pb =  cbind( rbind( vcov( m1 ), rep(0, (length(K)+1)  ) ) ,
                       c( rep(0, (length(K)+1) ),  2 * sig_m ^ 4 / ( nrow(dat1) - ( length(K) + 1 ) ) ) )

      pbs = rmvnorm(1000,
                    mean  = c( sum$coef[1:(length(K)+1)], sig_m ^ 2 ),
                    sigma =  cov_pb )

      pbs1  = pnorm( ( tau - pbs[,1] ) / sqrt( pbs[, (length(K)+2)] ) ) - pnorm( ( tau - ( pbs[,1] + pbs[,K] ) ) / sqrt( pbs[,(length(K)+2)] ) )

      NNT_PBS  = fun_g( pbs1 )

      av_ppbs   = 1 / length(K) * ( rowSums( pbs1 ) )

      NNT_UNPBS = fun_g(av_ppbs)

      # 95% NBS confidence intervals

      NNT_PBS    = NNT_PBS[is.finite(rowSums(NNT_PBS)),]

      NNT_UNPBS  = NNT_UNPBS[is.finite(NNT_UNPBS)]

      for( j in 1:length(K) ) {

        cipbs_l[j] = max( NNT_K[j] - 1.96 * sd(NNT_PBS[,j], na.rm = T), 1 )
        cipbs_r[j] = NNT_K[j]      + 1.96 * sd(NNT_PBS[,j], na.rm = T)

      }

      ########################
      ### MARGINAL NNT CIs ###
      ########################

      grad_ps = 1/length(K) * c( sum(deriv_psb0),
                                 deriv_psbk,
                                 sum(deriv_pssig) )

      grad_ml = - NNT_MLE ^ 2  * grad_ps


      cov_par = c( vcov( m1 )[1, 1:length(unique(x))], 0 )

      for( j in 2:length(unique(x)) ){

        cov_par   =   rbind( cov_par, c( vcov( m1 )[j, 1:length(unique(x))], 0 ) )

      }

      cov_bun  =  rbind( cov_par,
                         c( rep(0, length(unique(x)) ), 2 * sig_m ^ 4 / ( nrow(dat1) - ( length(unique(x)) )  ) ) )

      citr_mlel   =  max( fun_g( av_ps + 1.96  * sqrt( grad_ps %*% cov_bun %*% grad_ps ) ), 1 )
      citr_mler   =  fun_g( av_ps - 1.96  * sqrt(  grad_ps %*% cov_bun %*% grad_ps ) )

      cidl_mlel   =  max( NNT_MLE - 1.96  * sqrt( grad_ml %*% cov_bun %*% grad_ml ), 1 )
      cidl_mler   =  NNT_MLE      + 1.96  * sqrt( grad_ml %*% cov_bun %*% grad_ml )

      cinbs_mlel  = quantile(results$t[,(length(K)+1)], probs = c(0.025))
      cinbs_mler  = quantile(results$t[,(length(K)+1)], probs = c(0.975))

      cipbs_mlel   =  max( NNT_MLE - 1.96  * sd( NNT_UNPBS, na.rm = T ) , 1 )
      cipbs_mler   =  NNT_MLE      + 1.96  * sd( NNT_UNPBS, na.rm = T )

      ### MARGINAL NNT L CIs ###

      nntl = nnt_l( type      = "laupacis",
                    treat     = dat1[ x != base, 1 ],
                    control   = dat1[ x == base, 1 ],
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


      for( j in 1:length(K) ) {

        out_list[[j+2]] =  t(c(NNT_K[j],
                               citr_l[j],  citr_r[j],
                               cidl_l[j],  cidl_r[j],
                               cinbs_l[j], cinbs_r[j],
                               cipbs_l[j], cipbs_r[j]))

      }

      out_list1 = matrix( unlist( out_list ), ncol = 9, byrow = T )
      colnames(out_list1) = c("NNT",
                              "CI TR L",  "CI TR U",
                              "CI DL L",  "CI DL U",
                              "CI NBS L", "CI NBS U",
                              "CI PBS L", "CI PBS U"  )

      row_names = NULL

      for(j in 1:length(K)){

        row_names[j] = paste(c("NNT", "(", sort( unique( x[x != base ] ) )[j], ")" ), sep = "", collapse = "")

      }

      rownames(out_list1) = c("NNT L", "NNT MLE", row_names)

      return( out_list1 )

    }

###########################################################################################################################
##########################
### SUCCESS = DECREASE ###
##########################

    if( decrease == T ) {

      #####################
      ### ADJUSTED NNTs ###
      #####################

      ### ESTIMATION ###
      ps    =  - pnorm( ( tau - sum$coef[1] ) / sig_m ) + pnorm( ( tau - ( sum$coef[1] + sum$coef[K] ) ) / sig_m )

      NNT_K = fun_g( ps )

      #####################
      ### MARGINAL NNTs ###
      #####################

      ### ANOVA NNT ###
      av_ps    = 1 / ( length(unique(x)) - 1 ) * ( sum(ps) )

      NNT_MLE  = fun_g(av_ps)

      #########################
      ### ADJUSTED NNT CIs  ###
      #########################

      ### PARTIAL DERIVATIVES ###

      deriv_psb0  = - dnorm( (tau - sum$coef[1]) / sig_m )  * ( - 1 / sig_m ) +
        dnorm( (tau - sum$coef[1] - sum$coef[K]) / sig_m )  * ( - 1 / sig_m )

      deriv_psbk  =  - dnorm( (tau - sum$coef[1] - sum$coef[K] )  / sig_m  ) * ( - 1 / sig_m )

      deriv_pssig = - dnorm( (tau - sum$coef[1]) / sig_m  ) * ( - ( tau - sum$coef[1]) / sig_m ^ 2 ) +
        dnorm( (tau - sum$coef[1] - sum$coef[K] )  / sig_m  )  * ( - (tau - sum$coef[1] - sum$coef[K] ) / sig_m ^ 2 )

      cov_b   = list()
      grad    = list()

      citr_l  = numeric()
      citr_r  = numeric()

      cidl_l  = numeric()
      cidl_r  = numeric()

      cinbs_l = numeric()
      cinbs_r = numeric()

      cipbs_l = numeric()
      cipbs_r = numeric()

      ### TR & DL CIs ###
      for( j in 1:length(K) ){

        grad[[j]]   = c( deriv_psb0[j], deriv_psbk[j], deriv_pssig[j] )

        cov_b[[j]]  =  rbind( c( vcov( m1 )[1,     c(1, (j+1))], 0 ),
                              c( vcov( m1 )[(j+1), c(1, (j+1))], 0 ),
                              c( 0, 0, 2 * sig_m ^ 4 / ( nrow(dat1) - ( length(K) + 1 ) ) ) )

        citr_l[j]   =  max( fun_g( ps[j]  + 1.96  * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] ) ), 1)
        citr_r[j]   =  fun_g( ps[j] - 1.96  * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] ) )

        cidl_l[j]   =  max( NNT_K[j] - 1.96  *  NNT_K[j] ^ 2 * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] ), 1 )
        cidl_r[j]   =  NNT_K[j]      + 1.96  *  NNT_K[j] ^ 2 * sqrt( grad[[j]] %*% cov_b[[j]] %*% grad[[j]] )

      }

      ### NBS CIs ###

      # function to obtain NNT(K) from the data
      nnt_1  = function(data, indices) {

        d       = data[indices,]

        fit     = aov(y ~ relevel( factor( x ), ref = base ), data = d)

        sumb    = summary.lm(fit)

        psb     =  - pnorm( ( tau - sumb$coef[1] ) / sumb$sigma ) + pnorm( ( tau - ( sumb$coef[1] + sumb$coef[K] ) ) / sumb$sigma )

        NNT_KB  = fun_g( psb )

        av_psb  = 1 / length(K) * sum( psb )

        NNT_UNB = fun_g(av_psb)

        return(c( NNT_KB, NNT_UNB ))
      }

      # bootstrapping with 1000 replications
      results <- boot(data      = dat1,
                      statistic = nnt_1,
                      R         = 1000)


      # 95% NBS confidence intervals

      for( j in 1:length(K) ) {

        cinbs_l[j] = quantile(results$t[,j], probs = c(0.025))
        cinbs_r[j] = quantile(results$t[,j], probs = c(0.975))

      }

      ### Parametric Bootstrap 95% CI ###
  #    library(mvtnorm)


      cov_pb =  cbind( rbind( vcov( m1 ), rep(0, (length(K)+1)  ) ) ,
                       c( rep(0, (length(K)+1) ),  2 * sig_m ^ 4 / ( nrow(dat1) - ( length(K) + 1 ) ) ) )

      pbs = rmvnorm(1000,
                    mean  = c( sum$coef[1:(length(K)+1)], sig_m ^ 2 ),
                    sigma =  cov_pb )

      pbs1  = - pnorm( ( tau - pbs[,1] ) / sqrt( pbs[, (length(K)+2)] ) ) + pnorm( ( tau - ( pbs[,1] + pbs[,K] ) ) / sqrt( pbs[,(length(K)+2)] ) )

      NNT_PBS  = fun_g( pbs1 )

      av_ppbs   = 1 / length(K) * ( rowSums( pbs1 ) )

      NNT_UNPBS = fun_g(av_ppbs)

      # 95% NBS confidence intervals

      NNT_PBS    = NNT_PBS[is.finite(rowSums(NNT_PBS)),]

      NNT_UNPBS  = NNT_UNPBS[is.finite(NNT_UNPBS)]

      for( j in 1:length(K) ) {

        cipbs_l[j] = max( NNT_K[j] - 1.96 * sd(NNT_PBS[,j], na.rm = T), 1 )
        cipbs_r[j] = NNT_K[j]      + 1.96 * sd(NNT_PBS[,j], na.rm = T)

      }

      ########################
      ### MARGINAL NNT CIs ###
      ########################

      grad_ps = 1/length(K) * c( sum(deriv_psb0),
                                 deriv_psbk,
                                 sum(deriv_pssig) )

      grad_ml = - NNT_MLE ^ 2  * grad_ps


      cov_par = c( vcov( m1 )[1, 1:length(unique(x))], 0 )

      for( j in 2:length(unique(x)) ){

        cov_par   =   rbind( cov_par, c( vcov( m1 )[j, 1:length(unique(x))], 0 ) )

      }

      cov_bun  =  rbind( cov_par,
                         c( rep(0, length(unique(x)) ), 2 * sig_m ^ 4 / ( nrow(dat1) - ( length(unique(x)) )  ) ) )

      citr_mlel   =  max( fun_g( av_ps + 1.96  * sqrt( grad_ps %*% cov_bun %*% grad_ps ) ), 1 )
      citr_mler   =  fun_g( av_ps - 1.96  * sqrt(  grad_ps %*% cov_bun %*% grad_ps ) )

      cidl_mlel   =  max( NNT_MLE - 1.96  * sqrt( grad_ml %*% cov_bun %*% grad_ml ), 1 )
      cidl_mler   =  NNT_MLE      + 1.96  * sqrt( grad_ml %*% cov_bun %*% grad_ml )

      cinbs_mlel  = quantile(results$t[,(length(K)+1)], probs = c(0.025))
      cinbs_mler  = quantile(results$t[,(length(K)+1)], probs = c(0.975))

      cipbs_mlel   =  max( NNT_MLE - 1.96  * sd( NNT_UNPBS, na.rm = T ) , 1 )
      cipbs_mler   =  NNT_MLE      + 1.96  * sd( NNT_UNPBS, na.rm = T )

      ### MARGINAL NNT L CIs ###

      nntl = nnt_l( type      = "laupacis",
                    treat     = dat1[ x != base, 1 ],
                    control   = dat1[ x == base, 1 ],
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


      for( j in 1:length(K) ) {

        out_list[[j+2]] =  t(c(NNT_K[j],
                               citr_l[j],  citr_r[j],
                               cidl_l[j],  cidl_r[j],
                               cinbs_l[j], cinbs_r[j],
                               cipbs_l[j], cipbs_r[j]))

      }

      out_list1 = matrix( unlist( out_list ), ncol = 9, byrow = T )
      colnames(out_list1) = c("NNT",
                              "CI TR L",  "CI TR U",
                              "CI DL L",  "CI DL U",
                              "CI NBS L", "CI NBS U",
                              "CI PBS L", "CI PBS U"  )

      row_names = NULL

      for(j in 1:length(K)){

        row_names[j] = paste(c("NNT", "(", sort( unique( x[x != base ] ) )[j], ")" ), sep = "", collapse = "")

      }

      rownames(out_list1) = c("NNT L", "NNT MLE", row_names)

      return( out_list1 )

    } }

