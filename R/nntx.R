#' @title Laupacis' adjusted and marginal NNT calculator in linear and generalized regression models.
#'
#' @description Calculates Laupacis' unadjusted, adjusted and marginal NNT.
#' Takes a data-set suitable for a regression analysis and returns
#' the estimated adjusted NNT(x), the estimated unadjusted, and the estimated marginal NNT
#' given the explanatory variable, and a specified model.
#' @param model         specification of the regression model; anova for the one-way ANOVA model,
#'  linreg for the linear regression,
#'  and logreg for the logistic regression with the logit link-function.
#' @param response      vector of the response variable (i.e., the dependent variable).
#' @param x             vector of the explanatory variable.
#' @param cutoff        the MCID threshold. This argument is suitable for continuous response variables,
#'  namely for ANOVA and linear regression.
#' @param decrease      TRUE or FALSE. Indicates whether the MCID change is decrease in the response variable
#' @param group         allocated arm variable where 1 corresponds to the treatment arm, and 0 to the control
#' arm. Suitable for linear and logistic regression.
#' @param adj           value that the NNT need to be adjusted for. The default value is mean of x.
#' @param base          control group of the x variable in the one-way ANOVA model.
#' @param data          data frame that contains the required variables for the computations.
#'
#' @return The estimated unadjusted, marginal and adjusted NNT with their corresponding 95 percent confidence intervals
#' given a specified model, and adjusted for a specified value of the explanatory variables.
#' @examples
#' data(anova_data)
#'
#' ### SUCCESS = INCREASE
#' nnt_x(         model     = "anova",
#'                response  = anova_data$y,
#'                x         = anova_data$gr,
#'                cutoff    = 0,
#'                base      = 1,
#'                decrease  = FALSE,
#'                data      = anova_data)
#'
#' ### SUCCESS = DECREASE
#' nnt_x(         model     = "anova",
#'                response  = anova_data$y,
#'                x         = anova_data$gr,
#'                cutoff    = 2,
#'                base      = 4,
#'                decrease  = TRUE,
#'                data      = anova_data)
#'
#' data(linreg_data)
#'
#' ### SUCCESS = INCREASE
#' nnt_x(         model    = "linreg",
#'                response = linreg_data$y,
#'                x        = linreg_data$x_var,
#'                cutoff   = 3,
#'                group    = linreg_data$gr,
#'                decrease = FALSE,
#'                adj      = 2.6,
#'                data     = linreg_data )
#'
#' ### SUCCESS = DECREASE
#' inv_data = data.frame( y = linreg_data$y, x_var = linreg_data$x_var, gr = 1 - linreg_data$gr )
#'
#' nnt_x(  model    = "linreg",
#'        response  = inv_data$y,
#'         x        = inv_data$x_var,
#'         cutoff   = 3,
#'         group    = inv_data$gr,
#'         decrease = TRUE,
#'         adj      = 2.6,
#'         data     = inv_data )
#'
#' data(logreg_data)
#'
#' nnt_x( model    = "logreg",
#'        response = logreg_data$y,
#'        x        = logreg_data$x_var,
#'        group    = logreg_data$gr,
#'        adj      = 1.5,
#'        data     = logreg_data )
#'
#' @export
#'
nnt_x <- function( model,          # regression model; anova, linreg, logreg
                   response,       # vector of the response variable
                   x,              # vector of the explanatory variable
                   cutoff,         # the MCID
                   base,           # control group of the x variable in the one-way ANOVA model
                   group,          # the allocated arm variable where 1 corresponds to the treatment arm, and 0 to the control arm
                   adj,            # value that the NNT need to be adjusted for. The default value is mean(x)
                   decrease,       # whether the MCID change is decrease in the response variable (TRUE\FALSE)
                   data ) {        # the analyzed data-set


  if( !(model %in% c("anova", "linreg", "logreg") ) ) warning("type can be 'anova', 'linreg' or 'logreg' only")
  if( !is.numeric(response) )                         warning("response vector must be a numeric vector")
#  if( !is.numeric(cutoff) | length(cutoff) > 1 )      warning("cutoff must be a scalar")
#  if( !(decrease %in% c(T, F, TRUE, FALSE) ) )        warning("decrease must be TRUE or FALSE")
#  if( model %in% c("anova") & !(base %in% unique(x)) ) warning("base must be one the levels of x")
#  if( model == "logreg" )                             warning("if `cutoff` and/or `decrease` arguments were set they are ignored")
#  if( model == "anova" & ncol(x) > 1 )                warning("in one-way-ANOVA `x` must be column of factors")
#  if( model == "anova" & !(is.factor(x))  )           warning("in one-way-ANOVA `x` must be factor")
#  if( model == "logreg" & unique(response) > 2 )      warning("in logistic regression the response variable must be binary vector of 1 and 0")
#  if( length(response) != nrow(x) )                   warning("response vector and number of rows of `x` must be the same")

#  library(nntcalc)

  if( model == "anova" ){

    df = nnt_anova( response  = response,
                    x         = x,
                    cutoff    = cutoff,
                    base      = base,
                    decrease  = decrease,
                    data      = data )
    }

  if( model == "linreg" ){

  df = nnt_lm( response = response,
               x        = x,
               cutoff   = cutoff,
               group    = group,
               decrease = decrease,
               adj      = adj,
               data     = data )
  }

  if( model == "logreg" ){

  df =  nnt_logreg( response = response,
                    x        = x,
                    group    = group,
                    adj      = adj,
                    data     = data )

  }

  return(df)

  }

