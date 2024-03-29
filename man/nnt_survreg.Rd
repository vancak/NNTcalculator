% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nnt_survreg.R
\name{nnt_survreg}
\alias{nnt_survreg}
\title{NNT-KM and NNT-COX calculator}
\usage{
nnt_survreg(response, status, x, group, adj, time.point, data)
}
\arguments{
\item{response}{vector of the response variable; times of the events/censoring}

\item{status}{column that contains 0/1 indicator, where 1 is failure time, and 0 is censoring time.}

\item{x}{vector of the explanatory variable.}

\item{group}{allocated arm variable where 1 corresponds to the treatment arm, and 0 to the control arm.}

\item{adj}{the x value that the NNT(y|x) need to be conditioned on (adjusted for). The default value is the mean of x.}

\item{time.point}{the fixed time point y for NNT(y|x), and NNT(y)}

\item{data}{analyzed data frame that contains the required variables for the computations.}
}
\value{
The estimated unadjusted, harmonic mean (marginal) and conditinoal (adjusted) NNTs with their corresponding 95 percent confidence intervals.
}
\description{
Calculates Laupacis' type conditional (adjusted) and harmonic mean (marginal) NNT(y|x) for survival data.
Takes a data-set suitable for a survival analysis, a fixed time point y,
and an explanatory variables x,  and returns
the estimated conditional (adjusted) NNT(y|x), and the estimated unadjusted and harmonic mean (marginal) NNT(y).
In the Cox-model-based estimators, the baseline hazard is estimated by the Breslow's nonparametric MLE.
In the unadjusted estimator of NNT(y), the survival probabilities are estimated by the Kaplan-Meier nonparametric MLE.
}
\examples{
data(survreg_data)

nnt_survreg( response = survreg_data$stop,
             status   = survreg_data$status,
                   x  = survreg_data$x.1,
               group  = survreg_data$x,
                  adj = -1,
           time.point = 0.5,
                 data = survreg_data )
}
\references{
Vancak, V., Goldberg, Y., & Levine, S. Z. (2021). Guidelines to understand and compute the number needed to treat. Evidence-Based Mental Health, 24(4), 131-136.

Therneau, T. M., & Lumley, T. (2014). Package survival. Survival analysis Published on CRAN, 2, 3.

Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. Journal of the American statistical association, 53(282), 457-481.

Cox, D. R., & Oakes, D. (1984). Analysis of survival data (Vol. 21). CRC Press.

Cox, D. R. (1972). Regression models and life-tables. Journal of the Royal Statistical Society: Series B (Methodological), 34(2), 187-202.
}
