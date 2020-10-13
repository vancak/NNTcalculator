#' Data sample for the survival analysis example
#' @name survreg_data
#' @author V. Vancak
#' @docType data
#' @format A data frame with 1000 rows, and 4 columns. Suitable for the Cox regression model.
#' The file was generated using the \code{simple.surv.sim} function from the \code{survsim} package.
#' \describe{
#' \item{status}{is the status of each subject, where 1 corresponds to event, and 0 to censoring}
#' \item{stop}{is the time of the event/censoring}
#' \item{x}{is the allocated arm of each subject where 1 corresponds to treatment, and 0 to control}
#' \item{x.1}{is continuous explanatory variable}
#' }
#' @references Moriña, D., Navarro, A., & Soler, M. D. M. (2018). Package \code{survsim}.
#'
#' @keywords datasets
#'
#' @examples
#' data(survreg_data)
NULL
