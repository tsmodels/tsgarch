#' @rawNamespace useDynLib(tsgarch, .registration=TRUE); useDynLib(tsgarch_TMBExports)
#' @keywords internal
#' @importFrom TMB MakeADFun sdreport
#' @import tsmethods
#' @import data.table
#' @import methods
#' @importFrom stats nlminb na.omit printCoefmat ar coef dnorm integrate lm logLik residuals fitted pnorm var vcov confint qnorm sigma AIC BIC nobs simulate predict runif qqline qqplot
#' @importFrom zoo na.fill coredata index is.zoo `coredata<-` as.zoo
#' @importFrom xts xts as.xts is.xts merge.xts
#' @importFrom tsaux sampling_frequency check_xreg check_newxreg future_dates
#' @importFrom sandwich estfun bwNeweyWest vcovHAC vcovOPG bread
#' @importFrom numDeriv jacobian
#' @importFrom nloptr nloptr
#' @importFrom flextable flextable as_flextable set_caption italic fontsize separate_header add_footer_row add_footer_lines append_chunks as_chunk as_equation as_paragraph compose colformat_double set_header_labels padding bold align autofit hline width
#' @importFrom dplyr rename left_join
#' @importFrom generics glance tidy
#' @importFrom tibble as_tibble tibble
#' @importFrom graphics grid layout lines par
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom progressr handlers progressor
#' @importFrom car linearHypothesis
#' @importFrom tsdistributions ddist rdist qdist pdist
#' @importFrom Rdpack reprompt
#' @importFrom utils tail head data
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
