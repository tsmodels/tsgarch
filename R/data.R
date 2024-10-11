#' Deutschemark/British pound Exchange Rate
#'
#' The Bollerslev-Ghysel benchmark dataset. The variables in the data set are
#' the daily percentage nominal returns computed as 100 \deqn{ln(Pt) - ln(Pt-1)},
#' where Pt is the bilateral Deutschemark/British pound rate constructed from
#' the corresponding U.S. dollar rates, and a dummy variable that takes the
#' value of 1 on Mondays and other days following no trading in the Deutschemark
#' or British pound/ U.S. dollar market during regular European trading hours,
#' and 0 otherwise.
#' The data spans the period from 1984-01-03 through 1991-12-31, but exact dates
#' are not known as this dataset did not provide an index.
#' This dataset is included as it is used for the GARCH benchmark.
#' @format ## `dmbp`
#' A data.frame containing 2x1974 observations
#' \describe{
#'   \item{rate}{The exchange rate}
#'   \item{monday}{Dummy indicator (see descriptiom)}
#' }
#' @source Journal of Business & Economic Statistics Data Archive
"dmbp"


#' Japanese NIKKEI Stock Index
#'
#' The daily log returns in percent of the NIKKEI stock index spanning the
#' period 1984-01-04 to 2000-12-22. In the original dataset there was a
#' duplicate date on 2000-08-31 (but with a different value). Therefore,
#' in order to correct this we have moved up the duplicate 2000-08-31 to
#' become 2000-09-01, and the 2000-09-01 to 2000-09-02. Since the next date
#' after this was 2000-09-04, no further adjustments were made. These changes
#' preserve the original data in the order they appeared, with a minimal
#' adjustment only to the index which has no impact on estimation, but avoiding
#' internal warnings which arise on checks to the index.
#' This dataset is included as it is used for the APARCH benchmark.
#' @format ## `nikkei`
#' A data.frame containing 4246 observations in 2 columns:
#' \describe{
#'   \item{index}{The string date in YYYY-MM-DD format.}
#'   \item{value}{The daily log returns}
#' }
#' @source Journal of Applied Econometrics Data Archive 2003-v18.6/giot-laurent from
#' the paper \dQuote{Value-at-Risk for Long and Short Trading Positions} by Giot, Pierre
#' and Sebastien Laurent, 2003, \emph{Journal of Applied Econometrics}, 18(6), pp. 641--664.
"nikkei"
