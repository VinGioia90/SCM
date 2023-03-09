\name{GEF14_d4}
\alias{GEF14_d4}
\docType{data}
\title{ GEFCom14 Electricity demand data}
\description{Data set on electricity demand from USA. The data are from the electricity load forecasting track of the GEFCom2014 challenge (Hong et al., 2016). Here, only the hourly loads (in MW) from 17 p.m. to 20 p.m., spanning the period from 2005/01/02 to 2011/11/30, are considered.
}
\usage{data("GEF14_d4")}
\format{
  A data frame with 2520 observations on the following 21 variables.
  \describe{
    \item{\code{load_h17}}{hourly loads at 17 p.m.}
    \item{\code{load_h18}}{hourly loads at 18 p.m.}
    \item{\code{load_h19}}{hourly loads at 19 p.m.}
    \item{\code{load_h20}}{hourly loads at 20 p.m.}
    \item{\code{load24_h17}}{hourly loads at 17 p.m. of the previous day}
    \item{\code{load24_h18}}{hourly loads at 18 p.m. of the previous day}
    \item{\code{load24_h19}}{hourly loads at 19 p.m. of the previous day}
    \item{\code{load24_h20}}{hourly loads at 20 p.m. of the previous day}
    \item{\code{temp_h17}}{temperature at 17 p.m.}
    \item{\code{temp_h18}}{temperature at 18 p.m.}
    \item{\code{temp_h19}}{temperature at 19 p.m.}
    \item{\code{temp_h20}}{temperature at 20 p.m.}
    \item{\code{temp95_h17}}{exponentially smoothed temperature at 17 p.m.}
    \item{\code{temp95_h18}}{exponentially smoothed temperature at 18 p.m.}}
    \item{\code{temp95_h19}}{exponentially smoothed temperature at 19 p.m.}}
    \item{\code{temp95_h20}}{exponentially smoothed temperature at 20 p.m.}}
    \item{\code{doy}}{day of the year}
    \item{\code{dow}}{day of the week: a factor with levels \code{0} \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}}
    \item{\code{month}}{month }
    \item{\code{year}}{year}
    \item{\code{progtime}}{progressive time}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Hong, T., Pinson, P., Fan, S., Zareipour, H., Troccoli, A., and Hyndman, R.J. (2016). Probabilistic energy forecasting: global energy forecasting competition 2014 and beyond. International Journal of Forecasting, 32, 896 -- 913.
}
\examples{
data(GEF14_d4)
str(GEF14_d4)
}
\keyword{datasets}
