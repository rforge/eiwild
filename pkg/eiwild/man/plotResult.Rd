\name{plotResult}
\alias{plotResult}
\title{plots result of the beta table}
\usage{
plotResult(x, abs = FALSE, bgColors = c("white", "steelblue", 10),
  rectOpts = list(border = grey(0.7)), cellOpts = list(cex = 1),
  dimNameOpts = list(col = grey(0.2), cex = 0.8))
}
\arguments{
  \item{x}{\code{matrix} with results}

  \item{abs}{\code{TRUE} if values are not between [0,1]
  (Default is \code{FALSE})}

  \item{bgColors}{vector with 3 elements: starting colour
  of new colour palette, ending colour of new colour
  palette, length of colour palette}

  \item{rectOpts}{named \code{list} of options for
  rectangle}

  \item{cellOpts}{named \code{list} of options for inner
  cell text}

  \item{dimNameOpts}{named \code{list} of options for
  dimnames}
}
\description{
!!! Work in progress. not finished !!! only works for
betatables
}
\examples{
\dontrun{
# loading some fake election data
data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
set.seed(1234)
res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"), sample=1000, thinning=2, burnin=100,verbose=100)

tabs <- summary(res)
plotResult(round(tabs$relative,3))
plotResult(tabs$absolut, abs=TRUE)
bal <- getBalance(tabs$absolut, which=c("c","GRUN_2"))
plotResult(bal, abs=TRUE)

plotResult(round(tabs$relative,3), bgColors=c("white", "darkorange", 9))
plotResult(round(tabs$relative,3), bgColors=c("white", "darkorange", 5))

# ugly ;)
plotResult(round(tabs$relative,3), bgColors=c("blue", "red", 5))
#' }
}
\seealso{
\code{\link[eiwild]{getBalance}}
}

