\name{trimEiwild}
\alias{trimEiwild}
\title{Trims eiwild-object}
\usage{
trimEiwild(obj, burnin = 0, thinning = 1, sample = NULL)
}
\arguments{
  \item{obj}{object of type eiwild}

  \item{burnin}{number of draws to be cut away from the
  beginning of the chain. default=0}

  \item{thinning}{number specifying the thinning interval.
  default=1}

  \item{sample}{specifies sample size after burn-in and
  thinning (default is \code{NULL})}
}
\value{
eiwild-object
}
\description{
Trims eiwild-object with burn-in and thinning
}
\examples{
\dontrun{
# loading some fake election data
data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
set.seed(1234)
res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                sample=1000, thinning=2, burnin=100,verbose=100)

res
res2<- trimEiwild(res, burnin=100, thinning=3)
res2
}
}

