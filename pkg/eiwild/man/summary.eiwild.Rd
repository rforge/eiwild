\name{summary.eiwild}
\alias{summary.eiwild}
\title{summary method for \code{eiwild}-object}
\usage{
\method{summary}{eiwild}(object, cred = 0.95, ...)
}
\arguments{
  \item{object}{\code{eiwild}-object}

  \item{cred}{percentage for credibility interval of alphas
  and cellCounts}

  \item{...}{no function. included for S3 generic/method
  consistency}
}
\value{
tables and matrices \itemize{ \item \code{relative:} global
beta values calculated with \code{cellCounts} \item
\code{absolut:} \code{cellCounts} mean \item
\code{alphaMeans:} Means of \code{alphaDraws} \item
\code{relativeCol:} proportions with \code{colSum=1} \item
\code{countsCred:} Credibility Interval of length
\code{cred} for \code{cellCounts} \item \code{alphaCred:}
Credibility Interval of length \code{cred} for
\code{alphaDraws} \item \code{realtiveCred:} Credibility
Interval of length \code{cred} for global beta values
calculated with \code{cellCounts} }
}
\description{
summary method for \code{eiwild}-object
}
\examples{
\dontrun{
# loading some fake election data
data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
set.seed(1234)
res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"), sample=1000, thinning=2, burnin=100,verbose=100)

res
sumRes <- summary(res)
sumRes
str(sumRes)
}
}

