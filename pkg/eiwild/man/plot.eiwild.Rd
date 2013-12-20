\name{plot.eiwild}
\alias{plot.eiwild}
\title{Plot diagnostics for an eiwild Object}
\usage{
\method{plot}{eiwild}(x, whichPlot = NULL, whichParam = "cellCounts",
  whichCell = NULL, layout = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{eiwild}}

  \item{whichPlot}{which type of plot. see Details for more
  information}

  \item{whichParam}{which parameter should be plotted
  \code{"alphaDraws"} or \code{"cellCounts"}}

  \item{whichCell}{which cell to plot}

  \item{layout}{logical if automatic layout of plot should
  be made with help of \code{par=mfrow()}}

  \item{...}{further graphical parameters given to
  correspondent coda function}
}
\description{
plot diagnostics for an eiwild Object containing
autocorrelation plots, density plots and trace using
functions of coda package
}
\details{
\code{plot.eiwild} uses the plot diagnostic functions of
the coda package. The default is \code{NULL} which passes
arguments to \code{plot.mcmc}: \itemize{ \item \code{1}
passes arguments to \code{\link[coda]{traceplot}} \item
\code{2} passes arguments to
\code{\link[coda]{autocorr.plot}} \item \code{3} passes
arguments to \code{\link[coda]{densplot}} \item \code{4}
calculates and plots rolling mean }
}
\examples{
\dontrun{
# loading some fake election data
data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
set.seed(1234)
res <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                sample=1000, thinning=2, burnin=100,verbose=100)

plot(res, whichPlot=1)
plot(res, whichPlot=2)
plot(res, whichPlot=3)
plot(res, whichPlot=4)

plot(res, whichPlot=1, whichCell=c(1,4,5))
plot(res, whichPlot=2, whichCell=c(1,4,5))
plot(res, whichPlot=3, whichCell=c(1,4,5))
plot(res, whichPlot=4, whichCell=c(1,4,5))

plot(res, whichPlot=1, whichCell=c(1))
plot(res, whichPlot=2, whichCell=c(1))
plot(res, whichPlot=3, whichCell=c(1))
plot(res, whichPlot=4, whichCell=c(1))

plot(res, whichPlot=1, whichParam="alphaDraws")
plot(res, whichPlot=2, whichParam="alphaDraws")
plot(res, whichPlot=3, whichParam="alphaDraws")
plot(res, whichPlot=4, whichParam="alphaDraws")

par(mfrow=c(2,2))
plot(res, whichPlot=1, whichCell=1, layout=FALSE)
plot(res, whichPlot=2, whichCell=1, layout=FALSE)
plot(res, whichPlot=3, whichCell=1, layout=FALSE)
plot(res, whichPlot=4, whichCell=1, layout=FALSE)
}
}
\seealso{
\code{\link[coda]{plot.mcmc}} \code{\link[coda]{mcmc}}
}

