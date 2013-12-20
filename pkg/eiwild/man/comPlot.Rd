\name{comPlot}
\alias{comPlot}
\title{Plot Diagnostics to compare eiwild objects}
\usage{
comPlot(eiList, whichCell, whichPlot, whichParam = "cellCounts",
  rollLim = NULL, rollCol = NULL, ...)
}
\arguments{
  \item{eiList}{list of \code{eiwild} objects}

  \item{whichCell}{which cell to plot}

  \item{whichPlot}{which type of plot. see Details for more
  information}

  \item{whichParam}{which parameter should be plotted
  \code{"alphaDraws"} or \code{"cellCounts"}}

  \item{rollLim}{specifying other \code{ylim}-values for
  rolling mean plot (Default=\code{NULL})}

  \item{rollCol}{specifying other \code{col}-values for
  rolling mean plot (Default=\code{NULL})}

  \item{...}{arguments given to corresponding coda
  function}
}
\description{
This function uses plot diagnostics (see
\code{\link[eiwild]{plot.eiwild}}) to compare different
eiwild objects using functions of coda package
}
\details{
\code{whichPlot} controls the plot diagnostic to run:
\itemize{ \item \code{1} passes arguments to
\code{\link[coda]{traceplot}} \item \code{2} passes
arguments to \code{\link[lattice]{densityplot}} \item
\code{3} calculates Running Mean with
\code{eiwild:::rollMean} \item \code{4} passes arguments to
\code{\link[coda]{gelman.plot}}.  Output of
\code{\link[coda]{gelman.diag}} will be title of this plot.
}
}
\examples{
\dontrun{
# loading some fake election data
data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
set.seed(1234)
out1 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                sample=1000, thinning=2, burnin=100, verbose=100)
out2 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                sample=1000, thinning=2, burnin=100, verbose=100)
out3 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                sample=1000, thinning=2, burnin=100, verbose=100)
out4 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                sample=1000, thinning=2, burnin=100, verbose=100)

eiList <- list(out1, out2, out3, out4)

comPlot(eiList, whichCell=1, whichPlot=1)
comPlot(eiList, whichCell="counts.CSU_1.CSU_2", whichPlot=1)
comPlot(eiList, whichCell=1, whichPlot=1, smooth=TRUE)

comPlot(eiList, whichCell=1, whichPlot=2)

comPlot(eiList, whichCell=1, whichPlot=3)

comPlot(eiList, whichCell=1, whichPlot=4)
comPlot(eiList, whichCell=1, whichPlot=4)
comPlot(eiList, 1, 3, whichParam="alphaDraws")

comPlot(eiList, "alpha.CSU_1.CSU_2", 3, whichParam="alphaDraws")
}
}
\seealso{
\code{\link[coda]{mcmc}} \code{\link[eiwild]{plot.eiwild}}
\code{\link[eiwild]{indAggEi}}
}

