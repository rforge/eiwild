\name{prioriPlot}
\alias{prioriPlot}
\title{plots \eqn{\beta^{rc}} given hyperpriori parameters}
\usage{
prioriPlot(pars, which, cols, alphaSample = 10000, betaSample = 300,
  plot = TRUE, ...)
}
\arguments{
  \item{pars}{list of parameters for hyperpriori. if
  which="gamma" then parameter has to be a list with shape
  and rate as parameters if which="expo" then parameter has
  to be a list with only lam}

  \item{which}{specified priori. "gamma" or "expo"}

  \item{cols}{integer specifying how many columns the
  RxC-Table should have}

  \item{alphaSample}{integer specifying the number of times
  new alpha-values are drawn}

  \item{betaSample}{integer specifying the number of times
  betas will be drawn for each alpha-value}

  \item{plot}{logical TRUE/FALSE if histogram should be
  plotted}

  \item{...}{additional arguments for "hist" function}
}
\value{
nested \code{list} with each element containing another
\code{list}. First level are rows and second level are
columns per row.
}
\description{
prioriPlot simulates the \eqn{beta_i^rc} values at first
level given specific parameters at hyperpriori level
}
\details{
Calculation is made via the marginal beta distribution

function structure: \itemize{ \item \code{"gamma"} choose
one parameter for every alpha_rc-parameter or a two
matrices of parameters specifying lambda's for every
alpha_rc-parameter \item \code{"expo"} choose one parameter
for every alpha_rc-parameter or a one matrix of parameters
specifying lambda's for every alpha_rc-parameter }
}
\examples{
\dontrun{
test1 <- prioriPlot(list(shape=4,rate=2), "gamma",cols=4)
str(test1)

pars <- list(shape=matrix(1:9,3,3),rate=matrix(9:1,3,3))
test2 <- prioriPlot(pars, "gamma",breaks=100)
test3 <- prioriPlot(list(shape=8,rate=2),"gamma",breaks=100,cols=3)

pars4 <- list(shape=matrix(c(6,6,6),1,3), rate=matrix(c(4,4,4),1,3))
test4 <- prioriPlot(pars4, "gamma",breaks=100)

pars5 <- list(lam=2)
test5 <- prioriPlot(pars5, "expo",cols=4, breaks=100)

pars6 <- list(lam=matrix(1:9,3,3)/100)
test6 <- prioriPlot(pars6, "expo", breaks=25, col=grey(0.8))

# example for 3x4-table
set.seed(568)
pars7 <- list(shape=matrix(sample(1:20,12), 3,4), rate=matrix(sample(1:20,12),3,4))
test7 <- prioriPlot(pars7, "gamma",breaks=50)
}
}

