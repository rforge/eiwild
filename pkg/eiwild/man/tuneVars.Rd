\name{tuneVars}
\alias{tuneVars}
\title{Generate Variance parameters for proposal density}
\usage{
tuneVars(form, aggr, indi = NULL, IDCols = c("ID"), whichPriori = "gamma",
  prioriPars = list(shape = 4, rate = 2), accRat = c(0.4, 0.6),
  minProp = 0.7, maxiter = 20, sample = 10000, verbose = 10000,
  verboseTune = TRUE, improv = NULL, betaVars = NULL, alphaVars = NULL,
  startValsAlpha = NULL, startValsBeta = NULL, seed = NULL)
}
\arguments{
  \item{form}{formula in this format
  cbind(column_1,column_2,
  ...,column_c)~cbind(row_1,row_2,...,row_r))}

  \item{aggr}{data.frame with aggregate data. One district
  per line and one column giving one ID per district. (see
  Details)}

  \item{indi}{data.frame with individual data. One district
  per line and one column giving one ID per district. (see
  Details) If no individual data are present it defaults to
  NULL}

  \item{IDCols}{vector of length 2 (or 1) giving the
  columnnames or numbers of ID column.  First element for
  \code{aggr} and second element for \code{indi}.}

  \item{whichPriori}{character string defining the
  hyperpriori. default="gamma"}

  \item{prioriPars}{vector giving the parameters of the
  hyperpriori}

  \item{accRat}{vector with two elements describing the
  wished range of the acceptance ratios}

  \item{minProp}{numeric between 0 and 1 describing the
  percentage of parameters to have the wished acceptance
  ratios (\code{accRat}). \code{maxiter} will be the
  maximum iteration}

  \item{maxiter}{numeric how many times the algorithm
  should run maximum. If NULL tuning will run until minProp
  is reached}

  \item{sample}{the sample size to be drawn each tuning
  run.}

  \item{verbose}{an integer specifying whether the progress
  of the sampler is printed to the screen (defaults to 0).
  If verbose is greater than 0, the iteration number is
  printed to the screen every verboseth iteration}

  \item{verboseTune}{logical if tuning iteration should be
  printed (default=TRUE)}

  \item{improv}{numeric vector with 2 elements c(a,b).
  standard deviation will be calculated with the last b
  percentages of parameters to have the wished acceptance
  ratio.  If standard deviation is lower than a, than
  tuning is finished. Default is NULL}

  \item{betaVars}{array of dimensions (rows, columns,
  districts) giving variance of proposal density for
  betavalues}

  \item{alphaVars}{matrix of dimensions (rows, columns)
  giving variance of proposal density for alphavalues.}

  \item{startValsAlpha}{matrix with
  dimension=c(rows,columns) giving the starting values for
  alpha. If \code{NULL} random numbers of rdirichlet with
  prioriPars will be drawn}

  \item{startValsBeta}{array with
  dimension=c(rows,columns,districts) giving the starting
  values of beta If \code{NULL} random multinomial numbers
  with startValsAlpha or prioriPars will be draws}

  \item{seed}{Default is NULL. Can be given the "seed"
  attribute of an eiwild-object to reproduce an
  eiwild-object}
}
\value{
A list containing matrices of variance parameters for the
proposal densities
}
\description{
\code{tuneVars} generates variance parameters dependent on
the acceptance ratio of drawn parameters in the Metropolis
algorithm used in \code{indAggEi}. The target acceptance
rate will be given by parameter \code{accRat}
}
\details{
\code{indi} is a districts x [(r*c)+1] data.frame
containing one district per line. One column gives the ID
of the districts. This will we connected to the ID column
in the \var{aggr}-data.frame. The rest of one line in
\var{indi} is every row beside the nex row. For example a
2x3 ecological Inference problem with \var{formula}
\code{cbind(col1,col2,col3) ~ cbind(row1,row2)} will have
the row format : \code{[ID, row1.col1, row1.col2,
row1.col3, row2.col1, row2.col2, row2.col3]}

It is important that the formula names correspond to the
exact column number in the \var{indi}-data.frame.

The \var{aggr} data.frame can have more columns than the
names given in \var{formula} as long as the colnames exist

Priorities for finishing of tuning are as follows: If
\var{improv} isn't specified: \var{minProp} and
\var{maxiter} are checked. If \var{improv} is specified: 1)
\var{improv} is checked, 2) \var{minProp} and \var{maxiter}
are checked.
}
\examples{
\dontrun{
data(topleveldat)
out1 <- tuneVars(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),sample=10000, verbose=11000)
out2 <- tuneVars(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"), sample=10000, verbose=11000,
                 maxiter=NULL, improv=c(0.01,5))
out3 <- tuneVars(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"), sample=10000, verbose=11000,
                 maxiter=NULL, accRat=c(0.45,0.55), improv=c(0.01,5))
str(out3)
out4 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                 betaVars=out1$betaVars, alphaVars=out1$alphaVars,
                 sample=10000,thinning=1,burnin=100, verbose=1000)
out4
}
}
\seealso{
\code{\link[eiwild]{convertEiData}},
\code{\link[eiwild]{runMBayes}}, \code{\link[coda]{mcmc}}
\code{\link[eiwild]{tuneVars}},
\code{\link[eiwild]{indAggEi}}
}

