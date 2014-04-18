\name{indAggEi}
\alias{class-eiwild}
\alias{indAggEi}
\title{Ecological Inference including Individual Data}
\usage{
indAggEi(form, aggr, indi = NULL, IDCols = c("ID"), whichPriori = "gamma",
  prioriPars = list(shape = 4, rate = 2), startValsAlpha = NULL,
  startValsBeta = NULL, betaVars = NULL, alphaVars = NULL, sample,
  burnin = 0, thinning = 1, verbose = 1, retBeta = FALSE, seed = NULL)
}
\arguments{
  \item{form}{\code{formula} in this format
  \code{cbind(column_1,column_2,
  ...,column_c)~cbind(row_1,row_2,...,row_r))}}

  \item{aggr}{\code{data.frame} with aggregate data.  One
  district per line and one column giving one ID per
  district. (see Details)}

  \item{indi}{\code{data.frame} with individual data.  One
  district per line and one column giving one ID per
  district. (see Details) If no individual data are present
  it defaults to \code{NULL}}

  \item{IDCols}{vector of length 2 (or 1) giving the
  column-names or numbers of ID column.  First element for
  \code{aggr} and second element for \code{indi}.}

  \item{whichPriori}{character string specifying the chosen
  hyperpriori.  Options are \code{"gamma"} or \code{"expo"}
  (see Details)}

  \item{prioriPars}{vector or matrix of parameters for the
  specified hyperpriori in \var{whichPriori}}

  \item{startValsAlpha}{matrix with dimension
  \code{c(rows,columns)} giving the starting values for
  alpha.  If \code{NULL} random numbers of rdirichlet with
  chosen hyperpriori will be chosen.}

  \item{startValsBeta}{array with dimension
  \code{c(rows,columns,districts)} giving the starting
  values of beta If \code{NULL} random multinomial numbers
  will be chosen.}

  \item{sample}{the sample size to be saved in output.
  Total length of chain will be \var{burnin} \code{+}
  \var{sample} \code{*} \var{thinning}}

  \item{burnin}{number of draws to be cut away from the
  beginning of the Markov-Chain. \code{default=0}}

  \item{thinning}{number specifying the thinning interval.
  \code{default=1}}

  \item{verbose}{an integer specifying whether the progress
  of the sampler is printed to the screen (defaults to 0).
  If verbose is greater than 0, the iteration number is
  printed to the screen every verboseth iteration}

  \item{betaVars}{array-object with dimensions \code{(rows,
  columns-1, districts)} giving variance of proposal
  density for \eqn{\beta}-values}

  \item{alphaVars}{matrix of dimensions \code{(rows,
  columns)} giving variance of proposal density for
  \eqn{\alpha}-values.}

  \item{retBeta}{logical \code{TRUE} if estimated
  \eqn{\beta}-parameters should be returned.  With large
  number of precincts there can be problems with memory}

  \item{seed}{Default is \code{NULL}. Can be given the
  \code{"seed"}-attribute of an \code{eiwild}-object to
  reproduce an \code{eiwild}-object}
}
\value{
object of class \code{"eiwild"} which is a nested
\code{list} with elements: \itemize{ \item \code{draws}
\itemize{ \item \code{alphaDraws}
\code{\link[coda]{mcmc}}-object \item \code{cellCounts}
\code{\link[coda]{mcmc}}-object \item \code{betaDraws}
\code{\link[coda]{mcmc}}-object \item \code{betaAcc}
\code{"numeric"} with Acceptance ratios \item
\code{alphaAcc} \code{"numeric"} with Acceptance ratios
\item \code{alphaVars} \code{matrix} with variances for
proposal density \item \code{betaVars} \code{array} with
variances for proposal density } \item \code{rowdf}
original aggregate data \item \code{coldf} original
aggregate data } An Attribute called \code{"seed"} is also
saved to reproduce the eiwild-object
}
\description{
\code{indAggEi} (Individual and Aggregate Ecological
Inference) Calculating ecological Inference including
Individual data using the two functions
\code{\link[eiwild]{convertEiData}} and
\code{\link[eiwild]{runMBayes}}
}
\details{
\code{indi} is a \var{districts x [(r*c)+1]}
\code{data.frame} containing one district per line. One
column gives the ID of the districts which will be
connection to the ID column in the \code{aggr}-data.frame.

For example a 2x3 ecological Inference problem with
\code{formula} \code{cbind(col1,col2,col3) ~
cbind(row1,row2)} will have the row format : \code{[ID,
row1.col1, row1.col2, row1.col3, row2.col1, row2.col2,
row2.col3]}

It is important that the \code{formula} names correspond to
the exact column number in the \var{indi}-data.frame.

The \var{aggr} data.frame can have more columns than the
names given in \var{formula} as long as the colnames exist

The \code{whichPriori}-parameter has the options
\code{"gamma"} or \code{"expo"} and corresponding
\code{prioriPars}-parameters in a \code{"list"}: \itemize{
\item \code{"expo"} and \code{numeric} list-element called
\code{"lam"} corresponding to: \eqn{\alpha_{rc} \sim
Exp(\lambda)} \item \code{"expo"} and \code{matrix}
list-element called \code{"lam"} corresponding to:
\eqn{\alpha_{rc} \sim Exp(\lambda_{rc})} \item
\code{"gamma"} and two \code{numeric} list-element called
\code{"shape"} and \code{"rate"} corresponding to:
\eqn{\alpha_{rc} \sim Gamma(\lambda_1, \lambda_2)} \item
\code{"gamma"} and two \code{matrix} list-element called
\code{"shape"} and \code{"rate"} corresponding to:
\eqn{\alpha_{rc} \sim Gamma(\lambda_1^{rc},
\lambda_2^{rc})} }

The \code{"seed"} attribute is generated by the
\code{\link{.Random.seed}}-function.
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
summary(res)

# with individual alpha-hyperpriori-parameters
hypMat <- list(shape = matrix(c(30,4,4,4,
                                4,30,4,4,
                                4,4,30,4), nrow=3, ncol=4, byrow=TRUE),
               rate = matrix(c(1,2,2,2,
                               2,1,2,2,
                               2,2,1,2), nrow=3, ncol=4, byrow=TRUE))
set.seed(12345)
res2 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
                 sample=1000, thinning=2, burnin=100, verbose=100,
                 prioriPars=hypMat, whichPriori="gamma")
}
}
\seealso{
\code{\link[eiwild]{convertEiData}},
\code{\link[eiwild]{runMBayes}}, \code{\link[coda]{mcmc}}
\code{\link[eiwild]{tuneVars}}
}

