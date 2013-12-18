\name{betaCheck}
\alias{betaCheck}
\title{betaCheck function}
\usage{
betaCheck(arr, r, c, prec)
}
\arguments{
  \item{arr}{array with beta values. 1st dimension: rows,
  2nd dimension columns, 3rd dimension precincts}

  \item{r}{number of rows in the RxC-table}

  \item{c}{number of cols in the RxC-table}

  \item{prec}{number of precincts}
}
\value{
betaarry if everything is okay. \code{Error}-message if
something fails.
}
\description{
Given starting values for the beta values on the first
level of the ecological inference model, \code{betaCheck}
checks if parameters are correctly specified
}
\examples{
\dontrun{
# right beta array
beta1 <- rep(c(0,0.25,0.75), each=3)
beta2 <- rep(beta1, 3)
betaRight <- array(beta2, dim=c(3,3,3))
betaCheck(betaRight, 3,3,3)

# wrong beta array
betaWrong <- array(1:27, dim=c(3,3,3))
betaCheck(betaWrong)
}
}

