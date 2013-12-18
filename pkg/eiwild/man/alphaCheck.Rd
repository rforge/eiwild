\name{alphaCheck}
\alias{alphaCheck}
\title{alphaCheck function}
\usage{
alphaCheck(mat, r, c)
}
\arguments{
  \item{mat}{matrix with alpha values}

  \item{r}{number of rows in the RxC-table}

  \item{c}{number of cols in the RxC-table}
}
\value{
alphamatrix if everything is okay. \code{Error}-message if
something fails.
}
\description{
Given starting values for the alpha values of the dirichlet
distribution of the second level of the ecological
inference model \code{alphaCheck} checks if parameters are
correctly specified
}
\examples{
\dontrun{

# right alpha matrix
alphaRight <- matrix(1:9, 3,3)
alphaCheck(alphaRight, 3, 3)
# return is right alpha matrix

# wrong alpha matrix
alphaWrong <- matrix(0:8, 3,3)
alphaCheck(alphaWrong, 3,3)
alphaCheck(alphaRight, 3,4)
}
}

