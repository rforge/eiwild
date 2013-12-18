\name{getBalance}
\alias{getBalance}
\title{get profit and loss of partys}
\usage{
getBalance(x, rnd = 1, zero = TRUE, which = NULL)
}
\arguments{
  \item{x}{table of ecological inference (see
  \code{\link{summary.eiwild}})}

  \item{rnd}{rounding of values (default is \code{1})}

  \item{zero}{replace negative values with \code{0}
  (default=\code{TRUE})}

  \item{which}{if table isn't square it has to be a vector
  giving 2 arguments (see Details)}
}
\value{
table with balance BUT (!!!) order of rows or cols maybe
changed
}
\description{
get profit and loss of partys
}
\details{
if table isn't square the row or column not to be
calculated has to be given in \code{wich}.

First element has to be \code{"r"} for row or \code{"c"}
for column.

2nd element has to give name of row or column.
}
\seealso{
\code{\link{summary.eiwild}}
}

