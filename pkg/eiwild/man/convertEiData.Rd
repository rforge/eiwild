\name{convertEiData}
\alias{convertEiData}
\title{Extracting important values for ecological Inference}
\usage{
convertEiData(form, aggr, indi = NULL, IDCols = c("ID"))
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

  \item{IDCols}{vector of length 2 giving the column-names
  or numbers of ID column}
}
\value{
\code{list} with components needed for the Metropolis
algorithm in \code{\link[eiwild]{runMBayes}}
}
\description{
Extracting important values for calculation of the
ecological Inference with the
\code{\link[eiwild]{runMBayes}}-function
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
}
\examples{
\dontrun{
# loading some fake election data
data(topleveldat)
form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)

conv <- convertEiData(form=form, aggr=aggr, indi=indi, IDCols=c("ID", "ID"))

str(conv)
}
}
\seealso{
\code{\link[eiwild]{runMBayes}}, \code{\link[coda]{mcmc}}
\code{\link[eiwild]{tuneVars}}
}

