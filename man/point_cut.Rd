% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10_survival_cut_cox_km.R
\name{point_cut}
\alias{point_cut}
\title{point_cut}
\usage{
point_cut(exprSet_hub, meta)
}
\arguments{
\item{exprSet_hub}{a tumor expression set for hubgenes}

\item{meta}{meta data corresponds to expression set}
}
\value{
a vector with cutpoint for genes
}
\description{
calculate cut point for multiple genes
}
\examples{
point_cut(exprSet_hub1,meta1)
}
\seealso{
\code{\link{surv_KM}};\code{\link{surv_cox}}
}
\author{
Xiaojie Sun
}
