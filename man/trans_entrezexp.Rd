% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_tcga_group_trans_exp.R
\name{trans_entrezexp}
\alias{trans_entrezexp}
\title{trans_entrezexp}
\usage{
trans_entrezexp(entrezexp, species = "human")
}
\arguments{
\item{entrezexp}{expression set with entrezid as rownames}

\item{species}{choose human or mouse, or rat, default: human}
}
\value{
a transformed expression set with symbol
}
\description{
transform rownames of expression set from "entrez" to"symbol",according to the bitr function.
}
\examples{
exp = matrix(rnorm(200),ncol = 10)
rownames(exp) = c("79691", "56271", "8662", "10394", "55630", "159162", "23541",
                  "79723", "54413", "22927", "92342", "23787", "5550", "8924",
                  "55274", "866", "8844", "353299", "587", "1473")
colnames(exp) = paste0("s",1:10)
if(requireNamespace("org.Hs.eg.db",quietly = TRUE)){
exp2 = trans_entrezexp(exp)
}else{
    warning("Package \"org.Hs.eg.db\" needed for this function to work.
        Please install it by BiocManager::install('org.Hs.eg.db')",call. = FALSE)
}
}
\seealso{
\code{\link{trans_exp}}
}
\author{
Xiaojie Sun
}
