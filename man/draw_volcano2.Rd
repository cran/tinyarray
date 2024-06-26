% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6_multi_deg_all.R
\name{draw_volcano2}
\alias{draw_volcano2}
\title{draw_volcano2}
\usage{
draw_volcano2(deg, pkg = 4, lab, ...)
}
\arguments{
\item{deg}{a data.frame created by Differential analysis}

\item{pkg}{a integer ,means which Differential analysis packages you used,we support three packages by now, 1,2,3,4 respectively means "DESeq2","edgeR","limma(voom)","limma"}

\item{lab}{label for  x axis in volcano plot, if NA , x axis names by package}

\item{...}{other parameters from draw_volcano}
}
\value{
one or more volcano plot
}
\description{
print one or more volcano plot for Differential analysis result in data.frame fomat.
}
\examples{
\dontrun{
if(requireNamespace("Biobase",quietly = TRUE)&
   requireNamespace("AnnoProbe",quietly = TRUE)){
#two group
gse = "GSE42872"
geo = geo_download(gse,destdir=tempdir())
group_list = rep(c("A","B"),each = 3)
ids = AnnoProbe::idmap('GPL6244',destdir = tempdir())
deg = get_deg(geo$exp,group_list,ids)
draw_volcano2(deg)
#multigroup
gse = "GSE474"
geo = geo_download(gse,destdir=tempdir())
geo$exp[1:4,1:4]
geo$exp=log2(geo$exp+1)
group_list=ifelse(stringr::str_detect(geo$pd$title,"MObese"),"MObese",
ifelse(stringr::str_detect(geo$pd$title,"NonObese"),"NonObese","Obese"))
group_list=factor(group_list,levels = c("NonObese","Obese","MObese"))
find_anno(geo$gpl)
ids <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
deg = multi_deg(geo$exp,group_list,ids,adjust = FALSE,entriz = FALSE)
draw_volcano2(deg)
draw_volcano2(deg,color = c("darkgreen","grey","darkred"))
}else{
  if(!requireNamespace("AnnoProbe",quietly = TRUE)) {
    warning("Package 'AnnoProbe' needed for this function to work.
         Please install it by install.packages('AnnoProbe')",call. = FALSE)
  }
  if(!requireNamespace("Biobase",quietly = TRUE)) {
    warning("Package 'Biobase' needed for this function to work.
         Please install it by BiocManager::install('Biobase')",call. = FALSE)
  }
}
}
}
\seealso{
\code{\link{geo_download}};\code{\link{draw_volcano}};\code{\link{draw_venn}}
}
\author{
Xiaojie Sun
}
