##' geo_parser
##'
##' download gse data and get informations
##'
##' @param gse gse assession number
##' @param destdir	 The destination directory for data downloads
##' @importFrom stringr str_remove_all
##' @importFrom stringr str_remove
##' @importFrom stringr str_starts
##' @importFrom stringr str_split
##' @importFrom stringr str_detect
##' @importFrom stringr str_extract
##' @return an ExpressionSet object list
##' @author Xiaojie Sun
##' @export
##' @examples
##' \dontrun{
##' if(requireNamespace("GEOquery",quietly = TRUE)){
##'   gse = "GSE42872"
##'   a = geo_download(gse,destdir=tempdir())
##' }else{
##'   print("Package 'GEOquery' needed for this function to work.
##'          Please install it by BiocManager::install('GEOquery')")
##' }
##' }
##' @seealso
##' \code{\link{find_anno}}

geo_parser <- function(gse,destdir = getwd()) {
  series_matrix_file = paste0(gse,'_series_matrix.txt.gz')
  if(!file.exists(series_matrix_file)){
    if(!requireNamespace("GEOquery",quietly = TRUE)) {
      stop("Package \"GEOquery\" needed for this function to work.
         Please install it by BiocManager::install('GEOquery')",call. = FALSE)
    }
    tryCatch({GEOquery::getGEO(gse, destdir = destdir,getGPL = FALSE)
    },error = function(e){555})
  }
  if(!file.exists(series_matrix_file)){
    stop('Unable to download the series matrix file. Please download the series
         matrix file from GEO, place it in your working directory, and then
         rerun this function.')
  }
  if(!requireNamespace("Biobase",quietly = TRUE)) {
    stop("Package \"Biobase\" needed for this function to work.
         Please install it by BiocManager::install('Biobase')",call. = FALSE)
  }else{
    lines <- readLines(series_matrix_file)
    pd <- str_remove_all(lines, "^!|\"") |>
      str_split("\t", simplify = TRUE)
    k <- pd[,1] |> str_starts("Sample_")
    pd <- pd[k, ]
    rownames(pd) <- pd[,1] |> str_remove("Sample_")
    pd <- t(pd[,-1])
    rownames(pd) = pd[,2]
    metaData <- data.frame(labelDescription=colnames(pd))
    pd <- Biobase::AnnotatedDataFrame(data = data.frame(pd),varMetadata = metaData)
    exp <- utils::read.delim(series_matrix_file,
                      comment.char = "!",
                      row.names = 1) |>
      as.matrix()
    gpl = str_extract(lines[which(str_detect(lines,'GPL'))[1]],"GPL\\d+")
    return(list(Biobase::ExpressionSet(assayData = exp,phenoData = pd,annotation = gpl)))
  }
}
