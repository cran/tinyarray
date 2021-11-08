##' geo_download
##'
##' download gse data and get informations
##'
##' @param gse gse assession number
##' @param by_annopbrobe getGEO or geoChina
##' @param simpd get simplified pdata,drop out columns with all same values
##' @param colon_remove whether to remove duplicated columns with colons
##' @param destdir	 The destination directory for data downloads.
##' @return a list with exp,pd and gpl
##' @author Xiaojie Sun
##' @importFrom dplyr arrange
##' @importFrom dplyr filter
##' @importFrom dplyr %>%
##' @export
##' @examples
##' \donttest{
##' gse = "GSE42872"
##' a = geo_download(gse,by_annopbrobe = FALSE,destdir=tempdir())
##' }
##' @seealso
##' \code{\link{find_anno}}

geo_download <-  function(gse,by_annopbrobe = TRUE,
                          simpd = TRUE,colon_remove = FALSE,
                          destdir = getwd()){
  if(!requireNamespace("Biobase",quietly = TRUE)) {
    stop("you must install Biobase first by BiocManger::install('Biobase')",call. = FALSE)
  }
  if((!by_annopbrobe) & !requireNamespace("GEOquery",quietly = TRUE)) {
    stop("you must install GEOquery first by BiocManger::install('GEOquery')",call. = FALSE)
  }
  if((by_annopbrobe) & !requireNamespace("AnnoProbe",quietly = TRUE)) {
    stop("you must install AnnoProbe first by install.packages('AnnoProbe')",call. = FALSE)
  }
  if(by_annopbrobe){
    if(!file.exists(paste0(destdir,"/",gse,"_eSet.Rdata"))){
      eSet <- AnnoProbe::geoChina(gse, destdir = destdir)
    }else{
      suppressWarnings(load(paste0(destdir,"/",gse,"_eSet.Rdata")))
      eSet = gset
      rm(gset)
    }

  }else{
    eSet <- GEOquery::getGEO(gse,
                   destdir = destdir,
                   getGPL = FALSE)
  }
  exp <- Biobase::exprs(eSet[[1]])
  pd <- Biobase::pData(eSet[[1]])
  if(simpd){
    colname <- vector("character")
    count <- vector("integer")
    for (i in 1:ncol(pd)) {
      colname[i] = colnames(pd)[[i]]
      count[i] = nrow(pd[!duplicated(pd[, i]), ])
    }
    df <- data.frame(colname, count,stringsAsFactors = FALSE) %>% arrange(desc(count)) %>% dplyr::filter(count >1)
    pd <-  pd[,df$colname]
    pd <- pd[,!(colnames(pd) %in% c("geo_accession", "supplementary_file"))]
    if(colon_remove){
      pd = pd[,!apply(pd, 2, function(x){all(str_detect(x,": |https|www")|is.na(x)|x=="")})]
      colnames(pd) = str_remove(colnames(pd),":ch1")
    }
  }
  p1 = identical(rownames(pd),colnames(exp))
  p2 = all(rownames(pd) %in% colnames(exp) & colnames(exp) %in% rownames(pd))
  if(!p1) {
    exp = exp[,match(rownames(pd),colnames(exp))]
    if(!p2) {
      exp = exp[,intersect(rownames(pd),colnames(exp))]
      pd = pd[intersect(rownames(pd),colnames(exp)),]
    }
  }
  gpl <- eSet[[1]]@annotation
  re = list(exp=exp,pd=pd,gpl=gpl)
  if(is.null(dim(exp)) | nrow(exp)==0){
    warning("exp is empty")
  } else if (any(is.na(exp)|is.nan(exp))) {
    warning("NA or NAN values detected")
  }else if (any(exp<0)) {
    warning("nagtive values detected")
  } else{
    message(paste(nrow(exp),"probes,",
                  ncol(exp),"samples",
                  "from",min(exp),
                  "to",max(exp)))}
  return(re)
}

##' find annotation package or files
##'
##' find gpl annotation package or files
##'
##' @param gpl a gpl accession
##' @param install whether to install and library the package
##' @param update whether to update the package
##' @return a list with deg data.frame, volcano plot and a list with DEGs.
##' @author Xiaojie Sun
##' @importFrom stringr str_remove_all
##' @importFrom stringr str_to_upper
##' @importFrom BiocManager install
##' @export
##' @examples
##' find_anno("GPL570")
##' @seealso
##' \code{\link{geo_download}}

find_anno <-function(gpl,install = FALSE,update = FALSE){
  gpl = str_to_upper(gpl)
  if(!any(pkg_all$gpl==gpl)) {
    # R包不可用
    if(gpl %in% setdiff(exists_anno_list,pkg_all$gpl)){
      # 只有idmap可用
      ml1 = str_remove_all(paste0("`ids <- AnnoProbe::idmap\\(","\\'",gpl,"\\'","\\)`"),"\\\\")
      print(paste0("no annotation packages avliable,please use ",ml1))
    }else{
      # R包和idmap都不可用
      print("no annotation avliable in Bioconductor and AnnoProbe")
    }
  }else {
    qz = pkg_all$bioc_package[pkg_all$gpl== gpl]
    ml1 = str_remove_all(paste0("`ids <- AnnoProbe::idmap\\(","\\'",gpl,"\\'","\\)`"),"\\\\")
    ml2 = str_remove_all(paste0("`library\\(",qz,".db","\\)",";","ids <- toTable\\(",qz,"SYMBOL\\)`"),"\\\\")
    if(install){
      if(!suppressMessages(requireNamespace(paste0(qz,".db")))){
        BiocManager::install(paste0(qz,".db"),update = update)
        suppressMessages(requireNamespace(paste0(qz,".db")))
      }
    }
    if(!(gpl %in% exists_anno_list)) {
      #仅有R包可用
      print(paste0(ml2," is avaliable"))
    }else {
      #idmap和R包都可用
      print(paste0(ml2," and ",ml1 ," are both avaliable"))
    }
  }
}
utils::globalVariables(c("pkg_all","exists_anno_list","gset"))