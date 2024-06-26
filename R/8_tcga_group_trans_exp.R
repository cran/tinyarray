##' make_tcga_group
##'
##' make tcga group for given tcga expression matrix
##'
##' @inheritParams trans_exp
##' @importFrom stringr str_starts
##' @importFrom stringr str_sub
##' @export
##' @return a group factor with normal and tumor ,correspond to colnames for expression matrix
##' @author Xiaojie Sun
##' @examples
##' k = make_tcga_group(exp_hub1);table(k)
##' @seealso
##' \code{\link{sam_filter}};\code{\link{match_exp_cl}}

make_tcga_group <- function(exp){
  k1 = stringr::str_starts(colnames(exp),"TCGA")
  if(!any(k1))stop("no TCGA samples detected,please check it")
  k2 = suppressWarnings(as.numeric(stringr::str_sub(colnames(exp),14,15))<10)
  group_list = ifelse(k1&k2,"tumor","normal")
  group_list = factor(group_list,levels = c("normal","tumor"))
  return(group_list)
}

##' trans_exp
##'
##' transform rownames of TCGA or TCGA_Gtex expression set from gdc or xena,from ensembl id to gene symbol
##'
##' @param exp TCGA or TCGA_Gtex expression set from gdc or xena
##' @param mrna_only only keep mrna rows in result
##' @param lncrna_only only keep lncrna rows in result
##' @param gtex logical,whether including Gtex data
##' @return a transformed expression set with symbol
##' @author Xiaojie Sun
##' @importFrom stringr str_detect
##' @importFrom stringr str_remove
##' @importFrom dplyr inner_join
##' @export
##' @examples
##' exp = matrix(rnorm(1000),ncol = 10)
##' rownames(exp) = sample(mRNA_annov23$gene_id,100)
##' colnames(exp) = c(paste0("TCGA",1:5),paste0("GTEX",1:5))
##' k  = trans_exp(exp)
##' @seealso
##' \code{\link{trans_array}}

trans_exp = function(exp,mrna_only = FALSE,
                     lncrna_only = FALSE,gtex = FALSE){
  k00 = any(str_detect(colnames(exp),"TCGA"))
  if(!k00)warning("this expression set probably not from TCGA,please ensure it")
  k0 = any(str_detect(colnames(exp),"GTEX"))
  kd = any(str_detect(rownames(exp),"\\."))
  if((!(k0|gtex))){
    lanno = lnc_anno
    manno = mRNA_anno
  }else if(k00){
    lanno = lnc_annov23
    manno = mRNA_annov23
  }
  if(!kd){
    lanno$gene_id = str_remove(lanno$gene_id,"\\.\\d*")
    manno$gene_id = str_remove(manno$gene_id,"\\.\\d*")
  }
  n1 = sum(rownames(exp) %in% manno$gene_id)
  k1 = length(n1)/nrow(exp)< 0.25 & length(n1)<5000
  n2 = sum(rownames(exp) %in% lanno$gene_id)
  k2 = length(n2)/nrow(exp)< 0.25 & length(n2)<5000
  mRNA_exp = exp[rownames(exp) %in% manno$gene_id,]
  tmp = data.frame(gene_id = rownames(exp))
  x = dplyr::inner_join(tmp,manno,by = "gene_id")
  mRNA_exp = mRNA_exp[!duplicated(x$gene_name),]
  x = x[!duplicated(x$gene_name),]
  rownames(mRNA_exp) = x$gene_name
  lnc_exp = exp[rownames(exp) %in% lanno$gene_id,]
  tmp = data.frame(gene_id = rownames(exp))
  x = dplyr::inner_join(tmp,lanno,by = "gene_id")
  lnc_exp = lnc_exp[!duplicated(x$gene_name),]
  x = x[!duplicated(x$gene_name),]
  rownames(lnc_exp) = x$gene_name
  message(paste0(nrow(mRNA_exp),
                 " of genes successfully mapping to mRNA,",
                 nrow(lnc_exp),
                 " of genes successfully mapping to lncRNA"))
  if(mrna_only){
    return(mRNA_exp)
  }else if(lncrna_only){
      return(lnc_exp)
  }else{
    expa  = rbind(mRNA_exp,lnc_exp)
    k = !duplicated(rownames(expa))
    expa = expa[k,]
    return(expa)
    message(paste0(sum(!k),
                   " of duplicaterd genes removed"))
    }
}
utils::globalVariables(c("lnc_anno","mRNA_anno","lnc_annov23","mRNA_annov23"))

##' trans_array
##'
##' transform rownames for microarray or rnaseq expression matrix
##'
##' @param exp  microarray expression matrix with probe_id as rownames
##' @param ids data.frame  with original rownames and new rownames
##' @param from colname for original rownames
##' @param to colname for new rownames
##' @return a transformed expression set with new rownames
##' @author Xiaojie Sun
##' @export
##' @examples
##' exp = matrix(1:50,nrow = 10)
##' rownames(exp) = paste0("g",1:10)
##' ids = data.frame(probe_id = paste0("g",1:10),
##'                 symbol = paste0("G",c(1:9,9)))
##' trans_array(exp,ids)
##' @seealso
##' \code{\link{trans_exp}}

trans_array = function(exp,ids,from = "probe_id",
                       to = "symbol"){
  if(!is.character(ids[,from])) ids[,from] = as.character(ids[,from])
  a = intersect(rownames(exp),ids[,from])
  ids = ids[!duplicated(ids[,to]),]
  exp = exp[rownames(exp) %in% ids[,from],]
  ids = ids[ids[,from]%in% rownames(exp),]
  exp = exp[ids[,from],]
  rownames(exp)=ids[,to]
  message(paste0(nrow(exp)," rownames transformed after duplicate rows removed"))
  return(exp)
}

##' sam_filter
##'
##' drop duplicated samples from the same patients
##'
##' @param exp TCGA or TCGA_Gtex expression set from gdc or xena
##' @return a transformed expression set without duplicated samples
##' @author Xiaojie Sun
##' @export
##' @examples
##' cod[1:4,1:4]
##' dim(cod)
##' cod2 = sam_filter(cod)
##' dim(cod2)
##' g = make_tcga_group(cod);table(g)
##' library(stringr)
##' table(!duplicated(str_sub(colnames(cod[,g=="tumor"]),1,12)))
##' @seealso
##' \code{\link{make_tcga_group}};\code{\link{match_exp_cl}}

sam_filter = function(exp){
  exp = exp[,order(colnames(exp))]
  n1 = ncol(exp)
  group = make_tcga_group(exp)
  exptumor = exp[,group == "tumor"]

  expnormol = exp[,group == "normal"]
  exptumor = exptumor[,!duplicated(str_sub(colnames(exptumor),1,12))]
  expnormol = expnormol[,!duplicated(str_sub(colnames(expnormol),1,12))]

  exp = cbind(exptumor,expnormol)
  message(paste("filtered",n1-ncol(exp),"samples."))
  return(exp)
}


##' match_exp_cl
##'
##' match exp and clinical data from TCGA
##'
##' @param exp TCGA  expression set
##' @param cl TCGA clinical data.frame
##' @param id_column which column contains patient ids, column number or colnmn name.
##' @param sample_centric logical,deault T,keep all samples from the same patients.if FALSE,keep only one tumor sample for one patient.
##' @return a transformed clinical data.frame with sample ids.
##' @author Xiaojie Sun
##' @export
##' @examples
##' a = match_exp_cl(exp_hub1,meta1[,2:4],"X_PATIENT")
##' exp_matched = a[[1]]
##' cl_matched = a[[2]]
##' b = match_exp_cl(exp_hub1,meta1[,2:4],"X_PATIENT",sample_centric = FALSE)
##' exp_matched = b[[1]]
##' cl_matched = b[[2]]
##' @seealso
##' \code{\link{make_tcga_group}};\code{\link{sam_filter}}

match_exp_cl = function(exp,cl,id_column = "id",sample_centric = TRUE){
  colnames(cl)[colnames(cl)==id_column] = "id"
  cl = cl[cl$id %in% substr(colnames(exp),1,12),]
  exp = exp[,substr(colnames(exp),1,12) %in% cl$id]
  patient <- substr(colnames(exp),1,12)
  if(nrow(cl)==0) stop("your exp or cl doesn't match,please check them.")
  da = data.frame(sample_id = colnames(exp),
                  id = patient)
  cl = merge(da,cl,by ="id",all.y = TRUE)
  cl = cl[match(colnames(exp),cl$sample_id),]
  if(!sample_centric) {
    Group = make_tcga_group(exp)
    exp = exp[,Group=="tumor"]
    cl = cl[Group=="tumor",]
    cl = cl[order(colnames(exp)),]
    exp = exp[,sort(colnames(exp))]
    exp = exp[,!duplicated(cl$id)]
    cl = cl[!duplicated(cl$id),]
  }
  k = identical(colnames(exp),cl$sample_id)
  if(k)message("match successfully")
  rownames(cl) = cl$sample_id
  compiler::setCompilerOptions(suppressAll = TRUE)
  return(list(exp_matched = exp,
              cl_matched = cl))
  message("New version of tinyarray canceled global assigning inside the package,
          please obtain exp_matched and cl_matched by split this list result.")
}

##' trans_exp_new
##'
##' transform rownames of expression set from "ensembl" to"symbol",according to the new information from ensembl database.
##'
##' @param exp expression set with ensembl as rownames
##' @param mrna_only only keep mrna rows in result
##' @param lncrna_only only keep lncrna rows in result
##' @param species choose human or mouse, or rat, default: human
##' @return a transformed expression set with symbol
##' @author Xiaojie Sun
##' @importFrom stringr str_split
##' @export
##' @examples
##' exp = matrix(rnorm(1000),ncol = 10)
##' rownames(exp) = sample(mRNA_annov23$gene_id,100)
##' colnames(exp) = c(paste0("TCGA",1:5),paste0("GTEX",1:5))
##' if(requireNamespace("AnnoProbe")){
##' k  = trans_exp_new(exp)
##' }else{
##'   warning("Package \"AnnoProbe\" needed for this function to work.
##'          Please install it by install.packages('AnnoProbe')")
##' }
##' @seealso
##' \code{\link{trans_exp}}
trans_exp_new = function(exp,mrna_only = FALSE,
                         lncrna_only = FALSE,
                         species = "human"){
  if(is.data.frame(exp)){exp = as.matrix(exp)}
  if(!requireNamespace("AnnoProbe"))stop("Package \"AnnoProbe\" needed for this function to work.
         Please install it by install.packages('AnnoProbe')",call. = FALSE)
  rownames(exp) = str_split(rownames(exp),"\\.",simplify = T)[,1]
  re = AnnoProbe::annoGene(rownames(exp),ID_type = "ENSEMBL",species = species)
  if(mrna_only){
    re = re[re$biotypes=="protein_coding",]
  }else if(lncrna_only){
    re = re[re$biotypes=="lncRNA",]
  }else{
    re = re
  }
  exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")
  return(exp)
}

##' trans_entrezexp
##'
##' transform rownames of expression set from "entrez" to"symbol",according to the bitr function.
##'
##' @param entrezexp expression set with entrezid as rownames
##' @param species choose human or mouse, or rat, default: human
##' @importFrom clusterProfiler bitr
##' @return a transformed expression set with symbol
##' @author Xiaojie Sun
##' @export
##' @examples
##' exp = matrix(rnorm(200),ncol = 10)
##' rownames(exp) = c("79691", "56271", "8662", "10394", "55630", "159162", "23541",
##'                   "79723", "54413", "22927", "92342", "23787", "5550", "8924",
##'                   "55274", "866", "8844", "353299", "587", "1473")
##' colnames(exp) = paste0("s",1:10)
##' if(requireNamespace("org.Hs.eg.db",quietly = TRUE)){
##' exp2 = trans_entrezexp(exp)
##' }else{
##'     warning("Package \"org.Hs.eg.db\" needed for this function to work.
##'         Please install it by BiocManager::install('org.Hs.eg.db')",call. = FALSE)
##' }
##' @seealso
##' \code{\link{trans_exp}}
trans_entrezexp = function(entrezexp,species="human"){
  if(species == "human"){
    if(!requireNamespace("org.Hs.eg.db",quietly = TRUE)) {
      stop("Package \"org.Hs.eg.db\" needed for this function to work.
         Please install it by BiocManager::install('org.Hs.eg.db')",call. = FALSE)
    }
    or = org.Hs.eg.db::org.Hs.eg.db
  }
  if(species == "mouse"){
    if(!requireNamespace("org.Mm.eg.db",quietly = TRUE)) {
      stop("Package \"org.Mm.eg.db\" needed for this function to work.
         Please install it by BiocManager::install('org.Mm.eg.db')",call. = FALSE)
    }
    or = org.Mm.eg.db::org.Mm.eg.db
  }
  if(species == "rat"){
    if(!requireNamespace("org.Rn.eg.db",quietly = TRUE)) {
      stop("Package \"org.Rn.eg.db\" needed for this function to work.
         Please install it by BiocManager::install('org.Rn.eg.db')",call. = FALSE)
    }
    or = org.Rn.eg.db::org.Rn.eg.db
  }
  if(is.data.frame(entrezexp)){exp = as.matrix(entrezexp)}
  re = bitr(rownames(entrezexp),fromType = "ENTREZID",toType = "SYMBOL",OrgDb = or)
  exp = trans_array(entrezexp,ids = re,from = "ENTREZID",to = "SYMBOL")
  return(exp)
}
