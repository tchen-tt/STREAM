#' @title checkcellphonedb
#' @param cellphonedbPath
#' 
checkcellphonedb<-function(cellphonedbPath=NULL,...){
 if(is.null(cellphonedbPath)){
   checkCellphonedb <- system("which cellphonedb", intern = TRUE)
   if (length(checkCellphonedb) == 0) {
     system("pip3 install cellphonedb")
     #pyPath <- system("which python3", intern = TRUE)
     #stop(paste0("No Cellphonedb in",pyPath,"\n Install cellphonedb first. \n More details in https://github.com/Teichlab/cellphonedb"))
   }
   cellphonedbPath=system("which cellphonedb", intern = TRUE)
 }
  if(!file.exists("/Users/xiaoyihan/miniconda3/bin/cellphonedb")){
    stop("No cellphonedb found in the provided path.")
  }
  return(cellphonedbPath)
}

