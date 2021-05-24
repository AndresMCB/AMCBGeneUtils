changeMatrixIds <- function(GEMatrix
                           , from=c("Ensembl.ID","HGNC.ID","HGNC.symbol","NCBI.ID")
                           , to=c("HGNC.symbol","NCBI.ID","Ensembl.ID","HGNC.ID")){
  # change colnames of GEMatrix "from" one format "to" another

  GEMatrix <- as.matrix(GEMatrix)

  ##  default case, change colnames from Ensembl.ID to HGNC.symbol
  if(length(from)>1)
    from<-from[1]
  if(length(to)>1)
    to<-to[1]

  genes <- changeGeneId(colnames(GEMatrix),from = from, to = to)
  index <- which(is.na(genes[[2]]))

  if(length(index)>0){
    message(paste0("------ ",toString(length(index))
                   ," genes cannot be mapped ------"))
    GEMatrix <- GEMatrix[,-index,drop=F]
    colnames(GEMatrix) <- genes[[2]][-index]
  }
  else{
    GEMatrix <- GEMatrix[,,drop=F]
    colnames(GEMatrix) <- genes[[2]]
  }

  return(GEMatrix)

}
