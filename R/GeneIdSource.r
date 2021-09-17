GeneIdSource <- function(geneIDs){
  if(!require(stringr))
    install.packages("stringr")

  library(stringr)

  # make geneIDs data structure a char array
  geneIDs <- sapply(geneIDs, as.character, simplify = T)

  n <- length(geneIDs)
  pattern <- regex("ENSG0", ignore_case = T)
  aux <- sum(str_detect(geneIDs, pattern = pattern))
  if(aux == n){
    Source <- "Ensembl.ID"
  }else
  {
    pattern <- regex("HGNC:", ignore_case = T)
    aux <- sum(str_detect(geneIDs, pattern = pattern))
    if(aux == n){
      Source <- "HGNC.ID"
    }else
    {
      pattern <- "[:alpha:]"
      aux <- sum(str_detect(geneIDs, pattern = pattern))
      if(aux == 0){
        Source <- "NCBI.ID"
      }else
      {
        Source <- "HGNC.symbol"
      }
    }
  }
  return(Source)
}





