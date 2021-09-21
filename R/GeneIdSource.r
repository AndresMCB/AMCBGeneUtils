GeneIdSource <- function(geneIDs){
  if(!require(stringr))
    install.packages("stringr")

  library(stringr)

  # make geneIDs data structure a char array
  geneIDs <- sapply(geneIDs, as.character, simplify = T)

  n <- length(geneIDs)
  pattern <- regex("ENSG0", ignore_case = T)
  aux <- sum(str_detect(geneIDs, pattern = pattern))
  if(aux > n/2){
    Source <- "Ensembl.ID"
  }else
  {
    pattern <- regex("HGNC:", ignore_case = T)
    aux <- sum(str_detect(geneIDs, pattern = pattern))
    if(aux  > n/2){
      Source <- "HGNC.ID"
    }else
    {
      pattern <- "[:alpha:]"
      aux <- sum(str_detect(geneIDs, pattern = pattern))
      if(aux < n/2){
        Source <- "NCBI.ID"
      }else
      {
        message("Gene IDs do not coorespond to Ensembl.ID, HGNC.ID nor NCBI.ID")
        message("Assuming HGNC.symbol")
        Source <- "HGNC.symbol"
      }
    }
  }
  return(Source)
}





