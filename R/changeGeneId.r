changeGeneId <- function(IDs, from = "HGNC.symbol"
                          , to=c("Ensembl.ID","HGNC.ID","HGNC.symbol","NCBI.ID")){
  if (!requireNamespace("tidyverse"))
    install.packages('tidyverse')
  library(tidyverse)
  HGNC <- AMCBGeneUtils::HGNC

  IDs <- str_trim(IDs)
  IDs.in <- IDs
  res <- matrix(nrow = length(IDs), ncol = length(to)+1)

  # selecting proper name in file from https://www.genenames.org/download/custom/
  case <- which(c("Ensembl.ID","HGNC.ID","HGNC.symbol","NCBI.ID")%in%from)
  if (length(case)<1)
    stop("incorrect \"from\" label, correct labels are \"HGNC.ID\",\"HGNC.symbol\",\"NCBI.ID\" ")

  if(length(case)>1)
    warning("more than 1 \from\" label detected, taking the 1st one")

  if(case==3){ #if input is HGNC.symbol, we first update those them
    IDs.in <- IDs
    IDs <- updateGeneSymbol(IDs = IDs)[,2]
  }


  input <- switch(case[1],"Ensembl.ID.supplied.by.Ensembl."
                 ,"HGNC.ID"
                 ,"Approved.symbol"
                 ,"NCBI.Gene.ID.supplied.by.NCBI.")

  case <- which(c("Ensembl.ID","HGNC.ID","HGNC.symbol","NCBI.ID")%in%to)
  if (length(case)<1)
    stop("incorrect \"to\" label, correct labels are \"HGNC.ID\",\"HGNC.symbol\",\"NCBI.ID\" ")

  output <- sapply(case
                   ,function(a){switch(a,"Ensembl.ID.supplied.by.Ensembl."
                                            ,"HGNC.ID"
                                            ,"Approved.symbol"
                                            ,"NCBI.Gene.ID.supplied.by.NCBI.")
                     })

  index <- match(IDs,HGNC[[input]])
  res[,1] <- IDs.in
  res[,2:ncol(res)] <- as.matrix(HGNC[index,output])
  colnames(res) <- c("IDs",to)

  res<-as_tibble(res)
  res[res==""] <- NA
  return(res)

}
