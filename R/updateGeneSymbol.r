updateGeneSymbol <- function(IDs){
  if (!requireNamespace("tidyverse"))
    install.packages('tidyverse')
  library(tidyverse)
  HGNC <- AMCBGeneUtils::HGNC

  IDs <- str_trim(IDs)
  res <- matrix(nrow = length(IDs), ncol = 2)
  colnames(res) <- c("IDs","HGNC.symbol")
  res[,1] <- IDs
  # tidying up string data
  aux <- HGNC%>%transmute(Approved.symbol
                      , allNames = str_c(Approved.symbol,Previous.symbols
                                           ,Alias.symbols,sep = ","))%>%
    mutate(allNames=str_replace_all(allNames," ",""))%>%
    mutate(allNames=str_replace_all(allNames,",+"," "))%>%
    mutate(allNames=str_trim(allNames))


  b <- apply(aux, MARGIN = 1
             , function(x){
               temp <- t(str_split(x[2], " ", simplify = T))
               temp <- cbind(rep(temp[1,1],ncol(temp)),temp)
             })
  b <- do.call(what = rbind,args = b)
  colnames(b) <- c("HGNC.symbol","symbol")

  index <- match(res[,1],as.matrix(b[,2]))
  res[,2] <- b[index,1]
  return(res)

}
