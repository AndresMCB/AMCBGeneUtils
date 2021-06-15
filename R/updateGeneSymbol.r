updateGeneSymbol <- function(IDs){
  if (!requireNamespace("tidyverse"))
    install.packages('tidyverse')
  library(tidyverse)
  HGNC.Approved <- AMCBGeneUtils::HGNC%>%
    filter(Status=="Approved")

  IDs <- str_trim(IDs)
  res <- matrix(nrow = length(IDs), ncol = 2)
  colnames(res) <- c("IDs","HGNC.symbol")

  # Tidying up string data for previous symbols
  previous <- HGNC.Approved%>%
    transmute(Approved.symbol
                      , Previous.symbols = str_c(Previous.symbols
                                           ,sep = ","))%>%
    mutate(Previous.symbols=str_replace_all(Previous.symbols," ",""))%>%
    mutate(Previous.symbols=str_replace_all(Previous.symbols,",+"," "))%>%
    mutate(Previous.symbols=str_trim(Previous.symbols))%>%
    filter(Previous.symbols!="")

  aux <- apply(previous, MARGIN = 1
             , function(x){
               temp <- t(str_split(x[2], " ", simplify = T))
               temp <- cbind(rep(x[1, drop=T],ncol(temp)),temp)
             })

  previous <- do.call(what = rbind,args = aux)
  colnames(previous) <- c("HGNC.symbol","IDs")

  # Tidying up string data for Alias
  Alias <- HGNC.Approved%>%
    transmute(Approved.symbol
              , Alias.symbols = str_c(Alias.symbols
                                         ,sep = ","))%>%
    mutate(Alias.symbols=str_replace_all(Alias.symbols," ",""))%>%
    mutate(Alias.symbols=str_replace_all(Alias.symbols,",+"," "))%>%
    mutate(Alias.symbols=str_trim(Alias.symbols))%>%
    filter(Alias.symbols!="")

  aux <- apply(Alias, MARGIN = 1
                , function(x){
                  temp <- t(str_split(x[2], " ", simplify = T))
                  temp <- cbind(rep(x[1, drop=T],ncol(temp)),temp)
                })

  Alias <- do.call(what = rbind,args = aux)
  colnames(Alias) <- c("HGNC.symbol","IDs")


  # Find those IDs already in HGNC$Approved.symbol
  res[,1] <- IDs
  index <- res[,1]%in%HGNC.Approved$Approved.symbol
  res[index,2] <- res[index,1]

  # Check if IDs are in previous
  temp <- left_join(as_tibble(res[!index,1, drop=F])
                    ,as_tibble(previous), by = "IDs")

  if(sum(duplicated(temp$IDs))>0){
    message("There are duplicated gene symbols. Keeping the first ocurrence.")
    message("duplicated genes:")
    print(temp$IDs[duplicated(temp$IDs)])
    temp <- temp[!duplicated(temp$IDs),]
  }
  res[!index,2] <- temp$HGNC.symbol

  # Check if IDs are in Alias
  index <- is.na(res[,2])
  if(sum(index)>0){
    temp <- left_join(as_tibble(res[index,1, drop=F])
                      ,as_tibble(Alias), by = "IDs")

    if(sum(duplicated(temp$IDs))>0){
      message("There are duplicated gene symbols. Keeping the first ocurrence.")
      message("duplicated genes:")
      print(temp$IDs[duplicated(temp$IDs)])
      temp <- temp[!duplicated(temp$IDs),]
    }
    res[index,2] <- temp$HGNC.symbol
  }



  aux <- res[!index,2]
  aux %in% previous[,2]



  #index <- match(res[,1],as.matrix(b[,2]))
  #res[,2] <- b[index,1]
  return(res)

}
