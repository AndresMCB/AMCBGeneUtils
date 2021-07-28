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
  previous[,2] <-toupper(previous[,2])

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
  Alias[,2] <- toupper(Alias[,2])
  Alias <- as_tibble(Alias)%>%
    filter(!IDs%in%previous[,"IDs"])

  previous <- rbind(previous,Alias)

  rm(Alias, aux)
  gc()

  # Find those IDs already in HGNC$Approved.symbol
  res[,1] <- IDs
  index <- tolower(res[,1])%in%tolower(HGNC.Approved$Approved.symbol)
  res[index,2] <- res[index,1]

  # Check if IDs are in previous
  temp <- as_tibble(res)%>%
    filter(is.na(HGNC.symbol))%>%
    mutate(IDs = toupper(IDs))

  aux <- left_join(temp,previous, by = "IDs")

  if(nrow(aux) > nrow(temp)){
    warning(paste0("The following genes symbols are duplicated \n"
                   ,"and/or maps to a more than one updated symbol.\n"
                   ,na.omit(aux[duplicated(aux[,1]),1])
                   ,"\nThis can leads to an incorrect updated symbol.\n"
                   ,"Please check your input gene IDs"))
  }

  aux <- aux %>% distinct(IDs, .keep_all = T)%>%
    transmute(IDs,HGNC.symbol = HGNC.symbol.y)

  temp <- left_join(temp,aux, by ="IDs")%>%
    transmute(IDs,HGNC.symbol = HGNC.symbol.y)

  res[!index,2] <- temp$HGNC.symbol



  #index <- match(res[,1],as.matrix(b[,2]))
  #res[,2] <- b[index,1]
  return(res)

}
