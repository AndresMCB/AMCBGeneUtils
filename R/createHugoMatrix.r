
createHugoMatrix <- function(HGNCcustom=NULL, file = NULL){
  if (!requireNamespace("tidyverse"))
    install.packages('tidyverse')
  library(tidyverse)

  if(is.null(HGNCcustom)){
    HGNCdir <- system.file("extdata", "HGNCcustom.txt"
                       , package = "AMCBGeneUtils", mustWork = TRUE)
    HGNC <- read.delim(HGNCdir)
  }

  HGNC <- as_tibble(HGNC)
  #removing any possible white space at beginning/end of each data
  HGNC <- HGNC%>%mutate_all(str_trim)

  if(!is.null(file)){
    save(HGNC, file = file)
  }else{
    file <- paste0(getwd(),"/data/HGNC.rda")
    save(HGNC, file = file)
  }
  message("HGNCcustom.rda file created in:")
  message(file)

}


