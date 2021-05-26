rankByMut <- function(genesIds, project = "BRCA"){

  switch(project,
    BRCA={
      load(system.file("extdata"
                       , "BRCA.MAF.hg38.rda"
                       , package = "AMCBGeneUtils"
                       , mustWork = TRUE))
      rank <- sapply(genesIds
                     ,function(x,mutations){
                       res <- sum(mutations%in%x)
                     }
                     , BRCA.MAF.hg38$Gene)
    },
    GBM={

    }
  )

  rank <- as.data.frame(rank)
  row.names(rank) <- genesIds
  return(rank)
}
