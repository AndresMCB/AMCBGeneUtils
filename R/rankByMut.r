rankByMut <- function(genesIds, project = "BRCA"){

  switch(project,
    BRCA={
      load(system.file("extdata"
                       , "TCGA_BRCA_TP.rda"
                       , package = "AMCBGeneUtils"
                       , mustWork = TRUE))
      rank <- sapply(genesIds
                     ,function(x,mutations){
                       res <- sum(mutations%in%x)
                     }
                     , TCGA_BRCA_TP$MAF.hg38$Gene)
    },
    GBM={

    }
  )

  rank <- as.data.frame(rank)
  row.names(rank) <- genesIds
  return(rank)
}
