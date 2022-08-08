rankByMut <- function(genesIds, project = "BRCA", patient_IDs=NULL){

  switch(project,
    BRCA={
    x = load(system.file("extdata"
                       , "BRCA.MAF.hg38.rda"
                       , package = "AMCBGeneUtils"
                       , mustWork = TRUE))

    MAF.hg38 = get(x)
    rm(x)

    },
    GBM={

    }
  )

  if (!is.null(patient_IDs)){
    MAF.hg38 <- MAF.hg38%>%
      mutate(patient=str_trunc(Tumor_Sample_Barcode,12, "right", ellipsis = ""),.before = 1)%>%
      dplyr::filter(patient%in%patient_IDs)
  }

  aux <- MAF.hg38%>%
    dplyr::filter(Gene%in%genesIds)%>%
    dplyr::count(Ensembl.ID=Gene, name = "rank")

  output <- data.frame("Ensembl.ID"=genesIds)
  output <- left_join(x=output,y=aux,by="Ensembl.ID")
  output[is.na(output$rank),"rank"] <- 0
  output <- output%>%
    arrange(desc(rank))

  row.names(output) <- output$Ensembl.ID
  return(output)
}
