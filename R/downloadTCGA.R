downloadTCGA <- function(project="TCGA-BRCA"
                         , vital_status=c("Alive","Dead")
                         , shortLetterCode = c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC"
                                               ,"TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB"
                                               ,"CELL","XP","XCL")
                         , IncludeMut = TRUE
                         , clinicalSupp=TRUE
                         , normCounts = TRUE){

  code2type <- function(shortLetterCode){
    sample.type <- character(0)
    for (i in 1:length(shortLetterCode)) {
      switch(shortLetterCode[i]
             ,"TP"={sample.type<-c(sample.type,"Primary Tumor")}
             ,"TR"={sample.type<-c(sample.type,"Recurrent Tumor")}
             ,"TB"={sample.type<-c(sample.type,"Primary Blood Derived Cancer - Peripheral Blood")}
             ,"TRBM"={sample.type<-c(sample.type,"Recurrent Blood Derived Cancer - Bone Marrow")}
             ,"TAP"={sample.type<-c(sample.type,"Additional - New Primary")}
             ,"TM"={sample.type<-c(sample.type,"Metastatic")}
             ,"TAM"={sample.type<-c(sample.type,"Additional Metastatic")}
             ,"THOC"={sample.type<-c(sample.type,"Human Tumor Original Cells")}
             ,"TBM"={sample.type<-c(sample.type,"Primary Blood Derived Cancer - Bone Marrow")}
             ,"NB"={sample.type<-c(sample.type,"Blood Derived Normal")}
             ,"NT"={sample.type<-c(sample.type,"Solid Tissue Normal")}
             ,"NBC"={sample.type<-c(sample.type,"Buccal Cell Normal")}
             ,"NEBV"={sample.type<-c(sample.type,"EBV Immortalized Normal")}
             ,"NBM"={sample.type<-c(sample.type,"Bone Marrow Normal")}
             ,"CELLC"={sample.type<-c(sample.type,"Control Analyte")}
             ,"TRB"={sample.type<-c(sample.type,"Recurrent Blood Derived Cancer - Peripheral Blood")}
             ,"CELL"={sample.type<-c(sample.type,"Cell Lines")}
             ,"XP"={sample.type<-c(sample.type,"Primary Xenograft Tissue")}
             ,"XCL"={sample.type<-c(sample.type,"Cell Line Derived Xenograft Tissue")}
             ,{message(paste(shortLetterCode[i],
                             "code NOT FOUND.")
             )})
    }#end for
    return(sample.type)
  }#end code2type

  #####  ------- Main Code -------  #####

  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
    BiocManager::install("TCGAbiolinks")
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    BiocManager::install("SummarizedExperiment")
  if (!requireNamespace("EDASeq", quietly = TRUE))
    BiocManager::install("EDASeq")

  if (!requireNamespace("DT"))
    install.packages('DT')
  if (!requireNamespace("dplyr"))
    install.packages('dplyr')

  library(TCGAbiolinks)
  library(dplyr)
  library(DT)
  library(SummarizedExperiment)
  library(EDASeq)

  project <- toupper(project)
  vital_status <- tolower(vital_status)
  vital_status <- paste(toupper(substr(vital_status, 1, 1))
                          , substr(vital_status, 2, nchar(vital_status)), sep="")

  shortLetterCode<-toupper(shortLetterCode)

  #####   #####
  ##
  ## Legacy = TRUE:access to an unmodified copy of data that uses as
  ##             references GRCh37 (hg19) and GRCh36 (hg18).
  ## Legacy = FALSE: access data harmonized against GRCh38 (hg38), which provides
  ##             methods to the standardization of biospecimen and clinical data.

  sample.type = code2type(shortLetterCode)
  if(length(project)>1){
    message("Please retrieve one project at a time. Retrieving the first project.")
    project<-project[1]
  }

  GDCquery <- GDCquery(project = project
                         ,data.category = "Transcriptome Profiling"
                         ,data.type = "Gene Expression Quantification"
                         ,workflow.type = "HTSeq - Counts"
                         ,sample.type = sample.type
                         ,legacy = FALSE)


  ## Download data (saved in GDCdata folder by default)
  GDCdownload(GDCquery, method = "api", files.per.chunk = 10)
  data <- GDCprepare(GDCquery)

  ## keep data with selected vital status
  index <- which(data$vital_status==vital_status)
  data <- data[,index]
  res <- list()

  res$clinical<-as.data.frame(colData(data))

  ## download mutational data (MUSE pipeline)
  ## substring(project,6) removes "TCGA-" from project

  if(IncludeMut)
    res$MAF.hg38 <- GDCquery_Maf(substring(project,6), pipelines = "muse")

  if(clinicalSupp){
    query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Clinical",
                      data.type = "Clinical Supplement",
                      data.format = "BCR Biotab")
    GDCdownload(query)
    res$clinicalSup <- GDCprepare(query)
  }

  res$GenomicData <- as.data.frame(rowRanges(data))
  Counts <- assay(data)
  res$HTSeq_Counts <-as.data.frame(t(Counts))
  if(normCounts){
    dataNorm <- TCGAanalyze_Normalization(t(HTSeq_Counts)
                                          , geneInfoHT
                                          , method = "gcContent")
    res$HTSeq_Norm_Counts <- as.data.frame(t(dataNorm))
  }

  return(res)
}
