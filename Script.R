#calling the libraries
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)

#using GDCQuery for creating a query
GDCProjects <- getGDCprojects()
getProjectSummary('TCGA-LUAD')
getProjectSummary('TCGA-LUSC')


#query using GDCquery function
query_LUAD <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling")
queryoutput_luad <- getResults(query_BRCA)

query_LUSC <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling")
queryoutput_lusc <- getResults(query_LUSC)

listsamples<- c ("TCGA-78-7166-01A-12R-2066-07","TCGA-38-4630-01A-01R-1206-07",
                 "TCGA-55-1592-01A-01R-0946-07","TCGA-73-4670-01A-01R-1206-07",
                 "TCGA-44-7661-01A-11R-2066-07", "TCGA-50-5932-11A-01R-1755-07",
                 "TCGA-49-6742-11A-01R-1858-07", "TCGA-44-6147-11A-01R-1858-07",
                 "TCGA-55-6979-11A-01R-1949-07", "TCGA-50-5931-11A-01R-1858-07",
                 "TCGA-37-A5EN-01A-21R-A26W-07", "TCGA-85-A4QR-01A-11R-A24Z-07",
                 "TCGA-21-1071-01A-01R-0692-07", "TCGA-37-3783-01A-01R-1201-07",
                 "TCGA-66-2765-01A-01R-0851-07", "TCGA-56-7582-11A-01R-2045-07",
                 "TCGA-22-5482-11A-01R-1635-07", "TCGA-58-8386-11A-01R-2296-07",
                 "TCGA-56-7222-11A-01R-2045-07", "TCGA-77-7142-11A-01R-2045-07"
)

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(
  project = c("TCGA-LUAD","TCGA-LUSC"), 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  barcode = listsamples
)

GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
NSCLC.Rnaseq.SE <- GDCprepare(query)
NSCLCMatrix <- assay (NSCLC.Rnaseq.SE, "unstranded")

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
NSCLC.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(NSCLC.Rnaseq.SE)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(
  tabDF = NSCLC.RNAseq_CorOutliers, 
  geneInfo =  geneInfoHT
)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)


# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)

# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(Diff_Expressed_ECM_Genes_Metadata)
library(EnsDb.Hsapiens.v86)

geneID_NSCLC <- ensembldb::select(EnsDb.Hsapiens.v86, 
                                 keys= Genelist, 
                                 keytype = "GENEID", 
                                 columns = c("SYMBOL","GENEID"))

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  Genelist)

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP), 
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = Genelist, 
  nBar = 25
)

# selection of normal samples "NT" 
group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
# selection of normal samples "TP" 
group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

# Principal Component Analysis plot for ntop selected DEGs
pca <- TCGAvisualize_PCA(
  dataFilt = dataFilt,
  dataDEGsFiltLevel = DEGsFiltLevel,
  ntopgenes = 200, 
  group1 = group1,
  group2 =  group2
)
Diff_Expressed_ECM_Genes_Metadata <- Diff_Expressed_ECM_Genes_Metadata %>% remove_rownames %>% column_to_rownames(var="Gene Name")
x<- Diff_Expressed_ECM_Genes_Metadata$logFC
y<- Diff_Expressed_ECM_Genes_Metadata$FDR

TCGAVisualize_volcano(
  x,
  y,
  filename = "volcano_ECM.pdf",
  ylab = expression(paste(-Log[10], " (FDR corrected -P values)")),
  xlab = NULL,
  names = rownames (Diff_Expressed_ECM_Genes_Metadata),
  title = "Volcano plot",
  legend = NULL,
  label = NULL,
  xlim = NULL,
  ylim = NULL,
  color = c("black", "red", "green"),
  names.fill = TRUE,
  show.names = "significant",
  x.cut = 0,
  y.cut = 10^-5,
  height = 5,
  width = 10,
  highlight = NULL,
  highlight.color = "orange",
  names.size = 4,
  dpi = 300
)

#pathway analysis
install.packages("pathfindR")
