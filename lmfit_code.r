## try http:// if https:// URLs are not supported
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
library("limma")

#rm(sc_1HT_GEX_data_filter)
lmfit<-function(AE){
  design <-rbind(matrix(wt),as.matrix(lt))
  AE=apply(AE,c(1,2),as.numeric)
  fit = lmFit(AE,design)
  fit = eBayes(fit)
  topTable(fit)
  
  a=rownames(AE[which(fit$F.p.value<="0.05"),])
  b=fit$F.p.value[which(fit$F.p.value<="0.05")]
  c=AE[which(fit$F.p.value<="0.05"),]
  d=matrix(nrow=nrow(c),ncol=1)
  d=log2(rowMeans(c[,1:15])/rowMeans(c[,16:22]))
  selrow=which(abs(d)>=0.584)
  result=cbind(a,b,c,d)[selrow,]
  
  colnames(result)=c("geneid","adjust.p",rep("case",15),rep("control",8),"fold change")
  head(result)
  return(result)
}

library("DESeq2")
#combined_mtx
PERFORM_DESEQ2 = function(data_mtx,control,case,output_dir){
  print(case)
  ctrl_mtx = data_mtx[,colnames(data_mtx) == control]
  case_mtx = data_mtx[,colnames(data_mtx) == case]
  combined_mtx = cbind(ctrl_mtx,case_mtx)
  
  combined_mtx = apply(combined_mtx,c(1,2),as.numeric)
  
  coldata = data.frame(condition = c(rep(control,ncol(ctrl_mtx)),rep(case,ncol(case_mtx))))
  combined_mtx_int = data.frame(apply(combined_mtx, 2, function(x) as.integer(as.character(x))))
  dim(combined_mtx_int)
  dds <- DESeqDataSetFromMatrix(countData = combined_mtx_int,
                                colData = coldata,
                                design = ~ condition)
  #featureData <- data.frame(gene=rownames(combined_mtx))
  #mcols(dds) <- DataFrame(mcols(dds), featureData)
  
  dds <- DESeq(dds)
  res <- results(dds)
  rownames(res) = data_mtx$Name
  res['symbol'] = rownames(res)
  res_filtered = res[!is.na(res$pvalue) & res$pvalue<=0.05,]
  write.table(res_filtered,paste0(output_dir,case,"_sig_DE_Results.txt"),sep="\t",col.names = NA,quote=F)
}

#path='C:\\Users\\yzlco\\Desktop\\PAGER_COVID19\\casestudy\\input'
getwd()
path = '/home/zongliang/DEseq2/input/'
setwd(path)

#############################
### cbioPortal_melanoma_analysis link: ###

path = '/home/zongliang/PAGER/PAGER3/cBioPortal/'

meta_pat = read.csv(paste0(path,"data_clinical_patient.txt"),sep="\t")
meta_pat$Treatment.Response
dict_sampleID2res_status = as.matrix(meta_pat$Treatment.Response)
rownames(dict_sampleID2res_status)=meta_pat$X.Patient.Identifier

mtx = read.csv(paste0(path,"data_RNA_Seq_expression_median.txt"),sep="\t")
mtx
status = vector()
samples=colnames(mtx)[-1]

colnames(mtx)=c("Name",dict_sampleID2res_status[samples,])
mtx
mtx$Name[isNA(mtx$Name)]

### gene symbol retrieval ###
require('org.Hs.eg.db')
BiocManager::install("GeneAnswers")
entrezID2geneSymbol = GeneAnswers::getSymbols(as.character(mtx$Name), 'org.Hs.eg.db')#as.vector(mtx$Name)
# optional to do the average expression of the same gene symbol
duplicate_names = table(entrezID2geneSymbol)[table(entrezID2geneSymbol)>1]
duplicate_names

mtx$Name = make.unique(entrezID2geneSymbol)
### end gene symbol retrieval ###

mtx
table(colnames(mtx))

### do deSeq2 analysis ###
output_dir = path
control = c("response")
cases = c("nonresponse")
for (case in cases){
  PERFORM_DESEQ2(mtx,control,case,output_dir)
}
### end do deSeq2 analysis ###

##########################

#### 
gender_info = c("Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: female","Sex: male","Sex: female","Sex: male","Sex: female","Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: female","Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: female","Sex: female","Sex: male","Sex: female","Sex: male","Sex: male","Sex: male","Sex: male","Sex: female","Sex: female","Sex: male","Sex: male","Sex: male","Sex: female","Sex: male","Sex: male","Sex: female","Sex: female","Sex: male","Sex: female","Sex: male","Sex: male","Sex: female","Sex: male","Sex: male","Sex: female","Sex: male","Sex: female","Sex: female","Sex: female","Sex: female","Sex: male","Sex: male","Sex: male","Sex: female","Sex: female","Sex: female","Sex: male","Sex: male","Sex: female","Sex: male","Sex: male","Sex: female","Sex: male","Sex: male","Sex: female","Sex: male","Sex: female","Sex: male","Sex: female","Sex: female","Sex: male","Sex: male","Sex: male","Sex: male","Sex: female","Sex: female","Sex: male","Sex: female","Sex: male","Sex: female","Sex: female","Sex: female","Sex: female","Sex: female","Sex: female","Sex: male","Sex: male","Sex: male","Sex: male","Sex: male","Sex: female","Sex: male","Sex: male","Sex: male","Sex: female","Sex: male","Sex: male","Sex: male","Sex: female","Sex: female","Sex: male","Sex: female","Sex: female","Sex: male","Sex: female","Sex: male","Sex: male","Sex: female","Sex: unknown","Sex: female","Sex: female","Sex: male","Sex: female","Sex: female","Sex: male","Sex: male","Sex: female","Sex: female","Sex: male","Sex: male")
gender_info = gsub("Sex: ","",gender_info)
ICU_info = c("icu: no","icu: no","icu: no","icu: no","icu: no","icu: no","icu: no","icu: yes","icu: yes","icu: yes","icu: no","icu: yes","icu: no","icu: yes","icu: yes","icu: no","icu: no","icu: no","icu: no","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: no","icu: no","icu: no","icu: no","icu: no","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: no","icu: no","icu: no","icu: yes","icu: no","icu: no","icu: no","icu: no","icu: no","icu: no","icu: yes","icu: yes","icu: yes","icu: no","icu: no","icu: yes","icu: yes","icu: no","icu: no","icu: yes","icu: no","icu: no","icu: yes","icu: no","icu: yes","icu: yes","icu: no","icu: no","icu: yes","icu: no","icu: no","icu: no","icu: yes","icu: no","icu: yes","icu: no","icu: no","icu: no","icu: yes","icu: yes","icu: no","icu: no","icu: no","icu: yes","icu: no","icu: yes","icu: no","icu: yes","icu: no","icu: no","icu: no","icu: yes","icu: yes","icu: no","icu: no","icu: yes","icu: yes","icu: yes","icu: no","icu: yes","icu: no","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: yes","icu: no","icu: no","icu: no","icu: no","icu: no","icu: yes")
ICU_info = gsub("icu: ","",ICU_info)
disease_info = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27","C28","C29","C30","C31","C32","C33","C34","C35","C36","C37","C38","C39","C40","C41","C42","C43","C44","C45","C46","C47","C48","C49","C50","C51","C52","C53","C55","C56","C57","C58","C59","C60","C61","C62","C63","C64","C65","C66","C67","C68","C69","C70","C71","C72","C73","C74","C75","C76","C77","C78","C79","C80","C82","C83","C84","C85","C86","C87","C89","C90","C91","C92","C93","C94","C95","C96","C97","C98","C99","C100","C101","C102","C103","NC1","NC2","NC3","NC4","NC5","NC6","NC7","NC8","NC9","NC10","NC11","NC12","NC13","NC14","NC15","NC16","NC17","NC18","NC19","NC20","NC21","NC22","NC23","NC24","NC25","NC26")
disease_info = gsub("[0-9]+$","",disease_info)
header = paste0(gender_info,"_",ICU_info,"_",disease_info)
table(header)

# GSE156063, GSE149689 (single), GSE154104 mice, GSE154567 (single cell)

# GSE151803
GSE151803_mtx = read.delim("GSE151803.txt",sep="\t")
colnames(GSE151803_mtx)
new_names = gsub("_[0-9]+$","",colnames(GSE151803_mtx))
new_names = gsub("(_[0-9]+)_","",new_names)
new_names = gsub("\\.", "_", new_names)
colnames(GSE151803_mtx) = new_names
GSE151803_mtx$Name = GSE151803_mtx$GeneID
table(colnames(GSE151803_mtx))
##Donor 1
control = c("D1_Mock")
cases = c("D1_SARS_CoV_2")
for (case in cases){
  PERFORM_DESEQ2(GSE151803_mtx,control,case,output_dir)
}

##Pancreas
control = c("Pancreas_mock")
cases = c("Pancreas_SARSCoV2")
for (case in cases){
  PERFORM_DESEQ2(GSE151803_mtx,control,case,output_dir)
}

##Donor 1 R2
control = c("D1_MockR2")
cases = c("D1_SARS_CoV_2R2")
for (case in cases){
  PERFORM_DESEQ2(GSE151803_mtx,control,case,output_dir)
}



##Donor 2 R2
control = c("D2_MockR2")
cases = c("D2_SARS_CoV_2R2")
for (case in cases){
  PERFORM_DESEQ2(GSE151803_mtx,control,case,output_dir)
}



#### import data from GEO database:  

file_names = c("disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: SC2","disease state: other virus","disease state: other virus","disease state: no virus","disease state: other virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: other virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: SC2","disease state: SC2","disease state: SC2","disease state: no virus","disease state: other virus","disease state: no virus","disease state: other virus","disease state: other virus","disease state: other virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: other virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: other virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: other virus","disease state: other virus","disease state: no virus","disease state: other virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: other virus","disease state: SC2","disease state: no virus","disease state: other virus","disease state: other virus","disease state: other virus","disease state: no virus","disease state: no virus","disease state: other virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: other virus","disease state: no virus","disease state: other virus","disease state: no virus","disease state: other virus","disease state: no virus","disease state: other virus","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: other virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: no virus","disease state: other virus","disease state: other virus","disease state: no virus","disease state: SC2","disease state: other virus","disease state: other virus","disease state: other virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: no virus","disease state: other virus","disease state: no virus","disease state: other virus","disease state: SC2","disease state: no virus","disease state: SC2","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: SC2","disease state: no virus","disease state: other virus","disease state: no virus","disease state: SC2","disease state: other virus","disease state: SC2","disease state: other virus","disease state: no virus","disease state: SC2","disease state: SC2","disease state: no virus","disease state: no virus","disease state: SC2","disease state: no virus","disease state: SC2","disease state: other virus","disease state: other virus","disease state: no virus","disease state: no virus","disease state: other virus","disease state: other virus","disease state: SC2","disease state: no virus","disease state: other virus","disease state: SC2","disease state: other virus","disease state: no virus","disease state: SC2","disease state: other virus","disease state: other virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: no virus","disease state: SC2","disease state: SC2","disease state: no virus","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: no virus","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: SC2","disease state: no virus")


#### import data: GSE152418 ###
GSE152418_mtx = read.delim("GSE152418_p20047_Study1_RawCounts.txt",sep="\t")
GSE152418_mtx[1:5,1:5]
### biomaRt gene symbol retrieval 2 issues. 1. old version dplyr "filter_" not compatible to getBM function. 2. too slow
#require("biomaRt")
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#library(tidyverse)
#library(dplyr)
#G_list <- getBM(filters = "ensembl_gene_id", 
#                attributes = c("ensembl_gene_id", "description"),
#                values = GSE152418_mtx$ENSEMBLID, mart = mart,useCache = FALSE)

### EnsDb gene symbol retrieval ###
BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)

# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- as.vector(GSE152418_mtx$ENSEMBLID)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
head(geneIDs1)
GSE152418_mtx[1:5,1:5]
GSE152418_mtx = merge(GSE152418_mtx,geneIDs1,by.x="ENSEMBLID",by.y="GENEID")

new_names = gsub(".+(nC[oO]V).+","\\1",colnames(GSE152418_mtx))
new_names = gsub("nCOV","nCoV",new_names)
new_names = gsub(".+(_).+","\\1",new_names)
new_names[new_names == "nCoV"]=c("severity: Convalescent","severity: Moderate","severity: Severe","severity: Severe","severity: ICU","severity: Severe","severity: Moderate","severity: Severe","severity: ICU","severity: Severe","severity: Moderate","severity: ICU","severity: Severe","severity: ICU","severity: Severe","severity: Severe","severity: Moderate")
new_names = gsub("severity: ","",new_names)
colnames(GSE152418_mtx) = new_names
#rownames(GSE152418_mtx) = make.unique(GSE152418_mtx$SYMBOL)
GSE152418_mtx$Name = GSE152418_mtx$SYMBOL
output_dir = "../output/"
table(new_names)
control = "_"
cases = c("Moderate","Severe","ICU")

for (case in cases){
  PERFORM_DESEQ2(GSE152418_mtx,control,case,output_dir)
}




### import data GSE147507 ###
data_mtx = read.delim("GSE147507_NormalizedCount_gsea.txt",sep="\t")
#colname_new = gsub("Series[[:digit:]]+_","",colnames(data_mtx))
colname_new = gsub("_[[:digit:]]+$","",colnames(data_mtx))

colnames(data_mtx) = colname_new




### process data Series1 ###
control = "NHBE_Mock"
cases = c("NHBE_SARS_CoV_2","NHBE_IAV","NHBE_IAVdNS","NHBE_IFNB_4h","NHBE_IFNB_6h","NHBE_IFNB_12h")

case = cases[1]
ctrl_mtx = data_mtx[,colnames(data_mtx) == control]
case_mtx = data_mtx[,colnames(data_mtx) == case]
combined_mtx = cbind(ctrl_mtx,case_mtx)
design <-c(rep(1,ncol(ctrl_mtx)),rep(2,ncol(case_mtx)))
#fit = lmFit(combined_mtx,design)
#fit = eBayes(fit)
#topTable(fit)

table(colnames(data_mtx))
sum(table(colnames(data_mtx)))-2
### DEseq 2 protocal ###
table(colnames(data_mtx))
output_dir = "../output/"
## Series1 ##
control = "Series1_NHBE_Mock"
cases = c("Series1_NHBE_SARS_CoV_2")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series2 ##
control = "Series2_A549_Mock"
cases = c("Series2_A549_SARS_CoV_2")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series3 ##
control = "Series3_A549_Mock"
cases = c("Series3_A549_RSV")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series4 ##
control = "Series4_A549_Mock"
cases = c("Series4_A549_IAV")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series5 ##
control = "Series5_A549_Mock"
cases = c("Series5_A549_SARS_CoV_2")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series6 ##
control = "Series6_A549_ACE2_Mock"
cases = c("Series6_A549_ACE2_SARS_CoV_2")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series7 ##
control = "Series7_Calu3_Mock"
cases = c("Series7_Calu3_SARS_CoV_2")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series8 ##
control = "Series8_A549_Mock"
cases = c("Series8_A549_RSV","Series8_A549_HPIV3")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series9 ##
control = "Series9_NHBE_Mock"
cases = c("Series9_NHBE_IAV","Series9_NHBE_IAVdNS1","Series9_NHBE_IFNB_4h","Series9_NHBE_IFNB_6h","Series9_NHBE_IFNB_12h")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

## Series16 ##
control = "Series16_A549_ACE2_Mock"
cases = c("Series16_A549_ACE2_SARS_CoV_2","Series16_A549_ACE2_SARS_CoV_2_Rux")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

control = "NHBE_Mock"
cases = c("NHBE_SARS_CoV_2","NHBE_IAV","NHBE_IAVdNS1","NHBE_IFNB_4h","NHBE_IFNB_6h","NHBE_IFNB_12h")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

control = "A549_Mock"
cases = c("A549_SARS_CoV_2","A549_HPIV3","A549_IAV","A549_RSV")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

control = "A549_ACE2_Mock"
cases = c("A549_ACE2_SARS_CoV_2","A549_ACE2_SARS_CoV_2_Rux")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

control = "Calu3_Mock"
cases = c("Calu3_SARS_CoV_2")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}

control = "HealthyLungBiopsy"
cases = c("COVID19Lung")
for (case in cases){
  PERFORM_DESEQ2(data_mtx,control,case,output_dir)
}


