if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install('biomaRt')
}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install('jsonlite')
}
if("httr" %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install('httr')
}


library(biomaRt)
library(jsonlite)
library(httr)

# pathFun is a function connected to PAGER api to perform hypergeometric test and retrieve enriched PAGs associated to a list of genes
# 1.markers: a list of gene symbols, 2.filename: provide the name of the file to store, 3.sourceSel: a list of sources refered to 'http://discovery.informatics.uab.edu/PAGER/index.php/search/advanced' data source
# 4.typeSel 'All', 'P', 'A' or 'G', 5.pCutoff: p-valuecutoff, 6. fdrCutoff: FDR cutoff, 7.upperCutoff: maximum size of genes in PAG
# 8.bottomCutoff: minimum size of genes in PAG
#markers=unique(as.matrix(targetFmt))
#typeSel='All'
#upperCutoff=2000
#bottomCutoff=1
#pCutoff=0.05
#fdrCutoff=0.05
# example: 

pathFun = function(dir = './', 
                   fileName = 'test.txt', 
                   markers = c('BRCA1','BRCA2'), 
                   type = "All",
                   minSize = 2, 
                   maxSize = 5000, 
                   similarity = 0.05, 
                   overlap = 1, 
                   nCoCo = 0, 
                   pValue = 0.05,
                   FDR = 0.05, 
                   Species = 'All', 
                   source = c('KEGG',"WikiPathway","BioCarta","MSigDB","NCI-Nature Curated","Reactome","Protein Lounge","Spike")){
  url = "http://discovery.informatics.uab.edu/PAGER/index.php/geneset/pagerapi"
  marker_str = NULL
  marker_str=paste0(as.matrix(markers),collapse = "%20")
  source_str = NULL
  source_str=paste0(as.matrix(source),collapse = "%20")
  #Other adjustable parameters can be refered to http://discovery.informatics.uab.edu/PAGER/index.php/search/advanced
  contnt <- list(
    genes = marker_str,
    type = type,
    ge = minSize,
    le = maxSize,
    sim = similarity,
    olap = overlap,
    cohesion = nCoCo,
    pvalue = pValue,
    FDR = FDR,
    organism = Species,
    source = source_str
  )
  req <- httr::POST(url,
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  PAGER=fromJSON(json)
  write.table(PAGER,paste0(dir,fileName), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}

# PAGERCOV_ANALYSIS is a function connected to PAGER api to perform hypergeometric test and retrieve enriched PAGs associated to a list of genes
# The input parameters are:
# 1.dir: the directory of output file.
# 2.filename: the name of the output file.
# 3.markers: a list of gene symbols.
# 4.type: a list of PAG types consisting of 'P', 'A' or 'G'
# 5.minSize: the allowed minimum size of PAG genes
# 6.maxSize: the allowed maximum size of PAG genes
# 7.similarity: the similarity score cutoff
# 8.overlap: the allowed minimum overlap genes
# 9.nCoCo: the minimum nCoCo score
# 10.pValue: p-value cutoff
# 11.FDR: false discovery rate 
# 12.source: a list of sources refered to 'http://discovery.informatics.uab.edu/PAGER-COV/index.php/pages/help'

PAGERCOV_ANALYSIS = function(dir = './', fileName = 'test.txt', markers = c('BRCA1','BRCA2'), type = c("P","A","G"),
                             minSize = 2, maxSize = 5000, similarity = 0.05, overlap = 1, nCoCo = 0, pValue = 0.05,
                             FDR = 0.05, Species = 'All', source = c('PubChem','PAGER-MSigDB','PAGER-GOA','PAGER-GOA_AGG',
                                                                     'PAGER-GAD','PAGER-Pfam','PAGER-GeneSigDB','PAGER-Protein Lounge',
                                                                     'PAGER-Spike','PAGER-PheWAS','PAGER-PharmGKB','PAGER-Reactome',
                                                                     'PAGER-BioCarta','PAGER-GWAS Catalog','PAGER-KEGG','PAGER-NGS Catalog',
                                                                     'PAGER-DSigDB','PAGER-GTEx','PAGER-NCI-Nature Curated','PAGER-WikiPathway',
                                                                     'Am J Respir Crit Care Med','Microbiology and Molecular Biology Reviews',
                                                                     'Zenodo','Mouse Genome Informatics Database','Nature Cell Discovery',
                                                                     'GenBank/UniProt','Cell','The Annual Review of Cell and Developmental Biology',
                                                                     'Nature','Drugbank','Nature Medicine','Cell Host and Microbe','bioRxiv'
                             )){
  url = "http://discovery.informatics.uab.edu/PAGER-COV/index.php/geneset/pagerapi"
  marker_str = NULL
  marker_str = paste0(as.matrix(markers),collapse = "%20")
  source_str = NULL
  source_str = paste0(as.matrix(source),collapse = "%20")
  contnt <- list(
    genes = marker_str,
    type = type,
    ge = minSize,
    le = maxSize,
    sim = similarity,
    olap = overlap,
    cohesion = nCoCo,
    pvalue = pValue,
    FDR = FDR,
    organism=Species,
    source=source_str
  )
  req <- httr::POST(url,
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  result = fromJSON(json)
  write.table(result,paste0(dir,fileName), row.names = F, quote =F,col.names = T,sep="\t")
  return(result)
}

# pathMember is a function connected to PAGER api to retrieve the membership of PAGs using a list of PAG IDs
# example:
# PAGmember=pathMember(dir = "./",PAG_result$GS_ID)
pathMember = function(dir = './',PAGER){
  PAGstring = paste(PAGER,collapse=",")
  contnt <- list(
    pag=PAGstring
  )
  req <- httr::POST("http://discovery.informatics.uab.edu/PAGER/index.php/geneset/get_members_by_ids/",
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  PAGINF=fromJSON(json)
  write.table(PAGINF,paste0(dir,"PAG_members.txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGINF)
}

# pathInt is a function connected to PAGER api to retrieve the m-type relationships of PAGs using a list of PAG IDs
# example:
# PAGint=pathInt(dir = './',PAG_result$GS_ID)
pathInt = function(dir = './',PAGER){
  
  PAGstring = paste(PAGER,collapse=",")
  
  contnt <- list(
    pag=PAGstring
  )
  req <- httr::POST("http://discovery.informatics.uab.edu/PAGER/index.php/pag_pag/inter_network_int_api/",
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  PAGINF=fromJSON(json)
  write.table(PAGINF,paste0(dir,"mtype.txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGINF)
}


# pathReg is a function connected to PAGER api to retrieve the r-type relationships of PAGs using a list of PAG IDs
# example:
# PAGreg=pathReg(dir = './',PAG_result$GS_ID)
pathReg = function(dir = './',PAGER){
  PAGstring = paste(PAGER,collapse=",")
  contnt <- list(
    pag=PAGstring
  )
  req <- httr::POST("http://discovery.informatics.uab.edu/PAGER/index.php/pag_pag/inter_network_reg_api/",
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  PAGINF=fromJSON(json)
  write.table(PAGINF,paste0(dir,"rtype.txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGINF)
}

# pagRankedGene is a function connected to PAGER api to retrieve RP-ranked genes with RP-score of the given PAG_IDs
# example:
# PAGRankedGenes=pagRankedGene(dir = './',PAG_result$GS_ID[1])
pagRankedGene = function(dir = './',pag_id){
  
  contnt <- list(
    pag=pag_id
  )
  req <- httr::POST("http://discovery.informatics.uab.edu/PAGER/index.php/geneset/pag_ranked_gene/",
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  PAGINF=fromJSON(json)
  write.table(PAGINF,paste0(dir,pag_id,"_pagRankedGene.txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGINF)
}

# pagGeneInt is a function connected to PAGER api to retrieve gene interaction network
# example:
# pagGeneInts=pagGeneInt(dir = './',PAG_result$GS_ID[1])
pagGeneInt = function(dir = './',pag_id){
  req <- httr::GET(paste0("http://discovery.informatics.uab.edu/PAGER/index.php/pag_mol_mol_map/interactions/",pag_id));
  json <- httr::content(req, as = "text")
  PAGGENEINT=fromJSON(json)
  write.table(PAGGENEINT,paste0(dir,pag_id,"_pagGeneInt.txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGGENEINT)
}

# pagGeneReg is a function connected to PAGER api to retrieve gene regulatory network
# example:
# pagGeneRegs=pagGeneReg(dir = './',PAG_result$GS_ID[1])
pagGeneReg = function(dir = './',pag_id){
  req <- httr::GET(paste0("http://discovery.informatics.uab.edu/PAGER/index.php/pag_mol_mol_map/regulations/",pag_id));
  json <- httr::content(req, as = "text")
  PAGGENEINT=fromJSON(json)
  write.table(PAGGENEINT,paste0(dir,pag_id,"_pagGeneReg.txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGGENEINT)
}

# path_NGSEA is a function connected to PAGER api to generate the network-based GSEA result
# example:
# gene_df = data.frame(symbol = c("BRCA1","BRCA2","RAD51","PGR"), log2FoldChange = c(0.6,0.3,0.1,0.1))
# pagInNetworkWGSEA=path_NGSEA(dir = './',gene_df,PAG_member)
path_NGSEA = function(dir = './',genes,PAGmember){
  geneExpStr=paste(paste0(genes$symbol,"\\t\\t",genes$log2FoldChange),collapse = "\\t\\t\\t")
  PAGsetsStr=paste(paste0(PAGmember$data$GS_ID,"\\t\\t",PAGmember$data$GENE_SYM),collapse = "\\t\\t\\t")
  contnt <- list(
    PAGsetsStr=PAGsetsStr,
    geneExpStr=geneExpStr
  )
  req <- httr::POST("http://discovery.informatics.uab.edu/PAGER/index.php/geneset/ngseaapi/",
                    body=contnt
  );
  
  json <- httr::content(req, as = "text")
  PAGER=fromJSON(json)
  write.table(PAGER,paste0(dir,"PAG_InNetworkWeightedgsea.txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}

# CellAnnotation is a function connected to PAGER api to perform hypergeometric test and retrieve enriched PAGs associated to a list of genes
# 1.markers: a list of gene symbols, 2.filename: provide the name of the file to store, 3.pCutoff: p-valuecutoff, 
# 4. fdrCutoff: FDR cutoff, 5.upperCutoff: maximum size of genes in PAG
# 6.bottomCutoff: minimum size of genes in PAG
# example: 

CellAnnotation = function(dir = './', 
                          fileName = 'test.txt', 
                          markers = c('KLRB1','KLRC1','KLRD1'), 
                          minSize = 2, 
                          maxSize = 5000, 
                          similarity = 0.05, 
                          overlap = 1, 
                          pValue = 0.05,
                          FDR = 0.05, 
                          Species = 'All'
){
  url = "http://discovery.informatics.uab.edu/PAGER/index.php/geneset/pagerapi"
  marker_str = NULL
  marker_str=paste0(as.matrix(markers),collapse = "%20")
  #Other adjustable parameters can be refered to http://discovery.informatics.uab.edu/PAGER/index.php/search/advanced
  contnt <- list(
    genes = marker_str,
    type = "All",
    ge = minSize,
    le = maxSize,
    sim = similarity,
    olap = overlap,
    cohesion = 0,
    pvalue = pValue,
    FDR = FDR,
    organism = Species,
    source = 'CellMarker'
  )
  req <- httr::POST(url,
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  PAGER=fromJSON(json)
  write.table(PAGER,paste0(dir,fileName), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}
