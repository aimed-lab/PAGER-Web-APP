if("biomaRt" %in% rownames(installed.packages()) == FALSE) {source("http://bioconductor.org/biocLite.R");biocLite('biomaRt')}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {source("http://bioconductor.org/biocLite.R");biocLite('jsonlite')}
library(biomaRt)
library(jsonlite)
library(httr)

# input parameters:
#

matchedGene = function(
  marker=c("CMM","ENO1","MTOR","MIIP","PLA2G2A"),
  exptag="yes",
  ppiCutoff=0.45,
  ppiExtended=0.9
){
  #markers=as.matrix(R_SltCnd)
  marker_str = paste0(as.matrix(marker),collapse = "%20")
  #print(marker_str)
  contnt <- list(
    genes=marker_str,
    extended=exptag,#'yes',
    ppi=ppiCutoff,
    ppiExtended=ppiExtended
  )
  #http://discovery.informatics.uab.edu/BEERE/index.php/search/matchedgene
  #http://127.0.0.1:8084/beere/index.php/search/matchedgene
  req <- httr::POST("http://discovery.informatics.uab.edu/BEERE/index.php/search/matchedgene",#"http://138.26.131.247:8800/list/",
                    body=contnt
  );
  
  json <- httr::content(req, as = "text")
  
  PAGER=fromJSON(json)
  #write.table(PAGER,paste0(fileName,"_",paste0(as.matrix(""),collapse = "_"),".txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}
#seedCat=c('BRCA1','BRCA2',"AF10",'MLLT3')
#PPI=matchedGene(seedCat,'no','0.45','0')
#PPI

# input parameters: 
# version "_V40","_V31","_V30"
# matchVersion "pattern", "fuzzy","exact"
#markers=UPlist
SEMMEDsearch = function(
  markers=c("Alzheimer's", "Vitamin A","Vitamin D","Vitamin E","Vitamin K"),
  fileName="./",
  version="_V40",
  matchVersion="fuzzy",
  exttag='yes',
  extended="yes",
  ppi=1,
  ppiExtended = 10,
  identifierSub = "term"
){
  #markers=as.matrix(R_SltCnd)
  marker_str = toJSON(unlist(markers), pretty=TRUE)
  #print(marker_str)
  contnt <- list(
    genes=marker_str,
    fuzzy=matchVersion,#'pattern'
    extended=exttag,#s'yes',
    ppi=ppi,
    identifierSub=identifierSub,
    version=version,
    ppiExtended=ppiExtended
  )
  #http://discovery.informatics.uab.edu/BEERE/index.php/search/SEMMEDsearch
  #http://127.0.0.1:8084/beere/index.php/search/SEMMEDsearch
  req <- httr::POST("http://discovery.informatics.uab.edu/BEERE/index.php/search/SEMMEDsearch",#"http://138.26.131.247:8800/list/",
                    body=contnt
  );
  
  json <- httr::content(req, as = "text")
  
  PAGER=fromJSON(json)
  
  #write.table(PAGER,paste0(fileName,"_",paste0(as.matrix(""),collapse = "_"),".txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}
#resMtxUP$SYM=resMtxUP$SYM
#markers=resMtx$matched

# ext_type: 'ant',''pagerank'
DDTsemmed = function(
  markers=c("Alzheimer's", "Vitamin A","Vitamin D","Vitamin E","Vitamin K"),
  fileName='./',
  ppiCut=0,
  ExtCut=10,
  version='_V40',
  exttag='yes',
  exp_type=c(''),
  predicate=c('')
  ){
  #markers=as.matrix(R_SltCnd)
  markers=unique(markers)[!unique(markers)%in%"-"]
  marker_str = toJSON(unlist(markers), pretty=TRUE)
  exp_type_str = toJSON(unlist(exp_type), pretty=TRUE)
  predicate_str = toJSON(unlist(predicate), pretty=TRUE)
  #print(predicate_str)
  contnt <- list(
    genes=marker_str,
    #fuzzy='pattern',
    extended=exttag,#'yes',
    ppi=ppiCut,
    #identifierSub="term",
    version=version,
    ppiExtended=ExtCut,
    exp_type=exp_type_str,
    predicateSet=predicate_str
  )
  
  #req <- httr::POST("http://discovery.informatics.uab.edu/BEERE/index.php/search/SEMMEDttRelNew",#"discovery.informatics.uab.edu http://138.26.131.247:8800/list/",
  #                  body=contnt
  #);
  #http://discovery.informatics.uab.edu/BEERE/index.php/search/SEMMEDttRelNew
  #http://127.0.0.1:8084/BEERE/index.php/search/SEMMEDttRelNew
  req <- httr::POST("http://discovery.informatics.uab.edu/BEERE/index.php/search/SEMMEDttRelNew",#"discovery.informatics.uab.edu http://138.26.131.247:8800/list/",
                    body=contnt
  );
  json <- httr::content(req, as = "text")
  
  PAGER=fromJSON(json)
  
  #write.table(PAGER,paste0(fileName,"_",paste0(as.matrix(""),collapse = "_"),".txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}


Ranksemmed = function(
  CbnInt=DataFrame(
    SYM_A = c("Vitamin A","Vitamin E"),
    SYM_B = c("Vitamin E","Alzheimer's Disease"),
    SCORE = c(31,5),
    PREDICATE = c('compared_with','DISRUPTS')
  ),
  markers=c("Alzheimer's", "Vitamin A","Vitamin D","Vitamin E","Vitamin K"),
  method='ant',
  sigma=0.8,
  iteration=10
){
  #markers=as.matrix(R_SltCnd)
  queryStr=paste(paste0(CbnInt$SYM_A,"\t\t",CbnInt$SYM_B,"\t\t",CbnInt$SCORE,"\t\t",CbnInt$PREDICATE),collapse = "\t\t\t")
  geneStr=paste(markers,collapse = "\t\t")
  #print(marker_str)
  contnt <- list(
    customPPI=queryStr,
    method=method,#'ant',
    sigma=sigma,#0.8,
    genes=geneStr,
    iteration=iteration,#10,
    extended="",
    ppi=""
  )
  #http://discovery.informatics.uab.edu/BEERE/index.php/search/geneRank
  req <- httr::POST("http://discovery.informatics.uab.edu/BEERE/index.php/search/geneRank",#"discovery.informatics.uab.edu http://138.26.131.247:8800/list/",
                    body=contnt
  );
  
  json <- httr::content(req, as = "text")
  
  PAGER=fromJSON(json)
  
  #write.table(PAGER,paste0(fileName,"_",paste0(as.matrix(""),collapse = "_"),".txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}

PMIDsemmed = function(
  SYM_A=SYM_A,
  SYM_B=SYM_B,
  predicate=predicate,
  version=version){
  #markers=as.matrix(R_SltCnd)
  
  #print(marker_str)
  contnt <- list(
    SYM_A=SYM_A,
    SYM_B=SYM_B,
    predicate=predicate,
    version=version
  )
  
  req <- httr::POST("http://discovery.informatics.uab.edu/BEERE/index.php/search/SEMMEDsearchPMID",#"discovery.informatics.uab.edu http://138.26.131.247:8800/list/",
                    body=contnt
  );
  
  json <- httr::content(req, as = "text")
  
  PAGER=fromJSON(json)
  
  #write.table(PAGER,paste0(fileName,"_",paste0(as.matrix(""),collapse = "_"),".txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)
}

# 01/12/21
Pubtator = function(SYM_A=SYM_A,SYM_B=SYM_B){

  #print(marker_str)
  contnt <- list(
    SYM_A=SYM_A,
    SYM_B=SYM_B
  )
  
  req <- httr::POST("http://discovery.informatics.uab.edu/BEERE/index.php/search/PPIPUBMED",#"discovery.informatics.uab.edu http://138.26.131.247:8800/list/",
                    body=contnt
  );
  
  json <- httr::content(req, as = "text")
  
  PAGER=fromJSON(json)
  
  #write.table(PAGER,paste0(fileName,"_",paste0(as.matrix(""),collapse = "_"),".txt"), row.names = F, quote =F,col.names = T,sep="\t")
  return(PAGER)  
}
####