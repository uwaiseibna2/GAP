library(ANN2)
library(dplyr)
library(gtools)
library(parallel)
library(qqman)
library(rjson)

path_source <- commandArgs(trailingOnly = TRUE)[1]
use_tree_features<-as.logical(commandArgs(trailingOnly = TRUE)[2])
start<-as.integer(commandArgs(trailingOnly = TRUE)[3])
end<-as.integer(commandArgs(trailingOnly = TRUE)[4])


path_tree<-"data-raw/tree-feature-all-species.csv"
path_corresponding_gene_file<-"data-raw/corresponding_genes.json"
path_corresponding_genes<-"data-raw/unknown_genes/"
path_to_result<-"data-raw/results/"
path_genes<-"data-raw/known_genes/"

num_species<-34


path_genes<-paste0(path_source,path_genes)
path_tree<-paste0(path_source,path_tree)
path_corresponding_genes<-paste0(path_source,path_corresponding_genes)
corresponding_genes<-fromJSON(file=paste0(path_source,path_corresponding_gene_file))
path_to_results<-paste0(path_source,path_to_result)


gene_info<-function(path){
  return((data.frame(read.table(path), row.names=2))$V1[1])
}
get_data<-function(path,gene,corresponding_genes,tree_flag,path_tree){
  data = data.frame(read.table(paste0(path,gene)), row.names=2)
  corresponding_gene<-names(corresponding_genes)[which(corresponding_genes==as.character(gene))]
  correspoding_data=data.frame(read.table(paste0(path_corresponding_genes,corresponding_gene)),row.names = 2)
  gene_name<-(data$V1)[1]
  data=rbind(data,correspoding_data)
  data = data.frame(lapply(data, function(x){gsub("-", 0, x)}))
  data = data.frame(lapply(data, function(x){gsub("A|T|G|C|a|t|g|c|N|n", 1, x)}))
  unique_cols <- sapply(data, function(x) length(unique(x))) == 1
  data<- data[, !unique_cols]
  data_<-data.matrix(data)[,-c(1)]
  data_<-scale(data_)
  target<-data[,1]
  pca<-prcomp(data_)
  d <- which(cumsum(pca$sdev^2/sum(pca$sdev^2)) >= 0.95)[1]
  PCs<-pca$x[,1:d]
  V<-pca$rotation
  pcs_d<-data.frame(PCs)
  data<-cbind(target,pcs_d)
  if(tree_flag==TRUE)
  {
    dendo_tree<-read.csv(path_tree,row.names = 1)
    de<-data.frame(dendo_tree)
    unique_cols <- sapply(de, function(x) length(unique(x))) == 1
    de<- de[, !unique_cols]
    data<-cbind(data,de)
  }
  data<-na.omit(data)
  
  return(data)
  
}
fit.model <-function(data,index,results,hl,use_tree_features){
  nFolds <- 34
  myFolds <- cut(seq(1, nrow(data)),
                 breaks = nFolds,
                 labels=FALSE)
  result_no<-index
  flag<-TRUE
  misclassified<-0
  for (i in 1:nFolds) {
    testObs  <- which(myFolds == i, arr.ind = TRUE)

    dfTest   <- data[ testObs, ]
    dfTrain  <- data[-testObs, ]

    train_X<-data.matrix(dfTrain[,-c(1)])

    train_y<-as.numeric(dfTrain$target)

    test_X<-data.matrix(dfTest[,-c(1)])
    test_y<-as.numeric(dfTest$target)
    hl<-(ifelse(is.na(hl),hl,c(hl)))
    if(!is.na(hl)){hl<-list(hl)}
    nn.train<- neuralnetwork(X = train_X, y = train_y,
                             hidden.layers = hl[0],
                             val.prop = 0,
                             optim.type = 'adam',
                             loss.type = "log",
                             batch.size = 33,
                             activ.functions = "relu",
                             standardize = FALSE,
                             learn.rates = 1e-2,
                             L1=results[index,]$L1,
                             L2=results[index,]$L2,
                             n.epochs = 500,
                             verbose = FALSE,
                             random.seed = 1)
    logit.prob <- predict(nn.train, newdata = test_X)
    pred.class <- ifelse(logit.prob$probabilities[1]>logit.prob$probabilities[2], 0, 1)
    if (is.na(pred.class))
    {
      flag<-FALSE
      print("NAN encountered")
      return(34)
    }
    tryCatch(
      expr = {
        if (pred.class!=test_y)
        {
          misclassified=misclassified+1
        }
      },
      error= function(e)
      {})

  }
  return (misclassified)
}

get_azero<-function(path_gene,gene,results,tree_features,path_tree){
  if(length((which(results$CV==0)))==0)
  {
    return(FALSE)
  }
  set<-which(results$CV==0)
  for (subset in set)
  {
    hl<-results[subset,]$HL
    hl<-ifelse(hl!= 0,hl, NA)
    miss<-fit.model(get_data(path_gene,gene,corresponding_genes,tree_features,path_tree),subset,results,hl,tree_features)
    if(miss==0)
    {
      return (TRUE)
      break
    }
    else
    {
      return (FALSE)
    }
  }
}
find_associated_genes<-function(path_to_result,start_range,end_range,path_genes,tree_feats,path_tree){
  expected_files <- sprintf("%d.csv", start_range:end_range)
  files_in_folder <- list.files(path = path_to_results)
  actual_gene_numbers <- as.numeric(sub("(\\d+)\\.csv", "\\1", files_in_folder))
  missing<- setdiff(expected_files, actual_gene_numbers)



  genes<-list()
  for (gene in start_range:end_range){
    if(gene %in% missing){next}
    results<-data.frame(read.csv(paste0(path_to_results,gene,'.csv')))
     if(get_azero(path_genes,gene,results,tree_feats,path_tree)==TRUE)
    {
      append(genes,gene_info(paste0(path_genes,gene)))
    }
  }
  if(length(genes) > 0){write.csv(genes,paste(path_to_result,'associated_genes.csv'),row.names = FALSE)}
  else {print("No associated genes found.")}


}
find_associated_genes(path_to_results,start,end,path_genes,use_tree_features,path_tree)
