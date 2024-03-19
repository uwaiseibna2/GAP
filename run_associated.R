library(ANN2)
library(dplyr)
library(parallel)
library(rjson)


path_source <- commandArgs(trailingOnly = TRUE)[1]
use_dendrogram_features<-as.logical(commandArgs(trailingOnly = TRUE)[2])
hidden_layers<-commandArgs(trailingOnly = TRUE)[3]
start<-as.integer(commandArgs(trailingOnly = TRUE)[4])
end<-as.integer(commandArgs(trailingOnly = TRUE)[5])


path_tree<-"data-raw/tree-feature-all-species.csv"
path_corresponding_gene_file<-"data-raw/corresponding_genes.json"
path_corresponding_genes<-"data-raw/unknown_genes/"
path_to_results<-"data-raw/results/"
path_genes<-"data-raw/known_genes/"


num_species<-34


path_tree<-paste0(path_source,path_tree)
path_genes<-paste0(path_source,path_genes)
corresponding_geneset<-fromJSON(file=paste0(path_source,path_corresponding_gene_file))
path_corresponding_genes<-paste0(path_source,path_corresponding_genes)
path_to_results<-paste0(path_source,path_to_results)

gene_info<-function(path){
  return((data.frame(read.table(path), row.names=2))$V1[1])
}
get_data<-function(path,gene,corresponding_genes,tree_flag,path_tree){
  data = data.frame(read.table(path), row.names=2)
  corresponding_gene<-names(corresponding_geneset)[which(corresponding_geneset==as.character(gene))]
  correspoding_data=data.frame(read.table(paste0(path_corresponding_genes,corresponding_gene)),row.names = 2)
  gene<-(data$V1)[1]
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
  
  if(use_dendrogram_features==TRUE)
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
train_ANN<-function(data,hl,num_species,path_to_result,num){
  nFolds <- num_species
  cv_error<-list()
  myFolds <- cut(seq(1, nrow(data)),
                 breaks = nFolds,
                 labels=FALSE)
  l1Vals = 10^seq(-4, 3, length.out = 100)
  alphas = seq(0, 1, 0.05)
  count<-0
  for( l in l1Vals)
  {
    for( alpha in alphas)
    {
      l1=alpha*l
      l2=(1-alpha)*l
      misclassified<-0
      flag<-TRUE
      for (i in 1:nFolds) {
        testObs  <- which(myFolds == i, arr.ind = TRUE)

        dfTest   <- data[ testObs, ]
        dfTrain  <- data[-testObs, ]

        train_X<-data.matrix(dfTrain[,-c(1)])

        train_y<-as.numeric(dfTrain$target)

        test_X<-data.matrix(dfTest[,-c(1)])
        test_y<-as.numeric(dfTest$target)


        nn.train<- neuralnetwork(X = train_X, y = train_y,
                                 hidden.layers = ifelse(hl != "0", as.numeric(strsplit(hl, ",")[[1]]), NA),
                                 val.prop = 0,
                                 optim.type = 'adam',
                                 loss.type = "log",
                                 batch.size = num_species-1,
                                 standardize = FALSE,
                                 learn.rates = 0.01,
                                 L1=l1,
                                 L2=l2,
                                 n.epochs = 500,
                                 verbose = TRUE,
                                 random.seed = 1)
        logit.prob <- predict(nn.train, newdata = test_X)
        pred.class <- ifelse(logit.prob$probabilities[1]>logit.prob$probabilities[2], 0, 1)
        if (is.na(pred.class))
        {
          flag<-FALSE
          break
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
      if(flag!=FALSE)
      {
        cv_error <- append(cv_error,list(list(L1=l1,L2 = l2,alpha=alpha,CV=(misclassified/num_species),HL=hl)))
      }
    }
    count<-count+1
    print(paste0('architecture ',as.character(count),' completed.'))
  }
  results<-(do.call(rbind, lapply(cv_error, data.frame)))
  results<-apply(results,2,as.character)
  write.csv(results,paste0(path_to_result,num,'.csv'),row.names = FALSE)
}


f <- function(i) {
  cat(paste0("Working on gene number: ", i, "\n"))
  data <- get_data(paste0(path_genes, i),i,path_corresponding_genes, use_dendrogram_features, path_tree)
  
  train_ANN(data, hidden_layers, num_species, path_to_results, i)
  cat(paste0("Finished working on gene number: ", i, "\n"))
}

run_associated_genes <- function(start, end) {
  cat("Starting parallel processing...\n")
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(ANN2))
  clusterEvalQ(cl, library(rjson))
  clusterExport(cl, list("corresponding_geneset","path_source","path_genes","path_corresponding_genes","path_corresponding_gene_file", "gene_info", "get_data", "path_to_results", "path_tree", "hidden_layers", "use_dendrogram_features", "train_ANN", "num_species"))
  clusterApply(cl, start:end, f)
  stopCluster(cl)
  cat("Parallel processing completed.\n")
}


run_associated_genes(start,end)
