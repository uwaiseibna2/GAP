
packages<-c('ANN2','dplyr','gtools','parallel','qqman','rjson')

for (package in packages){
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = 'http://cran.us.r-project.org')
  }
}
library(ANN2)
library(dplyr)
library(gtools)
library(parallel)
library(qqman)
library(rjson)

#modify the path to source code
path_source <- commandArgs(trailingOnly = TRUE)[1]
use_tree_features<-as.logical(commandArgs(trailingOnly = TRUE)[2])
hidden_layers<-as.integer(commandArgs(trailingOnly = TRUE)[3])
#modify GULO if you want to run different genes
gulo<-21038
path_tree<-"data-raw/tree-feature-all-species.csv"
path_corresponding_gene_file<-"data-raw/corresponding_genes.json"
path_corresponding_genes<-"data-raw/unknown_genes/"
path_to_results<-"data-raw/results/"
path_genes<-"data-raw/known_genes/"
path_gulo<-paste0(path_genes,gulo)

num_species<-34






path_gene<-paste0(path_source,path_gulo)
path_tree<-paste0(path_source,path_tree)

corresponding_genes<-fromJSON(file=paste0(path_source,path_corresponding_gene_file))

find_weights<-function(nn.params,use_tree_features,V,d){
  weights<-nn.params[[1]][1,]
  weights<-abs(weights)
  if(length(nn.params)>1)
  {
    for(i in 2:length(nn.params))
    {
      weights<-t(weights)%*% t(nn.params[[i]])
    }
  }
  
  if (use_tree_features)
  {
    pca_weights<- ((V[,1:d])%*%as.matrix(weights[1:d]) )
  }
  else
  {
    pca_weights<- ((V[,1:d]))%*%as.matrix(weights)}
  return (pca_weights)
  
}
plot_manhattan<-function(weights){
  betas<-as.numeric(weights)
  dataf<-data.frame(betas)
  p<-length(betas)
  d2 <- mahalanobis( x = as.matrix(betas), center = mean(betas), cov = var(betas) )
  f.stat <- (p - 1) * d2 / ( p-1 )
  neglog10pval <- -pf(f.stat, 1, p-1, lower.tail = FALSE, log.p = TRUE)
  dataf$P= 10^-(neglog10pval)
  dataf$P=p.adjust(dataf$P,method = "BH")
  dataf$position=seq(1,length(neglog10pval),1)
  dataf$BP<-(dataf$position)
  dataf$SNP=dataf$position
  
  
  # Assign 'CHR' values based on 'pos' ranges
  dataf$CHR[dataf$pos >= 1 & dataf$pos <= 2] <- 1
  dataf$CHR[dataf$pos >= 3 & dataf$pos <= 108] <- 2
  dataf$CHR[dataf$pos >= 109 & dataf$pos <= 233] <- 3
  dataf$CHR[dataf$pos >= 234 & dataf$pos <= 326] <- 4
  dataf$CHR[dataf$pos >= 327 & dataf$pos <= 420] <- 5
  dataf$CHR[dataf$pos >= 421 & dataf$pos <= 605] <- 6
  dataf$CHR[dataf$pos >= 606 & dataf$pos <= 711] <- 7
  dataf$CHR[dataf$pos >= 712 & dataf$pos <= 788] <- 8
  dataf$CHR[dataf$pos >= 789 & dataf$pos <= 938] <- 9
  dataf$CHR[dataf$pos >= 939 & dataf$pos <= 1104] <- 10
  dataf$CHR[dataf$pos >= 1105 & dataf$pos <= 1195] <- 11
  dataf$CHR[dataf$pos >= 1196 & dataf$pos <= 1323] <- 12
  manhattan(dataf[,c("BP","CHR","P","SNP")],logp = TRUE,col = c('blue4','grey'),suggestiveline = 1.302)
}
generate_hidden_layer_tuples <- function(n, x) {
  if (n == 0) {
    return(NA)
  } else {
    return(permutations(x, n))
  }}
gene_info<-function(path){
  return((data.frame(read.table(path), row.names=2))$V1[1])
}
get_data<-function(path,gulo,corresponding_genes,tree_flag,path_tree){
  data = data.frame(read.table(path), row.names=2)
  corresponding_gene<-names(corresponding_genes)[which(corresponding_genes==as.character(gulo))]
  correspoding_data=data.frame(read.table(paste0(path_source,path_corresponding_genes,corresponding_gene)),row.names = 2)
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
get_vd<-function(path){
  data = data.frame(read.table(path), row.names=2)
  corresponding_gene<-names(corresponding_genes)[which(corresponding_genes==as.character(gulo))]
  correspoding_data=data.frame(read.table(paste0(path_source,path_corresponding_genes,corresponding_gene)),row.names = 2)
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
  V<-V[,1:34]
  return (list(V,d))
}
train_NN<-function(data,hidden_layer,num_species,path_gene){
  nFolds <- num_species
  cv_error<-list()
  myFolds <- cut(seq(1, nrow(data)),
                 breaks = nFolds,
                 labels=FALSE)
  l1Vals = 10^seq(-4, 3, length.out = 100)
  alphas = seq(0, 1, 0.05)
  hl= generate_hidden_layer_tuples(hidden_layer,num_species-1)
  count<-0
  print(paste0('Gene ',gene_info(path_gene),' running: '))
  for( l in l1Vals)
  {
    for( alpha in alphas)
    {
      if(is.null(dim(hl) ))
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
                                   hidden.layers = NA,
                                   val.prop = 0,
                                   optim.type = 'adam',
                                   loss.type = "log",
                                   batch.size = num_species-1,
                                   standardize = FALSE,
                                   learn.rates = 0.01,
                                   L1=l1,
                                   L2=l2,
                                   n.epochs = 500,
                                   verbose = FALSE,
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
          cv_error <- append(cv_error,list(list(L1=l1,L2 = l2,alpha=alpha,CV=(misclassified/num_species))))
        }
        count<-count+1
        print(paste0('Running architecture ',as.character(count)))
      }
      else
      {
        for(j in 1:nrow(hl))
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
                                     hidden.layers = hl[j,],
                                     val.prop = 0,
                                     optim.type = 'adam',
                                     loss.type = "log",
                                     batch.size = num_species-1,
                                     standardize = FALSE,
                                     learn.rates = 0.01,
                                     L1=l1,
                                     L2=l2,
                                     n.epochs = 500,
                                     verbose = FALSE,
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
            cv_error <- append(cv_error,list(list(L1=l1,L2 = l2,alpha=alpha,layers=hl[j,],CV=(misclassified/num_species))))
          }
          count<-count+1
          print(paste0('Running architecture ',as.character(count)))
        }

      }

    }
  }
  print(paste0((gene_info(path_gene)),' completed'))
  return (do.call(rbind, lapply(cv_error, data.frame)))}
fit.model <-function(data,index,results,hl,use_tree_features){
  print(results[index,])
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


    nn.train<- neuralnetwork(X = train_X, y = train_y,
                             hidden.layers = hl,
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
  
  if(misclassified==0){

    vd<-get_vd(path_gene)
    manhattan_plot<-plot_manhattan(find_weights(nn.train$Rcpp_ANN$getParams()$weights,use_tree_features,vd[[1]],vd[[2]]))
    
    png(file = path_source)
    print(manhattan_plot)
    dev.off()
  }
  return (misclassified)
}
get_zero<-function(results)
{

  if(length((which(results$CV==0)))==0)
  {
    return(list(FALSE,hl))
  }
  set<-which(results$CV==0)
  for (subset in set)
  {
    if(is.null(results[subset,]$layers))
    {
      hl<-NA
    }
    else
    {

      hl<-results[subset,]$layers
    }
    miss<-fit.model(get_data(path_gene,gulo,corresponding_genes,use_tree_features,path_tree),subset,results,hl,use_tree_features)
    if(miss==0)
    {
      return (list(TRUE,hl))
    }
  }
}


get_gene_architecture<-function(path_gene,use_tree_features,hidden_layers){
  if(hidden_layers>0)
  {
  print(paste0('This might take a while as we will try all different permutations of ',hidden_layers,' hidden layer architectures with each layer having nodes ranging [1 to number of species]'))
  }
  x<-get_zero(train_NN(get_data(path_gene,gulo,corresponding_genes,use_tree_features,path_tree),hidden_layers,num_species,path_gene))
  if(x[[1]])
  {
    print(paste0('found architecture with CV error = 0 with architecture: ',ifelse(x[[2]]=='NA',0,x[[2]])))
    hidden_layers<-x[[2]]
  }
  else
  {
    print('Please try different hidden layer architectures!')
  }}



get_gene_architecture(path_gene,use_tree_features,hidden_layers)