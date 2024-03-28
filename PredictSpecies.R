packages<-c('ANN2','dplyr','gtools','parallel','qqman','rjson','BiocManager','tidyr')

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
library(tidyr)
if(!requireNamespace('Biostrings',quietly = TRUE)){
BiocManager::install("Biostrings")}
library(Biostrings)
read_fasta <- function(file_path) {
  # Read the sequences from the FASTA file
  fasta_sequences <- readDNAStringSet(file_path)
  return(fasta_sequences)
}
extract_sequence <- function(sequences, transcript_id) {
  matching_ids <- grep(transcript_id, names(sequences), value = TRUE)
  check_lamprey <- paste(transcript_id, 'lamprey')
  matching_ids <- matching_ids[matching_ids != check_lamprey]
  if (length(matching_ids) > 0) {
    sequence_list <- lapply(matching_ids, function(id) unlist(strsplit(as.character(sequences[[id]]), "")))
    return(list(matching_ids, sequence_list))
  } else {
    print(paste("Transcript ID", transcript_id, "not found."))
    return(NULL)
  }
}

path_source <- commandArgs(trailingOnly = TRUE)[1]
use_tree_features<-as.logical(commandArgs(trailingOnly = TRUE)[2])
transcript_id<-commandArgs(trailingOnly = TRUE)[3]
path_status<-commandArgs(trailingOnly = TRUE)[4]
path_file<-commandArgs(trailingOnly = TRUE)[5]
path_tree<-commandArgs(trailingOnly = TRUE)[6]

path_to_results<-"results/"




path_file<-paste0(path_source,path_file)
path_tree<-paste0(path_source,path_tree)
path_to_results<-paste0(path_source,path_to_results)
path_status<-paste0(path_source,path_status)


get_num_species<-function(path_status)
{
phylogeny<-read.table(path_status,header=1)
colnames(phylogeny)<-c('species_name','target')
return(sum(!is.na(phylogeny$target)))
}
num_species<-get_num_species(path_status)


find_weights<-function(nn.params,use_tree_features,V,d){
  weights<-nn.params[[1]][[1]][1,]
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
generate_hidden_layer_tuples <- function(n, x) {
  if (n == 0) {
    return(NA)
  } else {
    return(permutations(x, n))
  }}

get_data<-function(transcript_id,path_status,path_file,tree_flag,path_tree)
{
  sequence <- extract_sequence(read_fasta(path_file), transcript_id)
  sequence_list_chars <- lapply(sequence[[2]], function(seq) unlist(strsplit(as.character(seq), "")))
  sequence_df <- as.data.frame(do.call(rbind, sequence_list_chars))
  sequence_df<-cbind(sequence[[1]],sequence_df)
  sequence_df <- separate(sequence_df, 'sequence[[1]]', into = c("transcript_id", "species_name"), sep = " ")
  phylogeny<-read.table(path_status,header=1)
  colnames(phylogeny)<-c('species_name','target')
  merged<-merge(phylogeny,sequence_df,by='species_name')
  merged <- merged[, c(3, 1, 2, seq(4, ncol(merged)))]
  data<-merged[,-c(1,2,3)]
  data = data.frame(lapply(data, function(x){gsub("-", 0, x)}))
  data = data.frame(lapply(data, function(x){gsub("A|T|G|C|a|t|g|c|N|n", 1, x)}))
  unique_cols <- sapply(data, function(x) length(unique(x))) == 1
  data<- data[, !unique_cols]
  data<-data.matrix(data)
  data<-scale(data)
  pca<-prcomp(data)
  d <- which(cumsum(pca$sdev^2/sum(pca$sdev^2)) >= 0.95)[1]
  PCs<-pca$x[,1:d]
  V<-pca$rotation
  pcs_d<-data.frame(PCs)
  target<-merged[,c(3)]
  data<-cbind(target,pcs_d)
  if(tree_flag==TRUE)
  {
    dendo_tree<-read.csv(path_tree,header = 1)
    de<-data.frame(dendo_tree)
    unique_cols <- sapply(de, function(x) length(unique(x))) == 1
    de<- de[, !unique_cols]
    de<-de[match(de$species_name,merged$species_name),]
    de<-de[,-c(1)]
    data<-cbind(data,de)
  }
  all_data<-data[,-c(1)]
  data<-na.omit(data)
  return (list(data,merged[,c(1,2)],V,d,all_data))
  
}

train_NN<-function(data,hidden_layer,num_species){
  nFolds <- num_species
  cv_error<-list()
  myFolds <- cut(seq(1, nrow(data)),
                 breaks = nFolds,
                 labels=FALSE)
  l1Vals = 10^seq(-4, 3, length.out = 100)
  alphas = seq(0, 1, 0.05)
  hl= generate_hidden_layer_tuples(hidden_layer,num_species-1)
  print(paste0("hidden_layers: ",dim(hl)))
  count<-0
  print(paste0('Gene ',transcript_id,' running: '))
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
  print(paste0(transcript_id,' completed'))
  return (do.call(rbind, lapply(cv_error, data.frame)))}
fit.model <-function(data,index,results,hl,use_tree_features,all_data,species){
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
    predictions=predict(nn.train,all_data)
    predicted_pheno=cbind(species,predictions$predictions)
    print(predicted_pheno)
    weights = nn.train$Rcpp_ANN$getParams()$weights
    return(list(misclassified,predicted_pheno,weights))
  }
  return (misclassified)
}
get_zero<-function(results,data,all_data,species)
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
    miss<-fit.model(data,subset,results,hl,use_tree_features,all_data,species)

  }
  if(length(miss)>1)
  {
    return (list(TRUE,hl,miss[2],miss[3]))
  }
}


get_gene_architecture<-function(use_tree_features,hidden_layers){
  if(hidden_layers>0)
  {
  print(paste0('This might take a while as we will try all different permutations of ',hidden_layers,' hidden layer architectures with each layer having nodes ranging [1 to number of species]'))
  }
  dataset<-get_data(transcript_id, path_status,path_file,use_tree_features,path_tree)
  x<-get_zero(train_NN(dataset[[1]],hidden_layers,num_species),dataset[[1]],dataset[[5]],dataset[[2]])
  if(x[[1]])
  {
    print(paste0('found architecture with CV error = 0 with architecture: ',ifelse(is.na(x[[2]]),0,x[[2]])))
    if (!dir.exists(path_to_results)) 
    {
      dir.create(path_to_results)
    }
    write.csv(x[[3]],paste0(path_to_results,'predictions.csv'))
    weights<-find_weights(x[[4]],use_tree_features,dataset[[3]],dataset[[4]])
    writeLines(as.character(weights),paste0(path_to_results,'weights.txt'))
    return(TRUE)
  }
  else
    {
      return(FALSE)
    }
  }


for (hl in 0:3) {
  if(get_gene_architecture(use_tree_features,hl))
  {break}
}
