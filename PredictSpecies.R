packages<-c('ANN2','gtools','parallel','BiocManager','tidyr')
for (package in packages){
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = 'http://cran.us.r-project.org')
  }
}
library(ANN2)
library(gtools)
library(tidyr)
if(!requireNamespace('Biostrings',quietly = TRUE)){
BiocManager::install("Biostrings")}
library(Biostrings)

#get file paths
path_source <- commandArgs(trailingOnly = TRUE)[1]
use_tree_features<-as.logical(commandArgs(trailingOnly = TRUE)[2])
transcript_id<-commandArgs(trailingOnly = TRUE)[3]
path_status<-commandArgs(trailingOnly = TRUE)[4]
path_file<-commandArgs(trailingOnly = TRUE)[5]
path_tree<-commandArgs(trailingOnly = TRUE)[6]


#implement paths
path_to_results<-"results/"
path_file<-paste0(path_source,path_file)
path_tree<-paste0(path_source,path_tree)
path_to_results<-paste0(path_source,path_to_results)
path_status<-paste0(path_source,path_status)

get_num_species<-function(path_status){
  phylogeny<-read.table(path_status,header=1)
  colnames(phylogeny)<-c('species_name','target')
  return(sum(!is.na(phylogeny$target)))
}
#number of known species
num_species<-get_num_species(path_status)

read_fasta <- function(file_path) {
  # Read the sequences from the FASTA file
  fasta_sequences <- readDNAStringSet(file_path)
  return(fasta_sequences)
}
extract_sequence <- function(seqs, t_id) {
  matching_ids <- grep(t_id, names(seqs), value = TRUE)
  check_lamprey <- paste(t_id, 'lamprey')
  matching_ids <- matching_ids[matching_ids != check_lamprey]
  if (length(matching_ids) > 0) {
    sequence_list <- lapply(matching_ids, function(id) unlist(strsplit(as.character(seqs[[id]]), "")))
    return(list(matching_ids, sequence_list))
  } else {
    print(paste("Transcript ID", t_id, "not found."))
    return(NULL)
  }
}
find_weights<-function(params,tf,V,d){
  weights<-params[[1]][[1]][1,]
  weights<-abs(weights)
  if(length(params)>1)
  {
    for(i in 2:length(params))
    {
      weights<-t(weights)%*% t(params[[i]])
    }
  }
  
  if (tf)
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
get_data<-function(trid,phenoStatus,SeqFile,TreeFlag,TreeFile){
  sequence <- extract_sequence(read_fasta(SeqFile), trid)
  sequence_list_chars <- lapply(sequence[[2]], function(seq) unlist(strsplit(as.character(seq), "")))
  sequence_df <- as.data.frame(do.call(rbind, sequence_list_chars))
  sequence_df<-cbind(sequence[[1]],sequence_df)
  sequence_df <- separate(sequence_df, 'sequence[[1]]', into = c("trid", "species_name"), sep = " ")
  phylogeny<-read.table(phenoStatus,header=1)
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
  if(TreeFlag==TRUE)
  {
    dendo_tree<-read.csv(TreeFile,header = 1)
    de<-data.frame(dendo_tree)
    unique_cols <- sapply(de, function(x) length(unique(x))) == 1
    de<- de[, !unique_cols]
    de<-de[match(merged$species_name,de$species_name),]
    de<-de[,-c(1)]
    data<-cbind(data,de)
  }
  unknownData<-data.frame(data)[is.na(data$target), ]
  unknownData<-unknownData[,-c(1)]
  data<-na.omit(data)
  spcsl<-data.frame(merged)[is.na(merged$target),]
  return (list(data,spcsl[,c(1,2)],V,d,unknownData))
  
}
train_NN<-function(InputNN,num_hidden_layer,numSpcs){
  nFolds <- numSpcs
  cv_error<-list()
  myFolds <- cut(seq(1, nrow(InputNN)),
                 breaks = nFolds,
                 labels=FALSE)
  l1Vals = 10^seq(-4, 3, length.out = 100)
  alphas = seq(0, 1, 0.05)
  hl= generate_hidden_layer_tuples(num_hidden_layer,numSpcs-1)
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

          dfTest   <- InputNN[ testObs, ]
          dfTrain  <- InputNN[-testObs, ]

          train_X<-data.matrix(dfTrain[,-c(1)])

          train_y<-as.numeric(dfTrain$target)

          test_X<-data.matrix(dfTest[,-c(1)])
          test_y<-as.numeric(dfTest$target)


          nn.train<- neuralnetwork(X = train_X, y = train_y,
                                   hidden.layers = NA,
                                   val.prop = 0,
                                   optim.type = 'adam',
                                   loss.type = "log",
                                   batch.size = numSpcs-1,
                                   standardize = FALSE,
                                   learn.rates = 1e-2,
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
          cv_error <- append(cv_error,list(list(L1=l1,L2 = l2,alpha=alpha,CV=(misclassified/numSpcs))))
        }
        count<-count+1
        cat(paste0('\rRunning architecture ', as.character(count),'/',length(alphas)*length(l1Vals)))
        flush.console()
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

            dfTest   <- InputNN[ testObs, ]
            dfTrain  <- InputNN[-testObs, ]

            train_X<-data.matrix(dfTrain[,-c(1)])

            train_y<-as.numeric(dfTrain$target)

            test_X<-data.matrix(dfTest[,-c(1)])
            test_y<-as.numeric(dfTest$target)


            nn.train<- neuralnetwork(X = train_X, y = train_y,
                                     hidden.layers = hl[j,],
                                     val.prop = 0,
                                     optim.type = 'adam',
                                     loss.type = "log",
                                     batch.size = numSpcs-1,
                                     standardize = FALSE,
                                     learn.rates = 1e-2,
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
            cv_error <- append(cv_error,list(list(L1=l1,L2 = l2,alpha=alpha,layers=hl[j,],CV=(misclassified/numSpcs))))
          }
          count=count+1
          cat(paste0('\rRunning architecture ', as.character(count),'/',length(alphas)*length(l1Vals)*nrow(hl)))
          flush.console()
        }

      }

    }
  }
  print(paste0(transcript_id,' completed'))
  res<-do.call(rbind, lapply(cv_error, data.frame))
  write.csv(res,'/Users/uwaiseibna/Downloads/result-tree.csv')
  return (do.call(rbind, lapply(cv_error, data.frame)))}
fit.model <-function(InputVal,idx,ResNN,HiddenLayers,UnknownData,SpcsList){
  nFolds <- num_species
  myFolds <- cut(seq(1, nrow(InputVal)),
                 breaks = nFolds,
                 labels=FALSE)
  result_no<-idx
  flag<-TRUE
  misclassified<-0
  for (i in 1:nFolds) {
    testObs  <- which(myFolds == i, arr.ind = TRUE)

    dfTest   <- InputVal[ testObs, ]
    dfTrain  <- InputVal[-testObs, ]

    train_X<-data.matrix(dfTrain[,-c(1)])

    train_y<-as.numeric(dfTrain$target)

    test_X<-data.matrix(dfTest[,-c(1)])
    test_y<-as.numeric(dfTest$target)


    nn.train<- neuralnetwork(X = train_X, y = train_y,
                             hidden.layers = HiddenLayers,
                             val.prop = 0,
                             optim.type = 'adam',
                             loss.type = "log",
                             batch.size = num_species-1,
                             activ.functions = "relu",
                             standardize = FALSE,
                             learn.rates = 1e-2,
                             L1=ResNN[idx,]$L1,
                             L2=ResNN[idx,]$L2,
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
    misclassified=misclassified/num_species
    predictions=predict(nn.train,UnknownData)
    predicted_pheno=cbind(SpcsList,predictions$predictions)
    weights = nn.train$Rcpp_ANN$getParams()$weights
    return(list(misclassified,predicted_pheno,weights))
}
get_zero<-function(res,dtEx,unknownData,sps){

  set<-which(res$CV==min(res$CV))
  print(set)
  for (subset in set)
  {
    if(is.null(res[subset,]$layers))
    {
      hl<-NA
    }
    else
    {

      hl<-res[subset,]$layers
    }
    miss<-fit.model(dtEx,subset,res,hl,unknownData,sps)

  }
  return (list(miss[1],hl,miss[2],miss[3]))
}
get_gene_architecture<-function(tid,pst,pfi,ptr,utf,pre,hla,numsp){
  if(hla>0)
  {
    print(paste0('This might take a while as we will try all different permutations of ',hidden_layers,' hidden layer architectures with each layer having nodes ranging [1 to number of species]'))
  }
    dataset<-get_data(tid, pst,pfi,utf,ptr)
    x<-get_zero(train_NN(dataset[[1]],hla,numsp),dataset[[1]],dataset[[5]],dataset[[2]])
    if(x[[1]]==0)
    {
      print(paste0('found architecture with CV error = 0 with architecture: ',ifelse(is.na(x[[2]]),0,x[[2]])))
      if (!dir.exists(pre)) 
      {
        dir.create(pre)
      }
      write.csv(x[[3]],paste0(pre,'predictions.csv'))
      weights<-find_weights(x[[4]],utf,dataset[[3]],dataset[[4]])
      writeLines(as.character(weights),paste0(pre,'weights.txt'))
      return(TRUE)
    }
    else if (x[[1]]<mincv) 
    {
      mincv=x[[1]]
      print(paste0('found architecture with CV error =',x[[1]],' with architecture: ',ifelse(is.na(x[[2]]),0,x[[2]])))
      if (!dir.exists(pre)) 
      {
        dir.create(pre)
      }
      write.csv(x[[3]],paste0(pre,'predictions.csv'))
      weights<-find_weights(x[[4]],utf,dataset[[3]],dataset[[4]])
      writeLines(as.character(weights),paste0(pre,'weights.txt'))
      return(FALSE)
    }
  }


mincv=1
for (hl in 0:3) {
  if(get_gene_architecture(transcript_id,path_status,path_file,path_tree,use_tree_features,path_to_results,hl,num_species))
  {break}
}
