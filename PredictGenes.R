library(ANN2)
library(dplyr)
library(Biostrings)
library(tidyr)
library(utils)
library(parallel)
path_source <- commandArgs(trailingOnly = TRUE)[1]
use_dendrogram_features<-as.logical(commandArgs(trailingOnly = TRUE)[2])
path_status<-commandArgs(trailingOnly = TRUE)[3]
path_file<-commandArgs(trailingOnly = TRUE)[4]
path_tree<-commandArgs(trailingOnly = TRUE)[5]
transcript_list_path<-commandArgs(trailingOnly = TRUE)[6]
hidden_layers<-commandArgs(trailingOnly = TRUE)[7]
mincv<-as.numeric(commandArgs(trailingOnly = TRUE)[8])

loading_chars <- c("|", "/", "-", "\\")





path_to_results<-"results/associated.txt"



path_transcripts<-paste0(path_source,transcript_list_path)
path_tree<-paste0(path_source,path_tree)
path_file<-paste0(path_source,path_file)
path_to_results<-paste0(path_source,path_to_results)
path_status<-paste0(path_source,path_status)

get_num_species<-function(path_status){
  phylogeny<-read.table(path_status,header=1)
  colnames(phylogeny)<-c('species_name','target')
  return(sum(!is.na(phylogeny$target)))
}
num_species<-get_num_species(path_status)



if (!file.exists(path_to_results)) {
  file.create(path_to_results)
}
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
get_data<-function(transcript_id,path_status,path_file,tree_flag,path_tree){
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
fit.model <-function(data,index,results,hl){
  print(results[index,])
  nFolds <- num_species
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
                             hidden.layers = ifelse(hl != "0", as.numeric(strsplit(hl, ",")[[1]]), NA),
                             val.prop = 0,
                             optim.type = 'adam',
                             loss.type = "log",
                             batch.size = num_species-1,
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
      return(num_species)
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
  return (misclassified/num_species)
}
get_azero<-function(results,data){
  set<-which(results$CV<=mincv)
  for (subset in set)
  {
    hl<-results[subset,]$HL
    hl<-ifelse(hl!= 0,hl, NA)
    miss<-fit.model(data,subset,results,hl)
    if(miss<=mincv)
    {
      return (TRUE)
    }
    else
    {
      return (FALSE)
    }
  }
  return (FALSE)

}
train_ANN<-function(data,hl,num_species,gene){
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
        cv_error <- append(cv_error,list(list(L1=l1,L2 = l2,alpha=alpha,CV=(misclassified/num_species),HL=hl)))
      }
      count<-count+1
      cat(sprintf("\r%s", loading_chars[count %% length(loading_chars) + 1]))
      flush.console()
      Sys.sleep(0.1)
    }

  }
  if(get_azero(do.call(rbind, lapply(cv_error, data.frame)),data))
  {
    print(paste0(gene,' associated!'))
    if(!file.exists(path_to_results))
    {
      file.create(path_to_results)
    }
    writeLines(gene,path_to_results)
  }
  else 
  {
    print(paste0(gene,' not associated'))
  }
}


f <- function(i) {
  cat(paste0("Working on gene: ", i, "\n"))
  data <- get_data(i,path_status,path_file, use_dendrogram_features, path_tree)[[1]]
  train_ANN(data, hidden_layers, num_species,i)
}
run_associated_genes<-function(gene_list){
  cat("Starting parallel processing...\n")
  cat("Please do not close the window...\n")
  cl <- makeCluster(detectCores() - 1)
  clusterEvalQ(cl, library(ANN2))
  clusterEvalQ(cl, library(Biostrings))
  clusterEvalQ(cl, library(tidyr))
  clusterExport(cl, list("f", "read_fasta","loading_chars", "fit.model", "gene_list", 
                         "get_azero", "path_transcripts", "extract_sequence", 
                         "path_file", "path_status", "get_data", "path_to_results", 
                         "path_tree", "hidden_layers", "use_dendrogram_features", 
                         "train_ANN", "num_species","mincv"))
  clusterApply(cl, gene_list, f)
  stopCluster(cl)
  cat("Parallel processing completed.\n")
}

genes<-read.table(path_transcripts)
gene_list<-genes$V1
run_associated_genes(gene_list)
