packages<-c('ANN2','parallel','BiocManager','tidyr')
for (package in packages){
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = 'http://cran.us.r-project.org')
  }
}
library(ANN2)
library(parallel)
library(tidyr)
if(!requireNamespace('Biostrings',quietly = TRUE)){
BiocManager::install("Biostrings")}
library(Biostrings)

path_source <- commandArgs(trailingOnly = TRUE)[1]
use_dendrogram_features<-as.logical(commandArgs(trailingOnly = TRUE)[2])
path_status<-commandArgs(trailingOnly = TRUE)[3]
path_file<-commandArgs(trailingOnly = TRUE)[4]
path_phylo<-commandArgs(trailingOnly=TRUE)[5]

path_tree<-'input/tree-features.csv'
#implement paths
if(use_dendrogram_features)
{
if(length(commandArgs(trailingOnly = TRUE)) == 4)
  {
    stop("Path to Phylogeny not found")
  }
  if(length(commandArgs(trailingOnly = TRUE)) > 4)
  {
    if(path_phylo=='default')
    {
      path_tree<-'input/tree-features.csv'
    }
    else
    {
    system(paste("python3", "extract_tree_feats.py", shQuote(readLines(path_phylo, warn = FALSE))))
    path_tree<-'input/custom-tree-features.csv'
    }
  }
}



hidden_layers<-0
mincv<-0

loading_chars <- c("|", "/", "-", "\\")


path_order<-'input/order.txt'
path_to_results<-"results/associated.txt"



path_tree<-paste0(path_source,path_tree)
path_file<-paste0(path_source,path_file)
path_to_results<-paste0(path_source,path_to_results)
path_status<-paste0(path_source,path_status)
path_order<-paste0(path_source,path_order)
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
  order<-read.table(path_order)
  merged<-merged[match(order$V1,merged$species_name),]
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
    con <- file(path_to_results, open = "a")
    writeLines(gene, con)
    close(con)
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
  clusterExport(cl, list("f", "read_fasta","loading_chars", "fit.model", "ids", 
                         "get_azero", "extract_sequence", 
                         "path_file", "path_status", "get_data", "path_to_results", 
                         "path_tree", "hidden_layers", "use_dendrogram_features", 
                         "train_ANN", "num_species","mincv","path_order"))
  clusterApply(cl, gene_list, f)
  stopCluster(cl)
  cat("Parallel processing completed.\n")
}

get_unique_transcript_ids <- function(fasta_file) {
  unique_ids <- character()
  con <- file(fasta_file, "r")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (startsWith(line, ">")) {
      transcript_id <- strsplit(line, " ")[[1]][1]
      transcript_id <- substr(transcript_id, 2, nchar(transcript_id))
      unique_ids <- c(unique_ids, transcript_id)
    }
  }
  close(con)
  unique_ids <- unique(unique_ids)
  return(unique_ids)
}
ids<-as.list(get_unique_transcript_ids(path_file))
run_associated_genes(ids)
