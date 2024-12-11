get_num_species<-function(path_status){
  phylogeny<-read.table(path_status,header=1)
  colnames(phylogeny)<-c('species_name','target')
  return(sum(!is.na(phylogeny$target)))
}
check_all_dashes <- function(df) {
  df<-na.omit(df)
  target_1_rows <- df[df$target == 1, ]
  
  if (nrow(target_1_rows) == 0) {
    return(FALSE)
  }
  all_dashes <- apply(target_1_rows[, 4:ncol(target_1_rows)], 1, function(row) all(row == "-"))
  return(any(all_dashes))
}
read_fasta <- function(file_path) {

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
generate_hidden_layer_tuples <- function(n, x) {
  if (n == 0) {
    return(NA)
  } else {
    return(permutations(x, n))
  }}
get_data<-function(trid,phenoStatus,SeqFile,TreeFlag,TreeFile,po){
  sequence <- extract_sequence(read_fasta(SeqFile), trid)
  sequence_list_chars <- lapply(sequence[[2]], function(seq) unlist(strsplit(as.character(seq), "")))
  sequence_df <- as.data.frame(do.call(rbind, sequence_list_chars))
  sequence_df<-cbind(sequence[[1]],sequence_df)
  sequence_df <- separate(sequence_df, 'sequence[[1]]', into = c("trid", "species_name"), sep = " ")
  phylogeny<-read.table(phenoStatus,header=1)
  colnames(phylogeny)<-c('species_name','target')
  merged<-merge(phylogeny,sequence_df,by='species_name')
  merged <- merged[, c(3, 1, 2, seq(4, ncol(merged)))]
  order<-read.table(po)
  merged<-merged[match(order$V1,merged$species_name),]
  if(check_all_dashes(merged)) {print('x')}
  data<-merged[,-c(1,2,3)]
  data = data.frame(lapply(data, function(x){gsub("-", 0, x)}))
  data = data.frame(lapply(data, function(x){gsub("A|T|G|C|a|t|g|c|N|n", 1, x)}))
  unique_cols <- sapply(data, function(x) length(unique(x))) == 1
  data<- data[, !unique_cols]
  data <- data.frame(lapply(data, function(x) as.numeric(as.character(x))))
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
train_NNLOOCV <- function(tid, InputNN, num_hidden_layer, numSpcs, learn_rate = 1e-2, n_epochs = 500){

  nFolds <- numSpcs
  if (nFolds <= 1) stop("numSpcs should be greater than 1")
  
  myFolds <- cut(seq(1, nrow(InputNN)), breaks = nFolds, labels = FALSE)
  l1Vals <- 10^seq(-4, 3, length.out = 100)
  alphas <- seq(0, 1, 0.05)
  hl <- generate_hidden_layer_tuples(num_hidden_layer, numSpcs - 1)
  
  message(paste0('Gene ', tid, ' running: '))
  
  results <- list()
  total_iterations <- length(l1Vals) * length(alphas) * ifelse(is.null(dim(hl)), 1, nrow(hl))
  
  # Add progress bar
  withProgress(message = "Running training", value = 0, {
    iteration <- 0
    
    for (i in seq_along(l1Vals)) {
      for (j in seq_along(alphas)) {
        l <- l1Vals[i]
        alpha <- alphas[j]
        l1 <- alpha * l
        l2 <- (1 - alpha) * l
        
        hidden_layers <- if (is.null(dim(hl))) list(NA) else asplit(hl, 1)
        
        for (hidden_layer in hidden_layers) {
          misclassified <- 0
          valid_run <- TRUE
          
          for (fold in 1:nFolds) {
            testObs <- which(myFolds == fold, arr.ind = TRUE)
            dfTest <- InputNN[testObs, ]
            dfTrain <- InputNN[-testObs, ]
            
            train_X <- data.matrix(dfTrain[, -1])
            train_y <- as.numeric(dfTrain$target)
            test_X <- data.matrix(dfTest[, -1])
            test_y <- as.numeric(dfTest$target)
            
            nn.train <- tryCatch({
              neuralnetwork(X = train_X, y = train_y,
                            hidden.layers = hidden_layer,
                            val.prop = 0,
                            optim.type = 'adam',
                            loss.type = "log",
                            batch.size = numSpcs - 1,
                            standardize = FALSE,
                            learn.rates = learn_rate,
                            L1 = l1,
                            L2 = l2,
                            n.epochs = n_epochs,
                            verbose = FALSE,
                            random.seed = 1)
            }, error = function(e) {
              message("Error in neural network training: ", e$message)
              NULL
            })
            
            if (is.null(nn.train)) {
              valid_run <- FALSE
              break
            }
            
            logit.prob <- predict(nn.train, newdata = test_X)
            pred.class <- ifelse(logit.prob$probabilities[1] > logit.prob$probabilities[2], 0, 1)
            
            if (is.na(pred.class)) {
              valid_run <- FALSE
              break
            }
            
            if (pred.class != test_y) {
              misclassified <- misclassified + 1
            }
          }
          
          if (valid_run) {
            results <- c(results, list(list(L1 = l1, L2 = l2, alpha = alpha,
                                            layers = if (is.na(hidden_layer[1])) NA else paste(hidden_layer, collapse = ","),
                                            CV = misclassified / numSpcs)))
          }
          
          iteration <- iteration + 1
          incProgress(1 / total_iterations, detail = paste("Iteration", iteration, "of", total_iterations))
        }
      }
    }
  })
  
  message(paste0(tid, ' completed'))
  res <- do.call(rbind, lapply(results, data.frame))
  return(res)
}
train_NN5CV <- function(tid, InputNN, num_hidden_layer, numSpcs, learn_rate = 1e-2, n_epochs = 500) {
  nFolds <- 5
  myFolds <- cut(seq(1, nrow(InputNN)), breaks = nFolds, labels = FALSE)
  l1Vals <- 10^seq(-4, 2, length.out = 100)
  alphas <- seq(0, 1, 0.05)
  hl <- generate_hidden_layer_tuples(num_hidden_layer, numSpcs - 1)
  
  message(paste0('Gene ', tid, ' running: '))
  
  results <- list()
  total_iterations <- length(l1Vals) * length(alphas) * ifelse(is.null(dim(hl)), 1, nrow(hl))
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
  iteration <- 0  # Counter for progress bar
  
  for (i in seq_along(l1Vals)) {
    for (j in seq_along(alphas)) {
      l <- l1Vals[i]
      alpha <- alphas[j]
      l1 <- alpha * l
      l2 <- (1 - alpha) * l
      
      hidden_layers <- if(is.null(dim(hl))) list(NA) else asplit(hl, 1)
      
      for (hidden_layer in hidden_layers) {
        misclassified <- 0
        valid_run <- TRUE
        
        for (fold in 1:nFolds) {
          testObs <- which(myFolds == fold, arr.ind = TRUE)
          dfTest <- InputNN[testObs, ]
          dfTrain <- InputNN[-testObs, ]
          
          train_X <- data.matrix(dfTrain[, -1])
          train_y <- as.numeric(dfTrain$target)
          test_X <- data.matrix(dfTest[, -1])
          test_y <- as.numeric(dfTest$target)
          batch_size <- nrow(train_X)
          
          nn.train <- tryCatch({
            neuralnetwork(X = train_X, y = train_y,
                          hidden.layers = hidden_layer,
                          val.prop = 0,
                          optim.type = 'adam',
                          loss.type = "log",
                          batch.size = batch_size,
                          standardize = FALSE,
                          learn.rates = learn_rate,
                          L1 = l1,
                          L2 = l2,
                          n.epochs = n_epochs,
                          verbose = FALSE,
                          random.seed = 1)
          }, error = function(e) {
            message("Error in neural network training: ", e$message)
            NULL
          })
          
          if (is.null(nn.train)) {
            valid_run <- FALSE
            break
          }
          logit.prob <- predict(nn.train, newdata = test_X)
          pred.class <- apply(logit.prob$probabilities, 1, function(row) {
            if (all(is.nan(row))) {
              return(NA)  
            } else {
              return(which.max(unlist(row)) - 1)
            }
          })
          misclassified <- misclassified + sum(pred.class != test_y)
        }
        
        if (valid_run) {
          results <- c(results, list(list(L1 = l1, L2 = l2, alpha = alpha,
                                          layers = if(is.na(hidden_layer[1])) NA else paste(hidden_layer, collapse = ","),
                                          CV = misclassified / nFolds)))
        }
        
        iteration <- iteration + 1
        setTxtProgressBar(pb, iteration)
      }
    }
  }
  
  close(pb)
  cat("\n")  # New line after progress bar
  
  message(paste0(tid, ' completed'))
  
  res <- do.call(rbind, lapply(results, data.frame))
  return(res)
}
fit.model <- function(InputVal, idx, ResNN, HiddenLayers, UnknownData, SpcsList,num_species) {
  nFolds <- num_species
  myFolds <- cut(seq_len(nrow(InputVal)), breaks = nFolds, labels = FALSE)
  misclassified <- 0
  
  for (i in seq_len(nFolds)) {
    # Split data into train and test sets
    testObs <- which(myFolds == i, arr.ind = TRUE)
    dfTest <- InputVal[testObs, ]
    dfTrain <- InputVal[-testObs, ]
    
    train_X <- as.matrix(dfTrain[, -1])  # Assuming target is the first column
    train_y <- as.numeric(dfTrain$target)
    test_X <- as.matrix(dfTest[, -1])
    test_y <- as.numeric(dfTest$target)
    nn.train <- tryCatch({
      neuralnetwork(X = train_X, y = train_y,
                    hidden.layers = HiddenLayers,
                    val.prop = 0,
                    optim.type = 'adam',
                    loss.type = "log",
                    batch.size = num_species - 1,
                    activ.functions = "relu",
                    standardize = FALSE,
                    learn.rates = 1e-2,
                    L1 = ResNN[idx, "L1"],
                    L2 = ResNN[idx, "L2"],
                    n.epochs = 500,
                    verbose = FALSE,
                    random.seed = 1)
    }, error = function(e) {
      message("Error in neural network training: ", e$message)
      return(NULL)
    })
    
    if (is.null(nn.train)) {
      return(list(error = "Failed to train neural network", misclassification_rate = 1))
    }
    
    # Make predictions
    logit.prob <- predict(nn.train, newdata = test_X)
    pred.class <- ifelse(logit.prob$probabilities[1] > logit.prob$probabilities[2], 0, 1)
    
    if (any(is.na(pred.class))) {
      warning("NaN encountered in predictions")
      return(list(error = "NaN in predictions", misclassification_rate = 1))
    }
    
    misclassified <- misclassified + sum(pred.class != test_y)
  }
  
  misclassification_rate <- misclassified / nrow(InputVal)
  predictions <- predict(nn.train, UnknownData)
  predicted_pheno <- cbind(SpcsList, predictions$predictions)
  weights <- nn.train$Rcpp_ANN$getParams()$weights

  return(list(
    misclassification_rate = misclassification_rate,
    predicted_phenotypes = predicted_pheno[, -1],
    weights = weights
  ))
}
get_zero<-function(res,dtEx,unknownData,sps,numsp){

  set<-which(res$CV==min(res$CV))
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
    miss<-fit.model(dtEx,subset,res,hl,unknownData,sps,numsp)

  }
  return (list(miss[1],hl,miss[2],miss[3]))
}
find_weights <- function(params, tf, V, d) {
  weights <- params[[1]][[1]][1, ]
  weights <- abs(weights)
  if (length(params) > 1) {
    for (i in 2:length(params))
    {
      weights <- t(weights) %*% t(params[[i]])
    }
  }
  
  if (tf) {
    pca_weights <- ((V[, 1:d]) %*% as.matrix(weights[1:d]))
  } else {
    pca_weights <- ((V[, 1:d])) %*% as.matrix(weights)
  }
  return(pca_weights)
}
plot_manhattan <- function(weights) {
  betas <- as.numeric(weights)
  dataf <- data.frame(betas)
  p <- length(betas)
  d2 <- mahalanobis(x = as.matrix(betas), center = mean(betas), cov = var(betas))
  f.stat <- (p - 1) * d2 / (p - 1)
  neglog10pval <- -pf(f.stat, 1, p - 1, lower.tail = FALSE, log.p = TRUE)
  dataf$P <- 10^-(neglog10pval)
  dataf$P <- p.adjust(dataf$P, method = "BH")
  dataf$position <- seq(1, length(neglog10pval), 1)
  dataf$BP <- (dataf$position)
  dataf$SNP <- dataf$position
  return(dataf)
}
get_gene_architecture1 <- function(tid, pst, pfi, ptr, utf, pre, hla, numsp,numcv,mincv,po){
  if(hla>0)
  {
    print(paste0('This might take a while as we will try all different permutations of ',hla,' hidden layer architectures with each layer having nodes ranging [1 to number of species]'))
  }
    tryCatch({
      dataset <- get_data(tid, pst, pfi, utf, ptr,po)
    }, error = function(e) {
      cat("Error occurred due to incorrect custom phylogeny, please make sure the tree has correct species names and tree format! ", conditionMessage(e), "\n")
      stop("Execution halted due to error.")
    })
    if(numcv==1)
    {
        x<-get_zero(train_NNLOOCV(tid,dataset[[1]],hla,numsp),dataset[[1]],dataset[[5]],dataset[[2]],numsp)
    }
    if(numcv==5)
    {
      x<-get_zero(train_NN5CV(tid,dataset[[1]],hla,numsp),dataset[[1]],dataset[[5]],dataset[[2]],numsp)
    }
   if(x[[1]]==0)
    {
      print(paste0('found architecture with CV error = 0 with architecture: ',ifelse(is.na(x[[2]]),0,x[[2]])))
      df<-data.frame(x[[3]])
      colnames(df)<-c('species','predictions')
      print(df)
      #write.csv(df,paste0(pre,'predictions',tid,'.csv'),row.names=FALSE)
      #print("predictions stored in specified directory.")
      return(TRUE )
    }
    else if (x[[1]]<mincv) 
    {
      mincv=x[[1]]
      print(paste0('found architecture with CV error =',x[[1]],' with architecture: ',ifelse(is.na(x[[2]]),0,x[[2]])))
      df<-data.frame(x[[3]])
      colnames(df)<-c('species','predictions')
      print(df)
      #write.csv(df,paste0(pre,'predictions',tid,'.csv'),row.names=FALSE)
      #print("predictions stored in specified directory.")
      return(FALSE)
    }
}
get_gene_architecture2 <- function(tid, pst, pfi, ptr, utf, pre, hla, numsp,numcv,mincv,po) {
  if (hla > 0) {
    print(paste0("This might take a while as we will try all different permutations of ", hla, " hidden layer architectures with each layer having nodes ranging [1 to number of species]"))
  }
  tryCatch(
    {
      dataset <- get_data(tid, pst, pfi, utf, ptr, po)
    },
    error = function(e) {
      cat("Error occurred due to incorrect custom phylogeny, please make sure the tree has correct species names and tree format! ", conditionMessage(e), "\n")
      stop("Execution halted due to error.")
    }
  )
  if(numcv==1)
  {
    x<-get_zero(train_NNLOOCV(tid,dataset[[1]],hla,numsp),dataset[[1]],dataset[[5]],dataset[[2]],numsp)
  }
  if(numcv==5)
  {
    x<-get_zero(train_NN5CV(tid,dataset[[1]],hla,numsp),dataset[[1]],dataset[[5]],dataset[[2]],numsp)
  }
  if (x[[1]] == 0) {
    print(paste0("found architecture with CV error = 0 with architecture: ", ifelse(is.na(x[[2]]), 0, x[[2]])))
    weights <- find_weights(x[[4]], utf, dataset[[3]], dataset[[4]])
    df <- plot_manhattan(weights)
    df2<-df[, c("position", "P")]
    colnames(df2)<-c('position','p-value')
    #write.csv(df2, paste0(pre, "/PostionalPvals.csv"), row.names = FALSE)
    df_ <- df[df$P < 0.05, ]
    print("Positions: ")
    print(df_$position)
    return(TRUE)
  } else if (x[[1]] < mincv) {
    mincv <- x[[1]]
    print(paste0("found architecture with CV error =", x[[1]], " with architecture: ", ifelse(is.na(x[[2]]), 0, x[[2]])))
    #print(x[[3]])
    weights <- find_weights(x[[4]], utf, dataset[[3]], dataset[[4]])
    df <- plot_manhattan(weights)
    df2<-df[, c("position", "P")]
    colnames(df2)<-c('position','p-value')
    #write.csv(df2, paste0(pre, "/PostionalPvals.csv"), row.names = FALSE)
    df_ <- df[df$P < 0.05, ]
    print("Positions: ")
    print(df_$position)
    return(FALSE)
  }
}
fit.model_ <-function(data,index,results,hl,num_species){
  #print(results[index,])
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
get_azero<-function(results,data,num_species,mincv){
  set<-which(results$CV<=mincv)

  for (subset in set)
  {
    hl<-results[subset,]$HL
    hl<-ifelse(hl!= 0,hl, NA)
    miss<-fit.model_(data,subset,results,hl,num_species)
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
train_ANN <- function(data, hl, num_species, gene, path_to_results,mincv) {
  nFolds <- num_species
  cv_error <- list()
  myFolds <- cut(seq(1, nrow(data)),
                 breaks = nFolds,
                 labels = FALSE)
  l1Vals <- 10^seq(-4, 3, length.out = 100)
  alphas <- seq(0, 1, 0.05)
  count <- 0
  total_steps <- length(l1Vals) * length(alphas)
  
  # Progress bar
  withProgress(message = "Training Neural Network", value = 0, {
    for (l in l1Vals) {
      for (alpha in alphas) {
        l1 <- alpha * l
        l2 <- (1 - alpha) * l
        misclassified <- 0
        flag <- TRUE
        for (i in 1:nFolds) {
          testObs <- which(myFolds == i, arr.ind = TRUE)
          dfTest <- data[testObs, ]
          dfTrain <- data[-testObs, ]
          
          train_X <- data.matrix(dfTrain[, -c(1)])
          train_y <- as.numeric(dfTrain$target)
          
          test_X <- data.matrix(dfTest[, -c(1)])
          test_y <- as.numeric(dfTest$target)
          
          nn.train <- neuralnetwork(X = train_X, y = train_y,
                                    hidden.layers = ifelse(hl != "0", as.numeric(strsplit(hl, ",")[[1]]), NA),
                                    val.prop = 0,
                                    optim.type = 'adam',
                                    loss.type = "log",
                                    batch.size = num_species - 1,
                                    standardize = FALSE,
                                    learn.rates = 0.01,
                                    L1 = l1,
                                    L2 = l2,
                                    n.epochs = 500,
                                    verbose = FALSE,
                                    random.seed = 1)
          logit.prob <- predict(nn.train, newdata = test_X)
          pred.class <- ifelse(logit.prob$probabilities[1] > logit.prob$probabilities[2], 0, 1)
          if (is.na(pred.class)) {
            flag <- FALSE
            break
          }
          tryCatch(
            expr = {
              if (pred.class != test_y) {
                misclassified <- misclassified + 1
              }
            },
            error = function(e) {}
          )
        }
        if (flag != FALSE) {
          cv_error <- append(cv_error, list(list(L1 = l1, L2 = l2, alpha = alpha, CV = (misclassified / num_species), HL = hl)))
        }
        count <- count + 1
        
        # Update progress bar
        incProgress(1 / total_steps, detail = paste("Step", count, "of", total_steps))
        flush.console()
        Sys.sleep(0.1)
      }
    }
  })
  
  if (get_azero(do.call(rbind, lapply(cv_error, data.frame)), data, num_species,mincv)) {
    print(paste0(gene, " associated!"))
    con <- file(path_to_results, open = "a")
    writeLines(gene, con)
    close(con)
  } else {
    print(paste0(gene, " not associated"))
  }
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

PredictRegions<-function(gene_list,path_status,path_file, use_tree_features, path_tree,path_order, hidden_layers, num_species,path_to_results,mincv){
  for(i in gene_list)
  {
    cat(paste0("Working on gene: ", i, "\n"))
    data <- get_data(i,path_status,path_file, use_tree_features, path_tree,path_order)[[1]]
    train_ANN(data, hidden_layers, num_species,i,path_to_results,mincv)
  }
}
PredictSpecies<-function(transcript_id,path_status,path_file,path_tree,use_tree_features,path_to_results,num_species,numcv,mincv=1,po){
  for (hl in 0:3) {
  if(get_gene_architecture1(transcript_id,path_status,path_file,path_tree,use_tree_features,path_to_results,hl,num_species,numcv,mincv,po))
  return(TRUE)
  }
}
PositionalPvals<-function(transcript_id,path_status,path_file,path_tree,use_tree_features,path_to_results,num_species,numcv,mincv=1,po){
  for (hl in 0:3) {
    if(get_gene_architecture2(transcript_id,path_status,path_file,path_tree,use_tree_features,path_to_results,hl,num_species,numcv,mincv,po))
      return(TRUE)
  }
}