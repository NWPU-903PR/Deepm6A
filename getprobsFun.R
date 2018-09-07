

## select model --------------------------------------------------------------

.modelselect <- function(filepath, epochid) {
  
  mod_list <- list.files(path = filepath,
                         pattern = "hdf5$", full.names = F, recursive = F)
  
  mod_list <- strsplit(mod_list, "-")
  repid <- unlist(lapply(mod_list, "[[", 3))
  
  meth <- unlist(lapply(mod_list, "[[", 1))
  negsub <- unlist(lapply(mod_list, "[[", 2))
  repsub <- unlist(lapply(mod_list, "[[", 4))
  repsub <- strsplit(repsub, "_")
  
  los <- as.numeric(unlist(lapply(repsub, "[[", 2)))
  acc <- as.numeric(unlist(lapply(repsub, "[[", 3)))
  
  epochs <- unlist(lapply(repsub, "[[", 4))
  epochs <- lapply(strsplit(epochs, "[.]"), "[[", 1)
  epochs <- as.numeric(epochs)
  
  id <- (epochs == epochid)
  
  metr <- acc - 0.8*los + (1:length(acc))*10^-6
  
  los_ind <- paste(meth, negsub)
  acc_max <- tapply(metr[id], los_ind[id], max)
  model_id <- match(acc_max, metr)
  names(model_id) <- names(acc_max)
  
  mod_list <- list.files(path = filepath,
                         pattern = "hdf5$", full.names = T, recursive = F)
  mod_list <- mod_list[model_id]
  
  seq_model <- mod_list[1:(length(model_id)/2)]
  seqreads_model <- mod_list[(length(model_id)/2 + 1):length(model_id)]
  
  return(seqreads_model)
  
}



## sample input ----------------------------------------------------------------------------------

.sampleinput <- function(input) {
  
  library(keras)
  
  ## get cnn input
  cnn_input <- input$cnn_input
  cnn_x_test <- lapply(cnn_input, t)
  cnn_x_test <- array_reshape(unlist(cnn_x_test), 
                              c(length(cnn_x_test), ncol(cnn_x_test[[1]]), nrow(cnn_x_test[[1]]), 1))
  
  return(cnn_x_test)
}


## get model probs --------------------------------------------------------------------------------

.modelpredict <- function(model, test) {
  
  library(keras)
  model <- load_model_hdf5(model, compile = F)
  
  ## compile model
  model %>% compile(
    loss = loss_categorical_crossentropy,
    optimizer = optimizer_adadelta(),
    metrics = c('accuracy')
  )
  
  ## predict
  y <- model %>% predict(test)
  probs <- y[,2]
  
  return(probs)
}


.combineprobs <- function(x) {
  n <- length(x[[1]])
  y <- matrix(unlist(x), nrow = n)
  y <- rowMeans(y)
  y <- data.frame(predict = y)
  return(y)
}

## get predict result -----------------------------------------------------------------------------

getpredict <- function(input, filepath, outputfile) {
  
  library(exomePeak)
  
  mod_list <- .modelselect(filepath, 30)
  print("get input ...")
  cnn_input <- .sampleinput(input)
  
  print("get probs ...")
  probs <- lapply(as.list(mod_list), .modelpredict, cnn_input)
  pred_prob <- .combineprobs(probs)
  
  peakloci <- input$motif_A_loci
  
  meripp <- input$meripp
  anno <- mcols(meripp)
  gene <- anno$V4
  peakname <- gene[input$motif_ind]
  
  print("get motifs ...")
  motif <- input$motif_seq
  motif <- lapply(motif, "[", 49:53)
  
  motif <- DNAStringSet(motif)
  motif <- as.character(motif)
  
  print("get results ...")
  xls <- data.frame(chr = as.character(seqnames(peakloci)),
                    chromStart = as.numeric(start(peakloci)) - 1,
                    chromEnd = as.numeric(end(peakloci)),
                    name = peakname,
                    score = pred_prob$predict,
                    strand = as.character(strand(peakloci)),
                    motif = motif)
  
  # id <- unique(peakloci)
  # ind <- findOverlaps(peakloci, id)
  # 
  # probid <- pred_prob$predict + (1:length(pred_prob$predict)) * 10^-9
  # probs <- tapply(probid, subjectHits(ind), max)
  # ind <- match(probs, probid)
  # 
  # xls <- xls[ind,]
  # peakloci <- peakloci[ind]
  
  write.table(xls, paste(outputfile, "CandidateSingleBasePeak.xls", sep = "/"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  return(xls)
}
















