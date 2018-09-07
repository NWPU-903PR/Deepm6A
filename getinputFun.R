

## make input -----------------------------------------------------------------------------------


.makeInput <- function(meripp, 
                       tx_genes, 
                       motif,
                       bamIp, 
                       extend = 200, 
                       min_reads = 5) {
  
  ## get bin of peak
  print("Cutting peak to Single Base...")
  tic <- Sys.time()
  peakSB <- .makePeak2SB(meripp, tx_genes, extend = extend)
  print(Sys.time() - tic)
  
  
  ## get motif info
  print("Searching peak for motifs...")
  tic <- Sys.time()
  peakMF <- .makePeak2Motif(meripp, tx_genes, motif, extend = extend)
  print(Sys.time() - tic)
  
  
  ## get bam reads count
  print("Counting reads for single base...")
  tic <- Sys.time()
  readsip <- .bamReadsCount(unlist(peakSB), bamIp, tx_genes)
  print(Sys.time() - tic)
  
  readsbypeak_ip <- .makeInputMat(readsip$base_reads_count, peakSB)
  
  
  ## get motif reads
  print("Making input...")
  motif_A_loci <- peakMF$motif_A_loci
  motif_ind <- peakMF$motif_ind
  motif_peak_loci <- pmapToTranscripts(motif_A_loci, meripp[motif_ind], ignore.strand = T)
  
  ind <- rbind(motif_ind, start(motif_peak_loci))
  ind <- as.vector(ind)
  id <- rep(1:length(motif_peak_loci), each = 2)
  
  readsbymotif_ip <- tapply(ind, id, .getmotifreads, readsbypeak_ip, extend)
  readsbymotif_ip <- matrix(unlist(readsbymotif_ip), 
                            nrow = (extend*2 + 1), ncol = length(readsbymotif_ip))
 
  
  ## get reads sizefactor
  sf <- readsip$total_reads_count
  
  ## normalize reads count
  ind0 <- as.character(strand(motif_A_loci)) == "-"
  readsbymotif_ip <- .normalizeReads(readsbymotif_ip, sf, min_reads, ind0)
  
  
  ## get sequance encode
  motif_seq <- peakMF$motif_seq
  tic <- Sys.time()
  seq_input <- lapply(motif_seq, .seq2vec, extend)
  print(Sys.time() - tic)
  
  
  ## make input
  cnn_input <- list()
  for(i in 1:length(seq_input)) {
    input_mat <- seq_input[[i]]
    input_mat <- t(t(input_mat)*readsbymotif_ip[,i])
    cnn_input[[i]] <- input_mat
  }
  
  re <- list(cnn_input = cnn_input,
             readsbymotif_ip = readsbymotif_ip,
             motif_seq = motif_seq,
             motif_A_loci = motif_A_loci,
             motif_peak_loci = motif_peak_loci,
             motif_ind = motif_ind,
             readsip = readsip,
             peakSB = peakSB,
             meripp = meripp)
  return(re)
}


.makeInputMat <- function(x, peakSB) {
  readsbypeak <- split(x, names(unlist(peakSB)))
  readsbypeak <- readsbypeak[names(peakSB)]
  return(readsbypeak)
}

.getmotifreads <- function(x, reads, extend) {
  reads[[x[1]]][x[2]:(x[2] + 2*extend)]
}

.normalizeReads <- function(reads, sf, min_reads, ind) {
  reads[,ind] <- reads[nrow(reads):1, ind]
  reads <- round(reads/sf*10^8)
  reads[reads < min_reads] <- min_reads
  reads <- log(reads)
  return(reads)
}

.seq2vec <- function(x, extend) {
  library(Biostrings)
  a <- matchPattern("A", x)
  a <- start(a)
  
  t <- matchPattern("T", x)
  t <- start(t)
  
  c <- matchPattern("C", x)
  c <- start(c)
  
  g <- matchPattern("G", x)
  g <- start(g)
  
  m <- matrix(0, nrow = 4, ncol = (2*extend + 1))
  m[1, a] <- 1
  m[2, t] <- 1
  m[3, c] <- 1
  m[4, g] <- 1
  return(m)
}
