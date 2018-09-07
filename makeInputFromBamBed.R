

getinput <- function(filepath, ip_bam, peak_bed, outfilepath) {
  
  library(exomePeak)
  library(Guitar)
  
  source(paste(filepath, "makePeak2SB.R", sep = "/"))
  source(paste(filepath, "makePeak2Motif.R", sep = "/"))
  source(paste(filepath, "bamReadsCount.R", sep = "/"))
  source(paste(filepath, "getinputFun.R", sep = "/"))
  load(paste(filepath, "psuedoGene.RData", sep = "/"))
  
  
  ## get merippp
  meripp <- BED12toGRangesList(peak_bed)
  id <- countOverlaps(meripp, tx_genes, type = "within")
  meripp <- meripp[id!=0]
  
  motif <- c("GGACA", "GGACC", "GGACT", "AGACA", "AGACC", "AGACT", 
             "GAACA", "GAACC", "GAACT", "AAACA", "AAACC", "AAACT",
             "TGACA", "TGACC", "TGACT", "TAACA", "TAACC", "TAACT")
  
  input <- .makeInput(meripp, tx_genes, motif, ip_bam,
                      extend = 50, min_reads = 5)
  
  save(input, file = paste(outfilepath, "input.RData", sep = "/"))
  
  return(input)
}
