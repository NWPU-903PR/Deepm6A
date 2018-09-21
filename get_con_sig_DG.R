

## run rwr ----------------------------------------------------------------------------

sigDG_test <- function(filepath, n) {
  
  setwd(filepath)
  load("F:/MeRIPSinglebase/review1/result/Pro_Dis.RData")
  load("ag0.RData")
  source("F:/MeRIPSinglebase/metdb2/process/rwrhFun.R")
  mmpath <- "rwrmm.RData"
  
  x <- read.table("F:/MeRIPSinglebase/review1/result/gene_id1_new.txt")
  siggene <- as.character(x$V1)
  siggene <- siggene[is.element(siggene, id)]
  
  sigdis <- sigdis[is.element(sigdis$gene, id),]
  
  re <- rwrh(siggene, mmpath, sigdis$doid)
  
  re0 <- list()
  filepath <- "disease/randomm"
  for(i in 1:100) {
    mmpath <- paste(filepath, "/rwrmm", i, ".RData", sep = "")
    re0[[i]] <- rwrhr(siggene, mmpath, sigdis$doid)
    print(i)
  }
  
  randomdis <- lapply(re0, "[[", 2)
  randomdis <- lapply(randomdis, "[", 1:n)
  randomdis <- unlist(randomdis)
  
  topdis <- re$topdis[1:n]
  pval <- lapply(as.list(topdis), function(x, y){sum(is.element(y,x))/100}, randomdis)
  
  re$topdis[1:n][unlist(pval) < 0.05]
  
  topdis <- re$topdis[1:n][unlist(pval) < 0.05]
  
  topdisgene <- sigdis[is.element(sigdis$disname, topdis),]
  topdisgene <- tapply(as.character(topdisgene$gene), as.character(topdisgene$disname), c)
  
  topdisgene <- lapply(topdisgene, function(x){paste(x, collapse = ",")})
  
  topdisgene <- data.frame(dis = names(topdisgene), gene = unlist(topdisgene))
  
  topdisgene <- data.frame(dis = topdis, gene = topdisgene$gene[match(topdis, topdisgene$dis)])
  
  write.table(topdisgene, file = paste("disease/top_", n, "_disgene10.xls", sep = ""), 
              quote = F, col.names = T, row.names = F, sep = "\t")
  
  return(topdisgene)
}


## get con sig dis ----------------------------------------------------------------------------------------------

n <- 10

filepath <- "F:/MeRIPSinglebase/review1/data/random/biogrid"
biogrid <- sigDG_test(filepath, n)

filepath <- "F:/MeRIPSinglebase/review1/data/random/iRef"
iRef <- sigDG_test(filepath, n)

filepath <- "F:/MeRIPSinglebase/review1/data/random/hint"
hint <- sigDG_test(filepath, n)

filepath <- "F:/MeRIPSinglebase/review1/data/random/multinet"
multinet <- sigDG_test(filepath, n)

id <- c(as.character(biogrid$dis), as.character(iRef$dis), 
        as.character(hint$dis), as.character(multinet$dis))

ind <- tapply(id, id, length)

condis <- names(ind)[ind == 4]

condis <- biogrid[is.element(biogrid$dis, condis),]

write.table(condis, file = paste("F:/MeRIPSinglebase/review1/result/top_", n, "_disgene.xls", sep = ""), 
            quote = F, col.names = T, row.names = F, sep = "\t")







