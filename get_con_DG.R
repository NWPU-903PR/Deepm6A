

## run rwr ----------------------------------------------------------------------------

get_top_DG <- function(filepath, n) {
  setwd(filepath)
  load("Pro_Dis.RData")
  load("ag0.RData")
  
  source("F:/MeRIPSinglebase/metdb2/process/rwrhFun.R")
  mmpath <- paste(filepath, "/rwrmm.RData", sep = "")
  
  x <- read.table("F:/MeRIPSinglebase/review1/result/gene_id1_15.txt")
  siggene <- as.character(x$V1)
  siggene <- siggene[is.element(siggene, id)]
  
  sigdis <- sigdis[is.element(sigdis$gene, id),]
  
  re <- rwrh(siggene, mmpath, sigdis$doid)
  
  topdis <- re$topdis[1:n]
  
  sigdis <- sigdis[is.element(sigdis$disname, topdis),]
  disgenelist <- tapply(as.character(sigdis$gene), as.character(sigdis$disname), paste, collapse = ",")
  
  topdis <- data.frame(dis = names(disgenelist), gene = disgenelist)
  topgene <- re$topgene[1:n]
  
  re <- list(topdis, topgene)
  
  return(re)
}

## 4 net-----------------------------------------------------------------------------------------------------------------

n <- 30

filepath <- "F:/MeRIPSinglebase/review1/data/random/biogrid"
biogrid <- get_top_DG(filepath, n)

filepath <- "F:/MeRIPSinglebase/review1/data/random/iRef"
iRef <- get_top_DG(filepath, n)

filepath <- "F:/MeRIPSinglebase/review1/data/random/hint"
hint <- get_top_DG(filepath, n)

filepath <- "F:/MeRIPSinglebase/review1/data/random/multinet"
multinet <- get_top_DG(filepath, n)

id <- c(as.character(biogrid[[1]]$dis), as.character(iRef[[1]]$dis), 
        as.character(hint[[1]]$dis), as.character(multinet[[1]]$dis))

ind <- tapply(id, id, length)

condis <- names(ind)[ind == 4]

biogrid[[2]]
iRef[[2]]
hint[[2]]
multinet[[2]]








