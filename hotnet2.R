

network_scale <- 15

setwd("F:/MeRIPSinglebase/review1/result")
load("F:/hotnet2/data/ag.RData")
load("F:/MeRIPSinglebase/review1/data/geneCorrelation10.RData")

fc <- abs(corr_m$corr)
gene <- corr_m$name

fc0 <- fc[is.element(gene, id)]
gene <- as.character(gene[is.element(gene, id)])
names(fc0) <- gene

delt_random <- rep(0, 100)
delt_m <- numeric()
csize_m <- numeric()

library("igraph")
for(i in 1:100){
  tic <- Sys.time()
  load(paste("F:/MeRIPSinglebase/metdb2/data/rank_shanmu/rank", i, ".RData", sep=""))
  
  id <- colnames(rank)
  idfc <- fc0[match(id, gene)]
  w0 <- t(t(rank)*as.numeric(idfc))
  
  deltv <- quantile(w0, seq(0.999, 0.9995, by = 0.00001))
  csizev <- rep(0, length(deltv))
  for(h in 1:length(deltv)) {
    delt <- deltv[h]
    w <- w0
    w[w>delt] <- 1
    w[w<=delt] <- 0
    
    # w[(w-t(w)) != 0] <- 0
    diag(w) <- 0
    
    g <- graph_from_adjacency_matrix(w)
    clu <- components(g, mode = "strong")
    csizev[h] <- max(clu$csize)
  }
  delt_m <- rbind(delt_m, deltv)
  csize_m <- rbind(csize_m, csizev)
  
  delt_random[i] <- as.numeric(deltv[which.min(abs(csizev-network_scale))])
  print(i)
  print(Sys.time() - tic)
}

save(delt_random, file = paste("F:/MeRIPSinglebase/metdb2/data/network_delt", network_scale, ".RData", sep = ""))
save(delt_m, csize_m, file = "F:/MeRIPSinglebase/metdb2/data/network_deltm.RData")

## get network ------------------------------------------------------------------------------------

network_scale <- 15

load("F:/hotnet2/data/ag.RData")
load("F:/hotnet2/data/restart0.5.RData")
load("F:/MeRIPSinglebase/metdb2/data/geneCorrelation.RData")
load("F:/MeRIPSinglebase/metdb2/data/network_deltm.RData")

fc <- abs(corr_m$corr)
gene <- corr_m$name

fc0 <- fc[is.element(gene, id)]
gene <- as.character(gene[is.element(gene, id)])
names(fc0) <- gene


id0 <- rep(0, length(id))
id0[is.element(id, gene)] <- fc0[match(id[is.element(id, gene)], gene)]

rank <- t(t(rank)*id0)
row.names(rank) <- id
colnames(rank) <- id

w0 <- rank[id0 != 0, id0 != 0]

w <- w0

delt_random <- rep(0, 100)
for(i in 1:100){
  deltv <- delt_m[i,]
  csizev <- csize_m[i,]
  
  delt_random[i] <- as.numeric(deltv[which.min(abs(csizev-network_scale))])
}

delt <- median(delt_random)

w[w>delt] <- 1
w[w<=delt] <- 0
diag(w) <- 0

# w[(w-t(w)) != 0] <- 0
sum(w)/2

library("igraph")
g <- graph_from_adjacency_matrix(w)
clu <- components(g, mode = "strong")
max(clu$csize)

cg <- clu$membership
cg <- tapply(names(cg), cg, c)

id <- unlist(lapply(cg, length))
cg <- cg[id > 1]

save(cg, file = paste("F:/MeRIPSinglebase/metdb2/data/sig_gene_clus", network_scale, ".RData", sep = ""))

id <- unlist(cg)
ag <- w
ag[(ag-t(ag)) == -1] <- 1
ag[lower.tri(ag)] <- 0
id <- is.element(row.names(ag), id)
ag <- ag[id,id]
ed <- which(ag == 1, arr.ind = T)
gene_edge <- cbind(row.names(ag)[ed[,1]], colnames(ag)[ed[,2]])
gene_id <- unique(c(row.names(ag)[ed[,1]], colnames(ag)[ed[,2]]))

write.table(gene_edge, file = paste("F:/MeRIPSinglebase/metdb2/result/gene_edge", network_scale, ".txt", sep = ""), 
            quote = F, col.names = F, row.names = F)
write.table(gene_id, file = paste("F:/MeRIPSinglebase/metdb2/result/gene_id", network_scale, ".txt", sep = ""), 
            quote = F, col.names = F, row.names = F)

filepath <- paste("F:/MeRIPSinglebase/metdb2/result/sig_gene", network_scale, sep = "")
dir.create(filepath)
for(i in 1:length(cg)) {
  id <- cg[[i]]
  ag <- w
  ag[(ag-t(ag)) == -1] <- 1
  ag[lower.tri(ag)] <- 0
  id <- is.element(row.names(ag), id)
  ag <- ag[id,id]
  ed <- which(ag == 1, arr.ind = T)
  gene_edge <- cbind(row.names(ag)[ed[,1]], colnames(ag)[ed[,2]])
  gene_id <- unique(c(row.names(ag)[ed[,1]], colnames(ag)[ed[,2]]))
  
  write.table(gene_edge, file = paste(filepath, "/gene_edge", nrow(gene_edge), "_", i, ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
  write.table(gene_id, file = paste(filepath, "/gene_id", length(gene_id), "_", i,  ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
}


