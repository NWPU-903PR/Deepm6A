

setwd("F:/MeRIPSinglebase/review1/data/random/biogrid")


load("ag0.RData")
load("rank0.RData")
load("F:/MeRIPSinglebase/review1/data/geneCorrelation12.RData")
load("delt.RData")

sn <- (spval <= 0.05)
sn <- colSums(sn)

sn <- sn[order(sdelt)]
delt <- sdelt[order(sdelt)]
delt <- delt[which.max(sn)]

fc <- abs(corr_m$corr)
gene <- corr_m$name

fc0 <- fc[is.element(gene, id)]
gene <- as.character(gene[is.element(gene, id)])
names(fc0) <- gene

id0 <- fc0[match(id[is.element(id, gene)], gene)]

ind <- colnames(rank)
ind0 <- is.element(ind, gene)
rank <- rank[ind0,ind0]

rank <- t(t(rank)*id0)
row.names(rank) <- id[is.element(id, gene)]
colnames(rank) <- id[is.element(id, gene)]

w0 <- rank

w <- w0

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

summary(id)
save(cg, file = "sig_gene_clus_new.RData")

id <- unlist(cg)
ag <- w
ag[(ag-t(ag)) == -1] <- 1
ag[lower.tri(ag)] <- 0
id <- is.element(row.names(ag), id)
ag <- ag[id,id]
ed <- which(ag == 1, arr.ind = T)
gene_edge <- cbind(row.names(ag)[ed[,1]], colnames(ag)[ed[,2]])
gene_id <- unique(c(row.names(ag)[ed[,1]], colnames(ag)[ed[,2]]))

write.table(gene_edge, file = "gene_edge_new.txt", 
            quote = F, col.names = F, row.names = F)
write.table(gene_id, file = "gene_id_new.txt", 
            quote = F, col.names = F, row.names = F)

filepath <- paste("sig_gene", "new", sep = "")
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
  
  write.table(gene_edge, file = paste(filepath, "/gene_edge", length(gene_id), "_", nrow(gene_edge), "_", i, ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
  write.table(gene_id, file = paste(filepath, "/gene_id", length(gene_id), "_", i,  ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
}