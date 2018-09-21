

network_scale <- "_new"

setwd("F:/MeRIPSinglebase/review1/result")

x <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/hint/gene_id", network_scale, ".txt", sep = ""))
hint <- as.character(x$V1)

x <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/iRef/gene_id", network_scale, ".txt", sep = ""))
iRef <- as.character(x$V1)

x <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/multinet/gene_id", network_scale, ".txt", sep = ""))
multinet <- as.character(x$V1)

x <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/biogrid/gene_id", network_scale, ".txt", sep = ""))
biogrid <- as.character(x$V1)

sum(is.element(biogrid, hint))
sum(is.element(biogrid, iRef))
sum(is.element(biogrid, multinet))

sum(is.element(hint, iRef))
sum(is.element(hint, multinet))
sum(is.element(iRef, multinet))

hi <- hint[is.element(hint, iRef)]
him <- hi[is.element(hi, multinet)]
himb <- him[is.element(him, biogrid)]

## get con edge
id <- unique(c(hint, iRef, multinet, biogrid))

he <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/hint/gene_edge", network_scale, ".txt", sep = ""), 
                 header = F, stringsAsFactors = F)
ie <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/iRef/gene_edge", network_scale, ".txt", sep = ""),
                 header = F, stringsAsFactors = F)
me <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/multinet/gene_edge", network_scale, ".txt", sep = ""),
                 header = F, stringsAsFactors = F)
be <- read.table(paste("F:/MeRIPSinglebase/review1/data/random/biogrid/gene_edge", network_scale, ".txt", sep = ""),
                 header = F, stringsAsFactors = F)

himbe <- rbind(he, ie, me, be)

inde <- cbind(match(himbe$V1, id), match(himbe$V2, id))

ag <- matrix(0, length(id), length(id))
for(i in 1:nrow(inde)) {
  ag[inde[i,1], inde[i,2]] <- ag[inde[i,1], inde[i,2]] + 1
  ag[inde[i,2], inde[i,1]] <- ag[inde[i,2], inde[i,1]] + 1
}
colnames(ag) <- id


## get 4 consensus -----------------------------------------------------------------------------------------------------------

ag0 <- ag
ag0[ag0<4] <- 0
ind <- colSums(ag0) != 0

ag0 <- ag0[ind, ind]
id0 <- colnames(ag0)

w <- ag0
library("igraph")
g <- graph_from_adjacency_matrix(w)
clu <- components(g, mode = "strong")
max(clu$csize)

cg <- clu$membership
cg <- tapply(names(cg), cg, c)

id <- unlist(lapply(cg, length))
cg <- cg[id > 1]

summary(id)


ag0[lower.tri(ag0)] <- 0
ed <- which(ag0 != 0, arr.ind = T)
gene_edge <- cbind(id0[ed[,1]], id0[ed[,2]])
gene_id <- unique(c(id0[ed[,1]], id0[ed[,2]]))

write.table(gene_edge, file = paste("gene_edge", network_scale, ".txt", sep = ""), 
            quote = F, col.names = F, row.names = F)
write.table(gene_id, file = paste("gene_id", network_scale, ".txt", sep = ""), 
            quote = F, col.names = F, row.names = F)

filepath <- paste("sig_gene", network_scale, sep = "")
dir.create(filepath)
for(i in 1:length(cg)) {
  id <- cg[[i]]
  ag0 <- w
  ag0[lower.tri(ag0)] <- 0
  id <- is.element(colnames(ag0), id)
  ag0 <- ag0[id,id]
  ed <- which(ag0 != 0, arr.ind = T)
  gene_edge <- cbind(colnames(ag0)[ed[,1]], colnames(ag0)[ed[,2]])
  gene_id <- unique(c(colnames(ag0)[ed[,1]], colnames(ag0)[ed[,2]]))
  
  write.table(gene_edge, file = paste(filepath, "/gene_edge", length(gene_id), "_", nrow(gene_edge), "_", i, ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
  write.table(gene_id, file = paste(filepath, "/gene_id", length(gene_id), "_", i,  ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
}


## get 1 network -----------------------------------------------------------------------------------------------------------

ag0 <- ag
ag0[ag0<4] <- 0
ind <- colSums(ag0) != 0

ag0 <- ag0[ind, ind]
id0 <- colnames(ag0)

w <- ag0
library("igraph")
g <- graph_from_adjacency_matrix(w)
clu <- components(g, mode = "strong")
max(clu$csize)

cg <- clu$membership
cg <- tapply(names(cg), cg, c)

id <- unlist(lapply(cg, length))
cg <- cg[id > 1]

summary(id)


ag0[lower.tri(ag0)] <- 0
ed <- which(ag0 != 0, arr.ind = T)
gene_edge <- cbind(id0[ed[,1]], id0[ed[,2]])
gene_id <- unique(c(id0[ed[,1]], id0[ed[,2]]))

## edge n ---------------------------------------------------------------------------------------------

n <- 2

ag0 <- ag
ag0[ag0<n] <- 0
ind <- colSums(ag0) != 0

ag0 <- ag0[ind, ind]
id0 <- colnames(ag0)

w <- ag0
library("igraph")
g <- graph_from_adjacency_matrix(w)
clu <- components(g, mode = "strong")
max(clu$csize)

cg <- clu$membership
cg <- tapply(names(cg), cg, c)

id <- unlist(lapply(cg, length))
cg <- cg[id > 1]

summary(id)

filtercon <- function(x, gene_id) {
  x <- x[is.element(x, gene_id)]
  return(length(x))
}

id <- unlist(lapply(cg, filtercon, gene_id))
cg <- cg[id > 1]
id <- unique(unlist(cg))

ag0 <- ag0[is.element(colnames(ag0), id), is.element(colnames(ag0), id)]
id0 <- colnames(ag0)

ag0[lower.tri(ag0)] <- 0

ed <- which(ag0 != 0, arr.ind = T)
edw <- ag0[which(ag0 != 0, arr.ind = T)]
gene_edge1 <- cbind(id0[ed[,1]], id0[ed[,2]], edw)
gene_id1 <- unique(c(id0[ed[,1]], id0[ed[,2]]))

write.table(gene_edge1, file = paste("gene_edge", n, network_scale, ".txt", sep = ""), 
            quote = F, col.names = F, row.names = F)
write.table(gene_id1, file = paste("gene_id", n, network_scale, ".txt", sep = ""), 
            quote = F, col.names = F, row.names = F)


filepath <- paste("sig_gene", n, network_scale, sep = "")
dir.create(filepath)
for(i in 1:length(cg)) {
  
  id <- cg[[i]]
  ag0 <- w
  ag0[lower.tri(ag0)] <- 0
  id <- is.element(colnames(ag0), id)
  ag0 <- ag0[id,id]
  ed <- which(ag0 != 0, arr.ind = T)
  edw <- ag0[which(ag0 != 0, arr.ind = T)]
  gene_edge <- cbind(colnames(ag0)[ed[,1]], colnames(ag0)[ed[,2]], edw)
  gene_id <- unique(c(colnames(ag0)[ed[,1]], colnames(ag0)[ed[,2]]))
  
  write.table(gene_edge, file = paste(filepath, "/gene_edge", length(gene_id), "_", nrow(gene_edge), "_", i, ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
  write.table(gene_id, file = paste(filepath, "/gene_id", length(gene_id), "_", i,  ".txt", sep = ""), 
              quote = F, col.names = F, row.names = F)
}


























