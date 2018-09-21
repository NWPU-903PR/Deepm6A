

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


## omim disease gene -------------------------------------------------------------------------------

network_scale <- 15

load(paste("F:/MeRIPSinglebase/metdb2/data/sig_gene_clus", network_scale, ".RData", sep = ""))
load("F:/MeRIPSinglebase/metdb2/data/omim/ID_map.RData")
load("F:/MeRIPSinglebase/metdb2/data/T0ZSY/DDsim_allD.RData")
load("F:/MeRIPSinglebase/metdb2/data/omim/Did_Dname.RData")
omim <- read.csv("F:/MeRIPSinglebase/metdb2/data/omim/morbidmap.csv",
                 header = T, stringsAsFactors = F)

dis <- omim$Phenotype
dis <- strsplit(as.character(dis), "[(][1234]")
dis <- unlist(lapply(dis, "[", 1))

dis0 <- strsplit(dis, ",")
dis0 <- unlist(lapply(dis0, function(x){x[length(x)]}))

dis0 <- strsplit(dis0, " ")
dis0 <- unlist(lapply(dis0, function(x){x[2]}))
dis0 <- as.numeric(dis0)

dis0[dis0 < 10000] <- NA
dis0[is.na(dis0)] <- omim$MIM.Number[is.na(dis0)]

dgene <- omim$Gene.Symbols
dgene <- strsplit(dgene, ",")
dgene <- lapply(dgene, "[", 1)

dis0 <- paste("OMIM:", dis0, sep = "")

disgene <- data.frame(dis = dis, id = dis0, gene = unlist(dgene))
sigdis <- disgene[is.element(disgene$gene, unlist(cg)),]

topdisgene <- tapply(as.character(sigdis$gene), as.character(sigdis$dis), c)
topdisgene <- lapply(topdisgene, function(x){paste(unique(x), collapse = ",")})
topdisgene <- data.frame(dis = names(topdisgene), gene = unlist(topdisgene))

write.table(topdisgene, file = "F:/MeRIPSinglebase/metdb2/result/alldisgene.xls",
            quote = F, col.names = T, row.names = F, sep = "\t")

new <- data.frame(rbind(c("DOID:8761", "OMIM:606078"), 
                        c("DOID:8761", "OMIM:606077"),
                        c("DOID:9119", "OMIM:164690")))
names(new) <- names(ID_map)
ID_map <- rbind(ID_map, new)

id <- is.element(disgene$id, ID_map$OMIM)
disgene <- disgene[id,]

doid <- ID_map$DOID[match(disgene$id, ID_map$OMIM)]
disgene <- cbind(disgene, doid)

id <- is.element(disgene$doid, colnames(DDsim_allD))
disgene <- disgene[id,]

disname <- Did_Dname$Dname[match(disgene$doid, Did_Dname$Did)]
disgene <- cbind(disgene, disname)

sigdis <- disgene[is.element(disgene$gene, unlist(cg)),]

topdisgene <- tapply(as.character(sigdis$gene), as.character(sigdis$dis), c)
topdisgene <- lapply(topdisgene, function(x){paste(unique(x), collapse = ",")})
topdisgene <- data.frame(dis = names(topdisgene), gene = unlist(topdisgene))

write.table(topdisgene, file = "F:/MeRIPSinglebase/metdb2/result/sigdisgene.xls",
            quote = F, col.names = T, row.names = F, sep = "\t")

save(disgene, sigdis, file = "F:/MeRIPSinglebase/metdb2/data/Pro_Dis.RData")



sigdis[sigdis$gene == "PTCH1",]
sigdis[sigdis$gene == "SUFU",]

sigdis <- tapply(as.character(sigdis$gene), as.character(sigdis$dis), c)
sigdis[unlist(lapply(sigdis, length)) > 1]








