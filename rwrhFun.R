
rwrh <- function(siggene, mmpath, sigdis) {
  
  load(mmpath)
  load("F:/MeRIPSinglebase/metdb2/data/omim/Did_Dname.RData")
  
  library(Matrix)
  
  seed <- match(siggene, gene)
  y0g <- Matrix(0, ncol=1, nrow=length(gene))
  y1g <- y0g
  y0g[seed] <- 1/length(seed)*0.5
  yg <- y0g
  
  seed <- is.element(dis, sigdis)
  y0d <- Matrix(0, ncol=1, nrow=length(dis))
  y1d <- y0d
  y0d[seed] <- 1/length(seed)*0.5
  # y0d <- Matrix(1/length(seed)*0.5, ncol=1, nrow=length(dis))
  yd <- y0d
  
  y0 <- rbind(y0g, y0d)
  y <- rbind(yg, yd)
  y1 <- rbind(y1g, y1d)
  
  while (sum(abs(y-y1))>10^-10)
  {
    y1 <- y
    y <- 0.5*mm%*%y1+0.5*y0
    # print(sum(abs(y-y1)))
  }
  y <- as.matrix(y)
  
  yd <- y[(length(gene)+1):length(y)]
  id <- order(yd, decreasing = T)
  topdis <- Did_Dname$Dname[match(dis[id[1:50]], Did_Dname$Did)]
  
  print(topdis[1:5])
  
  yg <- y[1:length(gene)]
  id0 <- order(yg, decreasing = T)
  topgene <- gene[id0[1:500]]
  
  print(topgene[1:5])
  
  return(list(topgene = topgene, topdis = topdis))
}

## rwrh random ----------------------------------------------------------

rwrhr <- function(siggene, mmpath, sigdis) {
  
  load(mmpath)
  load("F:/MeRIPSinglebase/metdb2/data/omim/Did_Dname.RData")
  
  library(Matrix)
  
  seed <- match(siggene, gene)
  y0g <- Matrix(0, ncol=1, nrow=length(gene))
  y1g <- y0g
  y0g[seed] <- 1/length(seed)*0.5
  yg <- y0g
  
  seed <- is.element(dis, sigdis)
  y0d <- Matrix(0, ncol=1, nrow=length(dis))
  y1d <- y0d
  # y0d[seed] <- 1/length(seed)*0.5
  y0d <- Matrix(1/length(seed)*0.5, ncol=1, nrow=length(dis))
  yd <- y0d
  
  y0 <- rbind(y0g, y0d)
  y <- rbind(yg, yd)
  y1 <- rbind(y1g, y1d)
  
  while (sum(abs(y-y1))>10^-10)
  {
    y1 <- y
    y <- 0.5*mm%*%y1+0.5*y0
    # print(sum(abs(y-y1)))
  }
  y <- as.matrix(y)
  
  yd <- y[(length(gene)+1):length(y)]
  id <- order(yd, decreasing = T)
  topdis <- Did_Dname$Dname[match(dis[id[1:50]], Did_Dname$Did)]
  
  print(topdis[1:5])
  
  yg <- y[1:length(gene)]
  id0 <- order(yg, decreasing = T)
  topgene <- gene[id0[1:500]]
  
  print(topgene[1:5])
  
  return(list(topgene = topgene, topdis = topdis))
}






























