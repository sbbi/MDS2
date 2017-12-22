args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one arguments must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "outfile.txt"
}

library("igraph")

###graph_merged3 <- read.delim("simi_g.txt", header = F, as.is = T)

graph_merged3 <- read.delim(args[1], header = F, as.is = T)

colnames(graph_merged3) <- c("source", "target", "weights")
#graph_merged <- subset(graph_merged, weights >= median(graph_merged$weights))

graph_merged2 <- subset(graph_merged3, weights > 0)
graph_merged <- as.matrix(graph_merged2)

#g3 <- graph.edgelist(as.matrix(graph_merged[,1:2], drop = F), directed = F)

if(nrow(graph_merged) > 1){

  g3 <- graph.edgelist(as.matrix(graph_merged[,1:2], drop = F), directed = F)
  E(g3)$weight <- as.numeric(graph_merged[,3])
  #layout <- layout_with_lgl(g3)


### louvain
  lo3 <- cluster_louvain(g3, weights = E(g3)$weight)
  lo3_modules <- communities(lo3)

  n.obs <- sapply(lo3_modules, length)
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(lo3_modules, "[", i = seq.max))
}else{
  mat <- t(as.data.frame(unique(unlist(graph_merged3[,1:2]))))
}

write.table(mat, args[2], row.names = F, col.names = F, quote = F, na = "", sep = "\t")
#write.table(mat, "null.txt", row.names = F, col.names = F, quote = F, na = "", sep = "\t")