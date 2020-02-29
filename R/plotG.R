plotG  <- function(G, mod, vertex.size = 5, label = TRUE, main = "Network"){
	p <- ncol(G)
	path <-  abs(sign(G)) - diag(rep(1,p))
	adj <- graph.adjacency(path, mode = mod)
	if(label) id <- colnames(G) else id <- NA
	l <- layout.fruchterman.reingold(adj)
	plot(adj, vertex.size = vertex.size, main = main, vertex.label = id , vertex.label.dist = 1.3)
}
