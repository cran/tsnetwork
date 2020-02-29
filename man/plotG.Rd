\name{plotG}
\alias{plotG}

\title{ Graph visualization }

\description{
Creates a graph from adjacency matrices
 }
 
\usage{ plotG(G, mod, vertex.size = 5, label = TRUE, main = "Network") }

\arguments{
	\item{G}{ a pxp (adjacency) matrix}
	\item{mod}{specifies how a graph should be created from the supplied matrix. Possible values are: directed, undirected.}
	\item{vertex.size}{The size of vertices in the graph.}
	\item{label}{The vertex labels. The default is TRUE.}
	\item{main}{Overall title for the graph.}
	}
	
\value{
A graph
}

\author{
Pariya Behrouzi \cr
Maintainer: Pariya Behrouzi <pariya.behrouzi@gmail.com>
}

\seealso{
\code{\link{tsnetwork}}
}