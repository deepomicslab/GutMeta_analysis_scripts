library(getopt)
library(vegan)
library(ica)
library(Rtsne)
library(uwot)
library(ecodist)

command <- matrix(c(
  "input",  "i", 1, "character", "Input abundance",
  "output",  "o", 1, "character",  "Output file name of matrix output",
  "method",  "m", 2, "character",  "Method for dimension reduction [PCA/ICA/tSNE/UMAP/NMDS/PCoA]",
  "axis",  "r", 2, "integer",  "Axis number[10]",
  "dist", "d", 2, "character",  "Method for beta diversity calculation [bray/euclidean/jaccard/manhattan]",
  "title",  "t", 1, "character", "Taxnomy level as rowname",
  "help",   "h", 0, "logical",  "Help"),
  byrow=TRUE, ncol=5)
args = getopt(command)

if( !is.null(args$help) | is.null(args$input) | is.null(args$output) | is.null(args$title)){
  cat(paste(getopt(command, usage = T), "\n"))
  quit()
}

if( is.null(args$method)){
  args$method = "PCA"
}
if( is.null(args$dist)){
  args$dist = "bray"
}
if( is.null(args$axis)){
  args$axis = 10
}

ABD <- args$input
OUTPUT <- args$output
METHOD <- args$method
TAX <- args$title
AXIS_N <- args$axis
DIST_METHOD <- args$dist

reduce_dim <- function(x, method = "PCA", axis_n = 10, dist_method = "bray") {
  INDICES <- c("PCA", "ICA", "tSNE","UMAP","NMDS","PCoA")
  #  method <- match.arg(method, INDICES)
  DIST_INDICES <- c("bray","euclidean","jaccard","manhattan")
  if (method == "PCA") {
    tmp <- prcomp(x)
    result <- tmp$x[,1:axis_n]
  } else if (method == "ICA") {
    tmp <- icafast(x, axis_n)
    result <- tmp$Y
  } else if (method == "tSNE") {
    tmp <- Rtsne(x)
    result<-tmp$Y
  } else if (method == "UMAP") {
    result <- umap(x, n_components = axis_n)
  } else if (method == "PCA") {
    dist <- vegdist(x, method = dist_method)
    tmp <- pco(dist)
    result <- tmp$vectors[,1:axis_n]
  } else if (method == "NMDS") {
    tmp <- vegan::metaMDS(x, distance = dist_method, k = axis_n)
    result<-tmp$points
  }
  rownames(result)<-rownames(x)
  colnames(result)<-paste("Axis",seq(1,ncol(result)),sep = '')
  result
}

abd_data <- read.table(ABD,header=T,row.names = 1,sep = "\t")
result <- reduce_dim(t(abd_data), method = METHOD, axis_n = AXIS_N, dist_method = DIST_METHOD)


tmp<-data.frame(rownames(result))
names(tmp)<-TAX
write.table(data.frame(tmp,result),OUTPUT,sep = "\t",quote = F, row.names = F)
