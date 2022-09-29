library(getopt)
library(cluster)
library(clusterSim)
library(ade4)

command <- matrix(c(
  "input",  "i", 1, "character", "Input abundance",
  "group", "g", 1, "character",  "Sample group information",
  "output",  "o", 1, "character",  "Output file prefix",
  "method",  "m", 2, "character",  "Method for visualization [PCA/PCoA]",
  "maxk",  "k", 2, "integer",  "Maxinum number of clusters to find the optimal k. default:10",
  "title",  "t", 1, "character", "Taxnomy level as rowname",
  "help",   "h", 0, "logical",  "Help"),
  byrow=TRUE, ncol=5)
args = getopt(command)

if( !is.null(args$help) | is.null(args$input) | is.null(args$output) | is.null(args$title)){
  cat(paste(getopt(command, usage = T), "\n"))
  quit()
}

if( is.null(args$method)){
  args$method = "PCoA"
}
if( is.null(args$maxk)){
  args$maxk = 10
}

ABD <- args$input
OUTPUT <- args$output
METHOD <- args$method
MAX_K <- args$maxk
TAX <- args$title

OUTPUT1 <- paste(OUTPUT,".cluster.tsv",sep='')
OUTPUT2 <- paste(OUTPUT,".CH_index.tsv",sep='')
OUTPUT3 <- paste(OUTPUT,".",METHOD,".coordinates.tsv",sep='')

dist_JSD <- function(InMat, pseudocount = 0.0000000000001){
        kld <- function(x, y){ sum(x * log(x/y)) }
        jsd <- function(x,y){ sqrt(0.5*kld(x, (x+y)/2) + 0.5*kld(y, (x+y)/2)) }
        ncol1 <- length(colnames(InMat))
        colname <- colnames(InMat)
        resultMatrix <- matrix(0, ncol1, ncol1)
        InMat <- apply(InMat, 1:2, function(x) ifelse (x==0, pseudocount, x))
        for(i in 1:ncol1){
                for(j in 1:ncol1){
                        resultMatrix[i, j] <- jsd(as.vector(InMat[, i]), as.vector(InMat[, j]))
                        }
        }
        colname -> colnames(resultMatrix) -> rownames(resultMatrix)
        resultMatrix <- as.dist(resultMatrix)
        attr(resultMatrix, "method") <- "dist"
        return(resultMatrix)
}

rm_noise <- function(dataframe, percent=0.01, top=NULL){ ## cutoff 0.01%
        Matrix <- dataframe
        big.ones <- rowSums(Matrix)*100 / (sum(rowSums(Matrix))) > percent
        Matrix.1 <- Matrix[big.ones, ]
        return(Matrix.1)
}
pam_clustering <- function(x, k) { # x is a distance matrix and k the number of clusters
        cluster <- as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
        return(cluster)
}

get_opt_k <- function (data, dist, max_k = 10) {
	num.cluster <- NULL
	num.cluster[1]=0
	for(k in 2:max_k){
		cluster_temp <- pam_clustering(dist, k)
		num.cluster[k] <- index.G1(data, cluster_temp, d=dist, centrotypes="medoids")
	}
	CH_index <- data.frame(num.cluster)
	colnames(CH_index) <- "CH_index"
	return(CH_index)
}

abd_data <- read.table(ABD,header=T,row.names = 1,sep = "\t")
dist <- dist_JSD(abd_data)
CH_index = get_opt_k(t(abd_data),dist,MAX_K)
opt_k = which.max(as.matrix(CH_index))
cluster <- pam_clustering(dist, k=opt_k)
cluster.result <- cbind(colnames(abd_data),cluster)
colnames(cluster.result) <- c(TAX,"fac")
write.table(cluster.result, OUTPUT1, row.names = F, col.names = T, sep = "\t", quote = F)
write.table(data.frame("k"=rownames(CH_index),CH_index), OUTPUT2, row.names = F, col.names = T, sep = "\t", quote = F)


#visualization
if(METHOD == "PCoA") {
	dim_redu <- dudi.pco(dist, scannf=F, nf=10)
}else{
	dim_redu <- dudi.pca(t(abd_data), scale=F,scannf=F, nf=10)
}


colnames(dim_redu$li) <- paste("PC", 1:ncol(dim_redu$li), sep = "")
result = data.frame(type = rownames(dim_redu$li))
result = cbind(result, dim_redu$li)
names(result)[names(result)=="type"]=TAX
write.table(result, OUTPUT3, quote = F, sep = "\t", row.names = F)



