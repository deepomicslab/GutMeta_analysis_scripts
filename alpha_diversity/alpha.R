library(getopt)
library(vegan)

command <- matrix(c(
  "input",  "i", 1, "character", "Input abundance",
  "group", "g", 1, "character",  "Sample group information",
  "output",  "o", 1, "character",  "Output file name of diversity output",
  "pvalue",  "p", 1, "character",  "Output file name of pvalue output",
  "method",  "m", 2, "character",  "Method for alpha diversity calculation [shannon/simpson/invsimpson/ACE/Chao1/observedSpecies]",
  "reads",  "r", 2, "integer",  "Reads number for ACE/Chao1/observedSpecies",
  "title",  "t", 1, "character", "Taxnomy level as rowname",
  "testing_method",  "e", 1, "character", "Method for testing [wilcox.test/t.test/kruskal.test/aov]",
  "append",  "a", 1, "logical", "Append to output file or not [True/False]",
  "help",   "h", 0, "logical",  "Help"),
  byrow=TRUE, ncol=5)
args = getopt(command)

if( !is.null(args$help) | is.null(args$input) | is.null(args$group) | is.null(args$output) | is.null(args$pvalue) | is.null(args$title) | is.null(args$append) | is.null(args$testing_method)){
  cat(paste(getopt(command, usage = T), "\n"))
  quit()
}

if( is.null(args$method)){
  args$method = "shannon"
}
if( is.null(args$reads)){
  args$reads = 500000
}

ABD <- args$input
GROUP <- args$group
OUTPUT <- args$output
PVALUE <- args$pvalue
METHOD <- args$method
READS <- args$reads
TAX <- args$title
APPEND <- args$append
TESTING <- args$testing_method

# ABD <- "C:/Users/yiqijiang3/jupyter/MergedInput.abundance.g.tsv"
# GROUP <- "C:/Users/yiqijiang3/jupyter/MergedInput.group_info.tsv"
# OUTPUT <- "C:/Users/yiqijiang3/jupyter/alpha.tsv"
# PVALUE <- "C:/Users/yiqijiang3/jupyter/p.tsv"
# TAX <- "g"
# METHOD <- "shannon"
# READS = 500000
# APPEND <- TRUE



alpha_diversity <- function(x, method = "shannon", reads_n = 500000) {
  INDICES <- c("shannon", "simpson", "invsimpson","observedSpecies ","Chao1","ACE")
#  method <- match.arg(method, INDICES)
  
  tmp = matrix(as.integer((x/rowSums(x))*reads_n*2), nrow = nrow(x), ncol = ncol(x))
  tmp_flat = as.data.frame((rrarefy(tmp,reads_n)))
  for (method in INDICES) {
    if (method == "shannon" | method == "invsimpson" | method == "simpson") {
      result <- diversity(x, index = method)
    } else {
      if (method == "observedSpecies")
        result <- estimateR(tmp_flat)[1, ]
      if (method == "Chao1")
        result <- estimateR(tmp_flat)[2, ]
      if (method == "ACE")
        result <- estimateR(tmp_flat)[4, ]
    }
    if(exists("output"))
      output <- data.frame(output,result)
    else
      output <- result
  }
  output[is.na(output)] <- 0
  colnames(output) <- INDICES
  rownames(output) <- rownames(x)
  output
}


abd_data <- read.table(ABD,header=T,row.names = 1,sep = "\t")
result <- alpha_diversity(t(abd_data))
res <- result[METHOD]
colnames(res)<-c(TAX)
res$samp <- rownames(res)

phen <- read.table(GROUP,head=T)
colnames(phen) <- c("ID","Group")
phen_samp <- cbind("samp"=as.vector(phen[,1]), "state"=as.vector(phen[,2]))
groupname <- as.matrix(sort(unique(phen$Group)))


dat <- merge(phen_samp, res, by="samp",sort =F)
dat[,c(TAX)] <- as.numeric(as.vector(sprintf("%0.4f",dat[,c(TAX)])))
dat <- dat[complete.cases(dat),]

testing_p <- function(data, x, y, method = "t.test"){
  x <- data[,x]
  y <- data[,y]
  if (method == "t.test" || method == "wilcox.test" || method == "kruskal.test")
    pvalue = do.call(method, list(x ~ y, data = data))$p.value
  if (method == "aov"){
    tmp <- summary(aov(x ~ y, data = data))
    pvalue = tmp[[1]][["Pr(>F)"]][1]
  }
    
  pvalue 
}

pvalue = testing_p(dat, TAX, "state", method = TESTING)
tmp <- data.frame(matrix(c(pvalue),dimnames = list(c(TAX),c("pvalue"))))

if (APPEND) {
  write.table(tmp,PVALUE,append = T,sep = "\t",quote=F,col.names=F)
  write.table(t(res[,TAX,drop=F]),OUTPUT,append = T,sep = "\t",quote=F,col.names=F)
} else{
  write.table(tmp,PVALUE,append = F,sep = "\t",quote=F,col.names=NA)
  write.table(t(res[,TAX,drop=F]),OUTPUT,append = F,sep = "\t",quote=F,col.names=NA)
}
