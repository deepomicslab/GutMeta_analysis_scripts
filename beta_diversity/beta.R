library(getopt)
library(vegan)
library(GUniFrac)
library(dplyr)
library(reshape2)
library(ape)

command <- matrix(c(
  "input",  "i", 1, "character", "Input abundance",
  "group", "g", 1, "character",  "Sample group information",
  "output",  "o", 1, "character",  "Output file name of diversity output",
  "pvalue",  "p", 1, "character",  "Output file name of pvalue output",
  "method",  "m", 2, "character",  "Method for beta diversity calculation [bray/euclidean/jaccard/manhattan/weighted_unifrac/unweighted_unifrac]",
  "reads",  "r", 2, "integer",  "Reads number for unifrac",
  "title",  "t", 1, "character", "Taxnomy level as rowname",
  "testing_method",  "e", 1, "character", "Method for testing [wilcox.test/t.test/kruskal.test/aov]",
  "tree",  "c", 2, "character", "Tree in newick format for unifrac distance calculation",
  "help",   "h", 0, "logical",  "Help"),
  byrow=TRUE, ncol=5)
args = getopt(command)

if( !is.null(args$help) | is.null(args$input) | is.null(args$group) | is.null(args$output) | is.null(args$pvalue) | is.null(args$title) | is.null(args$testing_method)){
  cat(paste(getopt(command, usage = T), "\n"))
  quit()
}

if( is.null(args$method)){
  args$method = "bray"
}
if( is.null(args$reads)){
  args$reads = 500000
}
if( is.null(args$tree)){
  args$tree = NA
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
TREE <- args$tree

OUTPUT1 <- paste(OUTPUT,".matrix.tsv",sep='')
OUTPUT2 <- paste(OUTPUT,".list.txt",sep='')



beta_diversity <- function(x, method = "bray", reads_n = 500000, tree = NA) {
  INDICES <- c("bray", "euclidean", "jaccard","manhattan ","weighted_unifrac","unweighted_unifrac")
  #  method <- match.arg(method, INDICES)

  if (method == "bray" | method == "euclidean" | method == "jaccard" | method == "manhattan") {
    result <- as.matrix(vegdist(x, method = method))
  } else {
    tmp = matrix(as.integer((x/rowSums(x))*reads_n*2), nrow = nrow(x), ncol = ncol(x))
    tmp_flat = as.data.frame((rrarefy(tmp,reads_n)))
    colnames(tmp_flat)<-colnames(x)
    rownames(tmp_flat)<-rownames(x)
    tree <- read.tree(tree)
    unifracs <- suppressWarnings(GUniFrac(tmp_flat,tree)$unifracs)
    if (method == "weighted_unifrac"){
      result <- unifracs[, , "d_1"]
    }else{result <- unifracs[, , "d_UW"]}
  }
  result
}

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

abd_data <- read.table(ABD,header=T,row.names = 1,sep = "\t")
result <- beta_diversity(t(abd_data), method = METHOD, reads_n = READS, tree = TREE)
write.table(data.frame("beta_diverisity_matrix"=rownames(result),result),OUTPUT1,sep = "\t",quote = F, row.names = F)


res = result
res[lower.tri(res)] <- 0
test<-melt(res,value.name = "dist",varname = c("sample1", "sample2"))
filter<-test[test$dist>0,]


phen <- read.table(GROUP,head=T)
colnames(phen) <- c("ID","Group")
phen_samp <- cbind("samp"=as.vector(phen[,1]), "state"=as.vector(phen[,2]))
groups <- as.matrix(sort(unique(phen$Group)))

merged<-left_join(left_join(filter,phen,by=c("sample1"="ID")),phen,by=c("sample2"="ID"))
merged$group<-paste(merged$Group.x,merged$Group.y,sep = "-")

g1<-paste(groups[1],groups[2],sep = "-")
g2<-paste(groups[2],groups[1],sep = "-")
merged[merged==g1]<-g2

write(TAX,OUTPUT2,append = F)
for (g in levels(as.factor(merged$group))){
  tmp<-merged[merged$group==g,"dist"]
  line<-paste(t(unlist(tmp)),collapse = "\t")
  line<-paste(g,line,sep = "\t")
  write(line,OUTPUT2,append=T)
}

write(paste(TAX,'p-value',sep="\t"),PVALUE,append = F)

groups <-as.matrix(sort(unique(merged$group)))
for (i in 1:(nrow(groups)-1)){
    for (j in (i+1):nrow(groups)){
        tmp = merged[merged$group==groups[i]|merged$group==groups[j],]
        pvalue = testing_p(tmp,"dist","group", method = TESTING)
        line = paste(paste(groups[i],groups[j],sep=":"),pvalue,sep="\t")
	write(line,PVALUE,append=T)
    }
}


