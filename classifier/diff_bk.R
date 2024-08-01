library(getopt)
library(vegan)
library(dplyr)
library(reshape2)
library(ape)

command <- matrix(c(
  "input",  "i", 1, "character", "Input abundance",
  "group", "g", 1, "character",  "Sample group information",
  "output",  "o", 1, "character",  "Output prefix",
  "method",  "m", 1, "character", "Method for testing [wilcox.test/t.test/kruskal.test/aov]",
  "adjust",  "a", 2, "character",  "Method for adjusting p-values, default: fdr. [fdr/BH/hochberg/holm/bonferroni/hommel]",
  "title",  "t", 1, "character", "Taxnomy level as rowname",
  "pvalue_cutoff", "p", 2, "numeric", "Threshold of adjusted p-value as significant, defualt: 0.05",
  "mean_cutoff", "e", 2, "numeric", "Threshold of mean value to filter taxonomy, at least one group should have a mean relative abundance greater than this threshold, defualt: 0",
  "occ_cutoff", "c", 2, "numeric", "Threshold of occurence to filter taxonomy, the taxonomy should at least exist in the threshold ratio of samples of certain group, defualt: 0.1",
  "help",   "h", 0, "logical",  "Help"),
  byrow=TRUE, ncol=5)
args = getopt(command)

if( !is.null(args$help) | is.null(args$input) | is.null(args$group) | is.null(args$output) | is.null(args$title)){
  cat(paste(getopt(command, usage = T), "\n"))
  quit()
}

if( is.null(args$method)){
  args$method = "wilcox.test"
}
if( is.null(args$adjust)){
  args$adjust = "BH"
}
if( is.null(args$pvalue_cutoff)){
  args$pvalue_cutoff = 0.05
}
if( is.null(args$mean_cutoff)){
  args$mean_cutoff = 0
}
if( is.null(args$occ_cutoff)){
  args$occ_cutoff = 0.1
}

ABD <- args$input
GROUP <- args$group
OUTPUT <- args$output
METHOD <- args$method
TAX <- args$title
ADJUST <- args$adjust

MEAN <- args$mean_cutoff
OCC <- args$occ_cutoff
PVALUE <- args$pvalue_cutoff

OUTPUT1 <- paste(OUTPUT,".stat.tsv",sep='')
OUTPUT2 <- paste(OUTPUT,".stat.pass.tsv",sep='')
OUTPUT3 <- paste(OUTPUT,".abd.pass.tsv",sep='')


testing_p <- function( x, y, method = "t.test"){
  if (method == "t.test" || method == "wilcox.test" || method == "kruskal.test")
    pvalue = do.call(method, list(x ~ y))$p.value
  if (method == "aov"){
    tmp <- summary(aov(x ~ y))
    pvalue = tmp[[1]][["Pr(>F)"]][1]
  }
  pvalue
}

outfile <-file(OUTPUT1,"w")

phen <- read.table(GROUP,head=T,sep="\t")
phen_samp <- cbind("samp"=as.vector(phen[,1]), "state"=as.vector(phen[,2]))
groupname <- as.matrix(sort(unique(phen[,2])))

title <- TAX

for (i in 1:nrow(groupname))
{
        mean <- paste("mean(",groupname[i,1],")",sep="")
        sd <- paste("sd(",groupname[i,1],")",sep="")
        occ <- paste("occ-rate(",groupname[i,1],")",sep="")
        title <- paste(title,mean,sd,occ,sep="\t")
}

outmatrixname <- paste(title,"enriched","pvalue",sep="\t")
writeLines(outmatrixname, con=outfile, sep="\n")

abd_file <- file(ABD,"r")
line <- readLines(abd_file, n=1)

samp <- as.vector(unlist(strsplit(line, split="\t")))

while(length(line <- readLines(abd_file, n=1))) {
	line <- as.vector(unlist(strsplit(line, split="\t")))
	tmp <- cbind(samp=samp[-1], abund=line[-1])
	dat <- merge(phen_samp, tmp, by="samp",sort =F)
	dat$abund <- as.numeric(as.vector((dat$abund)))
	dat <- dat[complete.cases(dat),]
	mean.sd <- as.matrix(aggregate(dat$abund,by=list(dat$state),FUN=function(x)c(mean=sprintf("%0.9f",mean(x)),sd=sprintf("%0.9f",sd(x)))))
	Occ <- as.matrix(aggregate(dat$abund,by=list(dat$state),FUN=function(x)sum(x != 0)/length(x)))
	if ( (mean.sd[1,2] < MEAN || as.numeric(Occ[1,2]) < OCC ) && (mean.sd[2,2] < MEAN || as.numeric(Occ[2,2]) < OCC) ){next}
	test.pvalue <- testing_p(dat$abund, dat$state, method = METHOD)
	Rank <- rank(dat$abund)
	Rank_mean <- as.matrix(aggregate(Rank,by=list(dat$state),FUN=function(x)c(mean=mean(x))))
	enriched <- as.character(Rank_mean[which.max(Rank_mean[,2]),1])
	output <- paste(line[1],mean.sd[1,2],mean.sd[1,3],format(as.numeric(Occ[1,2]),digits=3),mean.sd[2,2],mean.sd[2,3],format(as.numeric(Occ[2,2]),digits=3),enriched,test.pvalue,sep="\t")
	writeLines(output, con=outfile, sep="\n")

}

close(abd_file)
close(outfile)
data <- read.table(OUTPUT1,head=T,check.names=F,sep="\t")
p_adjusted <- p.adjust(data$pvalue, method=ADJUST)
data$qvalue <- p_adjusted
data <- data %>% mutate_if(is.numeric, round, digits = 8)
data <- data %>% arrange(.[[8]], desc(.[[2]]))
write.table(format(data,scientific=FALSE),file=OUTPUT1,row.names=F, col.names=T, quote=F, sep="\t")


pass <- data %>% filter (qvalue < PVALUE)
write.table(pass, OUTPUT2, sep = "\t",quote = F, row.names = F)

order <- as.factor(pass[[1]])
abd_data <- read.table(ABD,header=T,sep = "\t")
pass_abd <- abd_data %>% filter (abd_data[[1]] %in% pass[[1]])
pass_abd <- pass_abd[match(pass[[1]],pass_abd[[1]]),]
write.table(pass_abd, OUTPUT3, sep = "\t",quote = F, row.names = F)
