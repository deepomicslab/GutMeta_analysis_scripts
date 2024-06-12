library(getopt)
library("randomForest")
library("pROC")

command <- matrix(c(
  "input",  "i", 1, "character", "Input abundance",
  "group", "g", 1, "character",  "Sample group information",
  "output",  "o", 1, "character",  "Output prefix.",
  "rep",  "r", 2, "integer",  "Replicates number for cross-validation. default = 5",
  "cv",  "c", 2, "integer",  "CV fold. default = 10",
  "scale",  "s", 1, "logical", "Log scale the abundance or not [True/False]",
  "help",   "h", 0, "logical",  "Help"),
  byrow=TRUE, ncol=5)
args = getopt(command)

if( !is.null(args$help) | is.null(args$input) | is.null(args$group) | is.null(args$output) ){
  cat(paste(getopt(command, usage = T), "\n"))
  quit()
}

if( is.null(args$rep)){
  args$rep = 5
}
if( is.null(args$cv)){
  args$cv = 10
}
if( is.null(args$scale)){
  args$scale = NA
}

ABD <- args$input
GROUP <- args$group
OUTPUT <- args$output
SCALE <- args$scale
REP_n <- args$rep
CV_n <- args$cv



OUTPUT1 <- paste(OUTPUT,".cross_validation_error.tsv",sep='')
OUTPUT2 <- paste(OUTPUT,".cross_validation_pick.txt",sep='')
OUTPUT3 <- paste(OUTPUT,".feature.importance.tsv",sep='')
OUTPUT4 <- paste(OUTPUT,".cross_validation.marker.predict.in.train.tsv",sep='')
OUTPUT5 <- paste(OUTPUT,".train.ci_result.tsv",sep='')
OUTPUT6 <- paste(OUTPUT,".train.auc_info.txt",sep='')

abd_data <- read.table(ABD,header=T,row.names = 1,sep = "\t")
phen <- read.table(GROUP,header=T,row.names = 1,sep = "\t")
colnames(phen) <- c("Group")
groups <- as.matrix(sort(unique(phen$Group)))
outcome = phen$Group
outcome <- sub(groups[1],"0",outcome)
outcome <- sub(groups[2],"1",outcome)
outcome <- as.factor(outcome)
X <- as.data.frame(t(abd_data))
X$outcome = outcome



rfcv1 <- function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5,mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE,...)
{
	classRF <- is.factor(trainy)
	n <- nrow(trainx)
	p <- ncol(trainx)
	if (scale == "log") 
	{
		k <- floor(log(p, base = 1/step))
		n.var <- round(p * step^(0:(k - 1)))
		same <- diff(n.var) == 0
		if (any(same))
			n.var <- n.var[-which(same)]
			if (!1 %in% n.var)
				n.var <- c(n.var, 1)
	}else {
		n.var <- seq(from = p, to = 1, by = step)
	}
	k <- length(n.var)
	cv.pred <- vector(k, mode = "list")
	for (i in 1:k) cv.pred[[i]] <- rep(0,length(trainy))
		if (classRF) 
		{
			f <- trainy
		}else {
			f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
		}
		nlvl <- table(f)
		idx <- numeric(n)
		for (i in 1:length(nlvl)) 
		{
			idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold,length = nlvl[i]))
		}
		res=list()
		accuracy=list()
		gini=list()
		for (i in 1:cv.fold) 
		{
			all.rf <- randomForest(trainx[idx != i, , drop = FALSE],trainy[idx != i],importance = TRUE)
			aa = predict(all.rf,trainx[idx == i, , drop = FALSE],type="prob")
			cv.pred[[1]][idx == i] <- as.numeric(aa[,2])
			impvar <- (1:p)[order(all.rf$importance[, 3], decreasing = TRUE)]
			acc_value <- as.numeric(importance(all.rf)[, 3])
			gini_value <- as.numeric(importance(all.rf)[, 4])
			res[[i]]=impvar
			accuracy[[i]]= acc_value
			gini[[i]] = gini_value
			for (j in 2:k) 
			{
				imp.idx <- impvar[1:n.var[j]]
				sub.rf <- randomForest(trainx[idx != i, imp.idx,drop = FALSE], trainy[idx != i])
				bb <- predict(sub.rf,trainx[idx ==i,imp.idx, drop = FALSE],type="prob")
				cv.pred[[j]][idx == i] <- as.numeric(bb[,2])
				if (recursive) 
				{
					impvar <- (1:length(imp.idx))[order(sub.rf$importance[,3], decreasing = TRUE)]
					acc_value <- as.numeric(importance(sub.ref)[, 3])
					gini_value <- as.numeric(importance(sub.ref)[, 4])
				}
				NULL
			}
			NULL
		}
		if (classRF) 
		{
			error.cv <- sapply(cv.pred, function(x) mean(factor(ifelse(x>0.5,1,0))!=trainy))
		}else {
			error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
		}
		names(error.cv) <- names(cv.pred) <- n.var
		list(n.var = n.var, error.cv = error.cv, predicted = cv.pred,res=res, acc=accuracy, gini=gini)
}

set.seed(999)
if (SCALE){
        result <- replicate(REP_n, rfcv1(X[,-ncol(X)], X$outcome, cv.fold=CV_n, scale = "log"), simplify=FALSE)
}else{
        result <- replicate(REP_n, rfcv1(X[,-ncol(X)], X$outcome, cv.fold=CV_n), simplify=FALSE)
}
error.cv <- sapply(result, "[[", "error.cv")
error.cv.cbm <- cbind(rowMeans(error.cv), error.cv)
cutoff <- min(error.cv.cbm[,1])+sd(error.cv.cbm[,1])
cutoff.num <- nrow(error.cv.cbm[error.cv.cbm[,1]<cutoff,])
optimal.set.feature.num <- as.numeric(rownames(error.cv.cbm[error.cv.cbm[,1]<cutoff,])[cutoff.num])
write.table(data.frame( "K"=rownames(error.cv.cbm),error.cv.cbm),OUTPUT1,quote=F,sep="\t",row.names=F)

tax_n = nrow(abd_data)
item_n = REP_n * CV_n
#----- pick marker by corossvalidation -----
k = 1
b <- matrix(0,ncol=tax_n,nrow=item_n)
for(i in 1:REP_n)
{
        for(j in 1:CV_n)
        {
                b[k,] <- result[[i]]$res[[j]]
                k = k+1
        }
}

mlg.list <- b[,1:CV_n]
list <- c()
k = 1
for(i in 1:CV_n)
{
        for(j in 1:item_n)
        {
                list[k] <- mlg.list[j,i]
                k = k+1
        }
}

mlg.sort <- as.matrix(table(list))
mlg.sort <- mlg.sort[rev(order(mlg.sort[,1])),]
pick <- as.numeric(names(head(mlg.sort,optimal.set.feature.num)))
tmp = X[,-ncol(X)]
mlg.pick <- colnames(tmp)[pick]
pick_line <- paste(t(mlg.pick),collapse=",")
cat(pick_line,file = OUTPUT2)

#----- importance -----
all.importance <- randomForest(outcome~., data = X,importance = TRUE)
imp <- data.frame(importance(all.importance))
colnames(imp) <- c("DecreaseAccuracy","DecreaseGini","MeanDecreaseAccuracy","MeanDecreaseGini")
out_imp <- imp[,c(3,4)]
write.table(data.frame("Taxonomy"=rownames(out_imp),out_imp),OUTPUT3,sep="\t",quote=F,row.names=F)

#----- train.set -----
train1 <- X[,c(pick,tax_n+1)] ##
set.seed(999)
train1 <- data.frame(train1)
train1.rf <- randomForest(outcome~., data = train1,importance = TRUE)
train1.pre <- predict(train1.rf,type="prob")
p.train <- train1.pre[,2]
combine <- as.data.frame(cbind(predict.value=as.matrix(p.train)[match(rownames(as.matrix(p.train)),rownames(as.matrix(phen)))],as.matrix(phen)))

write.table(data.frame("Sample"=rownames(combine),combine),OUTPUT4,sep="\t",quote=F,row.names=F)

#----- ROC in train -----
train.roc <- roc(outcome,p.train,percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
sens.ci <- as.data.frame(ci.se(train.roc, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c("sp","se.low","se.median","se.high"))
write.table(sens.ci,OUTPUT5,sep="\t",quote=F,col.names=T,row.names=F)

#----- AUC info -----
auc_info1 <- paste("AUC =",round(train.roc$ci[2],2))
auc_info2 <- paste("95% CI:",round(train.roc$ci[1],2),"-",round(train.roc$ci[3],2))

cat(auc_info1,file = OUTPUT6, sep = "\n")
cat(auc_info2,file = OUTPUT6, append = T)
