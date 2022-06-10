library(reshape2)
library(vegan)
library(randomForest)
library(pROC)

ROSMAP_clean <- na.omit(ROSMAP[DEG_union,])

ROSMAP_clinical_clean <- ROSMAP_clinical[,c("tissue", "msex", "race", "cogdx")]
colnames(ROSMAP_clinical_clean) = c("tissue", "sex", "race", "CDR")
ROSMAP_clinical_clean$sex <- gsub("1", "male", ROSMAP_clinical_clean$sex)
ROSMAP_clinical_clean$sex <- gsub("0", "female", ROSMAP_clinical_clean$sex)
ROSMAP_clinical_clean$race <- gsub("1", "W", ROSMAP_clinical_clean$race)
ROSMAP_clinical_clean$race <- gsub("2", "B", ROSMAP_clinical_clean$race)
ROSMAP_clinical_clean$race <- gsub("3", "W", ROSMAP_clinical_clean$race)
ROSMAP_clinical_clean$CDR <- gsub("1", "CTR", ROSMAP_clinical_clean$CDR)
ROSMAP_clinical_clean$CDR <- gsub("2", "MCI", ROSMAP_clinical_clean$CDR)
ROSMAP_clinical_clean$CDR <- gsub("3", "MCI", ROSMAP_clinical_clean$CDR)
ROSMAP_clinical_clean$CDR <- gsub("4", "AD", ROSMAP_clinical_clean$CDR)
ROSMAP_clinical_clean$CDR <- gsub("5", "AD", ROSMAP_clinical_clean$CDR)
ROSMAP_clinical_clean$CDR <- gsub("6", "Dementia", ROSMAP_clinical_clean$CDR)
ROSMAP_clinical_clean$status <- ROSMAP_clinical_clean$CDR

dim(ROSMAP_clean)
dim(ROSMAP_clinical_clean)

ROSMAP_clinical_clean$CDR[ROSMAP_clinical_clean$CDR == "AD"] <- 2
ROSMAP_clinical_clean$CDR[ROSMAP_clinical_clean$CDR == "CTR"] <- 0
ROSMAP_clinical_clean$CDR[ROSMAP_clinical_clean$CDR == "MCI"] <- 1

exclude_id_ROSMAP <- c(which(rowSums(is.na(ROSMAP_clinical_clean)) > 0), which(ROSMAP_clinical_clean$CDR == "Dementia"))
ROSMAP_clean <- ROSMAP_clean[,-exclude_id_ROSMAP]
ROSMAP_clinical_clean <- ROSMAP_clinical_clean[-exclude_id_ROSMAP,]

dim(ROSMAP_clean)
dim(ROSMAP_clinical_clean)

rfcv1 <- function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5, mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE, ...){
	classRF <- is.factor(trainy)
	n <- nrow(trainx)
	p <- ncol(trainx)
	if (scale == "log") {
		k <- floor(log(p, base = 1/step))
		n.var <- round(p * step^(0:(k - 1)))
		same <- diff(n.var) == 0
	if (any(same))
		n.var <- n.var[-which(same)]
	if (!1 %in% n.var)
		n.var <- c(n.var, 1)
	} else {
		n.var <- seq(from = p, to = 1, by = step)
	}
		k <- length(n.var)
		cv.pred <- vector(k, mode = "list")
	for (i in 1:k) cv.pred[[i]] <- rep(0,length(trainy))
		if (classRF) {
			f <- trainy
		} else {
			f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
		}
		nlvl <- table(f)
		idx <- numeric(n)

	for (i in 1:length(nlvl)) {
		idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold,
		length = nlvl[i]))
	}
	res = list()
	for (i in 1:cv.fold) {
		all.rf <- randomForest(trainx[idx != i, , drop = FALSE], trainy[idx != i], importance = TRUE)
		aa = predict(all.rf,trainx[idx == i, , drop = FALSE], type = "prob")
		cv.pred[[1]][idx == i] <- as.numeric(aa[,2])
		impvar <- (1:p)[order(all.rf$importance[, 3], decreasing = TRUE)]
		res[[i]] = impvar
		for (j in 2:k) {
			imp.idx <- impvar[1:n.var[j]]
			sub.rf <- randomForest(trainx[idx != i, imp.idx, drop = FALSE], trainy[idx != i])
			bb <- predict(sub.rf, trainx[idx == i, imp.idx, drop = FALSE], type = "prob")
			cv.pred[[j]][idx == i] <- as.numeric(bb[,2])
			if (recursive) {
				impvar <- (1:length(imp.idx))[order(sub.rf$importance[,3], decreasing = TRUE)]
			}
	NULL
	}
	NULL
}
	if (classRF) {
		error.cv <- sapply(cv.pred, function(x) mean(factor(ifelse(x > 0.5, 1, 0)) != trainy))
	}
	else {
		error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
	}
	names(error.cv) <- names(cv.pred) <- n.var
	list(n.var = n.var, error.cv = error.cv, predicted = cv.pred, res = res)
}


roc_list1 <- list()
roc_list2 <- list()
roc_list3 <- list()
roc_list4 <- list()

id1_used <- union(which(ROSMAP_clinical_clean$matched_MSBB_subtype == "control"), which(ROSMAP_clinical_clean$matched_MSBB_subtype == "typical")) 
id2_used <- union(which(ROSMAP_clinical_clean$matched_MSBB_subtype == "control"), which(ROSMAP_clinical_clean$matched_MSBB_subtype == "intermediate")) 
id3_used <- union(which(ROSMAP_clinical_clean$matched_MSBB_subtype == "control"), which(ROSMAP_clinical_clean$matched_MSBB_subtype == "atypical")) 
id4_used <- union(which(ROSMAP_clinical_clean$matched_MSBB_subtype == "control"), which(ROSMAP_clinical_clean$matched_MSBB_subtype == "other")) 

used_id1 <- ROSMAP_clinical_clean$matched_MSBB_subtype[id1_used]
used_id2 <- ROSMAP_clinical_clean$matched_MSBB_subtype[id2_used]
used_id3 <- ROSMAP_clinical_clean$matched_MSBB_subtype[id3_used]
used_id4 <- ROSMAP_clinical_clean$matched_MSBB_subtype[id4_used]

used_id1[used_id1 == "control"] <- 0
used_id2[used_id2 == "control"] <- 0
used_id3[used_id3 == "control"] <- 0
used_id4[used_id4 == "control"] <- 0
used_id1[used_id1 == "typical"] <- 1
used_id2[used_id2 == "intermediate"] <- 1
used_id3[used_id3 == "atypical"] <- 1
used_id4[used_id4 == "other"] <- 1

used_id1 <- factor(used_id1, levels = c(0, 1))
used_id2 <- factor(used_id2, levels = c(0, 1))
used_id3 <- factor(used_id3, levels = c(0, 1))
used_id4 <- factor(used_id4, levels = c(0, 1))

ROSMAP_clean_1 <- ROSMAP_clean[, id1_used]
ROSMAP_clean_2 <- ROSMAP_clean[, id2_used]
ROSMAP_clean_3 <- ROSMAP_clean[, id3_used]
ROSMAP_clean_4 <- ROSMAP_clean[, id4_used]

i = 99
set.seed(i)
test_id1 <- sample(ncol(ROSMAP_clean_1), 0.1 * ncol(ROSMAP_clean_1))
used_id1_train <- used_id1[-test_id1]
used_id1_test <- used_id1[test_id1]

test_id2 <- sample(ncol(ROSMAP_clean_2), 0.1 * ncol(ROSMAP_clean_2))
used_id2_train <- used_id2[-test_id2]
used_id2_test <- used_id2[test_id2]

test_id3 <- sample(ncol(ROSMAP_clean_3), 0.1 * ncol(ROSMAP_clean_3))
used_id3_train <- used_id3[-test_id3]
used_id3_test <- used_id3[test_id3]


# Typical vs Control
result1 <- replicate(10, rfcv1(t(ROSMAP_clean_1), used_id1, cv.fold = 10, step = 0.9), simplify = FALSE)
error1.cv <- sapply(result1, "[[", "error.cv")
Sys.time()

# Intermediate vs Control  ---------------------------------------------------------------------------
result2 <- replicate(10, rfcv1(t(ROSMAP_clean_2), used_id2, cv.fold = 10, step = 0.9), simplify = FALSE)
error2.cv <- sapply(result2, "[[", "error.cv")
Sys.time()

# Atypical vs Control  ---------------------------------------------------------------------------
result3 <- replicate(10, rfcv1(t(ROSMAP_clean_3), used_id3, cv.fold = 10, step = 0.9), simplify = FALSE)
error3.cv <- sapply(result3, "[[", "error.cv")
Sys.time()

matplot(result1[[1]]$n.var, cbind(rowMeans(error1.cv), error1.cv), type="l", lwd=c(2, rep(1, ncol(error1.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 20, col="pink", lwd = 2)

matplot(result2[[1]]$n.var, cbind(rowMeans(error2.cv), error2.cv), type="l", lwd=c(2, rep(1, ncol(error2.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 50, col="pink", lwd = 2)

matplot(result3[[1]]$n.var, cbind(rowMeans(error3.cv), error3.cv), type="l", lwd=c(2, rep(1, ncol(error2.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 30, col="pink", lwd = 2)

#------------------------------------------------------------------------
# Best model 
set.seed(i)
train.rf1 <- randomForest(t(ROSMAP_clean_1[,-test_id1]), used_id1_train, importance = TRUE, mtry = 20)	
train.rf1
set.seed(i)
train.rf2 <- randomForest(t(ROSMAP_clean_2[,-test_id2]), used_id2_train, importance = TRUE, mtry = 20)	
train.rf2
set.seed(i)
train.rf3 <- randomForest(t(ROSMAP_clean_3[,-test_id3]), used_id3_train, importance = TRUE, mtry = 20)	
train.rf3

# Prediction
ROSMAP.pre1 <- predict(train.rf1, type = "prob")
ROSMAP.pre2 <- predict(train.rf2, type = "prob")
ROSMAP.pre3 <- predict(train.rf3, type = "prob")

ROSMAP.pre1_res <- as.data.frame(cbind(ROSMAP.pre1, used_id1_train))
ROSMAP.pre2_res <- as.data.frame(cbind(ROSMAP.pre2, used_id2_train))
ROSMAP.pre3_res <- as.data.frame(cbind(ROSMAP.pre3, used_id3_train))

library(vioplot)
colnames(ROSMAP.pre1) = c("CTR", "Typical")
colnames(ROSMAP.pre2) = c("CTR", "Intermediate")
colnames(ROSMAP.pre3) = c("CTR", "Atypical")

#Fig C1
vioplot(ROSMAP.pre1, col = c(3,4), ylab = "Prob")

#Fig C2		 
level1 <- used_id1_train
level1 <- factor(level1, levels = c(0, 1))

roc1 <- roc(level1, ROSMAP.pre1[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc1 <- roc(level1, ROSMAP.pre1[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified = FALSE, plot = TRUE, percent = roc1$percent, col = 2)
sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc1,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc1$ci[2],2),"%"),
paste("95% CI:",round(roc1$ci[1],2),"%-",round(roc1$ci[3],2),"%")))
print("ROC for test dataset 1 d")

#Fig C3
roc1_test <- roc(used_id1_test, Test_roc1[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc1_test <- roc(used_id1_test, Test_roc1[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, percent = roc1_test$percent, col = 2)
sens.ci <- ci.se(roc1_test, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc1_test,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc1_test$ci[2],2),"%"),
paste("95% CI:",round(roc1_test$ci[1],2),"%-",round(roc1_test$ci[3],2),"%")))
print(Sys.time())


#Fig D1	
vioplot(ROSMAP.pre2, col = c(3,2), ylab = "Prob")

#Fig D2	
level2 <- used_id2_train
level2 <- factor(level2, levels = c(0, 1))

roc2 <- roc(level2, ROSMAP.pre2[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc2 <- roc(level2, ROSMAP.pre2[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified = FALSE, plot = TRUE, percent = roc2$percent, col = 2)
sens.ci <- ci.se(roc2, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc2,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc2$ci[2],2),"%"),
paste("95% CI:",round(roc2$ci[1],2),"%-",round(roc2$ci[3],2),"%")))
print("ROC for test dataset 2 d")

#Fig D3
roc2_test <- roc(used_id2_test, Test_roc2[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc2_test <- roc(used_id2_test, Test_roc2[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, percent = roc2_test$percent, col = 2)
sens.ci <- ci.se(roc2_test, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc2_test,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc2_test$ci[2],2),"%"),
paste("95% CI:",round(roc2_test$ci[1],2),"%-",round(roc2_test$ci[3],2),"%")))
print(Sys.time())

#Fig E1	
vioplot(ROSMAP.pre3, col = c(3,7), ylab = "Prob")

#Fig E2		
level3 <- used_id3_train
level3 <- factor(level3, levels = c(0, 1))

#Fig E3
roc3 <- roc(level3, ROSMAP.pre3[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc3 <- roc(level3, ROSMAP.pre3[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified = FALSE, plot = TRUE, percent = roc3$percent, col = 2)
sens.ci <- ci.se(roc3, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc3,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc3$ci[2],2),"%"),
paste("95% CI:",round(roc3$ci[1],2),"%-",round(roc3$ci[3],2),"%")))
print("ROC for test dataset 2 d")

#Fig4
roc3_test <- roc(used_id3_test, Test_roc3[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc3_test <- roc(used_id3_test, Test_roc3[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, percent = roc3_test$percent, col = 2)
sens.ci <- ci.se(roc3_test, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc3_test,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc3_test$ci[2],2),"%"),
paste("95% CI:",round(roc3_test$ci[1],2),"%-",round(roc3_test$ci[3],2),"%")))
print(Sys.time())

dev.off()


train.rf1_importance <- cbind(train.rf1$importance, train.rf1$importanceSD)
write.table(train.rf1_importance, file = "train.rf1_importance_ROSMAP.xls", quote = F, sep = "\t")

train.rf2_importance <- cbind(train.rf2$importance, train.rf2$importanceSD)
write.table(train.rf2_importance, file = "train.rf2_importance_ROSMAP.xls", quote = F, sep = "\t")

train.rf3_importance <- cbind(train.rf3$importance, train.rf3$importanceSD)
write.table(train.rf3_importance, file = "train.rf3_importance_ROSMAP.xls", quote = F, sep = "\t")

pdf("ROSMAP_importance.pdf", 6, 6)
plot(train.rf1_importance[,3], train.rf2_importance[,3], xlab = "Important in Typical", ylab = "Important in Intermediate", pch = 1, cex = 0, xlim = c(-0.0005, 0.0105))
text(train.rf1_importance[,3], train.rf2_importance[,3], label = rownames(train.rf1_importance))

plot(train.rf1_importance[,3], train.rf3_importance[,3], xlab = "Important in Typical", ylab = "Important in Atypical", pch = 1, cex = 0, xlim = c(-0.0005, 0.0105))
text(train.rf1_importance[,3], train.rf3_importance[,3], label = rownames(train.rf1_importance))

plot(train.rf2_importance[,3], train.rf3_importance[,3], xlab = "Important in Intermediate", ylab = "Important in Atypical", pch = 1, cex = 0, xlim = c(-0.0005, 0.0075))	#,
text(train.rf2_importance[,3], train.rf3_importance[,3], label = rownames(train.rf1_importance))
dev.off()