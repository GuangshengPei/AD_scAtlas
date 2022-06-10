library(reshape2)
library(vegan)
library(randomForest)
library(pROC)

MSBB_clean <- na.omit(MSBB[DEG_union,])
dim(MSBB_clean)

MSBB_clinical_clean <- MSBB_clinical[,c("tissue", "sex", "race", "CDR")]
MSBB_clinical_clean$status[MSBB_clinical_clean$CDR >= 1] <- "AD"
MSBB_clinical_clean$status[MSBB_clinical_clean$CDR == 0] <- "CTR"
MSBB_clinical_clean$status[MSBB_clinical_clean$CDR == 0.5] <- "MCI"
MSBB_clinical_clean$status[is.na(MSBB_clinical_clean$CDR)] <- "CTR"
MSBB_clinical_clean$CDR[MSBB_clinical_clean$CDR >= 1] <- 2
MSBB_clinical_clean$CDR[MSBB_clinical_clean$CDR == 0] <- 0
MSBB_clinical_clean$CDR[MSBB_clinical_clean$CDR == 0.5] <- 1
MSBB_clinical_clean$CDR[is.na(MSBB_clinical_clean$CDR)] <- 0

dim(MSBB_clean)

# two randomForest model for AD and MCI
dim(MSBB_clean)
i = 99
set.seed(i)
exclude_id1 <- which(MSBB_clinical_clean$CDR == 1)
exclude_id2 <- which(MSBB_clinical_clean$CDR == 2)
used_id1 <- MSBB_clinical_clean$CDR[-exclude_id1]
used_id2 <- MSBB_clinical_clean$CDR[-exclude_id2]
used_id1 <- factor(used_id1, levels = c(0, 2))
used_id2 <- factor(used_id2, levels = c(0, 1))

MSBB_clean_1 <- MSBB_clean[,-exclude_id1]
MSBB_clean_2 <- MSBB_clean[,-exclude_id2]

dim(MSBB_clean_1)
dim(MSBB_clean_2)

i = 99
set.seed(i)
test_id1 <- sample(ncol(MSBB_clean_1), 0.1 * ncol(MSBB_clean_1))
used_id1_train <- used_id1[-test_id1]
used_id1_test <- used_id1[test_id1]

test_id2 <- sample(ncol(MSBB_clean_2), 0.1 * ncol(MSBB_clean_2))
used_id2_train <- used_id2[-test_id2]
used_id2_test <- used_id2[test_id2]

used_id1 = as.vector(used_id1)
used_id1[used_id1 == 2] <- 1
used_id1 = factor(used_id1)

table(used_id1)
table(used_id2)

length(used_id1_train)
length(used_id2_train)

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

#------------------------------------------------
# AD vs Control
Sys.time()
result1 <- replicate(10, rfcv1(t(MSBB_clean_1), used_id1, cv.fold = 10, step = 0.9), simplify = FALSE)
error1.cv <- sapply(result1, "[[", "error.cv")
Sys.time()

matplot(result1[[1]]$n.var, cbind(rowMeans(error1.cv), error1.cv), type="l", lwd=c(2, rep(1, ncol(error1.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 25, col="pink", lwd = 2)

matplot(result1[[1]]$n.var, cbind(rowMeans(error1.cv), error1.cv), type="l", lwd=c(2, rep(1, ncol(error1.cv))), col=c(1,"grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey"), lty=1, log="x", xlab="Number of variables", ylab="CV Error")

save(result1, error1.cv, file = "result.temp1.rda")
# MCI vs Control
Sys.time()
result2 <- replicate(10, rfcv1(t(MSBB_clean_2), used_id2, cv.fold = 10, step = 0.9), simplify = FALSE)
error2.cv <- sapply(result2, "[[", "error.cv")
Sys.time()

matplot(result2[[1]]$n.var, cbind(rowMeans(error2.cv), error2.cv), type="l", lwd=c(2, rep(1, ncol(error2.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 25, col="pink", lwd = 2)

matplot(result2[[1]]$n.var, cbind(rowMeans(error2.cv), error2.cv), type="l", lwd=c(2, rep(1, ncol(error2.cv))), col=c(1,"grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey"), lty=1, log="x", xlab="Number of variables", ylab="CV Error")

save(result2, error2.cv, file = "result.temp2.rda")

# Importance analysis
load("result.temp1.rda")
load("result.temp2.rda")

matplot(result1[[1]]$n.var, cbind(rowMeans(error1.cv), error1.cv), type="l", lwd=c(2, rep(1, ncol(error1.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 15, col="pink", lwd = 2)

matplot(result2[[1]]$n.var, cbind(rowMeans(error2.cv), error2.cv), type="l", lwd=c(2, rep(1, ncol(error2.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 40, col="pink", lwd = 2)

roc_list1 <- list()
roc_list2 <- list()

train.rf1_importance <- cbind(train.rf1$importance, train.rf1$importanceSD)
write.table(train.rf1_importance, file = "train.rf1_importance.xls", quote = F, sep = "\t")

train.rf2_importance <- cbind(train.rf2$importance, train.rf2$importanceSD)
write.table(train.rf2_importance, file = "train.rf2_importance.xls", quote = F, sep = "\t")

pdf("MSBB_importance.pdf", 6, 6)
plot(train.rf1_importance[,3], train.rf2_importance[,3], xlab = "Important in AD", ylab = "Important in MCI", pch = 1, cex = 0, xlim = c(-0.0001, 0.0032))
text(train.rf1_importance[,3], train.rf2_importance[,3], label = rownames(train.rf1_importance))
dev.off()

cor.test(train.rf1_importance[,3], train.rf2_importance[,3])

# Best model #==================================================

train.rf1 <- randomForest(t(MSBB_clean_1[,-test_id1]), used_id1_train, importance = TRUE, mtry = 15)	
train.rf2 <- randomForest(t(MSBB_clean_2[,-test_id2]), used_id2_train, importance = TRUE, mtry = 40)	

# Model evaluation

pdf(paste0("Random_forrest_model.pdf"), 18, 6)
par(mfrow = c(2,6))

#Fig A1		Feature selection
matplot(result1[[1]]$n.var, cbind(rowMeans(error1.cv), error1.cv), type="l", lwd=c(2, rep(1, ncol(error1.cv))), col=c(1,"grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey", "grey"), lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 20, col="pink", lwd = 2)

MSBB.pre1 <- predict(train.rf1, type = "prob")
MSBB.pre2 <- predict(train.rf2, type = "prob")
MSBB.pre1_res <- as.data.frame(cbind(MSBB.pre1, used_id1_train))
MSBB.pre2_res <- as.data.frame(cbind(MSBB.pre2, used_id2_train))

library(vioplot)
colnames(MSBB.pre1) = c("CTR", "AD")
colnames(MSBB.pre2) = c("CTR", "MCI")

#Fig A2		Prediction
vioplot(MSBB.pre1, col = c(3,2), ylab = "Prob")

#Fig A3		Training AUC
level1 <- used_id1_train
level1 <- factor(level1, levels = c(0, 2))

roc1 <- roc(level1, MSBB.pre1[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc1 <- roc(level1, MSBB.pre1[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified = FALSE, plot = TRUE, percent = roc1$percent, col = 2)
sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc1,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc1$ci[2],2),"%"),
paste("95% CI:",round(roc1$ci[1],2),"%-",round(roc1$ci[3],2),"%")))
print("ROC for test dataset 1 d")

#Fig A4		Training comparison
col1 = MSBB.pre1_res[,3]
col1[col1 == "1"] <- "green"
col1[col1 == "2"] <- "red"
plot(rank(MSBB.pre1[,2], na.last = TRUE), MSBB.pre1[,2], col = col1, pch=16, xlab="Samples", ylab="Probability of AD")
abline(h=0.7)

#Fig A5		Test AUC
roc1_test <- roc(used_id1_test, Test_roc1[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc1_test <- roc(used_id1_test, Test_roc1[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, percent = roc1_test$percent, col = 2)
sens.ci <- ci.se(roc1_test, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc1_test,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc1_test$ci[2],2),"%"),
paste("95% CI:",round(roc1_test$ci[1],2),"%-",round(roc1_test$ci[3],2),"%")))
print(Sys.time())

#Fig A6		Test comparison
Test_roc1 <- predict(train.rf1, t(MSBB_clean_1[,test_id1]), type = "prob")
col1 = used_id1_test
col1 = as.vector(col1)
col1[col1 == "0"] <- "green"
col1[col1 == "1"] <- "red"
Test_roc1[,2] = Test_roc1[,2] + sample(nrow(Test_roc1))*1e-8
plot(rank(Test_roc1[,2], na.last = TRUE), Test_roc1[,2], col = col1, pch=16, xlab="Samples", ylab="Probability of AD")
abline(h=0.7)

# Control vs. MCI
#Fig B1		Feature selection
matplot(result2[[1]]$n.var, cbind(rowMeans(error2.cv), error2.cv), type="l", lwd=c(2, rep(1, ncol(error2.cv))), col=c(1,"grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey", "grey"), lty=1, log="x", xlab="Number of variables", ylab="CV Error")
abline(v = 50, col="pink", lwd = 2)

#Fig B2		Prediction
vioplot(MSBB.pre2, col = c(3,4), ylab = "Prob")

#Fig B3		Training AUC
level2 <- used_id2_train
level2 <- factor(level2, levels = c(0, 1))

roc2 <- roc(level2, MSBB.pre2[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc2 <- roc(level2, MSBB.pre2[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified = FALSE, plot = TRUE, percent = roc2$percent, col = 2)
sens.ci <- ci.se(roc2, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc2,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc2$ci[2],2),"%"),
paste("95% CI:",round(roc2$ci[1],2),"%-",round(roc2$ci[3],2),"%")))
print("ROC for test dataset 2 d")

#Fig B4		Training comparison
col2 = MSBB.pre2_res[,3]
col2[col2 == "1"] <- "green"
col2[col2 == "2"] <- "blue"
plot(rank(MSBB.pre2[,2], na.last = TRUE), MSBB.pre2[,2], col = col2, pch=16, xlab="Samples", ylab="Probability of AD")
abline(h=0.5)

#Fig B5		Test AUC
roc2_test <- roc(used_id2_test, Test_roc2[,2], percent = TRUE, partial.auc.correct=TRUE, ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
roc2_test <- roc(used_id2_test, Test_roc2[,2], ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, percent = roc2_test$percent, col = 2)
sens.ci <- ci.se(roc2_test, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc2_test,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc2_test$ci[2],2),"%"),
paste("95% CI:",round(roc2_test$ci[1],2),"%-",round(roc2_test$ci[3],2),"%")))
print(Sys.time())

#Fig B6		Test comparison
Test_roc2 <- predict(train.rf2, t(MSBB_clean_2[,test_id2]), type = "prob")
col2 = used_id2_test
col2 = as.vector(col2)
col2[col2 == "0"] <- "green"
col2[col2 == "1"] <- "blue"
Test_roc2[,2] = Test_roc2[,2] + sample(nrow(Test_roc2))*1e-8
plot(rank(Test_roc2[,2], na.last = TRUE), Test_roc2[,2], col = col2, pch=16, xlab="Samples", ylab="Probability of AD")
abline(h=0.5)

dev.off()