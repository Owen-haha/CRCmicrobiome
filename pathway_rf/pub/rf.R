library(data.table)
library(mlr3)
library(mlr3learners)
library(ggplot2)
library(mlr3viz)

set.seed(1234567)

# take first column as the rownames
dt.to.matrix <- function(x) {
  x <- as.data.frame(x)
  x <- x[!is.na(x[, 1]), ]
  rownames(x) <- x[,1]
  x <- as.matrix(x[,-1])
  x
}

fudan_guangzhou <- fread("/hwfssz1/ST_HEALTH/P18Z10200N0127/yqin/crczz/fudan/modelling/data_ml/pathway_fudan_guangzhou_n660.tsv")
pub_spe <- fread("/hwfssz1/ST_HEALTH/P18Z10200N0127/yqin/crczz/fudan/modelling/data_ml/pathway_pub_data_n1260.tsv")

# prepare data for mlr3
pub_spe <- dcast(pub_spe, study_condition + sampleID~Pathway, value.var = "abundance_cpm")
task_pub_spe <- as_task_classif(pub_spe[, -c("sampleID")], target = "study_condition", id = "publicData", positive = "CRC")

oFudan <- fudan_guangzhou[Group %in% c("oControl", "oCRC_Fudan")]
oFudan <- dcast(oFudan,  study_condition + sampleID~Pathway, value.var = "abundance_cpm")
yFudan <- fudan_guangzhou[Group %in% c("yControl", "yCRC_Fudan")]
yFudan <- dcast(yFudan,  study_condition + sampleID~Pathway, value.var = "abundance_cpm")

oGuangzhou <- fudan_guangzhou[Group == "oCRC_Guangzhou"] # N=293
oGuangzhou <- dcast(oGuangzhou,  study_condition + sampleID~Pathway, value.var = "abundance_cpm")
yGuangzhou <- fudan_guangzhou[Group == "yCRC_Guangzhou"] # N=167
yGuangzhou <- dcast(yGuangzhou,  study_condition + sampleID~Pathway, value.var = "abundance_cpm")

# default parameters of classif.ranger see https://mlr3learners.mlr-org.com/reference/mlr_learners_classif.ranger.html 
# trained on public data
rf_tree10000 <- lrn("classif.ranger", predict_type = "prob", num.trees = 10000, importance = "impurity")
rf_tree10000$train(task_pub_spe)

# make prediction
cat("\npub predict oFudan:\n")
predict_oFudan <- rf_tree10000$predict_newdata(oFudan)
predict_oFudan$score(msr("classif.auc")) # 0.8104

cat("\npub predict yFudan:\n")
predict_yFudan <- rf_tree10000$predict_newdata(yFudan)
predict_yFudan$score(msr("classif.auc")) # 0.7788

cat("\npub predict oGuangzhou:\n")
predict_oGuangzhou <- rf_tree10000$predict_newdata(oGuangzhou)
1- predict_oGuangzhou$score(msr("classif.ce")) # 0.894198

cat("\npub predict yGuangzhou:\n")
predict_yGuangzhou <- rf_tree10000$predict_newdata(yGuangzhou)
1- predict_yGuangzhou$score(msr("classif.ce")) # 0.8742515

cat("\nAUC of 10-fold cross validation for pub:\n")
cv_auc <- matrix(0, nrow = 10, ncol = 1)
for(i in 1:10){
   rf_cv <- resample(
     task = task_pub_spe,
     learner = rf_tree10000,
     resampling = rsmp("cv", folds = 10),
     store_models = FALSE )
  cv_auc[i, 1] <- rf_cv$aggregate(msr("classif.auc"))
}
cv_auc
summary(cv_auc[, 1])

# save prediction results
fwrite(cbind(oFudan$sampleID, as.data.table(predict_oFudan)), "pub_predict_oFD.txt", sep = "\t", quote = F)
fwrite(cbind(yFudan$sampleID, as.data.table(predict_yFudan)), "pub_predict_yFD.txt", sep = "\t", quote = F)
fwrite(cbind(oGuangzhou$sampleID, as.data.table(predict_oGuangzhou)), "pub_predict_oGZ.txt", sep = "\t", quote = F)
fwrite(cbind(yGuangzhou$sampleID, as.data.table(predict_yGuangzhou)), "pub_predict_yGZ.txt", sep = "\t", quote = F)

# visualization
roc_rf_cv <- autoplot(rf_cv, type = "roc")
ggsave(roc_rf_cv, file = "roc_cv_pub.png")

roc_oFudan <- autoplot(predict_oFudan, type = "roc")
ggsave(roc_oFudan, file = "roc_oFD.png")

roc_yFudan <- autoplot(predict_yFudan, type = "roc")
ggsave(roc_yFudan, file = "roc_yFD.png")

pdf("pub_predict_oGZ_distr.pdf")
hist(as.numeric(as.data.table(predict_oGuangzhou)$prob.CRC), ylab = "Number of samples", 
     xlab = "Predicted CRC probability", main = "oCRC_Guangzhou", breaks = 10, cex.lab = 1.2)
dev.off()

pdf("pub_predict_yGZ_distr.pdf")
hist(as.numeric(as.data.table(predict_yGuangzhou)$prob.CRC), ylab = "Number of samples", 
     xlab = "Predicted CRC probability", main = "yCRC_Guangzhou", breaks = 10, cex.lab = 1.2)
dev.off()

