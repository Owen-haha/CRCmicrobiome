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

fudan_guangzhou <- fread("/hwfssz1/ST_HEALTH/P18Z10200N0127/yqin/crczz/fudan/modelling/data_ml/level4ec_fudan_guangzhou_n660.tsv")
pub <- fread("/hwfssz1/ST_HEALTH/P18Z10200N0127/yqin/crczz/fudan/modelling/data_ml/level4ec_pub_data_n1260.tsv")

fudan_guangzhou[, ec_number:=paste0("EC", ec_number)]
pub[, ec_number:=paste0("EC", ec_number)]

# prepare data for mlr3
pub <- dcast(pub, study_condition + sampleID~ec_number, value.var = "abundance_cpm")

oFudan <- fudan_guangzhou[Group %in% c("oControl", "oCRC_Fudan")]
oFudan <- dcast(oFudan,  study_condition + sampleID~ec_number, value.var = "abundance_cpm")
yFudan <- fudan_guangzhou[Group %in% c("yControl", "yCRC_Fudan")]
yFudan <- dcast(yFudan,  study_condition + sampleID~ec_number, value.var = "abundance_cpm")

oGuangzhou <- fudan_guangzhou[Group == "oCRC_Guangzhou"] # N=293
oGuangzhou <- dcast(oGuangzhou,  study_condition + sampleID~ec_number, value.var = "abundance_cpm")
yGuangzhou <- fudan_guangzhou[Group == "yCRC_Guangzhou"] # N=167
yGuangzhou <- dcast(yGuangzhou,  study_condition + sampleID~ec_number, value.var = "abundance_cpm")

# train on 1262 public samples (600 CRC + 662 control) + 473 Guangzhou samples (all CRC)
# 1735 samples: 1073 CRC + 662 control
pub_gz <- rbind(pub, oGuangzhou)
pub_gz <- rbind(pub_gz, yGuangzhou)
task_pub_gz <- as_task_classif(pub_gz[, -c("sampleID")], target = "study_condition", id = "pub_gz", positive = "CRC")

# default parameters of classif.ranger see https://mlr3learners.mlr-org.com/reference/mlr_learners_classif.ranger.html 
# trained on public data
rf_tree10000 <- lrn("classif.ranger", predict_type = "prob", num.trees = 10000, importance = "impurity")
rf_tree10000$train(task_pub_gz)

# make prediction
cat("\npub_gz predict oFudan:\n")
predict_oFudan <- rf_tree10000$predict_newdata(oFudan)
predict_oFudan$score(msr("classif.auc")) 

cat("\npub_gz predict yFudan:\n")
predict_yFudan <- rf_tree10000$predict_newdata(yFudan)
predict_yFudan$score(msr("classif.auc")) 

cat("\nAUC of 10-fold cross validation for pub_gz:\n")
cv_auc <- matrix(0, nrow = 10, ncol = 1)
for(i in 1:10){
   rf_cv <- resample(
     task = task_pub_gz,
     learner = rf_tree10000,
     resampling = rsmp("cv", folds = 10),
     store_models = FALSE )
  cv_auc[i, 1] <- rf_cv$aggregate(msr("classif.auc"))
}
cv_auc
summary(cv_auc[, 1])

# save prediction results
fwrite(cbind(oFudan$sampleID, as.data.table(predict_oFudan)), "pub_gz_predict_oFD.txt", sep = "\t", quote = F)
fwrite(cbind(yFudan$sampleID, as.data.table(predict_yFudan)), "pub_gz_predict_yFD.txt", sep = "\t", quote = F)

# visualization
roc_cv <- autoplot(rf_cv, type = "roc")
ggsave(roc_cv, file = "roc_cv_pub_gz.png")

roc_oFudan <- autoplot(predict_oFudan, type = "roc")
ggsave(roc_oFudan, file = "roc_oFD.png")

roc_yFudan <- autoplot(predict_yFudan, type = "roc")
ggsave(roc_yFudan, file = "roc_yFD.png")

