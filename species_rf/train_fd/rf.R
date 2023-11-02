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

fudan_guangzhou <- fread("data_ml/species_fudan_guangzhou_n660.tsv")

# prepare data for mlr3
oFudan <- fudan_guangzhou[Group %in% c("oControl", "oCRC_Fudan")]
oFudan <- dcast(oFudan,  study_condition + sampleID~Taxon, value.var = "abundance_count")
oFudan_task <- as_task_classif(oFudan[, sampleID:=NULL], target = "study_condition", id = "oFudan", positive = "CRC")

yFudan <- fudan_guangzhou[Group %in% c("yControl", "yCRC_Fudan")]
yFudan <- dcast(yFudan,  study_condition + sampleID~Taxon, value.var = "abundance_count")
yFudan_task <- as_task_classif(yFudan[, sampleID:=NULL], target = "study_condition", id = "yFudan", positive = "CRC")

oGuangzhou <- fudan_guangzhou[Group == "oCRC_Guangzhou"] # N=293
oGuangzhou <- dcast(oGuangzhou,  study_condition + sampleID~Taxon, value.var = "abundance_count")
yGuangzhou <- fudan_guangzhou[Group == "yCRC_Guangzhou"] # N=167
yGuangzhou <- dcast(yGuangzhou,  study_condition + sampleID~Taxon, value.var = "abundance_count")

# default parameters of classif.ranger see https://mlr3learners.mlr-org.com/reference/mlr_learners_classif.ranger.html 
# trained on public data
rf_tree10000 <- lrn("classif.ranger", predict_type = "prob", num.trees = 10000, importance = "impurity")

# cross-validation in oFudan, repeat 100 times 10-fold cross validation to get a fair evaluation
cv_auc_oFudan <- matrix(0, nrow = 100, ncol = 1)
for(i in 1:100){
   rf_oFudan <- resample(
     task = oFudan_task,
     learner = rf_tree10000,
     resampling = rsmp("cv", folds = 10),
     store_models = FALSE )
  cv_auc_oFudan[i, 1] <- rf_oFudan$aggregate(msr("classif.auc"))
}
cat("\nAUC of 100 times 10-fold cross validation for oFD:\n ")
cv_auc_oFudan
summary(cv_auc_oFudan[, 1])

# trained on oFudan
rf_tree10000$train(oFudan_task)

cat("\noFD predict yFD:\n")
predict_yFudan <- rf_tree10000$predict_newdata(yFudan)
predict_yFudan$score(msr("classif.auc")) # 0.7724

cat("\noFD predict oGuangzhou:\n")
predict_oGuangzhou <- rf_tree10000$predict_newdata(oGuangzhou)
1 - predict_oGuangzhou$score(msr("classif.ce")) # 0.8122867

cat("\noFD predict yGuangzhou:\n")
predict_yGuangzhou <- rf_tree10000$predict_newdata(yGuangzhou)
1 - predict_yGuangzhou$score(msr("classif.ce")) # 0.754491

# save prediction results
fwrite(cbind(yFudan$sampleID, as.data.table(predict_yFudan)), "oFD_predict_yFD.txt", sep = "\t", quote = F)
fwrite(cbind(oGuangzhou$sampleID, as.data.table(predict_oGuangzhou)), "oFD_predict_oGZ.txt", sep = "\t", quote = F)
fwrite(cbind(yGuangzhou$sampleID, as.data.table(predict_yGuangzhou)), "oFD_predict_yGZ.txt", sep = "\t", quote = F)


# cross-validation in yFudan, repeat 100 times 10-fold cross validation to get a fair evaluation
cv_auc_yFudan <- matrix(0, nrow = 100, ncol = 1)
for(i in 1:100){
   rf_yFudan <- resample(
     task = yFudan_task,
     learner = rf_tree10000,
     resampling = rsmp("cv", folds = 10),
     store_models = FALSE )
  cv_auc_yFudan[i, 1] <- rf_yFudan$aggregate(msr("classif.auc"))
}
cat("\nAUC of 100 times 10-fold cross validation for yFD:\n ")
cv_auc_yFudan
summary(cv_auc_yFudan[, 1])

# trained on yFudan
rf_tree10000$train(yFudan_task)

cat("\nyFD predict oFD:\n")
predict_oFudan <- rf_tree10000$predict_newdata(oFudan)
predict_oFudan$score(msr("classif.auc")) # 0.7788

cat("\nyFD predict oGuangzhou:\n")
predict_oGuangzhou <- rf_tree10000$predict_newdata(oGuangzhou)
1 - predict_oGuangzhou$score(msr("classif.ce")) # 0.5972696

cat("\nyFD predict yGuangzhou:\n")
predict_yGuangzhou <- rf_tree10000$predict_newdata(yGuangzhou)
1 - predict_yGuangzhou$score(msr("classif.ce")) # 0.6227545

# save prediction results
fwrite(cbind(oFudan$sampleID, as.data.table(predict_oFudan)), "yFD_predict_oFD.txt", sep = "\t", quote = F)
fwrite(cbind(oGuangzhou$sampleID, as.data.table(predict_oGuangzhou)), "yFD_predict_oGZ.txt", sep = "\t", quote = F)
fwrite(cbind(yGuangzhou$sampleID, as.data.table(predict_yGuangzhou)), "yFD_predict_yGZ.txt", sep = "\t", quote = F)

