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
pub_spe <- fread("data_ml/species_pub_data_n1260.tsv")

# prepare data for mlr3
pub_spe <- dcast(pub_spe, study_condition + sampleID~Taxon, value.var = "abundance_count")

oFudan <- fudan_guangzhou[Group %in% c("oControl", "oCRC_Fudan")]
oFudan <- dcast(oFudan,  study_condition + sampleID~Taxon, value.var = "abundance_count")
yFudan <- fudan_guangzhou[Group %in% c("yControl", "yCRC_Fudan")]
yFudan <- dcast(yFudan,  study_condition + sampleID~Taxon, value.var = "abundance_count")

oGuangzhou <- fudan_guangzhou[Group == "oCRC_Guangzhou"] # N=293
oGuangzhou <- dcast(oGuangzhou,  study_condition + sampleID~Taxon, value.var = "abundance_count")
yGuangzhou <- fudan_guangzhou[Group == "yCRC_Guangzhou"] # N=167
yGuangzhou <- dcast(yGuangzhou,  study_condition + sampleID~Taxon, value.var = "abundance_count")

# train on 1262 public samples (600 CRC + 662 control) + 200 Fudan samples
# 1462 samples: 700 CRC + 762 control
pub_fd <- rbind(pub_spe, oFudan)
pub_fd <- rbind(pub_fd, yFudan)
task_pub_fd <- as_task_classif(pub_fd[, -c("sampleID")], target = "study_condition", id = "pub_fd", positive = "CRC")

# default parameters of classif.ranger see https://mlr3learners.mlr-org.com/reference/mlr_learners_classif.ranger.html 
# trained on public data
rf_tree10000 <- lrn("classif.ranger", predict_type = "prob", num.trees = 10000, importance = "impurity")
rf_tree10000$train(task_pub_fd)

cat("\npub predict oGuangzhou:\n")
predict_oGuangzhou <- rf_tree10000$predict_newdata(oGuangzhou)
1- predict_oGuangzhou$score(msr("classif.ce")) 

cat("\npub predict yGuangzhou:\n")
predict_yGuangzhou <- rf_tree10000$predict_newdata(yGuangzhou)
1- predict_yGuangzhou$score(msr("classif.ce")) 

cat("\nAUC of 10-fold cross validation for pub_fd:\n")
cv_auc <- matrix(0, nrow = 10, ncol = 1)
for(i in 1:10){
   rf_cv <- resample(
     task = task_pub_fd,
     learner = rf_tree10000,
     resampling = rsmp("cv", folds = 10),
     store_models = FALSE )
  cv_auc[i, 1] <- rf_cv$aggregate(msr("classif.auc"))
}
cv_auc
summary(cv_auc[, 1])

# save prediction results
fwrite(cbind(oGuangzhou$sampleID, as.data.table(predict_oGuangzhou)), "pub_fd_predict_oGZ.txt", sep = "\t", quote = F)
fwrite(cbind(yGuangzhou$sampleID, as.data.table(predict_yGuangzhou)), "pub_fd_predict_yGZ.txt", sep = "\t", quote = F)

# visualization
roc_cv <- autoplot(rf_cv, type = "roc")
ggsave(roc_cv, file = "roc_cv_pub.png")

pdf("pub_fd_predict_oGZ_distr.pdf")
hist(as.numeric(as.data.table(predict_oGuangzhou)$prob.CRC), ylab = "Number of samples", 
     xlab = "Predicted CRC probability", main = "oCRC_Guangzhou", breaks = 10, cex.lab = 1.2)
dev.off()

pdf("pub_fd_predict_yGZ_distr.pdf")
hist(as.numeric(as.data.table(predict_yGuangzhou)$prob.CRC), ylab = "Number of samples", 
     xlab = "Predicted CRC probability", main = "yCRC_Guangzhou", breaks = 10, cex.lab = 1.2)
dev.off()

