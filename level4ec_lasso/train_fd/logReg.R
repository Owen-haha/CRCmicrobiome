library(data.table)
library("SIAMCAT") # https://siamcat.embl.de/articles/SIAMCAT_vignette.html 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02306-1 

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

oGZ <- dcast(fudan_guangzhou[Group=="oCRC_Guangzhou"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
oGZ <- dt.to.matrix(oGZ)

yGZ <- dcast(fudan_guangzhou[Group=="yCRC_Guangzhou"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
yGZ <- dt.to.matrix(yGZ)

oFD <- dcast(fudan_guangzhou[Group=="oCRC_Fudan" | Group=="oControl"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
oFD <- dt.to.matrix(oFD)
oFD_phe <- unique(fudan_guangzhou[, .(sampleID, Group, study_condition)])
oFD_phe <- as.data.frame(dt.to.matrix(oFD_phe))
oFD <- siamcat(feat = oFD, meta = oFD_phe, label='study_condition', case='CRC')

yFD <- dcast(fudan_guangzhou[Group=="yCRC_Fudan" | Group=="yControl"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
yFD <- dt.to.matrix(yFD)
yFD_phe <- unique(fudan_guangzhou[, .(sampleID, Group, study_condition)])
yFD_phe <- as.data.frame(dt.to.matrix(yFD_phe))
yFD <- siamcat(feat = yFD, meta = yFD_phe, label='study_condition', case='CRC')

############# ############# ############# ############# ############# ############# ############# 
############# train on oFD + yFD, make prediction on oGZ, yGZ
fd <- dcast(fudan_guangzhou[!grep("Guangzhou", Group)], ec_number~sampleID, value.var = "abundance_ra") # sample in column
fd <- dt.to.matrix(fd)
fd_phe <- unique(fudan_guangzhou[!grep("Guangzhou", Group), .(sampleID, study_condition, Group)])
fd_phe <- as.data.frame(dt.to.matrix(fd_phe))

sc.obj <- siamcat(feat = fd, meta = fd_phe, label="study_condition", case="CRC")
sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 1e-06)
sc.obj <- normalize.features(sc.obj, norm.method = "log.std",  
                             norm.param = list(log.n0 = 1e-09, n.p = 2, sd.min.q = 0.1, norm.margin = 1))

sc.obj <-  create.data.split(sc.obj, num.folds = 10, num.resample = 10)
sc.obj <- train.model(sc.obj, method = "lasso")

# nested feature selection
# sc.obj <- train.model(sc.obj, method = "lasso", perform.fs = TRUE, param.fs = list(no_features = 100, method = "AUC", direction="absolute"))

sc.obj <- make.predictions(sc.obj)
sc.obj <-  evaluate.predictions(sc.obj)
cat("cross-validation of fd:\n")
sc.obj
#model.evaluation.plot(sc.obj) # Average AUC: 0.73

save(sc.obj, file = "fd_lasso_model.siamcat")

cat("fd_predict_oGZ:\n")
fd_predict_oGZ <- make.predictions(sc.obj, siamcat.holdout = siamcat(feat = oGZ))
fd_predict_oGZ <- pred_matrix(fd_predict_oGZ)
fd_predict_oGZ <- fd_predict_oGZ>0.5
summary(colSums(fd_predict_oGZ)/nrow(fd_predict_oGZ)) # mean 0.8135

cat("fd_predict_yGZ:\n")
fd_predict_yGZ <- make.predictions(sc.obj, siamcat.holdout = siamcat(feat = yGZ))
fd_predict_yGZ <- pred_matrix(fd_predict_yGZ)
fd_predict_yGZ <- fd_predict_yGZ>0.5
summary(colSums(fd_predict_yGZ)/nrow(fd_predict_yGZ)) # mean 0.8128

############# trained on oFD, make prediction on yFD and oGZ, yGZ 
oFD_model <- filter.features(oFD, filter.method = 'abundance', cutoff = 1e-06)
oFD_model <- normalize.features(oFD_model, norm.method = "log.std",  
                                norm.param = list(log.n0 = 1e-09, n.p = 2, sd.min.q = 0.1, norm.margin = 1))

oFD_model <-  create.data.split(oFD_model, num.folds = 10, num.resample = 100)
oFD_model <- train.model(oFD_model, method = "lasso")

# nested feature selection
# oFD_model <- train.model(oFD_model, method = "lasso", perform.fs = TRUE, param.fs = list(no_features = 100, method = "AUC", direction="absolute"))

cat("cross-validation of oFD_model:\n")
oFD_model <- make.predictions(oFD_model)
oFD_model <-  evaluate.predictions(oFD_model)
oFD_model
save(oFD_model, file = "oFD_lasso_model.siamcat")

cat("oFD_predict_yFD:\n")
oFD_predict_yFD <- make.predictions(oFD_model, siamcat.holdout = yFD)
oFD_predict_yFD <- evaluate.predictions(oFD_predict_yFD) # Average AUC: 0.792
oFD_predict_yFD

cat("oFD_predict_oGZ:\n")
oFD_predict_oGZ <- make.predictions(oFD_model, siamcat.holdout = siamcat(feat = oGZ))
oFD_predict_oGZ <- pred_matrix(oFD_predict_oGZ)
oFD_predict_oGZ <- oFD_predict_oGZ>0.5
summary(colSums(oFD_predict_oGZ)/nrow(oFD_predict_oGZ)) # mean 0.7051, median 0.7201

cat("oFD_predict_yGZ:\n")
oFD_predict_yGZ <- make.predictions(oFD_model, siamcat.holdout = siamcat(feat = yGZ))
oFD_predict_yGZ <- pred_matrix(oFD_predict_yGZ)
oFD_predict_yGZ <- oFD_predict_yGZ>0.5
summary(colSums(oFD_predict_yGZ)/nrow(oFD_predict_yGZ)) # mean 0.6573, median 0.6766


############# trained on yFD, make prediction on yFD and oGZ, yGZ 
yFD_model <- filter.features(yFD, filter.method = 'abundance', cutoff = 1e-06)
yFD_model <- normalize.features(yFD_model, norm.method = "log.std",  
                                norm.param = list(log.n0 = 1e-09, n.p = 2, sd.min.q = 0.1, norm.margin = 1))

yFD_model <-  create.data.split(yFD_model, num.folds = 10, num.resample = 100)
yFD_model <- train.model(yFD_model, method = "lasso")

# nested feature selection
# yFD_model <- train.model(yFD_model, method = "lasso", perform.fs = TRUE, param.fs = list(no_features = 100, method = "AUC", direction="absolute"))

cat("cross-validation of yFD_model:\n")
yFD_model <- make.predictions(yFD_model)
yFD_model <-  evaluate.predictions(yFD_model)
yFD_model
save(yFD_model, file = "yFD_lasso_model.siamcat")

cat("yFD_predict_oFD:\n")
yFD_predict_oFD <- make.predictions(yFD_model, siamcat.holdout = oFD)
yFD_predict_oFD <- evaluate.predictions(yFD_predict_oFD) # Average AUC: 0.774
yFD_predict_oFD

cat("yFD_predict_oGZ:\n")
yFD_predict_oGZ <- make.predictions(yFD_model, siamcat.holdout = siamcat(feat = oGZ))
yFD_predict_oGZ <- pred_matrix(yFD_predict_oGZ)
yFD_predict_oGZ <- yFD_predict_oGZ>0.5
summary(colSums(yFD_predict_oGZ)/nrow(yFD_predict_oGZ)) # mean 0.7368, median 0.7509

cat("yFD_predict_yGZ:\n")
yFD_predict_yGZ <- make.predictions(yFD_model, siamcat.holdout = siamcat(feat = yGZ))
yFD_predict_yGZ <- pred_matrix(yFD_predict_yGZ)
yFD_predict_yGZ <- yFD_predict_yGZ>0.5
summary(colSums(yFD_predict_yGZ)/nrow(yFD_predict_yGZ)) # mean 0.7302, median 0.7335

model.interpretation.plot(sc.obj, fn.plot = 'interpretation_fd_model.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')
model.interpretation.plot(yFD_model, fn.plot = 'interpretation_yFD_model.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')
model.interpretation.plot(oFD_model, fn.plot = 'interpretation_oFD_model.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

