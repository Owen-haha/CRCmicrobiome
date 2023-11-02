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
pub <- fread("/hwfssz1/ST_HEALTH/P18Z10200N0127/yqin/crczz/fudan/modelling/data_ml/level4ec_pub_data_n1260.tsv")

pub_ra <- dcast(pub, ec_number~sampleID, value.var = "abundance_ra") # sample in column
pub_ra <- dt.to.matrix(pub_ra)
pub_phe <- unique(pub[, .(sampleID, study_condition, Group=dataset_name)]) # sample in row

oFD <- dcast(fudan_guangzhou[Group=="oCRC_Fudan" | Group=="oControl"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
oFD <- dt.to.matrix(oFD)
oFD_phe <- unique(fudan_guangzhou[Group=="oCRC_Fudan" | Group=="oControl", .(sampleID, Group, study_condition)])
oFD_phe <- as.data.frame(dt.to.matrix(oFD_phe))
oFD <- siamcat(feat = oFD, meta = oFD_phe, label='study_condition', case='CRC')

yFD <- dcast(fudan_guangzhou[Group=="yCRC_Fudan" | Group=="yControl"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
yFD <- dt.to.matrix(yFD)
yFD_phe <- unique(fudan_guangzhou[Group=="yCRC_Fudan" | Group=="yControl", .(sampleID, Group, study_condition)])
yFD_phe <- as.data.frame(dt.to.matrix(yFD_phe))
yFD <- siamcat(feat = yFD, meta = yFD_phe, label='study_condition', case='CRC')

# pub + guangzhou
gz <- dcast(fudan_guangzhou[grep("Guangzhou", Group)], ec_number~sampleID, value.var = "abundance_ra") # sample in column
gz <- dt.to.matrix(gz)
gz_phe <- unique(fudan_guangzhou[grep("Guangzhou", Group), .(sampleID, study_condition, Group)])

pub_gz <- cbind(pub_ra, gz)
pub_gz_phe <- rbind(pub_phe, gz_phe)
pub_gz_phe <- as.data.frame(dt.to.matrix(pub_gz_phe))

sc.obj <- siamcat(feat=pub_gz, meta=pub_gz_phe, label="study_condition", case="CRC")
sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 1e-06)
sc.obj <- normalize.features(sc.obj, norm.method = "log.std",  
                             norm.param = list(log.n0 = 1e-09, n.p = 2, sd.min.q = 0.1, norm.margin = 1))

sc.obj <-  create.data.split(sc.obj, num.folds = 10, num.resample = 10)
sc.obj <- train.model(sc.obj, method = "lasso")

cat("cross validation of pub_gz:\n")
sc.obj <- make.predictions(sc.obj)
sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj) # Average AUC: 0.872
sc.obj

save(sc.obj, file = "pub_gz_lasso_model.siamcat")
model.interpretation.plot(sc.obj, fn.plot = 'interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

cat("pub_gz predict oFD:\n")
oFD <- make.predictions(sc.obj, siamcat.holdout = oFD)
oFD <-  evaluate.predictions(oFD)
oFD
model.evaluation.plot(oFD)
# recall rate
oFD_pred <- pred_matrix(oFD)
oFD_pred <- oFD_pred[grep("M_O", rownames(oFD_pred)), ]
oFD_pred <- oFD_pred>0.5
summary(colSums(oFD_pred)/nrow(oFD_pred)) # mean 0.8678

cat("pub_gz predict yFD:\n")
yFD <- make.predictions(sc.obj, siamcat.holdout = yFD)
yFD <-  evaluate.predictions(yFD)
yFD
model.evaluation.plot(yFD)
# recall rate
yFD_pred <- pred_matrix(yFD)
yFD_pred <- yFD_pred[grep("M_Y", rownames(yFD_pred)), ]
yFD_pred <- yFD_pred>0.5
summary(colSums(yFD_pred)/nrow(yFD_pred)) # mean 0.8466

