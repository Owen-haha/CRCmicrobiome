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

fudan_guangzhou <- fread("/hwfssz1/ST_HEALTH/P18Z10200N0127/yqin/crczz/fudan/modelling/data_ml/species_fudan_guangzhou_n660.tsv")
pub_spe <- fread("/hwfssz1/ST_HEALTH/P18Z10200N0127/yqin/crczz/fudan/modelling/data_ml/species_pub_data_n1260.tsv")

# prepare input data for siamacat
pub_spe_ra <- dcast(pub_spe, Taxon~sampleID, value.var = "abundance_ra") # sample in column
pub_spe_ra <- dt.to.matrix(pub_spe_ra)

pub_phe <- unique(pub_spe[, .(sampleID, study_condition, age, gender, tnm, dataset_name)]) # sample in row
pub_phe <- dt.to.matrix(pub_phe)
pub_phe <- as.data.frame(pub_phe)

# guangzhou
oGZ <- dcast(fudan_guangzhou[Group=="oCRC_Guangzhou"], Taxon~sampleID, value.var = "abundance_ra") # sample in column
oGZ <- dt.to.matrix(oGZ)

yGZ <- dcast(fudan_guangzhou[Group=="yCRC_Guangzhou"], Taxon~sampleID, value.var = "abundance_ra") # sample in column
yGZ <- dt.to.matrix(yGZ)

oFD <- dcast(fudan_guangzhou[Group=="oCRC_Fudan" | Group=="oControl"], Taxon~sampleID, value.var = "abundance_ra") # sample in column
oFD <- dt.to.matrix(oFD)
oFD_phe <- unique(fudan_guangzhou[, .(sampleID, Group, study_condition)])
oFD_phe <- as.data.frame(dt.to.matrix(oFD_phe))
oFD <- siamcat(feat = oFD, meta = oFD_phe, label='study_condition', case='CRC')

yFD <- dcast(fudan_guangzhou[Group=="yCRC_Fudan" | Group=="yControl"], Taxon~sampleID, value.var = "abundance_ra") # sample in column
yFD <- dt.to.matrix(yFD)
yFD_phe <- unique(fudan_guangzhou[, .(sampleID, Group, study_condition)])
yFD_phe <- as.data.frame(dt.to.matrix(yFD_phe))
yFD <- siamcat(feat = yFD, meta = yFD_phe, label='study_condition', case='CRC')

sc.obj <- siamcat(feat=pub_spe_ra, meta=pub_phe, label='study_condition', case='CRC')

sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 0.001)

sc.obj <- normalize.features(sc.obj, norm.method = "log.std",  
                             norm.param = list(log.n0 = 1e-06, n.p = 2, sd.min.q = 0.1, norm.margin = 1))

sc.obj <-  create.data.split(sc.obj, num.folds = 10, num.resample = 10)

sc.obj <- train.model(sc.obj, method = "lasso")

cat("cross validation of pub:\n")
sc.obj <- make.predictions(sc.obj)
# pred_matrix <- pred_matrix(sc.obj)
sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj)
sc.obj
save(sc.obj, file = "pub_lasso_model.siamcat")
model.interpretation.plot(sc.obj, fn.plot = 'interpretation_pub_lasso_model.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

cat("pub predict oGZ:\n")
predict_oGZ <- make.predictions(sc.obj, siamcat.holdout = siamcat(feat = oGZ))
predict_oGZ <- pred_matrix(predict_oGZ)
predict_oGZ <- predict_oGZ>0.5
summary(colSums(predict_oGZ)/nrow(predict_oGZ)) 

cat("pub predict yGZ:\n")
predict_yGZ <- make.predictions(sc.obj, siamcat.holdout = siamcat(feat = yGZ))
predict_yGZ <- pred_matrix(predict_yGZ)
predict_yGZ <- predict_yGZ>0.5
summary(colSums(predict_yGZ)/nrow(predict_yGZ)) 

cat("pub predict oFD:\n")
oFD <- make.predictions(sc.obj, siamcat.holdout = oFD)
oFD <-  evaluate.predictions(oFD)
model.evaluation.plot(oFD) # Average AUC: 0.788
oFD

cat("pub predict yFD:\n")
yFD <- make.predictions(sc.obj, siamcat.holdout = yFD)
yFD <-  evaluate.predictions(yFD)
model.evaluation.plot(yFD) # Average AUC: 0.796
yFD


