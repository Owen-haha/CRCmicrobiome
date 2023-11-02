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

oGZ <- dcast(fudan_guangzhou[Group=="oCRC_Guangzhou"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
oGZ <- dt.to.matrix(oGZ)

yGZ <- dcast(fudan_guangzhou[Group=="yCRC_Guangzhou"], ec_number~sampleID, value.var = "abundance_ra") # sample in column
yGZ <- dt.to.matrix(yGZ)

fd <- dcast(fudan_guangzhou[!grep("Guangzhou", Group)], ec_number~sampleID, value.var = "abundance_ra") # sample in column
fd <- dt.to.matrix(fd)
fd_phe <- unique(fudan_guangzhou[!grep("Guangzhou", Group), .(sampleID, study_condition, Group)])

############# trained on pub + Fudan, make prediction on oGZ, yGZ
pub_fd <- cbind(pub_ra, fd)
pub_fd_phe <- rbind(pub_phe, fd_phe)
pub_fd_phe <- as.data.frame(dt.to.matrix(pub_fd_phe))

sc.obj <- siamcat(feat=pub_fd, meta=pub_fd_phe, label="study_condition", case="CRC")
sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 1e-6)
sc.obj <- normalize.features(sc.obj, norm.method = "log.std",  
                             norm.param = list(log.n0 = 1e-09, n.p = 2, sd.min.q = 0.1, norm.margin = 1))

sc.obj <-  create.data.split(sc.obj, num.folds = 10, num.resample = 10)
sc.obj <- train.model(sc.obj, method = "lasso")

cat("cross validation of pub_fd:\n")
sc.obj <- make.predictions(sc.obj)
sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj) 
sc.obj

save(sc.obj, file = "pub_fd_lasso_model.siamcat")
model.interpretation.plot(sc.obj, fn.plot = 'interpretation_pub_fd_model.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

cat("pub_fd predict oGZ:\n")
predict_oGZ <- make.predictions(sc.obj, siamcat.holdout = siamcat(feat = oGZ))
predict_oGZ <- pred_matrix(predict_oGZ)
predict_oGZ <- predict_oGZ>0.5
summary(colSums(predict_oGZ)/nrow(predict_oGZ))

cat("pub_fd predict yGZ:\n")
predict_yGZ <- make.predictions(sc.obj, siamcat.holdout = siamcat(feat = yGZ))
predict_yGZ <- pred_matrix(predict_yGZ)
predict_yGZ <- predict_yGZ>0.5
summary(colSums(predict_yGZ)/nrow(predict_yGZ))

