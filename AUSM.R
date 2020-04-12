# title         : Zambia Thesis
# purpose       : Updating legasy soil map by using digital soil mapping techiques
# producers     : Jeff Chalwe
# address       : Moscow State Agricltural University 

list.of.packages <- c("plyr", "parallel", "doParallel", "plotKML", "raster", "rgdal", "ranger", "caret",
                      "snowfall", "doParallel", "parallelMap", "sp", "mlr", "Boruta", "entropy")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) installed.packages(new.packages, dependencies=TRUE)

library(sp)
library(raster)
library(ranger)
library(caret)
library(rgdal)
library(plotKML)
library(parallel)
library(snowfall)
library(doParallel)
library(parallelMap)
library(entropy)

crs <- CRS('+init=epsg:32735')
points <- read.table("./training_points/points.txt", header=T, sep="\t")
coordinates(points) <- ~X+Y
points@proj4string <- crs

grid.list <- list.files(path="./Covariats/stack_520m/", pattern="*.tif")
grid.list <- grid.list[-unlist(sapply(c("CD"), function(x) {
grep(x, grid.list)
}))]
sfInit(parallel=TRUE, cpus=parallel::detectCores())
sfExport("grid.list", "points")
sfLibrary(rgdal)
sfLibrary(raster)
# overlay points & grids in parallel
ov.cov <- snowfall::sfClusterApplyLB(grid.list, function(i) {
  raster::extract(raster(paste0("./Covariats/stack_520m/", i)), points)
})
snowfall::sfStop()

ov.cov <- data.frame(ov.cov)
names(ov.cov) <- grid.list
reg.matrix <- cbind(as.data.frame(points), ov.cov)
str(reg.matrix)

saveRDS(reg.matrix, file=".regression_matrix.rds")

# layers statistics
stat.ov <- lapply(reg.matrix[,paste0(grid.list)], function(i) {
  data.frame(t(as.vector(summary(i))))
})
stat.ov <- dplyr::bind_rows(stat.ov)
names(stat.ov) <- c("MIN", "Q1", "MEDIAN", "MEAN", "Q3", "MAX", "NA")
stat.ov$layer_name <- basename(grid.list)
str(stat.ov)
write.csv(stat.ov, "./Results/area_weighted/cov_stats.csv")

fm <- as.formula(paste("Class ~ ", paste(grid.list, collapse="+")))
fm
sel.na <- complete.cases(reg.matrix[,all.vars(fm)])
summary(sel.na)
# some 7 points have missing values

# random forest
fit <- ranger::ranger(fm, reg.matrix[sel.na, ], num.trees=150, seed=1,
                      probability=T,
                      importance="impurity")
fit

# variable importance
variables <- as.list(ranger::importance(fit))
print(t(data.frame(variables[order(unlist(variables), decreasing=TRUE)[1:13]])))
# variable importance plot with "vip" package
vip::vip(fit, num_features=18L, width=0.55, bar=F)

saveRDS(fit, "./ranger_fit.rds")

r.file = "./Results/area_weighted/Second_level/ranger_fit.txt"
cat("Results of model fitting 'randomForest':\n", file=r.file)
cat("\n", file=r.file, append=TRUE)
cat(paste("Variable:", all.vars(fm)), file=r.file, append=TRUE)
cat("\n", file=r.file, append=TRUE)
sink(file=r.file, append=TRUE, type="output")
cat("\n Random forest model:", file=r.file, append=TRUE)
print(fit)
cat("\n Variable importance:\n", file=r.file, append=TRUE)
print(t(data.frame(variables[order(unlist(variables), decreasing=TRUE)[1:17]])))
sink()

# make predictions on random forest model
stack <- raster::stack(paste0("./Covariats/stack_520m/", grid.list))
names(stack) <- basename(grid.list)
s <- as(as(stack, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
sel.pr <- complete.cases(s@data)
# we can produce predictions of probabilities per class by running
pred_soil <- predict(fit, s@data[sel.pr, ])

# most probable predominant class
s[sel.pr,"Class"] <- as.factor(apply(pred_soil$predictions,1,which.max))
s$Class <- as.numeric(as.character(s$Class))
tax <- attr(pred_soil$predictions, "dimnames")[[2]]
writeGDAL(s["Class"], "./Results/test.tif",
          type="Int32", mvFlag=-99999, catNames=list(tax))

# calculate uncertainty
# https://github.com/ISRICWorldSoil/SoilGrids250m/blob/master/grids/models/uncertainty_SoilGrids.R
# derive Scaled Shannon Entropy (100 is a maximum error; 0 is perfect prediction)
v <- unlist(apply(pred_soil$predictions, 1, FUN=function(x){entropy.empirical(unlist(x))}))
s[sel.pr,"SSI"] <- round(v/entropy.empirical(rep(1/length(tax), length(tax)))*100)
# spplot((s["SSI"]), col.regions=SAGA_pal[[1]])
writeGDAL(s["SSI"], "./Results/area_weighted/Second_level/RSG2_SSI_Uganda.tif",
          type="Int32", mvFlag=-99999)

# write probabilities GeoTiffs
s <- s[sel.pr, ]
s@data <- round(data.frame(pred_soil$predictions)*100)
for(j in 1:ncol(s)){
  out <- paste0("./probs/", tax[j], ".tif")
  writeGDAL(s[j], out, type="Int32", mvFlag=-99999)
}
spplot((s["Haplic_Ferralsols"]), col.regions=SAGA_pal[[10]])

# CROSS VALIDATION (CV)
# CV by using "caret" package
# compile CV settings
trControl <- trainControl(method="LGOCV", number=10, p=.70, classProbs=T,
                      allowParallel=TRUE)
# parameter tuning
tuneGrid <- expand.grid(.mtry=fit$mtry,
                        .splitrule=fit$splitrule,
                        .min.node.size=fit$min.node.size)
set.seed(1234)
UTax.rf <- caret::train(fm,
                           data=reg.matrix[sel.na, ],
                           method="ranger",
                           tuneGrid=tuneGrid,
                           trControl=trControl,
                           importance="impurity",
                           num.trees=150)
UTax.rf

png("./Results/area_weighted/Second_level/covariates_importance.png", width=4096,
    height=2160, units='px', res=300)
plot(varImp(UTax.rf))
dev.off()
saveRDS(UTax.rf, "./Results/area_weighted/Second_level/caret_fit.rds")

# detailed CV for ranger
# https://github.com/Envirometrix/LandGISmaps/blob/master/functions/LandGIS_functions.R
# remove classes with too little observations
reg.matrix$ID <- as.character(1:nrow(reg.matrix))
keep.b <- levels(reg.matrix$First_level)[table(reg.matrix$First_level) > 10]
reg.matrix.t <- reg.matrix[reg.matrix$First_level %in% keep.b, ]
reg.matrix.t <- reg.matrix.t[complete.cases(reg.matrix.t[,all.vars(fm)]),]
# run 10-fold cross validation
FL.CV <- cv_factor(fm, rmatrix=reg.matrix.t, nfold=10, idcol="ID")
FL.CV[["Cohen.Kappa"]]
FL.CV[["Classes"]]
FL.CV[["Purity"]]
saveRDS(FL.CV, file="./Results/area_weighted/First_level/ranger_cv.rds")
write.csv(FL.CV[["Confusion.matrix"]], "./Results/area_weighted/First_level/cv_Confusion.matrix.csv")
write.csv(FL.CV[["Classes"]], "./Results/area_weighted/First_level/cv_classes.csv")

# END