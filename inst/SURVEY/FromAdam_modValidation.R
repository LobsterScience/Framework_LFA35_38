fit_cv = sdmTMB_cv(CodWt~
                     s(lZ,k=5),
                   data=aT,
                   time='WOS', 
                   mesh=bspde, 
                   family=tweedie(link='log'),
                   spatial='on',
                   fold_ids = 'fold_id',
                   
                   spatialtemporal='ar1',
                   k_folds=5,
                   #constant_mesh=F
)




mae<- function(x,y){
  sum(abs(x-y))/length(x)
}

rmse = function(x,y){
  sqrt((sum(y-x)^2)/length(x))
  
}


fit_cvTT = sdmTMBcv_tntpreds(fit_cv)

fitTT = dplyr::bind_rows(fit_cvTT)
fitTT$sqR = fitTT$CodWt - fitTT$pred
with(subset(fitTT,tt=='train'),mae(as.numeric(CodWt),as.numeric(pred)))
with(subset(fitTT,tt=='test'),mae(as.numeric(CodWt),as.numeric(pred)))
with(subset(fitTT,tt=='train'),rmse(as.numeric(CodWt),as.numeric(pred)))
with(subset(fitTT,tt=='test'),rmse(as.numeric(CodWt),as.numeric(pred)))

require(ggplot2)

ggplot(fitTT,aes(sqR,after_stat(density))) + 
  geom_histogram() + facet_wrap(~tt) + xlab('Residuals')