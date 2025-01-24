library(splines)
library(MASS)
library(sdmTMB)
library(sdmTMBextra)
library(dplyr)
library(ggplot2)
require(bio.lobster)
require(bio.utilities)
require(lubridate)
require(devtools)
require(dplyr)
require(ggplot2)
require(INLA)
options(stringAsFactors=F)
require(sf)
la()
crs_utm20 <- 32620

fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)
path = fd
survey = readRDS( file='survey_data_all_combined.rds')
survey = st_as_sf(survey)

ns_coast =readRDS(file.path( bio.directory,"bio.lobster.data", "mapping_data","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs


sf_use_s2(FALSE) #needed for cropping
survey <- suppressWarnings(suppressMessages(
  st_crop(survey,
          c(xmin = -82, ymin = 4539, xmax = 383, ymax = 5200))))


# Project our survey data coordinates:


ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -70, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast <- st_transform(ns_coast, crs_utm20)
st_geometry(ns_coast) = st_geometry(ns_coast)
st_crs(ns_coast) <- crs_utm20

survey <- survey %>%   
  st_as_sf()

surv_utm_coords <- st_coordinates(survey)

survey$X1000 <- surv_utm_coords[,1] 
survey$Y1000 <- surv_utm_coords[,2] 
survey = subset(survey,!is.na(N))



ss = as_tibble(survey)
ss$depth.scaled = (ss$z-mean(ss$z))/sd(ss$z)
ss$depth.scaled2 = ss$depth.scaled * ss$depth.scaled
data = ss
data$year = data$YEAR
data = subset(data,year>=2000)

spde <- make_mesh(data, xy_cols = c("X1000", "Y1000"),
                  n_knots=350,type = "cutoff_search")
plot(spde)

# Add on the barrier mesh component:
mesh <- sdmTMBextra::add_barrier_mesh(
  spde, ns_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)
mesh_df_water <- mesh$mesh_sf[mesh$normal_triangles, ]
mesh_df_land <- mesh$mesh_sf[mesh$barrier_triangles, ]
ggplot(ns_coast) +
  geom_sf() +
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  geom_sf(data = mesh_df_land, size = 1, colour = "green")+
  coord_sf()



## make new columns to model survey effecst (0 in both columns is ILTS_summer)
data$ILTS_Fall <- ifelse(data$Survey=="ILTS_Fall", 1, 0)
data$MNH_Fall <- ifelse(data$Survey=="MNH_Fall", 1, 0)
data$MNH_NEFSC_Fall <- ifelse(data$Survey=="NEFSC_Fall", 1, 0)
data$RV_Summer <- ifelse(data$Survey=="RV_Summer", 1, 0)
k_folds <-5


##make the folds for model validation

data$IDS = "I"
data = cv_SpaceTimeFolds(data,idCol = 'IDS',nfolds=k_folds)
data$n = data$N

##model selection table
##set up a blank table

columns <- c("Model", "Formula", "Time-varying", "Family", "AIC", "rho", "Matern range", "Spatial SD", "Hessian_positive", "Sum loglik", "MAE_train","MAE_test", "RMSE_train", "RMSE_test" )
mod.select <- as.data.frame( matrix(data=NA, nrow =0, ncol=length(columns), byrow=TRUE))
colnames(mod.select) <- columns

##model selection table function

mod.select.fn <- function (){
  
  c<- as.data.frame( matrix(data=NA, nrow =1, ncol=length(columns), byrow=TRUE))
  colnames(c) <- columns
  c$Model <- mod.label
  c$Formula <-m$formula [1]
  c$"Time-varying" <- ifelse(is.null (m$time_varying), NA, paste(m$time_varying[1], m$time_varying[2])  )
  c$"Family" <- ifelse(m$family[1]=="tweedie", paste0(m$family[1], "(link = ", m$family[2], ")"), m$family["clean_name"])
  
  
  ##spatial model 
  c$AIC <- AIC (m)
  c$rho <- m$sd_report[[1]]["rho"]
  c$`Matern range` <- m$sd_report[[1]]["range"]
  c$`Spatial SD` <- m$sd_report[[1]]["sigma_O"]
  c$"Hessian_positive" <- m$pos_def_hessian
  
  
  c$"Sum loglik" <- m_cv$sum_loglik
  
  m_cvTT = sdmTMBcv_tntpreds(m_cv)
  fitTT = dplyr::bind_rows(m_cvTT)
  fitTT$sqR = fitTT$n - fitTT$pred
  c$MAE_test<-  with(fitTT[fitTT$tt=='test',],mae(as.numeric(n),as.numeric(pred)))
  c$MAE_train<-  with(fitTT[fitTT$tt=='train',],mae(as.numeric(n),as.numeric(pred)))
  c$RMSE_test <- with(fitTT[fitTT$tt=='test',],rmse(as.numeric(n),as.numeric(pred)))
  c$RMSE_train <- with(fitTT[fitTT$tt=='train',],rmse(as.numeric(n),as.numeric(pred)))
  
  
  return(c)
  
}

# ###what models to run
models <-c( "m3", "m4",'m6','m7', "m9", "m10") 
# ##adding spatial effects, time

if ("m3" %in% models) {
  
  mod.label <- "m3" 
  
  ##model
  
  m <- sdmTMB(
    data=data,
    formula = n  ~ depth.scaled + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1")
  
  
  ##model cross-validation
  
  m_cv <- sdmTMB_cv(
    data=data,
    formula = n  ~ depth.scaled + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh, 
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(path, "MAR_Lobster_ALL", "_", mod.label, "_", label,  ".rdata"))
  
  m3 <- m
  m3_cv <- m_cv
}



# ##adding spatial effects,  time-varying depth relationship

if ("m4" %in% models) {
  
  mod.label <- "m4" 
  
  m <- sdmTMB(  
    data=data,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall, 
    mesh = mesh,
    time_varying = ~ 1 + depth.scaled, #the 1 means an annual intercept for depth
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall, 
    mesh = mesh,
    time_varying = ~ 1 + depth.scaled, #the 1 means an annual intercept for depth
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(path, "MAR_Lobster_ALL", "_", mod.label, "_", label,  ".rdata"))
  
  m4 <- m
  m4_cv <- m_cv
  
}

##adding basis splines 
if ("m6" %in% models) {
  
  mod.label <- "m6" 
  
  # # Build B-spline basis functions 
  k <- 4
  knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
  bs <- bs(data$depth.scaled, knots = knots, intercept = FALSE)
  bs <- as.data.frame(bs)
  names(bs) <- paste0("bs", names(bs))
  data <- data[, setdiff(names(data), names(data)[grep("^bs[0-9]+", names(data))])]
  data_bs <- cbind(data, bs)
  
  m <- sdmTMB(
    data = data_bs,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4 + bs5,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = data_bs,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4 + bs5,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    fold_ids = 'fold_id',
    k_folds=k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  # Save model results:
  save(m, file = paste0(path, "MAR_Lobster_ALL", "_", mod.label, "_", label,  ".rdata"))
  m6 <-m
  m6_cv <-m_cv
  
}

if ("m7" %in% models) {
  
  mod.label <- "m7" 
  
  # # Build B-spline basis functions 
  k <- 3
  knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
  bs <- bs(data$depth.scaled, knots = knots, intercept = FALSE)
  bs <- as.data.frame(bs)
  names(bs) <- paste0("bs", names(bs))
  data <- data[, setdiff(names(data), names(data)[grep("^bs[0-9]+", names(data))])]
  data_bs <- cbind(data, bs)
  
  m <- sdmTMB(
    data = data_bs,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = data_bs,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall + RV_Summer,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    fold_ids = 'fold_id',
    k_folds=k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  # Save model results:
  save(m, file = paste0(path, "MAR_Lobster_ALL", "_", mod.label, "_", label,  ".rdata"))
  
  m7_cv <- m_cv
  m7 <-m
  
}

if ("m9" %in% models) {
  
  mod.label <- "m9" 
  
  
  
  m <- sdmTMB(  
    data=data,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall, 
    mesh = mesh,
    time_varying = ~ 1 + depth.scaled + depth.scaled2, #the 1 means an annual intercept for depth
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = n ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall, 
    mesh = mesh,
    time_varying = ~ 1 + depth.scaled + depth.scaled2, #the 1 means an annual intercept for depth
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  # Save model results:
  save(m, file = paste0(path, "MAR_Lobster_ALL", "_", mod.label, "_", label,  ".rdata"))
  m9 <-m
  m9_cv <- m_cv
  
}

if ("m10" %in% models) {
  
  mod.label <- "m10" 
  
  
  
  m <- sdmTMB(  
    data=data,
    formula = n ~ 0 + s(depth.scaled, year) + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall, 
    mesh = mesh,
    time_varying = ~ 1, #the 1 means an annual intercept for depth
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = n ~ 0 + s(depth.scaled, year) + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall + RV_Summer, 
    mesh = mesh,
    time_varying = ~ 1, #the 1 means an annual intercept for depth
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  # Save model results:
  save(m, file = paste0(path, "MAR_Lobster_ALL", "_", mod.label, "_", label,  ".rdata"))
  m10 <-m
  m10_cv <- m_cv
  
}



##extra step because formulas mess up saving to a csv
mod.select = data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)
write.csv(mod.select, file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS','allsurveys_model runs_sept28_final.csv'), row.names = FALSE)

# Save model results:
###model m6 best by criteria but residuals had a lot of spatial and temporal patterns and thus discarded.
model = m6

saveRDS(model, file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster model7 abund_setp29.rds"))
model = readRDS( file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster model7 abund_setp29.rds"))


data$resids <- residuals(model) # randomized quantile residuals
qqnorm(data$resids)
qqline(data$resids)

####maps of residuals
da = st_as_sf(data)
ggplot(subset(da,year<2024)) +
  geom_sf(aes(fill=resids,color=resids),size=1.3) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  facet_wrap(~year) +
  theme_test_adam()


#prediction surface
ff = readRDS(file.path( bio.directory,"bio.lobster.data", "mapping_data","bathy_by_lfa_37split.rds"))
ff = subset(ff,z<max(data$z) & z>5 & PID %in% c(34,35,36,38))
ff$LFA = ff$PID
ff$depth.scaled = (ff$z-mean(data$z))/sd(data$z)

k <-  3 #need to specify based on based model
knots <-   quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
basis <- as.data.frame(bs(ff$depth.scaled, knots = knots, intercept = FALSE))
names(basis) <- paste0("bs", names(basis))
ff <- cbind(ff, basis)

yy = unique(data$year)
o = list()
for(i in 1:length(yy)){
  f = ff
  f$year=yy[i]
  o[[i]]=f
}
ff = do.call(rbind,o)

sg = st_coordinates(ff)
ff$X1000 = sg[,1]
ff$Y1000 = sg[,2]
ff$SOURCE='ILTS'
f = as_tibble(ff)
f$MNH <- 0
f$NEFSC <- 0
f$ILTS_Fall <- f$MNH_Fall <- f$MNH_NEFSC_Fall <- 0

#g = predict(model,newdata=f,return_tmb_object = T)

f34 = subset(f,LFA==35)
g34 = predict(model,newdata=f34,return_tmb_object = T)
ind35 = get_index(g34,bias_correct = T)
g5 = ggplot(subset(ind35,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000) ')

f34 = subset(f,LFA==36)
g34 = predict(model,newdata=f34,return_tmb_object = T)
ind36 = get_index(g34,bias_correct = T)
g6 = ggplot(subset(ind36,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


f3538 = subset(f,LFA %in% 38)
g3538 = predict(model,newdata=f3538,return_tmb_object = T)
ind38 = get_index(g3538,bias_correct = T)
g8 = ggplot(subset(ind38,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


saveRDS(list(ind35,ind36,ind38),'IndicesFromFullComboModelSept29.rds')


###everywhere
g = predict(model,newdata=f,return_tmb_object = T)
h = g$data
h$pred = model$family$linkinv(h$est)
gsf = st_as_sf(h)


##features

sf_use_s2(FALSE) #needed for cropping

ns_coast =readRDS(file.path( bio.directory,'bio.lobster.data',"mapping_data","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast = st_transform(ns_coast,32620) 
st_geometry(ns_coast) <- st_geometry(ns_coast)/1000
st_crs(ns_coast) <- 32620



rL = readRDS(file.path( bio.directory,'bio.lobster.data',"mapping_data","LFAPolys37Split.rds"))
rL = rL[rL$PID %in% c(33:38),]
st_crs(rL) <- 4326
crs_utm20 <- 32620
#rL = rL[-which(!(st_is_valid(rL))),]
rL <- suppressWarnings(suppressMessages(
  st_crop(rL,
          c(xmin = -67.5, ymin = 42, xmax = -62.1, ymax = 46))))
rL <- st_transform(rL, crs_utm20)
st_geometry(rL) <- st_geometry(rL)/1000
st_crs(rL) <- 32620



###prediction surfaces


gsfr = subset(gsf,LFA %in% c(34,35,36,37,38) & Y1000>4850)

#mm = c(0.,max(quantile(gsfr$pred,.999)))
#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
gsfr$CommB = gsfr$pred
ggplot(subset(gsfr, year %in% c(2000:2023))) +
  geom_sf(aes(fill=CommB,color=CommB),size=2.1) + 
  scale_fill_viridis_c(trans='log') +
  scale_color_viridis_c(trans='log') +
  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(st_bbox(gsfr)$ymin,st_bbox(gsfr)$ymax),
           expand = FALSE)

#####################
##### simulations for defining LRPEIV

simsSurf = function(model, newD=f,nsims=100){
      g = predict(model,newdata=newD,nsim = nsims)
      gg = model$family$linkinv(g)
      gg1 = as.data.frame(gg)
      gg1$year = rownames(gg)
      ga = aggregate(.~year,data=gg1,FUN=sum) #each column is a biomass estimate for the LFA
          return(ga)    
      }

f_final = subset(f, LFA ==35)
b5 = simsSurf(model=model,newD=f_final)

f_final = subset(f, LFA ==36)
b6 = simsSurf(model=model,newD=f_final)

f_final = subset(f, LFA ==38)
b8 = simsSurf(model=model,newD=f_final)
