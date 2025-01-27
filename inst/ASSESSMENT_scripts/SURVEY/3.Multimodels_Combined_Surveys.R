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

fd=file.path(project.datadirectory('Assessment_LFA35_38'),'outputs','SURVEYS')
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
label='jan'
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
write.csv(mod.select, file = file.path(project.datadirectory('Assessment_LFA35_38'),'outputs','SURVEYS','allsurveys_model runs_jan17_2025.csv'), row.names = FALSE)

# Save model results:
###model m6 best by criteria but residuals had a lot of spatial and temporal patterns and thus discarded.
model = m7

saveRDS(model, file = file.path(project.datadirectory('Assessment_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster model7 abund_jan17.rds"))
model = readRDS( file = file.path(project.datadirectory('Assessment_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster model7 abund_jan17.rds"))


data$resids <- residuals(model) # randomized quantile residuals
qqnorm(data$resids)
qqline(data$resids)

####maps of residuals
da = st_as_sf(data)
ggplot(subset(da)) +
  geom_sf(aes(fill=resids,color=resids),size=1.3) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  facet_wrap(~year) +
  theme_test_adam()
