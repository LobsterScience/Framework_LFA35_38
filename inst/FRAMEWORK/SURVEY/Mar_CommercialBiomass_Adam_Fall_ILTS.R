library(splines)
library(MASS) ##for glm.nb function

#remotes::install_github("TobieSurette/gulf.spatial", dependencies = TRUE)

library(sdmTMB)
library(sdmTMBextra)

##these are used to make the mesh
library(dplyr)
library(ggplot2)
library(assertthat)
library(INLA)

library(sf)


require(sdmTMB)
require(bio.lobster)
require(bio.utilities)
require(lubridate)
require(devtools)
require(dplyr)
require(ggplot2)
require(INLA)
options(stringAsFactors=F)
require(PBSmapping)
require(SpatialHub)
require(sf)
la()
crs_utm20 <- 32620

y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = T)
y = subset(y,select=c(TRIP_ID,SET_NO, SET_LONG,SET_LAT,SET_DATE,SA_CORRECTED_PRORATED_N))    
y$Survey='ILTS'

survey = (y)
survey$YEAR = year(survey$SET_DATE)
survey = subset(survey, month(SET_DATE)>8 & YEAR>2018)
survey = st_as_sf(survey,coords = c('SET_LONG','SET_LAT'),crs=4326)

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs

survey$Legal_wt = survey$SA_CORRECTED_PRORATED_N

survey = st_transform(survey,crs_utm20)

surv_utm_coords = st_coordinates(survey)
survey$X1000 <- surv_utm_coords[,1] /1000
survey$Y1000 <- surv_utm_coords[,2] /1000
st_geometry(survey) <- NULL
survey = st_as_sf(survey,coords = c('X1000','Y1000'),crs=crs_utm20)

sf_use_s2(FALSE) #needed for cropping
survey <- suppressWarnings(suppressMessages(
  st_crop(survey,
          c(xmin = -82, ymin = 4539, xmax = 383, ymax = 5200))))


# Project our survey data coordinates:


ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast <- st_transform(ns_coast, crs_utm20)
st_geometry(ns_coast) = st_geometry(ns_coast)
st_crs(ns_coast) <- crs_utm20

survey <- survey %>%   
  st_as_sf()

surv_utm_coords <- st_coordinates(survey)

survey$X1000 <- surv_utm_coords[,1] 
survey$Y1000 <- surv_utm_coords[,2] 
survey = subset(survey,!is.na(Legal_wt))


#add in bathy 
#allocating depth to location
ba = readRDS('~/git/bio.lobster.data/mapping_data/bathymetrySF.rds')
ba = ba %>% st_as_sf() 
st_geometry(ba) = st_geometry(ba)/1000
st_crs(ba) = 32620

ss = st_nearest_feature(survey,ba)
ds = st_distance(survey,ba[ss,],by_element=T)
st_geometry(ba) = NULL
survey$z = ba$z[ss]
survey$z_dist = as.numeric(ds)

ss = as_tibble(survey)
ss$depth.scaled = (ss$z-mean(ss$z))/sd(ss$z)

spde <- make_mesh(ss, xy_cols = c("X1000", "Y1000"),
                  n_knots=150,type = "cutoff_search")
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

models <-c( "m3", "m4",  "m6", "m7", "m8", "m9")


#models <-c("m3")

##number of k folds for model validation. 4 = 80% of data to train, 20$ to test
k_folds <- 5


##make the folds for model validation
data=ss
data$IDS = "I"

data = cv_SpaceTimeFolds(data,idCol = 'IDS',nfolds=k_folds)
data$year=data$YEAR
data$legal_wt = data$Legal_wt
##model selection table
##set up a blank table

columns <- c("Model", "Formula", "Time-varying", "AIC", "rho", "Matern range", "Spatial SD", "Sum loglik", "MAE_train","MAE_test", "RMSE_train", "RMSE_test" )
mod.select <- as.data.frame( matrix(data=NA, nrow =0, ncol=length(columns), byrow=TRUE))
colnames(mod.select) <- columns

##model selection table function

mod.select.fn <- function (){

  c<- as.data.frame( matrix(data=NA, nrow =1, ncol=length(columns), byrow=TRUE))
  colnames(c) <- columns
  c$Model <- mod.label
  c$Formula <-m$formula
  c$"Time-varying" <- ifelse(is.null (m$time_varying), NA, paste(m$time_varying[1], m$time_varying[2])  )
  
  ##spatial model 
  c$AIC <- AIC (m)
  c$rho <- m$sd_report[[1]]["rho"]
  c$`Matern range` <- m$sd_report[[1]]["range"]
  c$`Spatial SD` <- m$sd_report[[1]]["sigma_O"]

  
  ##model validation 
  m_cvTT = sdmTMBcv_tntpreds(m_cv)
  
  c$"Sum loglik" <- m_cv$sum_loglik
  fitTT = dplyr::bind_rows(m_cvTT)
  fitTT$sqR = fitTT$legal_wt - fitTT$pred
  

  c$MAE_test<-  with(fitTT[fitTT$tt=='test',],mae(as.numeric(legal_wt),as.numeric(pred)))
  c$MAE_train<-  with(fitTT[fitTT$tt=='train',],mae(as.numeric(legal_wt),as.numeric(pred)))
  c$RMSE_test <- with(fitTT[fitTT$tt=='test',],rmse(as.numeric(legal_wt),as.numeric(pred)))
  c$RMSE_train <- with(fitTT[fitTT$tt=='train',],rmse(as.numeric(legal_wt),as.numeric(pred)))
  
  
return(c)

}

# ###MODELS

if ("m1" %in% models) {
 
  mod.label <- "m1" 
m <- sdmTMB(
  data=data,
  formula = legal_wt ~ depth.scaled,
  mesh = mesh, # can be omitted for a non-spatial model
  spatial = "off",
  family = nbinom2(link = "log")
)

m_cv <- sdmTMB_cv(
  data=data,
  formula = legal_wt ~ depth.scaled,
  mesh = mesh, # can be omitted for a non-spatial model
  spatial = "off",
  family = nbinom2(link = "log"),
  fold_ids = 'fold_id',
  k_folds = k_folds
)

c <-mod.select.fn()
mod.select <- rbind(mod.select, c)

summary(m)
m1 <- m

}

##adding spatial effects
if ("m2" %in% models) {
  
  mod.label <- "m1" 
  
m <- sdmTMB(
  data=data,
  formula = legal_wt ~ depth.scaled,
  mesh = mesh,
  spatial = "on",
  family = tweedie(link = "log")
)

m_cv <- sdmTMB_cv(
  data=data,
  formula = legal_wt ~ depth.scaled,
  mesh = mesh,
  spatial = "on",
  family = tweedie(link = "log"),
  fold_ids = 'fold_id',
  k_folds = k_folds
)

c <-mod.select.fn()
mod.select <- rbind(mod.select, c)

summary(m)

m2 <- m
}

# ##adding spatial effects, time

if ("m3" %in% models) {

  mod.label <- "m3" 
  
  ##model
  
  m <- sdmTMB(
  data=data,
  formula = legal_wt  ~ depth.scaled,
  mesh = mesh,
  spatial = "on",
  family = tweedie(link = "log"),
  time = "year",
  spatiotemporal = "ar1")
  
  
  ##model cross-validation

  m_cv <- sdmTMB_cv(
    data=data,
    formula = legal_wt  ~ depth.scaled,
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
m3 <- m

}

# ##adding spatial effects,  time-varying depth relationship

if ("m4" %in% models) {

  mod.label <- "m4" 
  
m <- sdmTMB(  
  data=data,
  formula = legal_wt ~ 0, 
  mesh = mesh,
  time_varying = ~ 1 + depth.scaled, #the 1 means an annual intercept for depth
  spatial = "on",
  family = tweedie(link = "log"),
  time = "year",
  spatiotemporal = "ar1"
)

m_cv <- sdmTMB_cv(  
  data=data,
  formula = legal_wt ~ 0, 
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
m4 <- m

}

##tried adding survey as a predictor to see if improves model ##not running well
##AIC and BIC are higher. 

if ("m5" %in% models) {

  mod.label <- "m5" 
  
m <- sdmTMB(
  data=data,
  formula = legal_wt ~ as.factor(survey),
  mesh = mesh,
  time_varying = ~ 1 + depth.scaled, #the 1 means an annual intercept for depth
  spatial = "on",
  family = tweedie(link = "log"),
  time = "year",
  spatiotemporal = "ar1"
)

##the cross validation is throwing an error
m_cv <- sdmTMB_cv(
  data=data,
  formula = legal_wt ~ as.factor(survey),
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

m5 <- m
}

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
        formula = legal_wt ~ 0,
        mesh = mesh,
        time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4 + bs5,
        spatial = "on",
        family = tweedie(link = "log"),
        time = "year",
        spatiotemporal = "ar1"
        )

m_cv <- sdmTMB_cv(
  data = data_bs,
  formula = legal_wt ~ 0,
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

m6 <-m

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
    formula = legal_wt ~ 0,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = data_bs,
    formula = legal_wt ~ 0,
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
  
  m7 <-m
  
}



if ("m8" %in% models) {
  
  mod.label <- "m8" 
  
  # # Build B-spline basis functions 
  k <- 2
  knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
  bs <- bs(data$depth.scaled, knots = knots, intercept = FALSE)
  bs <- as.data.frame(bs)
  names(bs) <- paste0("bs", names(bs))
  data <- data[, setdiff(names(data), names(data)[grep("^bs[0-9]+", names(data))])]
  data_bs <- cbind(data, bs)
  
  m <- sdmTMB(
    data = data_bs,
    formula = legal_wt ~ 0,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = data_bs,
    formula = legal_wt ~ 0,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3,
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
  
  m8 <-m
  
}


if ("m9" %in% models) {
  
  mod.label <- "m9" 
  
  # # Build B-spline basis functions 
  k <- 1
  knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
  bs <- bs(data$depth.scaled, knots = knots, intercept = FALSE)
  bs <- as.data.frame(bs)
  names(bs) <- paste0("bs", names(bs))
  data <- data[, setdiff(names(data), names(data)[grep("^bs[0-9]+", names(data))])]
  data_bs <- cbind(data, bs)
  
  m <- sdmTMB(
    data = data_bs,
    formula = legal_wt ~ 0,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = data_bs,
    formula = legal_wt ~ 0,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2,
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
  
  m9 <-m
  
}


##extra step because formulas mess up saving to a csv
mod.select = data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)

write.csv(mod.select, file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS','ilts_leg2_only.csv'), row.names = FALSE)

####predictions 



##make maps 

model <- m8 ##to say which model to use for maps


# Save model results:
saveRDS(model, file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "ILTS_Lobster model8 leg2 biomass.rds"))
model = readRDS( file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "ILTS_Lobster model8 leg2 biomass.rds"))

data=model$data
##make a QQplot

data$resids <- residuals(model) # randomized quantile residuals
qqnorm(data$resids)
qqline(data$resids)


head(data)


ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,z<max(ss$z))
ff$depth.scaled = (ff$z-mean(ff$z))/sd(ff$z)
k <-  2
knots <-   quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
basis <- as.data.frame(bs(ff$depth.scaled, knots = knots, intercept = FALSE))
names(basis) <- paste0("bs", names(basis))
ff <- cbind(ff, basis)

yy = unique(ss$YEAR)
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


f34 = subset(f,LFA==35)
g34 = predict(model,newdata=f34,return_tmb_object = T)
ind34 = get_index(g34,bias_correct = T)
ggplot(ind34,aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test()+labs(x='Year',y='Commercial Biomass (t)')


f3538 = subset(f,LFA %in% 36)
g3538 = predict(model,newdata=f3538,return_tmb_object = T)
ind3538 = get_index(g3538,bias_correct = T)
ggplot(ind3538,aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test()+labs(x='Year',y='Commercial Biomass (t)')

###everywhere
g = predict(model,newdata=f,return_tmb_object = T)
h = g$data
h$pred = model$family$linkinv(h$est)
gsf = st_as_sf(h)


##features

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast = st_transform(ns_coast,32620) 
st_geometry(ns_coast) <- st_geometry(ns_coast)/1000
st_crs(ns_coast) <- 32620




rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = rL[rL$LFA %in% c(33:38),]
st_crs(rL) <- 4326
crs_utm20 <- 32620
#rL = rL[-which(!(st_is_valid(rL))),]
rL <- suppressWarnings(suppressMessages(
  st_crop(rL,
          c(xmin = -67.5, ymin = 42, xmax = -62.1, ymax = 46))))
rL <- st_transform(rL, crs_utm20)
st_geometry(rL) <- st_geometry(rL)/1000
st_crs(rL) <- 32620


####maps of residuals
da = st_as_sf(data)
ggplot(subset(da,year>2012)) +
  geom_sf(aes(fill=resids,color=resids),size=1.3) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
    facet_wrap(~year) +
  theme_test_adam()






###prediction surfaces


gsfr = subset(gsf,LFA %in% c(35,36,37))

mm = c(0.,max(quantile(gsfr$pred,.999)))
gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
ggplot(subset(gsfr, year %in% c(2013:2023))) +
  geom_sf(aes(fill=CommB,color=CommB),size=2.1) + 
  scale_fill_viridis_c(trans='sqrt',limits=mm) +
  scale_color_viridis_c(trans='sqrt',limits=mm) +
  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(st_bbox(gsfr)$ymin,st_bbox(gsfr)$ymax),
           expand = FALSE)



