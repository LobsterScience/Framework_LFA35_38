library(splines)
library(MASS) ##for glm.nb function
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
survey = readRDS( file='survey_data_all_combined.rds')
survey = st_as_sf(survey)

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs

sf_use_s2(FALSE) #needed for cropping
survey <- suppressWarnings(suppressMessages(
  st_crop(survey,
          c(xmin = -82, ymin = 4539, xmax = 383, ymax = 5200))))


ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -70, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast <- st_transform(ns_coast, crs_utm20)
st_geometry(ns_coast) = st_geometry(ns_coast)
st_crs(ns_coast) <- crs_utm20

surv_utm_coords <- st_coordinates(survey)

survey$X1000 <- surv_utm_coords[,1] 
survey$Y1000 <- surv_utm_coords[,2] 
survey = subset(survey,!is.na(N))

survey = subset(survey,YEAR>1999)

ss = as_tibble(survey)
ss$depth.scaled = (ss$z-mean(ss$z))/sd(ss$z)
ss$depth.scaled2 = ss$depth.scaled * ss$depth.scaled
spde <- make_mesh(ss, xy_cols = c("X1000", "Y1000"),
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

data = ss
data$year = data$YEAR
## make new columns to model survey effecst (0 in both columns is ILTS_summer)
data$ILTS_Fall <- ifelse(data$Survey=="ILTS_Fall", 1, 0)
data$MNH_Fall <- ifelse(data$Survey=="MNH_Fall", 1, 0)
data$MNH_NEFSC_Fall <- ifelse(data$Survey=="NEFSC_Fall", 1, 0)
data$RV_Summer <- ifelse(data$Survey=="RV_Summer", 1, 0)


  ##########this section will change depending on teh best model from the Multimodel section in 3.
  # # Build B-spline basis functions ; needed to have time varying depth relationship
  k <- 3 
  knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
  bs <- bs(data$depth.scaled, knots = knots, intercept = FALSE)
  bs <- as.data.frame(bs)
  names(bs) <- paste0("bs", names(bs))
  data <- data[, setdiff(names(data), names(data)[grep("^bs[0-9]+", names(data))])]
  data_bs <- cbind(data, bs)
  
  
  model <- sdmTMB(
    data = data_bs,
    formula = N ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  


# Save model results:
saveRDS(model, file = file.path(project.datadirectory('Assessment_LFA35_38'),'outputs','SURVEYS', "bestModel_commLobster2000+feb142025.rds"))


data$resids <- residuals(model) # randomized quantile residuals
qqnorm(data$resids)
qqline(data$resids)

ff = readRDS(file.path(bio.directory,"bio.lobster.data","mapping_data",'bathy_byLFA_noLFA37_noMidas_noPass.rds')) #fixed to remove LFA 37 from all predictions Feb 11, 2025
ff = subset(ff,PID %in% 34:38)
ff = subset(ff,z<max(data$z) & z>5)
ff$depth.scaled = (ff$z-mean(data$z))/sd(data$z)

##k from model above
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


#generate indices for each area

f35 = subset(f,PID==35)
g35 = predict(model,newdata=f35,return_tmb_object = T)
ind35 = get_index(g35,bias_correct = T)
g5 = ggplot(subset(ind35),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000) ')

f36 = subset(f,PID==36)
g36 = predict(model,newdata=f36,return_tmb_object = T)
ind36 = get_index(g36,bias_correct = T) 
g6 = ggplot(subset(ind36),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


f38 = subset(f,PID %in% 38)
g38 = predict(model,newdata=f38,return_tmb_object = T)
ind38 = get_index(g38,bias_correct = T)
g8 = ggplot(subset(ind38),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


saveRDS(list(ind35,ind36,ind38),'IndicesFromFullComboModelFeb14_2000+.rds')
#b = readRDS('IndicesFromFullComboModelJan17_2000+.rds')
#ind351=b[[1]]
#ind361=b[[2]]
#ind381=b[[3]]

g5 = ggplot(subset(ind35,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000) ')
g6 = ggplot(subset(ind36,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')
g8 = ggplot(subset(ind38),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


###everywhere
      g = predict(model,newdata=f,return_tmb_object = T)
      h = g$data
      h$pred = model$family$linkinv(h$est)
      gsf = st_as_sf(h)


##features

sf_use_s2(FALSE) #needed for cropping

ns_coast =readRDS(file.path(bio.directory,'bio.lobster.data',"mapping_data","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast = st_transform(ns_coast,32620) 
st_geometry(ns_coast) <- st_geometry(ns_coast)/1000
st_crs(ns_coast) <- 32620



rL = readRDS(file.path(bio.directory,'bio.lobster.data',"mapping_data","LFAPolys37Split.rds"))
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


####maps of residuals
da = st_as_sf(data)
ggplot(subset(da)) +
  geom_sf(aes(fill=resids,color=resids),size=1.3) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
    facet_wrap(~year) +
  theme_test_adam()




gsfr = subset(gsf,PID %in% c(34,35,36,38) & Y1000>4850)

#mm = c(0.,max(quantile(gsfr$pred,.999)))
#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
gsfr$CommN = gsfr$pred
gb = ggplot(subset(gsfr, year %in% c(2023))) +
  geom_sf(aes(fill=CommN,color=CommN),size=2.1) + 
  scale_fill_viridis_c(trans='log') +
  scale_color_viridis_c(trans='log') +
#  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(st_bbox(gsfr)$ymin,st_bbox(gsfr)$ymax),
           expand = FALSE)


###se of predictions 
g = predict(model,newdata=f,se=T,nsim=500) #dec 18 2024 need to finish
gg = model$family$linkinv(g)

f$se = apply(gg,1,sd)
f$med = apply(gg,1,median)
f$me = apply(gg,1,mean)

gsfer = st_as_sf(f)
gsfer = subset(gsfer,PID %in% c(34,35,36,38)& Y1000 > 4850)

mm = c(0.,max(quantile(gsfer$se,.99)))
gsfer$PredictionError = gsfer$se
ggplot(subset(gsfer, year %in% c(1970:2023))) +
  geom_sf(aes(fill=PredictionError,color=PredictionError),size=2.1) + 
  scale_fill_viridis_c(trans='sqrt',limits=mm) +
  scale_color_viridis_c(trans='sqrt',limits=mm) +
  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(4850,st_bbox(gsfr)$ymax),
           expand = FALSE)


gsfer$CV_Prediction = gsfer$se/gsfer$me
mm = c(0.,max(quantile(gsfer$CV_Prediction,.995)))

ggplot(subset(gsfer, year %in% c(1970:2023))) +
  geom_sf(aes(fill=CV_Prediction,color=CV_Prediction),size=2.1) + 
  scale_fill_viridis_c(limits=mm) +
  scale_color_viridis_c(limits=mm) +
  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(st_bbox(gsfr)$ymin,st_bbox(gsfr)$ymax),
           expand = FALSE)



##########
#depth splines; not needed, but left if you want it
year=unique(data$year)
nd <- expand.grid(depth = seq(5, 500, by = 5), year = min(year):max(year))

nd$depth.scaled <- nd$depth - mean(model$data$z)/sd(model$data$z)
nd$survey.ID <- 0
k <- 
knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])

basis <- as.data.frame(splines::bs(nd$depth.scaled, knots = knots, intercept = FALSE))
names(basis) <- paste0("bs", names(basis))
nd <- cbind(nd, basis)
nd$ILTS_Fall <- nd$MNH_Fall <- nd$MNH_NEFSC_Fall <- nd$RV_Summer <- 0

p <- predict(model, newdata = nd, se_fit = TRUE, re_form = NA)
p$pred = model$family$linkinv(p$est)
p$lower <- p$est - 1.96 * p$est_se
p$upper <- p$est + 1.96 * p$est_se

ggplot(p,aes(x=depth,y=est,group=year,colour=year))+geom_line()+scale_colour_viridis_c()+theme_test(base_size = 14)+labs(x='Depth (m)',y='Response')

