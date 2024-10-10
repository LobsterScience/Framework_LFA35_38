library(splines)
library(MASS) ##for glm.nb function

#remotes::install_github("TobieSurette/gulf.spatial", dependencies = TRUE)

library(sdmTMB)
library(sdmTMBextra)

##these are used to make the mesh
library(dplyr)
library(ggplot2)
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

fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)

survey = readRDS( file='fall_survey_data_all_commBio_Sept25.rds')
survey = st_as_sf(survey)

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs




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
survey = subset(survey,!is.na(Commwt))



ss = as_tibble(survey)
ss$depth.scaled = (ss$z-mean(ss$z))/sd(ss$z)

spde <- make_mesh(ss, xy_cols = c("X1000", "Y1000"),
                  n_knots=300,type = "cutoff_search")
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
data$year = data$Year
data$MNH <- ifelse(data$Survey=="MNH", 1, 0)
data$NEFSC <- ifelse(data$Survey=="NEFSC", 1, 0)

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
    formula = Commwt ~ 0 + MNH + NEFSC,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4 + bs5,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
  




##make maps 

model <- m ##to say which model to use for maps
model.lab <- "m6" ##to label maps


# Save model results:
saveRDS(model, file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "fall_Lobster model6 biomass.rds"))
model = readRDS( file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "fall_Lobster model6 biomass.rds")) #sept 25 300 knots and winner


data$resids <- residuals(model) # randomized quantile residuals
qqnorm(data$resids)
qqline(data$resids)



ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,z<max(data$z))
ff$depth.scaled = (ff$z-mean(data$z))/sd(data$z)
k <-  4
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

g = predict(model,newdata=f,return_tmb_object = T)

f34 = subset(f,LFA==35)
g34 = predict(model,newdata=f34,return_tmb_object = T)
ind35 = get_index(g34,bias_correct = T)
g5 = ggplot(ind35,aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Biomass (t) ')

f34 = subset(f,LFA==36)
g34 = predict(model,newdata=f34,return_tmb_object = T)
ind36 = get_index(g34,bias_correct = T)
g6 = ggplot(ind36,aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Biomass (t)')


f3538 = subset(f,LFA %in% 38)
g3538 = predict(model,newdata=f3538,return_tmb_object = T)
ind38 = get_index(g3538,bias_correct = T)
g8 = ggplot(ind38,aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Biomass (t)')


###everywhere
g = predict(model,newdata=f,return_tmb_object = T)
h = g$data
h$pred = model$family$linkinv(h$est)
gsf = st_as_sf(h)


##features

sf_use_s2(FALSE) #needed for cropping

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
ggplot(subset(da,year<2024)) +
  geom_sf(aes(fill=resids,color=resids),size=1.3) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
    facet_wrap(~year) +
  theme_test_adam()






###prediction surfaces


gsfr = subset(gsf,LFA %in% c(34,35,36,37,38) & Y1000>4850)

#mm = c(0.,max(quantile(gsfr$pred,.999)))
#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
gsfr$CommB = gsfr$pred
ggplot(subset(gsfr, year %in% c(1970:2023))) +
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



###se of predictions 
g = predict(model,newdata=f,se=T,nsim=50)
gg = model$family$linkinv(g)

f$se = apply(gg,1,sd)
f$med = apply(gg,1,median)
f$me = apply(gg,1,mean)

gsfer = st_as_sf(f)
gsfer = subset(gsfer,LFA %in% c(34,35,36,37,38)& Y1000 > 4850)

mm = c(0.,max(quantile(gsfer$se,.99)))
#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
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


#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
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
#depth splines
year=unique(data$year)
nd <- expand.grid(depth = seq(5, 500, by = 5), year = min(year):max(year))

nd$depth.scaled <- nd$depth - mean(model$data$z)/sd(model$data$z)
nd$survey.ID <- 0
k <- 4
knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])

basis <- as.data.frame(splines::bs(nd$depth.scaled, knots = knots, intercept = FALSE))
names(basis) <- paste0("bs", names(basis))
nd <- cbind(nd, basis)
nd$MNH = nd$NEFSC = 0
p <- predict(model, newdata = nd, se_fit = TRUE, re_form = NA)
p$pred = model$family$linkinv(p$est)
p$lower <- p$est - 1.96 * p$est_se
p$upper <- p$est + 1.96 * p$est_se

ggplot(p,aes(x=depth,y=est,group=year,colour=year))+geom_line()+scale_colour_viridis_c()+theme_test(base_size = 14)+labs(x='Depth (m)',y='Response')


