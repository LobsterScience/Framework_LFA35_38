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

survey = readRDS( file='summer_fall_survey_data_all_N_Sept24.rds')
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
survey = subset(survey,!is.na(N))



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


  
  # # Build B-spline basis functions 
  k <- 3
  knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
  bs <- bs(data$depth.scaled, knots = knots, intercept = FALSE)
  bs <- as.data.frame(bs)
  names(bs) <- paste0("bs", names(bs))
  data <- data[, setdiff(names(data), names(data)[grep("^bs[0-9]+", names(data))])]
  data_bs <- cbind(data, bs)
  
  
  m3a <- sdmTMB(
    data = data_bs,
    formula = N ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
m3=m3a


##############################################knots2
spde <- make_mesh(ss, xy_cols = c("X1000", "Y1000"),
                  n_knots=500,type = "cutoff_search")
plot(spde)

# Add on the barrier mesh component:
mesh <- sdmTMBextra::add_barrier_mesh(
  spde, ns_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)



# # Build B-spline basis functions 

m3b <- sdmTMB(
  data = data_bs,
  formula = N ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
  mesh = mesh,
  time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
  spatial = "on",
  family = tweedie(link = "log"),
  time = "year",
  spatiotemporal = "ar1"
)


##############################################knots3
spde <- make_mesh(ss, xy_cols = c("X1000", "Y1000"),
                  n_knots=800,type = "cutoff_search")
plot(spde)

# Add on the barrier mesh component:
mesh <- sdmTMBextra::add_barrier_mesh(
  spde, ns_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)



# # Build B-spline basis functions 

m3c <- sdmTMB(
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
saveRDS(list(m3a,m3b,m3c), file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster model3 abund_knots.rds"))
model = readRDS( file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster model3 abund_knots.rds")) #sept 25 300 knots and winner
m3a = model[[1]]
m3b = model[[2]]
m3c = model[[3]]


ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,z<max(data$z) & z>5)
ff$depth.scaled = (ff$z-mean(data$z))/sd(data$z)
k <-  3
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
g34a = predict(m3a,newdata=f34,return_tmb_object = T)
ind35a = get_index(g34a,bias_correct = T)


g34b = predict(m3b,newdata=f34,return_tmb_object = T)
ind35b = get_index(g34b,bias_correct = T)

g34c = predict(m3c,newdata=f34,return_tmb_object = T)
ind35c = get_index(g34c,bias_correct = T)



f34 = subset(f,LFA==36)
g36a = predict(m3a,newdata=f34,return_tmb_object = T)
ind36a = get_index(g36a,bias_correct = T)


g36b = predict(m3b,newdata=f34,return_tmb_object = T)
ind36b = get_index(g36b,bias_correct = T)

g36c = predict(m3c,newdata=f34,return_tmb_object = T)
ind36c = get_index(g36c,bias_correct = T)

f34 = subset(f,LFA==38)
g38a = predict(m3a,newdata=f34,return_tmb_object = T)
ind38a = get_index(g38a,bias_correct = T)


g38b = predict(m3b,newdata=f34,return_tmb_object = T)
ind38b = get_index(g38b,bias_correct = T)

g38c = predict(m3c,newdata=f34,return_tmb_object = T)
ind38c = get_index(g38c,bias_correct = T)

saveRDS(list(ind35a,ind35b,ind35c,ind36a,ind36b,ind36c,ind38a,ind38b,ind38c),file='testingKnotsAllsurveys.rds')


with(ind38c,plot(year,est))
with(ind38b,lines(year,est))
with(ind38a,lines(year,est))

with(ind36c,plot(year,est))
with(ind36b,lines(year,est))
with(ind36a,lines(year,est))

with(ind35c,plot(year,est))
with(ind35b,lines(year,est))
with(ind35a,lines(year,est))
