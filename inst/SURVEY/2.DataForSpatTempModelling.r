
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
fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)
crs_utm20 <- 32620


###data in

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
    st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
    ns_coast <- suppressWarnings(suppressMessages(
      st_crop(ns_coast,
              c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))
    ns_coast <- st_transform(ns_coast, crs_utm20)
    st_geometry(ns_coast) = st_geometry(ns_coast)
    st_crs(ns_coast) <- crs_utm20

ba = readRDS('~/git/bio.lobster.data/mapping_data/bathymetrySF.rds')
    ba = ba %>% st_as_sf() 
    st_geometry(ba) = st_geometry(ba)/1000
    st_crs(ba) = crs_utm20
    
    G <- suppressWarnings(suppressMessages(
      st_crop(ba,
              c(xmin = -82, ymin = 4539, xmax = 383, ymax = 5200))))
    
    
    rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
    rL = st_as_sf(rL)
    st_crs(rL) <- 4326
    rL = st_transform(rL,crs_utm20) 
    st_geometry(rL) <- st_geometry(st_as_sf(rL$geometry/1000)) 
    st_crs(rL) <- crs_utm20
    
    bathy= ba
    LFApolys = rL

x = RV_sets()
      cr = sum(x$Legal_wt,na.rm=T)/sum(x$WEIGHT_KG,na.rm=T)
      i = which(is.na(x$Legal))
      x$Legal_wt[i] = x$WEIGHT_KG[i] * cr
      x = subset(x,month(DATE) %in% 6:8)
      x = subset(x,select=c(mission,setno,LONGITUDE,LATITUDE,DATE,Legal_wt))
      x$Survey = 'RV'
      
y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = T)
y = subset(y,select=c(TRIP_ID,SET_NO, SET_LONG,SET_LAT,SET_DATE,SA_CORRECTED_PRORATED_N))    
y$Survey='ILTS'
names(y) = names(x)

survey = rbind(x,y)
survey = subset(survey, month(DATE) %in% c(6:8))
survey$YEAR = year(survey$DATE)
survey = st_as_sf(survey,coords = c('LONGITUDE','LATITUDE'),crs=4326)

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

survey = subset(survey,!is.na(Legal_wt))


#survey and depth from poly
        ss = st_nearest_feature(survey,ba)
        ds = st_distance(survey,ba[ss,],by_element=T)
        st_geometry(ba) = NULL
        survey$z = ba$z[ss]
        survey$z_dist = as.numeric(ds)

#windsorize
i = quantile(survey$Legal_wt,0.999)
survey$Legal_wt[which(survey$Legal_wt>i)]=i #windsorize extreme catches


saveRDS(list(survey,ns_coast,bathy,LFApolys),file='full_package_survey_coast_bathy_LFApolys.rds')

