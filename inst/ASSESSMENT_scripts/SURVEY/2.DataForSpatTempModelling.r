require(tidyr)
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
fd=file.path(project.datadirectory('Assessment_LFA35_38'),'outputs','SURVEYS')
setwd(fd)
crs_utm20 <- 32620
sf_use_s2(FALSE) #needed for cropping

###data in

ns_coast =readRDS(file.path( bio.directory, "bio.lobster.data","mapping_data","CoastSF.rds"))
    st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
    ns_coast <- suppressWarnings(suppressMessages(
      st_crop(ns_coast,
              c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))
    ns_coast <- st_transform(ns_coast, crs_utm20)
    st_geometry(ns_coast) = st_geometry(ns_coast)
    st_crs(ns_coast) <- crs_utm20

ba = readRDS(file.path(bio.directory, "bio.lobster.data","mapping_data",'bathymetrySF.rds'))
    ba = ba %>% st_as_sf() 
    st_geometry(ba) = st_geometry(ba)/1000
    st_crs(ba) = crs_utm20
    
    rL = readRDS(file.path( bio.directory, "bio.lobster.data","mapping_data","LFAPolys37Split.rds"))
    rL = st_as_sf(rL)
    st_crs(rL) <- 4326
    rL = st_transform(rL,crs_utm20) 
    st_geometry(rL) <- st_geometry(st_as_sf(rL$geometry/1000)) 
    st_crs(rL) <- crs_utm20
    
    bathy= ba
    LFApolys = rL


#################################################################################################################
####summer recruit N and commercial N along with fall commercial N for the full model
#summer
y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(70,300),aggregate=T,species=2550,biomass = F)
y = subset(y,month(SET_DATE)<8)
y$Survey='ILTS_Summer'
y$id = paste(y$TRIP_ID,y$SET_NO,sep="-")
y = subset(y,select=c(id,YEAR,SET_LONG,SET_LAT,SA_CORRECTED_PRORATED_N,Survey))

x = RV_sets()
x = subset(x,month(DATE) %in% 6:8)
x$N = x$Recruit+x$Legal
x$N = x$N/x$OFFSET

x$id = paste(x$mission,x$setno,sep="-")
x = subset(x,select=c(id,YEAR,LONGITUDE,LATITUDE,N))
x$LONGITUDE = ifelse(x$LONGITUDE>0,x$LONGITUDE*-1,x$LONGITUDE )
x$Survey = 'RV_Summer'

names(y) = names(x)
suC = st_as_sf(rbind(y,x),coords=c('LONGITUDE','LATITUDE'),crs=4326)

#fall

      z = lobster.db('DMR.redo')
      z$N = z$Comm
      z$YEAR = z$Year
      z = subset(z,select=c(id,YEAR,N))
      z$Survey = 'MNH_Fall'
      
      ###NEFSC
      
      n = NEFSC_sets()
      n = subset(n,month(DATE)>8 & LONGITUDE> -68 & LATITUDE>42)
      n$N = (n$Legal)/n$OFFSET
      n$id = paste(n$MISSION,n$SETNO,sep="-")
      n = subset(n,select=c(id,YEAR,LONGITUDE,LATITUDE,N))
      n$Survey = 'NEFSC_Fall'
      n = st_as_sf(n,coords=c('LONGITUDE','LATITUDE'),crs=4326)

      #ILTS
      y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = F)
      y = subset(y,month(SET_DATE)>8 & LFA %ni% 'L34')    
      y$id = paste(y$TRIP_ID,y$SET_NO,sep="-")
      y$Survey = 'ILTS_Fall'
      y = subset(y,select=c(id,YEAR,SET_LONG,SET_LAT,SA_CORRECTED_PRORATED_N,Survey))
      
      ys = st_as_sf(y,coords=c('SET_LONG','SET_LAT'),crs=4326)
      ys$N = ys$SA_CORRECTED_PRORATED_N
      ys = subset(ys,select=c(id,YEAR,N,Survey))

      suCF =do.call(rbind,list(z,n,ys))
      
combSUFall = do.call('rbind',list(suC,suCF))

ggplot(combSUFall,aes(colour=Survey))+geom_sf()

combSUFall = st_transform(combSUFall,crs_utm20)

surv_utm_coords = st_coordinates(combSUFall)
combSUFall$X1000 <- surv_utm_coords[,1] /1000
combSUFall$Y1000 <- surv_utm_coords[,2] /1000
st_geometry(combSUFall) <- NULL
survey = st_as_sf(combSUFall,coords = c('X1000','Y1000'),crs=crs_utm20)

sf_use_s2(FALSE) #needed for cropping
#survey and depth from poly
ss = st_nearest_feature(survey,ba)
ds = st_distance(survey,ba[ss,],by_element=T)
st_geometry(ba) = NULL
survey$z = ba$z[ss]
survey$z_dist = as.numeric(ds)

saveRDS(survey,file='survey_data_all_combined_jan152026.rds')

b = aggregate(N~YEAR+Survey,data=survey,FUN=mean)
ggplot(b,aes(x=YEAR,y=N))+geom_point()+facet_wrap(~Survey,scales='free_y')
