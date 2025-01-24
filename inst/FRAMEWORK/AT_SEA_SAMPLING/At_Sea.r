# at sea sampling 

require(bio.lobster)
require(bio.utilities)
require(dplyr)
require(tidyr)
require(ggplot2)
require(devtools)
require(sf)


ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,LFA %in% c(35,36,37,38) & z>5)
st_geometry(ff) <- NULL
ff = st_as_sf(ff,coords = c('X1','Y1'),crs=4326)



rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = rL[rL$LFA %in% c(35:38),]

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))


lobster.db('atSea')
v = subset(atSea,LFA %in% c(35,36,38))
v$UID = paste(v$TRIPNO,v$LFA,v$TRAPNO,v$STRINGNO,sep="_")
v$year = lubridate:::year(v$STARTDATE)
v$mm = lubridate::month(v$STARTDATE)
v$SYEAR = v$year
v$SYEAR = ifelse(v$mm %in% 10:12,v$SYEAR-1,v$SYEAR)
v$Season = ifelse(v$mm %in% c(10,11,12,1),'Fall',ifelse(v$mm %in% 2:4,'Winter',ifelse(v$mm %in% 5:7,'Spring','Outside')))

vv = subset(v,!is.na(LATITUDE))
vv$P=1
v1 = aggregate(P~LFA+LATITUDE+LONGITUDE+SYEAR+UID+mm,data=subset(vv),FUN=sum)
v1$Season = ifelse(v1$mm %in% c(10,11,12,1),'Fall',ifelse(v1$mm %in% 2:4,'Winter',ifelse(v1$mm %in% 5:7,'Spring','Outside')))

v1s = st_as_sf(v1,coords=c('LONGITUDE','LATITUDE'),crs=4326)


ss = st_nearest_feature(v1s,ff)
ds = st_distance(v1s,ff[ss,],by_element=T)
st_geometry(ba) = NULL
v1s$z = ff$z[ss]

d8= subset(ff,LFA==38)
vl = st_bbox(d8)
ggplot(d8)+geom_sf(size=2.6,aes(fill=z,colour=z))+  scale_fill_viridis_c() +
  scale_color_viridis_c()+
  geom_sf(data=ns_coast)+
  theme_test_adam() +
 coord_sf(xlim=c(vl[1],vl[3]),ylim=c(vl[2],vl[4]))
  

d6= subset(ff,LFA==36)
vl = st_bbox(d6)
ggplot(d6)+geom_sf(size=2.6,aes(fill=z,colour=z))+  scale_fill_viridis_c() +
  scale_color_viridis_c()+
  geom_sf(data=ns_coast)+
  theme_test_adam() +
  coord_sf(xlim=c(vl[1],vl[3]),ylim=c(vl[2],vl[4]))




ggplot(v1s)+geom_sf()+facet_wrap(~SYEAR)
