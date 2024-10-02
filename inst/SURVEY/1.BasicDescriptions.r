##ILTS and RV survey overviews

require(bio.lobster)
require(devtools)
require(bio.utilities)
require(sf)
require(ggplot2)
require(tidyr)
theme_set(theme_test(base_size = 14))
fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)

p = ggLobsterMap('34-38',addGrids = F)


##RV Survey
x = RV_sets()
x = subset(x,month(x$DATE) %in% c(6,7,8))
xs = st_as_sf(x,coords = c('LONGITUDE','LATITUDE'),crs=4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326

xs = st_join(xs,rL,join=st_within)
xs = subset(xs,LFA %in% c(34,35,36,37,38))


xle = xs %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
#xle= na.zero(xle,cols='value')
xxa = aggregate(value~Length+YEAR,data=xle,FUN=mean)
ggplot(subset(xxa,YEAR>1998),aes(x=Length,y=value))+geom_bar(stat='identity')+facet_wrap(~YEAR,scales = 'free_y')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')

br = round(seq(min(xs$YEAR),max(xs$YEAR),length=5))
lab = c(paste(br[1],br[2],sep="-"),paste(br[2]+1,br[3],sep="-"),paste(br[3]+1,br[4],sep="-"),paste(br[4]+1,br[5],sep="-"))
xs$rYR=cut(xs$YEAR,breaks=br,labels=lab,include.lowest = T,include.highest=F)
#to check
aggregate(YEAR~rYR,data=xs,FUN=function(x)c(min(x),max(x)))

ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(xs))+theme_test(base_size=14)
# 
# 
# 
# ####proportions to extend time series
# #legal weight
# xs$propL = xs$Legal_wt/xs$WEIGHT_KG
# cr = mean(xs$propL,na.rm=T)
#  ggplot(subset(xs,YEAR>1998),aes(propL))+geom_density(color="darkblue", fill="lightblue")+facet_wrap(~YEAR)+xlab('Proportion Legal Lobster Weight')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =cr,col='red' ,linewidth=1.)
#  cw = ggplot(xs,aes(propL))+geom_density(color="darkblue", fill="lightblue")+xlab('Proportion Legal Lobster (wt)')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =cr,col='red' ,linewidth=1.)
# ##legal N
# ri = aggregate(value~mission+setno,data=subset(xle,Length %in% 82:300),FUN=sum)
# names(ri)[3] = 'CommT'
# ti = aggregate(value~mission+setno+YEAR,data=xle,FUN=sum)
# rti=merge(ri,ti)
# rti$propcomm = (rti$CommT/rti$value)
# rti = subset(rti,is.finite(propcomm))
# commMult = mean(rti$propcomm)
# 
# ggplot(subset(rti,YEAR>1998),aes(propcomm))+geom_density(color="darkblue", fill="lightblue")+facet_wrap(~YEAR)+xlab('Proportion Commercial Lobster (N)')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =commMult,col='red' ,linewidth=1.)
# cN = ggplot(rti,aes(propcomm))+geom_density(color="darkblue", fill="lightblue")+xlab('Proportion Commercial Lobster (n)')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =commMult,col='red' ,linewidth=1.)
# 
# rtim = merge(rti,subset(x,select=c(mission,setno,LONGITUDE,LATITUDE)))
# saveRDS(rtim,'rvproportions82_300.rds')
# 
# 
# 
# ##recruitment
# ri = aggregate(value~mission+setno,data=subset(xle,Length %in% 70:82),FUN=sum)
# names(ri)[3] = 'recT'
# rti=merge(ri,ti)
# rti$propRec = (rti$recT/rti$value)
# rti = subset(rti,is.finite(propRec))
# recMult = mean(rti$propRec)
# ggplot(subset(rti,YEAR>1998),aes(propRec))+geom_density(color="darkblue", fill="lightblue")+facet_wrap(~YEAR)+xlab('Proportion Recruit Lobster (N)')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =commMult,col='red' ,linewidth=1.)
# rN = ggplot(rti,aes(propRec))+geom_density(color="darkblue", fill="lightblue")+xlab('Proportion Recruit Lobster (n)')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =recMult,col='red' ,linewidth=1.)
# 
# 
# ##recruitment + comm
# ri = aggregate(value~mission+setno,data=subset(xle,Length %in% 70:300),FUN=sum)
# 
# names(ri)[3] = 'recT'
# rti=merge(ri,ti)
# rti$propRec = (rti$recT/rti$value)
# rti = subset(rti,is.finite(propRec))
# 
# recMult = mean(rti$propRec)
# ggplot(subset(rti,YEAR>1998),aes(propRec))+geom_density(color="darkblue", fill="lightblue")+facet_wrap(~YEAR)+xlab('Proportion Recruit + Commercial Lobster (N)')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =commMult,col='red' ,linewidth=1.)
# rNC = ggplot(rti,aes(propRec))+geom_density(color="darkblue", fill="lightblue")+xlab('Proportion Commerical + Recruit Lobster (n)')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept =recMult,col='red' ,linewidth=1.)
# 
# rtim = merge(rti,subset(x,select=c(mission,setno,LONGITUDE,LATITUDE)))
# saveRDS(rtim,'rvproportions70_300.rds')
# ggpubr::ggarrange(cw,cN,rN,rNC)
# 
# 
# ################################################################################################################################################
# #####ILTS


y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = F)
y = subset(y,select=c(TRIP_ID,SET_NO, SET_LONG,SET_LAT,SET_DATE,SA_CORRECTED_PRORATED_N))    

survey = y
survey = subset(survey, month(SET_DATE) %in% c(6:8))
survey$YEAR = year(survey$SET_DATE)
survey = st_as_sf(survey,coords = c('SET_LONG','SET_LAT'),crs=4326)
survey$CommN = survey$SA_CORRECTED_PRORATED_N

pt = subset(survey,YEAR==1996 )
pt$YEAR='1995-2012'
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(pt),fw='~YEAR')
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(survey,YEAR %in% 2013:2023 ))


yp = ILTS_ITQ_All_Data(redo_base_data = F,size=c(0,300),aggregate=F,species=2550,biomass = F,extend_ts = F)
yp1 = subset(yp, month(SET_DATE) %in% 6:8 & YEAR<2024 & YEAR>2012)

yp1 = aggregate(SA_CORRECTED_PRORATED_N~YEAR+FISH_LENGTH,data=yp1,FUN=mean)
ggplot(yp1,aes(x=FISH_LENGTH,y=SA_CORRECTED_PRORATED_N))+geom_bar(stat='identity')+facet_wrap(~YEAR)+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')


yp2 = subset(yp, month(SET_DATE)>8)

yp2 = aggregate(SA_CORRECTED_PRORATED_N~YEAR+LFA+FISH_LENGTH,data=yp2,FUN=mean)
ggplot(subset(yp2,YEAR>2014),aes(x=FISH_LENGTH,y=SA_CORRECTED_PRORATED_N))+geom_bar(stat='identity')+facet_wrap(~YEAR)+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')




#################################################################################################################################
###fall survey 
y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = F)
survey$YEAR = year(survey$SET_DATE)
survey = subset(y, month(SET_DATE) >8 & YEAR>2014)
survey = st_as_sf(survey,coords = c('SET_LONG','SET_LAT'),crs=4326)
survey$CommN = survey$SA_CORRECTED_PRORATED_N
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(survey,YEAR %in% 2019:2023 ),fw='~YEAR')


###############################################################################################################################

#including DMNR

z = lobster.db('DMR.redo')
br = round(seq(min(z$Year),max(z$Year),length=4))
lab = c(paste(br[1],br[2],sep="-"),paste(br[2]+1,br[3],sep="-"),paste(br[3]+1,br[4],sep="-"))
z$rYR=cut(z$Year,breaks=br,labels=lab,include.lowest = T,include.highest=F)
#to check
aggregate(Year~rYR,data=z,FUN=function(x)c(min(x),max(x)))

ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=z,fw='~rYR')

xle = z %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
#xle= na.zero(xle,cols='value')
xxa = aggregate(value~Length+Year,data=xle,FUN=mean)
ggplot(subset(xxa),aes(x=Length,y=value))+geom_bar(stat='identity')+facet_wrap(~Year,scales = 'free_y')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')



###NEFSC

n = NEFSC_sets()
n = subset(n,month(DATE)>8 & LONGITUDE> -68 & LATITUDE>42)
ns = st_as_sf(n,coords=c('LONGITUDE','LATITUDE'),crs=4326)
ggplot(ns)+geom_sf()+geom_sf(data=ns_coast,fill='black')+geom_sf(data=z,colour='red')+geom_sf(data=ys,colour='blue')



xle = n %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
#xle= na.zero(xle,cols='value')
xxa = aggregate(value~Length+YEAR,data=xle,FUN=mean)
ggplot(subset(xxa, YEAR %in% floor(seq(min(xxa$YEAR),max(xxa$YEAR),length=12))),aes(x=Length,y=value))+geom_bar(stat='identity')+facet_wrap(~YEAR,scales='free_y')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')

br = round(seq(min(ns$YEAR),max(ns$YEAR),length=9))
lab = c(paste(br[1],br[2],sep="-"),paste(br[2]+1,br[3],sep="-"),paste(br[3]+1,br[4],sep="-"),paste(br[4]+1,br[5],sep="-"),paste(br[5]+1,br[6],sep="-"),paste(br[6]+1,br[7],sep="-"),paste(br[7]+1,br[8],sep="-"),paste(br[8]+1,br[9],sep="-"))
ns$rYR=cut(ns$YEAR,breaks=br,labels=lab,include.lowest = T,include.highest=F)
#to check
aggregate(YEAR~rYR,data=ns,FUN=function(x)c(min(x),max(x)))

ggLobsterMap('custom',xlim=c(-68,-63.5),ylim=c(42,46),addGrids = F,addPoints = T,pts=subset(ns),fw='~rYR')+theme_test(base_size=14)


#######################

#commercial weight proportion

y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = T,extend_ts = F)
y = subset(y, tot>0 & month(y$SET_DATE)<8 & LFA %in% c('L34','L35','L36','L37','L38'),select=c(TRIP_ID,SET_NO,SET_LONG,SET_LAT,SA_CORRECTED_PRORATED_N,tot))

x = RV_sets()
x = subset(x,month(x$DATE) %in% c(6,7,8))
xs = st_as_sf(x,coords = c('LONGITUDE','LATITUDE'),crs=4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326

xs = st_join(xs,rL,join=st_within)
xs$X=st_coordinates(xs)[,1]
xs$Y=st_coordinates(xs)[,2]
st_geometry(xs) <- NULL
xs = subset(xs,LFA %in% c(34,35,36,37,38,41,40,33) &!is.na(Legal_wt)& X< -64,select=c(mission,setno,X,Y,Legal_wt,WEIGHT_KG))

names(y)=names(xs)
y$mission = as.character(y$mission)
y$setno = as.character(y$setno)

xy = bind_rows(y,xs)
saveRDS(xy,'proportions_wt_82_300.rds')


#commercial numb proportion

y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = F,extend_ts = F)
y = subset(y, tot>0 & month(y$SET_DATE)<8 & LFA %in% c('L34','L35','L36','L37','L38'),select=c(TRIP_ID,SET_NO,SET_LONG,SET_LAT,SA_CORRECTED_PRORATED_N,tot))

x = RV_sets()
x = subset(x,month(x$DATE) %in% c(6,7,8))

xs = st_as_sf(x,coords = c('LONGITUDE','LATITUDE'),crs=4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326

xs = st_join(xs,rL,join=st_within)
xs$X=st_coordinates(xs)[,1]
xs$Y=st_coordinates(xs)[,2]
st_geometry(xs) <- NULL
xs = subset(xs,LFA %in% c(34,35,36,37,38,41,40,33) &!is.na(Legal_wt)& X< -64,select=c(mission,setno,X,Y,Legal,Lobster))

names(y)=names(xs)
y$mission = as.character(y$mission)
y$setno = as.character(y$setno)

xy = bind_rows(y,xs)
saveRDS(xy,'proportions82_300.rds')


#recruit numb proportion

y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(70,82),aggregate=T,species=2550,biomass = F,extend_ts = F)
y = subset(y, tot>0 & month(y$SET_DATE)<8 & LFA %in% c('L34','L35','L36','L37','L38'),select=c(TRIP_ID,SET_NO,SET_LONG,SET_LAT,SA_CORRECTED_PRORATED_N,tot))

x = RV_sets()
x = subset(x,month(x$DATE) %in% c(6,7,8))

xle = x %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
#xle= na.zero(xle,cols='value')
xxa = aggregate(value~mission+setno+Lobster+LATITUDE+LONGITUDE,data=subset(xle,Length %in% 70:82),FUN=sum)

xs = st_as_sf(xxa,coords = c('LONGITUDE','LATITUDE'),crs=4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326

xs = st_join(xs,rL,join=st_within)
xs$X=st_coordinates(xs)[,1]
xs$Y=st_coordinates(xs)[,2]
st_geometry(xs) <- NULL
xs = subset(xs,LFA %in% c(34,35,36,37,38,41,40,33) & X< -64,select=c(mission,setno,X,Y,value,Lobster))
names(xs)[5] = 'Recruits'
names(y)=names(xs)
y$mission = as.character(y$mission)
y$setno = as.character(y$setno)

xy = bind_rows(y,xs)
saveRDS(xy,'proportions70_82.rds')


######################################
#recruit numb proportion

y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(70,300),aggregate=T,species=2550,biomass = F,extend_ts = F)
y = subset(y, tot>0 & month(y$SET_DATE)<8 & LFA %in% c('L34','L35','L36','L37','L38'),select=c(TRIP_ID,SET_NO,SET_LONG,SET_LAT,SA_CORRECTED_PRORATED_N,tot))

x = RV_sets()
x = subset(x,month(x$DATE) %in% c(6,7,8))

xle = x %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
#xle= na.zero(xle,cols='value')
xxa = aggregate(value~mission+setno+Lobster+LATITUDE+LONGITUDE,data=subset(xle,Length %in% 70:300),FUN=sum)

xs = st_as_sf(xxa,coords = c('LONGITUDE','LATITUDE'),crs=4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326

xs = st_join(xs,rL,join=st_within)
xs$X=st_coordinates(xs)[,1]
xs$Y=st_coordinates(xs)[,2]
st_geometry(xs) <- NULL
xs = subset(xs,LFA %in% c(34,35,36,37,38,41,40,33) & X< -64,select=c(mission,setno,X,Y,value,Lobster))
names(xs)[5] = 'RecruitsComm'
names(y)=names(xs)
y$mission = as.character(y$mission)
y$setno = as.character(y$setno)

xy = bind_rows(y,xs)
saveRDS(xy,'proportions70_300.rds')


