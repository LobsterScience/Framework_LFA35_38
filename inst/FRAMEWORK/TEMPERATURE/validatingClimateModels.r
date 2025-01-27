#validating current data from models
#this was done for Climate change modelling

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
require(sf)

setwd(file.path(bio.datadirectory,'bio.lobster','analysis','ClimateModelling'))
	g = lobster.db('temperature.data')
	g$yr = year(g$T_DATE)
	g$W = ceiling(yday(g$T_DATE)/366*25)
	g$DOY = yday(g$T_DATE)
	#g = aggregate(TEMP~T_UID+LON_DD+LAT_DD+W+yr,data=g,FUN=median)
	g = g %>% st_as_sf(coords=c('LON_DD','LAT_DD'),crs=4326) %>% st_transform(32620)
	g$yr = year(g$T_DATE)
	g$W = ceiling(yday(g$T_DATE)/366*25)
	st_geometry(g) = st_geometry(g)/1000
	st_crs(g) = 32620
	d = dir(file.path(bio.datadirectory,'bio.lobster','Temperature Data/QuinnBT-2022/'),full.names = T)
	dd = d[grep('Can',d)]
	ee = d[grep('Had',d)]
	gg = d[grep('MPI',d)]

yy = 2016:2022
oo = list()
m=0
for(i in 1:length(yy)){
		k = subset(g,yr == yy[i])
		ddd = dd[grep(paste(yy[i],"_",sep=""),dd)]
		te = read.table(ddd,header=T)
		te = te %>% st_as_sf(coords = c("Longitude",'Latitude'),crs=4326) %>% st_transform(crs_utm20)
		st_geometry(te) = st_geometry(te)/1000
		st_crs(te) = 32620

		eee = ee[grep(paste(yy[i],"_",sep=""),ee)]
		fe = read.table(eee,header=T)
		fe = fe %>% st_as_sf(coords = c("Longitude",'Latitude'),crs=4326) %>% st_transform(crs_utm20)
		st_geometry(fe) = st_geometry(fe)/1000
		st_crs(fe) = 32620

		ggg = gg[grep(paste(yy[i],"_",sep=""),gg)]
		de = read.table(ggg,header=T)
		de = de %>% st_as_sf(coords = c("Longitude",'Latitude'),crs=4326) %>% st_transform(crs_utm20)
		st_geometry(de) = st_geometry(de)/1000
		st_crs(de) = 32620

		ll = unique(k$W)	
		for(l in 1:length(ll)){
			m = m+1
			kk = subset(k,W ==ll[l])
			tt = subset(te,Time==ll[l],select=c(Depth_m,BottomTemp))
			tt = bio.utilities::rename.df(tt,c('Depth_m','BottomTemp'),c('CanZ','CanBT'))
     		
     		ff = subset(fe,Time==ll[l],select=c(Depth_m,BottomTemp))
			ff = bio.utilities::rename.df(ff,c('Depth_m','BottomTemp'),c('HadZ','HadBT'))

    		pp = subset(de,Time==ll[l],select=c(Depth_m,BottomTemp))
			pp = bio.utilities::rename.df(pp,c('Depth_m','BottomTemp'),c('MPIZ','MPIBT'))
  		
       	  	 ou = st_nearest_feature(kk,tt)
       	 ds = st_distance(kk,tt[ou,],by_element=T)
       	 st_geometry(tt) = NULL
       	 kk$CanZ = tt$CanZ[ou]
       	 kk$CanBT = tt$CanBT[ou]
       	 kk$CanDist = as.numeric(ds)
       	 
       	 ou = st_nearest_feature(kk,ff)
       	 st_geometry(ff) = NULL
       	 kk$HadZ = ff$HadZ[ou]
       	 kk$HadBT = ff$HadBT[ou]
       	 
       	 ou = st_nearest_feature(kk,pp)
       	 st_geometry(pp) = NULL
       	 kk$MPIZ = pp$MPIZ[ou]
       	 kk$MPIBT = pp$MPIBT[ou]
       	 
       	 oo[[m]] = kk
		}
	}

xx = st_as_sf(do.call(rbind,oo))
xx = subset(xx,CanDist<quantile(CanDist,0.90))
xx$Haddiff = xx$HadBT - xx$TEMP
xx$Candiff = xx$CanBT - xx$TEMP 
xx$MPIdiff = xx$MPIBT - xx$TEMP 

ggplot(subset(xx,yr==2016)) + 
  geom_sf(aes(fill=Haddiff,color=Haddiff)) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  facet_wrap(~W) +
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()


###glorys
 
s = file.path(project.datadirectory('bio.lobster'),'Temperature Data','GLORYS','SummaryFiles')
k = dir(s,full.names=T)
k = k[grep('ShelfBoF',k)]
k= k[c(23:27,29,30,31)]

ol = list()
m=0
for(i in 1:length(k)){
			h = readRDS(k[i])
			h$Date = as.Date(h$Date)
			y = unique(year(h$Date))
			h$W  = ceiling(yday(h$Date)/366*25)
			h = aggregate(bottomT~W+X+Y,data=h,FUN=median)
			h = st_as_sf(h,coords =c('X',"Y"),crs=4326)
			h = h %>% st_as_sf(coords=c('X','Y'),crs=4326) %>% st_transform(32620)
			st_geometry(h) = st_geometry(h)/1000
			st_crs(h) = 32620

			x1 = subset(xx,yr == y)
			uW = unique(x1$W)
					for(j in 1:length(uW)){
				m=m+1
			
						nn = subset(x1,W==uW[j])
						hh = subset(h,W==uW[j])
						ou = st_nearest_feature(nn,hh)
       	 ds = st_distance(nn,hh[ou,],by_element=T)
       	 st_geometry(hh) = NULL
       	 nn$GlT = hh$bottomT[ou]
       	 nn$GlD = as.numeric(ds)
			ol[[m]] = nn       	

			}
}

xx = st_as_sf(do.call(rbind,ol))
xx=subset(xx,GlD<quantile(GlD,0.9,na.rm=T))

xx$Gldiff = xx$GlT - xx$TEMP

ggplot(xx,aes(Haddiff)) + geom_density(aes(y=..density..))+ geom_density(data=aes(y=..density..))

ggplot(xx,aes(Candiff)) + geom_histogram(aes(y=..density..)) + facet_wrap(~yr)	

ggplot(xx,aes(MPIdiff)) + geom_histogram(aes(y=..density..)) + facet_wrap(~yr)	

ggplot(xx,aes(Gldiff)) + geom_histogram(aes(y=..density..)) + facet_wrap(~yr)	


aggregate(cbind(Haddiff,Candiff,MPIdiff,Gldiff)~yr,data=xx,FUN=summary)

saveRDS(xx,file='ModelsToObserved.rds')
xx = readRDS('ModelsToObserved.rds')

##bnam
s = bnamR(redo=F)
loc = st_as_sf(s[[1]],coords=c('X','Y'),crs=4326)
loc = st_transform(loc,32620)
st_geometry(loc) = st_geometry(loc)/1000
st_crs(loc) = 32620
yy = 2016:2018
bt = s[[2]][,-1]
t = s[[3]]
out = list()
m=0
for(i in 1:length(yy)){
	w = subset(xx,yr==yy[i])
	l = grep(as.character(yy[i]),t)
	kl = unique(w$W)
	for(n in 1:length(kl)){
		m=m+1
			ww = subset(w,W==kl[n])
			j = st_nearest_feature(ww,loc)
			tt = ceiling(kl[n]*25/50)
			tt = ifelse(tt==13,12,tt)
			ww$bnamt = bt[j,l[tt]]
			ww$distBN = as.numeric(st_distance(ww,loc[j,],by_element=T))
    	out[[m]] = ww  
		}
	}

xx=do.call(rbind,out)
xx = subset(xx,distBN<4)
xx$bnamdiff = xx$bnamt - xx$TEMP

saveRDS(xx,file='ModelsToObserved_BNAM.rds')
xx = readRDS('ModelsToObserved_BNAM.rds')


###########################################################################################
##prune to LFA 34-38
sd = file.path( project.datadirectory("Framework_LFA35_38"), "figures")
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326
rL = st_transform(rL,32620) 
st_geometry(rL) <- st_geometry(st_as_sf(rL$geometry/1000)) 
st_crs(rL) <- 32620


glx = st_join(xx,rL,join=st_within)
xx = subset(glx,LFA %in% 34:38)


 mad(xx$Haddiff)
2.11
 mad(xx$Candiff)
 1.92
 mad(xx$MPIdiff)
 1.56
 mad(xx$bnamdiff)
 2.44
 mad(xx$Gldiff)
 2.01

 png(file.path(sd,'Temp2Modelcomparisons.png'))
  plot(density(xx$Candiff),xlab='Model - Measured',main="",xlim=c(-8,8),ylim=c(0,.27),lwd=2)
 lines(density(xx$Haddiff),col='red')
	lines(density(xx$MPIdiff),col='blue')
	lines(density(xx$bnamdiff),col='orange',lwd=2)
	lines(density(xx$Gldiff),col='purple',lwd=2)	
	abline(v=0,lwd=2,lty=3)
	legend('topleft',c('Can','Had','MPI','BNAM','Glorys'),lty=c(1,1,1,1,1),col=c('black','red','blue','orange','purple'),bty='n')
dev.off()

xx$Z = round(xx$CanZ/5)*5

require(tidyr)
xxc= xx %>% select(W,Z,LFA,Candiff,MPIdiff,Haddiff,Gldiff,bnamdiff) %>% pivot_longer(cols = ends_with('diff'))

xs = aggregate(cbind(Candiff,MPIdiff,Haddiff,bnamdiff,Gldiff)~Z,data=xx,FUN=median)
png(file.path(sd,'Temp2DepthModelcomparisons.png'))

with(xs,{
	plot(Z,Candiff,type='l',col='black',ylab='Diff from Observed Temp',ylim=c(-5,5))
	lines(Z,MPIdiff,type='l',col='blue')
	lines(Z,Haddiff,type='l',col='red')
	lines(Z,bnamdiff,type='l',col='orange')
	lines(Z,Gldiff,type='l',col='purple')
	
	})
	abline(h=0,lwd=2)
	legend('topright',c('Can','MPI','Had','BNAM','Glorys'),lty=c(1,1,1,1),col=c('black','blue','red','orange','purple'),bty='n')
dev.off()

ggplot(subset(xxc,W %in% 1:26),aes(x=factor(W),y=value))+geom_boxplot()+facet_wrap(~name)
ggplot(subset(xxc,name=='Gldiff'),aes(x=factor(W),y=value))+geom_boxplot()+facet_wrap(~LFA)

ggplot(subset(xxc,W %in% 13:21),aes(x=factor(W),y=value,colour=factor(name)))+geom_boxplot()+facet_wrap(~LFA)


xs = aggregate(cbind(Candiff,MPIdiff,Haddiff,bnamdiff,Gldiff)~W,data=xx,FUN=median)

png(file.path(sd,'Temp2WeekModelcomparisons.png'))

with(xs,{
  plot(W,Candiff,type='l',col='black',ylab='Diff from Observed Temp',ylim=c(-5,5))
  lines(W,MPIdiff,type='l',col='blue')
  lines(W,Haddiff,type='l',col='red')
  lines(W,bnamdiff,type='l',col='orange')
  lines(W,Gldiff,type='l',col='purple')
  
})
abline(h=0,lwd=2)
legend('topright',c('Can','MPI','Had','BNAM','Glorys'),lty=c(1,1,1,1),col=c('black','blue','red','orange','purple'),bty='n')
dev.off()



ggplot(subset(xx))+
  geom_sf(aes(fill=Gldiff,color=Gldiff)) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  facet_wrap(~W)+
  geom_sf(data=subset(rL,LFA %in% 34:38),size=1,color='black',fill=NA ) 
  

#############
##glorys only by day

###glorys
 
s = file.path(project.datadirectory('bio.lobster'),'Temperature Data','GLORYS','SummaryFiles')
k = dir(s,full.names=T)
k = k[grep('ShelfBoF',k)]
#k= k[c(23:27,29,30,31)]

ol = list()
m=0
for(i in 1:length(k)){
			h = readRDS(k[i])
			h$Date = as.Date(h$Date)
			y = unique(year(h$Date))
			h$DOY = yday(h$Date)
		#	h$W  = ceiling(yday(h$Date)/366*25)
			h = st_as_sf(h,coords =c('X',"Y"),crs=4326)
			h = h %>% st_as_sf(coords=c('X','Y'),crs=4326) %>% st_transform(32620)
			st_geometry(h) = st_geometry(h)/1000
			st_crs(h) = 32620

			x1 = subset(g,yr == y)
			uW = unique(x1$DOY)
					for(j in 1:length(uW)){
				m=m+1
			
						nn = subset(x1,DOY==uW[j])
						hh = subset(h,DOY==uW[j])
						ou = st_nearest_feature(nn,hh)
       	 ds = st_distance(nn,hh[ou,],by_element=T)
       	 st_geometry(hh) = NULL
       	 nn$GlT = hh$bottomT[ou]
       	 nn$GlD = as.numeric(ds)
			ol[[m]] = nn       	

			}
}

xx = st_as_sf(do.call(rbind,ol))

xx=subset(xx,GlD<quantile(GlD,0.95,na.rm=T))

xx$Gldiff = xx$GlT - xx$TEMP


ggplot(xx,aes(Gldiff)) + geom_histogram(aes(y=..density..)) + facet_wrap(~yr)	


aggregate(cbind(Gldiff)~yr,data=xx,FUN=summary)



rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326
rL = st_transform(rL,32620) 
st_geometry(rL) <- st_geometry(st_as_sf(rL$geometry/1000)) 
st_crs(rL) <- 32620


glo = st_join(xx,rL,join=st_within)
glo$mon = month(glo$T_DATE)
sd = file.path( project.datadirectory("Framework_LFA35_38"),'gloDay.rds')
saveRDS(glo,sd)

#Aug 7 2024
glo = readRDS(sd)
boxplot(Gldiff~mon,data=subset(glo,LFA %in% c(35,36,38) & yr %in% 1993:2022))
boxplot(Gldiff~mon,data=subset(glo,LFA %in% c(35,36,38) & yr %in% 2020:2022))

glo$X = st_coordinates(glo)[,1]
glo$Y = st_coordinates(glo)[,2]
glo$Date = format(as.POSIXct(glo$T_DATE,format='%Y-%m/%d %H:%M:%S'),format='%Y-%m-%d')
gloa = aggregate(cbind(GlT, TEMP)~mon+yr++X+Y+LFA+Date,data=glo,FUN=median)
gloa = st_as_sf(gloa,coords = c('X','Y'),crs=32620)
gloa$Gldiff = gloa$GlT - gloa$TEMP

aggregate(Gldiff~LFA+yr,data=subset(gloa,LFA %in% c(35,36,38) & mon %in% c(10:12,1:6) & yr %in% 1993:2022),FUN=length)

boxplot(Gldiff~yr,data=subset(gloa,LFA %in% c(38,35,36) & mon %in% c(10:12,1:6) & yr %in% 1993:2022))
abline(h=0,col='red')

#boxplot(Gldiff~mon,data=subset(glo,LFA %in% c(34) & yr %in% 2015:2019))
#boxplot(Gldiff~mon,data=subset(glo,LFA %in% c(34) & yr %in% 2020:2022))#

#bias correction surface BoF GOM
sf_use_s2(FALSE) #needed for cropping

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
crs_utm20 <- 32620
sf_use_s2(FALSE) #needed for cropping

ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -56.5, ymax = 47.5))))

ns_coast <- st_transform(ns_coast, crs_utm20)

st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
crs_utm20 <- 32620

gloGOM <- subset(glo,LFA %in% c(33,34,35,36,38,37,41,40)) %>%   
  st_as_sf()
 gloGOM$fc = as.factor(ifelse(gloGOM$yr %in% 2020:2022,1,0))
    
    ba = readRDS('~/git/bio.lobster.data/mapping_data/bathymetrySF.rds')
ba = ba %>% st_as_sf() 
st_geometry(ba) = st_geometry(ba)/1000
st_crs(ba) = 32620

				 ss = st_nearest_feature(gloGOM,ba)
       	 ds = st_distance(gloGOM,ba[ss,],by_element=T)
       	 st_geometry(ba) = NULL
       	 gloGOM$z = ba$z[ss]
       	 gloGOM$z_dist = as.numeric(ds)


surv_utm_coords <- st_coordinates(gloGOM)

gloGOM$X1000 <- surv_utm_coords[,1] 
gloGOM$Y1000 <- surv_utm_coords[,2] 
gloGOM$lz = log(gloGOM$z)
gloGOM = subset(gloGOM,!is.na(lz))
spde <- make_mesh(as_tibble(gloGOM), xy_cols = c("X1000", "Y1000"),
                   n_knots=200,type = "cutoff_search")
plot(spde)

# Add on the barrier mesh component:
bspde <- add_barrier_mesh(
  spde, ns_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)


fitBias = sdmTMB(Gldiff~ fc+ s(lz),
             data=as_tibble(gloGOM),
             time='mon', 
             mesh=bspde,
             family=gaussian(link='identity'),
             spatial='on',
             spatiotemporal='ar1')

saveRDS(list(data=gloGOM,mesh=bspde,model=fitBias),file='ClaimteModelBias.rds')
ou = readRDS('ClaimteModelBias.rds')
Glsur = readRDS('GlorysPredictSurfaceMonth.rds')
x = Glsur
x



x = bio.utilities::rename.df(x,c('bottomT','yr'),c('BT','YEAR'))
x = subset(x,z>0)
x$lZ = log(x$z)
x$X1000 = st_coordinates(x)[,1]
x$Y1000 = st_coordinates(x)[,2]
x = subset(x,exp(lZ)<400)

x = as_tibble(subset(x,select=c(Q,YEAR,BT,X1000,Y1000,lZ)))
x$geometry=NULL

g = predict(fitpa,newdata=subset(x,YEAR>1999))
