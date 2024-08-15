#### comps of observations and GLORYsv12 on a date and location basis

require(bio.lobster)
require(bio.utilities)
require(devtools)
require(sf)
require(ggplot2)

setwd(file.path(project.datadirectory('Framework_LFA35_38')))

gr = readRDS('C:/Users/Cooka/Documents/git/bio.lobster.data/mapping_data/GridPolyLand.rds')
gr$grid_n = as.numeric(gr$grid)
gr = subset(gr,grid_n<300)
st_crs(gr) <- 4326
#gr = st_transform(gr,32620) 
#st_geometry(gr) <- st_geometry(st_as_sf(gr$geometry/1000)) 
#st_crs(gr) <- 32620
sf_use_s2(FALSE)

if(redo.obs){
  g = lobster.db('temperature.data')
  g$Date = as.Date(g$T_DATE)
  g = aggregate(TEMP~Date+LAT_DD+LON_DD,data=g,FUN=median)
  g = g %>% st_as_sf(coords=c('LON_DD','LAT_DD'),crs=4326)
  g = subset(g,year(Date)>2003)
  
  bb = st_as_sfc(st_bbox(gr))
#prune to outer bounds first
  cg = st_intersection(g,bb,join=st_within)
  g1 = st_join(cg,gr,join=st_within)
  ga2 = subset(g1,!is.na(grid_n))
  saveRDS(ga2,file.path('outputs','temperature','obstemps.rds'))
}
ga2 = readRDS(file.path('outputs','temperature','obstemps.rds'))
ga2$lon = st_coordinates(ga2)[,1]
ga2$lat = st_coordinates(ga2)[,2]
ga2 = aggregate(TEMP~date+lon+lat+GRIDNO,data=ga2,FUN=median)
out = list()

for(i in 1:length(ind)){
  k = read.csv(paste('GLORYS12v1_daily_bottomT_LFA33_AdamPoly',ind[i],'.csv',sep=""),header=F)
  names(k)=c('date',rep(paste('X',1:(ncol(k)-1),sep="")))
  k$date = as.Date(k$date,format = '%d-%b-%y')
  j = read.csv(paste('GLORYS12v1_daily_bottomT_LFA33_AdamPoly',ind[i],'_lonlat.csv',sep=""),sep=" ",header=T)
  v = list()
  for(z in 1:(nrow(j))){
    kk = k[,c(1,(z+1))]
    kk$lat = j[z,2]
    kk$lon = j[z,1]
    names(kk)[2]='BT'
    kk$EID= z
    kk$GridGroup = ind[i]
    v[[z]] = kk
  }      
  out[[i]] = do.call('rbind',v)
}

oo = do.call(rbind,out)
oos = st_as_sf(oo,coords = c('lon','lat'),crs=4326)
oos$GRIDNO = as.numeric(oos$GridGroup)
gas = st_as_sf(ga2,coords = c('lon','lat'),crs=4326)
gas$dist = gas$GL = NA
o = list()


for(i in 1:nrow(ga2)){
  p = gas[i,]
  l = subset(oos,date==p$date & GRIDNO==p$GRIDNO)
  distances <- st_distance(p,l)
  st_geometry(l) <-NULL
  gas[i,'GL']   <- l[which.min(distances),'BT']
  gas[i,'dist'] <- min(distances)
}


saveRDS(gas,'obstemps2GL.rds')
gas = readRDS('obstemps2GL.rds')
gas$mn = lubridate::month(gas$date)
gas$fm = ifelse(gas$mn %in% c(1:5,12),1,0)
gas$diff = gas$TEMP - gas$GL

gas = subset(gas,dist<5000)

saveRDS(gas,'obstemps2GL_filtered.rds')

ggplot(subset(gas,fm==1),aes(x=diff))+geom_histogram()+facet_wrap(~GRIDNO)
ggplot(subset(gas,fm==1),aes(x=diff,after_stat(density)))+geom_histogram()+geom_vline(aes(xintercept=0),col='red') +facet_wrap(~GRIDNO,scales = 'free_y')
gas$Gr=as.factor(as.character(paste('Gr',gas$GRIDNO,sep="")))

ggplot(subset(gas,fm==1),aes(x=factor(GRIDNO,levels=c('1','2','3','4','5','6','7','8','9','10')),y=diff))+
  geom_violin(width=1)+geom_boxplot(width=.1, color="grey", alpha=0.5)+geom_hline(aes(yintercept=0),col='red')+
  theme_bw()+xlab('Grid')+ylab('(Observed Temperature - GLORYSv12)')


gas = st_as_sf(gas)
g2 = st_read(file.path(project.datadirectory('bio.lobster'),'data','maps','GIS','LFA 33 Grid Grouping polygons_NAD83_region.shp'))
st_crs(g2) <- 4326
gg2 = st_centroid(g2)
gg2$X = st_coordinates(gg2)[,1]
gg2$Y = st_coordinates(gg2)[,2]

g2$area = st_area(g2)/1000000

v = list()
for(i in 1:length(ind)){
  j = read.csv(paste('GLORYS12v1_daily_bottomT_LFA33_AdamPoly',ind[i],'_lonlat.csv',sep=""),sep=" ",header=T)
  v[[i]] = j
}
vv = do.call(rbind, v)
vg = st_as_sf(vv,coords=c('longitude','latitude'),crs=4326)
require(ggtext)
ggplot(g2)+geom_sf()+geom_sf(data=subset(gas,fm==1),col='red')+geom_sf(data=g2,fill=NA,linewidth=1.3)+
  geom_sf(data=vg,col='blue')+geom_richtext(data=gg2,aes(x=X,y=Y,label=ID), label.padding = grid::unit(rep(0, 6), "pt"))+
  theme_bw()
