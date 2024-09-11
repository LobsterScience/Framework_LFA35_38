###########################
#start with an overview of all fishing year by year
require(bio.lobster)
require(bio.utilities)
require(dplyr)
require(devtools)
require(sf)
require(ggplot2)
require(mgcv)
require(ggeffects)
require(ggforce)
fd = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','CPUE')
setwd(fd)

gr = readRDS('C:/Users/Cooka/Documents/git/bio.lobster.data/mapping_data/GridPolyLand.rds')
st_crs(gr) <- 4326
gr = st_transform(gr,32620) 
st_geometry(gr) <- st_geometry(st_as_sf(gr$geometry/1000)) 
st_crs(gr) <- 32620

#need the temperature - catch relationship using IP/catchability/TempCatch.r output
#preT = readRDS(file=file.path(project.datadirectory('bio.lobster'),'analysis','ClimateModelling','tempCatchability.rds'))
#preT$temp = round(preT$Temperature,2)
#preT = aggregate(pred~temp,data=preT,FUN=mean)

aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>2005 & SYEAR<2024 & LFA %in% c(35,36,38))

#start with an overall CPUE trend

aa = split(aT,f=list(aT$LFA,aT$SYEAR))
cpue.lst<-list()
cpue.ann<- list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS')]
  names(tmp)<-c('time','catch','effort')
  tmp$date<-as.Date(tmp$time)
  first.day<-min(tmp$date)
  tmp$time<-julian(tmp$date,origin=first.day-1)
  g<-as.data.frame(biasCorrCPUE(tmp,by.time=T))
  g$lfa=unique(aa[[i]]$LFA)
  g$yr = unique(aa[[i]]$SYEAR)
  gl = aggregate(effort~time, data=tmp, FUN=sum)
  g = merge(g,gl,by.x='t',by.y='time')
  cpue.lst[[i]] <- g
  
  g<-as.data.frame(t(biasCorrCPUE(tmp,by.time=F)))
  g=data.frame(g,LFA=as.numeric(unique(aa[[i]]$LFA)),SYEAR=as.numeric(unique(aa[[i]]$SYEAR)))
  cpue.ann[[i]]=g
}

ca =as.data.frame(do.call(rbind,cpue.ann))
ggplot(subset(ca,LFA %in% c(35,36,38)),aes(x=SYEAR,y=unBCPUE))+geom_point()+geom_line()+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.3)+
  facet_wrap(~LFA)+xlab('Fishing Year (ending)')+ylab('unbiased CPUE')+theme_test()

cc =as.data.frame(do.call(rbind,cpue.lst))
cc$dyear = cc$yr+cc$t/365
#cc = merge(cc,preT)

ggplot(subset(cc,lfa %in% c(35) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()


ggplot(subset(cc,lfa %in% c(36) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()

ggplot(subset(cc,lfa %in% c(38) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()


saveRDS(list(ca,cc),file='unBIASED_CPUE.rds')


####logbook errors

a = lobster.db('process.logs')
b=lobster.db('process.logs.unfiltered')

aa = aggregate(SD_LOG_ID~SYEAR+LFA,data=a,FUN=function(x) length(unique(x)))
bb = aggregate(SD_LOG_ID~SYEAR+LFA,data=b,FUN=function(x) length(unique(x)))

names(bb)[3]='SD_LOG_UNF'

ab = merge(aa,bb)
ab$prop = 1-ab$SD_LOG_ID/ab$SD_LOG_UNF

ab=subset(ab,LFA %in% c(35,36,38) & SYEAR>2005 & SYEAR<2024)
agg = aggregate(prop~LFA,data=ab,FUN=mean)
ggplot(ab,aes(x=SYEAR,y=prop))+geom_line()+facet_wrap(~LFA)+theme_test(base_size = 14)+xlab('Fishing Season')+ylab('Proportion Logs Containing Errors')+geom_hline(data=agg,aes(yintercept=prop),color='red')




