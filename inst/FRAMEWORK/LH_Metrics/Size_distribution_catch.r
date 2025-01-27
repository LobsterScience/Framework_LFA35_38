# at sea sampling 
require(bio.lobster)
require(bio.utilities)
require(dplyr)
require(tidyr)
require(ggplot2)
require(devtools)
require(sf)


lobster.db('atSea')
v = subset(atSea,LFA %in% c(35,36,38))
v$UID = paste(v$TRIPNO,v$LFA,v$TRAPNO,v$STRINGNO,sep="_")
v$year = lubridate:::year(v$STARTDATE)
v$mm = lubridate::month(v$STARTDATE)
v$SYEAR = v$year
v$SYEAR = ifelse(v$mm %in% 10:12,v$SYEAR-1,v$SYEAR)
v$Season = ifelse(v$mm %in% c(10,11,12,1),'Fall',ifelse(v$mm %in% 2:4,'Winter',ifelse(v$mm %in% 5:7,'Spring','Outside')))

NTraps = aggregate(UID~LFA+SYEAR,data=v,FUN=function(x) length(unique(x)))
NTrapsSeason = aggregate(UID~LFA+SYEAR+Season,data=v,FUN=function(x) length(unique(x))) #most in fall or spring, separate indices

vF = subset(v,Season=='Fall')

vFLob = subset(v,Season=='Fall' & SPECIESCODE==2550)
ggplot(subset(vFLob,LFA==35),aes(x=CARLENGTH))+geom_histogram()+facet_wrap(~SYEAR,scales='free_y')+geom_vline(xintercept = 82.5,col='red')+theme_test


vfLobA = aggregate(CARLENGTH~LFA+SYEAR,data=vFLob,FUN=function(x) quantile(x,probs=c(0.05,.25,.5,.75,.95)))
ggplot(subset(vfLobA),aes(x=SYEAR,y=CARLENGTH[,3],ymin=CARLENGTH[,1],ymax=CARLENGTH[,5]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Carapace Length')+theme_test()
