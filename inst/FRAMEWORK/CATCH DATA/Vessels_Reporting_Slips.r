###
require(bio.lobster)
require(bio.utilities)
require(devtools)
require(tidyr)
require(ggplot2)

#slips

b = lobster.db('process_slips')

bap = aggregate(WT_LBS~LFA+SYEAR+CFV_NO,data=b,FUN=sum)
bapa = aggregate(WT_LBS~LFA+SYEAR,data=bap,FUN=function(x) quantile(x,c(0.15,0.25,0.5,0.75,0.95)))
bb = bapal =data.frame(LFA=bapa$LFA,SYEAR=bapa$SYEAR,minL = bapa$WT_LBS[,1])

for(i in 1:nrow(bb)){
    j = which(bap$LFA==bb$LFA[i])
    k = which(bap$SYEAR==bb$SYEAR[i])
    l = which(bap$WT_LBS<bb$minL[i])
    m= intersect(intersect(j,k),l)
    bap = bap[-m,]
    
  
}


ggplot(subset(bapa,LFA %in% c(35,36,38)),aes(x=SYEAR,y=WT_LBS[,3],ymin=WT_LBS[,1],ymax=WT_LBS[,5]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA,scales='free_y')+xlab('Fishing Season')+ylab('Landings LBS ')


b$mn = lubridate::month(b$Date)
bap = aggregate(WT_LBS~LFA+SYEAR+CFV_NO,data=subset(b, mn %in% 10:12)  ,FUN=sum)

ba = aggregate(CFV_NO~LFA+SYEAR,data=bap,FUN=function(x)length(unique(x)))

ggplot(subset(ba,LFA %in% c(35) & SYEAR>1988),aes(x=SYEAR,CFV_NO))+geom_line()+facet_wrap(~LFA,scales = 'free_y')
       