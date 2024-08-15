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

aT = readRDS('CPUETempDepth.rds')
aT = subset(aT,SYEAR>2005 & SYEAR<2023)

#start with an overall CPUE trend

aa = split(aT,f=list(aT$LFA,aT$SYEAR))
cpue.lst<-list()
cpue.ann<- list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS','temp')]
  names(tmp)<-c('time','catch','effort','temp')
  tmp$date<-as.Date(tmp$time)
  first.day<-min(tmp$date)
  tmp$time<-julian(tmp$date,origin=first.day-1)
  g<-as.data.frame(biasCorrCPUE(tmp,by.time=T))
  g$lfa=unique(aa[[i]]$LFA)
  g$yr = unique(aa[[i]]$SYEAR)
  te = aggregate(temp~time, data=tmp,FUN=mean)
  g = merge(g,te,by.x='t',by.y='time')
  gl = aggregate(effort~time, data=tmp, FUN=sum)
  g = merge(g,gl,by.x='t',by.y='time')
  cpue.lst[[i]] <- g
  
  g<-as.data.frame(t(biasCorrCPUE(tmp,by.time=F)))
  g=data.frame(g,LFA=as.numeric(unique(aa[[i]]$LFA)),SYEAR=as.numeric(unique(aa[[i]]$SYEAR)),temp=as.numeric(aggregate(temp~1, data=tmp,FUN=mean)))
  cpue.ann[[i]]=g
}

ca =as.data.frame(do.call(rbind,cpue.ann))
ca$temp = round(ca$temp,1)
ggplot(subset(ca,LFA %in% c(35,36,38)),aes(x=SYEAR,y=unBCPUE))+geom_point()+geom_line()+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.3)+
  facet_wrap(~LFA)+xlab('Fishing Year (ending)')+ylab('unbiased CPUE')+theme_test()

cc =as.data.frame(do.call(rbind,cpue.lst))
cc$dyear = cc$yr+cc$t/365
cc$temp = round(cc$temp,1)
#cc = merge(cc,preT)

ggplot(subset(cc,lfa %in% c(35) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()


ggplot(subset(cc,lfa %in% c(36) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()

ggplot(subset(cc,lfa %in% c(38) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()


saveRDS(list(ca,cc),file='unBIASED_CPUE.rds')


###beginning models

aa = split(aT,f=list(aT$LFA,aT$SYEAR))
oo = list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  first.day<-min(tmp$DATE_FISHED)
  tmp$time<-julian(tmp$DATE_FISHED,origin=first.day-1)
  oo[[i]] = tmp
}

aT = do.call(rbind,oo)

aT$fYear = as.factor(aT$SYEAR)
aT$leffort = log(aT$NUM_OF_TRAPS)
l38 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(aT,LFA==38),family = Gamma(link='log'),method='REML')
l38a = gam(WEIGHT_KG~s(time)+fYear+offset(leffort),data=subset(aT,LFA==38),family = Gamma(link='log'),method='REML')
l38b = gam(WEIGHT_KG~s(time)+s(time,by=fYear)+fYear+offset(leffort),data=subset(aT,LFA==38),family = Gamma(link='log'),method='REML')

saveRDS(list(l38,l38a,l38b),'first3CPUEmodels.rds')

aT$rLID = as.factor(aT$LICENCE_ID)
l38c = gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rLID,bs='re')+offset(leffort),data=subset(aT,LFA==38),family = Gamma(link='log')) #best by AIC
saveRDS(l38c,'CPUEmodel_LicenceRE.rds')


aT$rV = as.factor(aT$VESSEL_NAME)
l38d = gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rV,bs='re')+offset(leffort),data=subset(aT,LFA==38),family = Gamma(link='log')) #best by AIC
saveRDS(l38c,'CPUEmodel_vesselRE.rds')


#plot index
#these are conditional
pr1 = ggpredict(l38,terms='fYear',condition = c(leffort=0))
pr2 = ggpredict(l38a,terms=c('fYear'),condition = c(leffort=0, time=20))
pr3 = ggpredict(l38b,terms='fYear',condition = c(leffort=0,time=20))
pr4 = ggpredict(l38c,terms='fYear',condition = c(leffort=0,time=20))
pr5 = ggpredict(l38d,terms='fYear',condition = c(leffort=0,time=20))

pr1t = data.frame(x=pr1$x,predicted=pr1$predicted,Model='Base',conf.low=pr1$conf.low,conf.high=pr1$conf.high,se = pr1$std.error)
pr2t = data.frame(x=pr2$x,predicted=pr2$predicted,Model='BD',conf.low=pr2$conf.low,conf.high=pr2$conf.high,se = pr2$std.error)
pr3t = data.frame(x=pr3$x,predicted=pr3$predicted,Model='BxD',conf.low=pr3$conf.low,conf.high=pr3$conf.high,se = pr3$std.error)
pr4t = data.frame(x=pr4$x,predicted=pr4$predicted,Model='BxDRL',conf.low=pr4$conf.low,conf.high=pr4$conf.high,se = pr4$std.error)
pr5t = data.frame(x=pr5$x,predicted=pr5$predicted,Model='BxDRV',conf.low=pr5$conf.low,conf.high=pr5$conf.high,se = pr5$std.error)

pp = rbind(pr5, pr4t,pr3t,pr2t,pr1t)
pp$x = as.numeric(as.character(pp$x))

ggplot(pp, aes(x = x, y = predicted,colour=Model,fill=Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Year')+
  ylab('Standardized CPUE')+
  theme_test(base_size = 14)


# ##these are marginal means and have equal weighting across all days of the season... which is why they are lower
# require(ggeffects) #new as of 
# b = emmeans::emmeans(l38,~fYear,offset=0,type='response',data=subset(aT,LFA==38))
# bmm = as.data.frame(summary(b))[c('response', 'lower.CL','upper.CL')]
# bmm$Model = 'base'\
# bmm$Year = 2006:2022
# d = emmeans::emmeans(l38a,~fYear+time,offset=0,type='response',data=subset(aT,LFA==38),weights="proportional")
# dmm = as.data.frame(summary(d))[c('response', 'lower.CL','upper.CL')]
# dmm$Model = 'BD'
# dmm$Year = 2006:2022
# f = emmeans::emmeans(l38b,~fYear,offset=0,type='response',data=subset(aT,LFA==38))
# fmm = as.data.frame(summary(f))[c('response', 'lower.CL','upper.CL')]
# fmm$Model = 'BxD'
# fmm$Year = 2006:2022
# oo = do.call(rbind,list(bmm,dmm,fmm))

##
packageurl = "https://cran.r-project.org/src/contrib/Archive/ggeffects/ggeffects_1.3.4.tar.gz"

install.packages(packageurl, repos=NULL, type="source")

#turning conditional into marginal weighting by prop of total catch
pr2 = ggpredict(l38a,terms=c('fYear','time[1:235]'),condition = c(leffort=0), back_transform = F)
pr2tab = data.frame(fYear=pr2$x,predicted=pr2$predicted,Model='BD',se=pr2$std.error,time=pr2$group)

#weighted by catch -- not what we are looking for
# ind = aggregate(WEIGHT_KG~time+SYEAR,data=subset(aT,LFA ==38),FUN=sum)
# ind1 = aggregate(WEIGHT_KG~SYEAR,data=subset(aT,LFA ==38),FUN=sum)
# names(ind1)[2] = 'SumWt'
# ind = merge(ind,ind1)
# ind$prop = ind$WEIGHT_KG/ind$SumWt
# ind$fYear=as.factor(ind$SYEAR)
# ppr = merge(pr2ta,ind[,c('fYear','time','prop')])
# ppr$wpre = ppr$predicted*ppr$prop
# ppr$wvar = (ppr$se)^2*ppr$prop
# p2 = aggregate(cbind(wpre,wvar)~fYear+Model,data=ppr,FUN=sum)
# p2$wse = sqrt(p2$wvar)
# 
# pr3 = ggpredict(l38b,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
# pr3ta = data.frame(fYear=pr3$x,predicted=pr3$predicted,Model='BxD',se=pr3$std.error,time=pr3$group)
# pp3r = merge(pr3ta,ind[,c('fYear','time','prop')])
# pp3r$wpre = pp3r$predicted*pp3r$prop
# pp3r$wvar = (pp3r$se)^2*pp3r$prop
# p3 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp3r,FUN=sum)
# p3$wse = sqrt(p3$wvar)
# 
# plot(as.numeric(as.character(pr1t$x)),pr1t$predicted,type='b',ylim=c(0,4))
# lines(as.numeric(as.character(p2$fYear)),p2$wpre,type='b',col='red')
# lines(as.numeric(as.character(p3$fYear)),p3$wpre,type='b',col='blue')

##marginal weighted by effort

ind = aggregate(NUM_OF_TRAPS~time+SYEAR,data=subset(aT,LFA ==38),FUN=sum)
ind1 = aggregate(NUM_OF_TRAPS~SYEAR,data=subset(aT,LFA ==38),FUN=sum)
names(ind1)[2] = 'SumTraps'
ind = merge(ind,ind1)
ind$prop = ind$NUM_OF_TRAPS/ind$SumTraps
ind$fYear=as.factor(ind$SYEAR)
ppr = merge(pr2ta,ind[,c('fYear','time','prop')])
ppr$wpre = ppr$predicted*ppr$prop
ppr$wvar = (ppr$se)^2*ppr$prop
p2 = aggregate(cbind(wpre,wvar)~fYear+Model,data=ppr,FUN=sum)
p2$wse = sqrt(p2$wvar)

pr3 = ggpredict(l38b,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
pr3ta = data.frame(fYear=pr3$x,predicted=pr3$predicted,Model='BxD',se=pr3$std.error,time=pr3$group)
pp3r = merge(pr3ta,ind[,c('fYear','time','prop')])
pp3r$wpre = pp3r$predicted*pp3r$prop
pp3r$wvar = (pp3r$se)^2*pp3r$prop
p3 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp3r,FUN=sum)
p3$wse = sqrt(p3$wvar)

pr4 = ggpredict(l38c,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
pr4ta = data.frame(fYear=pr4$x,predicted=pr4$predicted,Model='BxDRL',se=pr4$std.error,time=pr4$group)
pp4r = merge(pr4ta,ind[,c('fYear','time','prop')])
pp4r$wpre = pp4r$predicted*pp4r$prop
pp4r$wvar = (pp4r$se)^2*pp4r$prop
p4 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp4r,FUN=sum)
p4$wse = sqrt(p4$wvar)

pr5 = ggpredict(l38d,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
pr5ta = data.frame(fYear=pr5$x,predicted=pr5$predicted,Model='BxDRV',se=pr5$std.error,time=pr5$group)
pp5r = merge(pr5ta,ind[,c('fYear','time','prop')])
pp5r$wpre = pp5r$predicted*pp5r$prop
pp5r$wvar = (pp5r$se)^2*pp5r$prop
p5 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp5r,FUN=sum)
p5$wse = sqrt(p5$wvar)



p1 = pr1t
names(p1)[c(1,2,6)] =c('fYear','wpre','wse')

pS = rbind(p1[,c(1,2,3,6)],p2[,c(1,2,3,5)],p3[,c(1,2,3,5)],p4[,c(1,2,3,5)],p5[,c(1,2,3,5)])
pS$Year = as.numeric(as.character(pS$fYear))

ggplot(pS, aes(x = Year, y = wpre ,colour=Model,fill=Model)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = wpre-wse, ymax = wpre+wse), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Year')+
  ylab('Standardized CPUE')+
  theme_test(base_size = 14)



####cross valiation
library(tidyverse)
library(modelr)
library(mgcv)

#100 samples
cv_df = crossv_mc(subset(aT,LFA==38), 50) 

#names samples
cv_df =
  cv_df |> 
  mutate(
    train = map(train, as_tibble),
    test = map(test, as_tibble))

# run
cv_res = 
  cv_df |> 
  mutate(
    sim_mod  = map(train, \(df) mgcv::gam(WEIGHT_KG~fYear+offset(leffort),data=df,family = Gamma(link='log'))),
    mod2  = map(train, \(df) mgcv::gam(WEIGHT_KG~fYear+s(time)+offset(leffort),data=df,family = Gamma(link='log'))), 
    mod3  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+offset(leffort),data=df,family = Gamma(link='log'))),
    mod4  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rLID,bs='re')+offset(leffort),data=df,family = Gamma(link='log'))), 
    mod5  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rV,bs='re')+offset(leffort),data=df,family = Gamma(link='log')))) |>
   mutate(
    rmse_M1 = map2_dbl(sim_mod, test, \(mod, df) rmse(model = mod, data = df)),
    rmse_M2 = map2_dbl(mod2, test, \(mod, df) rmse(model = mod, data = df)),
    rmse_M2 = map2_dbl(mod3, test, \(mod, df) rmse(model = mod, data = df)),
    rmse_M2 = map2_dbl(mod4, test, \(mod, df) rmse(model = mod, data = df)),
    rmse_M2 = map2_dbl(mod5, test, \(mod, df) rmse(model = mod, data = df)))

#plot cv
cv_res |> 
  select(starts_with("rmse")) |> 
  pivot_longer(
    everything(),
    names_to = "model", 
    values_to = "rmse",
    names_prefix = "rmse_") |> 
  mutate(model = fct_inorder(model)) |> 
  ggplot(aes(x = model, y = rmse)) + geom_violin()

####################################################
#not current analysis fall 2023; held for posterity
#####################################################
# ##Part 1, raw CPUE
# xx = subset(cc,lfa==38)
# l38P1 = gam(land~fYear+offset(leffort),data=xx,family = 'nb') #best by AIC
# pr1 = ggpredict(l38P1,terms='fYear',condition = c(leffort=0))
# summary(l38P1)
# AIC(l38P1)
# 
# plot(pr1)
# 
# ##Part 2, add in Day of Season
# xx = subset(cc,lfa==38)
# l38P2 = gam(land~fYear+s(t)+offset(leffort),data=xx,family = 'nb') #best by AIC
# pr2 = ggpredict(l38P2,terms=c('t','fYear'),condition = c(leffort=0))
# summary(l38P2)
# AIC(l38P2)
# ggplot(pr2, aes(x = x, y = predicted)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
#   facet_wrap(~group)+
#   xlab('Day of Season')+
#   ylab('CPUE')
# 
# ##Part 3, add in Day of Season and Temp
# xx = subset(cc,lfa==38)
# l38P3 = gam(land~fYear+s(t)+s(temp)+offset(leffort),data=xx,family = 'nb') #best by AIC
# pr3 = ggpredict(l38P3,terms=c('t','fYear'),condition = c(leffort=0,temp=8))
# summary(l38P3)
# AIC(l38P3)
# ggplot(pr3, aes(x = x, y = predicted)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
#   facet_wrap(~group)+
#   xlab('Day of Season')+
#   ylab('CPUE')
# 
# 
# ##Part 4, add in Day of Season and Temp and interaction
# xx = subset(cc,lfa==38)
# l38P4 = gam(land~fYear+s(t)+s(temp)+s(t,by=fYear)+offset(leffort),data=xx,family = 'nb') #best by AIC
# pr4 = ggpredict(l38P4,terms=c('t','fYear'),condition = c(leffort=0,temp=8))
# summary(l38P4)
# AIC(l38P4)
# ggplot(pr4, aes(x = x, y = predicted)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
#   facet_wrap(~group)+
#   xlab('Day of Season')+
#   ylab('CPUE')
# 
# l38cl = gam(land~fYear+s(t)+s(temp)+s(t,by=fYear)+offset(leffort),data=subset(cc,lfa==38),family = 'nb') #best by AIC
# 
# 
# ##
# l38Basel = gam(land~s(t)+s(temp)+offset(leffort),data=subset(cc,lfa==38),family = 'nb')
# 
# 
# require(ggeffects)
# mydft <- ggpredict(l38cl, terms = c('t'),condition=c(leffort=0,fYear=2016)) #
# ggplot(mydft, aes(x = x, y = predicted)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
#   xlab('Day of Season')+
#   ylab('CPUE')
# 
# mydf <- ggpredict(l38cl, terms = c('t','fYear'),condition=c(leffort=0,temp=8)) #
# plot(mydf,color='bw', ci=F)
# ggplot(mydf, aes(x = x, y = predicted)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
#   facet_wrap(~group)+
#   xlab('Day of Season')+
#   ylab('CPUE')
# 
# #residuals analysis
# bprd = ggpredict(l38Basel,terms=c('t'),condition=c(leffort=0,temp=8))
# bprdP = data.frame(baset = bprd$x,basepred = bprd$predicted)
# 
# mypred = data.frame(t=mydf$x,pred=mydf$predicted,yr=mydf$group)
# 
# combData = merge(mypred,bprdP,by.x='t',by.y='baset')
# combData$res = combData$pred - combData$basepred
# 
# ggplot(combData, aes(x = t, y = res)) +
#   geom_link2(aes(colour=after_stat(ifelse(y>0,'positive','negative'))),size=2,show.legend = F)+
#   facet_wrap(~yr)+
#   xlab('Day of Season')+
#   ylab('Residuals')+
#   labs(fill='Residual')+
#   geom_hline(yintercept=0,color='grey')
# 
# xf = lobster.db('seasonal.landings')
# xf$yr = as.numeric(substr(xf$SYEAR,6,9))
# 
# combDataL = merge(combData,xf[,c('yr','LFA38')])
# combDataL = combDataL[order(combDataL$LFA38,combDataL$t),]
# ggplot(combDataL, aes(x = t, y = res)) +
#   geom_link2(aes(colour=after_stat(ifelse(y>0,'positive','negative'))),size=2,show.legend = F)+
#   facet_wrap(~LFA38)+
#   xlab('Day of Season')+
#   ylab('Residuals')+
#   labs(fill='Residual')+
#   geom_hline(yintercept=0,color='grey')
# 
# 
# #what about temperatures for those years
# xxS = subset(xx,select=c(yr,temp,t))
# 
# combDataLT = merge(combDataL,xxS)
# combDataLT = combDataLT[order(combDataLT$LFA38,combDataLT$t),]
# avgT = aggregate(temp~t,data=combDataLT,FUN=mean)
# names(avgT)[2] = 'avgT'
# combDataLTa = merge(combDataLT,avgT)
# combDataLTa$tempDiff = combDataLTa$temp-combDataLTa$avgT
# 
# ggplot(combDataLTa, aes(x = t, y = tempDiff)) +
#   geom_link2(aes(colour=after_stat(ifelse(y>0,'positive','negative'))),size=2,show.legend = F)+
#   facet_wrap(~LFA38)+
#   xlab('Day of Season')+
#   ylab('Temp Diffs')+
#   geom_hline(yintercept=0,color='grey')
# 
# ##index 
# require(gratia)
# o = data_slice(l38P4,var1='fYear',offset=0)
# f  = fitted_values(l38P4,data=o,terms=c('(Intercept)','fYear'))
# ggplot(f,aes(x=fYear,y=fitted,group=1))+ geom_point()+ geom_errorbar(aes(ymin=lower,ymax=upper),width=.2)+stat_summary(fun.y=sum,geom='line')

############################
#spatial

aa = subset(aT,SYEAR>2005 & SYEAR<2023)
xx = split(aa,f=aa$SYEAR)
m=0
ou = list()
for(i in 1:length(xx)){
     w=   xx[[i]]
     y = unique(w$SYEAR)
     ww = unique(w$GRID_NUM)
    for(j in 1:length(ww)){
      m=m+1
      www = subset(w,GRID_NUM==ww[j])
      g = unique(www$GRID_NUM)
      uL = length(unique(www$LICENCE_ID))
      uLFA = length(unique(www$LFA))
      #grArea = as.numeric(sum(st_area(subset(gr,grid==ww[j]))))#these are wrong do no deal with holes right
      wa = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~GRID_NUM+DATE_FISHED+temp,data=www,FUN=sum)
      wa = wa[order(wa$DATE_FISHED),]
      wa$CPUE = wa$WEIGHT_KG / wa$NUM_OF_TRAPS
      wa$CUME = cumsum(wa$NUM_OF_TRAPS)
      wa$CUML = cumsum(wa$WEIGHT_KG)
      wa$rT = round(wa$temp,2)
      wa = merge(wa,preT,by.x='rT',by.y='temp')
      wa = wa[order(wa$DATE_FISHED),]
      wa$SCPUE = wa$pred*wa$CPUE
      wa$cpue=wa$SCPUE
      wa$effort = wa$NUM_OF_TRAPS
      wa$catch = wa$WEIGHT_KG
      wa$weight = wa$NUM_OF_TRAPS
      le = deluryLeslie(y=wa,estimate='leslie',method='robust',plot=F,weight = wa$weight)
      ou[[m]] = c(y,g,uL,uLFA,le$q,le$N0,max(wa$CUML),max(wa$CUME),nrow(www))
    }
  
}
out = as.data.frame(do.call(rbind,ou))
names(out) = c('Year','Grid','NLics','NLFAs','q','No','Land','Effort','Trips')
out$CPUE = out$Land/out$Effort

################ grouping grids
ou = list()
m=0
LFAgrid <- read.csv(file.path( project.datadirectory("bio.lobster"), "data","maps","GridPolys.csv"))
grL = (LFAgrid)
grL$area=grL$PID
attr(grL,'projection') <- 'LL'

#grs is the grids
#making the list of connections to grids
require(sp)
require(spdep)


  aa = subset(aT,SYEAR>2005 & SYEAR<2023)
  xx = split(aa,f=aa$SYEAR)
  m=0
  ou = list()
  for(i in 1:length(xx)){
    w=   xx[[i]]
    y = unique(w$SYEAR)
    ww = unique(w$GRID_NUM)
    grL1 = subset(grL,SID %in% ww)
    g1 = split(grL1,grL1$SID)
    nm = c()
    gp = list()
    for(k in 1:length(g1)){
      gp[[k]] = Polygons(list(Polygon(g1[[k]][,c('X','Y')])),unique(g1[[k]]$SID))
    }
    gpp = SpatialPolygons(gp,proj4string=CRS("+proj=longlat +datum=WGS84"))
    gpnb = poly2nb(gpp,row.names=names(gpp))
    names(gpnb)=names(gpp)
    
    for(j in 1:length(gpnb)){
      gr = names(gpnb)[gpnb[[j]]]
      ss = ifelse(length(gr)>=3,3,length(gr))
      gr = sample(gr,ss)
      m=m+1
      www = subset(w,GRID_NUM==gr)
      g = length(unique(www$GRID_NUM))
      gn = paste(gr,collapse="_")
      uL = length(unique(www$LICENCE_ID))
      uLFA = length(unique(www$LFA))
      bgf = subset(grL,SID %in%gr)
      attr(bgf,'projection') <-"LL"
      grArea = sum(calcArea(convUL(bgf))$area)
      wa = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~GRID_NUM+DATE_FISHED+temp,data=www,FUN=sum)
      wa = wa[order(wa$DATE_FISHED),]
      wa$CPUE = wa$WEIGHT_KG / wa$NUM_OF_TRAPS
      wa$CUME = cumsum(wa$NUM_OF_TRAPS)
      wa$CUML = cumsum(wa$WEIGHT_KG)
      wa$rT = round(wa$temp,2)
      wa = merge(wa,preT,by.x='rT',by.y='temp')
      wa = wa[order(wa$DATE_FISHED),]
      wa$SCPUE = wa$pred*wa$CPUE
      wa$cpue=wa$SCPUE
      wa$effort = wa$NUM_OF_TRAPS
      wa$catch = wa$WEIGHT_KG
      wa$weight = wa$NUM_OF_TRAPS
      le = deluryLeslie(y=wa,estimate='leslie',method='robust',plot=F,weight = wa$weight)
      ou[[m]] = c(y,gn,g,uL,uLFA,le$q,le$N0,max(wa$CUML),max(wa$CUME),nrow(www))
    }
    
  }
  outa = as.data.frame(do.call(rbind,ou))
  names(outa) = c('Year','Grid','NGrids','NLics','NLFAs','q','No','Land','Effort','Trips')
  outa$CPUE = outa$Land/outa$Effort
  outa = toNums(outa,c(1,3:10))
  outas = subset(outa,No>0)
  outas$Rexp = outas$Land/outas$No
  ggplot(outas,aes(x=Land,y=No))+geom_point()+facet_wrap(~Year,scales='free')
  ggplot(outas,aes(x=Year,y=Rexp))
  
  
head(outa)  
