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


###beginning models

aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>2005 & SYEAR<2024 & LFA %in% c(35,36,38))

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
w = readRDS('first3CPUEmodels.rds')
l38 = w[[1]]
l38a = w[[2]]
l38b = w[[3]]

aT$rLID = as.factor(aT$LICENCE_ID)
l38c = gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rLID,bs='re')+offset(leffort),data=subset(aT,LFA==38),family = Gamma(link='log'),method='REML') #random intercept 
saveRDS(l38c,'CPUEmodel_LicenceRE.rds')
l38c = readRDS('CPUEmodel_LicenceRE.rds')


aT$rV = as.factor(aT$VR_NUMBER)
l38d = gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rV,bs='re')+offset(leffort),data=subset(aT,LFA==38),family = Gamma(link='log'),method='REML') # random intercept
saveRDS(l38d,'CPUEmodel_vesselRE.rds')
l38d = readRDS('CPUEmodel_vesselRE.rds')


# #plot index
# #these are conditional
# pr1 = ggpredict(l38,terms='fYear',condition = c(leffort=0))
# pr2 = ggpredict(l38a,terms=c('fYear'),condition = c(leffort=0, time=20))
# pr3 = ggpredict(l38b,terms='fYear',condition = c(leffort=0,time=20))
# pr4 = ggpredict(l38c,terms='fYear',condition = c(leffort=0,time=20))
# pr5 = ggpredict(l38d,terms='fYear',condition = c(leffort=0,time=20))
# 
# pr1t = data.frame(x=pr1$x,predicted=pr1$predicted,Model='Base',conf.low=pr1$conf.low,conf.high=pr1$conf.high,se = pr1$std.error)
# pr2t = data.frame(x=pr2$x,predicted=pr2$predicted,Model='BD',conf.low=pr2$conf.low,conf.high=pr2$conf.high,se = pr2$std.error)
# pr3t = data.frame(x=pr3$x,predicted=pr3$predicted,Model='BxD',conf.low=pr3$conf.low,conf.high=pr3$conf.high,se = pr3$std.error)
# pr4t = data.frame(x=pr4$x,predicted=pr4$predicted,Model='BxDRL',conf.low=pr4$conf.low,conf.high=pr4$conf.high,se = pr4$std.error)
# pr5t = data.frame(x=pr5$x,predicted=pr5$predicted,Model='BxDRV',conf.low=pr5$conf.low,conf.high=pr5$conf.high,se = pr5$std.error)
# 
# pp = rbind(pr5t, pr4t,pr3t,pr2t,pr1t)
# pp$x = as.numeric(as.character(pp$x))
# 
# ggplot(pp, aes(x = x, y = predicted,colour=Model,fill=Model)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
#   #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
#   xlab('Year')+
#   ylab('Standardized CPUE')+
#   theme_test(base_size = 14)
# 
# 
# 
# 
# 
# ###emmeans which is the mean for the year assuming equal weights per day--- not really the case in observational studies
# b = emmeans::emmeans(l38,~fYear,offset=0,type='response',data=subset(aT,LFA==38))
# bmm = as.data.frame(summary(b))[c('response', 'lower.CL','upper.CL')]
# bmm$Model = 'Base'
# bmm$Year = 2006:2023
# d = emmeans::emmeans(l38a,~fYear,offset=0,type='response',data=subset(aT,LFA==38),at=list(time=1:234))
# dmm = as.data.frame(summary(d))[c('response', 'lower.CL','upper.CL')]
# dmm$Model = 'BD'
# dmm$Year = 2006:2023
# f = emmeans::emmeans(l38b,~fYear,offset=0,type='response',data=subset(aT,LFA==38))
# fmm = as.data.frame(summary(f))[c('response', 'lower.CL','upper.CL')]
# fmm$Model = 'BxD'
# fmm$Year = 2006:2023
# g = emmeans::emmeans(l38c,~fYear,offset=0,type='response',data=subset(aT,LFA==38))
# gmm = as.data.frame(summary(g))[c('response', 'lower.CL','upper.CL')]
# gmm$Model = 'BxDRL'
# gmm$Year = 2006:2023
# 
# h = emmeans::emmeans(l38d,~fYear,offset=0,type='response',data=subset(aT,LFA==38))
# hmm = as.data.frame(summary(h))[c('response', 'lower.CL','upper.CL')]
# hmm$Model = 'BxDRV'
# hmm$Year = 2006:2023
# 
# 
# 
# 
# oo = do.call(rbind,list(bmm,dmm,fmm,gmm,hmm))
# 
# #marginal mean, internally consistent
# 
# ggplot(oo, aes(x = Year, y = response ,colour=Model,fill=Model)) +
#   geom_line(linewidth=1) +
#   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width=0.2)+
#   #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
#   xlab('Year')+
#   ylab('Marginal Mean CPUE')+
#   theme_test(base_size = 14)
# 


# trip weightings


ind = aggregate(SD_LOG_ID~time+SYEAR,data=subset(aT,LFA ==38),FUN=function(x) length(unique(x)))
ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
names(ind1)[2] = 'SumTrips'
ind = merge(ind,ind1)
ind$prop = ind$SD_LOG_ID/ind$SumTrips
ind$fYear=as.factor(ind$SYEAR)


b = emmeans::emmeans(l38,~fYear,offset=0,type='response',data=subset(aT,LFA==38))
bmm = as.data.frame(summary(b))[c('response', 'SE')]
bmm$Model = 'Base'
bmm$SYEAR = 2006:2023
names(bmm)[1:2]=c('wemm','wse')

d = emmeans::emmeans(l38a,~fYear+time,offset=0,type='response',data=subset(aT,LFA==38),at=list(time=1:234))

dmm = as.data.frame(summary(d))
dmm = merge(dmm,ind)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){

      z = subset(dmm,Year==ii[i])
      usm = sum(z$SD_LOG_ID)
      dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
      dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
      dom[i,'SYEAR'] = ii[i]
      }
dom$Model = 'BD'
dmme=dom

                                   
f = emmeans::emmeans(l38b,~fYear+time,offset=0,type='response',data=subset(aT,LFA==38),at=list(time=1:234))
dmm = as.data.frame(summary(f))
dmm = merge(dmm,ind)
usm = sum(dmm$SD_LOG_ID)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){
  z = subset(dmm,Year==ii[i])
  usm = sum(z$SD_LOG_ID)
  dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[i,'SYEAR'] = ii[i]
}
dom$Model = 'BxD'
fmm=dom


g = emmeans::emmeans(l38c,~fYear+time,offset=0,type='response',data=subset(aT,LFA==38),at=list(time=1:234))
dmm = as.data.frame(summary(g))
dmm = merge(dmm,ind)
usm = sum(dmm$SD_LOG_ID)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){
  z = subset(dmm,Year==ii[i])
  usm = sum(z$SD_LOG_ID)
  dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[i,'SYEAR'] = ii[i]
}
dom$Model = 'BxDRL'
gmm=dom

h = emmeans::emmeans(l38d,~fYear+time,offset=0,type='response',data=subset(aT,LFA==38),at=list(time=1:234))
dmm = as.data.frame(summary(h))
dmm = merge(dmm,ind)
usm = sum(dmm$SD_LOG_ID)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){
  z = subset(dmm,Year==ii[i])
  usm = sum(z$SD_LOG_ID)
  dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[i,'SYEAR'] = ii[i]
}
dom$Model = 'BxDRV'
hmm=dom




vb = readRDS(file='unBIASED_CPUE.rds')
unb = subset(vb[[1]],LFA==38,select=c(unBCPUE, unBVar,SYEAR))
unb$unBVar = sqrt(unb$unBVar)
names(unb)=c('wemm','wse','SYEAR')
unb$Model='unbiased'
oo = do.call(rbind,list(unb,bmm,dmme,fmm,gmm,hmm))
write.csv(oo,'LFA38MultipleModelsCPUE.csv')

#marginal mean, internally consistent

ggplot(oo, aes(x = SYEAR, y = wemm ,colour=Model,fill=Model)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = wemm-wse, ymax = wemm+wse), width=0.2)+
  scale_color_viridis_d( option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)


# 
# 
# ##
# packageurl = "https://cran.r-project.org/src/contrib/Archive/ggeffects/ggeffects_1.3.4.tar.gz"
# 
# install.packages(packageurl, repos=NULL, type="source")
# 
# # #turning conditional into marginal weighting by prop of total catch
#  pr2 = ggpredict(l38a,terms=c('fYear','time[1:235]'),condition = c(leffort=0), back_transform = F)
#  pr2ta = data.frame(fYear=pr2$x,predicted=pr2$predicted,Model='BD',se=pr2$std.error,time=pr2$group)
# # 
# # 
# # ##marginal weighted by effort
# # 
# ind = aggregate(NUM_OF_TRAPS~time+SYEAR,data=subset(aT,LFA ==38),FUN=sum)
# ind1 = aggregate(NUM_OF_TRAPS~SYEAR,data=subset(aT,LFA ==38),FUN=sum)
# names(ind1)[2] = 'SumTraps'
# ind = merge(ind,ind1)
# ind$prop = ind$NUM_OF_TRAPS/ind$SumTraps
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
# pr4 = ggpredict(l38c,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
# pr4ta = data.frame(fYear=pr4$x,predicted=pr4$predicted,Model='BxDRL',se=pr4$std.error,time=pr4$group)
# pp4r = merge(pr4ta,ind[,c('fYear','time','prop')])
# pp4r$wpre = pp4r$predicted*pp4r$prop
# pp4r$wvar = (pp4r$se)^2*pp4r$prop
# p4 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp4r,FUN=sum)
# p4$wse = sqrt(p4$wvar)
# 
# pr5 = ggpredict(l38d,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
# pr5ta = data.frame(fYear=pr5$x,predicted=pr5$predicted,Model='BxDRV',se=pr5$std.error,time=pr5$group)
# pp5r = merge(pr5ta,ind[,c('fYear','time','prop')])
# pp5r$wpre = pp5r$predicted*pp5r$prop
# pp5r$wvar = (pp5r$se)^2*pp5r$prop
# p5 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp5r,FUN=sum)
# p5$wse = sqrt(p5$wvar)
# 
# 
# 
# p1 = pr1t
# names(p1)[c(1,2,6)] =c('fYear','wpre','wse')
# 
# pS = rbind(p1[,c(1,2,3,6)],p2[,c(1,2,3,5)],p3[,c(1,2,3,5)],p4[,c(1,2,3,5)],p5[,c(1,2,3,5)])
# pS$Year = as.numeric(as.character(pS$fYear))
# 
# ggplot(pS, aes(x = Year, y = wpre ,colour=Model,fill=Model)) +
#   geom_line(linewidth=1) +
#   geom_errorbar(aes(ymin = wpre-wse, ymax = wpre+wse), width=0.2)+
#   #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
#   xlab('Year')+
#   ylab('Standardized CPUE')+
#   theme_test(base_size = 14)
# 
# 
# 
# ####cross valiation
# library(tidyverse)
# library(modelr)
# library(mgcv)
# 
# #100 samples
# cv_df = crossv_mc(subset(aT,LFA==38), 50) 
# 
# #names samples
# cv_df =
#   cv_df |> 
#   mutate(
#     train = map(train, as_tibble),
#     test = map(test, as_tibble))
# 
# # run
# cv_res = 
#   cv_df |> 
#   mutate(
#     sim_mod  = map(train, \(df) mgcv::gam(WEIGHT_KG~fYear+offset(leffort),data=df,family = Gamma(link='log'))),
#     mod2  = map(train, \(df) mgcv::gam(WEIGHT_KG~fYear+s(time)+offset(leffort),data=df,family = Gamma(link='log'))), 
#     mod3  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+offset(leffort),data=df,family = Gamma(link='log'))),
#     mod4  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rLID,bs='re')+offset(leffort),data=df,family = Gamma(link='log'))), 
#     mod5  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rV,bs='re')+offset(leffort),data=df,family = Gamma(link='log')))) |>
#    mutate(
#     rmse_M1 = map2_dbl(sim_mod, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod2, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod3, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod4, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod5, test, \(mod, df) rmse(model = mod, data = df)))
# 
# #plot cv
# cv_res |> 
#   select(starts_with("rmse")) |> 
#   pivot_longer(
#     everything(),
#     names_to = "model", 
#     values_to = "rmse",
#     names_prefix = "rmse_") |> 
#   mutate(model = fct_inorder(model)) |> 
#   ggplot(aes(x = model, y = rmse)) + geom_violin()
# 
# ############################
# #spatial
# 
# aa = subset(aT,SYEAR>2005 & SYEAR<2023)
# xx = split(aa,f=aa$SYEAR)
# m=0
# ou = list()
# for(i in 1:length(xx)){
#      w=   xx[[i]]
#      y = unique(w$SYEAR)
#      ww = unique(w$GRID_NUM)
#     for(j in 1:length(ww)){
#       m=m+1
#       www = subset(w,GRID_NUM==ww[j])
#       g = unique(www$GRID_NUM)
#       uL = length(unique(www$LICENCE_ID))
#       uLFA = length(unique(www$LFA))
#       #grArea = as.numeric(sum(st_area(subset(gr,grid==ww[j]))))#these are wrong do no deal with holes right
#       wa = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~GRID_NUM+DATE_FISHED+temp,data=www,FUN=sum)
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$CPUE = wa$WEIGHT_KG / wa$NUM_OF_TRAPS
#       wa$CUME = cumsum(wa$NUM_OF_TRAPS)
#       wa$CUML = cumsum(wa$WEIGHT_KG)
#       wa$rT = round(wa$temp,2)
#       wa = merge(wa,preT,by.x='rT',by.y='temp')
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$SCPUE = wa$pred*wa$CPUE
#       wa$cpue=wa$SCPUE
#       wa$effort = wa$NUM_OF_TRAPS
#       wa$catch = wa$WEIGHT_KG
#       wa$weight = wa$NUM_OF_TRAPS
#       le = deluryLeslie(y=wa,estimate='leslie',method='robust',plot=F,weight = wa$weight)
#       ou[[m]] = c(y,g,uL,uLFA,le$q,le$N0,max(wa$CUML),max(wa$CUME),nrow(www))
#     }
#   
# }
# out = as.data.frame(do.call(rbind,ou))
# names(out) = c('Year','Grid','NLics','NLFAs','q','No','Land','Effort','Trips')
# out$CPUE = out$Land/out$Effort
# 
# ################ grouping grids
# ou = list()
# m=0
# LFAgrid <- read.csv(file.path( project.datadirectory("bio.lobster"), "data","maps","GridPolys.csv"))
# grL = (LFAgrid)
# grL$area=grL$PID
# attr(grL,'projection') <- 'LL'
# 
# #grs is the grids
# #making the list of connections to grids
# require(sp)
# require(spdep)
# 
# 
#   aa = subset(aT,SYEAR>2005 & SYEAR<2023)
#   xx = split(aa,f=aa$SYEAR)
#   m=0
#   ou = list()
#   for(i in 1:length(xx)){
#     w=   xx[[i]]
#     y = unique(w$SYEAR)
#     ww = unique(w$GRID_NUM)
#     grL1 = subset(grL,SID %in% ww)
#     g1 = split(grL1,grL1$SID)
#     nm = c()
#     gp = list()
#     for(k in 1:length(g1)){
#       gp[[k]] = Polygons(list(Polygon(g1[[k]][,c('X','Y')])),unique(g1[[k]]$SID))
#     }
#     gpp = SpatialPolygons(gp,proj4string=CRS("+proj=longlat +datum=WGS84"))
#     gpnb = poly2nb(gpp,row.names=names(gpp))
#     names(gpnb)=names(gpp)
#     
#     for(j in 1:length(gpnb)){
#       gr = names(gpnb)[gpnb[[j]]]
#       ss = ifelse(length(gr)>=3,3,length(gr))
#       gr = sample(gr,ss)
#       m=m+1
#       www = subset(w,GRID_NUM==gr)
#       g = length(unique(www$GRID_NUM))
#       gn = paste(gr,collapse="_")
#       uL = length(unique(www$LICENCE_ID))
#       uLFA = length(unique(www$LFA))
#       bgf = subset(grL,SID %in%gr)
#       attr(bgf,'projection') <-"LL"
#       grArea = sum(calcArea(convUL(bgf))$area)
#       wa = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~GRID_NUM+DATE_FISHED+temp,data=www,FUN=sum)
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$CPUE = wa$WEIGHT_KG / wa$NUM_OF_TRAPS
#       wa$CUME = cumsum(wa$NUM_OF_TRAPS)
#       wa$CUML = cumsum(wa$WEIGHT_KG)
#       wa$rT = round(wa$temp,2)
#       wa = merge(wa,preT,by.x='rT',by.y='temp')
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$SCPUE = wa$pred*wa$CPUE
#       wa$cpue=wa$SCPUE
#       wa$effort = wa$NUM_OF_TRAPS
#       wa$catch = wa$WEIGHT_KG
#       wa$weight = wa$NUM_OF_TRAPS
#       le = deluryLeslie(y=wa,estimate='leslie',method='robust',plot=F,weight = wa$weight)
#       ou[[m]] = c(y,gn,g,uL,uLFA,le$q,le$N0,max(wa$CUML),max(wa$CUME),nrow(www))
#     }
#     
#   }
#   outa = as.data.frame(do.call(rbind,ou))
#   names(outa) = c('Year','Grid','NGrids','NLics','NLFAs','q','No','Land','Effort','Trips')
#   outa$CPUE = outa$Land/outa$Effort
#   outa = toNums(outa,c(1,3:10))
#   outas = subset(outa,No>0)
#   outas$Rexp = outas$Land/outas$No
#   ggplot(outas,aes(x=Land,y=No))+geom_point()+facet_wrap(~Year,scales='free')
#   ggplot(outas,aes(x=Year,y=Rexp))
#   
#   
# head(outa)  
