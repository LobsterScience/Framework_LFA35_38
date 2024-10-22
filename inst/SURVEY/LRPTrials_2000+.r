### reference points and relative fishing mortality from final model


library(splines)
library(MASS) ##for glm.nb function

#remotes::install_github("TobieSurette/gulf.spatial", dependencies = TRUE)

library(sdmTMB)
library(sdmTMBextra)

##these are used to make the mesh
library(dplyr)
library(ggplot2)
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
crs_utm20 <- 32620

fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)
io = readRDS('IndicesFromFullComboModelOct92000+.rds')
ind35 = io[[1]]
ind36 = io[[2]]
ind38 = io[[3]]

ggplot(subset(ind35,year>=2000 &year<2024),aes(year,est,ymin=lwr,ymax=upr))+geom_point()+geom_ribbon(alpha=.25)+theme_test()+geom_path()

ggplot(subset(ind36,year>=2000 &year<2024),aes(year,est,ymin=lwr,ymax=upr))+geom_point()+geom_ribbon(alpha=.25)+theme_test()+geom_path()

ggplot(subset(ind38,year>=2000 &year<2024),aes(year,est,ymin=lwr,ymax=upr))+geom_point()+geom_ribbon(alpha=.25)+theme_test()+geom_path()


##converting number to biomass using fall length frequencies from all fall surveys in bay of fundy

yp = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=F,species=2550,biomass = F,extend_ts = F)
yp2 = subset(yp, month(SET_DATE)>8)

yp2 = aggregate(SA_CORRECTED_PRORATED_N~FISH_LENGTH,data=yp2,FUN=mean)
names(yp2)[2]='meanden'
yp2$Prop = yp2$meanden/sum(yp2$meanden)
yp2$Weight = lobLW(CL=yp2$FISH_LENGTH,sex = 1)
yp2s = subset(yp2,Prop>0)
mnW = sum(yp2s$Prop*yp2s$Weight)
ggplot(yp2,aes(x=FISH_LENGTH,y=Prop))+geom_bar(stat='identity')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)

#commercial abundance to weight
ind35$estB = ind35$est*mnW/1000
ind35$lwrB = ind35$lwr*mnW/1000
ind35$uprB = ind35$upr*mnW/1000

ind36$estB = ind36$est*mnW/1000
ind36$lwrB = ind36$lwr*mnW/1000
ind36$uprB = ind36$upr*mnW/1000

ind38$estB = ind38$est*mnW/1000
ind38$lwrB = ind38$lwr*mnW/1000
ind38$uprB = ind38$upr*mnW/1000





g = lobster.db('seasonal.landings')
g$year= as.numeric(substring(g$SYEAR,6,9))
ind35$year=ind35$year+1
g5 = merge(g[,c('year','LFA35')],ind35)
g5$est1000 = g5$est/1000

#beverton holt with gamma errors (Jiao et al 2004 CJFAS 61)
a  <- glm(LFA35~I(1/(est1000)),data=g5,family=Gamma(link='inverse'))
b <- coef(a)

bvh <- function(a,b,S) {
  a*S / (1+b*S)
}

dn <- data.frame(est1000=c(0,g5[order(g5$est1000),'est1000']))
plot(LFA35~est1000,data=g5)
plot(dn$est1000,bvh(a=1/b[2],b=b[1]*b[2],S=dn$est1000))
legend('topleft',lty=c(1,2,3),cex=.75,c('Beverton-Holt LN','Density Independent',"Beverton-Holt Gamma"),bty='n')





ggplot(g5,aes(x=est,y=LFA35))+geom_point()+



ggplot(g5,aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')
ggplot(g5,aes(year,LFA35))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)

#tremblay 2012
ggplot(g5,aes(year,LFA35))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)+geom_hline(yintercept =292,colour='red',linewidth=2 )

#if we used Bmsy proxies
k = mean(ind35$est[order(ind35$est,decreasing = T)][1:10])/1000
ggplot(subset(ind35,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.2,colour='red',linewidth=2)




ggplot(subset(g5,year<2024 & year>=1995),aes(year,y=LFA35/(estB/1000)))+geom_point()+
  geom_line()+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality')

ind36$year=ind36$year+1
g6 = merge(g[,c('year','LFA36')],ind36)
g6$relf = g6$LFA36/g6$estB*1000
ggplot(subset(g6,year<2024 & year>1995),aes(year,y=LFA36/(estB/1000)))+geom_point()+
  geom_line()+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality')

ggplot(g6,aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')
ggplot(g6,aes(year,LFA36))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)

#tremblay 2012
ggplot(g6,aes(year,LFA36))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)+geom_hline(yintercept =266,colour='red',linewidth=2 )

#if we used Bmsy proxies
k = mean(ind36$est[order(ind36$est,decreasing = T)][1:10])/1000
ggplot(subset(ind36,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.2,colour='red',linewidth=2)




ind38$year=ind38$year+1

g8 = merge(g[,c('year','LFA38')],ind38)
g8$relf = g8$LFA38/g8$estB*1000

ggplot(g8,aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')
ggplot(g8,aes(year,LFA38))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)

#tremblay 2012
ggplot(g8,aes(year,LFA38))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)+geom_hline(yintercept =259,colour='red',linewidth=2 )

#if we used Bmsy proxies
k = mean(ind38$est[order(ind38$est,decreasing = T)][1:10])/1000
ggplot(subset(ind38,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.2,colour='red',linewidth=2)


ggplot(subset(g8,year<2024 & year>1995),aes(year,y=LFA38/(estB/1000)))+geom_point()+
  geom_line()+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality')


################################################
################cpue to survey biomass hyperstablity plots

vb = readRDS(file=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','CPUE','unBIASED_CPUE.rds'))
vb = vb[[1]]

vb$year=vb$SYEAR
vb5 = merge(subset(vb,LFA==35),g5)
vb5$estBred =  vb5$estB/1000000

ggplot(vb5,aes(x=estBred,y=unBCPUE))+
  geom_path()+
  geom_point()+
  geom_text(data=subset(vb5,year %in% c(2006,2023)),aes(label=year,x=estBred,y=unBCPUE),position=position_nudge(y=c(-.1,.1)))+
  geom_smooth(method = "nls",formula = y ~ a*x^b,se=FALSE, method.args = list(start = list(a = -1, b = .9)),   data = vb5)+theme_test(base_size=14)+labs(x='Survey Commercial Biomass',y='CPUE')

bb = nls(unBCPUE~a*estBred^b,start=(list(a=1,b=.9)),data=vb5)
bo5 = nlsBooter(bb,vb5,Qlow1 = .025,Qhigh1 = 0.975)
b05a = data.frame(b=bo5$coefboot[,2],LFA=35)

vb$year=vb$SYEAR
vb6 = merge(subset(vb,LFA==36),g6)
vb6$estBred =  vb6$estB/1000000

ggplot(vb6,aes(x=estBred,y=unBCPUE))+
  geom_path()+
  geom_point()+
  geom_text(data=subset(vb6,year %in% c(2006,2023)),aes(label=year,x=estBred,y=unBCPUE),position=position_nudge(x=c(0,-.09),y=c(-.05,-0.01)))+
  geom_smooth(method = "nls",formula = y ~ a*x^b,se=FALSE, method.args = list(start = list(a = -1, b = .9)),   data = vb6)+theme_test(base_size=14)+labs(x='Survey Commercial Biomass',y='CPUE')

b6 = nls(unBCPUE~a*estBred^b,start=(list(a=1,b=.9)),data=vb6)
bo6 = nlsBooter(b6,vb6,Qlow1 = .025,Qhigh1 = 0.975)
b06a = data.frame(b=bo6$coefboot[,2],LFA=36)



vb$year=vb$SYEAR
vb8 = merge(subset(vb,LFA==38),g8)
vb8$estBred =  vb8$estB/1000000

ggplot(vb8,aes(x=estBred,y=unBCPUE))+
  geom_path()+
  geom_point()+
  geom_text(data=subset(vb8,year %in% c(2006,2023)),aes(label=year,x=estBred,y=unBCPUE),position=position_nudge(x=c(0,.02),y=c(-.05,-0.01)))+
  geom_smooth(method = "nls",formula = y ~ a*x^b,se=FALSE, method.args = list(start = list(a = -1, b = .9)),   data = vb8)+theme_test(base_size=14)+labs(x='Survey Commercial Biomass',y='CPUE')

b8 = nls(unBCPUE~a*estBred^b,start=(list(a=1,b=.9)),data=vb8)

bo8 = nlsBooter(b8,vb8,Qlow1 = .025,Qhigh1 = 0.975)
b08a = data.frame(b=bo8$coefboot[,2],LFA=38)

pars = do.call(rbind,list(b05a,b06a,b08a))

ggplot(pars,aes(x=b))+geom_bar(stat='density')+facet_wrap(~LFA)+theme_test(base_size = 14)+labs(x='Shape Parameter',y='Density')

########
####example hyperstability plot

x = seq(0,1,by=.01)
df = data.frame(x)
ggplot(df,aes(x))+
  stat_function(fun=function(x) x^.25)+
  stat_function(fun=function(x) x^.5)+
  stat_function(fun=function(x) x^1)+
  stat_function(fun=function(x) x^2)+
  stat_function(fun=function(x) x^4)+
  theme_test(base_size = 14,axis.text.x=element_blank())+
  labs(x='Biomass',y='CPUE')+
  annotate("text",x=.25,y=.9,label='Hyperstability')+
  annotate("text",x=.8,y=.1,label='Hyperdepletion')+
  annotate("text",x=.13,y=.67,label=expression(paste(beta,"= 0.25",sep=" ")))+
  annotate("text",x=.23,y=.53,label=expression(paste(beta,"= 0.5",sep=" ")))+
  annotate("text",x=.33,y=.39,label=expression(paste(beta,"= 1",sep=" ")))+
  annotate("text",x=.43,y=.25,label=expression(paste(beta,"= 2",sep=" ")))+
  annotate("text",x=.53,y=.11,label=expression(paste(beta,"= 4",sep=" ")))
  
