### reference points and relative fishing mortality from final model


library(splines)
library(MASS) ##for glm.nb function
library(sdmTMB)
library(sdmTMBextra)
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
require(brms)

crs_utm20 <- 32620

fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)
io = readRDS('IndicesFromFullComboModelJan282000+.rds')
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


###BMSY plots + Brecover

#LFA 35
k = mean(ind35$est[order(ind35$est,decreasing = T)][1:3])/1000
br = min(ind35$est)/1000
g5 = ggplot(subset(ind35,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.4,colour='red',linewidth=1.1) + labs(title='LFA 35') +geom_hline(yintercept = br,colour='blue',linewidth=1.1)

#36
k = mean(ind36$est[order(ind36$est,decreasing = T)][1:3])/1000
br = min(ind36$est)/1000
(ind36$est[which(ind36==2023)]/1000-br) / br*100
g6 = ggplot(subset(ind36,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.4,colour='red',linewidth=1.1)+ labs(title='LFA 36')+geom_hline(yintercept = br,colour='blue',linewidth=1.1)

#38
br = min(ind38$est)/1000
k = mean(ind38$est[order(ind38$est,decreasing = T)][1:3])/1000*.4
g8 = ggplot(subset(ind38,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.4,colour='red',linewidth=1.1)+ labs(title='LFA 38')+geom_hline(yintercept = br,colour='blue',linewidth=1.1)

(ind38$est[which(ind38==2023)]/1000-br) / br*100

gridExtra::grid.arrange(g5,g6,ncol=2,g8, layout_matrix=rbind(c(1,2),c(3,4)))



g = lobster.db('landings_by_vessel')
g$year = g$SYEAR

g8 = aggregate(MT~year,data=subset(g,LFA %in% c('38B')),FUN=sum)
gg = aggregate(MT~year+LFA,data=subset(g,LFA %in% c(38,'38B')),FUN=sum)

gl = aggregate(VR_NUMBER~year+LFA,data=subset(g,LFA %in% c(38,'38B')),FUN=function(x) length(unique(x)))
glg = merge(gl,gg)
glg$MEAN_T = glg$MT/glg$VR_NUMBER
glg = aggregate(MEAN_T~year,data=glg,FUN=sum)



##lfa38
ind38$oyr = ind38$year
ind38$year=ind38$oyr+1

g5 = merge(glg,ind38)
g5$est1000 = g5$estB/1000
mL = mean(subset(glg,year %in% 1990:1995)$MEAN_T)
g5 = subset(g5,year<2024)
g3 <- ggplot(g5, aes(est1000, MEAN_T,xmin=lwrB/1000,xmax=uprB/1000)) + geom_point(size=2)+geom_errorbarh(data=g5, aes(xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test(base_size = 14)+labs(x='Survey Biomass',y='Landings per Vessel')

###linear
g5$sE = g5$est1000
g5 = subset(g5,year<2024)
model <- brm(MEAN_T~sE,data=g5,family=gaussian,
             prior = c(
               set_prior("normal(0,5)", class = "b"),
               set_prior("normal(0, 5)", class = "Intercept")
             ),
             chains = 4,
             iter = 2000,cores=4
             
)

###dist of params
ggplot(as.data.frame(model),aes(b_Intercept,b_sE))+
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()
mod = as.data.frame(model)
###lines
g3 <- ggplot(g5, aes(est1000, MEAN_T,xmin=lwrB/1000,xmax=uprB/1000)) + geom_point(size=2)+geom_errorbarh(data=g5, aes(xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()


xs = data.frame(sE=seq(0,max(g5$est1000),length.out=20))
xp = as.data.frame(model)

jnk = list()
for(i in 1:nrow(xp)){
  xs$Pred = (xp[i,'b_Intercept']+xp[i,'b_sE']*xs$sE)
  xs$i = i
  jnk[[i]] = xs
}

kk = bind_rows(jnk)
bk = aggregate(Pred~sE,data=kk, FUN=function(x) quantile(x,c( 0.025,0.975)))

###LRP from base model
ests = seq(0,300,by=.1)
#distribution of LRP  
LRP=c()
for(i in 1 :nrow(xp)) {
  y = (xp[i,'b_Intercept']+xp[i,'b_sE']*ests)
  LRP = c(LRP, ests[which.min(abs(y-mL))])
}     



i=0
gsa <- function(row) {
  mm = findMoments(as.numeric(row['lwrB'])/1000,as.numeric(row['uprB'])/1000,dist='lnorm',plot=F)
  rlnorm(1, mean = mm$xbar, sd = sqrt(mm$sig2))
}

junk=list()
iter=100
for(i in 1:iter){
  g5$sE <- apply(g5, 1, gsa)
  um = update(model, newdata=g5)
  junk[[i]]  = as.data.frame(um)
}

dd = bind_rows(junk)
LFA35Params =ggplot(dd,aes(b_Intercept,b_sE))+
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_test(base_size=14)+ labs(x='Intercept',y='Slope')


jnk = list()
for(i in 1:nrow(dd)){
  xs$Pred = (dd[i,'b_Intercept']+dd[i,'b_sE']*xs$sE)
  xs$i = i
  jnk[[i]] = xs
}

tos = bind_rows(jnk)



####lrp

LRP_D <- data.frame(est1000=seq(0,300,by=.2))

ests = seq(-1000,100,by=1)
#distribution of LRP  
LRP_jp=c()
LRP_buff = c()
#ddthin = dd[seq(1,nrow(dd),by=50),]
for(i in 1 :nrow(dd)) {
  y = dd[i,'b_Intercept']+dd[i,'b_sE']*ests
  LRP_jp = c(LRP_jp, ests[which.min(abs(y-mL))])
  LRP_buff = c(LRP_buff, ests[which.min(abs(y-mL*1.25))])
}     





HistLRPS_35 = ggplot(as.data.frame(LRP_jp),aes(LRP_jp))+geom_density()+geom_density(data=data.frame(LRP_buff),aes(LRP_buff),colour='red')+
  geom_density(data=data.frame(LRP),aes(LRP),colour='blue')+geom_vline(xintercept = c(mean(LRP_jp),mean(LRP_buff),mean(LRP)),colour=c('black','red','blue'),linewidth=1, linetype='dotted')+
  labs(x='LRP',y='Frequency') +theme_test(base_size = 14)

gridExtra::grid.arrange(LFA35Params,HistLRPS_35)

mean(LRP); quantile(LRP, c(0.025,0.975))      
mean(LRP_jp); quantile(LRP_jp, c(0.025,0.975))     
mean(LRP_buff); quantile(LRP_buff, c(0.025,0.975))  


b = aggregate(Pred~sE,data=tos, FUN=function(x) quantile(x,c( 0.025,0.975)))


#EIV
e1 = ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
  geom_ribbon(data=b, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
  labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)+
  geom_vline(xintercept=LRP_jp,colour='grey',alpha=.01)+geom_hline(yintercept=mL,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP_jp),colour='black')

#Buff
e2 = ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
  geom_ribbon(data=b, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
  labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)+
  geom_vline(xintercept=LRP_buff,colour='red',alpha=.005)+geom_hline(yintercept=mL*1.25,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP_buff),colour='black')

#base
e3 = ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
  geom_ribbon(data=bk, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
  labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)+
  geom_vline(xintercept=LRP,colour='blue',alpha=.005)+geom_hline(yintercept=mL,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP),colour='black')


gridExtra::grid.arrange(e1,e2,e3,ncol=2,layout_matrix=rbind(c(1,2),c(3,4)))
