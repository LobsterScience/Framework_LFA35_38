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
k = mean(ind38$est[order(ind38$est,decreasing = T)][1:3])/1000
g8 = ggplot(subset(ind38,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.4,colour='red',linewidth=1.1)+ labs(title='LFA 38')+geom_hline(yintercept = br,colour='blue',linewidth=1.1)

(ind38$est[which(ind38==2023)]/1000-br) / br*100

gridExtra::grid.arrange(g5,g6,ncol=2,g8, layout_matrix=rbind(c(1,2),c(3,4)))



g = lobster.db('landings_by_vessel')
g$year = g$SYEAR
g$LFA = ifelse(g$LFA=='38B',38,g$LFA)
glo = aggregate(MT~LFA+year,data=g,FUN=function(x) quantile(x,0.25)); names(glo)[3]='lCI'
gup = aggregate(MT~LFA+year,data=g,FUN=function(x) quantile(x,0.75)); names(gup)[3]='uCI'
gme = aggregate(MT~LFA+year,data=g,FUN=mean); names(gme)[3]='MEAN_T'

g = merge(merge(glo,gup),gme)


##lfa35
g5 = subset(g,LFA==35)
nd35$oyr = ind35$year
ind35$year=ind35$oyr+1

g5 = merge(g5,ind35)
g5$est1000 = g5$estB/1000
mL = mean(subset(g,LFA==35 &year %in% 1990:1995)$MEAN_T)

g3 <- ggplot(g5, aes(est1000, MEAN_T,ymin=lCI,ymax=uCI,xmin=lwrB/1000,xmax=uprB/1000)) + geom_point(size=2)+geom_errorbar(data=g5, aes(est1000, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test(base_size = 14)+labs(x='Survey Biomass',y='Landings per Vessel')

###hyperbolic
g5$sE = g5$est1000
g5 = subset(g5,year<2024)
model <- brm(MEAN_T~I(1/(sE)),data=g5,family=gaussian(link='inverse'),
             prior = c(
               set_prior("normal(15, 6)", class = "b"),
               set_prior("normal(0, 10)", class = "Intercept")
             ),
             chains = 4,
             iter = 2000, cores=4
)

plot(model, pars=c('b'))

xs = data.frame(sE=seq(0,max(g5$est1000),length.out=20))
xp = as.data.frame(model) #posterior parameters

###using the 
jnk = list()
for(i in 1:nrow(xp)){
  xs$Pred = 1/(xp[i,'b_Intercept']+xp[i,'b_I1DsE']/xs$sE)
  xs$i = i
  jnk[[i]] = xs
}

###LRP from base model
ests = seq(0,300,by=.1)
#distribution of LRP  
LRP=c()
for(i in 1 :nrow(xp)) {
  y = 1/(xp[i,'b_Intercept']+xp[i,'b_I1DsE']/ests)
  LRP = c(LRP, ests[which.min(abs(y-mL))])
}     


kk = bind_rows(jnk)
bk = aggregate(Pred~sE,data=kk, FUN=function(x) quantile(x,c( 0.025,0.975)))

#plot each line from the posterior params
ggplot(g5, aes(est1000, MEAN_T)) + geom_point(size=2)+geom_errorbar(data=g5, aes(est1000, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
  geom_line(data=kk,aes(x=sE,y=Pred,group=i),colour='blue',alpha=.01) +labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)


#Function for obtaining biomass estimates from teh distributions
gsa <- function(row) {
  mm = findMoments(as.numeric(row['lwrB'])/1000,as.numeric(row['uprB'])/1000,dist='lnorm',plot=F)
  rlnorm(1, mean = mm$xbar, sd = sqrt(mm$sig2))
}

##re run the model with each draw from teh distribution to generate EIV
junk=list()
iter=100
for(i in 1:iter){
  g5$sE <- apply(g5, 1, gsa)
  um = update(model, newdata=g5)
  junk[[i]]  = as.data.frame(um)
}

dd = bind_rows(junk)
LFA35Params = ggplot(dd,aes(b_Intercept,b_I1DsE))+
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  scale_x_continuous(limits=c(0.017,.058))+
  theme_test(base_size=14)+ labs(x='Intercept',y='Slope')+theme(legend.position='none')

#generate predictions from each pair of slope and response
jnk = list()
for(i in 1:nrow(dd)){
  xs$Pred = 1/(dd[i,'b_Intercept']+dd[i,'b_I1DsE']/xs$sE)
  xs$i = i
  jnk[[i]] = xs
}

tos = bind_rows(jnk)

####lrp

LRP_D <- data.frame(est1000=seq(0,300,by=.2))

invG = function(x,b0,b1){
  1/(b0+(b1*1/x))
}
ests = seq(-400,3000,by=.1)
#distribution of LRP  
LRP_jp=c()
LRP_buff = c()

for(i in 1 :nrow(dd)) {
  y = 1/(dd[i,'b_Intercept']+dd[i,'b_I1DsE']/ests)
  LRP_jp = c(LRP_jp, ests[which.min(abs(y-mL))])
  LRP_buff = c(LRP_buff, ests[which.min(abs(y-mL*1.25))])
}     

HistLRPS_35 = ggplot(as.data.frame(LRP_jp),aes(LRP_jp))+geom_density()+geom_density(data=data.frame(LRP_buff),aes(LRP_buff),colour='red')+
  geom_density(data=data.frame(LRP),aes(LRP),colour='blue')+geom_vline(xintercept = c(mean(LRP_jp),mean(LRP_buff),mean(LRP)),colour=c('black','red','blue'),linewidth=1, linetype='dotted')+
  labs(x='LRP',y='Frequency') +theme_test(base_size = 14)


b = aggregate(Pred~sE,data=tos, FUN=function(x) quantile(x,c( 0.025,0.975)))

gridExtra::grid.arrange(LFA35Params,HistLRPS_35)

#EIV
e1 = ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbar(data=g5, aes(est1000,y=MEAN_T, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
  geom_ribbon(data=b, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
  labs(x='Survey Commercial Biomass',y='Landings per Vessel')+theme_test(base_size = 14)+
  geom_vline(xintercept=LRP_jp,colour='grey',alpha=.01)+geom_hline(yintercept=mL,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP_jp),colour='black')+
  scale_x_continuous(limits= c(0,4500))

#Buff
e2 = ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbar(data=g5, aes(est1000,y=MEAN_T, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
  geom_ribbon(data=b, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
  labs(x='Survey Commercial Biomass',y='Landings per Vessel')+theme_test(base_size = 14)+
  geom_vline(xintercept=LRP_buff,colour='red',alpha=.005)+geom_hline(yintercept=mL*1.25,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP_buff),colour='black')+scale_x_continuous(limits= c(0,4500))

#base
e3 = ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbar(data=g5, aes(est1000,y=MEAN_T, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
  geom_ribbon(data=bk, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
  labs(x='Survey Commercial Biomass',y='Landings per Vessel')+theme_test(base_size = 14)+
  geom_vline(xintercept=LRP,colour='blue',alpha=.005)+geom_hline(yintercept=mL,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP),colour='black')+scale_x_continuous(limits= c(0,4500))

gridExtra::grid.arrange(e1,e2,e3,ncol=2,layout_matrix=rbind(c(1,2),c(3,4)))

mean(LRP); quantile(LRP, c(0.025,0.975))      
mean(LRP_jp); quantile(LRP_jp, c(0.025,0.975))     
mean(LRP_buff); quantile(LRP_buff, c(0.025,0.975))  
###################
#SP model

#LFA 38
require(bio.lobster)
require(devtools)
require(lubridate)
la()

fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)
io = readRDS('IndicesFromFullComboModelJan282000+.rds')
ind35 = io[[1]]



##converting number to biomass using fall length frequencies from all fall surveys in bay of fundy

yp = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=F,species=2550,biomass = F,extend_ts = F)
yp2 = subset(yp, month(SET_DATE)>8)

yp2 = aggregate(SA_CORRECTED_PRORATED_N~FISH_LENGTH,data=yp2,FUN=mean)
names(yp2)[2]='meanden'
yp2$Prop = yp2$meanden/sum(yp2$meanden)
yp2$Weight = lobLW(CL=yp2$FISH_LENGTH,sex = 1)
yp2s = subset(yp2,Prop>0)
mnW = sum(yp2s$Prop*yp2s$Weight)

ind35$estB = ind35$est*mnW/1000
ind35$lwrB = ind35$lwr*mnW/1000
ind35$uprB = ind35$upr*mnW/1000


land = lobster.db('seasonal.landings')
land$year = substr(land$SYEAR,6,9)
ind35$year = ind35$year+1
Da = merge(ind35,land[,c('year','LFA35')])

input1 <-  list()
input1$N <- nrow(Da)
input1$I <- Da$estB/1000
input1$C <- Da$LFA35



#Add in Priors
Ku2<-max(input1$I)

#rs <- find.moments(lo=0.1,up=0.8,l.perc=0.05,u.perc=0.95,dist='norm')
k2s <- findMoments(lo=1/Ku2,up=1/(Ku2*5),l.perc=0.05,u.perc=0.95,dist='lnorm')
K2s <<- findMoments(lo=Ku2,up=(Ku2*5),l.perc=0.05,u.perc=0.95,dist='lnorm')

input1$r.a<-9 ; input1$r.b <- 1/0.993347893 #informative prior from mcallister work
input1$k2.a<-k2s[[1]] ; input1$k2.b <- k2s[[2]]
input1$q.a<-0.01 ; input1$q.b <- 1
#Process = a b observation=c d
input1$a0<-0.000001;  input1$b0<-5; 
input1$c0<-0.000001; input1$d0<-5;



spBUGS3 <- function(input=input1,inits=inits,n = 10000, burn = 2000, thin = 10, debug = F, wd=getwd()){
  require(R2OpenBUGS)
  
  #	Initial values for WinBUGS one for each chain
  inits <- list(list(P=rep(0.75,times=input$N),
                     r2=0.8,
                     k2=0.003,
                     q2=0.4,
                     P.p = matrix(0.1,ncol=4,nrow=6)
  ),
  list(P=rep(0.8,times=input$N),
       r2=0.9,
       k2=0.005,
       q2=0.4,
       P.p = matrix(0.1,ncol=4,nrow=6)
       
  ))
  
  # Parameters to be returned 
  parameters <- c('Tau21','Sigma21','r2','K2','q2','B','Imean','P.res','BMSP2','BMSPc2','MSP2','MSPc2','FMSP2','FMSPc2')
  
  
  ## Call to WinBUGS ##
  tmp <- bugs(data = input, inits, parameters.to.save = parameters, model.file = file.path("C:",'Users',"CookA","Documents","git","bio.lobster","inst","SPmodel","BoFSp", "sp2000p.txt"), 
              n.chains = 2, n.iter = n, n.burnin = burn, n.thin = thin,  debug = debug)
  return(tmp)
} 

a<- spBUGS3(debug=F)


#########plots
op<-par(mfrow=c(3,2),mar=c(4,4,3,1))
# r
plot(density(rnorm(10000,input1$r.a,sqrt(input1$r.b))),xlab='r',xlim=c(0,1.5),lty=2,main='')
#lines(density(a$sims.list$r1),col='red',lty=1)
lines(density(a$sims.list$r2),col='blue',lty=1)
quantile(a$sims.list$r,c(0.025,0.5,0.975))

#K
plot(density(rlnorm(10000,K2s[[1]],sqrt(K2s[[2]]))),xlab='K',lty=2,main='')
#lines(density(rlnorm(10000,K2s[[1]],sqrt(K2s[[2]]))),lty=2)
#lines(density(a$sims.list$K1),col='red',lty=1)
lines(density(a$sims.list$K2),col='blue',lty=1)
quantile(a$sims.list$K2,c(0.025,0.5,0.975))
#q
plot(1,1,type='n',xlim=c(input1$q.a/1.1,input1$q.b*1.1),xlab='q',main='',ylab='Density')
polygon(x=c(input1$q.a,input1$q.a,input1$q.b,input1$q.b),y=c(0,dunif(c(input1$q.a,input1$q.b),input1$q.a,input1$q.b),0))
#lines(density(a$sims.list$q1),col='red',lty=1)
lines(density(a$sims.list$q2),col='blue',lty=1)
quantile(a$sims.list$q2,c(0.025,0.5,0.975))
#process errors
plot(1,1,type='n',xlim=c(0,3),xlab='Process Variance',main='',ylab='Density',ylim=c(0,0.6))
polygon(x=c(0.25^2,0.25^2,5,5),y=c(0,0.4,0.4,0),border='black',lty=2)
#lines(density(sqrt(a$sims.list$Sigma2)),col='red',lty=1,lwd=1)
lines(density(sqrt(a$sims.list$Sigma21)),col='blue',lty=1,lwd=1)

plot(1,1,type='n',xlim=c(0,1),xlab='Observation Variance',main='',ylab='Density',ylim=c(0,0.6))
polygon(x=c(0.25^2,0.25^2,1,1),y=c(0,0.4,0.4,0),border='black',lty=2)
#lines(density(sqrt(a$sims.list$Tau2)),col='red',lty=1,lwd=1)
lines(density(sqrt(a$sims.list$Tau21)),col='blue',lty=1,lwd=1)


plot(Da$year,a$median$P.res,type='h',ylab='Residuals',xlab='Year')
abline(h=0,col='red')

graphics.off()

##############################

dat1<-data.frame(years=Da$year,Biomass=0,Landings=0)
dat1$Biomass<-a$mean$B
dat1$Landings<-input1$C

dat1$ASP <- 0
for(i in 1:(nrow(dat1)-1)) {
  dat1[i,'ASP'] <- dat1[i+1,'Biomass'] - dat1[i,'Biomass'] +dat1[i,'Landings'] 
}

#MSY V B plots
plotArrows <- function (Y,X) {
  #//trace arrows around x,y points from start to end
  xl <- c(min(X),max(X))
  yl <- c(min(Y),max(Y))
  if(Y[length(Y)]==0) {Y<-Y[-length(Y)]}
  if(length(Y) !=length(X)) {X<-X[-length(X)]}
  n <- length(Y)
  cl <- rainbow(n)
  plot(1,1,ylim=yl,xlim=xl,type='n',xlab='Biomass',ylab = 'Annual Surplus Production')
  
  
  for(i in 1:n) {
    arrows(x0=X[i],x1=X[i+1],y0=Y[i],y1=Y[i+1],angle=45,length=0.1,lwd=2.5,col=cl[i])
  }
  points(X,Y,pch=16,cex=1,col='black')
  points(X[1],Y[1],pch=17,cex=2,col='black')
  points(X[n],Y[n],pch=17,cex=2,col='black')
}

plotArrows(dat1$ASP,dat1$Biomass)



#MSY plots with 40 and 80 BMSY
fit <- lower <- upper <- numeric()
N<-(Da$year[1:nrow(Da)])
n<-length(N)
for(i in 1:length(N)) {
  fit[i]<-median(a$sims.list$B[,i])
  lower[i]<-quantile(a$sims.list$B[,i],0.025)
  upper[i]<-quantile(a$sims.list$B[,i],0.975)
}
yl <- c(min(c(fit,lower,upper)),max(c(fit,lower,upper)))
plot(1,1, type='n',ylab='Biomass',xlab='Year',ylim=yl,xlim=c(min(N),max(N)))
#polygon(c(1992,2011,2011,1992),c(0.8*a$median$BMSP,0.8*a$median$BMSP,0.4*a$median$BMSP,0.4*a$median$BMSP),density=10,col='orange')
lines(N,fit)
lines(N,lower,lty=2)
lines(N,upper,lty=2)

polygon(c(1976,2024,2024,1976),c(0.4*a$median$BMSP2,0.4*a$median$BMSP2,0.4*a$median$BMSP2,0.4*a$median$BMSP2),density=10,col='blue')
#legend('topright',fill=c('orange','blue'),c('Det','Sto'),bty='n',cex=0.8,density=20)


par(mfrow=c(2,2))

plot(density(a$sims.list$BMSP1),main='',xlab='BMSP')
lines(density(a$sims.list$BMSP2),main='',xlab='BMSP')

plot(density(a$sims.list$MSP*10^af),main='',xlab='MSP')
lines(density(a$sims.list$MSPc*10^af),col='red')

plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
legend(0.2,0.7,c('Det','Sto'),lty=c(1,1),col=c('black','red'),lwd=2,bty='n')

dev.off()
}
#---------------------------

fit <- lower <- upper <- numeric()
N<-(Da$year[1:nrow(Da)])
n<-length(N)
for(i in 1:length(N)) {
  fit[i]<-median(a$sims.list$B[,i])
  lower[i]<-quantile(a$sims.list$B[,i],0.25)
  upper[i]<-quantile(a$sims.list$B[,i],0.75)
}
Land=input1$C
dd = data.frame(fit,lower,upper,N,Land)

ggplot(dd,aes(x=N,y=fit,ymin=lower,ymax=upper))+geom_point()+geom_ribbon(alpha=.25)+geom_path()+theme_test(base_size = 14)+labs(x='Year',y='Biomass')+geom_hline(yintercept=a$median$BMSP2*0.4, colour='blue',linewidth=1.2)
quantile(a$sims.list$BMSP2,c(0.025,0.5,0.975))*.4

dd$lwrF = U2F(dd$Land/dd$lower)
dd$uprF = U2F(dd$Land/dd$upper)
dd$F = U2F(dd$Land/dd$fit)


ggplot(dd,aes(x=N,y=F,ymin=lwrF,ymax=uprF))+geom_point()+geom_ribbon(alpha=.25)+geom_path()+theme_test(base_size = 14)+labs(x='Year',y='Fishing Mortality')+geom_hline(yintercept=a$median$FMSP2, colour='blue',linewidth=1.2)

