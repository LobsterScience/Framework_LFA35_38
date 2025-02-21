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

fd=file.path(project.datadirectory('Assessment_LFA35_38'),'outputs','SURVEYS')
setwd(fd)
io = readRDS('IndicesFromFullComboModelFeb14_2000+.rds')
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
ggplot(yp2,aes(x=FISH_LENGTH,y=Prop))+geom_bar(stat='identity')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)

#commercial abundance to weight
ind35$estB = ind35$est*mnW/1000
ind35$lwrB = ind35$lwr*mnW/1000
ind35$uprB = ind35$upr*mnW/1000



g = lobster.db('landings_by_vessel.redo')
g$year = g$SYEAR

glo = aggregate(MT~LFA+year,data=g,FUN=function(x) quantile(x,0.25)); names(glo)[3]='lCI'
gup = aggregate(MT~LFA+year,data=g,FUN=function(x) quantile(x,0.75)); names(gup)[3]='uCI'
gme = aggregate(MT~LFA+year,data=g,FUN=mean); names(gme)[3]='MEAN_T'

g = merge(merge(glo,gup),gme)

g = subset(g,LFA==35)

gl = lobster.db('seasonal.landings')
gl$year= as.numeric(substring(gl$SYEAR,6,9))
gl = subset(gl,select=c(year,LFA35))

gg = merge(gl,g)
maxl =max(gg$LFA35)
maxv = max(gg$MEAN_T)
ggplot(subset(gg), aes(x = year)) +
  geom_bar(aes(y =LFA35), stat = "identity", width = 0.75,alpha=.25) +
  geom_point(data = subset(gg), aes(y = MEAN_T * (maxl /maxv)),size=3)+
  geom_line(data = subset(gg), aes(y = MEAN_T * (maxl /maxv)),linewidth=1.1) +
  labs(x = "Season Year", y = "Landings (t)") +
    scale_y_continuous(
    limits = c(0, 4500),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * (maxv / maxl), name = "Mean Landings Per Vessel (t)")) + 
    scale_x_continuous(breaks = seq(1990, 2024, by = 5), expand = c(0.01, 0)) + 
  guides(fill = "none", color = "none") 



#############relationship between mean landings per vessel and Commercial Biomass
ind35$oyr = ind35$year

ind35$year=ind35$oyr+1

g5 = merge(gg,ind35)
g5$est1000 = g5$estB/1000
with(g5,plot(est1000,MEAN_T))



        
##############################################################################################
################################## bayesian to explore the errors in variables ### USE THIS POST FRAMEWORK
###############################################################################################                            
                            require(brms)
                            
                            g5$sE = g5$est1000
                            g5 = subset(g5,year<2024) ### this is a fixed time period for determining LRP
                            model <- brm(MEAN_T~I(1/(sE)),data=g5,family=gaussian(link='inverse'),
                                          prior = c(
                                            set_prior("normal(15, 6)", class = "b"),
                                            set_prior("normal(0, 10)", class = "Intercept")
                                          ),
                                          chains = 4,
                                          iter = 2000, cores=4
                            )
                            
                            ###dist of params
                            ggplot(as.data.frame(model),aes(b_Intercept,b_I1DsE))+
                              geom_hex(bins = 70) +
                              scale_fill_continuous(type = "viridis") +
                              scale_x_continuous(limits=c(0.015,.04))+
                            theme_bw()
                            mod = as.data.frame(model)
         
                            
                     ##rerun based on EIV This step takes time
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
                            ggplot(dd,aes(b_Intercept,b_I1DsE))+
                              geom_hex(bins = 70) +
                              scale_fill_continuous(type = "viridis") +
                              scale_x_continuous(limits=c(0.01,.058))+
                              theme_test(base_size=14)+ labs(x='Intercept',y='Slope')
                            
                            
                            jnk = list()
                            for(i in 1:nrow(dd)){
                              xs$Pred = 1/(dd[i,'b_Intercept']+dd[i,'b_I1DsE']/xs$sE)
                              xs$i = i
                              jnk[[i]] = xs
                            }
                            
                            tos = bind_rows(jnk)
                            
                            
                            
                            ####lrp
                            
                            mL = mean(subset(g,year %in% 1990:1995)$MEAN_T)
                            LRP_D <- data.frame(est1000=seq(0,300,by=.2))
                            
                            invG = function(x,b0,b1){
                              1/(b0+(b1*1/x))
                            }
                            ests = seq(-400,3000,by=.1)
                            #distribution of LRP  
                            LRP_jp=c()
                            LRP_buff = c()
                            #ddthin = dd[seq(1,nrow(dd),by=50),]
                            for(i in 1 :nrow(dd)) {
                              y = 1/(dd[i,'b_Intercept']+dd[i,'b_I1DsE']/ests)
                              LRP_buff = c(LRP_buff, ests[which.min(abs(y-mL*1.25))])
                            }     
                            
                            b = aggregate(Pred~sE,data=tos, FUN=function(x) quantile(x,c( 0.025,0.975)))
                            
                       
                            #####buffer is the recommended                       
                            #Buff
                            ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbar(data=g5, aes(est1000,y=MEAN_T, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
                              geom_ribbon(data=b, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
                              labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)+
                              geom_vline(xintercept=LRP_buff,colour='red',alpha=.005)+geom_hline(yintercept=mL*1.25,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP_buff),colour='black')
                            
                            
      
mean(LRP_buff); quantile(LRP_buff, c(0.025,0.975)) 

#biomass and LRP
ggplot(subset(ind35),aes(x=year,y=estB/1000,ymin=lwrB/1000,ymax=uprB/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+
  theme_test(base_size = 14)+labs(x='Year',y='Commercial Biomass (x000) ')+
  geom_hline(yintercept=LRP_buff,colour='red',alpha=.005)+geom_hline(yintercept=mean(LRP_buff),colour='black')

gll = merge(ind35,gl)

gll$relF = gll$LFA35/(gll$estB/1000)
gll$relF_l = gll$LFA35/(gll$lwrB/1000)
gll$relF_u = gll$LFA35/(gll$uprB/1000)

#relf
ggplot(subset(gll),aes(x=year,y=relF,ymin=relF_l,ymax=relF_u))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality ')

saveRDS(list(ind35,LRP_buff,g5,b,mL,gll),file='LFA35Index_RPs_Feb18_2025.rds')

