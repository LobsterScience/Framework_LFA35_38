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
io = readRDS('IndicesFromFullComboModelJan192026_2000+.rds')
ind38 = io[[3]]

##converting number to biomass using fall length frequencies from all fall surveys in bay of fundy

yp = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=F,species=2550,biomass = F)
yp2 = subset(yp, month(SET_DATE)>8)

yp2 = aggregate(SA_CORRECTED_PRORATED_N~FISH_LENGTH,data=yp2,FUN=mean)
names(yp2)[2]='meanden'
yp2$Prop = yp2$meanden/sum(yp2$meanden)
yp2$Weight = lobLW(CL=yp2$FISH_LENGTH,sex = 1)
yp2s = subset(yp2,Prop>0)
mnW = sum(yp2s$Prop*yp2s$Weight)
ggplot(yp2,aes(x=FISH_LENGTH,y=Prop))+geom_bar(stat='identity')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)

#commercial abundance to weight
ind38$estB = ind38$est*mnW/1000
ind38$lwrB = ind38$lwr*mnW/1000
ind38$uprB = ind38$upr*mnW/1000

##LRP brecover
br = min(ind38$estB[which(ind38$year<2024)])/1000

#biomass and LRP
ggplot(subset(ind38),aes(x=year,y=estB/1000,ymin=lwrB/1000,ymax=uprB/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+
  theme_test(base_size = 14)+labs(x='Year',y='Commercial Biomass (x000) ')+
  geom_hline(yintercept=br,colour='red')



#rel F
ind38$year=ind38$year+1

g = lobster.db('seasonal.landings')
ga = lobster.db('annual.landings')
ga = subset(ga,select=c('YR','LFA38B'))
ga = ga[order(ga$YR),]
g$YR= as.numeric(substring(g$SYEAR,6,9))
g = subset(g,select=c(YR,LFA38))
ga = na.zero(ga, 'LFA38B')
gg = merge(ga,g)
gg$LFA38 = gg$LFA38+gg$LFA38B
gg$year = gg$YR
gll = merge(ind38,gg[,c('year','LFA38')])

gll$relF = gll$LFA38/(gll$estB/1000)
gll$relF_l = gll$LFA38/(gll$lwrB/1000)
gll$relF_u = gll$LFA38/(gll$uprB/1000)



#relf
ggplot(subset(gll),aes(x=year,y=relF,ymin=relF_l,ymax=relF_u))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality ')

saveRDS(list(ind38,br,gll),file='LFA38Index_RPs_Jan21_2026.rds')

