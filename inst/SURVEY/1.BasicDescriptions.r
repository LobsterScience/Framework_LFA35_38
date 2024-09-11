##ILTS and RV survey overviews

require(bio.lobster)
require(devtools)
require(bio.utilities)
require(sf)
require(ggplot2)
require(tidyr)
theme_set(theme_test(base_size = 14))

p = ggLobsterMap('34-38',addGrids = F)


##RV Survey
x = RV_sets()
x = subset(x,month(x$DATE) %in% c(6,7,8))
xle = x %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
xle= na.zero(xle,cols='value')
xxa = aggregate(value~Length+YEAR,data=xle,FUN=mean)
ggplot(subset(xxa,YEAR>1998),aes(x=Length,y=value))+geom_bar(stat='identity')+facet_wrap(~YEAR,scales = 'free_y')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')

xs = st_as_sf(x,coords = c('LONGITUDE','LATITUDE'),crs=4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326

xs = st_join(xs,rL,join=st_within)
xs = subset(xs,LFA %in% c(34,35,36,37,38))

ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(xs,YEAR<1982),fw='~YEAR')
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(xs,YEAR>1981 & YEAR<1994),fw='~YEAR')
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(xs,YEAR>1993 & YEAR<2006),fw='~YEAR')
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(xs,YEAR>2005 & YEAR<2018),fw='~YEAR')
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(xs,YEAR>2017 & YEAR<2024),fw='~YEAR')



xle = xs %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
xle= na.zero(xle,cols='value')
xxa = aggregate(value~Length+YEAR,data=xle,FUN=mean)
ggplot(subset(xxa,YEAR>1998),aes(x=Length,y=value))+geom_bar(stat='identity')+facet_wrap(~YEAR,scales = 'free_y')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')



################################################################################################################################################
#####ILTS


y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = F)
y = subset(y,select=c(TRIP_ID,SET_NO, SET_LONG,SET_LAT,SET_DATE,SA_CORRECTED_PRORATED_N))    

survey = y
survey = subset(survey, month(SET_DATE) %in% c(6:8))
survey$YEAR = year(survey$SET_DATE)
survey = st_as_sf(survey,coords = c('SET_LONG','SET_LAT'),crs=4326)
survey$CommN = survey$SA_CORRECTED_PRORATED_N

pt = subset(survey,YEAR==1996 )
pt$YEAR='1995-2012'
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(pt),fw='~YEAR')
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(survey,YEAR %in% 2013:2023 ),fw='~YEAR')


yp = ILTS_ITQ_All_Data(redo_base_data = F,size=c(0,300),aggregate=F,species=2550,biomass = F)
yp1 = subset(yp, month(SET_DATE) %in% 6:8)

yp1 = aggregate(SA_CORRECTED_PRORATED_N~YEAR+LFA+FISH_LENGTH,data=yp1,FUN=mean)
ggplot(yp1,aes(x=FISH_LENGTH,y=SA_CORRECTED_PRORATED_N))+geom_bar(stat='identity')+facet_wrap(~YEAR)+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')


yp2 = subset(yp, month(SET_DATE)>8)

yp2 = aggregate(SA_CORRECTED_PRORATED_N~YEAR+LFA+FISH_LENGTH,data=yp2,FUN=mean)
ggplot(subset(yp2,YEAR>2014),aes(x=FISH_LENGTH,y=SA_CORRECTED_PRORATED_N))+geom_bar(stat='identity')+facet_wrap(~YEAR)+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)+geom_vline(xintercept = 82.5,color='red')




#################################################################################################################################
###fall survey 
y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=T,species=2550,biomass = F)
survey$YEAR = year(survey$SET_DATE)
survey = subset(y, month(SET_DATE) >8 & YEAR>2014)
survey = st_as_sf(survey,coords = c('SET_LONG','SET_LAT'),crs=4326)
survey$CommN = survey$SA_CORRECTED_PRORATED_N
ggLobsterMap('34-38',addGrids = F,addPoints = T,pts=subset(survey,YEAR %in% 2013:2023 &Survey=='ILTS'),fw='~YEAR')




