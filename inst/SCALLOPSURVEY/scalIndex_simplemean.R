library(sf)
require(bio.lobster)
require(bio.utilities)
require(Mar.utils)
require(devtools)
require(ggplot2)
require(roxygen2)
require(geosphere)
require(lubridate)
require(boot)
require(dplyr)
require(tidyr)



###### Bring in Survey Data #######
# import bycatch data from scallop survey
lobster.db('scallop')
size_range=c(70,82)
sex=0:3
lfa = 35      ###### PICK THE LFA  
layerDir=file.path(project.datadirectory("bio.lobster"), "data","maps")
return_sets_object=F


#polygons LFA
#rL = read.csv(file.path(layerDir,"Polygons_LFA.csv"))
#rL$label=rL$LFA
#rdef=Mar.utils::df_to_sf(rL,lat.field = "Y",lon.field = "X")  ## Make Polygons an sf object
#rdef<-merge(rdef,unique(rL[,c("PID", "LFA")]), all.x=T, by="PID")
rL = readRDS(file.path(layerDir,"Polygons_LFA.rds"))
rdef = st_as_sf(rL)

#subset polygons for target LFA
lf = subset(rdef,LFA==lfa)

scallop.tows$Y = convert.dd.dddd(scallop.tows$START_LAT)
scallop.tows$X = convert.dd.dddd(scallop.tows$START_LONG)
scT = subset(scallop.tows,select=c('TOW_SEQ','TOW_DATE','X','Y'))

scC = subset(scallopSurv,select=c("TOW_SEQ",  "ABUNDANCE_STD","MEAS_VAL", "SEX_ID")) ## Abundance_STD is the number of lobster per standard tow - standardized tow is 800m long by 5.334m wide ( 800* 5.334m = 4267 METERS SQUARED)

scC = merge(scT,scC,all.x=T)
scC$YEAR = lubridate::year(scC$TOW_DATE)
#year subset
s = subset(scC,YEAR>=1999)
i = which(is.na(s$ABUNDANCE_STD))
s$ABUNDANCE_STD[i]=0

#### Pull out Recruits
s$recruits[s$MEAS_VAL >= size_range[1]& s$MEAS_VAL<=size_range[2]]<-"yes"
s$recruits[s$MEAS_VAL < size_range[1]& s$MEAS_VAL<size_range[2]]<-"no"
s$recruits[is.na(s$recruits)]<-"no"

## Recruits per tow where they exist
recruits<- s %>% group_by(TOW_SEQ) %>% filter(recruits=="yes") %>% summarize(recruit_abun_std=sum(ABUNDANCE_STD))  

##need all tows with lat/longs

scT$YEAR<- year(scT$TOW_DATE)        
scT = subset(scT,YEAR>=1999)

allandrecruit<-merge(scT,recruits, by=c("TOW_SEQ"), all.x = T) ## left join

allandrecruit$recruit_abun_std[is.na(allandrecruit$recruit_abun_std)]<-0



## To get recruit per m2 would be  recruit_abund_std / 4267 per tow ---- this results in VERY SMALL NUMBERS - the issue is scale not so much data?
allandrecruit$ABUNDANCE_STD_PRU<-(allandrecruit$recruit_abun_std/4267) ## scaled for swept area to report on m2 or km2 because  data is already in standard tows


sa<-allandrecruit


totS = st_as_sf(sa,coords = c('X','Y'),crs=st_crs(4326)) ### Scallop data with position of sampling
if(return_sets_object) {return(totS); stop()}
ss = st_join(totS,lf, join=st_within) ### Adding the merged scallop data with the LFA 
xx = subset(ss,LFA==lfa)


###95% Cis - Bootstrap on the sets within strata

mean_function <- function(data, indices) {
  d <- data[indices, ]  
  return(mean(d$ABUNDANCE_STD_PRU))
}

# Perform bootstrapping
boot_results <- xx %>%
  group_by(YEAR) %>%
  do(boot_out = boot(data = ., statistic = mean_function, R = 1000))

# Extract confidence intervals
ci_results <- boot_results %>%
  rowwise() %>%
  mutate(
    ci = list(boot.ci(boot_out, type = "perc")$percent),
    ci_lower = ci[4],  # 4th element is the lower bound
    ci_upper = ci[5]   # 5th element is the upper bound
  ) %>%
  select(YEAR, ci_lower, ci_upper)




x = aggregate(ABUNDANCE_STD_PRU~YEAR,data=xx,FUN=mean) ## mean abundance by year  
names(x) <- c('YEAR', 'Index')
x$Index<-x$Index*10^6  #mean number of recruits per km2

##Variance
#simpVar<- aggregate(ABUNDANCE_STD_PRU~YEAR,data=xx,FUN=sd)
#x$SD<-simpVar$ABUNDANCE_STD_PRU*10^6 
#x$ci_upper<-x$Index+x$SD
#x$ci_lower<-x$Index-x$SD



ci_results$ci_lower<-ci_results$ci_lower*10^6
ci_results$ci_upper<-ci_results$ci_upper*10^6

recxx <- as.data.frame(do.call(cbind, rmed(x$YEAR, x$Index)))
names(recxx) <- c('YEAR', 'Rmed')
out <- merge(x, recxx, by="YEAR")
sc <- merge(out, ci_results, by = "YEAR")  
sc$Index <- as.numeric(sc$Index)


### Stratified mean abundance or recruits standardized for swept area  as mean number per Km^2   


ggplot(data = sc, aes(x = YEAR,y = as.numeric(Index))) +
  geom_point(size=2, colour = '#003f5c')+
  geom_line(data=out,aes(x=YEAR,y=Rmed),colour='#f95d6a',lwd=1.25)+
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2)+
  labs(x = "Year", y = expression("Scallop Survey Recruit Abundance (Mean Number /  "~Km^2~")"))+
  theme_test()

#ggplot(data = x, aes(x = YEAR,y = as.numeric(Index))) +
#  geom_point(size=2, colour = '#003f5c')+
#  geom_line(data=out,aes(x=YEAR,y=Rmed),colour='#f95d6a',lwd=1.25)+
#  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2)+
#  labs(x = "Year", y = expression("Scallop Survey Recruit Abundance (Mean Number /  "~Km^2~")"))+
#  theme_test()



#### Plot the recruit index for each LFA     
theme_set(theme_test(base_size = 14))

recruitscal<-ggplot(data = sc, aes(x = YEAR, y = Index)) +
  geom_point(size = 2, colour = 'black') +
  geom_line(aes(y = Rmed), colour = '#f95d6a', lwd = 1.25) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(x = "Year", y = expression("Scallop Survey Recruit Abundance (Mean Number /  "~m^2~")"))+
  scale_x_continuous(limits = c(min(sc$YEAR), max(sc$YEAR) + 1))



## Change file name to match lfa 
#ggsave(filename = "scalRecruitsLFA35.png", plot = recruitscal, path = "C:/Users/HowseVJ/OneDrive - DFO-MPO/LFA 35-38 Framework Resources/Figures", width = 8, height = 6, units = "in", dpi = 300)


