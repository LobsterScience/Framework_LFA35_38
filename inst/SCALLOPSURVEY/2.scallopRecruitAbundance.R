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


#### Import Inshore Scallop Survey Regular Domain shapefiles#####
## J.Sameoto Sept 2024 
## Note SPA 2 (South West and North West Banks -- at mouth of BoF only irregularly surveyed -- Scal_Area == "SPA2"
# Find where tempfiles are stored
temp <- tempfile()
# Download this to the temp directory
download.file("https://raw.githubusercontent.com/Mar-scal/GIS_layers/master/inshore_boundaries/inshore_boundaries.zip", temp)
# Figure out what this file was saved as
temp2 <- tempfile()
# Unzip it
unzip(zipfile=temp, exdir=temp2)
# Now read in the shapefiles
strata <- st_read(paste0(temp2, "/inshore_survey_strata/Inshore_Survey_Strata_update_sept2024/Scallop_Strata.shp"))
#note coordinate reference system 
st_crs(strata)
#view spatially 
plot(strata)


###### Bring in Survey Data #######


    # import bycatch data from scallop survey
    lobster.db('scallop')
    size_range=c(70,82)
    sex=0:3
    lfa = 38             ###### PICK THE LFA  
    layerDir=file.path(project.datadirectory("bio.lobster"), "data","maps")
    return_sets_object=F
    
    
#load Strata and Polygons    
    sdef<-strata
    
    #polygons LFA
    rL = read.csv(file.path(layerDir,"Polygons_LFA.csv"))
    rL$label=rL$LFA
    rdef=Mar.utils::df_to_sf(rL,lat.field = "Y",lon.field = "X")  ## Make Polygons an sf object
    rdef<-merge(rdef,unique(rL[,c("PID", "LFA")]), all.x=T, by="PID")
  
    
 #intersect LFA and strata
    lf = subset(rdef,LFA==lfa)
    
##Convert to match projection
    sdef <- st_transform(sdef, st_crs(lf))   
#intersect strata with lfa    
    v = st_intersection(sdef,lf)
    v$total_area_r = st_area(v) #area within LFA
    j = which(v$STRATA_ID ==49)
    if(length(j)>0) v$total_area_r[j] = v$total_area[j] = 239000000
  
    scallop.tows$Y = convert.dd.dddd(scallop.tows$START_LAT)
    scallop.tows$X = convert.dd.dddd(scallop.tows$START_LONG)
    scT = subset(scallop.tows,select=c('TOW_SEQ','TOW_DATE','STRATA_ID','X','Y'))
    
    scC = subset(scallopSurv,select=c("TOW_SEQ",  "ABUNDANCE_STD","MEAS_VAL", "SEX_ID"))
    
    scC = merge(scT,scC,all.x=T)
    scC$YEAR = lubridate::year(scC$TOW_DATE)
    #year subset
    s = subset(scC,YEAR>=1999)
    i = which(is.na(s$ABUNDANCE_STD))
    s$ABUNDANCE_STD[i]=0
    
    #add a dummy column for lobsters meeting a criteria for adjusting catch
    s$dumm= s$dumm2 = 0
    
    #size subset 
    j = which(s$ABUNDANCE_STD==0) #give me the zeros
    k = which(s$MEAS_VAL >= size_range[1]& s$MEAS_VAL<=size_range[2])
    s$dumm[c(j,k)] = 1
    
    #sex subset
    k = which(s$SEX_ID %in% sex) 
    s$dumm2[c(j,k)] = 1
    
    
    s$ABUNDANCE_STD_PRU = s$ABUNDANCE_STD * s$dumm * s$dumm2
    
    sa = aggregate(cbind(ABUNDANCE_STD_PRU)~TOW_SEQ+STRATA_ID+TOW_DATE+X+Y+YEAR,data=s,FUN=sum)
    
    totS = st_as_sf(sa,coords = c('X','Y'),crs=st_crs(4326))
    if(return_sets_object) {return(totS); stop()}
    ss = st_join(totS,v, join=st_within)
    xx = subset(ss,LFA==lfa)
    xx$STRATA_ID = xx$STRATA_ID.y 
    xx$STRATA_ID.y = xx$STRATA_ID.x = xx$total_area = NULL
    
    x = aggregate(ABUNDANCE_STD_PRU~STRATA_ID+YEAR+total_area_r,data=xx,FUN=mean) ## why is it aggregated by strata here rename this
    xa = aggregate(total_area_r~YEAR,data=x,FUN=sum)
    names(xa)[2]='total_area'
    xxa =  merge(x,xa)
 
    
   ###BEFORE 
    #    sca = aggregate(ABUNDANCE_STD_PRU*total_area_r/total_area~YEAR,data=xxa,FUN=sum)
    # names(sca)=c('YEAR','Index')
    # recxx = as.data.frame(do.call(cbind,rmed(sca$YEAR,sca$Index)))
    # names(recxx)=c('YEAR','Rmed')
    # out = merge(sca,recxx)
    # 
    #  sc=out
    # 
    # 
    #  
    #  ####Plot figures - CHECK LFAS #### 
    #  recruitscal35<-ggplot(data = sc, aes(x = YEAR,y = as.numeric(Index))) +
    #    geom_point(size=2, colour = 'black')+
    #    geom_line(data=sc,aes(x=YEAR,y=Rmed),colour='#f95d6a',lwd=1.25)+
    #    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    #    labs(x = "Year", y =expression("Scallop Survey Recruit Abundance (Mean Number/ m "^ 2 * ")"))+
    #    scale_x_continuous(limits = c(min(sc$YEAR), max(sc$YEAR) + 1))
    #  
     
     
##Bootstrap of 95% CI     
     mean_function <- function(data, indices) {
       d <- data[indices, ]  
       return(sum(d$ABUNDANCE_STD_PRU * d$total_area_r / d$total_area))
     }
   
     boot_results <- xxa %>%
       group_by(YEAR) %>%
       do(boot_out = boot(data = ., statistic = mean_function, R = 1000))
   
     
     # Extract confidence intervals correctly
     ci_results <- boot_results %>%
       rowwise() %>%
       mutate(
         ci = list(boot.ci(boot_out, type = "perc")$percent),
         ci_lower = ci[4],  # 4th element is the lower bound
         ci_upper = ci[5]   # 5th element is the upper bound
       ) %>%
       select(YEAR, ci_lower, ci_upper)
     
     
     sca <- aggregate(ABUNDANCE_STD_PRU * total_area_r / total_area ~ YEAR, data = xxa, FUN = sum)
     names(sca) <- c('YEAR', 'Index')
     
     recxx <- as.data.frame(do.call(cbind, rmed(sca$YEAR, sca$Index)))
     names(recxx) <- c('YEAR', 'Rmed')
     out <- merge(sca, recxx, by="YEAR")
     sc <- merge(out, ci_results, by = "YEAR")  
     sc$Index <- as.numeric(sc$Index)
     
     
theme_set(theme_test(base_size = 14))
    
recruitscal<-ggplot(data = sc, aes(x = YEAR, y = Index)) +
       geom_point(size = 2, colour = 'black') +
       geom_line(aes(y = Rmed), colour = '#f95d6a', lwd = 1.25) +
       geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(x = "Year", y = expression("Scallop Survey Recruit Abundance (Mean Number /  "~m^2~")"))+
       scale_x_continuous(limits = c(min(sc$YEAR), max(sc$YEAR) + 1))
     
     

## Change file name to match lfa 
ggsave(filename = "scalRecruitsLFA35.png", plot = recruitscal, path = "C:/Users/HowseVJ/OneDrive - DFO-MPO/LFA 35-38 Framework Resources/Figures", width = 8, height = 6, units = "in", dpi = 300)


