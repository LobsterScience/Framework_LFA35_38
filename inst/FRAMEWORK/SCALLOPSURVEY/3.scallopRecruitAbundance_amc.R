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

####amc rewrote Jan 23 2025 for correct index and boot strapping
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
    a = lobster.db('scallop')
    
    size_range=c(70,82)
    lfa = 36            ###### PICK THE LFA  
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
  
    scallop.tows=a[[1]]
    scallopSurv = a[[2]]
    scallop.tows$Y = convert.dd.dddd(scallop.tows$START_LAT)
    scallop.tows$X = convert.dd.dddd(scallop.tows$START_LONG)
    scT = subset(scallop.tows,select=c('TOW_SEQ','TOW_DATE','STRATA_ID','X','Y'))
    
    scC = subset(scallopSurv, c(MEAS_VAL>=size_range[1] & MEAS_VAL<=size_range[2]),select=c("TOW_SEQ",  "ABUNDANCE_STD","MEAS_VAL", "SEX_ID") )
    scC = aggregate(ABUNDANCE_STD~TOW_SEQ,data=scC,FUN=sum)
    scC = merge(scT,scC,all.x=T)
    scC = bio.utilities::na.zero(scC,'ABUNDANCE_STD')
    scC$YEAR = lubridate::year(scC$TOW_DATE)
    
    #year subset
    sa = subset(scC,YEAR>=1999,select=c(TOW_SEQ,TOW_DATE,X,Y,ABUNDANCE_STD,YEAR))
    
    totS = st_as_sf(sa,coords = c('X','Y'),crs=st_crs(4326))
    
    
    ss = st_join(totS,v, join=st_within)
    xx = subset(ss,LFA==lfa)
    
    
    require(bio.survey)
    
    
    gr = unique(xx$YEAR)
    stra = subset(v,select=c(STRATA_ID,total_area_r))
    stra$NH = as.numeric(stra$total_area_r)
    stra$Strata = stra$STRATA_ID
    
    st_geometry(stra) <- NULL
    
    out = data.frame(yr=NA,LFA=NA,yst=NA,yst.se=NA,ci.yst.l=NA,ci.yst.u=NA,Nsets=NA,NsetswithLobster=NA)
    
    for(i in 1:length(gr)){
      n = subset(xx,YEAR==gr[i]) 
      n$Strata = n$STRATA_ID
      ui = unique(n$STRATA_ID)
      j = subset(stra,STRATA_ID %in% ui,)  
      st = list(Strata=j$Strata,NH=j$NH)
      
      sc = subset(n,select=c(TOW_SEQ,ABUNDANCE_STD,Strata))
      oldClass(sc) <- c("data.frame", "strata.data")
      sW = Stratify(sc,st,sc$ABUNDANCE_STD)
      out[i,1:3]=c(unique(n$YEAR),lfa,0)
      if(sum(unlist(sW$yhi))>0){
      ssW = summary(sW)
      bsW=summary(boot.strata(sW,method='BWR',nresamp=1000),ci.method='Percentile',gini=F)
      nt  = sum(sW$Nh)
      out[i,3:8] = c(ssW[[1]],ssW[[2]],bsW[[1]][1],bsW[[1]][2],sum(sW[['nh']]),sum(sW[['nhws']]))
      }
    }
    
    #g8 = ggplot(out,aes(x=yr,y=yst,ymin=ci.yst.l,ymax=ci.yst.u))+geom_point()+geom_line()+geom_errorbar(width=0)+ theme_test(base_size = 14)+labs(y='Stratified Mean Abundance',x='Year')
    #g5 = ggplot(out,aes(x=yr,y=yst,ymin=ci.yst.l,ymax=ci.yst.u))+geom_point()+geom_line()+geom_errorbar(width=0)+ theme_test(base_size = 14)+labs(y='Stratified Mean Abundance',x='Year')
    g6 = ggplot(out,aes(x=yr,y=yst,ymin=ci.yst.l,ymax=ci.yst.u))+geom_point()+geom_line()+geom_errorbar(width=0)+ theme_test(base_size = 14)+labs(y='Stratified Mean Abundance',x='Year')
    
    require(patchwork)
    (g5|g6)/g8 +plot_layout(heights=c(1,1),widths = c(1,1))
    
    