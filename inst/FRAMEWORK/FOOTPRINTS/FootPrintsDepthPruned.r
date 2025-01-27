require(devtools)
require(bio.lobster)
require(bio.utilities)
require(sf)
load_all('C:/Users/Cooka/Documents/git/bio.utilities')


la()
wd = ('C:\\Users\\Cooka\\Documents\\bio_data\\Framework_LFA35_38')
setwd(wd)


layerDir=file.path(project.datadirectory("bio.lobster"), "data","maps")

##from inst/IP/FootprintMapping/1.EstimatingTrapHaulsFromSlipsandSplit2Grids.r

Tot = readRDS('DiscretizedData/PrivacyScreened_TrapHauls_Landings_Trips_Gridand.rds')


#making plots of Tot

GrMap = readRDS(file.path(project.datadirectory('bio.lobster'),'data','maps','LFA27-38GridsPrunedtoDepth-sf.rds'))
coa = st_as_sf(readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastlineSF_NY_NL.rds")))

GrMap1 = GrMap
GrMap1$area = st_area(GrMap1)/1000000
GrMap1$V2 = paste(GrMap1$LFA, GrMap1$GRID_NO,sep="-")
st_geometry(GrMap1)<- NULL
gg = aggregate(area~LFA+GRID_NO,data=GrMap1,FUN=function(x) abs(sum(x)))

GrMap2 =merge(GrMap,gg)

gTot = merge(GrMap2,Tot,by.x=c('LFA','GRID_NO'),by.y=c('LFA','Grid'),all.x=T)
gTot = subset(gTot, LFA %in% c(35,36,38))
###################################################################################################################

makeNewFootPrints <- function(x,data,grmap=GrMap){
  LFA1=x 
  gTot=data
  gTot$LandingsT = gTot$Landings/1000
  g27 = subset(gTot,LFA==LFA1 & FishingYear>2018 )
        g27p = subset(gTot,LFA==LFA1 & FishingYear>2018 &PrivacyScreen==1)
        g27n = subset(gTot,LFA==LFA1 & FishingYear>2018 &PrivacyScreen==0)
        
        ok1 = ggplot(g27p,aes(fill=TrapHauls))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
          coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        
        
        width_inch <- 8.5 * 0.7  # 80% of the paper width
        height_inch <- 11 * 0.7  # 80% of the paper height
        
        # Save the plot with the specified dimensions
        nm = paste('LFA', LFA1,'- Logbook reported number of trap hauls by grid.png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 

        g27p$TrapHaulspkm2 = as.numeric(g27p$TrapHauls/g27p$area)
        
        ok1 = ggplot(g27p,aes(fill=TrapHaulspkm2))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
                  coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        
        
        # Save the plot with the specified dimensions
        nm = paste('LFA', LFA1,'- Logbook reported number of trap hauls (TH) by grid corrected by the area of the grid (km2).png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
        
        
          ok1 = ggplot(g27p,aes(fill=LandingsT))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
            geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
            coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
            scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
            scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
          
        
        
        # Save the plot with the specified dimensions
          nm = paste('LFA',LFA1,'- Logbook reported Landings (t) by grid.png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
          
        g27p$Landingspkm2 = as.numeric(g27p$LandingsT/g27p$area)
        
        ok1 = ggplot(g27p,aes(fill=Landingspkm2))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
                    coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        
        
        
        # Save the plot with the specified dimensions
        nm = paste('LFA',LFA1,'- Logbook reported Landings (t) by grid adjusted by area (km2).png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
        
        
        ok1 = ggplot(g27p,aes(fill=Trips))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
                    coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        
        
        
        # Save the plot with the specified dimensions
        nm = paste('LFA',LFA1,'- Logbook reported number of trips by grid.png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
          
        g27p$Tripspkm2 = as.numeric(g27p$Trips/g27p$area)
        
        ok1 = ggplot(g27p,aes(fill=Tripspkm2))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
          coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        
        
        # Save the plot with the specified dimensions
        nm = paste('LFA',LFA1,'- Logbook reported number of trips by grid adjusted by area (km2).png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
        
        ok1 = ggplot(g27p,aes(fill=NLics))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
          coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        
        
        # Save the plot with the specified dimensions
        nm = paste('LFA', LFA1,'- Logbook reported number of licences fishing by grid.png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
        
        g27p$Licencespkm2 = as.numeric(g27p$NLics/g27p$area)
        
        ok1 = ggplot(g27p,aes(fill=Licencespkm2))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
          coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
                             
        
        # Save the plot with the specified dimensions
        nm = paste('LFA',LFA1,'- Logbook reported number of licences fishing by grid corrected by area (km2).png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
        
        g27p$CPUE = g27p$Landings/g27p$TrapHauls
        ok1 = ggplot(g27p,aes(fill=CPUE))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~FishingYear)+
          geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
            coord_sf(xlim = c(st_bbox(g27)$xmin,st_bbox(g27)$xmax),
                   ylim = c(st_bbox(g27)$ymin,st_bbox(g27)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        
        nm = paste('LFA',LFA1,'- Logbook reported catch (kg) per unit effort (trap haul) by grid.png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
        
        
#Trends
        g27 = subset(gTot,LFA==LFA1 & FishingYear>2005 )
        g27p = subset(gTot,LFA==LFA1 & FishingYear>2005 &PrivacyScreen==1)
        g27n = subset(gTot,LFA==LFA1 & FishingYear>2005 &PrivacyScreen==0)
        
        g27$CPUE = g27$Landings/g27$TrapHauls
        
        ui = unique(g27$GRID_NO)
        out = list()
        j=0
        for(i in 1:length(ui)){
            b=subset(g27,GRID_NO==ui[i])
            b1 = subset(g27,GRID_NO==ui[i]& FishingYear==max(b$FishingYear))
            bp=aggregate(CPUE~FishingYear,data=b,FUN=mean)
            
            base = mean(bp$CPUE[which(bp$FishingYear %in% 2006:2010)])
            p4b = p3b = p2b = p1b = b1
            p1b$PercentDiff = percentDifference(c(mean(bp$CPUE[6:8]),base))
            p2b$PercentDiff = percentDifference(c(mean(bp$CPUE[9:11]),base))
            p3b$PercentDiff = percentDifference(c(mean(bp$CPUE[12:14]),base))
            p4b$PercentDiff = percentDifference(c(mean(bp$CPUE[15:17]),base))
            p1b$label = '2011-2013'
            p2b$label = '2014-2016'
            p3b$label = '2017-2019'
            p4b$label = '2020-2022'
            if(nrow(bp)<14){
              p1b$PercentDiff =p2b$PercentDiff = p3b$PercentDiff = p4b$PercentDiff = NA
              }
            bp = do.call(rbind,list(p1b,p2b,p3b,p4b))
            out[[i]] = bp
            
        }
        o = do.call(rbind,out)
      
        ok1 = ggplot(o,aes(fill=PercentDiff))+
          geom_sf() +
          scale_fill_distiller(trans='identity',palette='Spectral') +
          facet_wrap(~label)+
          #geom_sf(data=g27n,fill='white')+  
          geom_sf(data=coa,fill='grey')+
          geom_sf(data=subset(grmap,LFA==LFA1),fill=NA)+
          coord_sf(xlim = c(st_bbox(o)$xmin,st_bbox(o)$xmax),
                   ylim = c(st_bbox(o)$ymin,st_bbox(o)$ymax),
                   expand = FALSE)+
          scale_x_continuous(breaks = c(round(seq(st_bbox(g27)$xmin,st_bbox(g27)$xmax,length.out=2),2)))+
          scale_y_continuous(breaks = c(round(seq(st_bbox(g27)$ymin,st_bbox(g27)$ymax,length.out=2),2)))
        nm = paste('LFA',LFA1,'- Logbook reported grid estimates of percent difference in catch per unit effort (CPUE) relative to mean CPUE between 2006 and 2010.png',sep="")
        ggsave(nm, plot = ok1, width = width_inch, height = height_inch, units = "in") 
}



makeNewFootPrints(x='35',data=gTot)

r<-readRDS(file.path( layerDir,"GridPolysSF.rds"))

b=subset(r,LFA %in% c(35,36,38))
o=subset(GrMap,LFA %in% c(34,35,36,38))

ggplot(b)+
         geom_sf()+
  geom_sf(data=coa,fill='grey')+
  geom_sf(data=o,fill='red')+
  coord_sf(xlim = c(st_bbox(b)$xmin,st_bbox(b)$xmax),
           ylim = c(st_bbox(b)$ymin,st_bbox(b)$ymax),
           expand = FALSE)
       


gTot$CPUE = gTot$Landings/gTot$TrapHauls
g27p = subset(gTot, FishingYear%in%2005:2014)
g27p$Landings_t = g27p$Landings/1000
ok1 = ggplot(g27p,aes(fill=Landings_t))+
    geom_sf() +
  scale_fill_distiller(trans='identity',palette='Spectral') +
  facet_wrap(~FishingYear)+
#  geom_sf(data=g27n,fill='white')+  
  geom_sf(data=coa,fill='grey')+
  geom_sf(data=GrMap,fill=NA)+
  coord_sf(xlim = c(st_bbox(b)$xmin,st_bbox(b)$xmax),
           ylim = c(st_bbox(b)$ymin,st_bbox(b)$ymax),
           expand=F)+
  #scale_x_continuous(breaks = c(round(seq(st_bbox(g27p)$xmin,st_bbox(g27p)$xmax,length.out=2),2)))+
  #scale_y_continuous(breaks = c(round(seq(st_bbox(g27p)$ymin,st_bbox(g27p)$ymax,length.out=2),2)))+
  theme_test_adam()


g27p = subset(gTot, FishingYear%in%2015:2024)
g27p$Landings_t = g27p$Landings/1000
ok1 = ggplot(g27p,aes(fill=Landings_t))+
  geom_sf() +
  scale_fill_distiller(trans='identity',palette='Spectral') +
  facet_wrap(~FishingYear)+
  #  geom_sf(data=g27n,fill='white')+  
  geom_sf(data=coa,fill='grey')+
  geom_sf(data=GrMap,fill=NA)+
  coord_sf(xlim = c(st_bbox(b)$xmin,st_bbox(b)$xmax),
           ylim = c(st_bbox(b)$ymin,st_bbox(b)$ymax),
           expand=F)+
  #scale_x_continuous(breaks = c(round(seq(st_bbox(g27p)$xmin,st_bbox(g27p)$xmax,length.out=2),2)))+
  #scale_y_continuous(breaks = c(round(seq(st_bbox(g27p)$ymin,st_bbox(g27p)$ymax,length.out=2),2)))+
  theme_test_adam()

########################
###how many grid constitute 50% of landings --- this data is not privacy screened

gl = subset(gTot,LFA %in% c(35,36,38))
glT = aggregate(Landings~LFA+FishingYear,data=gl,FUN=sum)
names(glT)[3]='TL'
gll = merge(gl,glT)
gll$pL = gll$Landings/gll$TL

gls = split(gll,f=list(gll$LFA,gll$FishingYear))
junk = data.frame(LFA=rep(NA,length(gls)),YR=rep(NA,length(gls)),Prop50=rep(NA,length(gls)),Prop75=rep(NA,length(gls)),Prop95=rep(NA,length(gls)),Prop40=rep(NA,length(gls)),Prop30=rep(NA,length(gls)),ngrids=rep(NA,length(gls)),H=rep(NA,length(gls)), E=rep(NA,length(gls)))
for(i in 1: length(gls)){
      j = gls[[i]]
      j = j %>% arrange(desc(pL))
      j$co=cumsum(j$pL)
      junk$LFA[i] = unique(j$LFA)
      junk$YR[i] = unique(j$FishingYear)
      junk$Prop30[i] = which.min(abs(j$co-.30))
      junk$Prop40[i] = which.min(abs(j$co-.40))
      
            junk$Prop50[i] = which.min(abs(j$co-.50))
      junk$Prop75[i] = which.min(abs(j$co-.75))
      junk$Prop95[i] = which.min(abs(j$co-.95))
      
      #evenness
      junk$ngrids[i] = nrow(j)
      junk$H[i] = -sum(j$pL * log(j$pL))
      junk$E[i] = junk$H[i]/junk$ngrids[i]
}

ggplot(junk,aes(x=YR,y=H))+geom_point()+facet_wrap(~LFA)+theme_test(base_size = 14)





################################
###############################
###survey biomass and first part of the season landings

glK = readRDS('DiscretizedData/LoglandingsFirsttoJan.rds')
glsm = aggregate(WEIGHT_KG~LFA+SYEAR,data=glK,FUN=sum)
names(glsm)[3] = 'TL'

ll = merge(glK,glsm)







################
gp = subset(g27p,FishingYear==2022)
gl = subset(g27p,FishingYear==2023)

gl$geometry<- NULL

gg = merge(gp,gl[,c('LFA','GRID_NO','CPUE')],by=c('LFA','GRID_NO'))


percent_diff <- function(row) {
  row$geometry<- NULL
  
  abs_diff <- (as.numeric(row[1]) - as.numeric(row[2]))
  mean_val <- mean(as.numeric(row))
  percent_diff <- (abs_diff / mean_val) * 100
  return(percent_diff)
}

gg$percentChange =  apply(gg[,c('CPUE.y','CPUE.x')],1,percent_diff)


require(colorspace)
 ggplot(subset(gg,PrivacyScreen==1),aes(fill=percentChange))+
  geom_sf() +
  scale_fill_continuous_diverging(palette='Purple-Green') +
  #facet_wrap(~FishingYear)+
  #  geom_sf(data=g27n,fill='white')+  
  geom_sf(data=coa,fill='grey')+
  geom_sf(data=GrMap,fill=NA)+
  coord_sf(xlim = c(st_bbox(g27p)$xmin,st_bbox(g27p)$xmax),
           ylim = c(st_bbox(g27p)$ymin,st_bbox(g27p)$ymax),
           expand = FALSE)+
  scale_x_continuous(breaks = c(round(seq(st_bbox(g27p)$xmin,st_bbox(g27p)$xmax,length.out=2),2)))+
  scale_y_continuous(breaks = c(round(seq(st_bbox(g27p)$ymin,st_bbox(g27p)$ymax,length.out=2),2)))
