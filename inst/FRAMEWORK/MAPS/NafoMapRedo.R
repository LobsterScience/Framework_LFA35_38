require(ggplot2)
require(sf)
require(bio.lobster)
require(bio.utilities)
require(tidyr)
require(dplyr)
require(dbscan)
require(devtools)



###Maritimes Region with  NAFO
cents = readRDS(file.path(project.datadirectory("bio.lobster"), "data","maps","LFALabelsSF.rds"))
cents$geometry[cents$label == "36"] <- st_sfc(st_point(c(-66, 45.15)))
cents$geometry[cents$label == "35"] <- st_sfc(st_point(c(-64.85, 45.22)))


##### With Mike's rework of LFAs 

##LFAs
lobarea = st_read(file.path(project.datadirectory("bio.lobster"), "data","maps","lfaRework","LFA_rework_mmm.shp"))
# Define the bounding box
bbox <- st_bbox(c(xmin = -68, xmax = -54.75, ymin = 40, ymax = 48), crs = st_crs(lobarea))
# Clip the lobarea data
lobarea_clipped <- st_crop(lobarea, bbox)


##Area 38b
L38b<-read.csv("C:/Users/HowseVJ/Documents/bio.data/bio.lobster/data/maps/lfaRework/shape38b.csv")
sf38<-st_as_sf(L38b, coords=c("X","Y"),crs=4326)
polygon_38b <- st_coordinates(sf38)
polygon_38b <- rbind(polygon_38b, polygon_38b[1, ])# Ensure the polygon is closed by repeating the first point at the end
polygon_38b <- st_polygon(list(polygon_38b))# Create a POLYGON object
sf_38b <- st_sfc(polygon_38b, crs = 4326)# Convert to an sf object



p=ggLobsterMap(area="custom", ylim=c(40,48), xlim=c(-68,-54.7), addGrids=F,
                 fill.colours = 'grey',bathy=T,color=NA,addNAFO = T, addNAFOLabels = F, nafo=c('5Y','4X','5Z'),colourLFA = F,addLFAlines=F)+
  geom_sf(data=lobarea_clipped, fill=NA,color="black",linewidth=0.75)+
  geom_sf(data =sf_38b,color="black",linewidth=0.7)+
  geom_sf_text(data=cents, aes(label=label),family='sans',fontface="bold", size=3,col="black")+
  annotate("text",y=44.41, x=-67.12,label = "B",family='sans',fontface="bold",size=2.7,col="black")+
  annotate("text", y=42.5, x=-63,label = "41",family='sans',fontface="bold",size=3,col="black")+
  annotate("text", y=44, x=-67.75,label = "5Y",family='sans',fontface="bold",size=3.25,col="orange")+
  annotate("text", y=41, x=-67.6,label = "5Z",family='sans',fontface="bold",size=3.25,col="orange")+
  annotate("text", y=41.5, x=-64.5,label = "4X",family='sans',fontface="bold",size=3.25,col="orange")

 p 
 
 ggsave(filename = "UpdatedMap_Nafo.png", plot = p, path = "C:/Users/HowseVJ/OneDrive - DFO-MPO/LFA 35-38 Framework Resources/Figures", width = 12, height = 8, units = "in", dpi = 300)
 
