# data building

require(bio.survey)
require(bio.lobster)
require(bio.groundfish)
require(bio.utilities)
require(ggplot2)
require(dplyr)
require(cowplot)
require(magick)
library(ggplotify)


p=list()
assessment.year = 2024 ##Check Year
p$syr = 1989
p$yrs = p$syr:assessment.year



#######  Data Import   #######

#lobster.db("season.dates.redo")
#lobster.db( DS = "logs.redo",  p=p)   # Offshore logs monitoring documents
#logsInSeason = lobster.db("process.logs.redo")
#lobster.db( DS = 'annual.landings.redo', p=p) #static annual landings tabke needs to be updated by CDenton
#lobster.db( DS = 'seasonal.landings.redo', p=p) #static seasonal landings table needs to be updated by CDenton
#lobster.db( DS = 'annual.landings', p=p) #static annual landings tabke needs to be updated by CDenton
#lobster.db( DS = 'seasonal.landings', p=p) #static seasonal landings table needs to be updated by CDenton
#lobster.db('percent_reporting.redo')

#lobster.db('logs')
logs=lobster.db("process.logs")
land =lobster.db("annual.landings")
Sland = lobster.db('seasonal.landings')
percentReport=lobster.db('percent_reporting')

#Oct 2023 to july 2024
PR_lfa35<-percentReport[percentReport$YEARMTH >= 202310 & percentReport$YEARMTH <= 202407,c("YEARMTH","L35MISS","L35RECD","L35PERCENT")]
PR_lfa35 <- PR_lfa35[order(PR_lfa35$YEARMTH), ]



#RDS objects
L35data<-readRDS("C:/Users/HowseVJ/Documents/bio.data/Framework_LFA35_38/data_templateonly/LFA35Index_RPs_jan_20_2025.RDS")
#[1]  ind35
#[2] LRP_buff
#[3] g5
#[4] b
#[5] mL
#[6] gll  

######## SET UP ######## 
#plotting as per csasdown 4 panel plot
#add in the theme_csas
theme_csas <- function(base_size = 14, base_family = "", text_col = "grey20",
                       panel_border_col = "grey70") {
  half_line <- base_size / 2
  theme_light(base_size = base_size, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = text_col),
      strip.text.y = element_text(colour = text_col),
      axis.text = element_text(colour = text_col),
      axis.title = element_text(colour = text_col),
      legend.title = element_text(colour = text_col, size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = panel_border_col, linewidth = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.9), colour = text_col),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = text_col, size = rel(1)),
      plot.subtitle = element_text(colour = text_col, size = rel(.85))
    )
}

#########  CATCH AND EFFORT  #########

a=lobster.db('seasonal.landings')
a$yr= as.numeric(substr(a$SYEAR,6,9))
aaa = a
aaa <- aaa %>%
  mutate(SYEAR = ifelse(row_number() == n(), "2024-2025", SYEAR),
         yr = ifelse(row_number() == n(), 2025, yr))
b = lobster.db('process.logs')
b = subset(b,SYEAR %in% 2005:2025 & LFA =='35') 

aa = split(b,f=list(b$LFA,b$SYEAR))
cpue.lst<-list()

for(i in 1:length(aa)){
  tmp<-aa[[i]]
  tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS')]
  names(tmp)<-c('time','catch','effort')
  tmp$date<-as.Date(tmp$time)
  first.day<-min(tmp$date)
  tmp$time<-julian(tmp$date,origin=first.day-1)
  tmp$time = ceiling(tmp$time/7) #convert to week of season
  g<-as.data.frame(biasCorrCPUE(tmp,by.time=F))
  g$lfa=unique(aa[[i]]$LFA)
  g$yr = unique(aa[[i]]$SYEAR)
  g = t(g)[,1]
  cpue.lst[[i]] <- g
}

cc =as.data.frame(do.call(rbind,cpue.lst))
cc$CPUE = as.numeric(cc$`biasCorrCPUE(tmp, by.time = F)`)
cc$yr = as.numeric(cc$yr)
cp = as.data.frame(do.call(cbind,rmed(cc$yr,cc$CPUE)))
#effort
ef = merge(cc,aaa)
ef$Effort = ef$LFA35/(ef$CPUE)


ymax <- 5200
scaleright <- max(ef$Effort) / ymax
aaa <- head(aaa, -1)

##Landings description numbers
landef<-aaa[, c("SYEAR", "yr","LFA35")]
filtered_df <- landef[landef$SYEAR >= 1994 & landef$SYEAR <= 2024, ]
median_landings <- median(filtered_df$LFA35, na.rm = TRUE)




plot1_catchef <- ggplot(data = subset(aaa, yr > 1990), aes(x = yr, y = LFA35)) +
  geom_bar(stat = 'identity', aes(fill = ifelse(yr == 2024, 'gray66', 'black')), show.legend = FALSE) +
  geom_point(data = ef, aes(x = yr, y = Effort / scaleright, colour = ifelse(yr == 2024, 'grey66', 'black')), shape = 17, size = 2.5, show.legend = FALSE) +
  geom_line(data = ef, aes(x = yr, y = Effort / scaleright, linetype = ifelse(yr == 2024, 'solid', 'dashed')), colour = 'black', linewidth = 0.75, show.legend = FALSE) +
  geom_point(data = ef, aes(x = yr, y = Effort / scaleright, colour = ifelse(yr == 2024, 'grey66', 'black')), shape = 16, size = 2, show.legend = FALSE) +
  scale_y_continuous(name = 'Catch (t)', sec.axis = sec_axis(~ . * scaleright, name = 'Effort', breaks = seq(0, 2000, by = 250))) +
  labs(x = "Fishing Year (ending)") +
  scale_fill_identity() +
  scale_colour_identity() +
  theme_csas()

#########  UNBIASED CPUE  #########


aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>2005 & SYEAR<2024 & LFA %in% c(35,36,38))


aa = split(aT,f=list(aT$LFA,aT$SYEAR))
cpue.lst<-list()
cpue.ann<- list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS')]
  names(tmp)<-c('time','catch','effort')
  tmp$date<-as.Date(tmp$time)
  first.day<-min(tmp$date)
  tmp$time<-julian(tmp$date,origin=first.day-1)
  g<-as.data.frame(biasCorrCPUE(tmp,by.time=T))
  g$lfa=unique(aa[[i]]$LFA)
  g$yr = unique(aa[[i]]$SYEAR)
  gl = aggregate(effort~time, data=tmp, FUN=sum)
  g = merge(g,gl,by.x='t',by.y='time')
  cpue.lst[[i]] <- g
  
  g<-as.data.frame(t(biasCorrCPUE(tmp,by.time=F)))
  g=data.frame(g,LFA=as.numeric(unique(aa[[i]]$LFA)),SYEAR=as.numeric(unique(aa[[i]]$SYEAR)))
  cpue.ann[[i]]=g
}

ca =as.data.frame(do.call(rbind,cpue.ann))


plot1_cpue <- ggplot(subset(ca,LFA %in% c(35)),aes(x=SYEAR,y=unBCPUE))+
              geom_point()+
              geom_line()+
              geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.3)+
              xlab('Fishing Year (ending)')+
              ylab('Unbiased CPUE')+
              theme_csas()
  
  
#########  BIOMASS & LRP #########

biomass35<-L35data[[1]]
LRP_BUFF<-L35data[[2]] ## For the mean  to get the LRP buffer



plot1_biomass<-ggplot(biomass35,aes(x=oyr,y=estB/1000,ymin=lwrB/1000,ymax=uprB/1000))+
  geom_point()+geom_line()+
  geom_ribbon(alpha=.25)+
  labs(x='Year',y='Commercial Biomass (x000) ')+
  geom_hline(yintercept=LRP_BUFF,colour='red',alpha=.005)+
  geom_hline(yintercept=mean(LRP_BUFF),colour='black')+
  theme_csas()


#########  RELATIVE FISHING MORTALITY NO LRP ##########    ###############     FIX    THE   YEARS?? 
relF35<-L35data[[6]]

plot1_relF<-ggplot(relF35,aes(x=year,y=relF,ymin=relF_l,ymax=relF_u))+
  geom_point()+
  geom_line()+
  geom_ribbon(alpha=.25)+
  labs(x='Year',y='Relative Fishing Mortality ')+
  theme_csas()

cowplot::plot_grid( plot1_catchef,plot1_cpue, plot1_biomass,plot1_relF, ncol = 2, labels = "AUTO", align = "hv")




#### Plot two 


plot2_recruits <- ggplot() +
  labs(x = "Year", y = "Settler Density") +
  annotate("text", x = 2008, y = 4, label = "Data Not Available", size = 5, colour = 'grey') +
  scale_x_continuous(limits = c(1990, 2024)) +
  scale_y_continuous(limits = c(0, 8)) +
  theme_minimal() +
  theme_csas()


image_path <- "C:/Users/HowseVJ/OneDrive - DFO-MPO/LFA 35-38 Assessments- FSAR/Figures/BOF_RV_Temp.jpg"
plot2_temp <- image_read(image_path)
plot2_temp_grob <- as.grob(plot2_temp)
combined_plot <- plot_grid(plot2_recruits, plot2_temp_grob, ncol = 2, labels = "AUTO", align = "hv", axis = "tb")
