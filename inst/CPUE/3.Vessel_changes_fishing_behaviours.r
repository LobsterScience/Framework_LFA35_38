#vessel and fishing characteristics

require(bio.lobster)
require(ggplot2)
require(bio.utilities)
require(devtools)
require(ggpubr)
load_all('~/git/bio.utilities')

theme_set(theme_test(base_size = 14))

#cpue index
co = readRDS(file.path(project.datadirectory('Framework_LFA35_38'),'outputs','CPUE','unBIASED_CPUE.rds'))
cpu = co[[1]]


d = lobster.db('vessels.by.port')
d = subset(d,LFA %in% 35:38 & YR_FISHED>2004 & YR_FISHED<2024)
d = na.zero(d)
d1 = aggregate(cbind(GROSS_TONNAGE, BHP, LOA, BREADTH, DEPTH,YEAR_BUILT)~VR_NUMBER+LFA+YR_FISHED,data=d,FUN=min)
d = na.zero(d1,rev=T)
d$YR_FISHED = as.numeric(d$YR_FISHED)
d$Age_of_Vessel = d$YR_FISHED-d$YEAR_BUILT
d$YR_FISHED = as.character(d$YR_FISHED)
l = aggregate(Age_of_Vessel~LFA,data=d,FUN=median)
oo = aggregate(Age_of_Vessel~YR_FISHED+LFA,data=d,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))
ocp = merge(oo,cpu,by.x=c('LFA','YR_FISHED'),by.y=c('LFA','SYEAR'))

ggplot(oo,aes(x=as.numeric(YR_FISHED),y=Age_of_Vessel[,2],ymin=Age_of_Vessel[,1],ymax=Age_of_Vessel[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Age of Vessel')+geom_hline(data=l,aes(yintercept=Age_of_Vessel),color='red')

ggplot(ocp,aes(x=unBCPUE,y=Age_of_Vessel[,2],label=YR_FISHED))+geom_point()+geom_path()+facet_wrap(~LFA)+
  xlab('CPUE')+ylab('Age of Vessel')+geom_text(data=subset(ocp,YR_FISHED %in% c(2006,2023)),aes(label=YR_FISHED,x=unBCPUE,y=Age_of_Vessel[,2]))

# v = aggregate(Age_of_Vessel~YR_FISHED,data=subset(d,LFA==35),FUN=length)
# ag = median(subset(d,LFA==35)$Age_of_Vessel,na.rm=T)
# ggplot(subset(d,LFA==35),aes(YR_FISHED,Age_of_Vessel))+geom_boxplot()+ geom_text(data=v,aes(YR_FISHED,Inf,label=Age_of_Vessel),vjust=1,size=3) + theme_test()+ylab('Age of Vessel')+xlab('Year Fished')+geom_hline(yintercept=ag,col='red',linewidth=1.1)
# 
# 
# v = aggregate(Age_of_Vessel~YR_FISHED,data=subset(d,LFA==36),FUN=length)
# ag = median(subset(d,LFA==36)$Age_of_Vessel,na.rm=T)
# ggplot(subset(d,LFA==36),aes(YR_FISHED,Age_of_Vessel))+geom_boxplot()+ geom_text(data=v,aes(YR_FISHED,Inf,label=Age_of_Vessel),vjust=1,size=3)+ theme_test()+ylab('Age of Vessel')+xlab('Year Fished')+geom_hline(yintercept=ag,col='red',linewidth=1.1)
# 
# 
# v = aggregate(Age_of_Vessel~YR_FISHED,data=subset(d,LFA==38),FUN=length)
# ag = median(subset(d,LFA==38)$Age_of_Vessel,na.rm=T)
# ggplot(subset(d,LFA==38),aes(YR_FISHED,Age_of_Vessel))+geom_boxplot()+ geom_text(data=v,aes(YR_FISHED,Inf,label=Age_of_Vessel),vjust=1,size=3)+ theme_test()+ylab('Age of Vessel')+xlab('Year Fished')+geom_hline(yintercept=ag,col='red',linewidth=1.1)

x=ggplot(d,aes(YEAR_BUILT,BHP))+geom_point()+theme_test()+facet_wrap(~LFA)+ xlab('Year Vessel Manufactured')+ylab('Brake Horse Power')
y=ggplot(d,aes(YEAR_BUILT,LOA))+geom_point()+theme_test()+facet_wrap(~LFA)+ xlab('Year Vessel Manufactured')+ylab('Length Overall')+geom_hline(yintercept=65,colour='red')
z=ggplot(subset(d,BREADTH>10),aes(YEAR_BUILT,BREADTH))+geom_point()+theme_test()+facet_wrap(~LFA)+ xlab('Year Vessel Manufactured')+ylab('Vessel Breadth')

ggpubr::ggarrange(x+rremove("xlab"),y+rremove("xlab"),z,ncol = 1)
#Demographics on Lic
o = lobster.db('licence_ages')
o = subset(o, LFA %in% 35:38)
o$id = paste(o$LICENCE_ID,o$FIN,o$LFA)
ui = unique(o$id)

ou = list()
m=0
for(i in 1:length(ui)){
      ju = subset(o,id==ui[i])
      nr = nrow(ju)
      for(j in 1:nr){
        m=m+1
      bd = lubridate::year(ju$BIRTHDATE[j])
      sd = lubridate::year(ju$START_DATE[j])
      ed = lubridate::year(ju$END_DATE[j])
      if(ed==4444) ed=2024
      ou[[m]] = data.frame(fishingyrs=sd:ed,birthyear=bd,lfa=ju$LFA[j],lic=ju$LICENCE_ID[j])
      }
}

o = do.call(rbind,ou)
o$Age = o$fishingyrs-o$birthyear
o = subset(o,fishingyrs>1994 &fishingyrs<2024) #missing info prior to 1995
oo = aggregate(Age~lfa+fishingyrs,data=o,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))

o$ageflag = ifelse(is.na(o$Age),0,1)
ol = aggregate(ageflag~lfa+fishingyrs,data=o,FUN=sum)
ol1 = aggregate(ageflag~lfa+fishingyrs,data=o,FUN=length)

ggplot(oo,aes(x=fishingyrs,y=Age[,2],ymin=Age[,1],ymax=Age[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~lfa)+xlab('Fishing Season')+ylab('Age')


####stacking partnerships
o = lobster.db('licence_ages')
o = subset(o, LFA %in% 35:38)
o$id = paste(o$LICENCE_ID,o$FIN,o$LFA)
ui = unique(o$id)

ou = list()
m=0
for(i in 1:length(ui)){
  ju = subset(o,id==ui[i])
  nr = nrow(ju)
  for(j in 1:nr){
    m=m+1
    sd = lubridate::year(ju$START_DATE[j])
    ed = lubridate::year(ju$END_DATE[j])
    if(ed==4444 | ed>2024) ed=2024
    ou[[m]] = data.frame(fishingyrs=sd:ed,lfa=ju$LFA[j],lic=ju$LICENCE_ID[j],type=ju$LIC_SUBTYPE[j])
  }
}

o = do.call(rbind,ou)
o = subset(o,fishingyrs>1998 &fishingyrs<2024) #missing info prior to 1995
o$type = ifelse(o$type %in% c('STACKED','PARTNERSHIP A'),'Stacked-Partnership',o$type)
o=subset(o,type %ni% c('EDUCATIONAL','TEMPORARY'))
z = aggregate(lic~fishingyrs+lfa+type,data=o,FUN=length)
zp = aggregate(lic~fishingyrs+lfa,data=o,FUN=length)
names(zp)[3]='tot'
z = merge(z,zp)
z$Licence_Type=z$lic/z$tot
ggplot(z,aes(fill=type,y=lic,x=fishingyrs))+geom_bar(position = 'stack',stat='identity')+facet_wrap(~lfa)+xlab('Fishing Season')+ylab('Proportiof of Licence Category')
#Logbook Processing
a = lobster.db('process.logs.unfiltered')
a = subset(a, LFA %in% c( 35,36,38) & SYEAR <2024 & SYEAR>2005)



###issues with grids in LFA 35

#centgrid = read.csv(file.path( project.datadirectory("bio.lobster"), "data","maps","lfa27_38_centgrid.csv"))
#grid.key = with(centgrid,paste(LFA,GRID_NUM,sep='.'))
##a35wronggrid = subset(subset(a,LFA==35),!is.na(GRID_NUM)&paste(LFA,GRID_NUM,sep='.')%ni%grid.key)
#a35wronggrid$P = a35wronggrid$WEIGHT_KG/a35wronggrid$NUM_OF_TRAPS
#xx = aggregate(P~SYEAR+LFA,data=a35wronggrid,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))
#ggplot(xx,aes(x=SYEAR,y=P[,2],ymin=P[,1],ymax=P[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('CPUE')

f = lobster.db('season.dates')
b = lobster.db('community_code')
f$days = f$END_DATE-f$START_DATE

sd = aggregate(days~LFA+SYEAR,data=subset(f,LFA %in% 35:38),FUN=sum)
z=ggplot(subset(sd,SYEAR>2004),aes(SYEAR,days))+geom_point()+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Length of Season (days)')


a$DYR = lubridate::decimal_date(a$DATE_FISHED) - lubridate::year(a$DATE_FISHED)
a$WYR = ceiling(a$DYR*52)
a$DWYR = lubridate::year(a$DATE_FISHED) + a$WYR/52
a$P=1

#how many trips
x1 = aggregate(SD_LOG_ID~SYEAR+LICENCE_ID+LFA,data=subset(a,WEIGHT_KG>0 & NUM_OF_TRAPS>0),FUN=function(x) length(unique(x)))
x1$P=x1$SD_LOG_ID
x1 = merge(x1,sd)
x1$Pp = x1$P/as.numeric(x1$days)

 xx = aggregate(P~SYEAR+LFA,data=x1,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))
 xp = aggregate(Pp~SYEAR+LFA,data=x1,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))

 md = aggregate(Pp[,2]~LFA,data=xp,FUN=median)
names(md)[2]='medP'
   
#ggplot(xx,aes(x=SYEAR,y=P[,2],ymin=P[,1],ymax=P[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Number of Trips')
y=ggplot(xp,aes(x=SYEAR,y=Pp[,2],ymin=Pp[,1],ymax=Pp[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Number Trips / Season Length ')+geom_hline(data=md,aes(yintercept=medP),color='red')

ggpubr::ggarrange(z+rremove('xlab'),y,ncol=1)

###distribution of CPUE per trip see if that changes and is the reason for fewer trips
xa = aggregate(cbind(WEIGHT_KG, NUM_OF_TRAPS)~SYEAR+LFA+SD_LOG_ID,data=subset(a,NUM_OF_TRAPS>0 & WEIGHT_KG>0),FUN=sum)
xa$Pc = xa$WEIGHT_KG/xa$NUM_OF_TRAPS
xx = aggregate(Pc~SYEAR+LFA,data=xa,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))
ggplot(xx,aes(x=SYEAR,y=Pc[,2],ymin=Pc[,1],ymax=Pc[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Median CPUE Per Trip')+theme_test()

##trips and cpues

xg = merge(xp,xx)

ggplot(xg,aes(x=Pc[,2],y=Pp[,2],label=SYEAR))+geom_point()+geom_path()+facet_wrap(~LFA)+
  xlab('CPUE')+ylab('Number of Trips / Season Length')+geom_text(data=subset(xg,SYEAR %in% c(2006,2023)),aes(label=SYEAR,x=Pc[,2],y=Pp[,2]))+theme_test()





###vrns moving between LFAs

a = lobster.db('process.logs.unfiltered')
a = subset(a, LFA %in% c(34, 35,36,38,33) & SYEAR <2024 & SYEAR>2005)

xap = aggregate(LFA~SYEAR+VR_NUMBER,data=a,FUN=function(x) length(unique(x)))
xap = subset(xap,LFA>1)
xaa = aggregate(SD_LOG_ID~SYEAR+LFA+VR_NUMBER,data=subset(a,VR_NUMBER %in% xap$VR_NUMBER),FUN=function(x) length(unique(x)))
xaw = pivot_wider(xaa,names_from=LFA,values_from=SD_LOG_ID)
names(xaw)[3:7] = c('l33','l34','l35','l36','l38')

xawr = subset(xaw,!is.na(l35) )
xawr35 = na.zero(xawr)
xawr35$id = ifelse(apply(xawr35[,3:7],1,function(row) sum(row>3))>0,1,0)
xawr35a = aggregate(id~SYEAR,data=xawr35,FUN=sum)
xawr35a$LFA=35
xawr = subset(xaw,!is.na(l36) )
xawr36 = na.zero(xawr)
xawr36$id = ifelse(apply(xawr36[,3:7],1,function(row) sum(row>3))>0,1,0)
xawr36a = aggregate(id~SYEAR,data=xawr36,FUN=sum)
xawr36a$LFA=36
xawr = subset(xaw,!is.na(l38) )
xawr38 = na.zero(xawr)
xawr38$id = ifelse(apply(xawr38[,3:7],1,function(row) sum(row>3))>0,1,0)
xawr38a = aggregate(id~SYEAR,data=xawr38,FUN=sum)
xawr38a$LFA=38

xxx = do.call(rbind,list(xawr35a,xawr36a,xawr38a))
xxx$nlic = ifelse(xxx$LFA==35,94,ifelse(xxx$LFA==36,177,136))
xxx$propLicswitch = xxx$id/xxx$nlic

ggplot(xxx,aes(SYEAR,propLicswitch))+geom_point()+geom_path()+facet_wrap(~LFA)+
  xlab('Fishing Season')+ylab('Proportion of Vessels moving between LFAs within Season')


####################################################################################################################################
####merge vessel info and licence info into logs

av = merge(a,d,by.x=c('VR_NUMBER','SYEAR','LFA'),by.y=c('VR_NUMBER','YR_FISHED','LFA'))
avo = merge(av,o[,c('lic','fishingyrs','lfa','Age')],by.x=c('LICENCE_ID','SYEAR','LFA'),by.y=c('lic','fishingyrs','lfa'))

avo$bCat = ifelse(avo$BREADTH<20,"<20",ifelse(avo$BREADTH>=20&avo$BREADTH<25,"20-24",">25"))
vo = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR+bCat,data=avo,FUN=sum)
vo$CPUE=vo$WEIGHT_KG/vo$NUM_OF_TRAPS

vo$Fishing_Season =vo$SYEAR
vo$Breadth_Category = as.character(vo$bCat)

vol = split(vo,f=list(vo$LFA,vo$bCat))
junk = list()
for(i in 1:length(vol)){
        z = vol[[i]]
        j = which.max(z$CPUE)
        zp = z[j:nrow(z),]
        op = (lm(CPUE~Fishing_Season,data=zp))
        zp$pred = predict(op)
        zp$coef = coef(op)[[2]]
        junk[[i]]=zp
}
jj = do.call(rbind,junk)


ggplot(vo,aes(x=Fishing_Season,y=CPUE,group=Breadth_Category,colour=Breadth_Category))+scale_fill_discrete(labels=c('<20','20-24','>25'))+geom_line(size=1.3)+facet_wrap(~LFA)+
  geom_line(data=jj,aes(x=Fishing_Season,y=pred,group=Breadth_Category,colour=Breadth_Category),linetype='dotted',size=1.3)+xlab('Fishing Season')+ylab('Annual CPUE')




#how many grids per year are they fishing
a = lobster.db('process.logs')
a = subset(a, LFA %in% c(35, 36,38) & SYEAR <2024 & SYEAR>2005)


xg = aggregate(cbind(WEIGHT_KG, NUM_OF_TRAPS)~SYEAR+LICENCE_ID+LFA+GRID_NUM,data=a,FUN=sum)
xgg = aggregate(GRID_NUM~SYEAR+LICENCE_ID+LFA,data=subset(xg,GRID_NUM>0),FUN=length)
xgga = aggregate(GRID_NUM~SYEAR+LFA,data=xgg,FUN=mean)

ggplot(xgga,aes(SYEAR,GRID_NUM))+geom_line()+facet_wrap(~LFA)+theme_test()+ylab('Number of Reported Grids')+xlab('Fishing Season')


xvog = merge(xvo, xgg, by.x=c('SYEAR','LICENCE_ID','LFA', 'VR_NUMBER'),by.y=c('SYEAR','LICENCE_ID','LFA', 'VR_NUMBER'),all.x=T)
xvog$ageBoat = xvog$SYEAR - xvog$YEAR_BUILT

ggplot(subset(xvog, SYEAR==2019 & CPUE<8 & LFA != 28),aes(x=CPUE)) + geom_histogram(aes(y=..density..)) + facet_wrap(~LFA, scales='free_y')

ggplot(subset(xvog, SYEAR==2019),aes(x=BHP)) + geom_histogram() + facet_wrap(~LFA, scales='free_y')

ggplot(subset(xvog, SYEAR==2019),aes(x=WEIGHT_KG)) + geom_histogram() + facet_wrap(~LFA, scales='free_y')
ggplot(subset(xvog, SYEAR==2019),aes(x=LOA)) + geom_histogram() + facet_wrap(~LFA, scales='free_y')

ggplot(subset(xvog, SYEAR==2019),aes(x=ageBoat)) + geom_histogram() + facet_wrap(~LFA, scales='free_y')


ggplot(subset(xvog, SYEAR==2019& GRID_NUM<10),aes(x=GRID_NUM)) + geom_histogram() + facet_wrap(~LFA, scales='free_y')

##Statistics

(aggregate(GRID_NUM~LFA, data=subset(xvog,SYEAR==2019),FUN=function(x) c(mean(x),sd(x),length(x))))
(aggregate(Age~LFA, data=subset(xvog,SYEAR==2019),FUN=function(x) c(mean(x),sd(x),length(x))))
(aggregate(log10(WEIGHT_KG)~LFA, data=subset(xvog,SYEAR==2019),FUN=function(x) c(mean(x),sd(x),length(x))))
(aggregate(P~LFA, data=subset(xvog,SYEAR==2019),FUN=function(x) c(mean(x),sd(x),length(x))))
(aggregate(ageBoat~LFA, data=subset(xvog,SYEAR==2019),FUN=function(x) c(mean(x),sd(x),length(x))))
(aggregate(GROSS_TONNAGE~LFA, data=subset(xvog,SYEAR==2019),FUN=function(x) c(mean(x),sd(x),length(x))))

####value per trip
a = lobster.db('process.logs')
a = subset(a, LFA %in% c(35, 36,38) & SYEAR <2024 & SYEAR>2005)

xg = aggregate(cbind(WEIGHT_KG, NUM_OF_TRAPS)~SYEAR+LICENCE_ID+LFA+GRID_NUM,data=a,FUN=sum)
sp = lobster.db('process_slips')
sp$woy = lubridate::week(sp$Date)

sp = aggregate(PRICE~woy+LFA+SYEAR,data=sp,FUN=median)
ggplot(subset(sp,SYEAR>2006 & LFA %in% c(35,36,38)),aes(woy,PRICE, group=LFA,colour=LFA))+geom_path()+facet_wrap(~SYEAR)+theme_test()


a$woy = lubridate::week(a$DATE_FISHED)

asp = merge(a,sp)
asp$valuePerTrip = asp$WEIGHT_KG*2.20462*(asp$PRICE)
asp$valuePerTrap = asp$WEIGHT_KG*2.20462*(asp$PRICE)/asp$NUM_OF_TRAPS

oo = aggregate(valuePerTrip~LFA+SYEAR,data=asp,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))
oa = aggregate(valuePerTrap~LFA+SYEAR,data=asp,FUN=function(x) quantile(x,probs=c(0.25,.5,.75)))

ggplot(oo,aes(x=SYEAR,y=valuePerTrip[,2],ymin=valuePerTrip[,1],ymax=valuePerTrip[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Value Per Trip')+theme_test(base_size = 14)


ggplot(oa,aes(x=SYEAR,y=valuePerTrap[,2],ymin=valuePerTrap[,1],ymax=valuePerTrap[,3]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Value Per Trap')+theme_test(base_size = 14)




x1 = aggregate(WEIGHT_KG~SYEAR+LFA+WOS,data=subset(asp,WEIGHT_KG>0 & NUM_OF_TRAPS>0),FUN=sum)
x1$P=x1$WEIGHT_KG


xxx = split(x1,f=list(x1$LFA,x1$SYEAR))
junk = list()

for(i in 1:length(xxx)){
  o = xxx[[i]]
  o$prp = cumsum(o$P) / sum(o$P)
  junk[[i]] = o  
}
x1a = do.call(rbind,junk)
x1a$Fishing_Season =x1a$SYEAR
ggplot(subset(x1a),aes(x=WOS,y=prp,group=Fishing_Season,colour=Fishing_Season))+scale_colour_viridis_c(option='inferno')+geom_line()+facet_wrap(~LFA)+xlab('Week of Season')+ylab('Proportion of Total Landings')+theme_test(base_size = 14)


x1 = aggregate(valuePerTrip~SYEAR+LFA+WOS,data=subset(asp,WEIGHT_KG>0 & NUM_OF_TRAPS>0),FUN=sum)
x1$P=x1$valuePerTrip


xxx = split(x1,f=list(x1$LFA,x1$SYEAR))
junk = list()

for(i in 1:length(xxx)){
  o = xxx[[i]]
  o$prp = cumsum(o$P) / sum(o$P)
  junk[[i]] = o  
}
x1a = do.call(rbind,junk)
x1a$Fishing_Season =x1a$SYEAR
ggplot(x1a,aes(x=WOS,y=prp,group=Fishing_Season,colour=Fishing_Season))+scale_colour_viridis_c(option='inferno')+geom_line()+facet_wrap(~LFA)+xlab('Week of Season')+ylab('Proportion of Total Income')+theme_test(base_size = 14)



####timing of trips
a = lobster.db('process.logs')
a = subset(a, LFA %in% c( 35,36,38) & SYEAR <2024 & SYEAR>2005)

x1 = aggregate(SD_LOG_ID~SYEAR+LICENCE_ID+LFA+WOS,data=subset(a,WEIGHT_KG>0 & NUM_OF_TRAPS>0),FUN=function(x) length(unique(x)))
x1$P=x1$SD_LOG_ID

x1a = aggregate(P~SYEAR+LFA+WOS,data=x1,FUN=sum)

xxx = split(x1a,f=list(x1a$LFA,x1a$SYEAR))
junk = list()

for(i in 1:length(xxx)){
  o = xxx[[i]]
  o$prp = cumsum(o$P) / sum(o$P)
  junk[[i]] = o  
}
x1a = do.call(rbind,junk)
x1a$Fishing_Season =x1a$SYEAR
ggplot(x1a,aes(x=WOS,y=prp,group=Fishing_Season,colour=Fishing_Season))+scale_colour_viridis_c(option='inferno')+geom_line()+facet_wrap(~LFA)+xlab('Week of Season')+ylab('Proportion of Total Trips')+theme_test(base_size = 14)


x1aa = aggregate(P~LFA,data=x1a,FUN=max)
names(x1aa)[2]='Pmax'
x1a = merge(x1a,x1aa)
x1a$prop = x1a$P/x1a$Pmax

x1a$YR = as.numeric(x1a$SYEAR)
x1a$YR1 = x1a$YR+x1a$prop




#expand
lfa = c(35,36,38)
ouu = list()
for(i in 1:length(lfa)){
  m = subset(x1a,LFA==lfa[i])
  mi = min(m$WOS)
  mx = max(m$WOS)
  si = min(m$SYEAR)
  sx = max(m$SYEAR)
  da = data.frame(LFA=lfa[i],SYEAR=rep(si:sx,times=length(mi:mx)),WOS=rep(mi:mx,each=length(si:sx)))
  ouu[[i]] = da
}

v=do.call(rbind,ouu)
xx = merge(x1a,v,all.y=T)

ggplot(xx,aes(x=WOS,y=YR1,group=SYEAR))+geom_line()+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Number of Trips')+geom_hline(yintercept=2006:2024)+theme_test()



##cpue and landings
lan = lobster.db('seasonal.landings')
lan$SYEAR=substr(lan$SYEAR,6,9)
require(tidyr)
lan = lan %>% select(SYEAR,LFA35,LFA36,LFA38) %>% pivot_longer(starts_with('LFA'))
lan$LFA = substr(lan$name,4,6)

v = merge(cpu,lan)

ggplot(v,aes(x=unBCPUE,y=value))+geom_point()+geom_path()+facet_wrap(~LFA)
