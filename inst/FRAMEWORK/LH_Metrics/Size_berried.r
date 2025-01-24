# at sea sampling 
require(bio.lobster)
require(bio.utilities)
require(dplyr)
require(tidyr)
require(ggplot2)
require(devtools)
require(sf)


lobster.db('atSea')
v = subset(atSea,LFA %in% c(35,36,38))
v$UID = paste(v$TRIPNO,v$LFA,v$TRAPNO,v$STRINGNO,sep="_")
v$year = lubridate:::year(v$STARTDATE)
v$mm = lubridate::month(v$STARTDATE)
v$SYEAR = v$year
v$SYEAR = ifelse(v$mm %in% 10:12,v$SYEAR-1,v$SYEAR)
v$Season = ifelse(v$mm %in% c(10,11,12,1),'Fall',ifelse(v$mm %in% 2:4,'Winter',ifelse(v$mm %in% 5:7,'Spring','Outside')))

NTraps = aggregate(UID~LFA+SYEAR,data=v,FUN=function(x) length(unique(x)))
NTrapsSeason = aggregate(UID~LFA+SYEAR+Season,data=v,FUN=function(x) length(unique(x))) #most in fall or spring, separate indices

vF = subset(v,Season=='Fall')

vFLob = subset(v,Season=='Fall' & SPECIESCODE==2550)
ggplot(subset(vFLob,LFA==35),aes(x=CARLENGTH))+geom_histogram()+facet_wrap(~SYEAR,scales='free_y')+geom_vline(xintercept = 82.5,col='red')+theme_test


vfLobA = aggregate(CARLENGTH~LFA+SYEAR,data=vFLob,FUN=function(x) quantile(x,probs=c(0.05,.25,.5,.75,.95)))
ggplot(subset(vfLobA),aes(x=SYEAR,y=CARLENGTH[,3],ymin=CARLENGTH[,1],ymax=CARLENGTH[,5]))+geom_point()+geom_errorbar(width=0)+facet_wrap(~LFA)+xlab('Fishing Season')+ylab('Carapace Length')+theme_test()



##############################
require(bio.lobster)
require(bio.utilities)
require(lubridate)
lobster.db('atSea')
a = subset(atSea,SPECIESCODE==2550)
a = subset(atSea,SPECIESCODE==2550 & SEX %in% 2:3 & CARLENGTH>40 & CARLENGTH<130)
a$yr=year(a$STARTDATE)
a$SYEAR = ifelse(month(a$STARTDATE) %in% c(10,11,12) & a$LFA %in% c(35,36,38),a$yr+1,a$yr)

a$syr = round(a$SYEAR/5)*5
a$cl = round(a$CARLENGTH)

#remove likely errors by LFA
g  = aggregate(TRIPNO~cl+LFA+SYEAR,data=subset(a,SEX==3),FUN=length)
g = split(g,f=list(g$LFA,g$SYEAR))
g = rm.from.list(g)
out = data.frame(LFA=NA,minCL=NA,SYEAR=NA)
for(i in 1:length(g)){
  k = g[[i]]
  l = which(diff(k$cl)==1)[1]
  out[i,] = c(unique(k$LFA),k$cl[l],unique(k$SYEAR))  
}


out = toNums(out,cols=2:3)
i = which(out$LFA==33 & out$minCL<60)
out = out[-i,]
i = which(out$LFA==38 & out$minCL<70)
out = out[-i,]

ggplot(out,aes(x=SYEAR,y=minCL))+geom_point()+geom_smooth(se=F)+facet_wrap(~LFA)



re = aggregate(TRIPNO~LFA+syr+cl,data=subset(a, SEX %in% c(2,3)),FUN=length)
reb = aggregate(TRIPNO~LFA+syr+cl,data=subset(a,SEX==3),FUN=length)
reb$N = reb$TRIPNO
reb$TRIPNO=NULL
rea = merge(re,reb,all=T)
rea = na.zero(rea)
rea$pB = rea$N/rea$TRIPNO

rea$cl = floor(rea$cl/3)*3+1
rea = aggregate(cbind(TRIPNO,N)~LFA+syr+cl,data=rea,FUN=sum)
rea$pB = rea$N/rea$TRIPNO


##nonparametric smoother
# local averaging (cv span selection)
#http://users.stat.umn.edu/~helwig/notes/smooth-notes.html
locavg <- with(subset(rea,LFA==35 & syr==1980), supsmu(cl, pB))
# local regression (gcv span selection)
locreg <- with(subset(rea,LFA==35 & syr==1980), loess.gcv(cl, pB))
# kernel regression (gcv span selection)
kern <- with(subset(rea,LFA==35 & syr==1980), ksmooth.gcv(cl, pB))




with(subset(rea,LFA==35 & syr==1980),plot(cl,pB))
lines(locavg,lwd=2)
lines(locreg,lwd=2,col='red')
lines(kern,lwd=2,col='blue')
legend("topleft", c("supsmu", "loess", "ksmooth"),
       lty = 1:3, lwd = 2, col = c("black", "red", "blue"), bty = "n")



#nls following Lebris using proportions

m1 = formula(pB~1/(1+(exp(a+(b*cl))))) #general
m2 = formula(pB~1/(1+(exp(-a * (cl-b))))) # a is slope b is som50 
m3 = formula(pB~d/(1+(exp(-a * (cl-b))))) # a is slope b is som50 , d = asymptote which is fixed

require(sfsmisc)
lf = c('35','36','38')
li = list()
m=0
for(k in 1:length(lf)){
  Y = unique(subset(rea,LFA==lf[k],select=syr))[,1]
  for(i in 1:length(Y)){
    
    x = subset(rea,LFA==lf[k]&syr==Y[i] & cl<125)
    if(sum(x$TRIPNO)>20){
    locavg <- with(x, supsmu(cl, pB))
    locreg <- with(x, loess.gcv(cl, pB))
    kern <- with(x, ksmooth.gcv(cl, pB))
    #find local maxima to id peak
    
    
    with(x,plot(cl,pB))
    lines(locavg,lwd=2)
    lines(locreg,lwd=2,col='red')
    lines(kern,lwd=2,col='blue')
    legend("topleft", c("supsmu", "loess", "ksmooth"),
           lty = 1:3, lwd = 2, col = c("black", "red", "blue"), bty = "n")
    title(paste('N=',sum(x$N),'SS=',sum(x$TRIPNO)))
    sm = min(x$cl[x$pB>0])
    
    maxs=c()
    maxs[1] =  findMax(kern$x,D1ss(kern$x,kern$y))[1]
    maxs[2] = findMax(locreg$x,D1ss(locreg$x,locreg$y))[1]
    maxs[3] = findMax(locavg$x,D1ss(locavg$x,locavg$y))[1]
    maxs = median(maxs,na.rm=T)
        di= max(x$pB)
        i1=try(nls(m3,data=x,start=list(a=.2,b=(maxs),d=di),lower=list(a=-100,b=30,d=di),upper=list(a=100,b=130,d=di),weights=TRIPNO,algorithm = 'port'),silent=T)  #weighted by sample size
        if(class(i1)=='try-error') next
        m=m+1
        rms = rmse(x$pB,predict(i1))
        ma = mae(x$pB,predict(i1))
        li[[m]] = c(lf[k],Y[i],coef(i1)[1],coef(i1)[2],di,sm,rms,ma,sum(x$TRIPNO),sum(x$N),maxs)
      }
    }
  }
  

som = as.data.frame(do.call(rbind,li))
names(som) = c('LFA','SYR','Slope','Som50','max','minBerr','rmse','mae','ntot','nberr','smooth')
som = toNums(som,2:10)

write.csv(som,file.path(project.datadirectory('bio.lobster'),'data','SoMFromBerried_Oct2024_BOF.csv'))

require(ggplot2)
ggplot(som,aes(x=SYR,y=Som50)) +geom_point() + geom_line()+facet_wrap(~LFA)+theme_test()

ggplot(som,aes(x=SYR,y=minBerr)) +geom_point() + geom_line()+facet_wrap(~LFA)+theme_test()
