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

fd=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS')
setwd(fd)

survey = readRDS( file='summer_fall_survey_data_all_N_Sept24.rds')
survey = st_as_sf(survey)

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs




sf_use_s2(FALSE) #needed for cropping
survey <- suppressWarnings(suppressMessages(
  st_crop(survey,
          c(xmin = -82, ymin = 4539, xmax = 383, ymax = 5200))))


# Project our survey data coordinates:


ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast <- st_transform(ns_coast, crs_utm20)
st_geometry(ns_coast) = st_geometry(ns_coast)
st_crs(ns_coast) <- crs_utm20

survey <- survey %>%   
  st_as_sf()

surv_utm_coords <- st_coordinates(survey)

survey$X1000 <- surv_utm_coords[,1] 
survey$Y1000 <- surv_utm_coords[,2] 
survey = subset(survey,!is.na(N))

survey = subset(survey,YEAR>1999)


ss = as_tibble(survey)
ss$depth.scaled = (ss$z-mean(ss$z))/sd(ss$z)
ss$depth.scaled2 = ss$depth.scaled * ss$depth.scaled
spde <- make_mesh(ss, xy_cols = c("X1000", "Y1000"),
                  n_knots=350,type = "cutoff_search")
plot(spde)

# Add on the barrier mesh component:
mesh <- sdmTMBextra::add_barrier_mesh(
  spde, ns_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)
mesh_df_water <- mesh$mesh_sf[mesh$normal_triangles, ]
mesh_df_land <- mesh$mesh_sf[mesh$barrier_triangles, ]
ggplot(ns_coast) +
  geom_sf() +
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  geom_sf(data = mesh_df_land, size = 1, colour = "green")+
  coord_sf()

data = ss
data$year = data$YEAR
## make new columns to model survey effecst (0 in both columns is ILTS_summer)
data$ILTS_Fall <- ifelse(data$Survey=="ILTS_Fall", 1, 0)
data$MNH_Fall <- ifelse(data$Survey=="MNH_Fall", 1, 0)
data$MNH_NEFSC_Fall <- ifelse(data$Survey=="NEFSC_Fall", 1, 0)
data$RV_Summer <- ifelse(data$Survey=="RV_Summer", 1, 0)


  
  # # Build B-spline basis functions 
  k <- 3
  knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
  bs <- bs(data$depth.scaled, knots = knots, intercept = FALSE)
  bs <- as.data.frame(bs)
  names(bs) <- paste0("bs", names(bs))
  data <- data[, setdiff(names(data), names(data)[grep("^bs[0-9]+", names(data))])]
  data_bs <- cbind(data, bs)
  
  
  m3a <- sdmTMB(
    data = data_bs,
    formula = N ~ 0 + ILTS_Fall + MNH_Fall + MNH_NEFSC_Fall,
    mesh = mesh,
    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
    spatial = "on",
    family = tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1"
  )
  
m3=m3a


##make maps 

model <- m3 ##to say which model to use for maps
model.lab <- "m3" ##to label maps


# Save model results:
saveRDS(model, file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster2000+.rds"))
model = readRDS( file = file.path(project.datadirectory('Framework_LFA35_38'),'outputs','SURVEYS', "summer_fall_Lobster2000+.rds")) #sept 25 300 knots and winner


data$resids <- residuals(model) # randomized quantile residuals
qqnorm(data$resids)
qqline(data$resids)



ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,z<max(data$z) & z>5)
ff$depth.scaled = (ff$z-mean(data$z))/sd(data$z)
k <-  3
knots <-   quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])
basis <- as.data.frame(bs(ff$depth.scaled, knots = knots, intercept = FALSE))
names(basis) <- paste0("bs", names(basis))
ff <- cbind(ff, basis)

yy = unique(data$year)
o = list()
for(i in 1:length(yy)){
  f = ff
  f$year=yy[i]
  o[[i]]=f
}
ff = do.call(rbind,o)

sg = st_coordinates(ff)
ff$X1000 = sg[,1]
ff$Y1000 = sg[,2]
ff$SOURCE='ILTS'
f = as_tibble(ff)
f$MNH <- 0
f$NEFSC <- 0
f$ILTS_Fall <- f$MNH_Fall <- f$MNH_NEFSC_Fall <- 0

#g = predict(model,newdata=f,return_tmb_object = T)

f34 = subset(f,LFA==35)
g34 = predict(model,newdata=f34,return_tmb_object = T)
ind35 = get_index(g34,bias_correct = T)
eao35 = get_eao(g34)
g5 = ggplot(subset(ind35,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000) ')

f34 = subset(f,LFA==36)
g34 = predict(model,newdata=f34,return_tmb_object = T)
ind36 = get_index(g34,bias_correct = T)
eao36 = get_eao(g34)

g6 = ggplot(subset(ind36,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


f3538 = subset(f,LFA %in% 38)
g3538 = predict(model,newdata=f3538,return_tmb_object = T)
ind38 = get_index(g3538,bias_correct = T)
g8 = ggplot(subset(ind38,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


saveRDS(list(ind35,ind36,ind38),'IndicesFromFullComboModelJan282000+.rds')
b = readRDS('IndicesFromFullComboModelOct92000+.rds')
ind35a=b[[1]]
ind36=b[[2]]
ind38=b[[3]]

g5 = ggplot(subset(ind35,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000) ')
g6 = ggplot(subset(ind36,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')
g8 = ggplot(subset(ind38,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')



###everywhere
      g = predict(model,newdata=f,return_tmb_object = T)
      
      v = get_index(g,bias_correct = T)
      #effective area occupied (Thorsen et al 2016)
      voe = get_eao(g)
      names(v)[5:6] = c('logAbund','logAbundse')
   
      ##ggplot
      v1 = subset(v,select=c(year,est,lwr,upr))
      names(v1)[2:4] = c('Abundance','lwrAb','uprAb')
      
      voe1 = subset(voe,select=c(year,est,lwr,upr))
      names(voe1)[2:4] = c('Area','lwrArea','uprArea')
      
      vv = merge(v1,voe1)
      ggplot(subset(vv,year<2024),aes(x=Abundance,y=Area,xmin=lwrAb,xmax=uprAb,ymin=lwrArea,ymax=uprArea))+geom_point(size=2)+geom_errorbar(linetype='dotted')+geom_errorbarh(linetype='dotted')+theme_test(base_size = 14)+geom_path()
      
      
      ################################################################   
      ##can we detect a positive slope between abund and distr (DDHS), there are errors in variables so need to account for that--bayesian makes most sense to me
      df = merge(subset(voe,select=c(year,log_est)),subset(v,select=c(year,logAbund,logAbundse)))
     df$sa = df$logAbund
     df = subset(df,year<2024)
     ##base model
     
     model <- brm(
        formula = bf(log_est ~ sa),  
        data = df,
        family = gaussian(),
        prior = c(
          set_prior("normal(0, 10)", class = "b"),
          set_prior("normal(0, 1)", class = "Intercept")
        ),
        chains = 4,
        iter = 2000
      )
      dr =  as.data.frame(model)
     ##sample across errors in x
     gsa <- function(row) {
       (rnorm(1, mean = row['logAbund'], sd = row['logAbundse']))
     }
     
     junk=list()
     iter=100
     for(i in 1:iter){
       df$sa <- apply(df, 1, gsa)
       um = update(model, newdata=df)
       junk[[i]]  = as.data.frame(um)
     }
     
     drd = bind_rows(junk)
     
     ggplot(drd,aes(b_sa))+
     geom_density() +
       theme_test(base_size = 14)+xlab('slope- log(Area)~log(Abundance)')
      
      h = g$data
      h$pred = model$family$linkinv(h$est)
      gsf = st_as_sf(h)


##features

sf_use_s2(FALSE) #needed for cropping

ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -64, ymax = 47.5))))

ns_coast = st_transform(ns_coast,32620) 
st_geometry(ns_coast) <- st_geometry(ns_coast)/1000
st_crs(ns_coast) <- 32620



rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = rL[rL$LFA %in% c(33:38),]
st_crs(rL) <- 4326
crs_utm20 <- 32620
#rL = rL[-which(!(st_is_valid(rL))),]
rL <- suppressWarnings(suppressMessages(
  st_crop(rL,
          c(xmin = -67.5, ymin = 42, xmax = -62.1, ymax = 46))))
rL <- st_transform(rL, crs_utm20)
st_geometry(rL) <- st_geometry(rL)/1000
st_crs(rL) <- 32620


####maps of residuals
da = st_as_sf(data)
ggplot(subset(da,year<2024)) +
  geom_sf(aes(fill=resids,color=resids),size=1.3) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
    facet_wrap(~year) +
  theme_test_adam()






###prediction surfaces


gsfr = subset(gsf,LFA %in% c(34,35,36,37,38) & Y1000>4850)

#mm = c(0.,max(quantile(gsfr$pred,.999)))
#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
gsfr$CommN = gsfr$pred
gb = ggplot(subset(gsfr, year %in% c(2000:2023))) +
  geom_sf(aes(fill=CommN,color=CommN),size=2.1) + 
  scale_fill_viridis_c(trans='log') +
  scale_color_viridis_c(trans='log') +
  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(st_bbox(gsfr)$ymin,st_bbox(gsfr)$ymax),
           expand = FALSE)




j1 = gb + transition_time(as.integer(year))+labs(title="Year: {frame_time}")
animate(j1,duration=25,height=7,width=7,units="in",res=250)
anim_save(filename = 'ModelM7anim.gif')

###se of predictions 
g = predict(model,newdata=f,se=T,nsim=50)
##gNsim = predict(model,newdata=f,se=F,nsim=3000)
gg = model$family$linkinv(g)

f$se = apply(gg,1,sd)
f$med = apply(gg,1,median)
f$me = apply(gg,1,mean)

gsfer = st_as_sf(f)
gsfer = subset(gsfer,LFA %in% c(34,35,36,37,38)& Y1000 > 4850)

mm = c(0.,max(quantile(gsfer$se,.99)))
#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
gsfer$PredictionError = gsfer$se
ggplot(subset(gsfer, year %in% c(1970:2023))) +
  geom_sf(aes(fill=PredictionError,color=PredictionError),size=2.1) + 
  scale_fill_viridis_c(trans='sqrt',limits=mm) +
  scale_color_viridis_c(trans='sqrt',limits=mm) +
  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(4850,st_bbox(gsfr)$ymax),
           expand = FALSE)


#gsfr$CommB = ifelse(gsfr$pred>mm[2],mm[2],gsfr$pred)
gsfer$CV_Prediction = gsfer$se/gsfer$me
mm = c(0.,max(quantile(gsfer$CV_Prediction,.995)))

ggplot(subset(gsfer, year %in% c(1970:2023))) +
  geom_sf(aes(fill=CV_Prediction,color=CV_Prediction),size=2.1) + 
  scale_fill_viridis_c(limits=mm) +
  scale_color_viridis_c(limits=mm) +
  facet_wrap(~year) +
  geom_sf(data=ns_coast,fill='black')+
  geom_sf(data=rL,colour='black',fill=NA)+
  
  theme_test_adam()+
  coord_sf(xlim = c(st_bbox(ns_coast)$xmin,st_bbox(gsfr)$xmax),
           ylim = c(st_bbox(gsfr)$ymin,st_bbox(gsfr)$ymax),
           expand = FALSE)



##########
#depth splines
year=unique(data$year)
nd <- expand.grid(depth = seq(5, 500, by = 5), year = min(year):max(year))

nd$depth.scaled <- nd$depth - mean(model$data$z)/sd(model$data$z)
nd$survey.ID <- 0
k <- 
knots <- quantile(data$depth.scaled, p = seq(0, 1, len = k)[-c(1,k)])

basis <- as.data.frame(splines::bs(nd$depth.scaled, knots = knots, intercept = FALSE))
names(basis) <- paste0("bs", names(basis))
nd <- cbind(nd, basis)
nd$ILTS_Fall <- nd$MNH_Fall <- nd$MNH_NEFSC_Fall <- nd$RV_Summer <- 0

p <- predict(model, newdata = nd, se_fit = TRUE, re_form = NA)
p$pred = model$family$linkinv(p$est)
p$lower <- p$est - 1.96 * p$est_se
p$upper <- p$est + 1.96 * p$est_se

ggplot(p,aes(x=depth,y=est,group=year,colour=year))+geom_line()+scale_colour_viridis_c()+theme_test(base_size = 14)+labs(x='Depth (m)',y='Response')


#####################
##### terminal year simulations
#kk = subset(ind35,year %in% 1995:1999)
#LRP = LRP35 = mean(kk$est)


lp = readRDS('IndicesFromFullComboModelSept26.rds')
ind35 = lp[[1]]
ind36 = lp[[2]] 
ind38 = lp[[3]]

#38
LRP38 = mean(ind38$est[order(ind38$est)][1])
ggplot(subset(ind38,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')+geom_hline(yintercept = LRP38/1000,col='red')

wb = 'https://cran.r-project.org/src/contrib/Archive/bcp/bcp_3.0.0.tar.gz'
install.packages(wb, repos=NULL, type="source")

#require(bcp)
#ao = ind38[,c('year','est')]
#b = bcp((ao$est),w0 = 0.2, p0 = 0.2)
#dataf = data.frame(y=ao$year,pr = b$posterior.prob,fi=fitted(b)[,1],da=ao$est)

#K = mean(subset(dataf,y %in% 2015:2018,da)[,1])
#LRP = K*.2
#ggplot(dataf,aes(x=y,y=da/1000))+geom_point()+geom_line(data=dataf,aes(x=y,y=fi/1000))+geom_vline(xintercept=subset(dataf,pr>.5,select=y)[,1],colour='red',linetype='dotted',linewidth=1)+theme_test(base_size = 14)+labs(x='Year',y='Abundance')+geom_hline(yintercept=LRP/1000,colour='blue',lwd=1.2)

#36
LRP36 = mean(ind36$est[order(ind36$est)][1])
ggplot(subset(ind36,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')+geom_hline(yintercept = LRP36/1000,col='red')

#require(bcp)
#ao = ind36[which(ind36$year<2024),c('year','est')]
#b = bcp((ao$est),w0 = 0.29, p0 = 0.1)
#dataf = data.frame(y=ao$year,pr = b$posterior.prob,fi=fitted(b)[,1],da=ao$est)

#K = mean(subset(dataf,y %in% 2015:2018,da)[,1])
#LRP = K*.2
#ggplot(dataf,aes(x=y,y=da/1000))+geom_point()+geom_line(data=dataf,aes(x=y,y=fi/1000))+geom_vline(xintercept=subset(dataf,pr>.91,select=y)[,1],colour='red',linetype='dotted',linewidth=1)+theme_test(base_size = 14)+labs(x='Year',y='Abundance')+geom_hline(yintercept=LRP/1000,colour='blue',lwd=1.2)




#35
LRP35 = mean(ind35$est[order(ind35$est)][1])
ggplot(subset(ind35,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')+geom_hline(yintercept = LRP35/1000,col='red')

#ao = ind35[which(ind35$year<2024),c('year','est')]
#b = bcp((ao$est),w0 = 0.29, p0 = 0.1)
#dataf = data.frame(y=ao$year,pr = b$posterior.prob,fi=fitted(b)[,1],da=ao$est)

#K = mean(subset(dataf,y %in% 2014:2018,da)[,1])
#LRP = K*.2
#ggplot(dataf,aes(x=y,y=da/1000))+geom_point()+geom_line(data=dataf,aes(x=y,y=fi/1000))+geom_vline(xintercept=subset(dataf,pr>.8,select=y)[,1],colour='red',linetype='dotted',linewidth=1)+theme_test(base_size = 14)+labs(x='Year',y='Abundance')+geom_hline(yintercept=LRP/1000,colour='blue',lwd=1.2)

#require(bcp)


##presentation plots
#35 

        newD = f_final = subset(f,year %in% 1970:2023 & LFA ==35)
        g = predict(model,newdata=newD,nsim = 1000)
        gg = as.data.frame(model$family$linkinv(g))
        nn = as.data.frame(cbind(gg,yr=newD$year))
        na = aggregate(.~yr,data=nn,FUN=sum)
        ll = data.frame(x=c(0,LRP35/1000),y=c(.25,0))
        op = ou=list()
        for(i in 1:nrow(na)){
          m = t(na[i,-1])[,1]
          g5h = hist(m,breaks=seq(min(floor(m)),max(ceiling(m)),length.out=1000),plot=F)
          x = cumsum(g5h$density)
          x=x/max(x)
          ou[[i]] = data.frame(yr = unique(na[i,'yr']),mids=g5h$mids,x=x)  
        op[[i]] = data.frame(yr = unique(na[i,'yr']), probLLRP=x[which.min(abs(g5h$mids-LRP35))])
            }
        oo = do.call(rbind,ou)
        op = do.call(rbind,op)
        ggplot(op,aes(yr,probLLRP))+geom_point()+theme_test(base_size = 14)+labs(x='Year',y='Proportion of Distribution Below LRP')+geom_hline(yintercept = 0.25,col='red')
        
        j = ggplot(oo,aes(x=mids/1000,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+
          geom_vline(xintercept = LRP35/1000,colour='red',linetype='solid',linewidth=1.3)+
          geom_segment(data=ll,aes(x=x[1],xend=x[2],y=y[1],yend=y[1]),linetype='solid',linewidth=1.3,colour='red')
        
        j1 = j + transition_time(as.integer(yr))+labs(title="Year: {frame_time}")+view_follow(fixed_y = T,fixed_x = c(0,NA))
        animate(j1,duration=25,height=7,width=7,units="in",res=250)
        anim_save(filename = 'LFA35LRPanim.gif')

cc = b5$data
cc$xD = c(0,cc$x[1:(nrow(cc)-1)] -cc$x[2:(nrow(cc))] )*-1
cc$midsR = floor(cc$mids/100)*100+50
cca = aggregate(xD~midsR,data=cc,FUN=sum)
ggplot(cca, aes(x=midsR,y=xD))+geom_bar(stat='identity')+labs(x='Terminal Year Abundance (x1000)',y='Density')+theme_test()+geom_vline(xintercept = LRP35/1000,color='red',linewidth=1.5)
ggplot(cc,aes(x=mids,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+geom_vline(xintercept = LRP35/1000,colour='red',linetype='solid',linewidth=1.3)
ggplot(cc,aes(x=mids/100,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+geom_vline(xintercept = LRP35/1000,colour='red',linetype='solid',linewidth=1.3)

##36
newD = f_final = subset(f,year %in% 1970:2023 & LFA ==36)
g = predict(model,newdata=newD,nsim = 1000)
gg = as.data.frame(model$family$linkinv(g))
nn = as.data.frame(cbind(gg,yr=newD$year))
na = aggregate(.~yr,data=nn,FUN=sum)
ll = data.frame(x=c(0,LRP36/1000),y=c(.25,0))
op = ou=list()
for(i in 1:nrow(na)){
  m = t(na[i,-1])[,1]
  g5h = hist(m,breaks=seq(min(floor(m)),max(ceiling(m)),length.out=1000),plot=F)
  x = cumsum(g5h$density)
  x=x/max(x)
  ou[[i]] = data.frame(yr = unique(na[i,'yr']),mids=g5h$mids,x=x)  
  op[[i]] = data.frame(yr = unique(na[i,'yr']), probLLRP=x[which.min(abs(g5h$mids-LRP36))])
}
oo = do.call(rbind,ou)
op = do.call(rbind,op)
ggplot(op,aes(yr,probLLRP))+geom_point()+theme_test(base_size = 14)+labs(x='Year',y='Proportion of Distribution Below LRP')+geom_hline(yintercept = 0.25,col='red')

j = ggplot(oo,aes(x=mids/1000,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+
  geom_vline(xintercept = LRP36/1000,colour='red',linetype='solid',linewidth=1.3)+
  geom_segment(data=ll,aes(x=x[1],xend=x[2],y=y[1],yend=y[1]),linetype='solid',linewidth=1.3,colour='red')

j1 = j + transition_time(as.integer(yr))+labs(title="Year: {frame_time}")+view_follow(fixed_y = T,fixed_x = c(0,NA))
animate(j1,duration=25,height=7,width=7,units="in",res=250)
anim_save(filename = 'LFA36LRPanim.gif')



####38

newD = f_final = subset(f,year %in% 1970:2023 & LFA ==38)
g = predict(model,newdata=newD,nsim = 1000)
gg = as.data.frame(model$family$linkinv(g))
nn = as.data.frame(cbind(gg,yr=newD$year))
na = aggregate(.~yr,data=nn,FUN=sum)
ll = data.frame(x=c(0,LRP38/1000),y=c(.25,0))
op = ou=list()
for(i in 1:nrow(na)){
  m = t(na[i,-1])[,1]
  g5h = hist(m,breaks=seq(min(floor(m)),max(ceiling(m)),length.out=1000),plot=F)
  x = cumsum(g5h$density)
  x=x/max(x)
  ou[[i]] = data.frame(yr = unique(na[i,'yr']),mids=g5h$mids,x=x)  
  op[[i]] = data.frame(yr = unique(na[i,'yr']), probLLRP=x[which.min(abs(g5h$mids-LRP38))])
}
oo = do.call(rbind,ou)
op = do.call(rbind,op)
ggplot(op,aes(yr,probLLRP))+geom_point()+theme_test(base_size = 14)+labs(x='Year',y='Proportion of Distribution Below LRP')+geom_hline(yintercept = 0.25,col='red')

j = ggplot(oo,aes(x=mids/1000,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+
  geom_vline(xintercept = LRP38/1000,colour='red',linetype='solid',linewidth=1.3)+
  geom_segment(data=ll,aes(x=x[1],xend=x[2],y=y[1],yend=y[1]),linetype='solid',linewidth=1.3,colour='red')

j1 = j + transition_time(as.integer(yr))+labs(title="Year: {frame_time}")+view_follow(fixed_y = T,fixed_x = c(0,NA))
animate(j1,duration=25,height=7,width=7,units="in",res=250)
anim_save(filename = 'LFA38LRPanim.gif')




cc = b5$data
cc$xD = c(0,cc$x[1:(nrow(cc)-1)] -cc$x[2:(nrow(cc))] )*-1
cc$midsR = floor(cc$mids/100)*100+50
cca = aggregate(xD~midsR,data=cc,FUN=sum)
ggplot(cca, aes(x=midsR,y=xD))+geom_bar(stat='identity')+labs(x='Terminal Year Abundance (x1000)',y='Density')+theme_test()+geom_vline(xintercept = LRP35/1000,color='red',linewidth=1.5)
ggplot(cc,aes(x=mids,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+geom_vline(xintercept = LRP35/1000,colour='red',linetype='solid',linewidth=1.3)
ggplot(cc,aes(x=mids/100,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+geom_vline(xintercept = LRP35/1000,colour='red',linetype='solid',linewidth=1.3)




####################
#terminal year LRPs comps

terminalYearPlot = function(model, newD=f_final,nsims=1000,LRP){
  g = predict(model,newdata=newD,nsim = nsims)
  gg = model$family$linkinv(g)
  g5s = colSums(gg)/1000
  g5h = hist(g5s,breaks=seq(min(floor(g5s)),max(ceiling(g5s)),length.out=1000))
  x = cumsum(g5h$density)
  x=x/max(x)
  ests = data.frame(mids=g5h$mids,x=x)  
  v = ggplot(ests,aes(x=mids,y=x))+geom_line(linewidth=1.1)+theme_test(base_size = 14)+labs(x='Abundance',y='Cumulative Distribution ')+geom_vline(xintercept = LRP/1000,colour='red',linetype='dotted',linewidth=1.3)
  return(v)    
}

f_final = subset(f,year==2023 & LFA ==35)
LRP35 = mean(ind35$est[order(ind35$est)][1:10])
b5 = terminalYearPlot(model=model,newD=f_final, LRP=LRP35)




f_final = subset(f,year==2023 & LFA ==36)
LRP36 = mean(ind36$est[order(ind36$est)][1:10])
b6 = terminalYearPlot(model=model,newD=f_final, LRP=LRP35)

f_final = subset(f,year==2023 & LFA ==38)
LRP38 = mean(ind38$est[order(ind38$est)][1:10])
b8 = terminalYearPlot(model=model,newD=f_final, LRP=LRP38)
