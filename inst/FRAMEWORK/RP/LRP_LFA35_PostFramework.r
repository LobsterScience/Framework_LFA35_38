### reference points and relative fishing mortality from final model


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
io = readRDS('IndicesFromFullComboModelOct92000+.rds')
ind35 = io[[1]]

ioF = readRDS('IndicesFromFullComboModelSept26.rds')
ind35F = ioF[[1]]


ggplot(subset(ind35,year>=2000 &year<2024),aes(year,est,ymin=lwr,ymax=upr))+geom_point()+geom_ribbon(alpha=.25)+geom_path()+
    geom_ribbon(data=subset(ind35F,year<2024),aes(year, est, ymin=lwr,ymax=upr),alpha=.25,fill='red')+geom_point(data=subset(ind35F,year<2024),aes(year, est),pch=1)+geom_path(data=subset(ind35F,year<2024),aes(year, est),linetype='dashed')+
  theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')


##converting number to biomass using fall length frequencies from all fall surveys in bay of fundy

yp = ILTS_ITQ_All_Data(redo_base_data = F,size=c(82,300),aggregate=F,species=2550,biomass = F,extend_ts = F)
yp2 = subset(yp, month(SET_DATE)>8)

yp2 = aggregate(SA_CORRECTED_PRORATED_N~FISH_LENGTH,data=yp2,FUN=mean)
names(yp2)[2]='meanden'
yp2$Prop = yp2$meanden/sum(yp2$meanden)
yp2$Weight = lobLW(CL=yp2$FISH_LENGTH,sex = 1)
yp2s = subset(yp2,Prop>0)
mnW = sum(yp2s$Prop*yp2s$Weight)
ggplot(yp2,aes(x=FISH_LENGTH,y=Prop))+geom_bar(stat='identity')+xlab('Carapace Length')+ylab('Density') +theme_test(base_size = 14)

#commercial abundance to weight
ind35$estB = ind35$est*mnW/1000
ind35$lwrB = ind35$lwr*mnW/1000
ind35$uprB = ind35$upr*mnW/1000


ind35F$estB = ind35F$est*mnW/1000
ind35F$lwrB = ind35F$lwr*mnW/1000
ind35F$uprB = ind35F$upr*mnW/1000

##################################################################################################################
####LRP Plots estimates from analyses below
#%%%40% BMSY

mx <- ind35 %>%
  summarise(mean_top3 = mean(sort(estB, decreasing = TRUE)[1:3]))

BMSYLRP = mx*.4/1000
LRPEIV = 44.5
LRPEIVBuff = 57.0

ggplot(subset(ind35,year>=2000 &year<2024),aes(year,estB/1000,ymin=lwrB/1000,ymax=uprB/1000))+geom_point()+geom_ribbon(alpha=.25)+geom_path()+
  geom_ribbon(data=subset(ind35,year<2024),aes(year, estB/1000, ymin=lwrB/1000,ymax=uprB/1000),alpha=.25,fill='red')+geom_point(data=subset(ind35,year<2024),aes(year, estB/1000),pch=1)+geom_path(data=subset(ind35,year<2024),aes(year, estB/1000),linetype='dashed')+
  theme_test(base_size = 14)+labs(x='Year',y='Commercial Biomass') +
  #geom_hline(yintercept = as.numeric(BMSYLRP),col='blue',linewidth=1.3)+
  #geom_hline(yintercept = as.numeric(LRPEIV),col='green',linewidth=1.3)+
  geom_hline(yintercept = as.numeric(LRPEIVBuff),col='purple',linewidth=1.3)



(subset(ind35,year==2023)$estB/1000 -  min(ind35$estB/1000))/ min(ind35$estB/1000)
(subset(ind35,year==2023)$estB/1000 -  LRPEIV)/ LRPEIV
(subset(ind35,year==2023)$estB/1000 -  LRPEIVBuff)/ LRPEIVBuff
###total landings and landings by vessel

g = lobster.db('landings_by_vessel')
g$year = g$SYEAR
g$LFA = ifelse(g$LFA=='38B',38,g$LFA)
glo = aggregate(MT~LFA+year,data=g,FUN=function(x) quantile(x,0.25)); names(glo)[3]='lCI'
gup = aggregate(MT~LFA+year,data=g,FUN=function(x) quantile(x,0.75)); names(gup)[3]='uCI'
gme = aggregate(MT~LFA+year,data=g,FUN=mean); names(gme)[3]='MEAN_T'

g = merge(merge(glo,gup),gme)

g = subset(g,LFA==35)

#all lfa plots
ggplot(g,aes(x=year,y=MEAN_T,ymin=lCI,ymax=uCI))+geom_point()+geom_line()+geom_errorbar(alpha=.5, width=0)+facet_wrap(~LFA)+theme_test(base_size=14)+labs(y='Landings Per Vessel', x='Year')


gl = lobster.db('seasonal.landings')
gl$year= as.numeric(substring(gl$SYEAR,6,9))
gl = subset(gl,select=c(year,LFA35))

gg = merge(gl,g)
maxl =max(gg$LFA35)
maxv = max(gg$MEAN_T)
ggplot(subset(gg, year<2024), aes(x = year)) +
  geom_bar(aes(y =LFA35), stat = "identity", width = 0.75,alpha=.25) +
  geom_point(data = subset(gg,year<2024), aes(y = MEAN_T * (maxl /maxv)),size=3)+
  geom_line(data = subset(gg,year<2024), aes(y = MEAN_T * (maxl /maxv)),linewidth=1.1) +
  labs(x = "Season Year", y = "Landings (t)") +
    scale_y_continuous(
    limits = c(0, 4500),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * (maxv / maxl), name = "Mean Landings Per Vessel (t)")) + 
    scale_x_continuous(breaks = seq(1990, 2023, by = 5), expand = c(0.01, 0)) + 
  guides(fill = "none", color = "none") 



#############relationship between mean landings per vessel and Commercial Biomass
ind35$oyr = ind35$year

ind35$year=ind35$oyr+1

g5 = merge(gg,ind35)
g5$est1000 = g5$estB/1000
with(g5,plot(est1000,MEAN_T))


#######################################################################
#########################MLE#############################################
            a  <- glm(MEAN_T~I(1/(est1000)),data=g5,family=gaussian(link='inverse'))
            ac = coef(a)
            dn <- data.frame(est1000=c(0,g5[order(g5$est1000),'est1000']))
            dn$Preds  =predict(a,newdata = dn,type = 'response')
            plot(MEAN_T~est1000,data=g5)
            lines(dn$est1000,dn$Preds)
            
            g3 <- ggplot(g5, aes(est1000, MEAN_T)) + geom_point()+theme_test()
            g4 = g3 + geom_smooth(method = "glm", formula = y ~ I(1/x), method.args=list(family = gaussian(link = "inverse")) , fill = "blue", alpha = 0.2,fullrange=T) +labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)
            
            
            
            mL = mean(subset(g,year %in% 1990:1995)$MEAN_T)
            LRP_D <- data.frame(est1000=seq(0,300,by=.2))
            LRP_D$Preds  =predict(a,newdata = LRP_D,type = 'response')
            
            LRP = LRP_D$est1000[which.min(abs(LRP_D$Preds-mL))]
            
            g4+geom_vline(xintercept=LRP,colour='red')+geom_hline(yintercept=mL,colour='black',linetype='dashed')
            
            ###################
            #### liklihood surface for inverse gaussian
            
                            log_likelihood <- function(beta0, beta1, X,Y) {
                              # Compute linear predictor and mean (mu)
                              eta <- beta0 + beta1 * 1/X
                              mu <- 1 / eta
                              
                              # Compute residuals and variance
                              residuals <- Y - mu
                              sigma2 <- var(residuals)
                              
                              # Calculate log-likelihood
                              n <- length(Y)
                              log_likelihood <- -n / 2 * log(2 * pi) - n / 2 * log(sigma2) - sum(residuals^2) / (2 * sigma2)
                              
                              return(log_likelihood)
                            }
                            beta0_seq <- seq(0.003, .2, length.out = 300)  # Avoid 0 to prevent division by zero
                            beta1_seq <- seq(-20,100, length.out = 600)
                            likelihood_matrix <- matrix(NA, nrow = length(beta0_seq), ncol = length(beta1_seq))
                            
                            # Compute log-likelihood for each combination of parameters
                            for (i in 1:length(beta0_seq)) {
                              for (j in 1:length(beta1_seq)) {
                                likelihood_matrix[i, j] <- log_likelihood(beta0_seq[i], beta1_seq[j], X=g5$est1000,Y=g5$MEAN_T)
                              }
                            }
                            
                            # Convert the matrix to a data frame for plotting
                            likelihood_df <- expand.grid(Beta0 = beta0_seq, Beta1 = beta1_seq)
                            likelihood_df$LogLikelihood <- as.vector(likelihood_matrix)
                            
                            mis = which(likelihood_df$LogLikelihood<quantile(likelihood_df$LogLikelihood,0.05))
                            
                            likelihood_df1 = likelihood_df
                            likelihood_df1$LogLikelihood[mis] <- NA
                            ggplot(likelihood_df1, aes(x = Beta0, y = Beta1, z = LogLikelihood)) +
                                geom_contour_filled(aes(color = ..level..), bins = 100) +
                              geom_contour(colour='white')+
                              labs(x = "Intercept (Beta0)", y = "Slope (Beta1)", color = "Log-Likelihood") +
                              theme(legend.position='none')+
                            geom_point(aes(x=ac[1],y=ac[2]))
            
                            
            ############ what is the distribution around the LRP
                      ###sample from the joint precision matrix
                            require(MASS)
                            # Get the coefficients and the design matrix
                            beta_hat <- coef(a)
                            X_matrix <- model.matrix(a)
                            
                            # Calculate the Hessian (negative of the second derivative of the log-likelihood)
                            # Use the `vcov` function to get the covariance matrix
                            vcov_matrix <- vcov(a)
                            
                            # The joint precision matrix is the inverse of the covariance matrix
                            precision_matrix <- solve(vcov_matrix)
                            
                            # Sample from the multivariate normal distribution
                            # Number of samples to draw
                            n_samples <- 1000
                            
                            # Generate samples from the multivariate normal using the mean (beta_hat) and the precision matrix
                            samples <- mvrnorm(n_samples, mu = beta_hat, Sigma = solve(precision_matrix))
                            
                            invG = function(x,b0,b1){
                              1/(b0+(b1*1/x))
                            }
                            ests = seq(10,80,by=.1)
                          #distribution of LRP  
                            LRP_jp=c()
                            for(i in 1 :nrow(samples)) {
                              y = invG(x=ests,b1=samples[i,2],b0=samples[i,1])
                              LRP_jp = c(LRP_jp, ests[which.min(abs(y-mL))])
                            }     
                            
                            ggplot(as.data.frame(LRP_jp),aes(LRP_jp))+geom_histogram()+labs(x='LRP',y='Frequency')+geom_vline(xintercept = LRP,col='red',linewidth=1.4)
            
                            
                            g4+geom_vline(xintercept=LRP_jp,colour='red',alpha=.01)+geom_hline(yintercept=mL,colour='black',linetype='dashed')+geom_vline(xintercept=LRP,colour='black')
                
###########################################################
########## bayesian to explore the errors in variables
                            require(brms)
                            
                            g5$sE = g5$est1000
                            g5 = subset(g5,year<2024)
                            model <- brm(MEAN_T~I(1/(sE)),data=g5,family=gaussian(link='inverse'),
                                          prior = c(
                                            set_prior("normal(15, 6)", class = "b"),
                                            set_prior("normal(0, 10)", class = "Intercept")
                                          ),
                                          chains = 4,
                                          iter = 2000, cores=4
                            )
                            
                            ###dist of params
                            ggplot(as.data.frame(model),aes(b_Intercept,b_I1DsE))+
                              geom_hex(bins = 70) +
                              scale_fill_continuous(type = "viridis") +
                              scale_x_continuous(limits=c(0.015,.04))+
                            theme_bw()
                            mod = as.data.frame(model)
                            ###lines
                            g3 <- ggplot(g5, aes(est1000, MEAN_T,ymin=lCI,ymax=uCI,xmin=lwrB/1000,xmax=uprB/1000)) + geom_point(size=2)+geom_errorbar(data=g5, aes(est1000, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test(base_size = 14)+labs(x='Survey Biomass',y='Landings per Vessel')
                            
                            
                            xs = data.frame(sE=seq(0,max(g5$est1000),length.out=20))
                            xp = as.data.frame(model)
                            
                            
                            jnk = list()
                            for(i in 1:nrow(xp)){
                              xs$Pred = 1/(xp[i,'b_Intercept']+xp[i,'b_I1DsE']/xs$sE)
                              xs$i = i
                              jnk[[i]] = xs
                            }
                            
                            ###LRP from base model
                            ests = seq(0,300,by=.1)
                            #distribution of LRP  
                            LRP=c()
                            for(i in 1 :nrow(xp)) {
                              y = 1/(xp[i,'b_Intercept']+xp[i,'b_I1DsE']/ests)
                              LRP = c(LRP, ests[which.min(abs(y-mL))])
                            }     
                            
                            
                            kk = bind_rows(jnk)
                            bk = aggregate(Pred~sE,data=kk, FUN=function(x) quantile(x,c( 0.025,0.975)))
                            
                            
                            ggplot(g5, aes(est1000, MEAN_T)) + geom_point(size=2)+geom_errorbar(data=g5, aes(est1000, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
                              geom_line(data=kk,aes(x=sE,y=Pred,group=i),colour='blue',alpha=.01) +labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)
                            
                            gsa <- function(row) {
                              mm = findMoments(as.numeric(row['lwrB'])/1000,as.numeric(row['uprB'])/1000,dist='lnorm')
                              rlnorm(1, mean = mm$xbar, sd = sqrt(mm$sig2))
                            }
                            
                            junk=list()
                            iter=100
                            for(i in 1:iter){
                              g5$sE <- apply(g5, 1, gsa)
                              um = update(model, newdata=g5)
                              junk[[i]]  = as.data.frame(um)
                            }
                            
                            dd = bind_rows(junk)
                            ggplot(dd,aes(b_Intercept,b_I1DsE))+
                              geom_hex(bins = 70) +
                              scale_fill_continuous(type = "viridis") +
                              scale_x_continuous(limits=c(0.017,.058))+
                              theme_test(base_size=14)+ labs(x='Intercept',y='Slope')
                            
                            
                            jnk = list()
                            for(i in 1:nrow(dd)){
                              xs$Pred = 1/(dd[i,'b_Intercept']+dd[i,'b_I1DsE']/xs$sE)
                              xs$i = i
                              jnk[[i]] = xs
                            }
                            
                            tos = bind_rows(jnk)
                            
                            
                            
                            ####lrp
                            
                            mL = mean(subset(g,year %in% 1990:1995)$MEAN_T)
                            LRP_D <- data.frame(est1000=seq(0,300,by=.2))
                            
                            invG = function(x,b0,b1){
                              1/(b0+(b1*1/x))
                            }
                            ests = seq(-400,3000,by=.1)
                            #distribution of LRP  
                            LRP_jp=c()
                            LRP_buff = c()
                            #ddthin = dd[seq(1,nrow(dd),by=50),]
                            for(i in 1 :nrow(dd)) {
                              y = 1/(dd[i,'b_Intercept']+dd[i,'b_I1DsE']/ests)
                              LRP_jp = c(LRP_jp, ests[which.min(abs(y-mL))])
                              LRP_buff = c(LRP_buff, ests[which.min(abs(y-mL*1.25))])
                            }     
                            
                            ggplot(as.data.frame(LRP_jp),aes(LRP_jp))+geom_density()+geom_density(data=data.frame(LRP_buff),aes(LRP_buff),colour='red')+
                              geom_density(data=data.frame(LRP),aes(LRP),colour='blue')+geom_vline(xintercept = c(mean(LRP_jp),mean(LRP_buff),mean(LRP)),colour=c('black','red','blue'),linewidth=1, linetype='dotted')+
                              labs(x='LRP',y='Frequency') +theme_test(base_size = 14)
                            
                            
                            b = aggregate(Pred~sE,data=tos, FUN=function(x) quantile(x,c( 0.025,0.975)))
                            
                            
                            #EIV
                            ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbar(data=g5, aes(est1000,y=MEAN_T, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
                              geom_ribbon(data=b, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
                              labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)+
                              geom_vline(xintercept=LRP_jp,colour='grey',alpha=.01)+geom_hline(yintercept=mL,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP_jp),colour='black')+
                              scale_x_continuous(limits= c(0,4500))
                            
                            #Buff
                            ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbar(data=g5, aes(est1000,y=MEAN_T, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
                              geom_ribbon(data=b, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
                              labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)+
                              geom_vline(xintercept=LRP_buff,colour='red',alpha=.005)+geom_hline(yintercept=mL*1.25,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP_buff),colour='black')
                            
                            #Buff
                            ggplot() + geom_point(data=g5, aes(x=est1000,y= MEAN_T),size=2)+geom_errorbar(data=g5, aes(est1000,y=MEAN_T, ymin=lCI,ymax=uCI),linetype='dotted')+geom_errorbarh(data=g5, aes(y=MEAN_T,xmin= lwrB/1000,xmax=uprB/1000),linetype='dotted')+theme_test()+
                              geom_ribbon(data=bk, aes(x=sE, ymin=Pred[,1],ymax=Pred[,2]),fill='blue', colour='blue', alpha=0.2) +
                              labs(x='Survey Commercial Biomass',y='Mean Landings per Vessel')+theme_test(base_size = 14)+
                              geom_vline(xintercept=LRP,colour='blue',alpha=.005)+geom_hline(yintercept=mL,colour='black',linetype='dashed')+geom_vline(xintercept=mean(LRP),colour='black')
                            
      
mean(LRP); quantile(LRP, c(0.025,0.975))      
mean(LRP_jp); quantile(LRP_jp, c(0.025,0.975))     
mean(LRP_buff); quantile(LRP_buff, c(0.025,0.975))  
#####################################################################################################################################################################
####################################################
ggplot(g5,aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')
ggplot(g5,aes(year,LFA35))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)

#tremblay 2012
ggplot(g5,aes(year,LFA35))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)+geom_hline(yintercept =292,colour='red',linewidth=2 )

#if we used Bmsy proxies
k = mean(ind35$est[order(ind35$est,decreasing = T)][1:10])/1000
ggplot(subset(ind35,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.2,colour='red',linewidth=2)




ggplot(subset(g5,year<2024 & year>=1995),aes(year,y=LFA35/(estB/1000)))+geom_point()+
  geom_line()+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality')

ind36$year=ind36$year+1
g6 = merge(g[,c('year','LFA36')],ind36)
g6$relf = g6$LFA36/g6$estB*1000
ggplot(subset(g6,year<2024 & year>1995),aes(year,y=LFA36/(estB/1000)))+geom_point()+
  geom_line()+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality')

ggplot(g6,aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')
ggplot(g6,aes(year,LFA36))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)

#tremblay 2012
ggplot(g6,aes(year,LFA36))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)+geom_hline(yintercept =266,colour='red',linewidth=2 )

#if we used Bmsy proxies
k = mean(ind36$est[order(ind36$est,decreasing = T)][1:10])/1000
ggplot(subset(ind36,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.2,colour='red',linewidth=2)




ind38$year=ind38$year+1

g8 = merge(g[,c('year','LFA38')],ind38)
g8$relf = g8$LFA38/g8$estB*1000

ggplot(g8,aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')
ggplot(g8,aes(year,LFA38))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)

#tremblay 2012
ggplot(g8,aes(year,LFA38))+geom_bar(stat='identity')+labs(x='Year',y='Landings')+theme_test(base_size=14)+geom_hline(yintercept =259,colour='red',linewidth=2 )

#if we used Bmsy proxies
k = mean(ind38$est[order(ind38$est,decreasing = T)][1:10])/1000
ggplot(subset(ind38,year<2024),aes(year,est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance')+
  geom_hline(yintercept = k*.2,colour='red',linewidth=2)


ggplot(subset(g8,year<2024 & year>1995),aes(year,y=LFA38/(estB/1000)))+geom_point()+
  geom_line()+
  theme_test(base_size = 14)+labs(x='Year',y='Relative Fishing Mortality')


################################################
################cpue to survey biomass hyperstablity plots

vb = readRDS(file=file.path(project.datadirectory('Framework_LFA35_38'),'outputs','CPUE','unBIASED_CPUE.rds'))
vb = vb[[1]]

vb$year=vb$SYEAR
vb5 = merge(subset(vb,LFA==35),g5)
vb5$estBred =  vb5$estB/1000000

ggplot(vb5,aes(x=estBred,y=unBCPUE))+
  geom_path()+
  geom_point()+
  geom_text(data=subset(vb5,year %in% c(2006,2023)),aes(label=year,x=estBred,y=unBCPUE),position=position_nudge(y=c(-.1,.1)))+
  geom_smooth(method = "nls",formula = y ~ a*x^b,se=FALSE, method.args = list(start = list(a = -1, b = .9)),   data = vb5)+theme_test(base_size=14)+labs(x='Survey Commercial Biomass',y='CPUE')

bb = nls(unBCPUE~a*estBred^b,start=(list(a=1,b=.9)),data=vb5)
bo5 = nlsBooter(bb,vb5,Qlow1 = .025,Qhigh1 = 0.975)
b05a = data.frame(b=bo5$coefboot[,2],LFA=35)

vb$year=vb$SYEAR
vb6 = merge(subset(vb,LFA==36),g6)
vb6$estBred =  vb6$estB/1000000

ggplot(vb6,aes(x=estBred,y=unBCPUE))+
  geom_path()+
  geom_point()+
  geom_text(data=subset(vb6,year %in% c(2006,2023)),aes(label=year,x=estBred,y=unBCPUE),position=position_nudge(x=c(0,-.09),y=c(-.05,-0.01)))+
  geom_smooth(method = "nls",formula = y ~ a*x^b,se=FALSE, method.args = list(start = list(a = -1, b = .9)),   data = vb6)+theme_test(base_size=14)+labs(x='Survey Commercial Biomass',y='CPUE')

b6 = nls(unBCPUE~a*estBred^b,start=(list(a=1,b=.9)),data=vb6)
bo6 = nlsBooter(b6,vb6,Qlow1 = .025,Qhigh1 = 0.975)
b06a = data.frame(b=bo6$coefboot[,2],LFA=36)



vb$year=vb$SYEAR
vb8 = merge(subset(vb,LFA==38),g8)
vb8$estBred =  vb8$estB/1000000

ggplot(vb8,aes(x=estBred,y=unBCPUE))+
  geom_path()+
  geom_point()+
  geom_text(data=subset(vb8,year %in% c(2006,2023)),aes(label=year,x=estBred,y=unBCPUE),position=position_nudge(x=c(0,.02),y=c(-.05,-0.01)))+
  geom_smooth(method = "nls",formula = y ~ a*x^b,se=FALSE, method.args = list(start = list(a = -1, b = .9)),   data = vb8)+theme_test(base_size=14)+labs(x='Survey Commercial Biomass',y='CPUE')

b8 = nls(unBCPUE~a*estBred^b,start=(list(a=1,b=.9)),data=vb8)

bo8 = nlsBooter(b8,vb8,Qlow1 = .025,Qhigh1 = 0.975)
b08a = data.frame(b=bo8$coefboot[,2],LFA=38)

pars = do.call(rbind,list(b05a,b06a,b08a))

ggplot(pars,aes(x=b))+geom_bar(stat='density')+facet_wrap(~LFA)+theme_test(base_size = 14)+labs(x='Shape Parameter',y='Density')

########
####example hyperstability plot

x = seq(0,1,by=.01)
df = data.frame(x)
ggplot(df,aes(x))+
  stat_function(fun=function(x) x^.25)+
  stat_function(fun=function(x) x^.5)+
  stat_function(fun=function(x) x^1)+
  stat_function(fun=function(x) x^2)+
  stat_function(fun=function(x) x^4)+
  theme_test(base_size = 14,axis.text.x=element_blank())+
  labs(x='Biomass',y='CPUE')+
  annotate("text",x=.25,y=.9,label='Hyperstability')+
  annotate("text",x=.8,y=.1,label='Hyperdepletion')+
  annotate("text",x=.13,y=.67,label=expression(paste(beta,"= 0.25",sep=" ")))+
  annotate("text",x=.23,y=.53,label=expression(paste(beta,"= 0.5",sep=" ")))+
  annotate("text",x=.33,y=.39,label=expression(paste(beta,"= 1",sep=" ")))+
  annotate("text",x=.43,y=.25,label=expression(paste(beta,"= 2",sep=" ")))+
  annotate("text",x=.53,y=.11,label=expression(paste(beta,"= 4",sep=" ")))
  
