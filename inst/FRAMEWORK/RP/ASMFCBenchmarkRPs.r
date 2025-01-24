#ASMFC Reference Points 2020

am = read.csv(file.path(project.datadirectory('Framework_LFA35_38'),'ASMFCData','AbundanceTrends.csv'))
rps = data.frame(Stock=c(rep('GOMGBK',2),rep('SNE',2)), RP=rep(c('Abundance_Limit','Abundance_Threshold'), times=2),Level=c(125,89,NA,20))



ggplot(am,aes(x=Year,y=ReferenceAbundance_million))+
  geom_line()+
  facet_wrap(~Stock,scales = 'free_y')+
  theme_test(base_size = 14)+
 # geom_hline(data=subset(rps,RP=='Abundance_Limit'),aes(yintercept = Level),colour='blue',linetype='dashed',linewidth=1.2)+
  geom_hline(data=subset(rps,RP=='Abundance_Threshold'),aes(yintercept = Level),colour='red',linetype='dashed',linewidth=1.2)+
  labs(x='Year',y='Reference Abundance')

  

mx <- am %>%
  group_by(Stock) %>%
  summarise(mean_top3 = mean(sort(ReferenceAbundance_million, decreasing = TRUE)[1:3]))

mxr = merge(mx,rps)

mxr$Prop = mxr$Level/mxr$mean_top3
