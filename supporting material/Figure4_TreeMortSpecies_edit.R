##########################################################################
#Figure 4: Species-Level Tree Mortality
#########################################################################
#Tori Meakem
#10/26/2016

#To run code:
#############
#Load R files: scbi.stem1.rdata, scbi.stem2.rdata, allmort.rdata
setwd("C:/Users/MeakemV/Dropbox (Smithsonian)/Tree Mortality/Data/Data Analysis_Tori")
load("C:/Users/MeakemV/Dropbox (Smithsonian)/Tree Mortality/Data/Data Analysis_Tori/scbi.stem1.rdata")
load("C:/Users/MeakemV/Dropbox (Smithsonian)/Tree Mortality/Data/Data Analysis_Tori/scbi.stem2.rdata")
load("C:/Users/MeakemV/Dropbox (Smithsonian)/Tree Mortality/Data/Data Analysis_Tori/allmort.rdata")

#Highlight and run all functions in the "Function List" (below)

#Run lines directly below
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(kSamples)
library(grid)
runplot()
dev.off()

#############################################
runplot=function(){
hectar=25.6 #Plot size in hectares

#Identify dead trees, apply filters
mind12<-indcalcmort12(scbi.stem1,scbi.stem2,100,rounddown=F)
mind23<-indcalcmort23(allmort)
mind34<-indcalcmort34(allmort)


#Only select species with more than 100 individuals
tab12<-table(mind12$species);tab12<-data.frame(tab12[tab12>=100]);tab12$Species<-rownames(tab12)
tab23<-table(mind23$species);tab23<-data.frame(tab23[tab23>=100]);tab23$Species<-rownames(tab23)
tab34<-table(mind34$species);tab34<-data.frame(tab34[tab34>=100]);tab34$Species<-rownames(tab34)
mind12<-mind12[mind12$species%in%tab12$Species,]
mind23<-mind23[mind23$species%in%tab23$Species,]
mind34<-mind34[mind34$species%in%tab34$Species,]

#Calculate mortality rate
##################################
mortrate12<-mszmn(mind12$mort,mind12$species,mind12$timeint)
mortrate23<-mszmn(mind23$mort,mind23$species,mind23$timeint)
mortrate34<-mszmn(mind34$mort,mind34$species,mind34$timeint)
mortrate12$Census<-"2008-2013";mortrate23$Census<-"2013-2014";mortrate34$Census<-"2014-2015"
speciesmort<-rbind(mortrate12,mortrate23,mortrate34)
speciesmort<-speciesmort[!is.na(speciesmort$mortrate)&speciesmort$mortrate!=Inf,]
row.names(speciesmort)<-seq(nrow(speciesmort))

#Calculate CIS: Normal Approximation
species.morts<-speciesmort[speciesmort$ndead>5,]
species.morts$mrtrate<-100*(1-((species.morts$ninit-species.morts$ndead)/species.morts$ninit))
#mortrate12<-100*(1-((speciesmort$ninit-speciesmort$ndead)/speciesmort$ninit)^(1/speciesmort$timeint))
ntot<-nrow(species.morts)
for (i in 1:ntot){
  prop<-prop.test(species.morts$ndead[i],species.morts$ninit[i],species.morts$mrtrate[i]/100)
  ci1<-prop$conf.int[1]*species.morts$ninit[i]
  ci2<-prop$conf.int[2]*species.morts$ninit[i]
  ci.lo<-100*(1-((species.morts$ninit[i]-ci1)/species.morts$ninit[i])^(1/species.morts$timeint[i]))
  ci.hi<-100*(1-((species.morts$ninit[i]-ci2)/species.morts$ninit[i])^(1/species.morts$timeint[i]))
  prop=data.frame(Species=species.morts$Species[i],Census=species.morts$Census[i],ci.lo=ci.lo,ci.hi=ci.hi)
  if (i==1)
    props=prop
  else
    props=rbind(props,prop)}
#If ndead <= 5: Exact Binomial Test
species_morts<-speciesmort[speciesmort$ndead<=5&speciesmort$ndead>0,]
species_morts$mrtrate<-100*(1-((species_morts$ninit-species_morts$ndead)/species_morts$ninit))
ntots<-nrow(species_morts)
for (i in 1:ntots){
binom<-binom.test(species_morts$ndead[i],species_morts$ninit[i],species_morts$mrtrate[i]/100)
ci1<-binom$conf.int[1]*species_morts$ninit[i]
ci2<-binom$conf.int[2]*species_morts$ninit[i]
ci.lo<-100*(1-((species_morts$ninit[i]-ci1)/species_morts$ninit[i])^(1/species_morts$timeint[i]))
ci.hi<-100*(1-((species_morts$ninit[i]-ci2)/species_morts$ninit[i])^(1/species_morts$timeint[i]))
Binom=data.frame(Species=species_morts$Species[i],Census=species_morts$Census[i],ci.lo=ci.lo,ci.hi=ci.hi)
if (i==1)
  binoms=Binom
else
  binoms=rbind(binoms,Binom)}
allcis<-rbind(props,binoms)
speciesmort<-merge(speciesmort,allcis,all.x=T,all.y=T)

#Formatting; fixes bar width for species that lack records for all 3 years (ex: tiam)
cen12<-speciesmort[speciesmort$Census=="2008-2013",]
cen23<-speciesmort[speciesmort$Census=="2013-2014",]
cen34<-speciesmort[speciesmort$Census=="2014-2015",]
cen12$Census<-NULL;cen23$Census<-NULL;cen34$Census<-NULL
cen12$ninit<-NULL;cen23$ninit<-NULL;cen34$ninit<-NULL
cen12$ndead<-NULL;cen23$ndead<-NULL;cen34$ndead<-NULL
cen12$timeint<-NULL;cen23$timeint<-NULL;cen34$timeint<-NULL
colnames(cen12)<-c("Species","mortrate12","ci.lo12","ci.hi12")
colnames(cen23)<-c("Species","mortrate23","ci.lo23","ci.hi23")
colnames(cen34)<-c("Species","mortrate34","ci.lo34","ci.hi34")

speciesmorts<-merge(cen12,cen23,all.x=T,all.y=T)
speciesmorts<-merge(speciesmorts,cen34,all.x=T,all.y=T)
morts12<-data.frame(mortrate=speciesmorts$mortrate12,ci.lo=speciesmorts$ci.lo12,ci.hi=speciesmorts$ci.hi12,Species=speciesmorts$Species,Census="2008-2013")
morts23<-data.frame(mortrate=speciesmorts$mortrate23,ci.lo=speciesmorts$ci.lo23,ci.hi=speciesmorts$ci.hi23,Species=speciesmorts$Species,Census="2013-2014")
morts34<-data.frame(mortrate=speciesmorts$mortrate34,ci.lo=speciesmorts$ci.lo34,ci.hi=speciesmorts$ci.hi34,Species=speciesmorts$Species,Census="2014-2015")
speciesmort<-rbind(morts12,morts23,morts34)
speciesmort <- with(speciesmort, speciesmort[order(Census, -as.numeric(mortrate)), ])
speciesmort$Species <- factor(speciesmort$Species, levels=speciesmort$Species) #Warning: levels deprecated - 
dodge <- position_dodge(width=.75)
##################################

#Plot mortality rate by species
c1<-ggplot(data=speciesmort,aes(x=Species,y=speciesmort$mortrate)) + geom_bar(stat="identity",position=dodge,aes(fill=Census))+ geom_errorbar(aes(ymin=ci.lo,ymax=ci.hi,fill=Census),width=0.75,size=.2,position=dodge) + labs(x = expression("Species"), y = expression(Mortality~Rate~("%"~y^{-1}))) + scale_fill_manual(values=c("darkblue","forestgreen","darkorange"),name="")+scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 45,size=10,hjust=1,vjust=1),text=element_text(size=12),axis.text.y=element_text(size=10),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.position=c(.75,.75),legend.key.size=unit(.4,"cm")) #+annotate("text",x=1,y=30,label="a)",size=8) 


#Biomass Mortality
####################################
#subset dead trees
mortagb12<-mind12[mind12$mort==1&!is.na(mind12$mort),]
mortagb23<-mind23[mind23$mort==1&!is.na(mind23$mort),]
mortagb34<-mind34[mind34$mort==1&!is.na(mind34$mort),]
#Calculate M per species
biomassmort12=eszmn(mortagb12$agb,mortagb12$species,hectar)
biomassmort23=eszmn(mortagb23$agb,mortagb23$species,hectar)
biomassmort34=eszmn(mortagb34$agb,mortagb34$species,hectar)
biomassmort12$Census<-"2008-2013"
biomassmort23$Census<-"2013-2014"
biomassmort34$Census<-"2014-2015"
speciesagbmort<-rbind(biomassmort12,biomassmort23,biomassmort34)
speciesagbmort<-speciesagbmort[!is.na(speciesagbmort$agb)&speciesagbmort$agb!=Inf,]

#Reformat data for the graph (ggplot)
cens12<-speciesagbmort[speciesagbmort$Census=="2008-2013",]
cens23<-speciesagbmort[speciesagbmort$Census=="2013-2014",]
cens34<-speciesagbmort[speciesagbmort$Census=="2014-2015",]
colnames(cens12)[1]<-"agb12"
colnames(cens23)[1]<-"agb23"
colnames(cens34)[1]<-"agb34"
cens12$Census<-NULL;cens23$Census<-NULL;cens34$Census<-NULL
speciesagbmorts<-merge(cens12,cens23,all.x=T,all.y=T)
speciesagbmorts<-merge(speciesagbmorts,cens34,all.x=T,all.y=T)
speciesagbmorts[is.na(speciesagbmorts)]<-0 #Fix later in code somehow
mortsagb12<-data.frame(agb=speciesagbmorts$agb12,Species=speciesagbmorts$Species,Census="2008-2013")
mortsagb23<-data.frame(agb=speciesagbmorts$agb23,Species=speciesagbmorts$Species,Census="2013-2014")
mortsagb34<-data.frame(agb=speciesagbmorts$agb34,Species=speciesagbmorts$Species,Census="2014-2015")
speciesagbmort<-rbind(mortsagb12,mortsagb23,mortsagb34)
speciesagbmort <- with(speciesagbmort, speciesagbmort[order(Census, -as.numeric(agb)), ])
speciesagbmort$Species <- factor(speciesagbmort$Species, levels=speciesagbmort$Species) #Warning: levels deprecated 
####################################

#Plot biomass mortality by species
c2<-ggplot(data=speciesagbmort,aes(x=Species,y=speciesagbmort$agb)) + geom_bar(stat="identity",position="dodge",aes(fill=Census)) + labs(x = expression("Species"), y = expression(Biomass~Mortality~(Mg~ha^{-1}~y^{-1})))+ scale_fill_manual(values=c("darkblue","forestgreen","darkorange")) +scale_y_continuous(expand = c(0,0))+ theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1,size=10),text=element_text(size=12),axis.text.y=element_text(size=10),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.position="none") #+annotate("text",x=2,y=2.5,label="b)",size=8) 

#Return csv file with results
colnames(tab12)[1]<-"N0_2008-2013"
colnames(tab23)[1]<-"N0_2013-2014"
colnames(tab34)[1]<-"N0_2014-2015"
tbz<-merge(tab12,tab23)
tbz<-merge(tbz,tab34)
mortality_rates<-merge(speciesmorts,speciesagbmorts)
mortality_rates<-merge(mortality_rates,tbz)
colnames(mortality_rates)<-c("Species","mortrate2008_2013","ci.lo2008_2013","ci.hi2008_2013","mortrate2013_2014","ci.lo2013_2014","ci.hi2013_2014","mortrate2014_2015","ci.lo2014_2015","ci.hi2014_2015","agbmort2008_2013","agbmort2013_2014","agbmort2014_2015","N0_2008-2013","N0_2013-2014","N0_2014-2015")
write.csv(mortality_rates,"mortality_rates.csv")

##Calculate relationship between m and M for each census
rel12<-lm(mortality_rates$mortrate2008_2013 ~ mortality_rates$agbmort2008_2013)
rel23<-lm(mortality_rates$mortrate2013_2014 ~ mortality_rates$agbmort2013_2014)
rel34<-lm(mortality_rates$mortrate2014_2015 ~ mortality_rates$agbmort2014_2015)
relsum12<-summary(rel12)
relsum23<-summary(rel23)
relsum34<-summary(rel34)
rp.12<-pf(relsum12$fstatistic[1], relsum12$fstatistic[2], relsum12$fstatistic[3],lower.tail = FALSE)
rp.23<-pf(relsum23$fstatistic[1], relsum23$fstatistic[2], relsum23$fstatistic[3],lower.tail = FALSE)
rp.34<-pf(relsum34$fstatistic[1], relsum34$fstatistic[2], relsum34$fstatistic[3], lower.tail = FALSE)
relstats<-data.frame(Census=c("2008-13","2013-14","2014-15"),R2=c(relsum12$r.squared,relsum23$r.squared,relsum34$r.squared),p.value=c(rp.12,rp.23,rp.34))
write.csv(relstats,"relative_stats.csv")

#Insert Community-wide means here
allmortmean12<-3.4
allmortmean23<-1.9
allmortmean34<-2.9

#Histograms:
#Mortality
all.n.mort<-read.csv("Histcounts.mort.csv",header=T)
c3<-ggplot(data=all.n.mort,aes(x=midcut,y=n,fill=Census))+geom_bar(stat="identity",position="dodge",width=1.5)+ labs(y = expression("n Species"), x = expression(Mortality ~ Rate ~("%"~y^{-1})))+ scale_fill_manual(name="\nCensus Period\n",values=c("darkblue","forestgreen","darkorange"))+ scale_color_manual(name="\nCensus Period\n",values=c("darkblue","forestgreen","darkorange"))+scale_y_continuous(expand = c(0,0),limits=c(0,15))+scale_x_continuous(breaks=c(seq(0,20,2)),labels=c(0,"",4,"",8,"",12,"",16,"",20))+theme(text=element_text(size=12),axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.position="none")+geom_vline(xintercept=c(allmortmean12,allmortmean23,allmortmean34),linetype="dashed",col=c("darkblue","forestgreen","darkorange"))

#Biomass mortality
all.n.agb<-read.csv("Histcounts.agb.csv",header=T)
c4<-ggplot(data=all.n.agb,aes(x=midcut,y=n,fill=Census))+geom_bar(stat="identity",position="dodge",width=0.075)+ labs(y = expression("n Species"), x = expression(Biomass~ Mortality ~(Mg~ha^{-1}~y^{-1})))+ scale_fill_manual(name="\nCensus Period\n",values=c("darkblue","forestgreen","darkorange"))+ scale_color_manual(name="\nCensus Period\n",values=c("darkblue","forestgreen","darkorange"))+scale_y_continuous(expand = c(0,0),limits=c(0,15))+scale_x_continuous(breaks=c(seq(0,1,.1)),labels=c(0,"",.2,"",.4,"",.6,"",.8,"",1))+theme(text=element_text(size=12),axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.position="none")#+geom_vline(xintercept=c(magb12,magb23,magb34),linetype="dashed",col=c("darkblue","forestgreen","darkorange"))

#All Plots in one figure:
windowsFonts(A=windowsFont("Arial")) #Change font to Arial

tiff(file = "Figure4.tiff", width = 7.5, height = 7.5, units = "in", res=600, compression = "lzw",family ="A",pointsize=12)

plot_grid(c1,c2,c3,c4,labels=c("A","B","C","D"),label_size=12) #Color
#dev.off()

}#end runall
############################################

################################################
indcalcmort12=function(full1,full2,mindbh,rounddown=F)
################################################
{
  if(rounddown) {
    sm=(full1$dbh<55 | full2$dbh<55)
    full1$dbh[sm]=5*floor(full1$dbh[sm]/5)
    full2$dbh[sm]=5*floor(full2$dbh[sm]/5)
  }
  inc=full1$dbh>=mindbh&!is.na(full1$dbh)&!full1$sp!=full2$sp
  mort=ifelse(full2$dbh<0|is.na(full2$dbh),1,0)
  timeint=(full2$date-full1$date)/365.25
  agb=(full1$agb)/timeint #flip to make positive for dead trees
  mind=data.frame(species=full1$sp[inc],gx=full1$gx[inc],gy=full1$gy[inc],dinit=full1$dbh[inc],mort=mort[inc],timeint=timeint[inc],agb=agb[inc],tag=full1$tag[inc],stemID=full1$stemID[inc],status=full1$status[inc],dbh2=full2$dbh[inc])
  return(mind)
} # end indcalcmort12

################################################
indcalcmort23=function(allmort)
  ################################################
{
  inc=allmort$status.2013=="Live"
  mort=ifelse(allmort$status.2014=="Dead",1,0)
  timeint=(allmort$date.2014-allmort$date.2013)/365.25
  agb=(allmort$agb.2013)/timeint #flip to make positive for dead trees
  mind=data.frame(species=allmort$sp[inc],gx=allmort$gx[inc],gy=allmort$gy[inc],dinit=allmort$dbh.2013[inc],mort=mort[inc],timeint=timeint[inc],agb=agb[inc],tag=allmort$tag[inc],stemID=allmort$stemID[inc],status=allmort$status.2013[inc])
  return(mind)
} # end indcalcmort23

################################################
indcalcmort34=function(allmort)
  ################################################
{
  
  inc=allmort$status.2014=="Live"
  mort=ifelse(allmort$status.2015=="Dead",1,0)
  timeint=(allmort$date.2015-allmort$date.2014)/365.25
  #agb=allmort$agb.ctfs.2013/timeint
  #Multi-Year AGB Problema Corregida:
  best.agb<-(pmax(allmort$agb.2013,allmort$agb.2015,na.rm=T))/timeint
  mind=data.frame(species=allmort$sp[inc],gx=allmort$gx[inc],gy=allmort$gy[inc],dinit=allmort$dbh.2013[inc],mort=mort[inc],timeint=timeint[inc],agb=best.agb[inc],tag=allmort$tag[inc],stemID=allmort$stemID[inc],status=allmort$status.2014[inc])
  return(mind)
} # end indcalcmort34


################################################
mszmn=function(x,species,timeint)
  # function for calculating mean mortality rate by size intervals
################################################
{
  inc=!is.na(x) #Tori added is.na(x)
  x=x[inc]
  species=species[inc]
  timeint=timeint[inc]
  tab=table(species) #move to before !is(na)?
  nclass=length(tab)
  ninit=tapply(x,species,length)
  ndead=tapply(x,species,sum)
  mntime=tapply(timeint,species,mean,na.rm=T) #Tori added na.rm=T
  mortrate=100*(1-((ninit-ndead)/ninit)^(1/mntime))
  #mortrate=100*((log(ninit)-log(ninit-ndead))/mntime) #Helene's
  names(mortrate)=paste(rownames(ninit))
  mortrate<-data.frame(mortrate)
  mortrate$Species <- rownames(mortrate)
  mortrate$ninit<-ninit
  mortrate$ndead<-ndead
  mortrate$timeint<-mntime
  row.names(mortrate)<-seq(nrow(mortrate))
  return(mortrate)
} # end mszmn
######################################################

################################################
eszmn=function(x,species,hectar)
  ################################################
{
  inc=!is.na(x)
  x=x[inc]
  species=species[inc]
  tab=table(species) #move to before !is(na)?
  nclass=length(tab)
  xmn=tapply(x,species,sum,na.rm=T)
  agb=xmn/hectar
  agb<-data.frame(agb)
  agb$Species<-rownames(agb)
  row.names(agb)<-seq(nrow(agb))
  return(agb)
} # end eszmn
######################################################








