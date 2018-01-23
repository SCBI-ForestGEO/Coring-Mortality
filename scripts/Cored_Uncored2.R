##########################################################################
#Mortality rates of trees cored in 2009 vs uncored trees 
#based on mortality censuses in 2013,14,15,16 and 17
#For use with SCBI ForestGEO grid
#Author: Ryan Helcoski - 12/14/2017
##########################################################################

rm(list = ls())

library(dplyr)
library (kimisc)

#############################################################
##### creating wd and loading data
# Set up working directory ####
setwd("C:/Users/helcoskiR/Dropbox (Smithsonian)/Tree Cores/literature/effect of coring on mortality/Data/R_project")

# INPUT DATA LOCATION ####
Input_data_location <- "INPUT_FILES/"

# OUTPUT DATA LOCATION ####
Output_data_location <- "OUTPUT_FILES/"

# load data ####
# load 2017 mortality data + trees cored and uncored. ####
#File is csv of 2017 mortality with one extra column, Cored_2009. A 1 indicates a tree was cored a 0 indicates it was not
cmortall <- read.csv(paste0(Input_data_location, "Mortality_Survey_2017+core2009.csv"))

 

############################################################
## adding columns for genus and size ####
#### remove unnecessary columns #
cmortall[14:34] <-NULL
cmortall[5:6] <-NULL


###subet of only genus included in csub = acer, carya, fagus, fraxinus, juglans, liriodendron, nyssa, quercus, tilia, ulmus
cmortall <- subset(cmortall, sp!= "aial" & sp!= "caca" & sp!= "amar" & sp!= "ceca"& sp!= "astr" & sp!= "casp"& sp!= "ceca"& sp!= "ceoc"& sp!= "coal"& sp!= "cofl"& sp!= "divi"& sp!= "frsp"& sp!= "havi"& sp!= "libe"& sp!= "pato"& sp!= "pipu"& sp!= "pist"& sp!= "pivi"& sp!= "ploc"& sp!= "prav"& sp!= "prse"& sp!= "qusp"& sp!= "rops"& sp!= "saal"& sp!= "ulsp"& sp!= "unk"& sp!= "vipr")

#create column for genus ####
cmortall$genus <- "NA"
unique(cmortall$sp)
#acer
cmortall$genus[which(cmortall$sp == "acru" | cmortall$sp == "acne" | cmortall$sp == "acpl")] <- "acer"
#carya
cmortall$genus[which(cmortall$sp == "caco" | cmortall$sp == "cagl" | cmortall$sp == "caovl" | cmortall$sp == "caovl" | cmortall$sp == "cade"| cmortall$sp == "cato")] <- "carya"
#fagus
cmortall$genus[which(cmortall$sp == "fagr" )] <- "fagus"
#fraxinus
cmortall$genus[which(cmortall$sp == "fram" | cmortall$sp == "frni" |cmortall$sp == "frpe" )] <- "fraxinus"
#juglans
cmortall$genus[which(cmortall$sp == "juci" | cmortall$sp == "juni")] <- "juglans"
#liriodendron
cmortall$genus[which(cmortall$sp == "litu" )] <- "liriodendron"
#nyssa
cmortall$genus[which(cmortall$sp == "nysy" )] <- "nyssa"
#quercus
cmortall$genus[which(cmortall$sp == "qual" | cmortall$sp == "quco" | cmortall$sp == "qufa" | cmortall$sp == "qumi" | cmortall$sp == "qumu"| cmortall$sp == "qupr" | cmortall$sp == "quru" | cmortall$sp == "quve")] <- "quercus"
#tilia
cmortall$genus[which(cmortall$sp == "tiam" )] <- "tillia"
#ulmus
cmortall$genus[which(cmortall$sp == "ulam" | cmortall$sp == "ulru")] <- "ulmus"
#remove any extra "NA"
cmortall <- subset(cmortall, genus!="NA")
table(cmortall$genus)

#create column for size bin
cmortall$size_bin <- "NA"
cmortall$size_bin [which(cmortall$dbh.2013 >= 350)] <- 2
cmortall$size_bin [which(cmortall$dbh.2013 < 350)] <- 1
table(cmortall$size_bin)

#remove all fraxinus (and others?) smaller than 100dbh
cmortall <- cmortall[!(cmortall$dbh.2013 < 100) , ]
summary(cmortall$new.status)

#create new column status.2017, where new.status A and AU = A, DC,DS,DP = D
cmortall$status.2017 <- "NA"
cmortall[cmortall$new.status == 'A' | cmortall$new.status == 'AU'| cmortall$new.status =='A ', "status.2017"] <- 'A'
cmortall[cmortall$new.status == 'PD' | cmortall$new.status == 'DS'| cmortall$new.status =='DC', "status.2017"] <- 'D'
table(cmortall$status.2017)

##############################################################################
# creating 2 DFs ####
# DF1, only trees cored, minus pinus and robinia ####
# DF A filter out only those cored
csub <- subset(cmortall, Cored_2009 == 1)

#pivot table, split into genus and size bin. Create new column for desired number of noncored needed - multiplier is 2x
csub = group_by(csub, genus, size_bin)
csub.table2 = summarise(csub, cored= sum(Cored_2009) , noncored_needed =sum((Cored_2009)*2))
csub.table2

#create tally column, each stem is just 1
csub$tally1 <- 1


# DF2, uncorred trees with only possible trees to comapre to ####
# subset of only noncored
ncsub <- subset(cmortall, Cored_2009 == 0)
#create tally column, each stem is just 1
ncsub$tally1 <- 1

#pivot table. Create new column to 
ncsub = group_by(ncsub, genus, size_bin)
ncsub.table2 = summarise(ncsub, noncored_total = sum(tally1))
ncsub.table2

#merge pivot tables
allsub.table <- merge(csub.table2, ncsub.table2, by.x = c("genus" , "size_bin"), by.y = c("genus" , "size_bin" ))
allsub.table
t6 <- allsub.table[!(allsub.table$noncored_needed > allsub.table$noncored_total),]
t6

#ulmus and fagus size class 2 (>= 350mm dbh) have been removed
#remove ulmus and fagus >= 350 from DF1 and DF2
csub2 <- csub[!(csub$genus == "fagus" & csub$size_bin == 2), ]
csub2 <- csub2[!(csub2$genus == "ulmus" & csub2$size_bin == 2), ]
csub3 <- csub2
table(csub3$genus)
csub5 = group_by(csub3, genus, size_bin)
csub5 = summarise(csub5, total= sum(tally1))
write.csv(csub5, "Cored_Species_Pivot.csv")
csub6 = group_by(csub3, size_bin)
csub6 = summarise(csub6, total= sum(tally1))
write.csv(csub6, "Cored_Total_Pivot.csv")




ncsub2 <- ncsub[!(ncsub$genus == "fagus" & ncsub$size_bin == 2), ]
ncsub2 <- ncsub2[!(ncsub2$genus == "ulmus" & ncsub2$size_bin == 2), ]
table(ncsub2$genus)


#randomly selecting desired numbers of noncored df
t1 <- sample.rows(subset(ncsub2, genus == "acer" & size_bin == 1), 50)
t2 <- sample.rows(subset(ncsub2, genus == "acer" & size_bin == 2), 16)
t3 <- sample.rows(subset(ncsub2, genus == "carya" & size_bin == 1), 396)
t4 <- sample.rows(subset(ncsub2, genus == "carya" & size_bin == 2), 116)
t5 <- sample.rows(subset(ncsub2, genus == "fagus" & size_bin == 1), 172)
t6 <- sample.rows(subset(ncsub2, genus == "fraxinus" & size_bin == 1), 34)
t7 <- sample.rows(subset(ncsub2, genus == "fraxinus" & size_bin == 2), 30)
t8 <- sample.rows(subset(ncsub2, genus == "juglans" & size_bin == 1), 18)
t9 <- sample.rows(subset(ncsub2, genus == "juglans" & size_bin == 2), 52)
t10 <- sample.rows(subset(ncsub2, genus == "liriodendron" & size_bin == 1), 90)
t11 <- sample.rows(subset(ncsub2, genus == "liriodendron" & size_bin == 2), 90)
t12 <- sample.rows(subset(ncsub2, genus == "nyssa" & size_bin == 1), 98)
t13 <- sample.rows(subset(ncsub2, genus == "nyssa" & size_bin == 2), 14)
t14 <- sample.rows(subset(ncsub2, genus == "quercus" & size_bin == 1), 168)
t15 <- sample.rows(subset(ncsub2, genus == "quercus" & size_bin == 2), 386)
t16 <- sample.rows(subset(ncsub2, genus == "tillia" & size_bin == 1), 54)
t17 <- sample.rows(subset(ncsub2, genus == "tillia" & size_bin == 2), 18)
t18 <- sample.rows(subset(ncsub2, genus == "ulmus" & size_bin == 1), 68)

#r bind, view
ncsub3 <- rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18)
ncsub3
ncsub5 = group_by(ncsub3, genus, size_bin)
ncsub5 = summarise(ncsub5, total= sum(tally1))
write.csv(ncsub5, "Non-cored_Species_Pivot.csv")
ncsub6 = group_by(ncsub3, size_bin)
ncsub6 = summarise(ncsub6, total= sum(tally1))
write.csv(ncsub6, "Non-Cored_Total_Pivot.csv")


#############################################################################
#Pivot tables of DF1 and DF2, looking for total cored and uncored and total from each group that died by 2017####

# cored trees
csub4 = group_by(csub3, genus, size_bin, status.2017)
pvcsub4 = summarise(csub4, total = sum(tally1))
write.csv(pvcsub4, "Table_Cored_Tree_Mortality.csv")

csubtotal = group_by(csub3, size_bin, status.2017)
csubtotal = summarise(csubtotal, total = sum(tally1))
write.csv(csubtotal, "Total_Cored_Tree_Mortality.csv")


# noncored trees
ncsub4 = group_by(ncsub3, genus, size_bin, status.2017)
pvncsub4 = summarise(ncsub4, total = sum(tally1))
write.csv(pvncsub4, "Table_Non-Cored_Tree_Mortality.csv")

ncsubtotal = group_by(ncsub3, size_bin, status.2017)
ncsubtotal = summarise(ncsubtotal, total = sum(tally1))
write.csv(ncsubtotal, "Total_Non-Cored_Tree_Mortality.csv")



##########################################################################
# annual mortality rate, m = (1- (Nt/N0)^(1/t))x100 #### 

# load in table of more precise t
t<- read.csv("New_t_cored_uncored.csv")
# uncorred using trees alive in 2008

# cored
#total anual mortality
csubalive <- csubtotal[!(csubtotal$status.2017 == "D") , ]
c_mort_rate <- (1-((sum(csubalive$total))/(sum(csub6$total)))^(1/t$cored_t)) * 100
#total by size class
csub6$mort_rate <- (1-((csubalive$total)/(csub6$total))^(1/t$cored_t)) * 100
#total by species and size class
csublivespecies <- pvcsub4[!(pvcsub4$status.2017 == "D") , ]
csub5 <- merge(csub5, csublivespecies, by.x=c("genus", "size_bin"), by.y =c("genus", "size_bin") )
csub5$mort_rate <- (1-((csub5$total.y)/(csub5$total.x))^(1/t$cored_t)) * 100



#uncored
#total anual mortality
ncsubalive <- ncsubtotal[!(ncsubtotal$status.2017 == "D") , ]
nc_mort_rate <- (1-((sum(ncsubalive$total))/(sum(ncsub6$total)))^(1/t$uncored_t)) * 100
total_mort <- cbind(nc_mort_rate, c_mort_rate)
#total by size class
ncsub6$mort_rate <- (1-((ncsubalive$total)/(ncsub6$total))^(1/t$uncored_t)) * 100
#total by species and size class
ncsublivespecies <- pvncsub4[!(pvncsub4$status.2017 == "D") , ]
ncsub5 <- merge(ncsub5, ncsublivespecies, by.x=c("genus", "size_bin"), by.y =c("genus", "size_bin") )
ncsub5$mort_rate <- (1-((ncsub5$total.y)/(ncsub5$total.x))^(1/t$uncored_t)) * 100



##########################################################################
# clean up data, prep for 95% CI ####
#remove unecessary columns
csub5 <- csub5[, -c(1,4,7)]
ncsub5 <- ncsub5[, -c(1,4,6)]
csub6 <- csub6[, -c(1,2)]
ncsub6 <- ncsub6[, -c(1,2)]

#ninitial, nfinal
csub5<- rename(csub5, n_initial=total.x , n_final=total.y)
ncsub5<- rename(ncsub5, n_initial=total.x , n_final=total.y)

ncsubdead <- ncsubtotal[!(ncsubtotal$status.2017 == "A") , ]
ncsub6 <- merge(ncsub6, ncsubdead, by.x="size_bin", by.y="size_bin")
ncsub6$n_final <- ncsub6$total.x - ncsub6$total.y
ncsub6 <- ncsub6[, -c(4,5)]
ncsub6 <- rename(ncsub6, n_initial=total.x, n_dead=total.y)

csubdead <- csubtotal[!(csubtotal$status.2017 == "A") , ]
csub6 <- merge(csub6, csubdead, by.x="size_bin", by.y="size_bin")
csub6$n_final <- csub6$total.x - csub6$total.y
csub6 <- csub6[, -c(4,5)]
csub6 <- rename(csub6, n_initial=total.x, n_dead=total.y)



#bring back ndead
csub5$n_dead <- (csub5$n_initial - csub5$n_final)
ncsub5$n_dead <- (ncsub5$n_initial - ncsub5$n_final)
csub6$n_dead <- (csub6$n_initial - csub6$n_final)
ncsub6$n_dead <- (ncsub6$n_initial - ncsub6$n_final)


#t
csub5$t <- t$cored_t
ncsub5$t <- t$uncored_t
csub6$t <- t$cored_t
ncsub6$t <- t$uncored_t


#fixing total_mort
c_initial <- sum(csub6$n_initial)
c_dead <- sum(csub6$n_dead)
c_mort <- total_mort$c_mort_rate
c_final <- sum(csub6$n_final)
c = cbind(c_initial, c_dead, c_mort, c_final)
nc_initial <- sum(ncsub6$n_initial)
nc_dead <- sum(ncsub6$n_dead)
nc_mort <- total_mort$nc_mort_rate
nc_final <- sum(ncsub6$n_final)
nc = cbind(nc_initial, nc_dead, nc_mort, nc_final)
total_mort = rbind(c,nc)
total_mort <- data.frame(total_mort)
total_mort <- rename(total_mort, n_initial= c_initial, n_dead=c_dead, mort_rate=c_mort, n_final=c_final)
type = c("cored" , "uncored")
t2 = c(t$cored_t, t$uncored_t)
t4 = c("total", "total")
total_mort$t <- t2
total_mort$genus <- type
total_mort$size_bin <- t4

total_mort$type <- NULL


#####################################################################
# creating large dfs ####

#cored
#fix data
csub5$status.2017 <- NULL
size = c("<35 dbh", ">=35 dbh")
csub6$genus <- size
csub6 <- rename(csub6, mort_rate = cmort)
cored = rbind.data.frame(csub5, csub6, total_mort)
cored = cored[-22,]
write.csv(cored, "Cored_mortality.csv")


#uncored
ncsub5$status.2017 <- NULL
size = c("<35 dbh", ">=35 dbh")
ncsub6$genus <- size
uncored = rbind.data.frame(ncsub5, ncsub6, total_mort)
uncored = uncored[-21,]
write.csv(uncored, "Uncored_mortality.csv")



#####################################################################
# 95% CI #####
#test with uncored species
#Calculate CIS: Normal Approximation
uncored2<-uncored[uncored$n_dead>5,]
ntot<-nrow(uncored2)
for (i in 1:ntot){
  prop<-prop.test(uncored2$n_dead[i],uncored2$n_initial[i],uncored2$mort_rate[i]/100)
  ci1<-prop$conf.int[1]*uncored2$n_initial[i]
  ci2<-prop$conf.int[2]*uncored2$n_initial[i]
  ci.lo<-100*(1-((uncored2$n_initial[i]-ci1)/uncored2$n_initial[i])^(1/uncored2$t[i]))
  ci.hi<-100*(1-((uncored2$n_initial[i]-ci2)/uncored2$n_initial[i])^(1/uncored2$t[i]))
  prop=data.frame(genus=uncored2$genus[i],   size_bin = uncored2$size_bin[i],n_initial=uncored2$n_initial[i], mort_rate=uncored2$mort_rate[i], n_final = uncored2$n_final[i], n_dead = uncored2$n_dead[i]   ,    ci.lo=ci.lo,ci.hi=ci.hi)
  if (i==1)
    props=prop
  else
    props=rbind(props,prop)}

#If ndead <= 5: Exact Binomial Test
uncored_2<-uncored[uncored$n_dead<=5&uncored$n_dead>0,]
ntots<-nrow(uncored_2)
for (i in 1:ntots){
  binom<-binom.test(uncored_2$n_dead[i],uncored_2$n_initial[i],uncored_2$mort_rate[i]/100)
  ci1<-binom$conf.int[1]*uncored_2$n_initial[i]
  ci2<-binom$conf.int[2]*uncored_2$n_initial[i]
  ci.lo<-100*(1-((uncored_2$n_initial[i]-ci1)/uncored_2$n_initial[i])^(1/uncored_2$t[i]))
  ci.hi<-100*(1-((uncored_2$n_initial[i]-ci2)/uncored_2$n_initial[i])^(1/uncored_2$t[i]))
  Binom=data.frame(genus=uncored_2$genus[i],size_bin = uncored_2$size_bin[i],n_initial=uncored_2$n_initial[i], mort_rate=uncored_2$mort_rate[i], n_final = uncored_2$n_final[i], n_dead = uncored_2$n_dead[i],ci.lo=ci.lo,ci.hi=ci.hi)
  if (i==1)
    binoms=Binom
  else
    binoms=rbind(binoms,Binom)}

alluncored<-rbind(props,binoms)
a2<-merge(uncored,alluncored,all.x=T,all.y=T)

write.csv(a2, "Uncored_mortality+CI.csv")





