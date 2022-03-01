
pbs.all<-read.csv("All_SNPs_Chemosensory_Project_QCPass.csv",header=T, stringsAsFactors = T)


### Farming Populations as Derived (Hunter-Gatherers as Outgroup)####
##Subset into PBS comparisons and filter monomorphic SNPs between focal pairs
#Twa vs Bakiga, outgroup Agta
ug1<-subset(pbs.all,pbs.all$Twa.Bakiga<0)
ug2<-subset(pbs.all,pbs.all$Twa.Bakiga>0)
ug3<-subset(pbs.all,pbs.all$Twa.Bakiga==0)
twa.pbs<-rbind(ug1,ug2,ug3)
remove(ug1,ug2,ug3)

#Sua vs Bakiga, outgroup Agta
ug1<-subset(pbs.all,pbs.all$Sua.Bakiga<0)
ug2<-subset(pbs.all,pbs.all$Sua.Bakiga>0)
ug3<-subset(pbs.all,pbs.all$Sua.Bakiga==0)
sua.pbs<-rbind(ug1,ug2,ug3)
remove(ug1,ug2,ug3)

#Agta vs Manobo, outgroup Twa
ph1<-subset(pbs.all,pbs.all$Agta.Manobo<0)
ph2<-subset(pbs.all,pbs.all$Agta.Manobo>0)
ph3<-subset(pbs.all,pbs.all$Agta.Manobo==0)
agta.pbs<-rbind(ph1,ph2,ph3)
remove(ph1,ph2,ph3)

#Mamanwa vs Manobo, outgroup Twa
ph1<-subset(pbs.all,pbs.all$Mamanwa.Manobo<0)
ph2<-subset(pbs.all,pbs.all$Mamanwa.Manobo>0)
ph3<-subset(pbs.all,pbs.all$Mamanwa.Manobo==0)
mamanwa.pbs<-rbind(ph1,ph2,ph3)
remove(ph1,ph2,ph3)

#change NAN to 0 in outgroup. ##outgroups are Agta and Twa
sua.pbs[,44][is.na(sua.pbs[,44])] <- 0  ##Agta.Sua,
sua.pbs[,42][is.na(sua.pbs[,42])] <- 0  ##Agta.Bakiga2
twa.pbs[,45][is.na(twa.pbs[,45])] <- 0  ##Agta.Twa,
twa.pbs[,42][is.na(twa.pbs[,42])] <- 0  ##Agta.Bakiga2
agta.pbs[,45][is.na(agta.pbs[,45])] <- 0  ##Agta.Twa,
agta.pbs[,48][is.na(agta.pbs[,48])] <- 0  ##Manobo.Twa
mamanwa.pbs[,47:48][is.na(mamanwa.pbs[,47:48])] <- 0  ##Mamanwa.Twa & Manobo.Twa

##Calculate PBS ###
#Bakiga v Sua, outgroup Agta
sua.pbs$t.ab<--log(1-sua.pbs$Sua.Bakiga) #HG vs AG
sua.pbs$t.ao<--log(1-sua.pbs$Agta.Sua) #HG vs outgroup
sua.pbs$t.bo<--log(1-sua.pbs$Agta.Bakiga) #AG vs outgroup

sua.pbs$PBS<-((sua.pbs$t.ab + sua.pbs$t.bo - sua.pbs$t.ao)/2)
sua.pbs$PBS[which(sua.pbs$PBS<0)]=0 # set to 0 negative values, for convenience

#Bakiga v Twa, outgroup Agta
twa.pbs$t.ab<--log(1-twa.pbs$Twa.Bakiga) #HG vs AG
twa.pbs$t.ao<--log(1-twa.pbs$Agta.Twa) #HG vs outgroup
twa.pbs$t.bo<--log(1-twa.pbs$Agta.Bakiga) #AG vs outgroup

twa.pbs$PBS<-((twa.pbs$t.ab + twa.pbs$t.bo - twa.pbs$t.ao)/2)
twa.pbs$PBS[which(twa.pbs$PBS<0)]=0 # set to 0 negative values, for convenience

#Mamanwa v Manobo, outgroup Twa
mamanwa.pbs$t.ab<--log(1-mamanwa.pbs$Mamanwa.Manobo) #HG vs AG
mamanwa.pbs$t.ao<--log(1-mamanwa.pbs$Mamanwa.Twa) #HG vs outgroup
mamanwa.pbs$t.bo<--log(1-mamanwa.pbs$Manobo.Twa) #AG vs outgroup

mamanwa.pbs$PBS<-((mamanwa.pbs$t.ab + mamanwa.pbs$t.bo - mamanwa.pbs$t.ao)/2)
mamanwa.pbs$PBS[which(mamanwa.pbs$PBS<0)]=0 # set to 0 negative values, for convenience

#Agta v Manobo, outgroup Twa
agta.pbs$t.ab<--log(1-agta.pbs$Agta.Manobo) #HG vs AG
agta.pbs$t.ao<--log(1-agta.pbs$Agta.Manobo) #HG vs outgroup
agta.pbs$t.bo<--log(1-agta.pbs$Agta.Manobo) #AG vs outgroup

agta.pbs$PBS<-((agta.pbs$t.ab + agta.pbs$t.bo - agta.pbs$t.ao)/2)
agta.pbs$PBS[which(agta.pbs$PBS<0)]=0 # set to 0 negative values, for convenience

###### PBS for combined HG populations per region #####
##Subset into PBS comparisons and filter monomorphic SNPs between focal pairs

#Uganda HG vs Bakiga, outgroup Agta
ug1<-subset(pbs.all,pbs.all$Bakiga.UgandaHG<0)
ug2<-subset(pbs.all,pbs.all$Bakiga.UgandaHG>0)
ug3<-subset(pbs.all,pbs.all$Bakiga.UgandaHG==0)
uganda.pbs<-rbind(ug1,ug2,ug3)
remove(ug1,ug2,ug3)

#Philipines HG vs Manobo, outgroup Twa
ph1<-subset(pbs.all,pbs.all$Manobo.PhilHG<0)
ph2<-subset(pbs.all,pbs.all$Manobo.PhilHG>0)
ph3<-subset(pbs.all,pbs.all$Manobo.PhilHG==0)
phil.pbs<-rbind(ph1,ph2,ph3)
remove(ph1,ph2,ph3)

#change NAN to 0 in outgroup. ##outgroups are Agta and Twa
phil.pbs[,55][is.na(phil.pbs[,55])] <- 0  ##Twa.PhilHG,
phil.pbs[,48][is.na(phil.pbs[,48])] <- 0  ##Manobo.Twa
uganda.pbs[,57][is.na(uganda.pbs[,57])] <- 0  ##Agta.UgandaHG,
uganda.pbs[,42][is.na(uganda.pbs[,42])] <- 0  ##Agta.Bakiga

##Calculate PBS ###
#Bakiga v UgandaHG, outgroup Agta
uganda.pbs$t.ab<--log(1-uganda.pbs$Bakiga.UgandaHG) #HG vs AG
uganda.pbs$t.ao<--log(1-uganda.pbs$Agta.UgandaHG) #HG vs outgroup
uganda.pbs$t.bo<--log(1-uganda.pbs$Agta.Bakiga) #AG vs outgroup

uganda.pbs$PBS<-((uganda.pbs$t.ab + uganda.pbs$t.bo - uganda.pbs$t.ao)/2)
uganda.pbs$PBS[which(uganda.pbs$PBS<0)]=0 # set to 0 negative values, for convenience

#Philipines HG v Manobo, outgroup Twa
phil.pbs$t.ab<--log(1-phil.pbs$Manobo.PhilHG) #HG vs AG
phil.pbs$t.ao<--log(1-phil.pbs$Twa.PhilHG) #HG vs outgroup
phil.pbs$t.bo<--log(1-phil.pbs$Manobo.Twa) #AG vs outgroup

phil.pbs$PBS<-((phil.pbs$t.ab + phil.pbs$t.bo - phil.pbs$t.ao)/2)
phil.pbs$PBS[which(phil.pbs$PBS<0)]=0 # set to 0 negative values, for convenience

#####empirical pvalue for Farmers (HG outgroup) ####

#Manobo vs. Agta
agta.neutral<-subset(agta.pbs,Analysis_Set=="neutral_SNPs")
agta.exonic<-subset(agta.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(agta.neutral$PBS)
agta.exonic$PBS.p<-1-neutral.Pdist(agta.exonic$PBS)
agta.exonic$PBS.fdr<-p.adjust(agta.exonic$PBS.p,method="fdr")
agta.exonic.PBS.outliers<-subset(agta.exonic,agta.exonic$PBS.fdr<=0.05)

#Manobo vs. Mamanwa
mam.neutral<-subset(mamanwa.pbs,Analysis_Set=="neutral_SNPs")
mam.exonic<-subset(mamanwa.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(mam.neutral$PBS)
mam.exonic$PBS.p<-1-neutral.Pdist(mam.exonic$PBS)
mam.exonic$PBS.fdr<-p.adjust(mam.exonic$PBS.p,method="fdr")
mam.exonic.PBS.outliers<-subset(mam.exonic,mam.exonic$PBS.fdr<=0.05)

#BaKiga vs. Sua
sua.neutral<-subset(sua.pbs,Analysis_Set=="neutral_SNPs")
sua.exonic<-subset(sua.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(sua.neutral$PBS)
sua.exonic$PBS.p<-1-neutral.Pdist(sua.exonic$PBS)
sua.exonic$PBS.fdr<-p.adjust(sua.exonic$PBS.p,method="fdr")
sua.exonic.PBS.outliers<-subset(sua.exonic,sua.exonic$PBS.fdr<=0.05)

#BaKiga vs. Twa
twa.neutral<-subset(twa.pbs,Analysis_Set=="neutral_SNPs")
twa.exonic<-subset(twa.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(twa.neutral$PBS)
twa.exonic$PBS.p<-1-neutral.Pdist(twa.exonic$PBS)
twa.exonic$PBS.fdr<-p.adjust(twa.exonic$PBS.p,method="fdr")
twa.exonic.PBS.outliers<-subset(twa.exonic,twa.exonic$PBS.fdr<=0.05)

#Manobo vs. Philippines HG
phil.neutral<-subset(phil.pbs,Analysis_Set=="neutral_SNPs")
phil.exonic<-subset(phil.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(phil.neutral$PBS)
phil.exonic$PBS.p<-1-neutral.Pdist(phil.exonic$PBS)
phil.exonic$PBS.fdr<-p.adjust(phil.exonic$PBS.p,method="fdr")
phil.exonic.PBS.outliers<-subset(phil.exonic,phil.exonic$PBS.fdr<=0.05)

#BaKiga vs. Uganda HG
uganda.neutral<-subset(uganda.pbs,Analysis_Set=="neutral_SNPs")
uganda.exonic<-subset(uganda.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(uganda.neutral$PBS)
uganda.exonic$PBS.p<-1-neutral.Pdist(uganda.exonic$PBS)
uganda.exonic$PBS.fdr<-p.adjust(uganda.exonic$PBS.p,method="fdr")
uganda.exonic.PBS.outliers<-subset(uganda.exonic,uganda.exonic$PBS.fdr<=0.05)


##### Hunter-Gatherer as Derived Group (Farming as Outgroup) #####

#change NAN to 0 in outgroup. ##outgroups are BaKiga and Manobo
sua.pbs[,53][is.na(sua.pbs[,53])] <- 0  ##Sua.Manobo
sua.pbs[,52][is.na(sua.pbs[,52])] <- 0  ##Manobo.Bakiga

twa.pbs[,48][is.na(twa.pbs[,48])] <- 0  ##Twa.Manobo,
twa.pbs[,52][is.na(twa.pbs[,52])] <- 0  ##Manobo.Bakiga

agta.pbs[,42][is.na(agta.pbs[,42])] <- 0  ##Agta.Bakiga,
agta.pbs[,52][is.na(agta.pbs[,52])] <- 0  ##Manobo.Bakiga

mamanwa.pbs[,51][is.na(mamanwa.pbs[,51])] <- 0  ##Mamanwa Bakiga
mamanwa.pbs[,52][is.na(mamanwa.pbs[,52])] <- 0  ##Manobo Bakiga

uganda.pbs[,59][is.na(uganda.pbs[,59])] <- 0  ##Manobo.UgandaHG
uganda.pbs[,52][is.na(uganda.pbs[,52])] <- 0  ##Manobo Bakiga

phil.pbs[,58][is.na(phil.pbs[,58])] <- 0  ##Bakiga.PhilUG
phil.pbs[,52][is.na(phil.pbs[,52])] <- 0  ##Manobo Bakiga

#BaSua v Bakiga, outgroup Manobo
sua.pbs$t.ab<--log(1-sua.pbs$Sua.Bakiga) #HG vs AG
sua.pbs$t.ao<--log(1-sua.pbs$Sua.Manobo) #HG vs outgroup
sua.pbs$t.bo<--log(1-sua.pbs$Manobo.Bakiga) #AG vs outgroup

sua.pbs$PBS.HG<-((sua.pbs$t.ab + sua.pbs$t.ao - sua.pbs$t.bo)/2)
sua.pbs$PBS.HG[which(sua.pbs$PBS.HG<0)]=0 # set to 0 negative values, for convenience

#Twa v Bakiga, outgroup Manobo
twa.pbs$t.ab<--log(1-twa.pbs$Twa.Bakiga) #HG vs AG
twa.pbs$t.ao<--log(1-twa.pbs$Manobo.Twa) #HG vs outgroup
twa.pbs$t.bo<--log(1-twa.pbs$Manobo.Bakiga) #AG vs outgroup

twa.pbs$PBS.HG<-((twa.pbs$t.ab + twa.pbs$t.ao - twa.pbs$t.bo)/2)
twa.pbs$PBS.HG[which(twa.pbs$PBS.HG<0)]=0 # set to 0 negative values, for convenience

#Mamanwa v Manobo, outgroup BaKiga
mamanwa.pbs$t.ab<--log(1-mamanwa.pbs$Mamanwa.Manobo) #HG vs AG
mamanwa.pbs$t.ao<--log(1-mamanwa.pbs$Mamanwa.Bakiga) #HG vs outgroup
mamanwa.pbs$t.bo<--log(1-mamanwa.pbs$Manobo.Bakiga) #AG vs outgroup

mamanwa.pbs$PBS.HG<-((mamanwa.pbs$t.ab + mamanwa.pbs$t.ao - mamanwa.pbs$t.bo)/2)
mamanwa.pbs$PBS.HG[which(mamanwa.pbs$PBS.HG<0)]=0 # set to 0 negative values, for convenience

#Agta v Manobo, outgroup BaTwa
agta.pbs$t.ab<--log(1-agta.pbs$Agta.Manobo) #HG vs AG
agta.pbs$t.ao<--log(1-agta.pbs$Agta.Bakiga) #HG vs outgroup
agta.pbs$t.bo<--log(1-agta.pbs$Manobo.Bakiga) #AG vs outgroup

agta.pbs$PBS.HG<-((agta.pbs$t.ab + agta.pbs$t.ao - agta.pbs$t.bo)/2)
agta.pbs$PBS.HG[which(agta.pbs$PBS.HG<0)]=0 # set to 0 negative values, for convenience

#UgandaHG v Bakiga, outgroup Manobo
uganda.pbs$t.ab<--log(1-uganda.pbs$Bakiga.UgandaHG) #HG vs AG
uganda.pbs$t.ao<--log(1-uganda.pbs$Manobo.UgandaHG) #HG vs outgroup
uganda.pbs$t.bo<--log(1-uganda.pbs$Manobo.Bakiga) #AG vs outgroup

uganda.pbs$PBS.HG<-((uganda.pbs$t.ab + uganda.pbs$t.ao - uganda.pbs$t.bo)/2)
uganda.pbs$PBS.HG[which(uganda.pbs$PBS.HG<0)]=0 # set to 0 negative values, for convenience

#Phillipines HG v Manobo, outgroup BaTwa
phil.pbs$t.ab<--log(1-phil.pbs$Manobo.PhilHG) #HG vs AG
phil.pbs$t.ao<--log(1-phil.pbs$Bakiga.PhilHG) #HG vs outgroup
phil.pbs$t.bo<--log(1-phil.pbs$Manobo.Bakiga) #AG vs outgroup

phil.pbs$PBS.HG<-((phil.pbs$t.ab + phil.pbs$t.ao - phil.pbs$t.bo)/2)
phil.pbs$PBS.HG[which(phil.pbs$PBS.HG<0)]=0 # set to 0 negative values, for convenience

#####empirical pvalue for Hunter-Gatherers ####
#Agta vs. Manobo
agta.neutral<-subset(agta.pbs,Analysis_Set=="neutral_SNPs")
agta.exonic<-subset(agta.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(agta.neutral$PBS.HG)
agta.exonic$PBS.p<-1-neutral.Pdist(agta.exonic$PBS.HG)
agta.exonic$PBS.fdr<-p.adjust(agta.exonic$PBS.p,method="fdr")
agta.exonic.PBS.HG.outliers<-subset(agta.exonic,agta.exonic$PBS.fdr<=0.05)

#Mamanwa vs. Manobo
mam.neutral<-subset(mamanwa.pbs,Analysis_Set=="neutral_SNPs")
mam.exonic<-subset(mamanwa.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(mam.neutral$PBS.HG)
mam.exonic$PBS.p<-1-neutral.Pdist(mam.exonic$PBS.HG)
mam.exonic$PBS.fdr<-p.adjust(mam.exonic$PBS.p,method="fdr")
mam.exonic.PBS.HG.outliers<-subset(mam.exonic,mam.exonic$PBS.fdr<=0.05)

#Sua vs. BaKiga
sua.neutral<-subset(sua.pbs,Analysis_Set=="neutral_SNPs")
sua.exonic<-subset(sua.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(sua.neutral$PBS.HG)
sua.exonic$PBS.p<-1-neutral.Pdist(sua.exonic$PBS.HG)
sua.exonic$PBS.fdr<-p.adjust(sua.exonic$PBS.p,method="fdr")
sua.exonic.PBS.HG.outliers<-subset(sua.exonic,sua.exonic$PBS.fdr<=0.05)

#Twa vs. BaKiga
twa.neutral<-subset(twa.pbs,Analysis_Set=="neutral_SNPs")
twa.exonic<-subset(twa.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(twa.neutral$PBS.HG)
twa.exonic$PBS.p<-1-neutral.Pdist(twa.exonic$PBS.HG)
twa.exonic$PBS.fdr<-p.adjust(twa.exonic$PBS.p,method="fdr")
twa.exonic.PBS.HG.outliers<-subset(twa.exonic,twa.exonic$PBS.fdr<=0.05)

#Philippines HG vs Manobo
phil.neutral<-subset(phil.pbs,Analysis_Set=="neutral_SNPs")
phil.exonic<-subset(phil.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(phil.neutral$PBS.HG)
phil.exonic$PBS.p<-1-neutral.Pdist(phil.exonic$PBS.HG)
phil.exonic$PBS.fdr<-p.adjust(phil.exonic$PBS.p,method="fdr")
phil.exonic.PBS.HG.outliers<-subset(phil.exonic,phil.exonic$PBS.fdr<=0.05)

#Uganda HG vs. 
uganda.neutral<-subset(uganda.pbs,Analysis_Set=="neutral_SNPs")
uganda.exonic<-subset(uganda.pbs,Analysis_Set=="sensory_exonic")
neutral.Pdist=ecdf(uganda.neutral$PBS.HG)
uganda.exonic$PBS.p<-1-neutral.Pdist(uganda.exonic$PBS.HG)
uganda.exonic$PBS.fdr<-p.adjust(uganda.exonic$PBS.p,method="fdr")
uganda.exonic.PBS.HG.outliers<-subset(uganda.exonic,uganda.exonic$PBS.fdr<=0.05)

