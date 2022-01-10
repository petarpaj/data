setwd("~/Dropbox/Gokcumen_Lab/Projects/Amanda_coll")

library (ggplot2)
library (gplots)
library (RColorBrewer)
library(dplyr)

### New Tajima's d


#ENTIRE DATASET - UPDATED MAY24 (FOR SWITCHING RESCUE'S TO LOF)
a <- read.csv("LOF_Dec2020_updated_may24.csv")
head(a)

#Read-depth analysis to filter out "bad" regions

rd <- read.csv("read-dept.csv")
rd <- subset (rd, rd$filter=="include")
a <- subset (a, a$SYMBOL%in%rd$gene)

unique (a$Consequence)

a <- subset (a, a$Consequence!="missense*")
a <- subset (a, a$Consequence!="non_coding_transcript_exon*")

head (a)

head (a)

#### Genotype frequencies with RESCUE's are switched
a <- read.csv("Table_LOF_variants_allele&genotype.csv")

require("ggrepel")
png ("LOF2.png",  width = 10, height = 6, units = 'in', res = 300)

df <- subset (a, a$Ugan_freq>0 | a$Phil_freq>0)
ggplot(df, aes(x=Phil_freq, y=Ugan_freq)) + 
  theme_classic()+ 
  theme(text = element_text(size=20))+
  geom_jitter(size=2,  aes(colour=factor(Consequence))) +
  scale_color_manual(values=c("#00AFBB", "light gray", "#E7B800", "#FC4E07"))+
  geom_text_repel(
    data = subset(a, Ugan_freq>0.2 | Phil_freq > 0.2),
    aes(label = SYMBOL),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
dev.off()

#00AFBB", "#E7B800", "#FC4E07"

########## All heatmap with RESCUE's are switched

#same data from above with allele sums

a <- read.csv("Table_LOF_variants_allele&genotype.csv")
mat_data <- a

row.names(mat_data) <- paste(mat_data$SYMBOL, mat_data$Type, mat_data$POS)
#get rid of singletons

mat_data <- subset (mat_data, mat_data$Ugan_freq>0.05 | mat_data$Phil_freq>0.05)
names(a)
mat_data <- as.matrix(mat_data [,16:148])


#my_palette <- rev(brewer.pal (2,"RdBu"))

png ("heatmap_all.png",  width = 10, height = 7, units = 'in', res = 300)

heatmap.2 (mat_data, density.info="none", trace="none", scale="none",
           sepwidth=c(0.0001,0.0001), col=c("light gray", "#E7B800", "#FC4E07"), sepcolor="black", colsep=1:ncol(mat_data), rowsep=1:nrow(mat_data),       
           Rowv = T, Colv=F, cexRow=0.5, cexCol=1, main="")

dev.off()



#### homozygote gene comparisions with RESCUE switched

a <- read.csv ("LOF_INDIVIDUAL_Dec2020_updated_may24.csv")

a$Population <- factor(a$Population, levels = c("Bakiga", "Basua", "Batwa", "Manobo", "Agta", "Manan"))

unique (a$Population)
p <- ggplot(a, aes(Population, Homozygous_genes))


png ("LOF1.png",  width = 10, height = 6, units = 'in', res = 300)

p + geom_boxplot(aes(fill=factor(Farmer)), alpha=0.2) + 
  scale_fill_manual(values=c( "#E7B800", "#FC4E07")) +   
  scale_color_manual(values=c("#E7B800", "#FC4E07"))  +
  theme_light()+ theme(text = element_text(size=20))+
  geom_jitter(width=0.1, size=3, aes(col=factor(Farmer)))

dev.off()


##### statistical significance - rescue swithced

a <- read.csv ("LOF_INDIVIDUAL_Dec2020_updated_may24.csv")

afr <- subset (a, a$Country=="Uganda")
asia <- subset (a, a$Country=="Phil")

mean(afr$Homozygous_genes)
mean(asia$Homozygous_genes)

stdev

wilcox.test(afr$Homozygous_genes, asia$Homozygous_genes)

hunt_afr <- subset (afr, afr$Farmer=="HG")
agr_afr <- subset (afr, afr$Farmer=="Farm")

mean(hunt_afr$Homozygous_genes)
mean(agr_afr$Homozygous_genes)

wilcox.test(agr_afr$Homozygous, hunt_afr$Homozygous)

hunt_asia <- subset (asia, asia$Farmer=="HG")
agr_asia <- subset (asia, asia$Farmer=="Farm")

mean (hunt_asia$Homozygous_genes)
mean(agr_asia$Homozygous_genes)

wilcox.test(hunt_asia$Homozygous_genes, agr_asia$Homozygous_genes)


unique (a$Subsistence)
head(a)


##### PBS DISTRIBUTION

setwd("~/Dropbox/Gokcumen_Lab/Projects/Amanda_coll")

a <- read.csv("PBS_Bakiga_basua.csv")
b <- read.csv ("PBS_Bakiga_batwa.csv")
a <- a [1:5364, ]
b <- b [1:5384,]

a$pop <- "bak_basua"
b$pop <- "bak_batwa"



df <- rbind (a, b)
df$cat <- paste (df$pop, df$Analysis_Set)

#RIDGE PLOT OF THE NEUTRAL DISTRIBUTION
library(ggridges)
head(df)


df$cat <- factor(df$cat, levels = c("bak_basua sensory_exonic",
                                    "bak_basua neutral_SNPs",
                                    "bak_batwa sensory_exonic",
                                    "bak_batwa neutral_SNPs"  ))


#pbs_agr
png ("PBS_agr.png",  width = 7, height = 5, units = 'in', res = 300)
ggplot(df, aes(x = AGR.PBS, y = cat, fill = cat)) +
  geom_density_ridges_gradient( scale = 0.5, size = 0.7, rel_min_height = 0.01,
                                jittered_points=TRUE, position = "raincloud",
                                alpha = 0.2, point_size=1, point_color="#666666",
                                quantile_lines = TRUE, quantiles = 1,  vline_size = 1, vline_color = "red"
  ) + 
  scale_fill_manual(values = (alpha(c("#7570B3", "#7570B3",  "#D95F02", "#D95F02", "#1B9E77", "#666666"), 0.5)))+
  theme_bw()

dev.off()

#pbs_hg
png ("PBS_hg.png",  width = 7, height = 5, units = 'in', res = 300)
ggplot(df, aes(x = HG.PBS, y = cat, fill = cat)) +
  geom_density_ridges_gradient( scale = 0.5, size = 0.7, rel_min_height = 0.01,
                                jittered_points=TRUE, position = "raincloud",
                                alpha = 0.2, point_size=1, point_color="#666666",
                                quantile_lines = TRUE, quantiles = 1,  vline_size = 1, vline_color = "red"
  ) + 
  scale_fill_manual(values = (alpha(c("#7570B3", "#7570B3",  "#D95F02", "#D95F02", "#1B9E77", "#666666"), 0.5)))+
  theme_bw()

dev.off()

#quantiles and specific snp information

unique (df$cat)

quantile (subset (df$HG.PBS, df$cat=="bak_basua neutral_SNPs"), 0.99)
quantile (subset (df$HG.PBS, df$cat=="bak_batwa neutral_SNPs"), 0.99)

quantile (subset (df$AGR.PBS, df$cat=="bak_basua neutral_SNPs"), 0.99)
quantile (subset (df$AGR.PBS, df$cat=="bak_batwa neutral_SNPs"), 0.99)

a <- subset (df, df$Existing_variation=="rs2961144" |
               df$Existing_variation=="rs2227264")
a




###################older stuff not used in the last version


########## LOF analysis - frequencies

a <- read.csv("Table_LOF_variants_allele&genotype.csv")
head(a)
require("ggrepel")
png ("LOF2.png",  width = 8, height = 6, units = 'in', res = 300)

ggplot(a, aes(x=Hom_Philip, y=Hom_Uganda)) + 
  theme_classic()+ 
  geom_jitter(size=2,  aes(colour=Hom_Philip/Hom_Uganda)) +
  geom_text_repel(
    data = subset(a, Hom_Philip + Hom_Uganda > 0.1),
    aes(label = SYMBOL),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
dev.off()


a <- read.csv("LOF_Dec2020_filteringForHomozygous.csv")

rd <- read.csv("read-dept.csv")
rd <- subset (rd, rd$filter=="include")
a <- subset (a, a$SYMBOL%in%rd$gene)

a<- subset (a, a$Simp_conseq!="Missense/UTR")


write.csv (a, "Table_LOF_variants.csv")



p <- ggplot(a, aes(Simp_conseq, Homozgygous))

png ("LOF1.png",  width = 6, height = 4, units = 'in', res = 300)

p + geom_boxplot(aes(fill=factor(Simp_conseq)), alpha=0.2) + 
 scale_fill_manual(values=c("#67001F", "#053061")) +   
  scale_color_manual(values=c("#67001F", "#053061"))  +
  theme_light()+
  geom_jitter(width=0.1, size=3, aes(col=factor(Simp_conseq)))

dev.off()


########## LOF analysis

a <- read.csv("LOF_Dec2020_filteringForHomozygous.csv")

a <- subset (a, a$Simp_conseq!="Missense/UTR" & a$Homozgygous>0)

require("ggrepel")
png ("LOF2.png",  width = 8, height = 6, units = 'in', res = 300)


ggplot(a, aes(x=Hom_Philip, y=Hom_Uganda)) + 
  theme_classic()+ 
  geom_jitter(size=2,  aes(colour=Simp_conseq)) +
  geom_text_repel(
    data = subset(a, Homozgygous > 15),
    aes(label = SYMBOL),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
dev.off()

########## LOF analysis with individuals

a <- read.csv("lof_homozygous.csv")

head(a)
a$Populatoin <- factor(a$Populatoin, levels = c("Bakiga", "Basua", "Batwa", "Manobo", "Agta", "Manan"))

png ("LOF2.png",  width = 10, height = 7, units = 'in', res = 300)

p <- ggplot(a, aes(Populatoin, Homozygous))

p + geom_boxplot(aes(fill=factor(Subsistence)), alpha=0.2) + 
  scale_fill_manual(values=c("#67001F", "#053061")) +   
  scale_color_manual(values=c("#67001F", "#053061"))  +
  theme_light()+
  geom_jitter(width=0.1, size=3, aes(col=factor(Subsistence)))

dev.off()

afr <- subset (a, a$Continent=="Africa")
asia <- subset (a, a$Continent=="Asia")

mean(afr$Homozygous)
mean(asia$Homozygous)

wilcox.test(afr$Homozygous, asia$Homozygous,)

hunt_afr <- subset (afr, afr$Subsistence=="Hunter-Gatherer")
agr_afr <- subset (afr, afr$Subsistence=="Farmer")

mean(hunt_afr$Homozygous)
mean(agr_afr$Homozygous)


wilcox.test(agr_afr$Homozygous, hunt_afr$Homozygous, alternative = "less")

hunt_asia <- subset (asia, asia$Subsistence=="Hunter-Gatherer")
agr_asia <- subset (asia, asia$Subsistence=="Farmer")

mean (hunt_asia$Homozygous)
mean(agr_asia$Homozygous)

wilcox.test(hunt_asia$Homozygous, agr_asia$Homozygous)


unique (a$Subsistence)
head(a)


########## INDIVIDUALS

a <- read.csv("number_individuals.csv")

head(a)
a$Populatoin <- factor(a$Populatoin, levels = c("Bakiga", "Basua", "Batwa", "Manobo", "Agta", "Manan"))

p <- ggplot(a, aes(Populatoin, Value))
p + geom_violin(aes(fill=factor(Continent)), col=NA, alpha=0.2) + 
  facet_grid(Genotype~., scales="free") +  scale_fill_manual(values=c("#67001F", "#053061")) +   
  scale_color_manual(values=c("#67001F", "#053061"))  +
  theme_light()+
  geom_jitter(width=0.1, aes(col=factor(Continent)))


########## TAJIMA'S D
setwd("~/Dropbox/Gokcumen_Lab/Projects/Amanda_coll")
df <- read.csv ("ALL_TajD.csv")
head (df)


a <- subset (df, df$N_SNPs>2 & df$Type=="SENS")

library (ggplot2)
library (RColorBrewer)

require("ggrepel")

a$Population <- factor(a$Population, levels = c("BaKiga",  "BaTwa",   "BaSua" , "Manobo" , "Mamanwa" , "Agta"   ))

png ("Tajima.png",  width = 10, height = 7, units = 'in', res = 300)

p <- ggplot(a, aes(Population, TajD))
p + geom_boxplot(aes(fill=factor(Subsistence)), col="black", alpha=0.2) + 
  theme_light()+
  geom_jitter(width=0.1, aes(col=factor(Subsistence)))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  geom_text_repel(
    data = subset(a, N_SNPs > 3 & (TajD>3.5)),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
dev.off()

a <- subset (df, df$Type=="NEU" & df$N_SNPs>2 & df$Population=="Manobo")
mean (a$TajD)
