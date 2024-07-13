###Analysis code for 
###Shechonge et al. 
###Environmental DNA metabarcoding details the spatial structure of a diverse tropical fish assemblage in a major East African river system

library(vegan)
library(ape)
library(ggplot2)
library(ggpubr)
library(DMwR2)
library(plyr)
library(ggplot2)
library(raster)
library(reshape2)
library(indicspecies)
library(tibble)
library(Hmisc)
library(lme4)
library(lmerTest)

### Testing for an effect of sampling volume on the recovered number of species
WaterFilteredTest <- read.table("Water_filtered.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
WaterFiltered_Species <- lmer(Total_species ~ Volume_water + (1|eventID), data = WaterFilteredTest)
summary(WaterFiltered_Species)

WaterFiltered_Reads <- lmer(log10(Assigned_reads) ~ Volume_water + (1|eventID), data = WaterFilteredTest)
summary(WaterFiltered_Reads)

WaterFiltered_Reads_Species <- lmer(Total_species ~ Volume_water + (1|eventID), data = WaterFilteredTest)
summary(WaterFiltered_Reads_Species)

WaterFiltered_Reads_Species <- lmer(Total_species ~ log10(Assigned_reads) + (1|eventID), data = WaterFilteredTest)
summary(WaterFiltered_Reads_Species)

### Load Ifakara samples only - 50 in total.

IfakaraData <- read.table("SpeciesData_Ifakara.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
IfakaraData <- (IfakaraData[2:67])

###Make species accumulation curve, and the dataframe for plotting.

acc <- specaccum(IfakaraData, "coleman")
data <- data.frame(Samples=acc$sites, Richness=acc$richness, SD=acc$sd)
data$min <- data$Richness - data$SD
data$max <- data$Richness + data$SD
data

###plot species accumulation curve

acc_plot <- ggplot(data, aes(x=Samples, y=Richness)) + geom_ribbon(aes(ymin=min, ymax=max), fill='lightblue') + theme_classic() + geom_point() + geom_line()
acc_plot <- acc_plot + labs(x ="Samples", y = "Number of species")
acc_plot

### Load All samples only - 174 samples in total and 66 species

SpeciesMatrix<- read.table("SpeciesData.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
HabitatList <- read.table("Habitat.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
SpeciesMatrix_Hab <- cbind(HabitatList$Habitat,SpeciesMatrix)
names(SpeciesMatrix_Hab)[1] <- "Habitat"

#Make a Hellinger matrix

SpeciesMatrixHell = decostand(x=SpeciesMatrix_Hab[,3:68], method="hellinger")
SpeciesMatrixHell <- cbind(SpeciesMatrix_Hab[,1:2],SpeciesMatrixHell)

#Principal Coordinate Analyses (PCOA) on the data matrices

SpeciesMatrixHell.D <- vegdist(SpeciesMatrixHell[,3:68], "hellinger")
PCOA_SpeciesMatrixHell<- pcoa(SpeciesMatrixHell.D, correction="none", rn=NULL)
PCOA_SpeciesMatrixHell

#Make the scores files, add habitat column

PCOAscoresHell <- PCOA_SpeciesMatrixHell$vectors
PCOAscoresHell <- as.data.frame(PCOAscoresHell)
PCOAscoresHell <- cbind(HabitatList,PCOAscoresHell)

#Plot the PCOA scores (manually enter the variance captured on axis labels)

PlotHell <- ggscatter(PCOAscoresHell, x = "Axis.1", y = "Axis.2", color = "Habitat",
                      palette = c("#999999", "#E69F00", "#56B4E9", "black", "red","#CCCC00"),
                      ellipse = FALSE, ellipse.type = "convex", mean.point = FALSE,
                      star.plot = TRUE) +
  labs(x ="PCOA1 (27.4% of variation)", y = "PCOA2 (14.7% of variation)") +
  theme(legend.position = "bottom") 
PlotHell

###Test for differences among sites, and among habitats

#Stat Test to test for differences (we will narrow this to only those wth enough samples)
GlobalHellTest <- adonis2(SpeciesMatrixHell.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = HabitatList)
GlobalHellTest

###Pairwise test of differences among habitats

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Great Ruaha River" |  Habitat =="Kilombero River (tributory high)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_GRR_KRH <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_GRR_KRH

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Great Ruaha River" |  Habitat =="Mansi", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_GRR_MAN <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_GRR_MAN

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Great Ruaha River" |  Habitat =="Kilombero River (plain)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_GRR_KRP <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_GRR_KRP

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Great Ruaha River" |  Habitat =="Kilombero River (tributory low)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_GRR_KRL <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_GRR_KRL

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Great Ruaha River" |  Habitat =="Rufiji River", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_GRR_RUR <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_GRR_RUR

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Kilombero River (tributory high)" |  Habitat =="Mansi", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_KRH_MAN <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_KRH_MAN

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Kilombero River (tributory high)" |  Habitat =="Kilombero River (plain)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_KRH_KRP <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_KRH_KRP

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Kilombero River (tributory high)" |  Habitat =="Kilombero River (tributory low)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_KRH_KRL <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_KRH_KRL

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Kilombero River (tributory high)" |  Habitat =="Rufiji River", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_KRH_RUR <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_KRH_RUR

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Mansi" |  Habitat =="Kilombero River (plain)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_MAN_KRP <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_MAN_KRP

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Mansi" |  Habitat =="Kilombero River (tributory low)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_MAN_KRL <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_MAN_KRL

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Mansi" |  Habitat =="Rufiji River", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_MAN_RUR <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_MAN_RUR

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Kilombero River (plain)" |  Habitat =="Kilombero River (tributory low)", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_KRP_KRL <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_KRP_KRL

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Kilombero River (plain)" |  Habitat =="Rufiji River", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_KRP_RUR <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_KRP_RUR

SpeciesMatrixHell_subset <- subset (SpeciesMatrixHell, Habitat == "Kilombero River (tributory low)" |  Habitat =="Rufiji River", select=c(1:68))
SpeciesMatrixHell_subset.D <- vegdist(SpeciesMatrixHell_subset[,3:68], "hellinger")
HellTest_KRL_RUR <- adonis2(SpeciesMatrixHell_subset.D ~ Habitat/eventID, strata=NULL, permutations = 10000, data = SpeciesMatrixHell_subset)
HellTest_KRL_RUR

###RDA to identify the key environmental variables that predict differences in communities
###Read in the Environmental File, which includes NAs for missing data

Environmental <- read.table("Environmental.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

Environmental_corr_1<-rcorr(as.matrix(Environmental[,4:11]))
Environmental_corr_1

Environmental_Imputed <- knnImputation(Environmental[4:11])
Environmental_Site <- cbind(Environmental[1:1], Environmental_Imputed)

Environmental_corr_2<-rcorr(as.matrix(Environmental_Site[,2:9]))
Environmental_corr_2

Environmental_Site_174 <- merge(Environmental_Site,SpeciesMatrixHell,by="eventID")

### Execute the RDA, for each variable independently.

elev_rda = rda(Environmental_Site_174[11:76] ~ Elevation, data = Environmental_Site_174)
summary(elev_rda)
anova(elev_rda, permutations = 10000)

DO_rda = rda(Environmental_Site_174[11:76] ~ DO, data = Environmental_Site_174)
summary(DO_rda)
anova(DO_rda, permutations = 10000)

Flow_rda = rda(Environmental_Site_174[11:76] ~ Still_Flowing, data = Environmental_Site_174)
summary(Flow_rda)
anova(Flow_rda, permutations = 10000)

Cond_rda = rda(Environmental_Site_174[11:76] ~ Conductivity, data = Environmental_Site_174)
summary(Cond_rda)
anova(Cond_rda, permutations = 10000)

pH_rda = rda(Environmental_Site_174[11:76] ~ pH, data = Environmental_Site_174)
summary(pH_rda)
anova(pH_rda, permutations = 10000)

sal_rda = rda(Environmental_Site_174[11:76] ~ Salinity, data = Environmental_Site_174)
summary(sal_rda)
anova(sal_rda, permutations = 10000)

temp_rda = rda(Environmental_Site_174[11:76] ~ Temperature, data = Environmental_Site_174)
summary(temp_rda)
anova(temp_rda, permutations = 10000)

turb_rda = rda(Environmental_Site_174[11:76] ~ Temperature, data = Environmental_Site_174)
summary(turb_rda)
anova(turb_rda, permutations = 10000)

####Plot PCOA1 x Elevation

PCOAscoresHell <- cbind(PCOAscoresHell,Environmental_Site_174[2:9])
PlotHell_Elev <- ggscatter(PCOAscoresHell, x = "Elevation", y = "Axis.1", color = "Habitat",
                           palette = c("#999999", "#E69F00", "#56B4E9", "black", "red","#CCCC00"),
                           ellipse = FALSE, ellipse.type = "convex", mean.point = FALSE,
                           star.plot = FALSE) +
  labs(x ="Elevation (m)", y = "PCOA1 (27.4% of variation)") +
  theme(legend.position = "bottom") 
PlotHell_Elev

####Plot PCOA1 x Temperature

PlotHell_Temp <- ggscatter(PCOAscoresHell, x = "Temperature", y = "Axis.1", color = "Habitat",
                           palette = c("#999999", "#E69F00", "#56B4E9", "black", "red","#CCCC00"),
                           ellipse = FALSE, ellipse.type = "convex", mean.point = FALSE,
                           star.plot = FALSE) +
  labs(x ="Temperature (°C)", y = "PCOA1 (27.4% of variation)") +
  theme(legend.position = "bottom") 
PlotHell_Temp

####Plot aggregate figure

ggarrange(acc_plot,PlotHell,PlotHell_Elev,PlotHell_Temp,
          labels = c("a", "b", "c","d"), align = c("hv"), font.label = list(size = 14),
          ncol = 2, nrow = 2, legend = "bottom",  common.legend = TRUE)

#save as 9 x 7

##Indicator Species (core analysis)

Abund <- SpeciesMatrix[,2:ncol(SpeciesMatrix)]
Group <- PCOAscoresHell$Habitat
indval <- multipatt(Abund, Group, control = how(nperm=999)) 

##Indicator Species (make table if IndVal scores)

IndVal_Index <- indval$A*indval$B*100
IndVal_Index <- IndVal_Index[,1:6]
IndVal_Index <- as.data.frame(IndVal_Index)
IndVal_Index <- tibble::rownames_to_column(IndVal_Index, "Species")

##Indicator Species (make plot of IndVal scores)

IndVal_Index_Long = melt(IndVal_Index , id = c("Species"))
colnames(IndVal_Index_Long)[2] <- "Habitat"
colnames(IndVal_Index_Long)[3] <- "IndVal"

#Pick your colours:
  
  colours = c("#999999", "#E69F00", "#56B4E9", "black", "red","#CCCC00")
  

##Bubble plot, use geom_point and scale size to relative abundance

level_order <- c('Amphilius chalei',
                 'Amphilius crassus',
                 'Anguilla bengalensis',
                 'Anguilla mossambica',
                 'Chiloglanis sp. mbarali river',
                 'Engraulicypris brevianalis',
                 'Enteromius atkinsoni',
                 'Enteromius kerstenii',
                 'Labeo cylindricus',
                 'Mastacembelus frenatus',
                 'Parakneria tanzaniae',
                 'Zaireichthys sp. ruaha',
                 'Atopochilus vogti',
                 'Chiloglanis sp. mzombe river sp.A',
                 'Clarias gariepinus',
                 'Clarias theodorae',
                 'Enteromius apleurogramma',
                 'Enteromius luikae',
                 'Enteromius radiatus',
                 'Enteromius zanzibaricus',
                 'Haplochromis vanheusdeni',
                 'Labeobarbus macrolepis',
                 'Labeobarbus oxyrhynchus',
                 'Opsaridium loveridgii',
                 'Oreochromis leucostictus',
                 'Bagrus orientalis',
                 'Brycinus lateralis',
                 'Ctenopoma muriei',
                 'Cyphomyrus discorhynchus',
                 'Enteromius lineomaculatus',
                 'Enteromius macrotaenia',
                 'Hydrocynus tanzaniae',
                 'Lacustricola kongoranensis',
                 'Marcusenius livingstonii',
                 'Mormyrus longirostris',
                 'Nothobranchius kilomberoensis',
                 'Pareutropius longifilis',
                 'Pollimyrus castelnaui',
                 'Pseudocrenilabrus sp. kilombero',
                 'Hemigrammopetersius barnardi',
                 'Schilbe moebiusii',
                 'Synodontis sp. utete',
                 'Anguilla marmorata',
                 'Astatotilapia sp. ruaha blue',
                 'Astatotilapia sp. ruaha redcheek',
                 'Oreochromis sp. mtera',
                 'Alestes stuhlmannii',
                 'Anguilla bicolor',
                 'Brycinus affinis',
                 'Brycinus imberi',
                 'Citharinus congicus',
                 'Distichodus petersii',
                 'Eleotris klunzingerii',
                 'Glossogobius giuris',
                 'Labeo congoro',
                 'Megalops cyprinoides',
                 'Petrocephalus catostoma',
                 'Synodontis rufigiensis',
                 'Synodontis rukwaensis',
                 'Astatotilapia gigliolii',
                 'Astatotilapia sp. rufiji blue',
                 'Coptodon rendalli',
                 'Marcusenius macrolepidotus',
                 'Nannaethiops sp. kilombero',
                 'Oreochromis niloticus',
                 'Oreochromis urolepis')

level_order2 <- c('Kilombero River (tributory high)',
                  'Kilombero River (tributory low)',
                  'Kilombero River (plain)',
                  'Great Ruaha River',
                  'Rufiji River',
                  'Mansi')

IndValPlot = ggplot(IndVal_Index_Long, aes(x = factor(Habitat, level = level_order2), y = factor(Species, level = level_order))) + 
  geom_point(aes(size = IndVal, fill = Habitat), alpha = 0.7, shape = 21) +
  scale_size_continuous(limits = c(0.000000001, 100)) + scale_y_discrete(limits = rev) +
  labs( x= "", y = "", size = "Indicator Value", fill = "") +
  scale_fill_manual(values = colours, guide = "none") + 
  theme_classic() + theme(legend.position='right') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
IndValPlot

#Save as 13 x 7, missing values are zeros, and they should be missing.

#Associating geographic coordinates with sample read counts for the 49 sites
SpeciesMatrix49 <- SpeciesMatrix[(1:67)]
library(data.table)
library(dplyr)
DT <- as.data.table(SpeciesMatrix49 )
SpeciesMatrix49 <- DT[, lapply(.SD,sum), by = "eventID"]
eventID_Lat_Long <- read.table("eventID_Lat_Long.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
SpeciesMatrix49 <- left_join(eventID_Lat_Long, SpeciesMatrix49, by = join_by(eventID == eventID))
write.table (SpeciesMatrix49, file="SpeciesMatrix49.txt", sep = ";")

#End of code

#Core Rare Analsus

library(dplyr)

CoreRare_Ifakara <- as.data.frame(t(IfakaraData))
CoreRare_Ifakara$AverageReads <- rowMeans(CoreRare_Ifakara[,1:50])
CoreRare_Ifakara$NumberSamples <-rowSums(CoreRare_Ifakara[,1:50] > 0)
CoreRare_Ifakara <- subset(CoreRare_Ifakara, !CoreRare_Ifakara$NumberSamples==0)
fitmodel <- lm(log10(AverageReads) ~ poly(NumberSamples, 3, raw=TRUE),data=CoreRare_Ifakara)
summary(fitmodel)

Plot <- ggplot(CoreRare_Ifakara, aes(NumberSamples,log10(AverageReads))) + 
  stat_smooth(method="lm",
              formula=y ~ poly(x, 3, raw=TRUE),colour="black", se=FALSE)+
  geom_point() +
  theme_classic()+
  labs(x ="Number of Samples", y = "Mean reads per sample (log10 transformed)")
Plot

CoreRareBlank<- read.table("Core_Rare_Blank.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
CoreRareBlank$Prediction<- predict(fitmodel,newdata=CoreRareBlank)
CoreRareBlank <- CoreRareBlank %>%mutate(Calc = Prediction - lag(Prediction))

#Output theinflexion point
CoreRareBlank$NumberSamples[which.min(CoreRareBlank$Calc[1:50])]

#AssignCoreRare and clean up
CoreRare_Ifakara$Group = as.numeric(CoreRare_Ifakara$NumberSamples)
CoreRare_Ifakara$Group <- ifelse(CoreRare_Ifakara$NumberSamples>=32,"Core","Rare")
CoreRare_Ifakara$Species <- rownames(CoreRare_Ifakara)
CoreRare_Ifakara <- CoreRare_Ifakara[51:54]

#

