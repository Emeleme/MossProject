################################################################################
##### MOSS PROJECT SCRIPT #####
################################################################################

##### Data manipulation #####
pacman::p_load(openxlsx,tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,
               phytools,ggtree, treeio,ggimage, broom)


#### Uploading data ####
Moss_data<-read.xlsx("Data/Moss_Data_Full.xlsx", sheet="Dry_data")
Moss_complete<-read.xlsx("Data/Moss_Data_Full.xlsx", sheet="Complete") 
Moss_separate<-read.xlsx("Data/Moss_Data_Full.xlsx", sheet="Separate") 
Moss_tree_raw<- read.tree("Data/Concat.treefile")
Moss_tips <- read.csv("Data/tips.txt", h=F)


#### Tree prunning ####
#Are all species in the tips file in the tree?
table(Moss_tips$V1 %in% Moss_tree_raw$tip.label)

#Are all tree tips in the dataset?
table(Moss_tree_raw$tip.label %in% Moss_tips$V1)

#Which are missing?
missingtips<-Moss_tree_raw$tip.label[Moss_tree_raw$tip.label %in% Moss_tips$V1 == "FALSE"]

#Remove these species from the tree using drop.tip
Moss_tree<-drop.tip(Moss_tree_raw, missingtips)

#Remove all the indicators from the sequences to leave only the genus name
Moss_tree$tip.label <- sub(".*\\.", "", Moss_tree$tip.label)
Moss_tree$tip.label <- sub("_.*", "", Moss_tree$tip.label)

#Convert to ultrametric tree
Moss_tree<- chronos(Moss_tree)

plot(Moss_tree)

#### Calculating rates ####
#We use the rates from min 30 to min 120 because 0-30min is water dropping
rate_times <- c(30, 60, 90, 120)

Filter_complete <- Moss_complete %>%
  filter(Time %in% rate_times)
Filter_separate <- Moss_separate %>%
  filter(Time %in% rate_times)
  
Moss_rates_complete <- Filter_complete %>%
  group_by(Genus) %>%
  do(tidy(lm(Weight ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  select(Genus, estimate)
  
Moss_rates_separate <- Filter_separate %>%
  group_by(Genus) %>%
  do(tidy(lm(Weight ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  select(Genus, estimate)

Moss_rates <- merge(Moss_rates_complete, Moss_rates_separate, by="Genus")
colnames(Moss_rates) <- c("Genus", "Rate_complete", "Rate_separate")

Moss_rates$Rate_complete<- abs(Moss_rates$Rate_complete)
Moss_rates$Rate_separate<- abs(Moss_rates$Rate_separate)

#### Calculate modes of evolution ####

MossRC_label<-as.numeric(Moss_rates$Rate_complete[match(Moss_tree$tip.label,Moss_rates$Genus)])
names(MossRC_label)<-Moss_tree$tip.label
MossSC_label<-as.numeric(Moss_rates$Rate_separate[match(Moss_tree$tip.label,Moss_rates$Genus)])
names(MossSC_label)<-Moss_tree$tip.label

BM_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "BM")
OU_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "OU", bounds=list(alpha=c(0,10000)))
EB_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "EB")

BM_C$opt$aic
OU_C$opt$aic
EB_C$opt$aic

##BM versus OU
BMvsOU_C <- -2*(BM_C$opt$lnL - OU_C$opt$lnL)
BMvsOU_C
pchisq(BMvsOU_C, df=1,lower.tail = FALSE)

##BM versus EB
BMvsEB_C <- -2*(EB_C$opt$lnL - BM_C$opt$lnL)
BMvsEB_C
pchisq(BMvsEB_C, df=1,lower.tail = FALSE)

##OU versus EB
OUvsEB_C <- -2*(EB_C$opt$lnL - OU_C$opt$lnL)
OUvsEB_C 
pchisq(OUvsEB_C, df=1,lower.tail = FALSE)

#OU is better than both EB and BM. The values are restrained to ne value in the model
#Use alpha when doing the phylogenetic approaches: 1387.735201
