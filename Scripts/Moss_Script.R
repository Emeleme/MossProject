################################################################################
##### MOSS PROJECT SCRIPT #####
################################################################################

#usethis::use_github()

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

#Delete node.labels 
Moss_tree$node.label <- NULL

plot(Moss_tree)
nodelabels()

#### Calculating rates ID ####
#We use the rates from min 30 to min 120 because 0-30min is water dropping
rate_times <- c(30, 60, 90, 120)

Filter_complete <- Moss_complete %>%
  filter(Time %in% rate_times)
Filter_separate <- Moss_separate %>%
  filter(Time %in% rate_times)
  
Moss_rates_complete_ID <- Filter_complete %>%
  group_by(ID, Genus) %>%
  do(tidy(lm(Weight ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  select(Genus, estimate)
  
Moss_rates_separate_ID <- Filter_separate %>%
  group_by(ID, Genus) %>%
  do(tidy(lm(Weight ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  select(Genus, estimate)

Moss_rates_ID <- merge(Moss_rates_complete_ID, Moss_rates_separate_ID, 
                       by=c("Genus", "ID"))
colnames(Moss_rates_ID) <- c("Genus", "ID", "Rate_complete", "Rate_separate")

Moss_rates_ID$Rate_complete<- abs(Moss_rates_ID$Rate_complete)
Moss_rates_ID$Rate_separate<- abs(Moss_rates_ID$Rate_separate)

#### Calculating rates _Genus ####
#We use the rates from min 30 to min 120 because 0-30min is water dropping
Moss_rates_complete_Genus <- Filter_complete %>%
  group_by(Genus) %>%
  do(tidy(lm(Weight ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  select(Genus, estimate)

Moss_rates_separate_Genus <- Filter_separate %>%
  group_by(Genus) %>%
  do(tidy(lm(Weight ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  select(Genus, estimate)

Moss_rates_Genus <- merge(Moss_rates_complete_Genus, Moss_rates_separate_Genus, 
                       by="Genus")
colnames(Moss_rates_Genus) <- c("Genus", "Rate_complete", "Rate_separate")

Moss_rates_Genus$Rate_complete<- abs(Moss_rates_Genus$Rate_complete)
Moss_rates_Genus$Rate_separate<- abs(Moss_rates_Genus$Rate_separate)

#### Calculate the differences between initial and end wetness ID ####
diff_times_immediate <- c(0, 30)

Filter_complete_i <- Moss_complete %>%
  filter(Time %in% diff_times_immediate)
Filter_separate_i <- Moss_separate %>%
  filter(Time %in% diff_times_immediate)

Diff_complete_i <- Filter_complete_i %>%
  group_by(ID) %>%
  spread(Time, Weight) %>%  
  mutate(Immediate_diff_complete = `0` - `30`) 

Diff_separate_i <- Filter_separate_i %>%
  group_by(ID) %>%
  spread(Time, Weight) %>%  
  mutate(Immediate_diff_separated = `0` - `30`)

diff_times_final <- c(30, 660)

Filter_complete_f <- Moss_complete %>%
  filter(Time %in% diff_times_final)
Filter_separate_f <- Moss_separate %>%
  filter(Time %in% diff_times_final)

Diff_complete_f <- Filter_complete_f %>%
  group_by(ID) %>%
  spread(Time, Weight) %>%  
  mutate(Final_diff_complete = `30` - `660`) 

Diff_separate_f <- Filter_separate_f %>%
  group_by(ID) %>%
  spread(Time, Weight) %>%  
  mutate(Final_diff_separated = `30` - `660`)

Moss_diff_c_ID <- merge(Diff_complete_i,Diff_complete_f, by=c("ID", "Genus"))
Moss_diff_c_ID <- select(Moss_diff_c_ID, -c("0","30.x","30.y","660"))

Moss_diff_s_ID <- merge(Diff_separate_i,Diff_separate_f, by=c("ID", "Genus"))
Moss_diff_s_ID <- select(Moss_diff_s_ID, -c("0","30.x","30.y","660"))

#### Calculate the differences between initial and end wetness _Genus ####

Moss_diff_c_Genus <- Moss_diff_c_ID %>% 
  group_by(Genus) %>% 
  summarise(mean_imm_diff=mean(Immediate_diff_complete),
            mean_final_diff=mean(Final_diff_complete))

Moss_diff_s_Genus <- Moss_diff_s_ID %>% 
  group_by(Genus) %>% 
  summarise(mean_imm_diff=mean(Immediate_diff_separated),
            mean_final_diff=mean(Final_diff_separated))

#### Make the complete dataframe for all values ####

Moss_ID_data <- reduce(list(Moss_data, Moss_rates_ID, Moss_diff_c_ID, Moss_diff_s_ID), 
                       left_join, by = c("ID", "Genus"))

#### Calculate modes of evolution Genus ####

MossRC_label<-as.numeric(Moss_rates_Genus$Rate_complete[match(
  Moss_tree$tip.label,Moss_rates_Genus$Genus)])
names(MossRC_label)<-Moss_tree$tip.label
MossSC_label<-as.numeric(Moss_rates_Genus$Rate_separate[match(
  Moss_tree$tip.label,Moss_rates_Genus$Genus)])
names(MossSC_label)<-Moss_tree$tip.label

BM_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "BM")
OU_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "OU", 
                    bounds=list(alpha=c(0,10000)))
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

#### MCMCglmm ####

inv.mossphylo<-inverseA(Moss_tree,nodes="TIPS",scale=TRUE)

p1 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m1<-MCMCglmm(Rate_complete ~ Substrate, random = ~Genus, 
             ginverse=list(Genus=inv.mossphylo$Ainv), 
             family ="poisson", data = Moss_ID_data, 
             prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
summary(m1)

