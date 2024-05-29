################################################################################
##### MOSS PROJECT SCRIPT #####
################################################################################

#usethis::use_github()

##### Data manipulation #####
pacman::p_load(openxlsx,tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,
               phytools,ggtree, treeio,ggimage, broom, ggally)


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

#### MCMCglmm z - Rate complete ~ substrate ####

Moss_ID_data$z_Rate_complete <- scale(Moss_ID_data$Rate_complete)
inv.mossphylo<-inverseA(Moss_tree,nodes="TIPS",scale=TRUE)

p1 = list(B=list(mu=rep(0,4), V=diag(4)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m1a<-MCMCglmm(z_Rate_complete ~ Substrate, random = ~Genus,
             ginverse=list(Genus=inv.mossphylo$Ainv),
             family ="gaussian", data = Moss_ID_data,
             prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
m1b<-MCMCglmm(z_Rate_complete ~ Substrate, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
m1c<-MCMCglmm(z_Rate_complete ~ Substrate, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))
#And random effects
plot(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))

gelman.diag(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))

summary(m1a)

boxplot(Moss_ID_data$z_Rate_complete~Moss_ID_data$Substrate)
boxplot(Moss_ID_data$z_Rate_separate~Moss_ID_data$Substrate)
#### MCMCglmm z - Rate separate ~ substrate  ####

Moss_ID_data$z_Rate_separate <- scale(Moss_ID_data$Rate_separate)
inv.mossphylo<-inverseA(Moss_tree,nodes="TIPS",scale=TRUE)

p2 = list(B=list(mu=rep(0,4), V=diag(4)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m2a<-MCMCglmm(z_Rate_separate ~ Substrate, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p2, nitt=110000, burnin=10000, thin=100,verbose=F)
m2b<-MCMCglmm(z_Rate_separate ~ Substrate, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p2, nitt=110000, burnin=10000, thin=100,verbose=F)
m2c<-MCMCglmm(z_Rate_separate ~ Substrate, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p2, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(m2a$Sol,m2b$Sol,m2c$Sol))
#And random effects
plot(mcmc.list(m2a$Sol,m2b$Sol,m2c$Sol))

gelman.diag(mcmc.list(m2a$Sol,m2b$Sol,m2c$Sol))

summary(m2a)


#### MCMCglmm z - Rate complete ~ Environment ####

Moss_ID_data$z_Rate_complete <- scale(Moss_ID_data$Rate_complete)
inv.mossphylo<-inverseA(Moss_tree,nodes="TIPS",scale=TRUE)

p3 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m3a<-MCMCglmm(z_Rate_complete ~ Environment, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p3, nitt=110000, burnin=10000, thin=100,verbose=F)
m3b<-MCMCglmm(z_Rate_complete ~ Environment, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p3, nitt=110000, burnin=10000, thin=100,verbose=F)
m3c<-MCMCglmm(z_Rate_complete ~ Environment, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p3, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(m3a$Sol,m3b$Sol,m3c$Sol))
#And random effects
plot(mcmc.list(m3a$Sol,m3b$Sol,m3c$Sol))

gelman.diag(mcmc.list(m3a$Sol,m3b$Sol,m3c$Sol))

summary(m3a)

#### MCMCglmm z - Rate separate ~ Environment  ####

Moss_ID_data$z_Rate_separate <- scale(Moss_ID_data$Rate_separate)
inv.mossphylo<-inverseA(Moss_tree,nodes="TIPS",scale=TRUE)

p4 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m4a<-MCMCglmm(z_Rate_separate ~ Environment, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p4, nitt=110000, burnin=10000, thin=100,verbose=F)
m4b<-MCMCglmm(z_Rate_separate ~ Environment, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p4, nitt=110000, burnin=10000, thin=100,verbose=F)
m4c<-MCMCglmm(z_Rate_separate ~ Environment, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p4, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(m4a$Sol,m4b$Sol,m4c$Sol))
#And random effects
plot(mcmc.list(m4a$Sol,m4b$Sol,m4c$Sol))

gelman.diag(mcmc.list(m4a$Sol,m4b$Sol,m4c$Sol))

summary(m4a)
#### Data with means per Genus ####

Moss_Genus_data <- Moss_ID_data %>% 
  group_by(Genus) %>% 
  summarise(mean_Rate_complete=mean(Rate_complete),
            mean_Rate_separate=mean(Rate_separate),
            mean_Immediate_diff_complete=mean(Immediate_diff_complete),
            mean_Final_diff_complete=mean(Final_diff_complete),
            mean_Immediate_diff_separated=mean(Immediate_diff_separated),
            mean_Final_diff_separated=mean(Final_diff_separated))

Moss_Genus_data <- as.data.frame(Moss_Genus_data)

#### MCMCglmm Ancestral reconstruction rate_complete ####

inv.mossphylo_all<-inverseA(Moss_tree,nodes="ALL",scale=TRUE)
# set priors
p5 = list(B=list(mu=rep(0,1), V=diag(1)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))
#model
m5a<-MCMCglmm(mean_Rate_complete ~ 1, random = ~Genus,
             ginverse=list(Genus=inv.mossphylo_all$Ainv),
             family ="gaussian", data = Moss_Genus_data, pr=TRUE,
             prior=p5, nitt=110000, burnin=10000, thin=100,verbose=F)
m5b<-MCMCglmm(mean_Rate_complete ~ 1, random = ~Genus,
             ginverse=list(Genus=inv.mossphylo_all$Ainv),
             family ="gaussian", data = Moss_Genus_data, pr=TRUE,
             prior=p5, nitt=110000, burnin=10000, thin=100,verbose=F)
m5c<-MCMCglmm(mean_Rate_complete ~ 1, random = ~Genus,
             ginverse=list(Genus=inv.mossphylo_all$Ainv),
             family ="gaussian", data = Moss_Genus_data, pr=TRUE,
             prior=p5, nitt=110000, burnin=10000, thin=100,verbose=F)

#This are the logit not transformed back of the values for each family 
#and each node
posterior.mode(m5a$Sol)

#Plot the fixed effects
plot(mcmc.list(m5a$Sol,m5b$Sol,m5c$Sol))
#And random effects
plot(mcmc.list(m5a$Sol,m5b$Sol,m5c$Sol))

gelman.diag(mcmc.list(m5a$Sol,m5b$Sol,m5c$Sol))

summary(m5a)

#BLUPS
#To get estimates of ancestral states add intercept to node BLUPs
blupsm5a<-data.frame(effect=colnames(m5a$Sol), estimate=posterior.mode(
  m5a$Sol[,'(Intercept)'])+posterior.mode(m5a$Sol),CI=HPDinterval(
    m5a$Sol[,'(Intercept)']+m5a$Sol))

#We added the intercept value to all values and so "Node1" value is incorrect 
#- lets replace it
blupsm5a['(Intercept)','estimate']<-posterior.mode(m5a$Sol[,'(Intercept)'])
blupsm5a['(Intercept)',c('CI.lower','CI.upper')]<-HPDinterval(
  m5a$Sol[,'(Intercept)'])

#The intercept is the value of the root ("Node1") so lets rename it
blupsm5a$effect<-as.character(blupsm5a$effect)
blupsm5a['(Intercept)','effect']<-"Node1"

#and remove the text "treetip." from the effect
blupsm5a$effect<-gsub("Genus.","",blupsm5a$effect)

#match with phylogeny
#match the vitamin data and the insect tree data
mean_Rate_complete_phylo<-(Moss_Genus_data$mean_Rate_complete)[match(
  Moss_tree$tip.label,Moss_Genus_data$Genus)]

#plot them and add the node numbers to being able to do 4.10
plot(Moss_tree, cex=0.8, no.margin =T, label.offset = 0.7)
nodelabels(pch=21, cex=(abs(blupsm5a$estimate[1:Moss_tree$Nnode])*100))
tiplabels(pch=21,cex=mean_Rate_complete_phylo*100,bg="black")


#### check phylogenetic signal heredability rate_complete ####

m5aPhyloSig<-m5a$VCV[,'Genus']/((m5a$VCV[,'Genus']+m5a$VCV[,'units']))

#posterior.mode is the heredability value in this case
posterior.mode(m5aPhyloSig)
HPDinterval(m5aPhyloSig)

plot(m5aPhyloSig)

#### MCMCglmm Ancestral reconstruction rate_separate ####
inv.mossphylo_all<-inverseA(Moss_tree,nodes="ALL",scale=TRUE)
# set priors
p6 = list(B=list(mu=rep(0,1), V=diag(1)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))
#model
m6a<-MCMCglmm(mean_Rate_separate ~ 1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo_all$Ainv),
              family ="gaussian", data = Moss_Genus_data, pr=TRUE,
              prior=p6, nitt=110000, burnin=10000, thin=100,verbose=F)
m6b<-MCMCglmm(mean_Rate_separate ~ 1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo_all$Ainv),
              family ="gaussian", data = Moss_Genus_data, pr=TRUE,
              prior=p6, nitt=110000, burnin=10000, thin=100,verbose=F)
m6c<-MCMCglmm(mean_Rate_separate ~ 1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo_all$Ainv),
              family ="gaussian", data = Moss_Genus_data, pr=TRUE,
              prior=p6, nitt=110000, burnin=10000, thin=100,verbose=F)

#This are the logit not transformed back of the values for each family 
#and each node
posterior.mode(m6a$Sol)

#Plot the fixed effects
plot(mcmc.list(m6a$Sol,m6b$Sol,m6c$Sol))
#And random effects
plot(mcmc.list(m6a$Sol,m6b$Sol,m6c$Sol))

gelman.diag(mcmc.list(m6a$Sol,m6b$Sol,m6c$Sol))

summary(m6a)

#BLUPS
#To get estimates of ancestral states add intercept to node BLUPs
blupsm6a<-data.frame(effect=colnames(m6a$Sol), estimate=posterior.mode(
  m6a$Sol[,'(Intercept)'])+posterior.mode(m6a$Sol),CI=HPDinterval(
    m6a$Sol[,'(Intercept)']+m6a$Sol))

#We added the intercept value to all values and so "Node1" value is incorrect 
#- lets replace it
blupsm6a['(Intercept)','estimate']<-posterior.mode(m6a$Sol[,'(Intercept)'])
blupsm6a['(Intercept)',c('CI.lower','CI.upper')]<-HPDinterval(
  m6a$Sol[,'(Intercept)'])

#The intercept is the value of the root ("Node1") so lets rename it
blupsm6a$effect<-as.character(blupsm6a$effect)
blupsm6a['(Intercept)','effect']<-"Node1"

#and remove the text "treetip." from the effect
blupsm6a$effect<-gsub("Genus.","",blupsm6a$effect)

#match with phylogeny
#match the vitamin data and the insect tree data
mean_Rate_separate_phylo<-(Moss_Genus_data$mean_Rate_separate)[match(
  Moss_tree$tip.label,Moss_Genus_data$Genus)]

#plot them and add the node numbers to being able to do 4.10
plot(Moss_tree, cex=0.8, no.margin =T, label.offset = 0.7)
nodelabels(pch=21, cex=(abs(blupsm6a$estimate[1:Moss_tree$Nnode])*100))
tiplabels(pch=21,cex=mean_Rate_separate_phylo*100,bg="black")

#### check phylogenetic signal heredability rate_separate ####

m6aPhyloSig<-m6a$VCV[,'Genus']/((m6a$VCV[,'Genus']+m6a$VCV[,'units']))

#posterior.mode is the heredability value in this case
posterior.mode(m6aPhyloSig)
HPDinterval(m6aPhyloSig)

plot(m6aPhyloSig)

#### MCMCglmm Ancestral reconstruction substrate ####
inv.mossphylo_all<-inverseA(Moss_tree,nodes="ALL",scale=TRUE)
# set priors
p7 = list(B=list(mu=rep(0,1), V=diag(1)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

#model
m7<-MCMCglmm(Substrate ~ 1, random = ~Genus,
             ginverse=list(Genus=inv.mossphylo_all$Ainv),
             family ="categorical", data = Moss_ID_data, pr=TRUE, 
             rcov = ~us(trait):units,
             prior=p7, nitt=110000, burnin=10000, thin=100,verbose=F)

#This are the logit not transformed back of the values for each family 
#and each node
posterior.mode(m7$Sol)

#BLUPS
#To get estimates of ancestral states add intercept to node BLUPs
blupsm7<-data.frame(effect=colnames(m7$Sol), estimate=posterior.mode(
  m7$Sol[,'(Intercept)'])+posterior.mode(m7$Sol),CI=HPDinterval(
    m7$Sol[,'(Intercept)']+m7$Sol))

#We added the intercept value to all values and so "Node1" value is incorrect 
#- lets replace it
blupsm7['(Intercept)','estimate']<-posterior.mode(m7$Sol[,'(Intercept)'])
blupsm7['(Intercept)',c('CI.lower','CI.upper')]<-HPDinterval(
  m7$Sol[,'(Intercept)'])

#The intercept is the value of the root ("Node1") so lets rename it
blupsm7$effect<-as.character(blupsm7$effect)
blupsm7['(Intercept)','effect']<-"Node1"

#and remove the text "treetip." from the effect
blupsm7$effect<-gsub("Genus.","",blupsm7$effect)

#match with phylogeny
#match the vitamin data and the insect tree data
Rate_complete_phylo<-(Moss_ID_data$Rate_separate)[match(
  Moss_tree$tip.label,Moss_ID_data$Genus)]

#plot them and add the node numbers to being able to do 4.10
plot(Moss_tree, cex=0.8, no.margin =T, label.offset = 0.7)
nodelabels(pch=21, cex=(abs(blupsm7$estimate[1:Moss_tree$Nnode])*200))
tiplabels(pch=21,cex=Rate_complete_phylo*200,bg="black")