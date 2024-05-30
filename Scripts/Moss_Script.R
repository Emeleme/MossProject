################################################################################
##### MOSS PROJECT SCRIPT #####
#Author: Maria Laura Mahecha Escobar
#Organization: Lund University
#Course: Evolution Methods and applications
#Date: 2024a05m30d
#Project description: 1. Are moss species (Bryophyta) growing in dry substrates
#                     better at water retention?
#                     2. Is there a phylogenetic signal in water retention of 
#                     mosses?
################################################################################

#usethis::use_github()

##### Data manipulation #####
pacman::p_load(openxlsx,tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,
               phytools,ggtree, treeio,ggimage, broom, ggally, wesanderson,
               paletteer, coda)


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

Moss_rates_ID$Rate_complete<- abs((Moss_rates_ID$Rate_complete)*1440)
Moss_rates_ID$Rate_separate<- abs((Moss_rates_ID$Rate_separate)*1440)

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

# Moss_diff_c_Genus <- Moss_diff_c_ID %>% 
#   group_by(Genus) %>% 
#   summarise(mean_imm_diff=mean(Immediate_diff_complete),
#             mean_final_diff=mean(Final_diff_complete))
# 
# Moss_diff_s_Genus <- Moss_diff_s_ID %>% 
#   group_by(Genus) %>% 
#   summarise(mean_imm_diff=mean(Immediate_diff_separated),
#             mean_final_diff=mean(Final_diff_separated))

#### Make the complete dataframe for all values ####

Moss_ID_data <- reduce(list(Moss_data, Moss_rates_ID, Moss_diff_c_ID, Moss_diff_s_ID), 
                       left_join, by = c("ID", "Genus"))

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


################################################################################
#### Calculate modes of evolution Genus ####
 
#Match tips with rate complete
MossRC_label<-as.numeric(Moss_Genus_data$mean_Rate_complete[match(
 Moss_tree$tip.label,Moss_Genus_data$Genus)])
names(MossRC_label)<-Moss_tree$tip.label

#Match tips with rate separate
MossSC_label<-as.numeric(Moss_Genus_data$mean_Rate_separate[match(
  Moss_tree$tip.label,Moss_Genus_data$Genus)])
names(MossSC_label)<-Moss_tree$tip.label

#Match tips with immediate loss complete
MossImC_label<-as.numeric(Moss_Genus_data$mean_Immediate_diff_complete[match(
  Moss_tree$tip.label,Moss_Genus_data$Genus)])
names(MossImC_label)<-Moss_tree$tip.label

#Match tips with immediate loss separate
MossImS_label<-as.numeric(Moss_Genus_data$mean_Immediate_diff_separated[match(
  Moss_tree$tip.label,Moss_Genus_data$Genus)])
names(MossImS_label)<-Moss_tree$tip.label

#### Compare modes of evolution for rates and immediate loss ####
#Rates of evolution function
modes_evol <- function(x) {  
  # Fit the models
  BM <- fitContinuous(phy = Moss_tree, dat = x, 
                      model = "BM")
  cat(paste0("Brownian motion (BM) sigma^2 value is: ", BM$opt$sigsq, "\n",
             ",AIC value is: ", BM$opt$aic, "\n",
             ",lnL is: ", BM$opt$lnL, "\n"))
  
  OU <- fitContinuous(phy = Moss_tree, dat = x, 
                      model = "OU")
  cat(paste0("Ornstein-Uhlenbeck (OU) alpha value is: ", OU$opt$alpha, "\n",
             ",AIC value is: ", OU$opt$aic, "\n",
             ",lnL is: ", OU$opt$lnL, "\n"))
  
  EB <- fitContinuous(phy = Moss_tree, dat = x, 
                      model = "EB")
  cat(paste0("Early burst (EB) a value is: ", EB$opt$a, "\n",
             ",AIC value is: ", EB$opt$aic, "\n",
             ",and lnL is: ", EB$opt$lnL, "\n"))
  
  # Function to compare models and print results
  compare_models <- function(model1, model2, model1_name, model2_name) {
    if (model1$opt$aic < model2$opt$aic) {
      larger_model <- model1_name
      smaller_model <- model2_name
      larger_lnL <- model1$opt$lnL
      smaller_lnL <- model2$opt$lnL
      larger_aic <- model1$opt$aic
      smaller_aic <- model2$opt$aic
    } else {
      larger_model <- model2_name
      smaller_model <- model1_name
      larger_lnL <- model2$opt$lnL
      smaller_lnL <- model1$opt$lnL
      larger_aic <- model2$opt$aic
      smaller_aic <- model1$opt$aic
    }
    
    LRT_stat <- -2 * (smaller_lnL - larger_lnL)
    p_value <- pchisq(LRT_stat, df = 1, lower.tail = FALSE)
    
    cat(paste0("Comparison: ", model1_name, " vs ", model2_name, "\n"))
    cat(paste0("Likelihood ratio test statistic: ", LRT_stat, "\n"))
    cat(paste0("p-value: ", p_value, "\n"))
    cat(paste0("smaller aic: ", smaller_model, 
               " (", smaller_aic, ") vs. ",larger_model, 
               " (", larger_aic, ")\n\n"))
  }
  
  # Print the comparisons
  cat(paste0("\n","...", "\n"))
  cat(paste0("Comparing models", "\n"))
  cat(paste0("...", "\n", "\n"))
  
  # BM vs OU
  compare_models(BM, OU, "BM", "OU")
  
  # BM vs EB
  compare_models(BM, EB, "BM", "EB")
  
  # OU vs EB
  compare_models(OU, EB, "OU", "EB")
}

#Calculate the modes of evolution of all the traits
rates_complete_evol<- modes_evol(MossRC_label)
#For rates complete the best model is OU with a value of 2.71828182845905
rates_separate_evol<- modes_evol(MossSC_label)
#For rates separate the best model is OU with a value of 2.71828182845905
immediate_complete_evol<- modes_evol(MossImC_label)
#For immediate loss complete the best model is OU with a value of 2.71828182845905
immediate_separate_evol<- modes_evol(MossImS_label)
#For immediate loss separate the best model is OU with a value of 2.71828182845905


#### Set the tree to an OU mode of evolution with value 2.71 ####

moss_tree_OU<-rescale(Moss_tree, "OU", alpha=2.71828182845905)

################################################################################
#### MCMCglmm z - Rate complete ~ substrate ####

Moss_ID_data$z_Rate_complete <- scale(Moss_ID_data$Rate_complete)

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

p1 = list(B=list(mu=rep(0,4), V=diag(4)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m1a<-MCMCglmm(z_Rate_complete ~ Substrate-1, random = ~Genus,
             ginverse=list(Genus=inv.mossphylo$Ainv),
             family ="gaussian", data = Moss_ID_data,
             prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
m1b<-MCMCglmm(z_Rate_complete ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
m1c<-MCMCglmm(z_Rate_complete ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))
#And random effects
plot(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))

gelman.diag(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))

summary(m1a)

# boxplot(Moss_ID_data$z_Rate_complete~Moss_ID_data$Substrate)
# boxplot(Moss_ID_data$z_Rate_separate~Moss_ID_data$Substrate)
#### MCMCglmm z - Rate separate ~ substrate  ####

Moss_ID_data$z_Rate_separate <- scale(Moss_ID_data$Rate_separate)

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

p2 = list(B=list(mu=rep(0,4), V=diag(4)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m2a<-MCMCglmm(z_Rate_separate ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p2, nitt=110000, burnin=10000, thin=100,verbose=F)
m2b<-MCMCglmm(z_Rate_separate ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p2, nitt=110000, burnin=10000, thin=100,verbose=F)
m2c<-MCMCglmm(z_Rate_separate ~ Substrate-1, random = ~Genus,
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

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

p3 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m3a<-MCMCglmm(z_Rate_complete ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p3, nitt=110000, burnin=10000, thin=100,verbose=F)
m3b<-MCMCglmm(z_Rate_complete ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p3, nitt=110000, burnin=10000, thin=100,verbose=F)
m3c<-MCMCglmm(z_Rate_complete ~ Environment-1, random = ~Genus,
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

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

p4 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

m4a<-MCMCglmm(z_Rate_separate ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p4, nitt=110000, burnin=10000, thin=100,verbose=F)
m4b<-MCMCglmm(z_Rate_separate ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p4, nitt=110000, burnin=10000, thin=100,verbose=F)
m4c<-MCMCglmm(z_Rate_separate ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=p4, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(m4a$Sol,m4b$Sol,m4c$Sol))
#And random effects
plot(mcmc.list(m4a$Sol,m4b$Sol,m4c$Sol))

gelman.diag(mcmc.list(m4a$Sol,m4b$Sol,m4c$Sol))

summary(m4a)



################################################################################
#### MCMCglmm Immediate loss complete ~ substrate ####

#Moss_ID_data$z_Rate_complete <- scale(Moss_ID_data$Immediate_diff_complete)

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

i_p1 = list(B=list(mu=rep(0,4), V=diag(4)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

i_m1a<-MCMCglmm(Immediate_diff_complete ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p1, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m1b<-MCMCglmm(Immediate_diff_complete ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p1, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m1c<-MCMCglmm(Immediate_diff_complete ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p1, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(i_m1a$Sol,i_m1b$Sol,i_m1c$Sol))
#And random effects
plot(mcmc.list(i_m1a$Sol,i_m1b$Sol,i_m1c$Sol))

gelman.diag(mcmc.list(i_m1a$Sol,i_m1b$Sol,i_m1c$Sol))

summary(i_m1a)

# boxplot(Moss_ID_data$z_Rate_complete~Moss_ID_data$Substrate)
# boxplot(Moss_ID_data$z_Rate_separate~Moss_ID_data$Substrate)
#### MCMCglmm Immediate loss separate ~ substrate  ####

#Moss_ID_data$z_Rate_separate <- scale(Moss_ID_data$Immediate_diff_separated)

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

i_p2 = list(B=list(mu=rep(0,4), V=diag(4)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

i_m2a<-MCMCglmm(Immediate_diff_separated ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p2, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m2b<-MCMCglmm(Immediate_diff_separated ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p2, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m2c<-MCMCglmm(Immediate_diff_separated ~ Substrate-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p2, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(i_m2a$Sol,i_m2b$Sol,i_m2c$Sol))
#And random effects
plot(mcmc.list(i_m2a$Sol,i_m2b$Sol,i_m2c$Sol))

gelman.diag(mcmc.list(i_m2a$Sol,i_m2b$Sol,i_m2c$Sol))

summary(i_m2a)


#### MCMCglmm Immediate loss complete ~ Environment ####

#Moss_ID_data$z_Rate_complete <- scale(Moss_ID_data$Rate_complete)

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

i_p3 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

i_m3a<-MCMCglmm(Immediate_diff_complete ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p3, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m3b<-MCMCglmm(Immediate_diff_complete ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p3, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m3c<-MCMCglmm(Immediate_diff_complete ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p3, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(i_m3a$Sol,i_m3b$Sol,i_m3c$Sol))
#And random effects
plot(mcmc.list(i_m3a$Sol,i_m3b$Sol,i_m3c$Sol))

gelman.diag(mcmc.list(i_m3a$Sol,i_m3b$Sol,i_m3c$Sol))

summary(i_m3a)

#### MCMCglmm Immediate loss separate ~ Environment  ####

#Moss_ID_data$z_Rate_separate <- scale(Moss_ID_data$Rate_separate)

#Inverse tree using OU
inv.mossphylo<-inverseA(moss_tree_OU,nodes="TIPS",scale=TRUE)

i_p4 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))

i_m4a<-MCMCglmm(Immediate_diff_separated ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p4, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m4b<-MCMCglmm(Immediate_diff_separated ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p4, nitt=110000, burnin=10000, thin=100,verbose=F)
i_m4c<-MCMCglmm(Immediate_diff_separated ~ Environment-1, random = ~Genus,
              ginverse=list(Genus=inv.mossphylo$Ainv),
              family ="gaussian", data = Moss_ID_data,
              prior=i_p4, nitt=110000, burnin=10000, thin=100,verbose=F)

#Plot the fixed effects
plot(mcmc.list(i_m4a$Sol,i_m4b$Sol,i_m4c$Sol))
#And random effects
plot(mcmc.list(i_m4a$Sol,i_m4b$Sol,i_m4c$Sol))

gelman.diag(mcmc.list(i_m4a$Sol,i_m4b$Sol,i_m4c$Sol))

summary(i_m4a)





################################################################################
#### Pairwise comparison immediate loss per Substrate & Environment ####
# Function to compute posterior modes and HPD intervals for specified comparisons
pairwise_comparison <- function(model, comparisons) {
  results <- data.frame(Comparison = character(),
                        PosteriorMode = numeric(),
                        LowerBoundary = numeric(),
                        UpperBoundary = numeric(),
                        stringsAsFactors = FALSE)
  
  for (comparison in comparisons) {
    var1 <- comparison[1]
    var2 <- comparison[2]
    diff <- model$Sol[, var1] - model$Sol[, var2]
    mode <- posterior.mode(diff)
    hpd <- HPDinterval(diff)
    
    results <- rbind(results, data.frame(Comparison = paste0(colnames(model$Sol)[var1], " vs ", colnames(model$Sol)[var2]),
                                         PosteriorMode = mode,
                                         LowerBoundary = hpd[1],
                                         UpperBoundary = hpd[2]))
  }
  
  return(results)
}
comparisons_substrate <- list(c(1, 2),  # Rock vs Sand
                    c(2, 3),  # Sand vs Soil
                    c(3, 4),  # Soil vs Wood
                    c(1, 3),  # Rock vs Soil
                    c(1, 4),  # Rock vs Wood
                    c(2, 4))  # Sand vs Wood

pairwise_table_m1a <- pairwise_comparison(m1a, comparisons_substrate)
pairwise_table_m2a <- pairwise_comparison(m2a, comparisons_substrate)
pairwise_table_m3a <- pairwise_comparison(m3a, comparisons = list(c(1,2)))
pairwise_table_m4a <- pairwise_comparison(m4a, comparisons = list(c(1,2)))

pairwise_table_i_m1a <- pairwise_comparison(i_m1a, comparisons_substrate)
pairwise_table_i_m2a <- pairwise_comparison(i_m2a, comparisons_substrate)
pairwise_table_i_m3a <- pairwise_comparison(i_m3a, comparisons = list(c(1,2)))
pairwise_table_i_m4a <- pairwise_comparison(i_m4a, comparisons = list(c(1,2)))
################################################################################
#### MCMCglmm Ancestral reconstruction rate_complete ####
#Inverse tree using OU
inv.mossphylo_all<-inverseA(moss_tree_OU,nodes="ALL",scale=TRUE)

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
posterior.mode(m5b$Sol)
posterior.mode(m5c$Sol)

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
  moss_tree_OU$tip.label,Moss_Genus_data$Genus)]

#plot them and add the node numbers to being able to do 4.10
plot(moss_tree_OU, cex=0.8, no.margin =T, label.offset = 0.02)
nodelabels(pch=21, cex=(blupsm5a$estimate[1:moss_tree_OU$Nnode]/2.5))
tiplabels(pch=21,cex=mean_Rate_complete_phylo/2.5,bg="black")


#### check phylogenetic signal heredability rate_complete ####

m5aPhyloSig<-m5a$VCV[,'Genus']/((m5a$VCV[,'Genus']+m5a$VCV[,'units']))

#posterior.mode is the heredability value in this case
posterior.mode(m5aPhyloSig)
HPDinterval(m5aPhyloSig)

plot(m5aPhyloSig)

#### MCMCglmm Ancestral reconstruction rate_separate ####
#Inverse tree using OU
inv.mossphylo_all<-inverseA(moss_tree_OU,nodes="ALL",scale=TRUE)

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
  moss_tree_OU$tip.label,Moss_Genus_data$Genus)]

#plot them and add the node numbers to being able to do 4.10
plot(moss_tree_OU, cex=0.8, no.margin =T, label.offset = 0.02)
nodelabels(pch=21, cex=(blupsm6a$estimate[1:moss_tree_OU$Nnode]/2.5))
tiplabels(pch=21,cex=mean_Rate_separate_phylo/2.5,bg="black")

#### check phylogenetic signal heredability rate_separate ####

m6aPhyloSig<-m6a$VCV[,'Genus']/((m6a$VCV[,'Genus']+m6a$VCV[,'units']))

#posterior.mode is the heredability value in this case
posterior.mode(m6aPhyloSig)
HPDinterval(m6aPhyloSig)

plot(m6aPhyloSig)

################################################################################
#### Graphics ####
wes_palette <- paletteer_d("tvthemes::Day")
selected_palette <- paletteer_d("tvthemes::gravityFalls")
pale_wes_palette <- scales::alpha(wes_palette, 0.5)


ggplot(data = Moss_ID_data) +
  geom_violin(aes(x = Substrate, y = Immediate_diff_separated, fill = Substrate), 
              trim = FALSE) +
  stat_summary(fun = mean, geom = "point", aes(
    x = Substrate, y = Immediate_diff_separated, group = Genus, col = Genus)) +
  stat_summary(fun = mean, geom = "line", aes(
    x = Substrate, y = Immediate_diff_separated, group = Genus, col = Genus)) +
  scale_fill_manual(values = pale_wes_palette) +
  scale_color_manual(values = selected_palette) +
  labs(title = "Violin Plot of Immediate Loss Separted by Substrate",
       x = "Substrate",
       y = "Immediate Loss Complete") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = Moss_ID_data) +
  geom_violin(aes(x = Substrate, y = Immediate_diff_complete, fill = Substrate), 
              trim = FALSE) +
  stat_summary(fun = mean, geom = "point", aes(
    x = Substrate, y = Immediate_diff_complete, group = Genus, col = Genus)) +
  stat_summary(fun = mean, geom = "line", aes(
    x = Substrate, y = Immediate_diff_complete, group = Genus, col = Genus)) +
  scale_fill_manual(values = pale_wes_palette) +
  scale_color_manual(values = selected_palette) +
  labs(title = "Violin Plot of Immediate Loss Separted by Substrate",
       x = "Substrate",
       y = "Immediate Loss Complete") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#### PGLS Ancestral reconstruction substrate ####
# inv.mossphylo_all<-inverseA(Moss_tree,nodes="ALL",scale=TRUE)
# # set priors
# p7 = list(B=list(mu=rep(0,1), V=diag(1)*1e+8), G=list(G1=list(V=1,nu=0.002)),
#           R=list(V=1,nu=0.002))
# 
# #model
# m7<-MCMCglmm(Substrate ~ 1, random = ~Genus,
#              ginverse=list(Genus=inv.mossphylo_all$Ainv),
#              family ="categorical", data = Moss_ID_data, pr=TRUE, 
#              rcov = ~us(trait):units,
#              prior=p7, nitt=110000, burnin=10000, thin=100,verbose=F)
# 
# #This are the logit not transformed back of the values for each family 
# #and each node
# posterior.mode(m7$Sol)
# 
# #BLUPS
# #To get estimates of ancestral states add intercept to node BLUPs
# blupsm7<-data.frame(effect=colnames(m7$Sol), estimate=posterior.mode(
#   m7$Sol[,'(Intercept)'])+posterior.mode(m7$Sol),CI=HPDinterval(
#     m7$Sol[,'(Intercept)']+m7$Sol))
# 
# #We added the intercept value to all values and so "Node1" value is incorrect 
# #- lets replace it
# blupsm7['(Intercept)','estimate']<-posterior.mode(m7$Sol[,'(Intercept)'])
# blupsm7['(Intercept)',c('CI.lower','CI.upper')]<-HPDinterval(
#   m7$Sol[,'(Intercept)'])
# 
# #The intercept is the value of the root ("Node1") so lets rename it
# blupsm7$effect<-as.character(blupsm7$effect)
# blupsm7['(Intercept)','effect']<-"Node1"
# 
# #and remove the text "treetip." from the effect
# blupsm7$effect<-gsub("Genus.","",blupsm7$effect)
# 
# #match with phylogeny
# #match the vitamin data and the insect tree data
# Rate_complete_phylo<-(Moss_ID_data$Rate_separate)[match(
#   Moss_tree$tip.label,Moss_ID_data$Genus)]
# 
# #plot them and add the node numbers to being able to do 4.10
# plot(Moss_tree, cex=0.8, no.margin =T, label.offset = 0.7)
# nodelabels(pch=21, cex=(abs(blupsm7$estimate[1:Moss_tree$Nnode])*200))
# tiplabels(pch=21,cex=Rate_complete_phylo*200,bg="black")

#### Calculate modes of evolution Genus ####

# MossRC_label<-as.numeric(Moss_rates_Genus$Rate_complete[match(
#   Moss_tree$tip.label,Moss_rates_Genus$Genus)])
# names(MossRC_label)<-Moss_tree$tip.label
# MossSC_label<-as.numeric(Moss_rates_Genus$Rate_separate[match(
#   Moss_tree$tip.label,Moss_rates_Genus$Genus)])
# names(MossSC_label)<-Moss_tree$tip.label
# 
# BM_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "BM")
# OU_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "OU", 
#                     bounds=list(alpha=c(0,10000)))
# EB_C<-fitContinuous(phy= Moss_tree, dat = MossRC_label, model = "EB")
# 
# BM_C$opt$aic
# OU_C$opt$aic
# EB_C$opt$aic
# 
# ##BM versus OU
# BMvsOU_C <- -2*(BM_C$opt$lnL - OU_C$opt$lnL)
# BMvsOU_C
# pchisq(BMvsOU_C, df=1,lower.tail = FALSE)
# 
# ##BM versus EB
# BMvsEB_C <- -2*(EB_C$opt$lnL - BM_C$opt$lnL)
# BMvsEB_C
# pchisq(BMvsEB_C, df=1,lower.tail = FALSE)
# 
# ##OU versus EB
# OUvsEB_C <- -2*(EB_C$opt$lnL - OU_C$opt$lnL)
# OUvsEB_C 
# pchisq(OUvsEB_C, df=1,lower.tail = FALSE)
# 
# #OU is better than both EB and BM. The values are restrained to ne value in the model
# #Use alpha when doing the phylogenetic approaches: 1387.735201
#### graph ####
# ggplot(data = Moss_ID_data) +
#   geom_boxplot(aes(x = Substrate, y = Immediate_diff_complete, fill = Substrate)) +
#   stat_summary(fun.y="mean", geom = "point", aes(x =Substrate, y = Immediate_diff_complete, group = Genus, col=Genus)) +
#   stat_summary(fun.y="mean", geom = "line", aes(x =Substrate, y = Immediate_diff_complete, group = Genus, col=Genus))

#### More unused code ####
# #posterior.mode(i_m1a$Sol[,"Substraterock"]-i_m1a$Sol[,"Substratesand"])
# #HPDinterval(i_m1a$Sol[,"Substraterock"]-i_m1a$Sol[,"Substratesand"])
# colnames(i_m1a$Sol)
#
# #Rock vs sand If the HPD interval is not
# posterior.mode(i_m1a$Sol[,1]-i_m1a$Sol[,2])
# HPDinterval(i_m1a$Sol[,1]-i_m1a$Sol[,2])
# 
# #Sand vs soil If the HPD interval is not
# posterior.mode(i_m1a$Sol[,2]-i_m1a$Sol[,3])
# HPDinterval(i_m1a$Sol[,2]-i_m1a$Sol[,3])
# 
# #Soil vs wood If the HPD interval is not
# posterior.mode(i_m1a$Sol[,3]-i_m1a$Sol[,4])
# HPDinterval(i_m1a$Sol[,3]-i_m1a$Sol[,4])
# 
# #Rock vs soil If the HPD interval is not
# posterior.mode(i_m1a$Sol[,1]-i_m1a$Sol[,3])
# HPDinterval(i_m1a$Sol[,1]-i_m1a$Sol[,3])
# 
# #Rock vs wood If the HPD interval is not
# posterior.mode(i_m1a$Sol[,1]-i_m1a$Sol[,4])
# HPDinterval(i_m1a$Sol[,1]-i_m1a$Sol[,4])
# 
# #Sand vs wood If the HPD interval is not
# posterior.mode(i_m1a$Sol[,2]-i_m1a$Sol[,4])
# HPDinterval(i_m1a$Sol[,2]-i_m1a$Sol[,4])