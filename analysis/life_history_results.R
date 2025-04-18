setwd("/Users/bell/Documents/GitHub/LTREB-life-history/analysis")
setwd("C:/Users/tm9/Dropbox/github/LTREB-life-history")
library(tidyverse)
library(vegan)
#library(factoextra)
##read in life history outputs
lifehistorypost<-read.csv("analysis/lifehistorypost.csv")

# PCA ---------------------------------------------------------------------
## need to pivot trait data from wide to long
## select traits that will go into PCA
lifehistorypost %>% 
  select(species,R0_em,R0_ep) %>% 
  pivot_longer(R0_em:R0_ep,names_to="endo",values_to="R0") %>% 
  mutate(endo=case_when(endo=="R0_em"~"S-",endo=="R0_ep"~"S+"))->R0
lifehistorypost %>% 
  select(species,G_em,G_ep) %>% 
  pivot_longer(G_em:G_ep,names_to="endo",values_to="G") %>% 
  mutate(endo=case_when(endo=="G_em"~"S-",endo=="G_ep"~"S+"))->G
lifehistorypost %>% 
  select(species,meanelexp_em,meanelexp_ep) %>% 
  pivot_longer(meanelexp_em:meanelexp_ep,names_to="endo",values_to="meanelexp") %>% 
  mutate(endo=case_when(endo=="meanelexp_em"~"S-",endo=="meanelexp_ep"~"S+"))->meanelexp
lifehistorypost %>% 
  select(species,matlifexp_em,matlifexp_ep) %>% 
  pivot_longer(matlifexp_em:matlifexp_ep,names_to="endo",values_to="matlifexp") %>% 
  mutate(endo=case_when(endo=="matlifexp_em"~"S-",endo=="matlifexp_ep"~"S+"))->matlifexp
lifehistorypost %>% 
  select(species,longevity_em,longevity_ep) %>% 
  pivot_longer(longevity_em:longevity_ep,names_to="endo",values_to="longevity") %>% 
  mutate(endo=case_when(endo=="longevity_em"~"S-",endo=="longevity_ep"~"S+"))->longevity
lifehistorypost %>% 
  select(species,entropyd_em,entropyd_ep) %>% 
  pivot_longer(entropyd_em:entropyd_ep,names_to="endo",values_to="entropyd") %>% 
  mutate(endo=case_when(endo=="entropyd_em"~"S-",endo=="entropyd_ep"~"S+"))->entropyd
lifehistorypost %>% 
  select(species,firstrepro_em,firstrepro_ep) %>% 
  pivot_longer(firstrepro_em:firstrepro_ep,names_to="endo",values_to="firstrepro") %>% 
  mutate(endo=case_when(endo=="firstrepro_em"~"S-",endo=="firstrepro_ep"~"S+"))->firstrepro

pca.dat<-bind_cols(R0,G$G,meanelexp$meanelexp,longevity$longevity,entropyd$entropyd, firstrepro$firstrepro)
names(pca.dat)<-c("Species","Endo","R0","Gen_time","Life_expect","Longevity","EntropyD","FirstRepro")
#row.names(pca.dat)<-c("AGPE-","AGPE+",
#                      "ELRI-","ELRI+",
#                      "ELVI-","ELVI+",
#                      "FESU-","FESU+",
#                      "POAL-","POAL+",
#                      "POAU-","POAU+",
#                      "POSY-","POSY+")

## first create a reference PCA from posterior means
pca.dat %>% 
  group_by(Species,Endo) %>% 
  summarise_all(mean) -> mean.pca.dat

## reference PCA
lifehistory_pca_mean<-prcomp(mean.pca.dat[,-(1:2)],scale=T)
lifehistory_pca_mean$rotation

## now align posterior samples with respect to reference

# From chatgpt:
# Align PCA loadings from a sample to a reference using Procrustes rotation
align_pca <- function(sample_loadings, ref_loadings) {
  proc <- procrustes(ref_loadings, sample_loadings, symmetric = TRUE)
  aligned <- proc$Yrot
  return(aligned)
}

##first posterior draw

plot(lifehistory_pca$x,pch=c(1,16))
arrows(lifehistory_pca$x[c(1,3,5,7,9,11,13),1],
       lifehistory_pca$x[c(1,3,5,7,9,11,13),2],
       lifehistory_pca$x[c(2,4,6,8,10,12,14),1],
       lifehistory_pca$x[c(2,4,6,8,10,12,14),2])

fviz_pca_var(lifehistory_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

groups<-as.factor(R0$endo)
fviz_pca_biplot(lifehistory_pca, repel = TRUE,
                col.var = "gray", # Variables color
                col.ind = groups, # color by groups
                palette = c("tomato",  "cornflowerblue")  # Individuals color
)+
  annotate("segment",x=lifehistory_pca$x[c(1,3,5,7,9,11,13),1], 
                   y=lifehistory_pca$x[c(1,3,5,7,9,11,13),2], 
                   xend=lifehistory_pca$x[c(2,4,6,8,10,12,14),1], 
                   yend=lifehistory_pca$x[c(2,4,6,8,10,12,14),2],
               arrow = arrow(length=unit(.4, 'cm')),col=alpha("cornflowerblue",0.25))+
  theme(legend.position = "none")

# posterior density plots of single traits ---------------------------------
## generate E+/E- contrasts
lifehistorypost$R0_diff<-lifehistorypost$R0_ep-lifehistorypost$R0_em
lifehistorypost$G_diff<-lifehistorypost$G_ep-lifehistorypost$G_em
lifehistorypost$lambda_diff<-lifehistorypost$lambda_ep-lifehistorypost$lambda_em
lifehistorypost$pRep_diff<-lifehistorypost$pRep_ep-lifehistorypost$pRep_em
lifehistorypost$La_diff<-lifehistorypost$La_ep-lifehistorypost$La_em
lifehistorypost$matlifexp_diff<-lifehistorypost$matlifexp_ep-lifehistorypost$matlifexp_em
lifehistorypost$meanelexp_diff<-lifehistorypost$meanelexp_ep-lifehistorypost$meanelexp_em
lifehistorypost$entropyd_diff<-lifehistorypost$entropyd_ep-lifehistorypost$entropyd_em
lifehistorypost$entropyk_diff<-lifehistorypost$entropyk_ep-lifehistorypost$entropyk_em
lifehistorypost$gini_diff<-lifehistorypost$gini_ep-lifehistorypost$gini_em
lifehistorypost$longevity_diff<-lifehistorypost$longevity_ep-lifehistorypost$longevity_em
lifehistorypost$firstrepro_diff<-lifehistorypost$firstrepro_ep-lifehistorypost$firstrepro_em

## some quick and dirty plots of life history effects
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(pr_epos_lambda = mean(lambda_diff>0),
            pr_epos_R0 = mean(R0_diff>0),
            pr_epos_G = mean(G_diff>0),
            pr_epos_meanelexp = mean(meanelexp_diff>0),
            pr_epos_longevity = mean(longevity_diff>0))

##lambda
ggplot(lifehistorypost,aes(lambda_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##R0
ggplot(lifehistorypost,aes(R0_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")


##G
ggplot(lifehistorypost,aes(G_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##pRep-- there is no variation here - not sure if this is calculated correctly
ggplot(lifehistorypost,aes(pRep_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##La--same
ggplot(lifehistorypost,aes(La_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##mean life expectancy
ggplot(lifehistorypost,aes(meanelexp_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##life expectancy from maturity
ggplot(lifehistorypost,aes(matlifexp_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##entropy D
ggplot(lifehistorypost,aes(entropyd_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

##entropy K
ggplot(lifehistorypost,aes(entropyk_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

## Gini -- no variation here
ggplot(lifehistorypost,aes(gini_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

## longevity
ggplot(lifehistorypost,aes(longevity_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

## longevity
ggplot(lifehistorypost,aes(firstrepro_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")

                