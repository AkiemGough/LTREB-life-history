setwd("/Users/bell/Documents/GitHub/LTREB-life-history/analysis")
setwd("C:/Users/tm9/Dropbox/github/LTREB-life-history")
library(tidyverse)
library(vegan)
library(factoextra)
library(ggrepel)
library(ggridges)
library(patchwork)
##read in life history outputs
lifehistorypost<-read.csv("analysis/lifehistorypost.csv")
## add posterior draw -- we used 100 samples for each species
lifehistorypost$draw<-rep(1:n_post,times=7) #n_post from vital rates script

# PCA ---------------------------------------------------------------------
## need to pivot trait data from wide to long
## select traits that will go into PCA
lifehistorypost %>% 
  select(draw,species,R0_em,R0_ep) %>% 
  pivot_longer(R0_em:R0_ep,names_to="endo",values_to="R0") %>% 
  mutate(endo=case_when(endo=="R0_em"~"S-",endo=="R0_ep"~"S+"))->R0
lifehistorypost %>% 
  select(draw,species,G_em,G_ep) %>% 
  pivot_longer(G_em:G_ep,names_to="endo",values_to="G") %>% 
  mutate(endo=case_when(endo=="G_em"~"S-",endo=="G_ep"~"S+"))->G
lifehistorypost %>% 
  select(draw,species,meanelexp_em,meanelexp_ep) %>% 
  pivot_longer(meanelexp_em:meanelexp_ep,names_to="endo",values_to="meanelexp") %>% 
  mutate(endo=case_when(endo=="meanelexp_em"~"S-",endo=="meanelexp_ep"~"S+"))->meanelexp
lifehistorypost %>% 
  select(draw,species,matlifexp_em,matlifexp_ep) %>% 
  pivot_longer(matlifexp_em:matlifexp_ep,names_to="endo",values_to="matlifexp") %>% 
  mutate(endo=case_when(endo=="matlifexp_em"~"S-",endo=="matlifexp_ep"~"S+"))->matlifexp
lifehistorypost %>% 
  select(draw,species,longevity_em,longevity_ep) %>% 
  pivot_longer(longevity_em:longevity_ep,names_to="endo",values_to="longevity") %>% 
  mutate(endo=case_when(endo=="longevity_em"~"S-",endo=="longevity_ep"~"S+"))->longevity
lifehistorypost %>% 
  select(draw,species,entropyd_em,entropyd_ep) %>% 
  pivot_longer(entropyd_em:entropyd_ep,names_to="endo",values_to="entropyd") %>% 
  mutate(endo=case_when(endo=="entropyd_em"~"S-",endo=="entropyd_ep"~"S+"))->entropyd
lifehistorypost %>% 
  select(draw,species,firstrepro_em,firstrepro_ep) %>% 
  pivot_longer(firstrepro_em:firstrepro_ep,names_to="endo",values_to="firstrepro") %>% 
  mutate(endo=case_when(endo=="firstrepro_em"~"S-",endo=="firstrepro_ep"~"S+"))->firstrepro

pca.dat<-bind_cols(R0,G$G,meanelexp$meanelexp,longevity$longevity,entropyd$entropyd, firstrepro$firstrepro)
names(pca.dat)<-c("Draw","Species","Endo","R0","Gen_time","Life_expect","Longevity","EntropyD","FirstRepro")

pca.dat %>% select(-Draw) %>% 
  group_by(Species,Endo) %>% 
  summarise_all(mean) -> mean.pca.dat

## reference PCA and loadings based on posterior mean
pca_ref <- prcomp(mean.pca.dat[,-(1:2)], scale. = TRUE)
ref_loadings <- pca_ref$rotation  # this is your fixed basis
ref_center <- pca_ref$center
ref_scale <- pca_ref$scale

## align each posterior draw to the reference loadings
## calculate shift along PC1 and PC2
posterior_shifts <- list()
for (d in unique(pca.dat$Draw)) {
  draw_df <- pca.dat %>% filter(Draw == d)
  # scale using reference
  mat <- draw_df %>% select(-c(Draw,Species,Endo)) %>% 
    scale(center = ref_center, scale = ref_scale)
  # project onto fixed PCA basis
  scores <- mat %*% ref_loadings[, 1:2]  # project to PC1 and PC2
  scores_df <- bind_cols(draw_df %>% select(Species, Endo), as.data.frame(scores))
  # compute shift vector for each species
  shift_df <- scores_df %>%
    group_by(Species) %>%
    summarise(
      dPC1 = PC1[Endo == "S+"] - PC1[Endo == "S-"],
      dPC2 = PC2[Endo == "S+"] - PC2[Endo == "S-"]
    ) %>%
    mutate(Draw = d)
  posterior_shifts[[length(posterior_shifts)+1]] <- shift_df
}
posterior_shift_df <- bind_rows(posterior_shifts)

posterior_shift_means<-posterior_shift_df %>% 
  group_by(Species) %>% 
  summarise(mean1 = mean(dPC1),
            mean2 = mean(dPC2))%>%
    mutate(direction1 = ifelse(mean1 >= 0, "positive", "negative"),
           direction2 = ifelse(mean2 >= 0, "positive", "negative"))

##Figure of shift in PC space
scores_df <- as_tibble(pca_ref$x[, 1:2], .name_repair = "unique") %>%
  mutate(Species = mean.pca.dat$Species,
         Endo = mean.pca.dat$Endo)
arrows_df <- scores_df %>%
  pivot_wider(names_from = Endo, values_from = c(PC1, PC2)) %>%
  mutate(
    x_start = `PC1_S-`,
    y_start = `PC2_S-`,
    x_end   = `PC1_S+`,
    y_end   = `PC2_S+`)
loadings <- as_tibble(pca_ref$rotation[, 1:2], rownames = "Trait")

ggplot() +
  geom_point(data = scores_df, aes(x = PC1, y = PC2, color = Endo), size = 4) +
  geom_segment(data = arrows_df,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.75,
               color = "black") +
  scale_color_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue")) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.6,
               color = "gray80") +
  geom_text_repel(data = loadings,
                  aes(x = PC1, y = PC2, label = Trait),
                  color = "gray80", size = 3) +
  geom_text_repel(
    data = arrows_df,
    aes(x = (x_start + x_end)/2, y = (y_start + y_end)/2, label = Species),
    size = 4,
    color = "black",
    max.overlaps = Inf
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.9),   # (x, y) in normalized plot coordinates
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )+
  labs(title = "A)",
       x = "PC1", y = "PC2", color = "Symbiont status")->PCA_plot

## Figure of posterior shifts
lims1 <- quantile(posterior_shift_df$dPC1, probs = c(0.01, 0.99))
ggplot(posterior_shift_df, aes(x = dPC1, y = Species, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "tomato",
    mid = "white",
    high = "cornflowerblue",
    midpoint = 0,
    limits = c(-5, 5),  # clip values outside this range
    oob = scales::squish,  # squish outliers into endpoint colors
    name = "Shift"
  ) +  coord_cartesian(xlim = lims1) + 
  theme_minimal(base_size = 13) + theme(legend.position = "none") +
  labs(title = "B)",
       x = expression(Delta~PC1), y = "Species")+ geom_vline(xintercept = 0)+
  geom_point(
    data = posterior_shift_means,
    aes(x = mean1, y = Species, color = direction1),  # now fill is mapped!
    inherit.aes = FALSE,
    shape = 16,
    size = 2
  ) +
  scale_color_manual(
    values = c("positive" = "cornflowerblue", "negative" = "tomato")
  )-> PC1_plot

lims2 <- quantile(posterior_shift_df$dPC2, probs = c(0.01, 0.99))
ggplot(posterior_shift_df, aes(x = dPC2, y = Species, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "tomato",
    mid = "white",
    high = "cornflowerblue",
    midpoint = 0,
    limits = c(-5, 5),  # clip values outside this range
    oob = scales::squish,  # squish outliers into endpoint colors
    name = "Shift"
  ) + coord_cartesian(xlim = lims1) + 
  theme_minimal(base_size = 13) + theme(legend.position = "none") +
  labs(title = "C)",
       x = expression(Delta~PC2), y = "Species")+ geom_vline(xintercept = 0)+
  geom_point(
    data = posterior_shift_means,
    aes(x = mean2, y = Species, color = direction2),  # now fill is mapped!
    inherit.aes = FALSE,
    shape = 16,
    size = 2
  ) +
  scale_color_manual(
    values = c("positive" = "cornflowerblue", "negative" = "tomato")
  )-> PC2_plot

PCA_combo <- PCA_plot / (PC1_plot | PC2_plot)

ggsave("manuscript/figures/pca_shift.jpg",
       plot = PCA_combo,
       width = 6,      # width in inches
       height = 8,      # height in inches
       dpi = 300,       # resolution
       units = "in")

ggplot(posterior_shift_df, aes(x = dPC1, y = Species, fill = factor(stat(quantile)))) +
stat_density_ridges(
  geom = "density_ridges_gradient",
  calc_ecdf = TRUE,
  quantiles = c(0.025, 0.975)
) +
  scale_fill_manual(
    name = "Probability", values = c("white", "gray", "white"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  )+
  theme_minimal(base_size = 13) +
  labs(
    title = "Posterior Distribution of PC1 Shifts (S+ − S−)",
    x = "Shift in PC1", y = "Species"
  ) + geom_vline(xintercept = 0)

## mean shift and probability of PC shift
posterior_shift_df %>% 
  group_by(Species) %>% 
  summarise(mean_shift1=mean(dPC1),
            mean_shift2=mean(dPC2),
            prob_shift1=mean(dPC1>0),
            prob_shift2=mean(dPC2>0))

groups<-as.factor(mean.pca.dat$Endo)
fviz_pca_biplot(pca_ref, repel = TRUE,
                col.var = "gray", # Variables color
                col.ind = groups, # color by groups
                palette = c("tomato",  "cornflowerblue")  # Individuals color
)+
  annotate("segment",x=mean.pca.dat$x[c(1,3,5,7,9,11,13),1], 
           y=mean.pca.dat$x[c(1,3,5,7,9,11,13),2], 
           xend=mean.pca.dat$x[c(2,4,6,8,10,12,14),1], 
           yend=mean.pca.dat$x[c(2,4,6,8,10,12,14),2],
           arrow = arrow(length=unit(.4, 'cm')),col=alpha("cornflowerblue",0.25))+
  theme(legend.position = "none")




lifehistory_pca<-prcomp(pca.dat[pca.dat$Draw == 5, -(1:3)], center = FALSE, scale. = FALSE)
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

## data frame for life history effects
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(`R0` = mean(R0_diff>0),
            `Generation time` = mean(G_diff>0),
            `Life expectancy` = mean(meanelexp_diff>0),
            `Longevity` = mean(longevity_diff>0),
            `Entropy` = mean(entropyd_diff>0),
            `Reproductive age` = mean(firstrepro_diff>0))->life_history_effects

life_history_effects %>%
  pivot_longer(cols = -species, names_to = "metric", values_to = "probpos") %>% 
  ggplot(aes(x = species, y = metric, fill = probpos)) +
  geom_tile() +
  scale_fill_gradientn(colors=c("tomato", "grey", "cornflowerblue")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "S+ Adv.e", x = "Species", y = "Metric", title = "Heatmap of Life History Metrics")

life_history_effects %>%
  pivot_longer(cols = -species, names_to = "metric", values_to = "probpos") %>% 
  ggplot(aes(x = species, y = metric, fill = probpos)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", probpos)), size = 3) +  # Add values to each tile
  scale_fill_gradient2(
    low = "tomato",
    mid = "grey",
    high = "cornflowerblue",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "Probability of\nS+ advantage"
  ) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Life history trait") -> LHtraits_heatmap

ggsave("manuscript/figures/LHtraits_heatmap.jpg",
       plot = LHtraits_heatmap,
       width = 10,      # width in inches
       height = 8,      # height in inches
       dpi = 300,       # resolution
       units = "in")

## effects on net fitness (lambda)

##lambda
ggplot(lifehistorypost,aes(lambda_diff,fill=species,col=species))+
  geom_density(alpha=0.1)+geom_vline(xintercept=0)+
  facet_wrap(~species,scales="free")


lifehistorypost %>% 
  group_by(species) %>% 
  summarise(mean_diff = mean(lambda_diff),
            prop_positive = mean(lambda_diff>0),
            prop_negative = mean(lambda_diff<0),
            x_pos = quantile(lambda_diff[lambda_diff > 0], 0.7, na.rm = TRUE),  
            x_neg = quantile(lambda_diff[lambda_diff < 0], 0.3, na.rm = TRUE)) %>%
              mutate(direction = ifelse(mean_diff >= 0, "positive", "negative"))->proportions_df

lambdadiffplot<-ggplot(lifehistorypost, aes(x = lambda_diff, y = species, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "tomato",
    mid = "white",
    high = "cornflowerblue",
    midpoint = 0,
    limits = c(-.65, .65),  # clip values outside this range
    oob = scales::squish,  # squish outliers into endpoint colors
    name = "Shift"
  ) +  
  geom_text(
    data = proportions_df,
    aes(
      x = x_neg,
      y = as.numeric(factor(species)) + 0.13,  # position just below the ridge
      label = paste0(round(prop_negative, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  # Add positive side (right tail) proportions
  geom_text(
    data = proportions_df,
    aes(
      x = x_pos,
      y = as.numeric(factor(species)) + 0.13,  # same vertical offset
      label = paste0(round(prop_positive, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_minimal(base_size = 13) + theme(legend.position = "none") +
  labs(x = expression("Endophyte effect on host fitness ("*Delta~lambda*")"), y = "Species")+ geom_vline(xintercept = 0)+ 
  geom_point(
    data = proportions_df,
    aes(x = mean_diff, y = species, color = direction),  # now fill is mapped!
    inherit.aes = FALSE,
    shape = 16,
    size = 4
  ) +
  scale_color_manual(
    values = c("positive" = "cornflowerblue", "negative" = "tomato")
  )

ggsave("manuscript/figures/lambdadiffplot.jpg",
       plot = lambdadiffplot,
       width = 6,      # width in inches
       height = 8,      # height in inches
       dpi = 300,       # resolution
       units = "in")

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

                