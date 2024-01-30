



#*************************************#
### Maria Angeles Perez-Navarro
### CREAF  05/2021
### statistical analyses part 1
#*************************************#

library(ggplot2)
library(nlme)
library(lme4)
library(lmerTest)
library(MASS)
library(visreg)
library("DHARMa")
library(tidyverse)
library(effects)
library(lsmeans)
library(viridis)
library(effectsize)
library(readr)
library(MuMIn)
library(parameters)
select=dplyr::select

## needed functions 
source('./niche_modelling/niche_analyses/niche_required_functions.R')



# 1. load datatables ------------------------------------------------------

categ_df <- read_delim("./statistical_analyses/output/Facilitation_ClimaticDisequilibrium_chisq.csv", 
                       delim=",", col_names=T)# verdu chisq analyses in excel


diseq_mode <-  read_delim("./statistical_analyses/output/disequilibrium_mode_centroid.csv", delim=";", col_names=T)
diseq_cent <-  read_delim("./statistical_analyses/output/disequilibrium_nurse.csv", delim=";", col_names=T)
diseq_90 <-  read_delim("./statistical_analyses/output/disequilibrium_perc90_in0.csv", delim=";", col_names=T)
diseq_95 <-  read_delim("./statistical_analyses/output/disequilibrium_perc95_in0.csv", delim=";", col_names=T)
expl <-  read_delim("./statistical_analyses/output/explanatory_vars.csv", delim=";", col_names=T)
change_abun <- read_delim("./statistical_analyses/output/change_abundances.csv", delim=";", col_names=T)

change_abun <- change_abun%>%
  rename(plot=community,
         plot_code=code)

diseq_mode <- diseq_mode%>% tibble(type="max_dens")
diseq_cent <- diseq_cent %>% tibble(type="centroid")
diseq_90 <- diseq_90 %>% tibble(type="90")
diseq_95 <- diseq_95 %>% tibble(type="95")


diseq <- diseq_95 %>% 
  dplyr::bind_rows(diseq_90) %>% 
  dplyr::bind_rows(diseq_cent) %>% 
  dplyr::bind_rows(diseq_mode) %>%
  right_join(expl, by=c("plot", "plot_code",
                        "nurse", "nurse_code"))%>%
  left_join(change_abun, by=c("plot", "plot_code"))



# 2. Chisq analyses and plots ---------------------------------------------


names(categ_df)
categ_df <- categ_df%>%
  mutate(facil_per=as.numeric(`% facilitated`),
         chisq= as.numeric(Chi2))%>%
  mutate(n=`open_obs (individuals)`+
           `canopy_obs (individuals)`)%>%
  filter(n>0)


nrow(categ_df)


chis_mix <- lmer(obs_minus_exp~facilitated_species+
                   (1|Community), data = categ_df)
summary(chis_mix)
emmeans::test(emmeans::emmeans(chis_mix, ~facilitated_species),
              adjust="tukey", null=0)


## add barplot with species x community facilitated or competing basing on verdu test

categ_sum <- categ_df%>%
  group_by(Interaction)%>%
  summarise(count=n())%>%
  mutate(canopy_int=case_when(
    Interaction=="Competition"~"less than expected",
    Interaction=="Facilitation" ~"more than expected",
    Interaction=="Neutral" ~ "expected by chance"
  ))

mycols <- colorRampPalette(viridis_pal(option="plasma")(10))(20)
scales::show_col(mycols)

chisq_barplot <- ggplot(data=categ_sum,aes(x=Interaction, y=count,
                                           fill=canopy_int))+
  geom_bar(stat = "identity")+
  #scale_fill_viridis_d(option="plasma")+
  scale_fill_manual(values=c("#FCCB25",
                             "#CA4778",
                             "#E87158"))+
  ylab("Species x Locality")+
  xlab("Recruit proportion under canopies")+
  scale_x_discrete(labels=c("Competition" = expression("less than \nexpected"),
                            "Facilitation" = expression("more than \n expected"),
                            "Neutral" = expression(" expected \nby chance")))+
  theme_bw()+
  theme(legend.position = "none",
        legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text=element_text(color="black", size=6),
        axis.title=element_text(color="black", size=7),
        axis.title.x = element_text(margin = margin(t = 11, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
        text=element_text(family="serif"))

chisq_barplot

# ggsave(paste0("./statistical_analyses/output/preliminar results/chisq_count.png"), 
#        plot=chisq_barplot, width = 8, height = 7 , dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/preliminar results/chisq_count.pdf"),
#        plot=chisq_barplot, width = 8, height = 7, units="cm")


# 3. Anova models comparing debt for facilitated and open communities -----


## see mean and median per group and method

diseq_sum <- diseq%>%
  group_by(type, nurse)%>%
  summarise(mean_cd=mean(cd),
            median_cd=median(cd),
            sd_cd=sd(cd))


method <- unique(diseq$type)
dif_estimate_df <- tibble(
  type=character(),
  estimate= numeric(),
  pvalue= numeric(),
  r2=numeric(),
  pvalue_paired_t=numeric()
)

for(i in 1:length(method)){
  diseq_type <- diseq%>%
    dplyr::filter(type==method[i])
  
  paired_t <- t.test(scale(cd)~nurse,
                     data=diseq_type,
                     paired=T, alterantive="two.sided")
  
  mixmod <- lmer(scale(cd)~nurse+ (1|plot), data=diseq_type)
  summary(mixmod)
  MuMIn::r.squaredGLMM(mixmod)
  dif_estimate_df[i, "type"] <- method[i]
  dif_estimate_df[i, "estimate"] <- round(summary(mixmod)$coefficients[2],3)
  dif_estimate_df[i, "pvalue"] <- round(summary(mixmod)$coefficients[10], 3)
  dif_estimate_df[i, "r2"] <- round(MuMIn::r.squaredGLMM(mixmod)[1], 3)
  dif_estimate_df[i, "pvalue_paired_t"] <- round(paired_t$p.value, 3)# paired t give exactly the same values as mixed effect model values
}


### plot boxplots anova mixed models all methods together ####

(total_boxplot <- diseq %>% dplyr::filter(type!="max_dens")%>% 
   ggplot(aes(x=type,y=cd,fill=nurse,color=nurse))+
   geom_violin(position = position_dodge(0.8),
               color="transparent",alpha=0.4,width=1)+
   geom_boxplot(position =  position_dodge(0.8),
                fill='transparent',width=0.4, lwd=0.2,
                outlier.size = 0.3)+
   scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                      labels=c("under canopies", "open"))+
   scale_fill_manual(values=c("#E69F00", "#56B4E9"), 
                     labels=c("under canopies", "open"))+
   labs(fill = "Community")+
   scale_x_discrete(labels=c("90" = "90 percentile", "95" = "95 percentile",
                             "centroid" = "Centroid"))+
   guides(color=F)+
   annotate("text", x =c(1:3), y = c(1.8, 2, 2.3),
            label = "*", size=2.5)+
   xlab("Approach")+
   ylab("Climatic disequilibrium")+
   theme_classic()+
   theme(legend.position = c(0.18, 0.9),
         legend.title = element_text(size = 7), 
         legend.text  = element_text(size = 6),
         legend.key.size = unit(3, "mm"),
         # panel.grid.major=element_blank(),
         # panel.grid.minor=element_blank(),
         axis.text=element_text(color="black", size=6),
         axis.title=element_text(color="black", size=7),
         text=element_text(family="serif")))

# ggsave(paste0("./statistical_analyses/output/preliminar results/boxplot nurse/all_methods.png"), 
#        plot=total_boxplot, width = 8, height = 6, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/preliminar results/boxplot nurse/all_methods.pdf"),
#        plot=total_boxplot, width = 8, height = 6, units="cm")
# 



### plot boxplots anova mixed models each method separately ####


diseq_change <- diseq%>%
  pivot_wider(id_cols = c(plot,plot_code, 
                          oc_x, oc_y, oc_bio01, oc_bio04, 
                          oc_bio05, oc_bio06, oc_bio12, 
                          lat, long, altitude, radiation,
                          canopy_percent, gap_percent, site, 
                          annual_PPET_mean, annual_PPET_sd,
                          annual_PPET_sum, type), 
              names_from = nurse, 
              values_from = c("cic_x", "cic_y","cd"))%>%
  mutate(cd_change=abs(cd_facil)-abs(cd_open),
         cd_change_rel=(cd_facil-cd_open)/abs(cd_facil),
         cd_sign=case_when(
           cd_change>=0 ~ "pos",
           cd_change<0 ~ "neg"
         ))

head(diseq)

method <- unique(diseq$type)


# * Figure 3  -------------------------------------------------------------
# create main text figure 3 when i = 1, ie., method=95 percentile


for(i in 1:length(method)){
  
  diseq_type <- diseq%>%
    dplyr::filter(type==method[i])
  
  diseq_change_type <- diseq_change%>%
    dplyr::filter(type==method[i])
  
  mixmod <- lmer(cd~nurse+ (1|plot), data=diseq_type)
  summary(mixmod)
  car::Anova(mixmod)
  MuMIn::r.squaredGLMM(mixmod)
  
  simulationOutput <- simulateResiduals(fittedModel = mixmod, n = 250)
  plot( simulationOutput)
  testZeroInflation(simulationOutput)
  
  c <- lsmeans::lsmeans(mixmod,  ~ nurse)
  lsmeans::contrast(c,  method="pairwise")
  effectsize::effectsize(mixmod)
  (param_tab <- parameters(mixmod))
  param_df <- param_tab%>%
    as.data.frame()%>%
    drop_na()
  (effect_size <- t_to_r(param_df$t, param_df$df_error))# cohen's d
  effect_size <- data.frame(cbind(param_df$Parameter, data.frame(effect_size)))
  names(effect_size)[1] <- "parameter" 
  
  ef_dis <- effects::Effect("nurse", mixmod)# only fixed effects
  eff_df_dis<- data.frame(ef_dis)
  plot(ef_dis)
  
  
  rsq <- round(r.squaredGLMM(mixmod)[[1]], 3)
  pvalue <- round(summary(mixmod)[[10]][2,5], 3)
  
  
  (ggboxpoint <- ggplot(data=diseq_type, aes(x=nurse, y=cd,
                                             fill=nurse,color=nurse))+ 
      geom_violin(position = position_dodge(0.8),
                  color="transparent",alpha=0.4,width=1)+
      geom_boxplot(position =  position_dodge(0.8),
                   fill='transparent',width=0.4, lwd=0.2,
                   outlier.size = 0.3)+
      scale_color_manual(values=c("#E69F00", "#56B4E9"))+
      scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
      geom_point(aes(group=plot), width=0.15, size=0.3,
                 position=position_dodge(0.2)) +
      geom_line(aes(group=plot), 
                color="#999999", lwd=0.2, 
                position=position_dodge(0.2)) +
      annotate("text", x=2.25, y=1.7, 
               label=paste0("R2= ",
                            signif(rsq, 5)), 
               size=2.5, family="serif")+
      annotate("text", x=2.25, y=1.5, 
               label=paste0("Pv=",
                            round(pvalue, 4)), 
               size=2.5, family="serif")+
      xlab(" ")+
      ylab("Climatic disequilibrium")+
      scale_x_discrete(labels=c("under \ncanopies","open"))+
      theme_linedraw()+
      theme(legend.position = "none",
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(color="black", size=8),
            axis.title=element_text(color="black", size=10),
            text=element_text(family="serif"))
  )
  
  
  # (ggbox <- ggplot(data=diseq_type, aes(x=nurse, y=cd,
  #                                      fill=nurse,color=nurse))+ 
  #   geom_violin(position = position_dodge(0.8),
  #               color="transparent",alpha=0.4,width=1)+
  #   geom_boxplot(position =  position_dodge(0.8),
  #                fill='transparent',width=0.4, lwd=0.2,
  #                outlier.size = 0.3)+
  #   scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  #   scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #   geom_line(aes(group=plot), 
  #             color="#999999", lwd=0.2) +
  #   annotate("text", x=2.25, y=2.4, 
  #            label=paste0("R2= ",
  #                         signif(rsq, 5)), 
  #            size=2.5, family="serif")+
  #   annotate("text", x=2.25, y=2.2, 
  #            label=paste0("Pv=",
  #                         round(pvalue, 4)), 
  #            size=2.5, family="serif")+
  #   xlab(" ")+
  #   ylab("Climatic disequilibrium")+
  #   theme_linedraw()+
  #   theme(legend.position = "none",
  #         panel.grid.major=element_blank(),
  #         panel.grid.minor=element_blank(),
  #         axis.text=element_text(color="black", size=8),
  #         axis.title=element_text(color="black", size=10),
  #         text=element_text(family="serif"))
  # )
  
  


  # ggsave(paste0("./statistical_analyses/output/boxplot nurse/", method[i], "_lines.png"), 
  #        plot=ggboxpoint, width = 8, height = 6, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/boxplot nurse/", method[i], "_lines.pdf"),
  #        plot=ggboxpoint, width = 8, height = 6, units="cm")
  
}




# 4. MIXED EFFECT MODEL changes in diseq vs environmental variables -------

i <- 1 #for diseq_type=95 perc

diseq_change <- diseq%>%
  pivot_wider(id_cols = c(plot,plot_code, 
                          oc_x, oc_y, oc_bio01, oc_bio04, 
                          oc_bio05, oc_bio06, oc_bio12, lat, long,
                          canopy_percent, gap_percent, site, 
                          annual_PPET_mean, annual_PPET_sd,
                          annual_PPET_sum, altitude, radiation,
                          type, names(change_abun)[3:19]), 
              names_from = nurse, 
              values_from = c("cic_x", "cic_y","cd"))%>%
  mutate(aridity=1/annual_PPET_sum,
         abs_oc_bio06=abs(oc_bio06-max(oc_bio06)),
         cd_dif=cd_facil-cd_open,# it is the same as cd change as 2d climatic debt is always positive
         cd_change=abs(cd_facil)-abs(cd_open),
         cd_ratio=abs(cd_facil)/abs(cd_open),# use abs values for 2d climatic debt not univariate
         cd_change_rel=(cd_facil-cd_open)/abs(cd_facil),
         ai=1/annual_PPET_sum,
         cd_sign=case_when(
           cd_change>=0 ~ "pos",
           cd_change<0 ~ "neg"
         ))

head(diseq_change)

method <- unique(diseq_change$type)
diseq_c <- diseq_change%>%
  dplyr::filter(type==method[i])

## def selected models

diseq <- diseq%>%
  mutate(ai=1/annual_PPET_sum)

# run i = 1 for results with percentile 95

for(i in 1:length(method)){
  

diseq_type <- diseq%>%
  dplyr::filter(type==method[i])

aeet_mod <- lmer(cd~nurse*#scale(annual_PPET_sum) + 
                   nurse*scale(ai)+
                   nurse*scale(radiation)+#nurse*annual_PPET_sd+
                   nurse*scale(oc_bio01)+
                   nurse*scale(canopy_percent)+
                   (1|plot_code), data=diseq_type)# run

summary(aeet_mod)
vif <- performance::check_collinearity(aeet_mod)#run it without interactions
vif

car::Anova(aeet_mod, type="II")
MuMIn::r.squaredGLMM(aeet_mod)
MuMIn::AICc(aeet_mod)

expl_vars <- c("radiation",
               #"annual_PPET_sum",
               "ai",
               "oc_bio01",
               "canopy_percent")

aeet_slopes <- expl_vars %>% 
  purrr::map(function(x){
    emmeans::emtrends(aeet_mod, "nurse", var = x) %>% 
      emmeans::test()%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
      as_tibble() %>% 
      rename(slope = 2) %>%
      mutate(variable = x,
             ci_upper= slope+1.96*SE,
             ci_lower= slope-1.96*SE)
  }) %>% 
  bind_rows()


aeet_pairs <- expl_vars %>% 
  purrr::map(function(x){
    emmeans::emtrends(aeet_mod, "nurse", var = x) %>% 
      pairs()%>%
      as_tibble() %>% 
      rename(trend_dif = 2,# add also facil open column in case they are reverse for some continuous variables
             pvalue_dif=6) %>%
      mutate(variable = x)%>%
      dplyr::select(variable,
                    trend_dif,
                    pvalue_dif)
  }) %>% 
  bind_rows()


aeet_slopes <- aeet_slopes%>%
  mutate(ef_size=effectsize::t_to_d(t.ratio, 
                                    nrow(aeet_mod@frame))$d,#use nrow(mod@frame)instead nrow(data_table)in case na are authomatically removed in the model
         ef_ci_lower=effectsize::t_to_d(t.ratio, 
                                        nrow(aeet_mod@frame))$CI_low,
         ef_ci_upper=effectsize::t_to_d(t.ratio, 
                                        nrow(aeet_mod@frame))$CI_high,
         r2m=MuMIn::r.squaredGLMM(aeet_mod)[1],
         r2c=MuMIn::r.squaredGLMM(aeet_mod)[2])

aeet_df<- aeet_slopes%>%
  left_join(aeet_pairs, by="variable")


aeet_df <- aeet_df%>%
    mutate(nurse=as.character(nurse),
           var_name=case_when(
      variable=="radiation"~"Radiation",
      variable=="ai"~"Aridity index",
      #variable=="annual_PPET_sum"~"PPET",
      variable=="oc_bio01"~"Annual mean \ntemperature  ",
      variable=="canopy_percent"~"Canopy percent"),
      group=case_when(
      nurse=="facil"~"Under canopies",
      nurse=="open"~"Open"),
      significance=case_when(
      p.value<0.05~"signif",
      p.value>0.05~"no signif")
  )%>%
  mutate(nursexsig = case_when(
           significance == "no signif"~"no signif",
           significance != "no signif"~ nurse))


# data.table::fwrite(aeet_df, paste0("./statistical_analyses/output/diseq change as response/model_df_",
#                                    method[i],".csv"),
#                    sep=";", dec=".",
#                    col.names=T, row.names=F)


aeet_df$var_name <- factor(aeet_df$var_name, levels=c("Annual mean \ntemperature  ",
                                                      "Aridity index",
                                                      "Canopy percent",
                                                      "Radiation"))

(effsize_plot <- aeet_df%>%
    mutate(sig = symnum(pvalue_dif,corr = FALSE,na = FALSE,
                        cutpoints = c(0,0.05, 1),
                        symbols = c("*", " "))
    ) %>% 
    #arrange(var_name)%>%
    ggplot(aes(x= ef_size, y = var_name,
               color= nurse
               #group=interaction(origin4,var_name)
    ))+
    geom_vline(xintercept= 0, 
               linetype=2, 
               color="grey30", 
               size=0.15)+
    geom_point(aes( 
      fill = nursexsig), 
      size=1.3, 
      stroke=0.4,
      position = position_dodge(0.3, preserve="total"),
      shape=21
    )+
    geom_errorbarh(aes(y = var_name, 
                       xmin = ef_ci_lower,
                       xmax = ef_ci_upper,
                       color=nurse),
                   height=0,
                   size=0.25,
                   position = position_dodge(0.3, preserve="total"),
                   # alpha=0.7
    )+
    scale_fill_manual(values= c("white",
                                "#E69F00",
                                "#56B4E9"
    ))+
    scale_color_manual(values=c("#E69F00", "#56B4E9"),
                       labels=c("Under canopy", "Open"))+
    geom_text(aes(x= max(ef_ci_upper)+0.05, y = var_name,label=sig),
              vjust=0.7, size=3, family="serif",
              color="black")+
    # geom_text(data=dat_text,
    #           mapping=aes(x=x, y=y, label=label),
    #           hjust=0, size=2.6, family="serif")+
    guides(fill=F)+
    xlab("Effect size")+
    ylab(" ")+
    theme_bw()+
    theme(axis.title.y=element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(color="black", size=8),
          axis.title.x = element_text(color="black", size=8),
          axis.ticks= element_line(size=0.4),
          legend.title= element_text(color="black", size=8),
          legend.text=element_text(color="black",size=7),
          legend.key.size = unit(0.4, "mm"),
          legend.key.height = unit(3.5, "mm"),
          legend.key.width = unit(1.2, "mm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          text=element_text(family="serif"),
          legend.box="vertical", 
          panel.border=element_rect(colour="black", size=0.25)
          
    ))


# ggsave(paste0("./statistical_analyses/output/diseq change as response/effectsize_plot_", method[i], "ai.png"), 
#        plot=effsize_plot, width = 8, height = 6, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/effectsize_plot_", method[i], "ai.pdf"),
#        plot=effsize_plot, width = 8, height = 6, units="cm")
# 

#**** REGRESSION PLOTS ********

# Figure 5 ----------------------------------------------------------------
#* main text figure 5. i=1 for centroid estimated as 95 percentile


mixmod <- aeet_mod

summary(mixmod)
MuMIn::r.squaredGLMM(mixmod)

mixmod_coeff <- as.data.frame(summary(mixmod)$coefficients)

pv <- round(mixmod_coeff$`Pr(>|t|)`[7], 3)
r2 <- round(MuMIn::r.squaredGLMM(mixmod)[1], 3)


ef <- effects::Effect(c("nurse","ai"), mixmod, 25)
eff_df<- data.frame(ef)

(nurse_ppet <-  ggplot()+
    geom_point(data=diseq_type, aes(x=ai, y=cd, color=nurse),
               size=1, shape=19)+
    # ggrepel::geom_label_repel(data=diseq_type,
    #                           aes(x=annual_PPET_sum, y=abs(cd),
    #                               label = plot),
    #                           box.padding   = 0.15,
    #                           point.padding = 0.15,
    #                           label.padding = 0.13,
    #                           segment.color = 'grey50',
    #                           max.overlaps = 50, size=1.3)+
    geom_line(data=eff_df, aes(x=ai,
                               y= fit,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower, ymax=upper, 
                                 x=ai, fill = nurse), alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    # annotate("text",  x= 0.43, y= max(diseq_type$cd)-0.05,
    #          label=paste0("Pvalue = ",pv), 
    #          size=2, hjust=0, vjust="top",  family="serif")+
    # annotate("text",  x= 0.43, y= max(diseq_type$cd)-0.18, 
    #          label=paste0("R^2 == ",r2), 
    #          parse=TRUE, size=2, hjust=0, vjust="top", family="serif")+
    ylab("Multivariate Climatic Disequilibrium")+
    xlab("Aridity Index (PET/P)")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title= element_text(color="black", size=7),
          legend.text=element_text(color="black",size=6),
          legend.key.size = unit(0.4, "mm"),
          legend.key.height = unit(3.5, "mm"),
          legend.key.width = unit(1.2, "mm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text=element_text(color="black", size=7),
      axis.title=element_text(color="black", size=8),
      text=element_text(family="serif"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      axis.ticks = element_line(colour = "black", size = 0.2))
)


# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_vs_nurse_ai_", method[i],".png"),
#        plot=nurse_ppet, width = 8, height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_vs_nurse_ai_", method[i],".pdf"),
#        plot=nurse_ppet, width = 8, height = 6.3, units="cm")


ef <- effects::Effect(c("nurse","canopy_percent"), mixmod, 25)# only fixed effects
eff_df<- data.frame(ef)

pv <- round(mixmod_coeff$`Pr(>|t|)`[10], 3)
r2 <- round(MuMIn::r.squaredGLMM(mixmod)[1], 3)

(nurse_gap <-  ggplot()+
    geom_point(data=diseq_type, aes(x=canopy_percent, y=cd, color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df, aes(x=canopy_percent,
                               y= fit,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower, ymax=upper, 
                                 x=canopy_percent, fill = nurse), alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    # annotate("text",  x= 30, y= max(diseq_type$cd)-0.05, 
    #          label=paste0("Pvalue = ",pv), 
    #          size=2, hjust=0, vjust="top",  family="serif")+
    # annotate("text",  x= 30, y= max(diseq_type$cd)-0.18, 
    #          label=paste0("R^2 == ",r2), 
    #          parse=TRUE, size=2, hjust=0, vjust="top", family="serif")+
    ylab("Climatic Disequilibrium")+
    xlab("Plot canopy percentage")+
    theme_bw()+
    theme(legend.position = "bottom",
      legend.title= element_text(color="black", size=7),
      legend.text=element_text(color="black",size=6),
      legend.key.size = unit(0.4, "mm"),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(1.2, "mm"),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text=element_text(color="black", size=7),
      axis.title=element_text(color="black", size=8),
      text=element_text(family="serif"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      axis.ticks = element_line(colour = "black", size = 0.2))
)


# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_vs_nurse_canopypercent_", method[i],".png"),
#        plot=nurse_gap, width = 8, height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_vs_nurse_canopypercent_", method[i],".pdf"),
#        plot=nurse_gap, width = 8, height = 6.3, units="cm")
# 
# 


ef <- effects::Effect(c("nurse","oc_bio01"), mixmod, 25)# only fixed effects
eff_df<- data.frame(ef)

pv <- round(mixmod_coeff$`Pr(>|t|)`[10], 3)
r2 <- round(MuMIn::r.squaredGLMM(mixmod)[1], 3)

(nurse_temp <-  ggplot()+
    geom_point(data=diseq_type, aes(x=oc_bio01, y=cd, color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df, aes(x=oc_bio01,
                               y= fit,color=nurse),
              size=0.5,alpha=1)+
    # geom_line(data=eff_df, aes(gap_percent,
    #                            lower,color=nurse),
    #           linetype=2, size=0.8,
    #           alpha=1)+
    # geom_line(data=eff_df, aes(gap_percent,
    #                            upper,color=nurse),
    #           linetype=2, size=0.8,
    #           alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower, ymax=upper, 
                                 x=oc_bio01, fill = nurse), alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    # annotate("text",  x= 30, y= max(diseq_type$cd)-0.05, 
    #          label=paste0("Pvalue = ",pv), 
    #          size=2, hjust=0, vjust="top",  family="serif")+
    # annotate("text",  x= 30, y= max(diseq_type$cd)-0.18, 
    #          label=paste0("R^2 == ",r2), 
    #          parse=TRUE, size=2, hjust=0, vjust="top", family="serif")+
    ylab("Climatic Disequilibrium")+
    xlab("Annual mean temperature")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title= element_text(color="black", size=7),
          legend.text=element_text(color="black",size=6),
          legend.key.size = unit(0.4, "mm"),
          legend.key.height = unit(3.5, "mm"),
          legend.key.width = unit(1.2, "mm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text=element_text(color="black", size=7),
      axis.title=element_text(color="black", size=8),
      text=element_text(family="serif"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      axis.ticks = element_line(colour = "black", size = 0.2))
)

# 
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_vs_nurse_annualtemp_",method[i],".png"),
#        plot=nurse_temp, width = 8, height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_vs_nurse_annualtemp_",method[i],".pdf"),
#        plot=nurse_temp, width = 8, height = 6.3, units="cm")
# 
# 



ef <- effects::Effect(c("nurse","radiation"), mixmod, 25)# only fixed effects
eff_df<- data.frame(ef)

pv <- round(mixmod_coeff$`Pr(>|t|)`[10], 3)
r2 <- round(MuMIn::r.squaredGLMM(mixmod)[1], 3)

(nurse_rad <-  ggplot()+
    geom_point(data=diseq_type, aes(x=radiation, y=cd, color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df, aes(x=radiation,
                               y= fit,color=nurse),
              size=0.5,alpha=1)+
    # geom_line(data=eff_df, aes(gap_percent,
    #                            lower,color=nurse),
    #           linetype=2, size=0.8,
    #           alpha=1)+
    # geom_line(data=eff_df, aes(gap_percent,
    #                            upper,color=nurse),
    #           linetype=2, size=0.8,
    #           alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower, ymax=upper, 
                                 x=radiation, fill = nurse), alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    # annotate("text",  x= 30, y= max(diseq_type$cd)-0.05, 
    #          label=paste0("Pvalue = ",pv), 
    #          size=2, hjust=0, vjust="top",  family="serif")+
    # annotate("text",  x= 30, y= max(diseq_type$cd)-0.18, 
    #          label=paste0("R^2 == ",r2), 
    #          parse=TRUE, size=2, hjust=0, vjust="top", family="serif")+
    ylab("Climatic Disequilibrium")+
    xlab("Radiation")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title= element_text(color="black", size=7),
          legend.text=element_text(color="black",size=6),
          legend.key.size = unit(0.4, "mm"),
          legend.key.height = unit(3.5, "mm"),
          legend.key.width = unit(1.2, "mm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text=element_text(color="black", size=7),
      axis.title=element_text(color="black", size=8),
      text=element_text(family="serif"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      axis.ticks = element_line(colour = "black", size = 0.2))
)


# ggsave(paste0("./statistical_analyses/output/preliminar results/diseq change as response/diseq_vs_nurse_radiation_", method[i], ".png"),
#        plot=nurse_rad, width = 8, height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/preliminar results/diseq change as response/diseq_vs_nurse_radiation_", method[i], ".pdf"),
#        plot=nurse_rad, width = 8, height = 6.3, units="cm")

}




# 5. MODEL ACCOUNTING SPATIAL AUTOCORRELATION -----------------------------


diseq_change <- diseq%>%
  pivot_wider(id_cols = c(plot,plot_code, 
                          oc_x, oc_y, oc_bio01, oc_bio04, 
                          oc_bio05, oc_bio06, oc_bio12, lat, long,
                          canopy_percent, gap_percent, site, 
                          annual_PPET_mean, annual_PPET_sd,
                          annual_PPET_sum, altitude, radiation,
                          type, names(change_abun)[3:19]), 
              names_from = nurse, 
              values_from = c("cic_x", "cic_y","cd"))%>%
  mutate(aridity=1/annual_PPET_sum,
         abs_oc_bio06=abs(oc_bio06-max(oc_bio06)),
         cd_dif=cd_facil-cd_open,# it is the same as cd change as 2d climatic debt is always positive
         cd_change=abs(cd_facil)-abs(cd_open),
         cd_ratio=abs(cd_facil)/abs(cd_open),# use abs values for 2d climatic debt not univariate
         cd_change_rel=(cd_facil-cd_open)/abs(cd_facil),
         ai=1/annual_PPET_sum,
         cd_sign=case_when(
           cd_change>=0 ~ "pos",
           cd_change<0 ~ "neg"
         ))

head(diseq_change)

method <- unique(diseq_change$type)
diseq_c <- diseq_change%>%
  dplyr::filter(type=="95")

## accounting for spatial autocorrelation using gls and cd_change to avoid using mixed effects

# check some advices here https://www.flutterbys.com.au/stats/tut/tut8.4a.html

aeet_mod0gls <- gls(cd_change~scale(ai )+ 
                      scale(radiation)+ scale(oc_bio01)+
                      scale(canopy_percent), data=diseq_c)
d <- diseq_c
d$res <- as.numeric(residuals(aeet_mod0gls))
ggplot(d, aes(long, lat, colour = res)) + 
  geom_point(size = 5) + scale_color_gradient2()

plot(Variogram(aeet_mod0gls, form= ~ lat+long, data=diseq_c))

aeet_mod1gls <- update(aeet_mod0gls, 
                       corr = corSpher(c(3, 0.9), 
                                       form = ~ lat + long, 
                                       nugget = TRUE))# range and nugget value given the semivariogram. It seems that spatial autorrelation is too little

aeet_mod2gls <- gls(cd_change~scale(annual_PPET_sum )+ 
                      scale(radiation)+ scale(oc_bio01)+
                      scale(canopy_percent), data=diseq_c, 
                    correlation = corSpatial(form = ~lat+ long, 
                                             type ="spherical", nugget = T))

aeet_mod3gls <- gls(cd_change~scale(annual_PPET_sum )+ 
                      scale(radiation)+ scale(oc_bio01)+
                      scale(canopy_percent), data=diseq_c, 
                    correlation = corSpatial(form = ~lat+ long, 
                                             type ="gaussian", nugget = T))

aeet_mod4gls <- gls(cd_change~scale(annual_PPET_sum )+ 
                      scale(radiation)+ scale(oc_bio01)+
                      scale(canopy_percent), data=diseq_c, 
                    correlation = corSpatial(form = ~lat+ long, 
                                             type ="exponential", nugget = T))

aeet_mod5gls <- gls(cd_change~scale(annual_PPET_sum )+ 
                      scale(radiation)+ scale(oc_bio01)+
                      scale(canopy_percent), data=diseq_c, 
                    correlation = corSpatial(form = ~lat+ long, 
                                             type ="rational", nugget = T))

summary(aeet_mod1gls)
str(summary(aeet_mod1gls))
rsq::rsq(aeet_mod1gls)
R2 <- cor(diseq_c$cd_change,predict(aeet_mod1))^2
R2
R2.1 <- 1 - with(diseq_c, (sum((cd_change-predict(aeet_mod1))^2)/sum((cd_change-mean(cd_change))^2)))
R2.1

summary(aeet_mod1)$coefficients


AIC(aeet_mod0gls,
    aeet_mod1gls,
    aeet_mod2gls,
    aeet_mod3gls,
    aeet_mod4gls,
    aeet_mod5gls)


anova(aeet_mod0gls,
    aeet_mod1gls,
    aeet_mod2gls,
    aeet_mod3gls,
    aeet_mod4gls,
    aeet_mod5gls)




##*************************************************************************
# 6. UNIVARIATE MODELS AND BARPLOTS ---------------------------------------
##*************************************************************************


## load multivariate

diseq_mode <-  read_delim("./statistical_analyses/output/disequilibrium_mode_centroid.csv", delim=";", col_names=T)
diseq_cent <-  read_delim("./statistical_analyses/output/disequilibrium_nurse.csv", delim=";", col_names=T)
diseq_90 <-  read_delim("./statistical_analyses/output/disequilibrium_perc90_in0.csv", delim=";", col_names=T)
diseq_95 <-  read_delim("./statistical_analyses/output/disequilibrium_perc95_in0.csv", delim=";", col_names=T)

diseq_mode <- diseq_mode%>% tibble(type="max_dens")
diseq_cent <- diseq_cent %>% tibble(type="centroid")
diseq_90 <- diseq_90 %>% tibble(type="90")
diseq_95 <- diseq_95 %>% tibble(type="95")


diseq_multi <- diseq_95 %>% 
  dplyr::bind_rows(diseq_90) %>% 
  dplyr::bind_rows(diseq_cent) %>% 
  dplyr::bind_rows(diseq_mode) 

diseq_multi <- diseq_multi%>%
  mutate(cd_axis1=cic_x-oc_x,
         cd_axis2=cic_y-oc_y)

## load univariate 
diseq_mode_u <-  read_delim("./statistical_analyses/output/disequilibrium_centroid_univariate.csv", delim=";", col_names=T)
diseq_cent_u <-  read_delim("./statistical_analyses/output/disequilibrium_mode_univariate.csv", delim=";", col_names=T)
diseq_90_u <-  read_delim("./statistical_analyses/output/disequilibrium_perc_90_in0_univariate.csv", delim=";", col_names=T)
diseq_95_u <-  read_delim("./statistical_analyses/output/disequilibrium_perc_95_in0_univariate.csv", delim=";", col_names=T)

diseq_mode_u <- diseq_mode_u %>% tibble(type="max_dens")
diseq_cent_u <- diseq_cent_u %>% tibble(type="centroid")
diseq_90_u <- diseq_90_u %>% tibble(type="90")
diseq_95_u <- diseq_95_u %>% tibble(type="95")


diseq_uni <- diseq_95_u %>% 
  dplyr::bind_rows(diseq_90_u) %>% 
  dplyr::bind_rows(diseq_cent_u) %>% 
  dplyr::bind_rows(diseq_mode_u) 

method <- unique(diseq_uni$type)

#i=1 for niche method 95 percentile

for(i in 1:length(method)){
  
  ## cd_abs is estimated but it is inacuarate as it is estimated as abs(cic)-abs(oc) and this make no sense
  
  multiv_total <- diseq_multi%>%
    dplyr::filter(type==method[i])
  
  multiv_cd0 <- multiv_total%>%
    dplyr::select(-c(oc_x, oc_y,
              cic_x, cic_y))%>%
    rename(multiv=cd,
           axis1=cd_axis1,
           axis2=cd_axis2)%>%
    pivot_longer(!c(plot,
                    plot_code,
                    nurse,
                    nurse_code,
                    type),names_to = "variable",
                 values_to ="cd")
  
  multiv_cic <- multiv_total%>%
    dplyr::select(-c(oc_x, oc_y,
              cd, cd_axis1,
              cd_axis2))%>%
    mutate(multiv= NA_integer_)%>%
    rename(axis1=cic_x,
           axis2=cic_y)%>%
    pivot_longer(!c(plot,
                    plot_code,
                    nurse,
                    nurse_code,
                    type),names_to = "variable",
                 values_to ="cic")
  
  multiv_oc <- multiv_total%>%
    dplyr::select(-c(cic_x, cic_y,
              cd, cd_axis1,
              cd_axis2))%>%
    mutate(multiv= NA_integer_)%>%
    rename(axis1=oc_x,
           axis2=oc_y)%>%
    pivot_longer(!c(plot,
                    plot_code,
                    nurse,
                    nurse_code,
                    type),names_to = "variable",
                 values_to ="oc")
  
  multiv_all <- purrr::reduce(list(multiv_oc,
                                   multiv_cic,
                                   multiv_cd0),
                              dplyr::left_join, by=c("plot",
                                                    "plot_code",
                                                    "nurse",
                                                    "nurse_code",
                                                    "type",
                                                    "variable"))
  
  multiv_cd <- multiv_all%>%
    mutate(cd_norm = case_when(
             variable!="multiv"~cd/abs(oc),
             variable=="multiv"~cd
           ),
           cd_abs=abs(cd))%>%
    dplyr::select(
      "plot",  "plot_code", "nurse", "nurse_code",
      "variable", "oc", "cic", "cd", "cd_norm", "cd_abs"  
    )
  
  univ_total <- diseq_uni%>%
    dplyr::filter(type==method[i])
  
  univ_oc <- univ_total%>%
    dplyr::select(1:23)%>%
    rename(bio01= oc_bio01,
           bio02= oc_bio02,
           bio03= oc_bio03,
           bio04= oc_bio04,
           bio05= oc_bio05,
           bio06= oc_bio06,
           bio07= oc_bio07,
           bio08= oc_bio08,
           bio09= oc_bio09,
           bio10= oc_bio10,
           bio11= oc_bio11,
           bio12= oc_bio12,
           bio13= oc_bio13,
           bio14= oc_bio14,
           bio15= oc_bio15,
           bio16= oc_bio16,
           bio17= oc_bio17,
           bio18= oc_bio18,
           bio19= oc_bio19)%>%
    tidyr::gather(key="variable", value="oc", 
                  -c(plot, plot_code, nurse, nurse_code))%>%
    dplyr::filter(variable %in% c("bio01", "bio04",
                                  "bio05", "bio06",
                                  "bio07", "bio10",
                                  "bio11", "bio12",
                                  "bio13", "bio14",
                                  "bio15", "bio16",
                                  "bio17"))
  univ_cic <- univ_total%>%
    dplyr::select(1:4,contains("cic"))%>%
    rename(bio01= cic_bio01,
           bio02= cic_bio02,
           bio03= cic_bio03,
           bio04= cic_bio04,
           bio05= cic_bio05,
           bio06= cic_bio06,
           bio07= cic_bio07,
           bio08= cic_bio08,
           bio09= cic_bio09,
           bio10= cic_bio10,
           bio11= cic_bio11,
           bio12= cic_bio12,
           bio13= cic_bio13,
           bio14= cic_bio14,
           bio15= cic_bio15,
           bio16= cic_bio16,
           bio17= cic_bio17,
           bio18= cic_bio18,
           bio19= cic_bio19)%>%
    tidyr::gather(key="variable", value="cic", 
                  -c(plot, plot_code, nurse, nurse_code))%>%
    dplyr::filter(variable %in% c("bio01", "bio04",
                                  "bio05", "bio06",
                                  "bio07", "bio10",
                                  "bio11", "bio12",
                                  "bio13", "bio14",
                                  "bio15", "bio16",
                                  "bio17"))
  
  univ_cd <- univ_total%>%
    dplyr::select(1:4, contains("cd_b"))%>%
    rename(bio01= cd_bio01,
           bio02= cd_bio02,
           bio03= cd_bio03,
           bio04= cd_bio04,
           bio05= cd_bio05,
           bio06= cd_bio06,
           bio07= cd_bio07,
           bio08= cd_bio08,
           bio09= cd_bio09,
           bio10= cd_bio10,
           bio11= cd_bio11,
           bio12= cd_bio12,
           bio13= cd_bio13,
           bio14= cd_bio14,
           bio15= cd_bio15,
           bio16= cd_bio16,
           bio17= cd_bio17,
           bio18= cd_bio18,
           bio19= cd_bio19)%>%
    tidyr::gather(key="variable", value="cd", 
                  -c(plot, plot_code, nurse, nurse_code))%>%
    dplyr::filter(variable %in% c("bio01", "bio04",
                                  "bio05", "bio06",
                                  "bio07", "bio10",
                                  "bio11", "bio12",
                                  "bio13", "bio14",
                                  "bio15", "bio16",
                                  "bio17"))
  
  univ_cd_abs <- univ_total%>%
    dplyr::select(1:4, contains("cd_abs_"))%>%
    rename(bio01= cd_abs_bio01,
           bio02= cd_abs_bio02,
           bio03= cd_abs_bio03,
           bio04= cd_abs_bio04,
           bio05= cd_abs_bio05,
           bio06= cd_abs_bio06,
           bio07= cd_abs_bio07,
           bio08= cd_abs_bio08,
           bio09= cd_abs_bio09,
           bio10= cd_abs_bio10,
           bio11= cd_abs_bio11,
           bio12= cd_abs_bio12,
           bio13= cd_abs_bio13,
           bio14= cd_abs_bio14,
           bio15= cd_abs_bio15,
           bio16= cd_abs_bio16,
           bio17= cd_abs_bio17,
           bio18= cd_abs_bio18,
           bio19= cd_abs_bio19)%>%
    tidyr::gather(key="variable", value="cd_abs", 
                  -c(plot, plot_code, nurse, nurse_code))%>%
    dplyr::filter(variable %in% c("bio01", "bio04",
                                  "bio05", "bio06",
                                  "bio07", "bio10",
                                  "bio11", "bio12",
                                  "bio13", "bio14",
                                  "bio15", "bio16",
                                  "bio17"))
  
  
  
  univ_cdoc <- univ_cd%>%
    left_join(univ_cd_abs, by=c("plot", "plot_code", 
                                "nurse", "nurse_code",
                                "variable"))%>%
    left_join(univ_oc, by= c("plot", "plot_code", 
                             "nurse", "nurse_code",
                             "variable"))%>%
    left_join(univ_cic, 
              by= c("plot", "plot_code", 
                    "nurse", "nurse_code",
                    "variable"))%>%
    mutate(cd_norm=cd/abs(oc))%>%
    dplyr::select(plot, plot_code, nurse, 
                  nurse_code, variable, 
                  oc,cic, cd, cd_norm, cd_abs)
  
  total_cd <- rbind(multiv_cd, univ_cdoc)
  
  total_cd <- total_cd%>%
    mutate(var_name=case_when(
      variable=="axis1"~"Axis1",
      variable=="axis2"~"Axis2",
      variable=="bio01"~"Annual \nTemp",
      variable=="bio04"~"Temp \nseason",
      variable=="bio05"~"Max \nTemp WM",
      variable=="bio06"~"Min \nTemp CM",
      variable=="bio07"~"Temp \nAnn range",
      variable=="bio10"~"Mean \nTemp WQ",
      variable=="bio11"~"Mean \nTemp CQ",
      variable=="bio12"~"Annual \nPrec",
      variable=="bio13"~"Prec WM",
      variable=="bio14"~"Prec DM",
      variable=="bio15"~"Prec \nSeason",
      variable=="bio16"~"Prec WQ",
      variable=="bio17"~"Prec DQ",
      variable=="multiv"~"Multiv"
    ),
    var_name1=case_when(
      variable=="axis1"~"Axis1",
      variable=="axis2"~"Axis2",
      variable=="bio01"~"Annual Temp",
      variable=="bio04"~"Temp season",
      variable=="bio05"~"Max Temp WM",
      variable=="bio06"~"Min Temp CM",
      variable=="bio07"~"Temp Ann range",
      variable=="bio10"~"Mean Temp WQ",
      variable=="bio11"~"Mean Temp CQ",
      variable=="bio12"~"Annual Prec",
      variable=="bio13"~"Prec WeM",
      variable=="bio14"~"Prec DM",
      variable=="bio15"~"Prec Season",
      variable=="bio16"~"Prec WeQ",
      variable=="bio17"~"Prec DQ",
      variable=="multiv"~"Multiv"
    ))
  
  total_sum <- total_cd%>%
    group_by(nurse, var_name1)%>%
    summarise(mean_cic=mean(cic),
              se_cic=sd(cic),
              mean_cd=mean(cd),
              se_cd= sd(cd),
              mean_cd_abs=mean(abs(cd)),
              se_cd_abs=sd(abs(cd)))%>%
    arrange(var_name1)
 
  # write.table(total_sum, paste0("./statistical_analyses/output/mean_cd_univariate_",
  #             method[i],".csv"), sep=";", dec=".", col.names=T, row.names=F, quote=F)
  # 
  
  biovariable <- sort(unique(total_cd$variable))
  
  cd_estimate_univ <- tibble()
  
  
  for(j in 1:length(biovariable)){
    
    diseq_bio <-  total_cd%>%
      dplyr::filter(variable==biovariable[j])
    
    # (ggboxpoint <- ggplot(data=diseq_bio, aes(x=nurse, y=cd,
    #                                            fill=nurse,color=nurse))+ 
    #     geom_violin(position = position_dodge(0.8),
    #                 color="transparent",alpha=0.4,width=1)+
    #     geom_boxplot(position =  position_dodge(0.8),
    #                  fill='transparent',width=0.4, lwd=0.2,
    #                  outlier.size = 0.3)+
    #     scale_color_manual(values=c("#E69F00", "#56B4E9"))+
    #     scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    #     geom_point(aes(group=plot), width=0.15, size=0.3,
    #                position=position_dodge(0.2)) +
    #     geom_line(aes(group=plot), 
    #               color="#999999", lwd=0.2, 
    #               position=position_dodge(0.2)) +
    #     xlab(" ")+
    #     ylab("Climatic disequilibrium")+
    #     scale_x_discrete(labels=c("under \ncanopies","open"))+
    #     theme_linedraw()+
    #     theme(legend.position = "none",
    #           panel.grid.major=element_blank(),
    #           panel.grid.minor=element_blank(),
    #           axis.text=element_text(color="black", size=8),
    #           axis.title=element_text(color="black", size=10),
    #           text=element_text(family="serif"))
    # )
    
    if(biovariable[j]!="multiv"){
      mixmod_cic <- lmer(cic~nurse+ (1|plot), data=diseq_bio)
      summary(mixmod_cic)
      MuMIn::r.squaredGLMM(mixmod_cic)
      
      simulationOutput_cic <- simulateResiduals(fittedModel = mixmod_cic, n = 250)
      plot( simulationOutput_cic)
      hist(resid(mixmod_cic))
    }else if(biovariable[j]=="multiv"){
      mixmod_cic <- NA
    }
    
    mixmod_cd <- lmer(cd~nurse+ (1|plot), data=diseq_bio)
    summary(mixmod_cd)
    MuMIn::r.squaredGLMM(mixmod_cd)
    
    simulationOutput_cd <- simulateResiduals(fittedModel = mixmod_cd, n = 250)
    plot( simulationOutput_cd)
    hist(resid(mixmod_cd))
    
    mixmod_cd_norm <- lmer(cd_norm~nurse+ (1|plot), data=diseq_bio)
    summary(mixmod_cd_norm)
    MuMIn::r.squaredGLMM(mixmod_cd_norm)
    
    simulationOutput_cd_norm <- simulateResiduals(fittedModel = mixmod_cd_norm, n = 250)
    plot( simulationOutput_cd_norm)
    hist(resid(mixmod_cd_norm))
    
    mixmod_cd_abs <- lmer(cd_abs~nurse+ (1|plot), data=diseq_bio)
    summary(mixmod_cd_abs)
    MuMIn::r.squaredGLMM(mixmod_cd_abs)
    
    simulationOutput_cd_abs <- simulateResiduals(fittedModel = mixmod_cd_abs, n = 250)
    plot( simulationOutput_cd_abs)
    hist(resid(mixmod_cd_abs))
    
    mixmod_cd_abs2 <- lmer(abs(cd)~nurse+ (1|plot), data=diseq_bio)
    summary(mixmod_cd_abs2)
    MuMIn::r.squaredGLMM(mixmod_cd_abs2)
    
    simulationOutput_cd_abs2 <- simulateResiduals(fittedModel = mixmod_cd_abs2, n = 250)
    plot( simulationOutput_cd_abs2)
    hist(resid(mixmod_cd_abs2))
    
    cd_estimate_univ[j, "row_number"] <- j
    cd_estimate_univ[j, "variable"] <- biovariable[j]
   
    if(biovariable[j]!="multiv"){
      cd_estimate_univ[j, "estimate_cic"] <- round(summary(mixmod_cic)$coefficients[2], 3)
      cd_estimate_univ[j, "se_cic"] <- round(summary(mixmod_cic)$coefficients[4], 3)
      cd_estimate_univ[j, "pvalue_cic"] <- round(summary(mixmod_cic)$coefficients[10], 3)
      cd_estimate_univ[j, "r2_cic"] <- round(MuMIn::r.squaredGLMM(mixmod_cic)[1], 3)
    }else if(biovariable[j]=="multiv"){
      cd_estimate_univ[j, "estimate_cic"] <- NA_real_
      cd_estimate_univ[j, "se_cic"] <- NA_real_
      cd_estimate_univ[j, "pvalue_cic"] <- 1
      cd_estimate_univ[j, "r2_cic"] <- NA_real_
    }
    
    cd_estimate_univ[j, "estimate_cd"] <- round(summary(mixmod_cd)$coefficients[2], 3)# negative estimates do not means that facil has lower disequilibrium but more negative disequilibrium values, ie. lower cic in relation to oc
    cd_estimate_univ[j, "se_cd"] <- round(summary(mixmod_cd)$coefficients[4], 3)
    cd_estimate_univ[j, "pvalue_cd"] <- round(summary(mixmod_cd)$coefficients[10], 3)
    cd_estimate_univ[j, "r2_cd"] <- round(MuMIn::r.squaredGLMM(mixmod_cd)[1], 3)
    
    cd_estimate_univ[j, "estimate_cd_norm"] <- round(summary(mixmod_cd_norm)$coefficients[2], 3)
    cd_estimate_univ[j, "se_cd_norm"] <- round(summary(mixmod_cd_norm)$coefficients[4], 3)
    cd_estimate_univ[j, "pvalue_cd_norm"] <- round(summary(mixmod_cd_norm)$coefficients[10], 3)
    cd_estimate_univ[j, "r2_cd_norm"] <- round(MuMIn::r.squaredGLMM(mixmod_cd_norm)[1], 3)
    
    cd_estimate_univ[j, "estimate_cd_abs"] <- round(summary(mixmod_cd_abs)$coefficients[2], 3)
    cd_estimate_univ[j, "se_cd_abs"] <- round(summary(mixmod_cd_abs)$coefficients[4], 3)
    cd_estimate_univ[j, "pvalue_cd_abs"] <- round(summary(mixmod_cd_abs)$coefficients[10], 3)
    cd_estimate_univ[j, "r2_cd_abs"] <- round(MuMIn::r.squaredGLMM(mixmod_cd_abs)[1], 3)
    
    cd_estimate_univ[j, "estimate_cd_abs2"] <- round(summary(mixmod_cd_abs2)$coefficients[2], 3)
    cd_estimate_univ[j, "se_cd_abs2"] <- round(summary(mixmod_cd_abs2)$coefficients[4], 3)
    cd_estimate_univ[j, "pvalue_cd_abs2"] <- round(summary(mixmod_cd_abs2)$coefficients[10], 3)
    cd_estimate_univ[j, "r2_cd_abs2"] <- round(MuMIn::r.squaredGLMM(mixmod_cd_abs2)[1], 3)
    
  }
  
  cd_mixed_univ <- cd_estimate_univ%>%
    mutate(signif_cic=case_when(pvalue_cic<=0.1&pvalue_cic>0.05~".",
                               pvalue_cic<=0.05~"*",
                               pvalue_cic>0.1~""),# should give same results as cd signif
           signif_cd=case_when(pvalue_cd<=0.1&pvalue_cd>0.05~".",
                               pvalue_cd<=0.05~"*",
                               pvalue_cd>0.1~""),
           signif_cd_norm=case_when(pvalue_cd_norm<=0.1&pvalue_cd_norm>0.05~".",
                                    pvalue_cd_norm<=0.05~"*",
                                    pvalue_cd_norm>0.1~""),
           signif_cd_abs=case_when(pvalue_cd_abs<=0.1&pvalue_cd_abs>0.05~".",
                                    pvalue_cd_abs<=0.05~"*",
                                    pvalue_cd_abs>0.1~""),
           signif_cd_abs2=case_when(pvalue_cd_abs2<=0.1&pvalue_cd_abs2>0.05~".",
                                   pvalue_cd_abs2<=0.05~"*",
                                   pvalue_cd_abs2>0.1~""))
  cd_mixed_univ

  # write.table(cd_estimate_univ, paste0("./statistical_analyses/output/estimate_univariate_",
  #             method[i],".csv"), sep=";", dec=".", col.names=T, row.names=F, quote=F)

  
  total_cd_mean <- total_cd%>%
    group_by(variable, nurse)%>%
    summarise(oc_mean=mean(oc),
              oc_sd=sd(oc),
              oc_median=median(oc),
              cd_mean=mean(cd),
              cic_median=median(cic),
              cic_se=se(cic),
              cic_mean=mean(cic),
              cd_median=median(cd),
              cd_se=se(cd),
              cd_abs_mean=mean(cd_abs),
              cd_abs_median=median(cd_abs),
              cd_abs_se=se(cd_abs),
              cd_abs2_mean=mean(abs(cd)),
              cd_abs2_median=median(abs(cd)),
              cd_abs2_se=se(abs(cd)),
              cd_norm_mean=mean(cd_norm),
              cd_norm_se=se(cd_norm),
              cd_norm_median=median(cd_norm))
  
  total_cd_mean <- total_cd_mean%>%
    mutate(var_name=case_when(
      variable=="axis1"~"Axis1",
      variable=="axis2"~"Axis2",
      variable=="bio01"~"Annual \nTemp",
      variable=="bio04"~"Temp \nseason",
      variable=="bio05"~"Max \nTemp WM",
      variable=="bio06"~"Min \nTemp CM",
      variable=="bio07"~"Temp \nAnn range",
      variable=="bio10"~"Mean \nTemp WQ",
      variable=="bio11"~"Mean \nTemp CQ",
      variable=="bio12"~"Annual \nPrec",
      variable=="bio13"~"Prec WM",
      variable=="bio14"~"Prec DM",
      variable=="bio15"~"Prec \nSeason",
      variable=="bio16"~"Prec WQ",
      variable=="bio17"~"Prec DQ",
      variable=="multiv"~"Multiv"
    ),
    var_name1=case_when(
      variable=="axis1"~"Axis1",
      variable=="axis2"~"Axis2",
      variable=="bio01"~"Annual Temp",
      variable=="bio04"~"Temp season",
      variable=="bio05"~"Max Temp WM",
      variable=="bio06"~"Min Temp CM",
      variable=="bio07"~"Temp Ann range",
      variable=="bio10"~"Mean Temp WQ",
      variable=="bio11"~"Mean Temp CQ",
      variable=="bio12"~"Annual Prec",
      variable=="bio13"~"Prec WeM",
      variable=="bio14"~"Prec DM",
      variable=="bio15"~"Prec Season",
      variable=="bio16"~"Prec WeQ",
      variable=="bio17"~"Prec DQ",
      variable=="multiv"~"Multiv"
    ))
  
  
  ## estimate differences t test 0
  
  univ_cdchange_facil <- total_cd%>%
    filter(nurse=="facil")%>%
    select(-nurse_code)%>%
    pivot_wider(names_from = nurse, 
                values_from=c(cic, cd, cd_norm, cd_abs))
  
  univ_cdchange_open<- total_cd%>%
    filter(nurse=="open")%>%
    select(-nurse_code)%>%
    pivot_wider(names_from = nurse, 
                values_from=c(cic, cd, cd_norm, cd_abs))
  
  univ_cdchange <- univ_cdchange_facil%>%
    left_join(univ_cdchange_open%>%
                select(plot, plot_code,
                       variable, cic_open,
                       cd_open,
                       cd_norm_open,
                       cd_abs_open), by=c("plot",
                                           "plot_code",
                                           "variable"))%>%
    mutate(cic_change=cic_facil-cic_open,
           cd_change_abs=cd_abs_facil - cd_abs_open,#incorrect
           cd_change_abs2=abs(cd_facil) - abs(cd_open),# is the disequilibrium higher in quantity between facil and open
           cd_change=cd_facil-cd_open,# is the disequlibrium different between facil and open ie. more or less negative or positive
           cd_change_norm=cd_norm_facil-cd_norm_open)%>%
    mutate(cd_change_combi=case_when(
      cd_change_abs2>=0~abs(cd_change)*1,# positive means cd facil is higher
      cd_change_abs2<0~abs(cd_change)*(-1)# negative means cd facil is lower
    ))# it combines the real magnitude, coming from cd_facil-cd_open with the real sign of difference ie. abs(cd_facil)-abs(cd_open)
  
  cd_estimate_change <- tibble()
  
  
  for(k in 1:length(biovariable)){# other way of view the same anova is check if difference between groups is different from 0
    
    diseq_bio_change <-  univ_cdchange%>%
      dplyr::filter(variable==biovariable[k])
    
    if(biovariable[k]!="multiv"){
      z0 <- t.test(diseq_bio_change$cic_change, mu=0)
    }else if(biovariable[k]=="multiv"){
        z0 <-NA_real_
      }
    
    z <- t.test(diseq_bio_change$cd_change_abs, mu=0)
    z1 <- t.test(diseq_bio_change$cd_change, mu=0)
    z2 <- t.test(diseq_bio_change$cd_change_norm, mu=0)#realmente no hace falta relativizar por el valor de oc si se hace independientemente para cada varaible
    z3 <- t.test(diseq_bio_change$cd_change_abs2, mu=0)
    z4 <- t.test(diseq_bio_change$cd_change_combi, mu=0)
    
    
    cd_estimate_change[k, "row_number"] <- k
    cd_estimate_change[k, "variable"] <- biovariable[k]
    
    if(biovariable[k]!="multiv"){
      
    cd_estimate_change[k, "estimate_cic"] <- z0$estimate
    cd_estimate_change[k, "se_cic"] <- z0$stderr
    cd_estimate_change[k, "pvalue_cic"] <- z0$p.value
    cd_estimate_change[k, "upper_ci_cic"] <- z0$conf.int[2]
    cd_estimate_change[k, "lower_ci_cic"] <- z0$conf.int[1]
   
    } else if( biovariable[k]=="multiv"){
      
    cd_estimate_change[k, "estimate_cic"] <- NA_real_
    cd_estimate_change[k, "se_cic"] <- NA_real_
    cd_estimate_change[k, "pvalue_cic"] <- 1# to simplify ulterior code but should be also NA
    cd_estimate_change[k, "upper_ci_cic"] <- NA_real_
    cd_estimate_change[k, "lower_ci_cic"] <- NA_real_
      
    }
    
    cd_estimate_change[k, "estimate_cd_abs"] <- z$estimate
    cd_estimate_change[k, "se_cd_abs"] <- z$stderr
    cd_estimate_change[k, "pvalue_cd_abs"] <- z$p.value
    cd_estimate_change[k, "upper_ci_cd_abs"] <- z$conf.int[2]
    cd_estimate_change[k, "lower_ci_cd_abs"] <- z$conf.int[1]
    
    cd_estimate_change[k, "estimate_cd_abs2"] <- z3$estimate
    cd_estimate_change[k, "se_cd_abs2"] <- z3$stderr
    cd_estimate_change[k, "pvalue_cd_abs2"] <- z3$p.value
    cd_estimate_change[k, "upper_ci_cd_abs2"] <- z3$conf.int[2]
    cd_estimate_change[k, "lower_ci_cd_abs2"] <- z3$conf.int[1]
    
    cd_estimate_change[k, "estimate_cd"] <- z1$estimate
    cd_estimate_change[k, "se_cd"] <- z1$stderr
    cd_estimate_change[k, "pvalue_cd"] <- z1$p.value
    cd_estimate_change[k, "upper_ci_cd"] <- z1$conf.int[2]
    cd_estimate_change[k, "lower_ci_cd"] <- z1$conf.int[1]
    
    cd_estimate_change[k, "estimate_cd_norm"] <- z2$estimate
    cd_estimate_change[k, "se_cd_norm"] <- z2$stderr
    cd_estimate_change[k, "pvalue_cd_norm"] <- z2$p.value
    cd_estimate_change[k, "upper_ci_cd_norm"] <- z2$conf.int[2]
    cd_estimate_change[k, "lower_ci_cd_norm"] <- z2$conf.int[1]
    
    cd_estimate_change[k, "estimate_cd_combi"] <- z4$estimate
    cd_estimate_change[k, "se_cd_combi"] <- z4$stderr
    cd_estimate_change[k, "pvalue_cd_combi"] <- z4$p.value
    cd_estimate_change[k, "upper_ci_cd_combi"] <- z4$conf.int[2]
    cd_estimate_change[k, "lower_ci_cd_combi"] <- z4$conf.int[1]
    
    cd_estimate_change[k, "var_name"] <- unique(diseq_bio_change$var_name1)
    
  }
  
  #cd_estimate_change same as cd_estimate_univ
  #**cd_norm realmente solo es interesante para plotear, y ni siquiera eso porque lo mejor es liberar ejes**
  
  signi_df <- univ_cdchange %>%
    mutate(cd_norm_change=cd_norm_facil-cd_norm_open)%>%
    group_by(variable)%>%
    summarise(mean_cic_change=mean(cic_change),
              median_cic_change=median(cic_change),
              se_cic_change=se(cic_change),
              mean_cd_change=mean(cd_change),
              median_cd_change=median(cd_change),
              se_cd_change=se(cd_change),
              mean_cd_norm_change=mean(cd_norm_change),
              median_cd_norm_change=median(cd_norm_change),
              se_cd_norm_change=se(cd_norm_change),
              mean_cd_change_abs=mean(cd_change_abs),
              median_cd_change_abs=median(cd_change_abs),
              se_cd_change_abs=se(cd_change_abs),
              mean_cd_change_abs2=mean(cd_change_abs2),
              median_cd_change_abs2=median(cd_change_abs2),
              se_cd_change_abs2=se(cd_change_abs2),
              mean_cd_change_combi=mean(cd_change_combi),
              median_cd_change_combi=median(cd_change_combi),
              se_cd_change_combi=se(cd_change_combi))%>%
    left_join(cd_estimate_univ%>%
                rename(pvalue_cic_mix= pvalue_cic,
                       pvalue_cd_mix=pvalue_cd,
                       pvalue_cd_norm_mix=pvalue_cd_norm,
                       pvalue_cd_abs_mix=pvalue_cd_abs,
                       pvalue_cd_abs2_mix=pvalue_cd_abs2),
              by="variable")%>%
    left_join(cd_estimate_change%>%
                select(variable, var_name,
                       upper_ci_cic, lower_ci_cic,
                       upper_ci_cd, lower_ci_cd,
                       upper_ci_cd_norm, lower_ci_cd_norm,
                       upper_ci_cd_abs, lower_ci_cd_abs,
                       upper_ci_cd_abs2, lower_ci_cd_abs2,
                       upper_ci_cd_combi, lower_ci_cd_combi,
                       pvalue_cic,
                       pvalue_cd, 
                       pvalue_cd_abs, 
                       pvalue_cd_abs2, 
                       pvalue_cd_norm, 
                       pvalue_cd_combi),
              by="variable")# estimate should change the sign
  
  signi_sum <- signi_df%>%
    select(variable, var_name,
           pvalue_cic,
           pvalue_cic_mix,#should be the same as pvalue cd
           pvalue_cd,
           pvalue_cd_mix,
          # pvalue_paired_cd,
           pvalue_cd_abs,
           pvalue_cd_abs_mix,
           pvalue_cd_abs2,
           pvalue_cd_abs2_mix,
           pvalue_cd_norm,
           pvalue_cd_norm_mix,
           pvalue_cd_combi)# mixmodels and t test actually give the same results
  
  almost_signif <- signi_sum%>%
    filter(pvalue_cd>0.05&# also pvalue_cd or pvalue_cd_abs, cd_norm no hace falta pq no es necesario normalizar
             pvalue_cd<0.1)%>%
    select(variable,
           var_name,
           pvalue_cd)
  
  signif <- signi_sum%>%
    filter(pvalue_cd<0.05)%>%
    select(variable,
           var_name,
           pvalue_cd)
  
  
  (uni_bar <- total_cd_mean%>%
    ggplot(aes(x=var_name, y=cd_mean,fill=nurse)) +
    geom_bar(position = "dodge", stat="identity")+
    geom_linerange(aes(ymin = cd_mean-cd_se,
                       ymax = cd_mean+cd_se),
                   position=position_dodge(width=0.9),
                   lwd=0.2)+
    viridis::scale_fill_viridis(discrete = TRUE,option = 'E')+
    #coord_flip()+ to make horizontal plot
    #geom_text(data=signi, mapping=aes(x=Inf,y=Inf,label=signi),hjust = 5, vjust = 1.5)+
    # annotate("text", x =almost_significant_var$row_number, y = c(-2,-0.3),
    #          label = ".", size=4)+
    # labs(fill = "Community")+
    ylab("Climatic Debt")+
    # xlab("")+
    #facet_wrap(~variable, scales="free_y",nrow=2)+
    xlab(" ")+
    theme_classic()+
    theme(#legend.position = c(0.12, 0.85),
      legend.title = element_text(size = 5.7),
      legend.text  = element_text(size = 5),
      legend.key.size = unit(2, "mm"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.title=element_text(color="black", size=7),
      axis.text.x =element_text(color="black", size=4),
      axis.text.y =element_text(color="black", size=5),
      text=element_text(family="serif"))
  )
  
  
  total_cd$var_name2 <- factor(total_cd$var_name1, 
                               levels = c("Annual Temp","Temp season", "Max Temp WM",
                                          "Min Temp CM","Temp Ann range","Mean Temp WQ",
                                          "Mean Temp CQ","Annual Prec","Prec WeM",
                                          "Prec DM","Prec Season","Prec WeQ", "Prec DQ","Multiv",
                                          "Axis1", "Axis2"),
                               labels = c("Annual Temp","Temp season", "Max Temp WM",
                                          "Min Temp CM","Temp Ann range","Mean Temp WQ",
                                          "Mean Temp CQ","Annual Prec","Prec WeM",
                                          "Prec DM","Prec Season","Prec WeQ", "Prec DQ","Multiv",
                                          "Axis1", "Axis2"))
  
  almost_significant_var <- cd_estimate_univ%>%
    dplyr::filter(pvalue_cd>0.05&
                    pvalue_cd<=0.1)
  
  significant_var <- cd_estimate_univ%>%
    dplyr::filter(pvalue_cd<=0.05)
  
  
  sig_labels <- data.frame(variable = significant_var$variable, 
                           label = rep("*", length(significant_var$variable)))
  a <- sig_labels 
  a$nurse <- "open"
  b <- sig_labels
  b$nurse <- "facil"
  sig_labels <- rbind(a, b)
  
  sig_labels <- sig_labels%>%
    left_join(total_cd%>%
                select(variable, var_name,
                       var_name1,
                       var_name2)%>%
                distinct(), by="variable")
  
  almsig_labels <- data.frame(variable = almost_significant_var$variable, 
                              label = rep(".", length(almost_significant_var$variable)))
  
  c <- almsig_labels 
  c$nurse <- "open"
  d <- almsig_labels
  d$nurse <- "facil"
  almsig_labels <- rbind(c, d)
  
  almsig_labels <- almsig_labels%>%
    left_join(total_cd%>%
                select(variable, var_name,
                       var_name1, var_name2)%>%
                distinct(), by="variable")
  
  y_signif <- total_cd%>%
    filter(variable%in%sig_labels$variable)%>%
    group_by(variable)%>%
    summarise(max_cd=max(cd),
              min_cd=min(cd))%>%
    mutate(cd_range=max_cd-min_cd,
           y=max_cd+0.15*cd_range)
  
  sig_labels <- sig_labels%>%
    right_join(y_signif, by="variable")
  
  y_almost <- total_cd%>%
    filter(variable%in%almsig_labels$variable)%>%
    group_by(variable)%>%
    summarise(max_cd=max(cd),
              min_cd=min(cd))%>%
    mutate(cd_range=max_cd-min_cd,
           y=max_cd+0.2*cd_range)
  
  almsig_labels <- almsig_labels%>%
    right_join(y_almost, by="variable")
  

  (gg_cd_boxpoint <- ggplot(data=total_cd,#%>%
                              # filter(variable%in%c("multiv",
                              #                      "axis1", "axis2",
                              #                      "bio01", "bio06",
                              #                      "bio12")),
                            aes(x=nurse, y=cd,
                                fill=nurse,color=nurse))+ 
      geom_violin(position = position_dodge(0.8),
                  color="transparent",alpha=0.4,width=1)+
      geom_boxplot(position =  position_dodge(0.8),
                   fill='transparent',width=0.4, lwd=0.4,
                   outlier.size = 0.3)+
      scale_color_manual(values=c("#E69F00", "#56B4E9"))+
      scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
      geom_jitter(width=0.15, size=0.3) +
      geom_hline(aes(yintercept = 0), linetype = "dashed", lwd=0.4) +
      geom_text(data=sig_labels,
                aes(x = 1.5, y = y,
                    label=label),
                col="black")+
      geom_text(data = almsig_labels,
                aes(x=1.5,
                    y=y,
                    label=label),
                col="black", size=6)+
      xlab(" ")+
      ylab("Climatic Disequilibrium")+
      facet_wrap(~var_name2, scales="free_y",
                 nrow=4)+#labeller=labeller(variable='var_name1')
      # facet_wrap(~variable, scales="free_y",
      #            nrow=2)+
      scale_x_discrete(labels=c("under \ncanopies", "open"))+
      theme_linedraw()+
      theme(strip.background =element_rect(fill="grey80"),
            strip.text = element_text(colour = 'black', size=6),
            legend.position = "none",
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(color="black", size=5),
            axis.title=element_text(color="black", size=7),
            text=element_text(family="serif"))
  )
  
  # 
  # ggsave(paste0("./statistical_analyses/output/boxplot nurse/univariate_cd_",
  #               method[i], "_with0_boxplot_names.png"),
  #        plot=gg_cd_boxpoint, width = 12, height = 16, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/boxplot nurse/univariate_cd_",
  #               method[i], "_with0_boxplot_names.pdf"),
  #        plot=gg_cd_boxpoint, width = 12, height = 16, units="cm")
  # 
  

  almost_significant_var <- cd_estimate_univ%>%
    dplyr::filter(pvalue_cic>0.05&
                    pvalue_cic<=0.1)
  
  significant_var <- cd_estimate_univ%>%
    dplyr::filter(pvalue_cic<=0.05)
  
  
  sig_labels <- data.frame(variable = significant_var$variable, 
                           label = rep("*", length(significant_var$variable)))
  a <- sig_labels 
  a$nurse <- "open"
  b <- sig_labels
  b$nurse <- "facil"
  sig_labels <- rbind(a, b)
  
  sig_labels <- sig_labels%>%
    left_join(total_cd%>%
                select(variable, var_name,
                       var_name1,
                       var_name2)%>%
                distinct(), by="variable")
  
  almsig_labels <- data.frame(variable = almost_significant_var$variable, 
                              label = rep(".", length(almost_significant_var$variable)))
  
  c <- almsig_labels 
  c$nurse <- "open"
  d <- almsig_labels
  d$nurse <- "facil"
  almsig_labels <- rbind(c, d)
  
  almsig_labels <- almsig_labels%>%
    left_join(total_cd%>%
                select(variable, var_name,
                       var_name1, var_name2)%>%
                distinct(), by="variable")
  
  y_signif <- total_cd%>%
    filter(variable%in%sig_labels$variable)%>%
    group_by(variable)%>%
    summarise(max_cic=max(cic),
              min_cic=min(cic))%>%
    mutate(cic_range=max_cic-min_cic,
           y=max_cic+0.1*cic_range)
  
  sig_labels <- sig_labels%>%
    right_join(y_signif, by="variable")
  
  y_almost <- total_cd%>%
    filter(variable%in%almsig_labels$variable)%>%
    group_by(variable)%>%
    summarise(max_cic=max(cic),
              min_cic=min(cic))%>%
    mutate(cic_range=max_cic-min_cic,
           y=max_cic+0.2*cic_range)
  
  almsig_labels <- almsig_labels%>%
    right_join(y_almost, by="variable")
  
  
  
  (gg_cic_boxpoint <- ggplot(data=total_cd%>%
                               filter(variable!="multiv"),
                             aes(x=nurse, y=cic,
                                fill=nurse,color=nurse))+ 
      geom_violin(position = position_dodge(0.8),
                  color="transparent",alpha=0.4,width=1)+
      geom_boxplot(position =  position_dodge(0.8),
                   fill='transparent',width=0.4, lwd=0.4,
                   outlier.size = 0.3)+
      scale_color_manual(values=c("#E69F00", "#56B4E9"))+
      scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
      geom_jitter(width=0.15, size=0.3) +
      #geom_hline(aes(yintercept = 0), linetype = "dashed", lwd=0.4) +
      geom_text(data=sig_labels,
                aes(x = 1.5, y = y,
                    label=label),
                col="black")+
      geom_text(data = almsig_labels,
                aes(x=1.5,
                    y=y,
                    label=label),
                col="black", size=6)+
      xlab(" ")+
      ylab("Community Inferred Climate")+
      facet_wrap(~var_name2, scales="free_y",
                 nrow=3)+#labeller=labeller(variable='var_name1')
      # facet_wrap(~variable, scales="free_y",
      #            nrow=2)+
      scale_x_discrete(labels=c("under \ncanopies", "open"))+
      theme_linedraw()+
      theme(strip.background =element_rect(fill="grey80"),
            strip.text = element_text(colour = 'black', size=6),
            legend.position = "none",
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(color="black", size=5),
            axis.title=element_text(color="black", size=7),
            text=element_text(family="serif"))
  )
  
  # 
  # ggsave(paste0("./statistical_analyses/output/boxplot nurse/univariate_cic_",
  #               method[i], ".png"),
  #        plot=gg_cic_boxpoint, width = 16, height = 13, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/boxplot nurse/univariate_cic_",
  #               method[i], ".pdf"),
  #        plot=gg_cic_boxpoint, width = 16, height = 13, units="cm")
  # 
  # 
  # 
  # 
  

# Figure 4 ----------------------------------------------------------------
# i=1


  signi_df$var_name <- factor(signi_df$var_name, 
                               levels = c("Annual Temp","Temp season", "Max Temp WM",
                                          "Min Temp CM","Temp Ann range","Mean Temp WQ",
                                          "Mean Temp CQ","Annual Prec","Prec WeM",
                                          "Prec DM","Prec Season","Prec WeQ", "Prec DQ","Multiv",
                                          "Axis1", "Axis2"))
  

  # linerrange plot comparing differences in cd_abs2== magnitude of dis

  
  (mis_linerrange <- signi_df %>%
      mutate(estimate_cd_s=case_when(
        variable=="bio04"~-estimate_cd_abs2/10,
        variable!="bio04"~-estimate_cd_abs2
      ),
      upper_ci_cd_s= case_when(
        variable=="bio04"~ upper_ci_cd_abs2/10,
        variable!="bio04"~ upper_ci_cd_abs2
      ),
      lower_ci_cd_s = case_when(
        variable=="bio04"~ lower_ci_cd_abs2/10,
        variable!="bio04"~ lower_ci_cd_abs2
      ))%>%
      filter(variable%in%c("multiv",
                           "axis1", "axis2",
                           "bio01", "bio06",
                           "bio12"))%>%
      ggplot(aes(y=estimate_cd_s,x=var_name, 
                 fill = pvalue_cd_abs2<=0.05))+
      coord_flip()+
      geom_linerange(aes(ymax = upper_ci_cd_s, 
                         ymin = lower_ci_cd_s),
                     size=0.3,
                     show.legend = FALSE)+
      geom_point(size=1.4,
                 shape=21,
                 stroke=0.3,
                 show.legend = FALSE)+
      geom_hline(yintercept = 0, linetype=3)+
      scale_fill_manual(name = 'P value < 0.05', 
                        values = setNames(c('white','black'),c(F, T))) +
      ylab("Absolute Climatic Disequilibrium")+
      #xlim(-0.1,1.15)+
      theme_classic()+
      theme(legend.position = "none",
            axis.title.y= element_blank(),
            axis.title.x= element_text(color="black",
                                       family = "serif",
                                       size=9),
            axis.text=element_text(color="black",
                                   family = "serif", size=8),
            text=element_text(family="serif"))
  )
  
  # ggsave(paste0("./statistical_analyses/output/change (facil_open)/univariate_",
  #               method[i], "_linerrange_mismatch_quantity.png"),
  #        plot=mis_linerrange, width = 8, height = 14, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/change (facil_open)/univariate_",
  #               method[i], "_linerrange_mismatch_quantity.pdf"),
  #        plot=mis_linerrange, width = 8, height = 14, units="cm")
  # 


  
  
}




