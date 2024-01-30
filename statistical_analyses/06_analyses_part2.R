

#******************************************#
## MARIA ANGELES PEREZ
## CREAF 10/2023
## statistical analyses part 2 
##*****************************************#


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


## 0. Load and prepare tables ####
#**************************************************************************************

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

# load extra explanatory variables

expl <-  read_delim("./statistical_analyses/output/explanatory_vars.csv", delim=";", col_names=T)

# join tables

method <- unique(diseq_uni$type)
i <- 1

## cd_abs is estimated but it is inacuarate as it is estimated as abs(cic)-abs(oc) and this make no sense

multiv_total <- diseq_multi%>%
  dplyr::filter(type==method[i])

multiv_cd0 <- multiv_total%>%
  select(-c(oc_x, oc_y,
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
  select(-c(oc_x, oc_y,
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
  select(-c(cic_x, cic_y,
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

expl <- expl%>%
  select(-c("plot",
            "plot_code",
            "nurse",
            "oc_bio01",
            "oc_bio04",
            "oc_bio05", 
            "oc_bio06",
            "oc_bio12",
            "community"))%>%
  mutate(ai=1/annual_PPET_sum)

univ_total1 <- univ_total%>%
  left_join(expl, by="nurse_code")

nrow(univ_total)  
nrow(univ_total1)

univ_total <- univ_total1


## 1. check differences between CD vs OC trends and uniform CIC and CD expectation ####
#**************************************************************************************

##*** 1.1 temperature model ***####

names(univ_total)
unique(univ_total$type)

mixmod_temp1 <- lmer(cd_bio01~nurse*oc_bio01+(1|plot), data=univ_total)
summary(mixmod_temp1)
MuMIn::r.squaredGLMM(mixmod_temp1)

temp1_slopes <- emmeans::emtrends(mixmod_temp1, "nurse", var = "oc_bio01") %>% 
  emmeans::test()%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio01",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)


temp1_pairs <- emmeans::emtrends(mixmod_temp1, "nurse", var = "oc_bio01") %>% 
  pairs()%>%
  as_tibble() %>% 
  rename(trend_dif = 2,# add also facil open column in case they are reverse for some continuous variables
         pvalue_dif=6) %>%
  mutate(variable = "oc_bio01")%>%
  dplyr::select(variable,
                trend_dif,
                pvalue_dif)


temp1_slopes <- temp1_slopes%>%
  mutate(ef_size=effectsize::t_to_d(t.ratio, 
                                    nrow(mixmod_temp1@frame))$d,#use nrow(mod@frame)instead nrow(data_table)in case na are authomatically removed in the model
         ef_ci_lower=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_temp1@frame))$CI_low,
         ef_ci_upper=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_temp1@frame))$CI_high,
         r2m=MuMIn::r.squaredGLMM(mixmod_temp1)[1],
         r2c=MuMIn::r.squaredGLMM(mixmod_temp1)[2])

temp1_df<- temp1_slopes%>%
  left_join(temp1_pairs, by="variable")%>%
  mutate(model="annual temp")



ef <- effects::Effect(c("nurse","oc_bio01"), mixmod_temp1, 25)# only fixed effects
eff_df<- data.frame(ef)

(nurse_temp1 <-  ggplot()+
    geom_point(data=univ_total, aes(x=oc_bio01, y=cd_bio01, 
                                    color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df, aes(x=oc_bio01,
                               y= fit,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower, ymax=upper, 
                                 x=oc_bio01, fill = nurse), alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    # ggrepel::geom_label_repel(data=univ_total,
    #                           aes(x=oc_bio01, y=cd_bio01,
    #                               label = plot_code),
    #                           box.padding   = 0.15,
    #                           point.padding = 0.15,
    #                           label.padding = 0.13,
    #                           segment.color = 'grey50',
    #                           max.overlaps = 50, size=1.3)+
    ylab("Climatic Disequilibrium \n Mean Annual Temperature")+
    xlab("Observed Mean Annual Temperature")+
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


# ggsave(paste0("../../output/diseq change as response/diseq_ave_temp", method[i], ".png"),
#        plot=nurse_temp1, width = 8, height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("../../output/diseq change as response/diseq_ave_temp", method[i], ".pdf"),
#        plot=nurse_temp1, width = 8, height = 6.3, units="cm")
# 


##*** 1.2 temperature model with uniform random CIC ***####

univ_random <- univ_total%>%
  dplyr::select(plot, plot_code,
                nurse, nurse_code,
                oc_bio01, cic_bio01,
                cd_bio01,type,
                lat, long,
                canopy_percent, gap_percent,
                site)%>%
  mutate(cic_bio01_random=mean(univ_total$cic_bio01),
         cd_bio01_random=cic_bio01_random-oc_bio01)

plot(univ_random$cd_bio01_random~univ_random$oc_bio01)

names(univ_random)
unique(univ_random$type)

lm_temp1 <- lm(cd_bio01_random~oc_bio01, data=univ_random)# not required to differentiate by nurse effect
summary(lm_temp1)

ef_random <- effects::Effect(c("oc_bio01"), lm_temp1, 25)# only fixed effects
eff_df_random<- data.frame(ef_random)

(random_temp1 <-  ggplot()+
    geom_point(data=univ_random, aes(x=oc_bio01, y=cd_bio01_random, 
                                     color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df_random, aes(x=oc_bio01,
                                      y= fit),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df_random, aes(ymin=lower, ymax=upper, 
                                        x=oc_bio01), alpha = 0.2,
                show.legend=F)+
    ylab("Climatic Disequilibrium \n Mean Annual Temperature")+
    xlab("Observed Mean Annual Temperature")+
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


##*** 1.3 test differences in temp for random slope -1 and add plot ***####

# Figure 6a ---------------------------------------------------------------

# remember method used to estimate niche centroid = perc95

emmeans::emtrends(mixmod_temp1, "nurse", var = "oc_bio01")%>%
  emmeans::test()

mod_slope1 <- emmeans::emtrends(mixmod_temp1, ~-1,var = "oc_bio01")%>%
  emmeans::test()

temp1_slopes <- emmeans::emtrends(mixmod_temp1, "nurse", var = "oc_bio01") %>% 
  emmeans::test(null=-1)%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio01",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)


(random_temp1 <-  ggplot()+
    geom_point(data=univ_total, aes(x=oc_bio01/10, 
                                    y=cd_bio01/10, 
                                    color=nurse),
               size=1, shape=19)+
    geom_hline(yintercept=0, linetype=2,
               color="grey40")+
    geom_line(data=eff_df_random, aes(x=oc_bio01/10,
                                      y= fit/10),
              color="grey40", linetype=2,
              size=0.5,alpha=1)+
    geom_line(data=eff_df, aes(x=oc_bio01/10,
                               y= fit/10,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower/10, ymax=upper/10, 
                                 x=oc_bio01/10, fill = nurse), 
                alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    # geom_abline(slope=1,
    #             intercept=mean(univ_total$oc_bio01),
    #             linetype="dashed")+
    xlim(c(6.5, 18))+
    labs(color = "Community")+
    # ggrepel::geom_label_repel(data=univ_total,
    #                           aes(x=oc_bio01, y=cd_bio01,
    #                               label = plot_code),
    #                           box.padding   = 0.15,
    #                           point.padding = 0.15,
    #                           label.padding = 0.13,
    #                           segment.color = 'grey50',
    #                           max.overlaps = 50, size=1.3)+
    ylab("Climatic Disequilibrium \n mean annual temperature")+
    xlab("Observed mean annual temperature")+
    #ggtitle("a)")+
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


# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_ave_temp",
#               method[i], "null_slope.png"),
#        plot=random_temp1, width = 8, 
#        height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_ave_temp",
#               method[i], "null_slope.pdf"),
#        plot=random_temp1, width = 8,
#        height = 6.3, units="cm")


##*** 1.4 temperature of coldest month model ***####


mixmod_temp6 <- lmer(cd_bio06~nurse*oc_bio06+ (1|plot), data=univ_total)
summary(mixmod_temp6)
MuMIn::r.squaredGLMM(mixmod_temp6)

temp6_slopes <- emmeans::emtrends(mixmod_temp6, "nurse", var = "oc_bio06") %>% 
  emmeans::test()%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio06",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)


temp6_pairs <- emmeans::emtrends(mixmod_temp6, "nurse", var = "oc_bio06") %>% 
  pairs()%>%
  as_tibble() %>% 
  rename(trend_dif = 2,# add also facil open column in case they are reverse for some continuous variables
         pvalue_dif=6) %>%
  mutate(variable = "oc_bio06")%>%
  dplyr::select(variable,
                trend_dif,
                pvalue_dif)


temp6_slopes <- temp6_slopes%>%
  mutate(ef_size=effectsize::t_to_d(t.ratio, 
                                    nrow(mixmod_temp6@frame))$d,#use nrow(mod@frame)instead nrow(data_table)in case na are authomatically removed in the model
         ef_ci_lower=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_temp6@frame))$CI_low,
         ef_ci_upper=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_temp6@frame))$CI_high,
         r2m=MuMIn::r.squaredGLMM(mixmod_temp6)[1],
         r2c=MuMIn::r.squaredGLMM(mixmod_temp6)[2])

temp6_df<- temp6_slopes%>%
  left_join(temp6_pairs, by="variable")%>%
  mutate(model="minimum temp")


ef <- effects::Effect(c("nurse","oc_bio06"), mixmod_temp6, 25)# only fixed effects
eff_df<- data.frame(ef)

(nurse_temp6 <-  ggplot()+
    geom_point(data=univ_total, aes(x=oc_bio06, y=cd_bio06, 
                                    color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df, aes(x=oc_bio06,
                               y= fit,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower, ymax=upper, 
                                 x=oc_bio06, fill = nurse), alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    # ggrepel::geom_label_repel(data=univ_total,
    #                           aes(x=oc_bio06, y=cd_bio06,
    #                               label = plot_code),
    #                           box.padding   = 0.15,
    #                           point.padding = 0.15,
    #                           label.padding = 0.13,
    #                           segment.color = 'grey50',
    #                           max.overlaps = 50, size=1.3)+
    ylab("Climatic disequilibrium\n miniumum temperature")+
    xlab("Observed minimum temperature")+
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


# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_min_temp", method[i], ".png"),
#        plot=nurse_temp6, width = 8, height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_min_temp", method[i], ".pdf"),
#        plot=nurse_temp6, width = 8, height = 6.3, units="cm")


##*** 1.5 temperature of coldest month model with uniform random CIC ***####

univ_random <- univ_total%>%
  dplyr::select(plot, plot_code,
                nurse, nurse_code,
                oc_bio06, cic_bio06,
                cd_bio06,type,
                lat, long,
                canopy_percent, gap_percent,
                site)%>%
  mutate(cic_bio06_random=mean(univ_total$cic_bio06),
         cd_bio06_random=cic_bio06_random-oc_bio06)

plot(univ_random$cd_bio06_random~univ_random$oc_bio06)

names(univ_random)
unique(univ_random$type)

lm_temp6 <- lm(cd_bio06_random~oc_bio06, data=univ_random)# not required to differentiate by nurse effect
summary(lm_temp6)

ef_random <- effects::Effect(c("oc_bio06"), lm_temp6, 25)# only fixed effects
eff_df_random<- data.frame(ef_random)

(random_temp <-  ggplot()+
    geom_point(data=univ_random, aes(x=oc_bio06, y=cd_bio06_random, 
                                     color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df_random, aes(x=oc_bio06,
                                      y= fit),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df_random, aes(ymin=lower, ymax=upper, 
                                        x=oc_bio06), alpha = 0.2,
                show.legend=F)+
    ylab("Climatic disequilibrium\n miniumum temperature")+
    xlab("Observed minimum temperature")+
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


##*** 1.6 test differences in temp of coldest month for random slope -1 and add plot ***####



emmeans::emtrends(mixmod_temp6, "nurse", var = "oc_bio06")%>%
  emmeans::test()

temp6_slopes <- emmeans::emtrends(mixmod_temp6, "nurse", var = "oc_bio06") %>% 
  emmeans::test(null=0)%>%#null=-1 #in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio06",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)


(random_temp6 <-  ggplot()+
    geom_point(data=univ_total, aes(x=oc_bio06/10, y=cd_bio06/10, 
                                    color=nurse),
               size=1, shape=19)+
    geom_hline(yintercept=0, linetype=2,
               color="grey40")+
    geom_line(data=eff_df_random, aes(x=oc_bio06/10,
                                      y= fit/10),
              color="grey40", linetype=2,
              size=0.5,alpha=1)+
    geom_line(data=eff_df, aes(x=oc_bio06/10,
                               y= fit/10,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower/10, ymax=upper/10, 
                                 x=oc_bio06/10, fill = nurse), alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    # ggrepel::geom_label_repel(data=univ_total,
    #                           aes(x=oc_bio01, y=cd_bio01,
    #                               label = plot_code),
    #                           box.padding   = 0.15,
    #                           point.padding = 0.15,
    #                           label.padding = 0.13,
    #                           segment.color = 'grey50',
    #                           max.overlaps = 50, size=1.3)+
    ylab("Climatic Disequilibrium\n miniumum temperature")+
    xlab("Observed minimum temperature")+
    #ggtitle("b)")+
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



ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_min_temp",
              method[i], "null_slope.png"),
       plot=random_temp6, width = 8, 
       height = 6.3, dpi = 300, units="cm")
ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_min_temp",
              method[i], "null_slope.pdf"),
       plot=random_temp6, width = 8,
       height = 6.3, units="cm")



##*** 1.7 precip vs Aridity index model ***####


mixmod_prec <- lmer(cd_bio12~nurse*ai+ (1|plot), data=univ_total)
summary(mixmod_prec)
MuMIn::r.squaredGLMM(mixmod_prec)


prec_slopes <- emmeans::emtrends(mixmod_prec, "nurse", var = "ai") %>% 
  emmeans::test()%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio12",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)


prec_pairs <- emmeans::emtrends(mixmod_prec, "nurse", var = "ai") %>% 
  pairs()%>%
  as_tibble() %>% 
  rename(trend_dif = 2,# add also facil open column in case they are reverse for some continuous variables
         pvalue_dif=6) %>%
  mutate(variable = "oc_bio12")%>%
  dplyr::select(variable,
                trend_dif,
                pvalue_dif)


prec_slopes <- prec_slopes%>%
  mutate(ef_size=effectsize::t_to_d(t.ratio, 
                                    nrow(mixmod_prec@frame))$d,#use nrow(mod@frame)instead nrow(data_table)in case na are authomatically removed in the model
         ef_ci_lower=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_prec@frame))$CI_low,
         ef_ci_upper=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_prec@frame))$CI_high,
         r2m=MuMIn::r.squaredGLMM(mixmod_prec)[1],
         r2c=MuMIn::r.squaredGLMM(mixmod_prec)[2])

prec_df<- prec_slopes%>%
  left_join(prec_pairs, by="variable")%>%
  mutate(model="precipitation")

ef <- effects::Effect(c("nurse","ai"), mixmod_prec, 25)# only fixed effects
eff_df<- data.frame(ef)


(nurse_prec <-  ggplot()+
    geom_point(data=univ_total, aes(x=ai, y=cd_bio12, 
                                    color=nurse),
               size=1, shape=19)+
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
    # ggrepel::geom_label_repel(data=univ_total,
    #                           aes(x=ai, y=cd_bio12,
    #                               label = plot_code),
    #                           box.padding   = 0.15,
    #                           point.padding = 0.15,
    #                           label.padding = 0.13,
    #                           segment.color = 'grey50',
    #                           max.overlaps = 50, size=1.3)+
    ylab("Climatic disequilibrium\n annual precipitation")+
    xlab("Aridity index (PET/P)")+
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




# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_prec_ai", method[i], ".png"),
#        plot=nurse_prec, width = 8, height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_prec_ai", method[i], ".pdf"),
#        plot=nurse_prec, width = 8, height = 6.3, units="cm")

##*** 1.8 precip vs Aridity index with uniform random CIC ***####


univ_random <- univ_total%>%
  dplyr::select(plot, plot_code,
                nurse, nurse_code,
                oc_bio12, cic_bio12,
                cd_bio12,type,
                lat, long,
                canopy_percent, gap_percent,
                site, ai)%>%
  mutate(cic_bio12_random=mean(univ_total$cic_bio12),
         cd_bio12_random=cic_bio12_random-oc_bio12)

plot(univ_random$cd_bio12_random~univ_random$ai)
# this is not a perfect straight line as precipitation and PET/P do 
# not perfectly correlate as temperature with itself in the previous models

names(univ_random)
unique(univ_random$type)

lm_prec <- lm(cd_bio12_random~ai, data=univ_random)# not required to differentiate by nurse effect
summary(lm_prec)

ef_random <- effects::Effect(c("ai"), lm_prec, 25)# only fixed effects
eff_df_random<- data.frame(ef_random)

(random_prec <-  ggplot()+
    geom_point(data=univ_random, aes(x=ai, y=cd_bio12_random, 
                                     color=nurse),
               size=1, shape=19)+
    geom_line(data=eff_df_random, aes(x=ai,
                                      y= fit),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df_random, aes(ymin=lower, ymax=upper, 
                                        x=ai), alpha = 0.2,
                show.legend=F)+
    ylab("Climatic disequilibrium\n annual precipitation")+
    xlab("Aridity index (PET/P)")+
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




##*** 1.9 test differences in precip for random slope x and add plot ***####


# Figure 6b ---------------------------------------------------------------



emmeans::emtrends(mixmod_prec, "nurse", var = "ai")%>%
  emmeans::test()

prec_slopes <- emmeans::emtrends(mixmod_prec, "nurse", var = "ai") %>% 
  emmeans::test(null=304.75)%>%#null=0 #in this case the intercept for random CIC community is 304.75
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "ai",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)


(random_prec <-  ggplot()+
    geom_point(data=univ_total, aes(x=ai, y=cd_bio12, 
                                    color=nurse),
               size=1, shape=19)+
    geom_hline(yintercept=0, linetype=2,
               color="grey40")+
    geom_line(data=eff_df_random, aes(x=ai,
                                      y= fit),
              size=0.5,alpha=1, color="grey40", linetype="dashed")+
    geom_ribbon(data=eff_df_random, aes(ymin=lower, ymax=upper,
                                        x=ai), alpha = 0.1,
                fill="grey40",
                show.legend=F)+# consider not to run
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
    # ggrepel::geom_label_repel(data=univ_total,
    #                           aes(x=ai, y=cd_bio12,
    #                               label = plot_code),
    #                           box.padding   = 0.15,
    #                           point.padding = 0.15,
    #                           label.padding = 0.13,
    #                           segment.color = 'grey50',
    #                           max.overlaps = 50, size=1.3)+
    ylab("Climatic Disequilibrium\n Annual Precipitation")+
    xlab("Aridity Index (PET/P)")+
    #ggtitle("c)")+
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


# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_prec",
#               method[i], "null_slope.png"),
#        plot=random_prec, width = 8, 
#        height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_prec",
#               method[i], "null_slope.pdf"),
#        plot=random_prec, width = 8,
#        height = 6.3, units="cm")



#*** 1.10 join all 3 panels ####


triplot <- gridExtra::grid.arrange(random_temp1, random_temp6,
                                   random_prec,
                                   ncol=2, widths = unit(c(8,8), "cm"),
                                   heights = unit(c(6.5, 7.4),'cm'))

ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_triplot",
              method[i], "_null_slope.png"),
       plot=triplot, width = 17, 
       height = 14.3, dpi = 300, units="cm")
ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_triplot",
              method[i], "_null_slope.pdf"),
       plot=triplot, width = 17,
       height = 14.3, units="cm")






biplot <- gridExtra::grid.arrange(random_temp1,
                                  random_prec,
                                  ncol=1, widths = unit(8, "cm"),
                                  heights = unit(c(6.5, 7.4),'cm'))

ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_biplot",
              method[i], "_null_slope_fig6.png"),
       plot=biplot, width = 17, 
       height = 14.3, dpi = 300, units="cm")
ggsave(paste0("./statistical_analyses/output/diseq change as response/diseq_biplot",
              method[i], "_null_slope_fig6.pdf"),
       plot=biplot, width = 17,
       height = 14.3, units="cm")



## 2. Check differences between CIC vs OC ####
#*********************************************

##*** 2.1 Average Temperature model ***####


mixmod_temp1 <- lmer(cic_bio01~nurse*oc_bio01+(1|plot), data=univ_total)
summary(mixmod_temp1)
MuMIn::r.squaredGLMM(mixmod_temp1)


emmeans::emtrends(mixmod_temp1, "nurse", var = "oc_bio01")%>%
  emmeans::test()

mod_slope1 <- emmeans::emtrends(mixmod_temp1, ~1,var = "oc_bio01")%>%
  emmeans::test()

temp1_slopes <- emmeans::emtrends(mixmod_temp1, "nurse", var = "oc_bio01") %>% 
  emmeans::test(null=1)%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio01",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)

ef <- effects::Effect(c("nurse","oc_bio01"), mixmod_temp1, 25)# only fixed effects
eff_df<- data.frame(ef)

(cic_temp1 <-  ggplot()+
    geom_point(data=univ_total, aes(x=oc_bio01/10, 
                                    y=cic_bio01/10, 
                                    color=nurse),
               size=1, shape=19)+
    geom_hline(yintercept=mean(univ_total$cic_bio01/10), 
               linetype=2,
               color="grey40")+
    geom_abline(intercept = 0, slope = 1,
                color="grey40", linetype=2)+
    geom_line(data=eff_df, aes(x=oc_bio01/10,
                               y= fit/10,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower/10, ymax=upper/10, 
                                 x=oc_bio01/10, fill = nurse), 
                alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    ylab("Community Inferred Climate \nmean annual temperature (ºC)")+
    xlab("Observed mean annual temperature")+
    #ggtitle("a)")+
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


# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_ave_temp",
#               method[i], "null_slope.png"),
#        plot=cic_temp1, width = 8, 
#        height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_ave_temp",
#               method[i], "null_slope.pdf"),
#        plot=cic_temp1, width = 8,
#        height = 6.3, units="cm")
# 
# 



##*** 2.2 Temperature of coldest month model ***####



mixmod_temp6 <- lmer(cic_bio06~nurse*oc_bio06+(1|plot), data=univ_total)
summary(mixmod_temp6)
MuMIn::r.squaredGLMM(mixmod_temp6)


emmeans::emtrends(mixmod_temp6, "nurse", var = "oc_bio06")%>%
  emmeans::test()

mod_slope1 <- emmeans::emtrends(mixmod_temp6, ~1,var = "oc_bio06")%>%
  emmeans::test()

temp1_slopes <- emmeans::emtrends(mixmod_temp6, "nurse", var = "oc_bio06") %>% 
  emmeans::test(null=1)%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio06",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)

ef <- effects::Effect(c("nurse","oc_bio06"), mixmod_temp6, 25)# only fixed effects
eff_df<- data.frame(ef)


(cic_temp6 <-  ggplot()+
    geom_point(data=univ_total, aes(x=oc_bio06/10, 
                                    y=cic_bio06/10, 
                                    color=nurse),
               size=1, shape=19)+
    geom_hline(yintercept=mean(univ_total$cic_bio06/10), 
               linetype=2,
               color="grey40")+
    geom_abline(intercept = 0, slope = 1,
                color="grey40", linetype=2)+
    geom_line(data=eff_df, aes(x=oc_bio06/10,
                               y= fit/10,color=nurse),
              size=0.5,alpha=1)+
    geom_ribbon(data=eff_df, aes(ymin=lower/10, ymax=upper/10, 
                                 x=oc_bio06/10, fill = nurse), 
                alpha = 0.2,
                show.legend=F)+
    scale_color_manual(values=c("#E69F00", "#56B4E9"), 
                       labels=c("Under canopy", "Open"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
    labs(color = "Community")+
    ylab("Community Inferred Climate \nminimum temperature (ºC)")+
    xlab("Observed minimum temperature")+
    #ggtitle("b)")+
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


# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_min_temp",
#               method[i], "null_slope.png"),
#        plot=cic_temp6, width = 8, 
#        height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_min_temp",
#               method[i], "null_slope.pdf"),
#        plot=cic_temp6, width = 8,
#        height = 6.3, units="cm")
# 


##*** 2.3 Precipitation model ***####

univ_random <- univ_total%>%
  dplyr::select(plot, plot_code,
                nurse, nurse_code,
                oc_bio12, cic_bio12,
                cd_bio12,type,
                lat, long,
                canopy_percent, gap_percent,
                site, ai)%>%
  mutate(cic_bio12_random=mean(univ_total$cic_bio12),
         cd_bio12_random=cic_bio12_random-oc_bio12)

names(univ_random)
unique(univ_random$type)

lm_prec <- lm(oc_bio12~ai, data=univ_random)# not required to differentiate by nurse effect
summary(lm_prec)

ef_random <- effects::Effect(c("ai"), lm_prec, 25)# only fixed effects
eff_df_random<- data.frame(ef_random)


mixmod_prec <- lmer(cic_bio12~nurse*ai+ (1|plot), data=univ_total)
summary(mixmod_prec)
MuMIn::r.squaredGLMM(mixmod_prec)


prec_slopes <- emmeans::emtrends(mixmod_prec, "nurse", var = "ai") %>% 
  emmeans::test(null=-304.75)%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
  as_tibble() %>% 
  rename(slope = 2) %>%
  mutate(variable = "oc_bio12",
         ci_upper= slope+1.96*SE,
         ci_lower= slope-1.96*SE)


prec_pairs <- emmeans::emtrends(mixmod_prec, "nurse", var = "ai") %>% 
  pairs()%>%
  as_tibble() %>% 
  rename(trend_dif = 2,# add also facil open column in case they are reverse for some continuous variables
         pvalue_dif=6) %>%
  mutate(variable = "oc_bio12")%>%
  dplyr::select(variable,
                trend_dif,
                pvalue_dif)


prec_slopes <- prec_slopes%>%
  mutate(ef_size=effectsize::t_to_d(t.ratio, 
                                    nrow(mixmod_prec@frame))$d,#use nrow(mod@frame)instead nrow(data_table)in case na are authomatically removed in the model
         ef_ci_lower=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_prec@frame))$CI_low,
         ef_ci_upper=effectsize::t_to_d(t.ratio, 
                                        nrow(mixmod_prec@frame))$CI_high,
         r2m=MuMIn::r.squaredGLMM(mixmod_prec)[1],
         r2c=MuMIn::r.squaredGLMM(mixmod_prec)[2])

prec_df<- prec_slopes%>%
  left_join(prec_pairs, by="variable")%>%
  mutate(model="precipitation")

ef <- effects::Effect(c("nurse","ai"), mixmod_prec, 25)# only fixed effects
eff_df<- data.frame(ef)


(cic_prec <-  ggplot()+
    geom_point(data=univ_total, aes(x=ai, y=cic_bio12, 
                                    color=nurse),
               size=1, shape=19)+
    geom_hline(yintercept=mean(univ_total$cic_bio12),
               linetype=2,
               color="grey40")+
    geom_line(data=eff_df_random, aes(x=ai,
                                      y= fit),
              size=0.5,alpha=1, color="grey40", linetype="dashed")+
    geom_ribbon(data=eff_df_random, aes(ymin=lower, ymax=upper,
                                        x=ai), alpha = 0.1,
                fill="grey40",
                show.legend=F)+# consider not to run
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
    ylab("Community Inferred Climate\nannual precipitation (mm)")+
    xlab("Aridity Index (PET/P)")+
    #ggtitle("c)")+
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


# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_prec",
#               method[i], "null_slope.png"),
#        plot=random_prec, width = 8, 
#        height = 6.3, dpi = 300, units="cm")
# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_prec",
#               method[i], "null_slope.pdf"),
#        plot=random_prec, width = 8,
#        height = 6.3, units="cm")




##*** 2.4 Join all 3 panels ***####

triplot <- gridExtra::grid.arrange(cic_temp1, cic_temp6,
                                   cic_prec,
                                   ncol=2, widths = unit(c(8,8), "cm"),
                                   heights = unit(c(6.5, 7.4),'cm'))
                                   
# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_triplot",
#               method[i], "_null_slope.png"),
#        plot=triplot, width = 17,
#        height = 14.3, dpi = 300, units="cm")
# 
# ggsave(paste0("./statistical_analyses/output/diseq change as response/cic_triplot",
#               method[i], "_null_slope.pdf"),
#        plot=triplot, width = 17,
#        height = 14.3, units="cm")        
# 


