



#*************************************#
### Maria Angeles Perez-Navarro
### CREAF  02/2021
### Community climatic diagrams univariate
#*************************************#


library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(cowplot)

## needed functions #

source('./niche_modelling/niche_analyses/niche_required_functions.R')

## create folder
dir.create("./statistical_analyses/output/")
dir.create("./statistical_analyses/output/climatic_diagrams/univariate/")


# 1. CIC, CD estimation with NICHE CENTROIDS -------------------------------

#* 1.1 prepare tables ------------------------------------------------------

table <-  read_delim("./niche_modelling/output/centroid_distances_univariate.csv", delim=";", col_names=T)
codes <- read_delim("./niche_modelling/output/perc95_distances_in0.csv", delim=";", col_names=T)
codes <- codes%>%
  dplyr::select(community, code)

table <- left_join(table, codes, by="community")%>%
  distinct()
head(table)
names(table)

table <- table%>%
  rename(plot= community,
         plot_code=code)%>%
  mutate(species_code= paste0(substring(genus, 1, 1),
                              substring(sp, 1,4)),
         nurse_code= paste0(plot_code, "_", nurse))

table_fq <- table %>%
  group_by(plot_code)%>%
  summarise(plot_n = sum(n))

table_fqq <- table%>%
  group_by(plot_code, nurse)%>%
  summarise(plot_nurse_n= sum(n))

table1 <- full_join(table, table_fq, by="plot_code")
table2 <- full_join(table1, table_fqq, by=c("plot_code", "nurse"))
table <- table2
rm(table1, table2)

table <- table%>%
  mutate(freq_plot = n/plot_n,
         freq_plot_nurse = n/plot_nurse_n)

#* 1.2 CIC, CD estimation---------------------------------------------------

# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

cic_df <-  table%>%
  dplyr::select(c(1:13, 71:77, 33:70))%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, names(table)[52:70])%>%
  distinct()%>%
  rename(oc_bio01= population_bio01, 
         oc_bio02= population_bio02,
         oc_bio03= population_bio03, 
         oc_bio04= population_bio04,
         oc_bio05= population_bio05, 
         oc_bio06= population_bio06,
         oc_bio07= population_bio07, 
         oc_bio08= population_bio08,
         oc_bio09= population_bio09, 
         oc_bio10= population_bio10,
         oc_bio11= population_bio11, 
         oc_bio12= population_bio12,
         oc_bio13= population_bio13, 
         oc_bio14= population_bio14,
         oc_bio15= population_bio15, 
         oc_bio16= population_bio16,
         oc_bio17= population_bio17, 
         oc_bio18= population_bio18,
         oc_bio19= population_bio19)

plot_list <- unique(cic_df$plot_code)
nurse_list <- unique(cic_df$nurse)


for (k in 1:19){
      
      biovariable <- if(k<=9){
        paste0("bio0", k )}else{
          paste0("bio", k )}
      
    for(i in 1:length(plot_list)) {
        
      table_plot <- table%>%
          filter(plot_code== plot_list[i])
        
        gg_list <- list()
        oc_list <- list()
        cic_list <- list()
        
        
      for(j in 1:length(nurse_list)){
      
      table_plot_n <- table_plot%>%
            filter(nurse== nurse_list[j])
          
      nurse_id <- unique(table_plot_n$nurse_code) 
      

      oc_bio <- cic_df[cic_df$nurse_code==nurse_id, paste0("oc_", biovariable)]
      cic_bio <- wt.mean(x=table_plot_n[, paste0("centroid_", biovariable)][[1]], # paste0("centroid_", biovariable)
                         wt=table_plot_n[,"freq_plot_nurse"][[1]])# freq_plot_nurse
      cic_df[cic_df$nurse_code==nurse_id, paste0("cic_", biovariable)] <- cic_bio
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_abs_", biovariable)] <- abs(cic_bio) - abs(oc_bio)
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_", biovariable)] <- cic_bio - oc_bio
      
      ggplot_x <- table_plot_n[, paste0("centroid_", biovariable)][[1]]
      x_range <- max(c(ggplot_x, oc_bio[[1]]))-min(c(ggplot_x, oc_bio[[1]]))
      gg_i <- ggplot(data= table_plot_n %>%
                       dplyr::mutate(freq_plot_nurse = case_when(
                         freq_plot_nurse == 0~NA_real_,
                         freq_plot_nurse != 0~freq_plot_nurse)),
                     aes(x=ggplot_x)) +
        geom_bar( aes( y=freq_plot_nurse,
                        fill=species_code), stat="identity",
                  position="dodge", alpha =0.85,  width =0.035*x_range )+#
        scale_fill_viridis_d(option="plasma")+
        labs(fill = "Species")+
        guides(fill=guide_legend(ncol=2)) +
        xlim(min(c(ggplot_x, oc_bio[[1]]))-2, 
             max(c(ggplot_x, oc_bio[[1]]))+2)+ 
        ylim(0, max(table_plot_n$freq_plot_nurse)+0.05)+ 
        # geom_segment(aes(x=oc_bio[[1]], y=0, 
        #                  xend=oc_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75, linetype="dashed")+
        # geom_segment(aes(x=cic_bio, y=0, 
        #                  xend=cic_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75)+
        xlab(paste0("Centroid ", biovariable))+
        ylab("Community relative abundance")+
        ggtitle(paste0(nurse_id, " ", biovariable))+
        # theme_minimal()+
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_text(size = 9), 
              legend.text  = element_text(size = 8),
              legend.key.size = unit(5, "mm"),
              panel.grid.minor=element_blank(),
              axis.text=element_text(color="black"),
              axis.title=element_text(color="black"),
              text=element_text(family="serif"))
      
      
      gg_list[[j]] <- gg_i
      names(gg_list)[j] <- paste0(nurse_id, "_", biovariable)
      oc_list[j] <- oc_bio[[1]]
      cic_list[j] <- cic_bio
      
    }
        ## add pairwise plots per locality
        
        legend <- cowplot::get_legend(gg_list[[1]])
        grid::grid.draw(legend)
        
        first <- gg_list[[1]]+
          geom_vline(xintercept = oc_list[[1]], linetype="dashed", 
                     color = "black", size=0.75)+
          geom_vline(xintercept = cic_list[[1]],
                     color = "black", size=0.75)+
          theme(legend.position = "none")
        
        second <- gg_list[[2]]+
          geom_vline(xintercept = oc_list[[2]], linetype="dashed", 
                     color = "black", size=0.75)+
          geom_vline(xintercept = cic_list[[2]],
                     color = "black", size=0.75)+
          theme(legend.position = "none")
        
        r <- rectGrob(gp=gpar(fill="white", col="white"))
        biplot <- gridExtra::grid.arrange(first, second, r, legend, 
                                          ncol=4, widths = unit(c(7.3,7.3,1,1.5), "cm"),
                                          heights = unit(7,'cm'))
        
        # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/",
        #               plot_list[i],"_", biovariable, "nurse.png"),
        #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
        # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/", 
        #               plot_list[i],"_", biovariable, "nurse.pdf"),
        #        plot=biplot, width = 19, height = 8, units="cm")
        # 
}

      
      
}


# write.table(cic_df, "./statistical_analyses/output/disequilibrium_centroid_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)



#* 1.3 add corresponding plots adding also multivariate CD ------------------------


multiv_total <-  read_delim( "./statistical_analyses/output/disequilibrium_nurse.csv", delim=";", col_names=T)
univ_total <-  read_delim( "./statistical_analyses/output/disequilibrium_centroid_univariate.csv", delim=";", col_names=T)

multiv_cd <- multiv_total%>%
  mutate(variable="multiv", 
         cd_norm = cd,
         oc= oc_x,
         cd_abs=cd)%>%
  dplyr::select(1:4, 10, 12, 9, 13, 11)


univ_oc <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("oc_b"))

names(univ_oc) <- names(univ_oc) %>%
  gsub("oc_b..", "bio", .) 

univ_oc <- univ_oc%>%
  tidyr::gather(key="variable", value="oc", 
         -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))


univ_cd <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_b"))

names(univ_cd) <- names(univ_cd) %>%
  gsub("cd_b..", "bio", .) 

univ_cd <- univ_cd%>%
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
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_abs_b"))

names(univ_cd_abs) <- names(univ_cd_abs) %>%
  gsub("cd_abs_b..", "bio", .) 

univ_cd_abs <- univ_cd_abs%>%
  tidyr::gather(key="variable", value="cd_abs", 
                -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))

univ_cdoc <- purrr::reduce(list(univ_oc, univ_cd,
                                univ_cd_abs),
                           dplyr::left_join,
                           by= c("plot", "plot_code", 
                                 "nurse", "nurse_code",
                                 "variable"))


univ_cdoc <- univ_cdoc%>%
  mutate(cd_norm=cd/oc)%>%
  dplyr::select(plot, plot_code, nurse, 
                nurse_code, variable, 
                oc, cd, cd_abs,
                cd_norm)

plot(univ_cdoc$cd, univ_cdoc$cd_abs)

total_cd <- rbind(multiv_cd, univ_cdoc)
plot_id <- unique(total_cd$plot_code)

for(i in 1:length(plot_id)){
  
  plot_cd <- total_cd%>%
    dplyr::filter(plot_code==plot_id[i])
  
  (uni_bar <- ggplot(data=plot_cd, 
                    aes(x=variable, y=cd_norm, fill=nurse)) +
    geom_bar(position = "dodge", stat="identity")+
    viridis::scale_fill_viridis(discrete = TRUE,option = 'E')+
    #♣coord_flip()+ to make horizontal plot
    labs(fill = "Community")+
    ylab("Climatic disequlibrium")+
    xlab("")+
    theme_minimal()+
    ggtitle(plot_id[i])+
    theme(plot.title = element_text(hjust = 0.5, color="black", size=9),
          legend.title = element_text(size = 5.7), 
          legend.text  = element_text(size = 5),
          legend.key.size = unit(2, "mm"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title=element_text(color="black", size=7),
          axis.text.x =element_text(color="black", size=4),
          axis.text.y =element_text(color="black", size=5),
          text=element_text(family="serif")))
  
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", plot_id[i], "_centroid.png"), 
  #        plot=uni_bar, width = 10, height = 6, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", plot_id[i], "_centroid.pdf"),
  #        plot=uni_bar, width = 8, height = 6, units="cm")
  
  
}



# 2. CIC, CD estimation with NICHE MODE -------------------------------


#* 2.1 prepare tables ------------------------------------------------------

table <-  read_delim("./niche_modelling/output/mode_centroid_distances_univariate.csv", delim=";", col_names=T)
codes <- read_delim("./niche_modelling/output/perc95_distances_in0.csv", delim=";", col_names=T)
codes <- codes%>%
  dplyr::select(community, code)

table <- left_join(table, codes, by="community")%>%
  distinct()
head(table)
names(table)
nrow(table)

table <- table%>%
  rename(plot= community,
         plot_code=code)%>%
  mutate(species_code= paste0(substring(genus, 1, 1),
                              substring(sp, 1,4)),
         nurse_code= paste0(plot_code, "_", nurse))

table_fq <- table %>%
  group_by(plot_code)%>%
  summarise(plot_n = sum(n))

table_fqq <- table%>%
  group_by(plot_code, nurse)%>%
  summarise(plot_nurse_n= sum(n))

table1 <- full_join(table, table_fq, by="plot_code")
table2 <- full_join(table1, table_fqq, by=c("plot_code", "nurse"))
table <- table2
rm(table1, table2)

table <- table%>%
  mutate(freq_plot = n/plot_n,
         freq_plot_nurse = n/plot_nurse_n)

#* 2.2 CIC, CD estimation---------------------------------------------------

# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

cic_df <-  table%>%
  dplyr::select(c(1:13, 71:77, 33:70))%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, names(table)[52:70])%>%
  distinct()%>%
  rename(oc_bio01= population_bio01, 
         oc_bio02= population_bio02,
         oc_bio03= population_bio03, 
         oc_bio04= population_bio04,
         oc_bio05= population_bio05, 
         oc_bio06= population_bio06,
         oc_bio07= population_bio07, 
         oc_bio08= population_bio08,
         oc_bio09= population_bio09, 
         oc_bio10= population_bio10,
         oc_bio11= population_bio11, 
         oc_bio12= population_bio12,
         oc_bio13= population_bio13, 
         oc_bio14= population_bio14,
         oc_bio15= population_bio15, 
         oc_bio16= population_bio16,
         oc_bio17= population_bio17, 
         oc_bio18= population_bio18,
         oc_bio19= population_bio19)

plot_list <- unique(cic_df$plot_code)
nurse_list <- unique(cic_df$nurse)

for (k in 1:19){
  
  biovariable <- if(k<=9){
    paste0("bio0", k )}else{
      paste0("bio", k )}
  
  for(i in 1:length(plot_list)) {
    
    table_plot <- table%>%
      filter(plot_code== plot_list[i])
    
    gg_list <- list()
    oc_list <- list()
    cic_list <- list()
    
    
    for(j in 1:length(nurse_list)){
      
      table_plot_n <- table_plot%>%
        filter(nurse== nurse_list[j])
      
      nurse_id <- unique(table_plot_n$nurse_code) 
      
      
      oc_bio <- cic_df[cic_df$nurse_code==nurse_id, paste0("oc_", biovariable)]
      cic_bio <- wt.mean(x=table_plot_n[, paste0("mode_", biovariable)][[1]], # paste0("mode_", biovariable)
                         wt=table_plot_n[,"freq_plot_nurse"][[1]])# freq_plot_nurse
      cic_df[cic_df$nurse_code==nurse_id, paste0("cic_", biovariable)] <- cic_bio
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_abs_", biovariable)] <- abs(cic_bio) - abs(oc_bio)
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_", biovariable)] <- cic_bio - oc_bio
      
      ggplot_x <- table_plot_n[, paste0("mode_", biovariable)][[1]]
      x_range <- max(c(ggplot_x, oc_bio[[1]]))-min(c(ggplot_x, oc_bio[[1]]))
      gg_i <- ggplot(data= table_plot_n %>%
                       dplyr::mutate(freq_plot_nurse = case_when(
                         freq_plot_nurse == 0~NA_real_,
                         freq_plot_nurse != 0~freq_plot_nurse)),
                     aes(x=ggplot_x)) +
        geom_bar( aes( y=freq_plot_nurse,
                       fill=species_code), stat="identity",
                  position="dodge", alpha =0.85,  width =0.035*x_range )+#
        scale_fill_viridis_d(option="plasma")+
        labs(fill = "Species")+
        guides(fill=guide_legend(ncol=2)) +
        xlim(min(c(ggplot_x, oc_bio[[1]]))-2, 
             max(c(ggplot_x, oc_bio[[1]]))+2)+ 
        ylim(0, max(table_plot_n$freq_plot_nurse)+0.05)+ 
        # geom_segment(aes(x=oc_bio[[1]], y=0, 
        #                  xend=oc_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75, linetype="dashed")+
        # geom_segment(aes(x=cic_bio, y=0, 
        #                  xend=cic_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75)+
        xlab(paste0("Centroid ", biovariable))+
        ylab("Species relative abundance")+
        ggtitle(paste0(nurse_id, " ", biovariable))+
        # theme_minimal()+
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_text(size = 9), 
              legend.text  = element_text(size = 8),
              legend.key.size = unit(5, "mm"),
              panel.grid.minor=element_blank(),
              axis.text=element_text(color="black"),
              axis.title=element_text(color="black"),
              text=element_text(family="serif"))
      
      
      gg_list[[j]] <- gg_i
      names(gg_list)[j] <- paste0(nurse_id, "_", biovariable)
      oc_list[j] <- oc_bio[[1]]
      cic_list[j] <- cic_bio
      
    }
    ## add pairwise plots per locality
    
    legend <- cowplot::get_legend(gg_list[[1]])
    grid::grid.draw(legend)
    
    first <- gg_list[[1]]+
      geom_vline(xintercept = oc_list[[1]], linetype="dashed", 
                 color = "black", size=0.75)+
      geom_vline(xintercept = cic_list[[1]],
                 color = "black", size=0.75)+
      theme(legend.position = "none")
    
    second <- gg_list[[2]]+
      geom_vline(xintercept = oc_list[[2]], linetype="dashed", 
                 color = "black", size=0.75)+
      geom_vline(xintercept = cic_list[[2]],
                 color = "black", size=0.75)+
      theme(legend.position = "none")
    
    r <- rectGrob(gp=gpar(fill="white", col="white"))
    biplot <- gridExtra::grid.arrange(first, second, r, legend, 
                                      ncol=4, widths = unit(c(7.3,7.3,1,1.5), "cm"),
                                      heights = unit(7,'cm'))
    
    # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/", 
    #               plot_list[i],"_", biovariable, "mode.png"), 
    #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
    # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/", 
    #               plot_list[i],"_", biovariable, "mode.pdf"),
    #        plot=biplot, width = 19, height = 8, units="cm")
    
  }
  
  
  
}


# write.table(cic_df, "./statistical_analyses/output/disequilibrium_mode_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)



#* 2.3 add corresponding plots adding also multivariate CD ------------------------


multiv_total <-  read_delim( "./statistical_analyses/output/disequilibrium_mode_centroid.csv", delim=";", col_names=T)
univ_total <-  read_delim( "./statistical_analyses/output/disequilibrium_mode_univariate.csv", delim=";", col_names=T)


multiv_cd <- multiv_total%>%
  mutate(variable="multiv", 
         cd_norm = cd,
         oc= oc_x,
         cd_abs=cd)%>%
  dplyr::select(1:4, 10, 12, 9, 13, 11)

univ_oc <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("oc_b"))

names(univ_oc) <- names(univ_oc) %>%
  gsub("oc_b..", "bio", .) 

univ_oc <- univ_oc%>%
  tidyr::gather(key="variable", value="oc", 
                -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))


univ_cd <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_b"))

names(univ_cd) <- names(univ_cd) %>%
  gsub("cd_b..", "bio", .) 

univ_cd <- univ_cd%>%
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
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_abs_b"))

names(univ_cd_abs) <- names(univ_cd_abs) %>%
  gsub("cd_abs_b..", "bio", .) 

univ_cd_abs <- univ_cd_abs%>%
  tidyr::gather(key="variable", value="cd_abs", 
                -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))

univ_cdoc <- purrr::reduce(list(univ_oc, univ_cd,
                                univ_cd_abs),
                           dplyr::left_join,
                           by= c("plot", "plot_code", 
                                 "nurse", "nurse_code",
                                 "variable"))


univ_cdoc <- univ_cdoc%>%
  mutate(cd_norm=cd/oc)%>%
  dplyr::select(plot, plot_code, nurse, 
                nurse_code, variable, 
                oc, cd, cd_abs,
                cd_norm)

plot(univ_cdoc$cd, univ_cdoc$cd_abs)

total_cd <- rbind(multiv_cd, univ_cdoc)
plot_id <- unique(total_cd$plot_code)

for(i in 1:length(plot_id)){
  
  plot_cd <- total_cd%>%
    dplyr::filter(plot_code==plot_id[i])
  
  (uni_bar <- ggplot(data=plot_cd, 
                    aes(x=variable, y=cd_norm, fill=nurse)) +
    geom_bar(position = "dodge", stat="identity")+
    viridis::scale_fill_viridis(discrete = TRUE,option = 'E')+
    #♣coord_flip()+ to make horizontal plot
    labs(fill = "Community")+
    ylab("Climatic disequlibrium")+
    xlab("")+
    theme_minimal()+
    ggtitle(plot_id[i])+
    theme(plot.title = element_text(hjust = 0.5, color="black", size=9),
          legend.title = element_text(size = 5.7), 
          legend.text  = element_text(size = 5),
          legend.key.size = unit(2, "mm"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title=element_text(color="black", size=7),
          axis.text.x =element_text(color="black", size=4),
          axis.text.y =element_text(color="black", size=5),
          text=element_text(family="serif")))
  
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", 
  #               plot_id[i], "_mode.png"), 
  #        plot=uni_bar, width = 10, height = 6, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", 
  #               plot_id[i], "_mode.pdf"),
  #        plot=uni_bar, width = 8, height = 6, units="cm")
  # 
  # 
}




# 3. CIC, CD estimation with NICHE PERCENTILE 95 ---------------------------

#* 3.1 prepare tables ------------------------------------------------------

table <-  read_delim("./niche_modelling/output/perc95_distances_in0_univariate.csv", delim=";", col_names=T)


table <- table%>%
  rename(plot= community,
         plot_code=code)%>%
  mutate(species_code= paste0(substring(genus, 1, 1),
                              substring(sp, 1,4)),
         nurse_code= paste0(plot_code, "_", nurse))

table_fq <- table %>%
  group_by(plot_code)%>%
  summarise(plot_n = sum(n))

table_fqq <- table%>%
  group_by(plot_code, nurse)%>%
  summarise(plot_nurse_n= sum(n))

table1 <- full_join(table, table_fq, by="plot_code")
table2 <- full_join(table1, table_fqq, by=c("plot_code", "nurse"))
table <- table2
rm(table1, table2)

table <- table%>%
  mutate(freq_plot = n/plot_n,
         freq_plot_nurse = n/plot_nurse_n)

#* 3.2 CIC, CD estimation --------------------------------------------------

# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

cic_df <-  table%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, names(table)[15:33])%>%
  distinct()%>%
  rename(oc_bio01= pop_bio01, 
         oc_bio02= pop_bio02,
         oc_bio03= pop_bio03, 
         oc_bio04= pop_bio04,
         oc_bio05= pop_bio05, 
         oc_bio06= pop_bio06,
         oc_bio07= pop_bio07, 
         oc_bio08= pop_bio08,
         oc_bio09= pop_bio09, 
         oc_bio10= pop_bio10,
         oc_bio11= pop_bio11, 
         oc_bio12= pop_bio12,
         oc_bio13= pop_bio13, 
         oc_bio14= pop_bio14,
         oc_bio15= pop_bio15, 
         oc_bio16= pop_bio16,
         oc_bio17= pop_bio17, 
         oc_bio18= pop_bio18,
         oc_bio19= pop_bio19)

plot_list <- unique(cic_df$plot_code)
nurse_list <- unique(cic_df$nurse)



for (k in 1:19){
  
  biovariable <- if(k<=9){
    paste0("bio0", k )}else{
      paste0("bio", k )}
  
  for(i in 1:length(plot_list)) {
    
    table_plot <- table%>%
      filter(plot_code== plot_list[i])
    
    gg_list <- list()
    oc_list <- list()
    cic_list <- list()
    
    
    for(j in 1:length(nurse_list)){
      
      table_plot_n <- table_plot%>%
        filter(nurse== nurse_list[j])
      
      nurse_id <- unique(table_plot_n$nurse_code) 
      
      
      oc_bio <- cic_df[cic_df$nurse_code==nurse_id, paste0("oc_", biovariable)]
      cic_bio <- wt.mean(x=table_plot_n[, paste0("closest_p95_", biovariable)][[1]], # paste0("mode_", biovariable)
                         wt=table_plot_n[,"freq_plot_nurse"][[1]])# freq_plot_nurse
     
      cic_df[cic_df$nurse_code==nurse_id, paste0("cic_", biovariable)] <- cic_bio
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_abs_", biovariable)] <- abs(cic_bio) - abs(oc_bio)
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_", biovariable)] <- cic_bio - oc_bio
      
      ggplot_x <- table_plot_n[, paste0("closest_p95_", biovariable)][[1]]
      x_range <- max(c(ggplot_x, oc_bio[[1]]))-min(c(ggplot_x, oc_bio[[1]]))
      gg_i <- ggplot(data= table_plot_n %>%
                       dplyr::mutate(freq_plot_nurse = case_when(
                         freq_plot_nurse == 0~NA_real_,
                         freq_plot_nurse != 0~freq_plot_nurse)),
                     aes(x=ggplot_x)) +
        geom_bar( aes( y=freq_plot_nurse,
                       fill=species_code), stat="identity",
                  position="dodge", alpha =0.85,  width =0.035*x_range )+#
        scale_fill_viridis_d(option="plasma")+
        labs(fill = "Species")+
        guides(fill=guide_legend(ncol=2)) +
        xlim(min(c(ggplot_x, oc_bio[[1]]))-2, 
             max(c(ggplot_x, oc_bio[[1]]))+2)+ 
        ylim(0, max(table_plot_n$freq_plot_nurse)+0.05)+ 
        # geom_segment(aes(x=oc_bio[[1]], y=0, 
        #                  xend=oc_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75, linetype="dashed")+
        # geom_segment(aes(x=cic_bio, y=0, 
        #                  xend=cic_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75)+
        xlab(paste0("Closest perc 95 ", biovariable))+
        ylab("Species relative abundance")+
        ggtitle(paste0(nurse_id, " ", biovariable))+
        # theme_minimal()+
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_text(size = 9), 
              legend.text  = element_text(size = 8),
              legend.key.size = unit(5, "mm"),
              panel.grid.minor=element_blank(),
              axis.text=element_text(color="black"),
              axis.title=element_text(color="black"),
              text=element_text(family="serif"))
      
      
      gg_list[[j]] <- gg_i
      names(gg_list)[j] <- paste0(nurse_id, "_", biovariable)
      oc_list[j] <- oc_bio[[1]]
      cic_list[j] <- cic_bio
      
    }
    ## add pairwise plots per locality
    
    legend <- cowplot::get_legend(gg_list[[1]])
    grid::grid.draw(legend)

    first <- gg_list[[1]]+
      geom_vline(xintercept = oc_list[[1]], linetype="dashed",
                 color = "black", size=0.75)+
      geom_vline(xintercept = cic_list[[1]],
                 color = "black", size=0.75)+
      theme(legend.position = "none")

    second <- gg_list[[2]]+
      geom_vline(xintercept = oc_list[[2]], linetype="dashed",
                 color = "black", size=0.75)+
      geom_vline(xintercept = cic_list[[2]],
                 color = "black", size=0.75)+
      theme(legend.position = "none")

    r <- rectGrob(gp=gpar(fill="white", col="white"))
    biplot <- gridExtra::grid.arrange(first, second, r, legend,
                                      ncol=4, widths = unit(c(7.3,7.3,1,1.5), "cm"),
                                      heights = unit(7,'cm'))

    # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/",
    #               plot_list[i],"_", biovariable, "perc_95_in0.png"),
    #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
    # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/",
    #               plot_list[i],"_", biovariable, "perc_95_in0.pdf"),
    #        plot=biplot, width = 19, height = 8, units="cm")

  }
  
  
  
}


# write.table(cic_df, "./statistical_analyses/output/disequilibrium_perc_95_in0_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


#* 3.3 add corresponding plots adding multivariate cd --------------------------


multiv_total <-  read_delim( "./statistical_analyses/output/disequilibrium_perc95_in0.csv", delim=";", col_names=T)
univ_total <-  read_delim( "./statistical_analyses/output/disequilibrium_perc_95_in0_univariate.csv", delim=";", col_names=T)

multiv_cd <- multiv_total%>%
  mutate(variable="multiv", 
         cd_norm = cd,
         oc= oc_x,
         cd_abs=cd)%>%
  dplyr::select(1:4, 10, 12, 9, 13, 11)


univ_oc <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("oc_b"))

names(univ_oc) <- names(univ_oc) %>%
  gsub("oc_b..", "bio", .) 

univ_oc <- univ_oc%>%
  tidyr::gather(key="variable", value="oc", 
                -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))


univ_cd <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_b"))

names(univ_cd) <- names(univ_cd) %>%
  gsub("cd_b..", "bio", .) 

univ_cd <- univ_cd%>%
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
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_abs_b"))

names(univ_cd_abs) <- names(univ_cd_abs) %>%
  gsub("cd_abs_b..", "bio", .) 

univ_cd_abs <- univ_cd_abs%>%
  tidyr::gather(key="variable", value="cd_abs", 
                -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))

univ_cdoc <- purrr::reduce(list(univ_oc, univ_cd,
                                univ_cd_abs),
                           dplyr::left_join,
                           by= c("plot", "plot_code", 
                                 "nurse", "nurse_code",
                                 "variable"))


univ_cdoc <- univ_cdoc%>%
  mutate(cd_norm=cd/oc)%>%
  dplyr::select(plot, plot_code, nurse, 
                nurse_code, variable, 
                oc, cd, cd_abs,
                cd_norm)

plot(univ_cdoc$cd, univ_cdoc$cd_abs)

total_cd <- rbind(multiv_cd, univ_cdoc)
plot_id <- unique(total_cd$plot_code)


for(i in 1:length(plot_id)){
  
  plot_cd <- total_cd%>%
    dplyr::filter(plot_code==plot_id[i])
  
  (uni_bar <- ggplot(data=plot_cd, 
                    aes(x=variable, y=cd_norm, fill=nurse)) +
    geom_bar(position = "dodge", stat="identity")+
    viridis::scale_fill_viridis(discrete = TRUE,option = 'E')+
    #♣coord_flip()+ to make horizontal plot
    labs(fill = "Community")+
    ylab("Climatic disequlibrium")+
    xlab("")+
    theme_minimal()+
    ggtitle(plot_id[i])+
    theme(plot.title = element_text(hjust = 0.5, color="black", size=9),
          legend.title = element_text(size = 5.7), 
          legend.text  = element_text(size = 5),
          legend.key.size = unit(2, "mm"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title=element_text(color="black", size=7),
          axis.text.x =element_text(color="black", size=4),
          axis.text.y =element_text(color="black", size=5),
          text=element_text(family="serif")))
  
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", 
  #               plot_id[i], "_perc_95_in0.png"), 
  #        plot=uni_bar, width = 10, height = 6, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", 
  #               plot_id[i], "_perc_95_in0.pdf"),
  #        plot=uni_bar, width = 8, height = 6, units="cm")
  
  
}




# 4. CIC, CD estimation with NICHE PERCENTILE 90 --------------------------


#* 4.1 prepare tables ------------------------------------------------------


table <-  read_delim("./niche_modelling/output/perc90_distances_in0_univariate.csv", delim=";", col_names=T)

table <- table%>%
  rename(plot= community,
         plot_code=code)%>%
  mutate(species_code= paste0(substring(genus, 1, 1),
                              substring(sp, 1,4)),
         nurse_code= paste0(plot_code, "_", nurse))

table_fq <- table %>%
  group_by(plot_code)%>%
  summarise(plot_n = sum(n))

table_fqq <- table%>%
  group_by(plot_code, nurse)%>%
  summarise(plot_nurse_n= sum(n))

table1 <- full_join(table, table_fq, by="plot_code")
table2 <- full_join(table1, table_fqq, by=c("plot_code", "nurse"))
table <- table2
rm(table1, table2)

table <- table%>%
  mutate(freq_plot = n/plot_n,
         freq_plot_nurse = n/plot_nurse_n)

#* 4.2 CIC, CD estimation -------------------------------------------------

# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

cic_df <-  table%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, names(table)[15:33])%>%
  distinct()%>%
  rename(oc_bio01= pop_bio01, 
         oc_bio02= pop_bio02,
         oc_bio03= pop_bio03, 
         oc_bio04= pop_bio04,
         oc_bio05= pop_bio05, 
         oc_bio06= pop_bio06,
         oc_bio07= pop_bio07, 
         oc_bio08= pop_bio08,
         oc_bio09= pop_bio09, 
         oc_bio10= pop_bio10,
         oc_bio11= pop_bio11, 
         oc_bio12= pop_bio12,
         oc_bio13= pop_bio13, 
         oc_bio14= pop_bio14,
         oc_bio15= pop_bio15, 
         oc_bio16= pop_bio16,
         oc_bio17= pop_bio17, 
         oc_bio18= pop_bio18,
         oc_bio19= pop_bio19)

plot_list <- unique(cic_df$plot_code)
nurse_list <- unique(cic_df$nurse)



for (k in 1:19){
  
  biovariable <- if(k<=9){
    paste0("bio0", k )}else{
      paste0("bio", k )}
  
  for(i in 1:length(plot_list)) {
    
    table_plot <- table%>%
      filter(plot_code== plot_list[i])
    
    gg_list <- list()
    oc_list <- list()
    cic_list <- list()
    
    
    for(j in 1:length(nurse_list)){
      
      table_plot_n <- table_plot%>%
        filter(nurse== nurse_list[j])
      
      nurse_id <- unique(table_plot_n$nurse_code) 
      
      
      oc_bio <- cic_df[cic_df$nurse_code==nurse_id, paste0("oc_", biovariable)]
      cic_bio <- wt.mean(x=table_plot_n[, paste0("closest_p90_", biovariable)][[1]], # paste0("mode_", biovariable)
                         wt=table_plot_n[,"freq_plot_nurse"][[1]])# freq_plot_nurse
      cic_df[cic_df$nurse_code==nurse_id, paste0("cic_", biovariable)] <- cic_bio
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_abs_", biovariable)] <- abs(cic_bio) - abs(oc_bio)
      cic_df[cic_df$nurse_code==nurse_id, paste0("cd_", biovariable)] <- cic_bio - oc_bio
      
      ggplot_x <- table_plot_n[, paste0("closest_p90_", biovariable)][[1]]
      x_range <- max(c(ggplot_x, oc_bio[[1]]))-min(c(ggplot_x, oc_bio[[1]]))
      gg_i <- ggplot(data= table_plot_n %>%
                       dplyr::mutate(freq_plot_nurse = case_when(
                         freq_plot_nurse == 0~NA_real_,
                         freq_plot_nurse != 0~freq_plot_nurse)),
                     aes(x=ggplot_x)) +
        geom_bar( aes( y=freq_plot_nurse,
                       fill=species_code), stat="identity",
                  position="dodge", alpha =0.85,  width =0.035*x_range )+#
        scale_fill_viridis_d(option="plasma")+
        labs(fill = "Species")+
        guides(fill=guide_legend(ncol=2)) +
        xlim(min(c(ggplot_x, oc_bio[[1]]))-2, 
             max(c(ggplot_x, oc_bio[[1]]))+2)+ 
        ylim(0, max(table_plot_n$freq_plot_nurse)+0.05)+ 
        # geom_segment(aes(x=oc_bio[[1]], y=0, 
        #                  xend=oc_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75, linetype="dashed")+
        # geom_segment(aes(x=cic_bio, y=0, 
        #                  xend=cic_bio[[1]], 
        #                  yend=max(table_plot_n$freq_plot_nurse)+0.05),
        #              color = "black", size=0.75)+
        xlab(paste0("Closest perc 90 ", biovariable))+
        ylab("Species relative abundance")+
        ggtitle(paste0(nurse_id, " ", biovariable))+
        # theme_minimal()+
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_text(size = 9), 
              legend.text  = element_text(size = 8),
              legend.key.size = unit(5, "mm"),
              panel.grid.minor=element_blank(),
              axis.text=element_text(color="black"),
              axis.title=element_text(color="black"),
              text=element_text(family="serif"))
      
      
      gg_list[[j]] <- gg_i
      names(gg_list)[j] <- paste0(nurse_id, "_", biovariable)
      oc_list[j] <- oc_bio[[1]]
      cic_list[j] <- cic_bio
      
    }
    ## add pairwise plots per locality
    
    legend <- cowplot::get_legend(gg_list[[1]])
    grid::grid.draw(legend)
    
    first <- gg_list[[1]]+
      geom_vline(xintercept = oc_list[[1]], linetype="dashed", 
                 color = "black", size=0.75)+
      geom_vline(xintercept = cic_list[[1]],
                 color = "black", size=0.75)+
      theme(legend.position = "none")
    
    second <- gg_list[[2]]+
      geom_vline(xintercept = oc_list[[2]], linetype="dashed", 
                 color = "black", size=0.75)+
      geom_vline(xintercept = cic_list[[2]],
                 color = "black", size=0.75)+
      theme(legend.position = "none")
    
    r <- rectGrob(gp=gpar(fill="white", col="white"))
    biplot <- gridExtra::grid.arrange(first, second, r, legend, 
                                      ncol=4, widths = unit(c(7.3,7.3,1,1.5), "cm"),
                                      heights = unit(7,'cm'))
    
    # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/", 
    #               plot_list[i],"_", biovariable, "perc_90_in0.png"), 
    #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
    # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/per_variable/", 
    #               plot_list[i],"_", biovariable, "perc_90_in0.pdf"),
    #        plot=biplot, width = 19, height = 8, units="cm")
    
  }
  
  
  
}


#write.table(cic_df, "./statistical_analyses/output/disequilibrium_perc_90_in0_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)



#* 4.3 add corresponding plots adding multivariate cd-------------------------------

multiv_total <-  read_delim( "./statistical_analyses/output/disequilibrium_perc90_in0.csv", delim=";", col_names=T)
univ_total <-  read_delim( "./statistical_analyses/output/disequilibrium_perc_90_in0_univariate.csv", delim=";", col_names=T)

multiv_cd <- multiv_total%>%
  mutate(variable="multiv", 
         cd_norm = cd,
         oc= oc_x,
         cd_abs=cd)%>%
  dplyr::select(1:4, 10, 12, 9, 13, 11)

univ_oc <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("oc_b"))

names(univ_oc) <- names(univ_oc) %>%
  gsub("oc_b..", "bio", .) 

univ_oc <- univ_oc%>%
  tidyr::gather(key="variable", value="oc", 
                -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))


univ_cd <- univ_total%>%
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_b"))

names(univ_cd) <- names(univ_cd) %>%
  gsub("cd_b..", "bio", .) 

univ_cd <- univ_cd%>%
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
  select(plot,
         plot_code,
         nurse,
         nurse_code,
         contains("cd_abs_b"))

names(univ_cd_abs) <- names(univ_cd_abs) %>%
  gsub("cd_abs_b..", "bio", .) 

univ_cd_abs <- univ_cd_abs%>%
  tidyr::gather(key="variable", value="cd_abs", 
                -c(plot, plot_code, nurse, nurse_code))%>%
  dplyr::filter(variable %in% c("bio01", "bio04",
                                "bio05", "bio06",
                                "bio07", "bio10",
                                "bio11", "bio12",
                                "bio13", "bio14",
                                "bio15", "bio16",
                                "bio17"))

univ_cdoc <- purrr::reduce(list(univ_oc, univ_cd,
                                univ_cd_abs),
                           dplyr::left_join,
                           by= c("plot", "plot_code", 
                                 "nurse", "nurse_code",
                                 "variable"))


univ_cdoc <- univ_cdoc%>%
  mutate(cd_norm=cd/oc)%>%
  dplyr::select(plot, plot_code, nurse, 
                nurse_code, variable, 
                oc, cd, cd_abs,
                cd_norm)

plot(univ_cdoc$cd, univ_cdoc$cd_abs)

total_cd <- rbind(multiv_cd, univ_cdoc)
plot_id <- unique(total_cd$plot_code)


for(i in 1:length(plot_id)){
  
  plot_cd <- total_cd%>%
    dplyr::filter(plot_code==plot_id[i])
  
  uni_bar <- ggplot(data=plot_cd, 
                    aes(x=variable, y=cd_norm, fill=nurse)) +
    geom_bar(position = "dodge", stat="identity")+
    viridis::scale_fill_viridis(discrete = TRUE,option = 'E')+
    #♣coord_flip()+ to make horizontal plot
    labs(fill = "Community")+
    ylab("Climatic disequlibrium")+
    xlab("")+
    theme_minimal()+
    ggtitle(plot_id[i])+
    theme(plot.title = element_text(hjust = 0.5, color="black", size=9),
          legend.title = element_text(size = 5.7), 
          legend.text  = element_text(size = 5),
          legend.key.size = unit(2, "mm"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title=element_text(color="black", size=7),
          axis.text.x =element_text(color="black", size=4),
          axis.text.y =element_text(color="black", size=5),
          text=element_text(family="serif"))
  
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", 
  #               plot_id[i], "_perc_90_in0.png"), 
  #        plot=uni_bar, width = 10, height = 6, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/univariate/all_variables/", 
  #               plot_id[i], "_perc_90_in0.pdf"),
  #        plot=uni_bar, width = 8, height = 6, units="cm")
  
  
}

