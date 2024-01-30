


#*************************************#
### Maria Angeles Perez-Navarro
### CREAF  02/2021
### Community climatic diagrams multivariate
#*************************************#


library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)

## needed functions #

source('./niche_modelling/niche_analyses/niche_required_functions.R')

## create folder
dir.create("./statistical_analyses/output/")
dir.create("./statistical_analyses/output/climatic_diagrams/multivariate/")

# 1. CIC, CD estimation with NICHE CENTROIDS -------------------------------


#* 1.1 prepare tables ------------------------------------------------------


table <-  read_delim("./niche_modelling/output/centroid_distances.csv", delim=";", col_names=T)
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

write.table(table%>%
              dplyr::select(plot, site, plot_code,
                            recruit_sp, species_code),
            "./statistical_analyses/output/table_codes.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


#* 1.2 CIC, CD estimation---------------------------------------------------


# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

centroid_df <-  table%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, Axis1, Axis2)%>%
  distinct()%>%
  rename(oc_x= Axis1, 
         oc_y= Axis2)

plot_list <- unique(centroid_df$plot_code)
nurse_list <- unique(centroid_df$nurse)

for(i in 1:length(plot_list)) {
  table_plot <- table%>%
    filter(plot_code== plot_list[i])
  
  gg_list <- list()
  
  
  for(j in 1:length(nurse_list)){
    table_plot_n <- table_plot%>%
      filter(nurse== nurse_list[j])
    
    nurse_id <- unique(table_plot_n$nurse_code)
    cic <- COGravity(x=table_plot_n$centroid_pca1, 
                     y=table_plot_n$centroid_pca2, 
                     wt=table_plot_n$freq_plot_nurse)
  
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_x"] <- cic[1]
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_y"] <- cic[3]
  
    centroid_df_plot <- centroid_df%>%
      filter(nurse_code== nurse_id)
    centroid_df[centroid_df$nurse_code==nurse_id, "cd"]<- GMINdistance(centroid_df_plot$cic_x, centroid_df_plot$cic_y,
                                                             centroid_df_plot$oc_x, centroid_df_plot$oc_y)[1]
    
    
    table_plot_n$species_code <- factor(table_plot_n$species_code)
    cols <- colorRampPalette(brewer.pal(12, "Spectral"))
    mypal <- cols(length(unique(table_plot_n$species_code)))
    
    gg_i <- ggplot(table_plot_n %>%
                     dplyr::mutate(freq_plot_nurse = case_when(
                       freq_plot_nurse == 0~NA_real_,
                       freq_plot_nurse != 0~freq_plot_nurse))) +
      geom_point( aes(x=centroid_pca1, y=centroid_pca2, 
                      color=species_code, size=freq_plot_nurse*10
      ), alpha =0.85)+ 
      #scale_color_manual(values = mypal)+
      scale_color_viridis_d(option="plasma")+
      labs(color = "Species")+
      guides(color=guide_legend(ncol=2)) +
      xlim(min(table_plot_n$centroid_pca1)-2, 
           max(table_plot_n$centroid_pca1)+2)+ 
      ylim(min(table_plot_n$centroid_pca2)-2, 
           max(table_plot_n$centroid_pca2)+2)+ 
      scale_size(range = c(2, 7))+
      scale_alpha_continuous("freq_plot_nurse", range=c(0.2, 0.9))+
      guides(size = FALSE,
             alpha = FALSE
             )+
      geom_point(data=centroid_df_plot, aes(x=cic_x, y=cic_y),
                 color="black", size=2, shape=19)+
      geom_point(data=centroid_df_plot, aes(x=oc_x, y=oc_y), 
                 color="black", size=2, shape=17)+
      xlab("Climatic Axis 1 (48.2%)")+
      ylab("Climatic Axis 2 (26.7%)")+
      ggtitle(nurse_id)+
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_text(size = 9), 
            legend.text  = element_text(size = 8),
            legend.key.size = unit(5, "mm"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(color="black"),
            text=element_text(family="serif"))
      
    gg_list[[j]] <- gg_i
  
}

  legend <- cowplot::get_legend(gg_list[[1]])
  grid::grid.draw(legend)
  
  first <- gg_list[[1]]+
    theme(legend.position = "none")
  
  second <- gg_list[[2]]+
    theme(legend.position = "none")
  
  r <- rectGrob(gp=gpar(fill="white", col="white"))
  biplot <- gridExtra::grid.arrange(first, second, r, legend, 
                                    ncol=4, widths = unit(c(6.75,6.75,1.25,1.5), "cm"),heights = unit(7,'cm'))
 
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/multivariate/", plot_list[i], "nurse.png"), 
  #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/multivariate/", plot_list[i], "nurse.pdf"),
  #        plot=biplot, width = 19, height = 8, units="cm")
  
   
}

# write.table(centroid_df, "./statistical_analyses/output/disequilibrium_nurse.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


# 2. CIC, CD estimation with NICHE MODE -------------------------------

#* 2.1 prepare tables ------------------------------------------------------

table <-  read_delim("./niche_modelling/output/centroid_distances.csv", delim=";", col_names=T)
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


#* 2.2 CIC, CD estimation---------------------------------------------------

# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

centroid_df <-  table%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, Axis1, Axis2)%>%
  distinct()%>%
  rename(oc_x= Axis1, 
         oc_y= Axis2)

plot_list <- unique(centroid_df$plot_code)
nurse_list <- unique(centroid_df$nurse)

for(i in 1:length(plot_list)) {
  table_plot <- table%>%
    filter(plot_code== plot_list[i])
  
  gg_list <- list()
  
  
  for(j in 1:length(nurse_list)){
    table_plot_n <- table_plot%>%
      filter(nurse== nurse_list[j])
    
    nurse_id <- unique(table_plot_n$nurse_code)
    cic <- COGravity(x=table_plot_n$mode_pca1, 
                               y=table_plot_n$mode_pca2, 
                               wt=table_plot_n$freq_plot_nurse)
    
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_x"] <- cic[1]
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_y"] <- cic[3]
    
    centroid_df_plot <- centroid_df%>%
      filter(nurse_code== nurse_id)
    centroid_df[centroid_df$nurse_code==nurse_id, "cd"]<- GMINdistance(centroid_df_plot$cic_x, centroid_df_plot$cic_y,
                                                                       centroid_df_plot$oc_x, centroid_df_plot$oc_y)[1]
    
    
    table_plot_n$species_code <- factor(table_plot_n$species_code)
    cols <- colorRampPalette(brewer.pal(12, "Spectral"))
    mypal <- cols(length(unique(table_plot_n$species_code)))
    
    gg_i <- ggplot(table_plot_n %>%
                     dplyr::mutate(freq_plot_nurse = case_when(
                       freq_plot_nurse == 0~NA_real_,
                       freq_plot_nurse != 0~freq_plot_nurse))) +
      geom_point( aes(x=mode_pca1, y=mode_pca2, 
                      color=species_code, size=freq_plot_nurse*10
      ), alpha =0.85)+ 
      #scale_color_manual(values = mypal)+
      scale_color_viridis_d(option="plasma")+
      labs(color = "Species")+
      guides(color=guide_legend(ncol=2)) +
      xlim(min(table_plot_n$mode_pca1)-2, 
           max(table_plot_n$mode_pca1)+2)+ 
      ylim(min(table_plot_n$mode_pca2)-2, 
           max(table_plot_n$mode_pca2)+2)+ 
      scale_size(range = c(2, 7))+
      scale_alpha_continuous("freq_plot_nurse", range=c(0.2, 0.9))+
      guides(size = FALSE,
             alpha = FALSE
      )+
      geom_point(data=centroid_df_plot, aes(x=cic_x, y=cic_y),
                 color="black", size=2, shape=19)+
      geom_point(data=centroid_df_plot, aes(x=oc_x, y=oc_y), 
                 color="black", size=2, shape=17)+
      xlab("Climatic Axis 1 (48.2%)")+
      ylab("Climatic Axis 2 (26.7%)")+
      ggtitle(nurse_id)+
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_text(size = 9), 
            legend.text  = element_text(size = 8),
            legend.key.size = unit(5, "mm"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(color="black"),
            text=element_text(family="serif"))
    
    gg_list[[j]] <- gg_i
    
  }
  
  legend <- cowplot::get_legend(gg_list[[1]])
  grid::grid.draw(legend)
  
  first <- gg_list[[1]]+
    theme(legend.position = "none")
  
  second <- gg_list[[2]]+
    theme(legend.position = "none")
  
  r <- rectGrob(gp=gpar(fill="white", col="white"))
  biplot <- gridExtra::grid.arrange(first, second, r, legend, 
                                    ncol=4, widths = unit(c(6.75,6.75,1.25,1.5), "cm"),heights = unit(7,'cm'))
  
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/multivariate/", plot_list[i], "mode_centroid.png"), 
  #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/multivariate/", plot_list[i], "mode_centroid.pdf"),
  #        plot=biplot, width = 19, height = 8, units="cm")
  # 
  
}

# write.table(centroid_df, "./statistical_analyses/output/disequilibrium_mode_centroid.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)



# 3. CIC, CD estimation with NICHE PERCENTILE 95 ---------------------------

#* 3.1 prepare tables ------------------------------------------------------

table <-  read_delim("./niche_modelling/output/perc95_distances_in0.csv", delim=";", col_names=T)

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

#* 3.2 CIC, CD estimation --------------------------------------------------

# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

centroid_df <-  table%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, Axis1, Axis2)%>%
  distinct()%>%
  rename(oc_x= Axis1, 
         oc_y= Axis2)

plot_list <- unique(centroid_df$plot_code)
nurse_list <- unique(centroid_df$nurse)

for(i in 1:length(plot_list)) {
  table_plot <- table%>%
    filter(plot_code== plot_list[i])
  
  gg_list <- list()
  
  
  for(j in 1:length(nurse_list)){
    table_plot_n <- table_plot%>%
      filter(nurse== nurse_list[j])
    
    nurse_id <- unique(table_plot_n$nurse_code)
    cic <- COGravity(x=table_plot_n$perc_95_pca1, 
                               y=table_plot_n$perc_95_pca2, 
                               wt=table_plot_n$freq_plot_nurse)
    
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_x"] <- cic[1]
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_y"] <- cic[3]
    
    centroid_df_plot <- centroid_df%>%
      filter(nurse_code== nurse_id)
    centroid_df[centroid_df$nurse_code==nurse_id, "cd"]<- GMINdistance(centroid_df_plot$cic_x, centroid_df_plot$cic_y,
                                                                       centroid_df_plot$oc_x, centroid_df_plot$oc_y)[1]
    
    
    table_plot_n$species_code <- factor(table_plot_n$species_code)
    cols <- colorRampPalette(brewer.pal(12, "Spectral"))
    mypal <- cols(length(unique(table_plot_n$species_code)))
    
    gg_i <- ggplot(table_plot_n %>%
                     dplyr::mutate(freq_plot_nurse = case_when(
                       freq_plot_nurse == 0~NA_real_,
                       freq_plot_nurse != 0~freq_plot_nurse))) +
      geom_point( aes(x=perc_95_pca1, y=perc_95_pca2, 
                      color=species_code, size=freq_plot_nurse*10
      ), alpha =0.85)+ 
      #scale_color_manual(values = mypal)+
      scale_color_viridis_d(option="plasma")+
      labs(color = "Species")+
      guides(color=guide_legend(ncol=2)) +
      xlim(min(table_plot_n$perc_95_pca1)-2, 
           max(table_plot_n$perc_95_pca1)+2)+ 
      ylim(min(table_plot_n$perc_95_pca2)-2, 
           max(table_plot_n$perc_95_pca2)+2)+ 
      scale_size(range = c(2, 7))+
      scale_alpha_continuous("freq_plot_nurse", range=c(0.2, 0.9))+
      guides(size = FALSE,
             alpha = FALSE
      )+
      geom_point(data=centroid_df_plot, aes(x=cic_x, y=cic_y),
                 color="black", size=2, shape=19)+
      geom_point(data=centroid_df_plot, aes(x=oc_x, y=oc_y), 
                 color="black", size=2, shape=17)+
      xlab("Climatic Axis 1 (48.2%)")+
      ylab("Climatic Axis 2 (26.7%)")+
      ggtitle(nurse_id)+
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_text(size = 9), 
            legend.text  = element_text(size = 8),
            legend.key.size = unit(5, "mm"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(color="black"),
            text=element_text(family="serif"))
    
    gg_list[[j]] <- gg_i
    
  }
  
  legend <- cowplot::get_legend(gg_list[[1]])
  grid::grid.draw(legend)
  
  first <- gg_list[[1]]+
    theme(legend.position = "none")
  
  second <- gg_list[[2]]+
    theme(legend.position = "none")
  
  r <- rectGrob(gp=gpar(fill="white", col="white"))
  biplot <- gridExtra::grid.arrange(first, second, r, legend, 
                                    ncol=4, widths = unit(c(6.75,6.75,1.25,1.5), "cm"),heights = unit(7,'cm'))
  
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/multivariate/", plot_list[i], "perc95_in0.png"), 
  #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/multivariate/", plot_list[i], "perc95_in0.pdf"),
  #        plot=biplot, width = 19, height = 8, units="cm")
  # 
  
}

# write.table(centroid_df, "./statistical_analyses/output/disequilibrium_perc95_in0.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)



# 4. CIC, CD estimation with NICHE PERCENTILE 90 --------------------------


#* 4.1 prepare tables ------------------------------------------------------


table <-  read_delim("./niche_modelling/output/perc90_distances_in0.csv", delim=";", col_names=T)
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


#* 4.2 CIC, CD estimation -------------------------------------------------

# where CIC= community inferred climate
# where OC= observed climate
# where CD= climatic disequilibrium

centroid_df <-  table%>%
  dplyr::select(plot, plot_code, nurse,
                nurse_code, Axis1, Axis2)%>%
  distinct()%>%
  rename(oc_x= Axis1, 
         oc_y= Axis2)

plot_list <- unique(centroid_df$plot_code)
nurse_list <- unique(centroid_df$nurse)

for(i in 1:length(plot_list)) {
  table_plot <- table%>%
    filter(plot_code== plot_list[i])
  
  gg_list <- list()
  
  
  for(j in 1:length(nurse_list)){
    table_plot_n <- table_plot%>%
      filter(nurse== nurse_list[j])
    
    nurse_id <- unique(table_plot_n$nurse_code)
    cic <- COGravity(x=table_plot_n$perc_90_pca1, 
                               y=table_plot_n$perc_90_pca2, 
                               wt=table_plot_n$freq_plot_nurse)
    
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_x"] <- cic[1]
    centroid_df[centroid_df$nurse_code==nurse_id, "cic_y"] <- cic[3]
    
    centroid_df_plot <- centroid_df%>%
      filter(nurse_code== nurse_id)
    centroid_df[centroid_df$nurse_code==nurse_id, "cd"]<- GMINdistance(centroid_df_plot$cic_x, centroid_df_plot$cic_y,
                                                                       centroid_df_plot$oc_x, centroid_df_plot$oc_y)[1]
    
    
    table_plot_n$species_code <- factor(table_plot_n$species_code)
    cols <- colorRampPalette(brewer.pal(12, "Spectral"))
    mypal <- cols(length(unique(table_plot_n$species_code)))
    
    gg_i <- ggplot(table_plot_n %>%
                     dplyr::mutate(freq_plot_nurse = case_when(
                       freq_plot_nurse == 0~NA_real_,
                       freq_plot_nurse != 0~freq_plot_nurse))) +
      geom_point( aes(x=perc_90_pca1, y=perc_90_pca2, 
                      color=species_code, size=freq_plot_nurse*10
      ), alpha =0.85)+ 
      #scale_color_manual(values = mypal)+
      scale_color_viridis_d(option="plasma")+
      labs(color = "Species")+
      guides(color=guide_legend(ncol=2)) +
      xlim(min(table_plot_n$perc_90_pca1)-2, 
           max(table_plot_n$perc_90_pca1)+2)+ 
      ylim(min(table_plot_n$perc_90_pca2)-2, 
           max(table_plot_n$perc_90_pca2)+2)+ 
      scale_size(range = c(2, 7))+
      scale_alpha_continuous("freq_plot_nurse", range=c(0.2, 0.9))+
      guides(size = FALSE,
             alpha = FALSE
      )+
      geom_point(data=centroid_df_plot, aes(x=cic_x, y=cic_y),
                 color="black", size=2, shape=19)+
      geom_point(data=centroid_df_plot, aes(x=oc_x, y=oc_y), 
                 color="black", size=2, shape=17)+
      xlab("Climatic Axis 1 (48.2%)")+
      ylab("Climatic Axis 2 (26.7%)")+
      ggtitle(nurse_id)+
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_text(size = 9), 
            legend.text  = element_text(size = 8),
            legend.key.size = unit(5, "mm"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(color="black"),
            text=element_text(family="serif"))
    
    gg_list[[j]] <- gg_i
    
  }
  
  legend <- cowplot::get_legend(gg_list[[1]])
  grid::grid.draw(legend)
  
  first <- gg_list[[1]]+
    theme(legend.position = "none")
  
  second <- gg_list[[2]]+
    theme(legend.position = "none")
  
  r <- rectGrob(gp=gpar(fill="white", col="white"))
  biplot <- gridExtra::grid.arrange(first, second, r, legend, 
                                    ncol=4, widths = unit(c(6.75,6.75,1.25,1.5), "cm"),
                                    heights = unit(7,'cm'))
  
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/all_sp/", plot_list[i], "perc90_in0.png"), 
  #        plot=biplot, width = 19, height = 8, dpi = 300, units="cm")
  # ggsave(paste0("./statistical_analyses/output/climatic_diagrams/all_sp/", plot_list[i], "perc90_in0.pdf"),
  #        plot=biplot, width = 19, height = 8, units="cm")
  # 

}

# write.table(centroid_df, "./statistical_analyses/output/disequilibrium_perc90_in0.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)




