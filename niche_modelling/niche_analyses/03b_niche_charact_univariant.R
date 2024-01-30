
#*************************************#
### Maria Angeles Perez-Navarro
### CREAF  12/2020
### univariate niche characterization 
#*************************************#

library(raster) 
library(dismo) 
library(sp)
library(tidyverse)
library(data.table)
library(ade4)
library(stringr)
library(RColorBrewer)
library(ks)
library(readr)
library(ggplot2)

# load functions
source('./niche_modelling/niche_analyses/niche_required_functions.R')

# 1. Load clean occurrences and climate  -------------------------------

#**** loading species clean occurrences with climate values  ****

# dir <- "./niche_modelling/occurrences/clean/"
# 
# list_sp <- list.files(path=dir,  pattern="*_clean.csv$", full.names=TRUE)
# list_sp <- sub('.*./niche_modelling/occurrences/clean/', '', list_sp)
# species_list <- str_remove(list_sp, "_clean.csv")

# table1 <- read_delim(paste0("./niche_modelling/occurrences/clean/",
#                             species_list[1], "_clean.csv"), delim=";", col_names=T)
# 
# for (i in 2:length(species_list)){
#   table <- read_delim(paste0("./niche_modelling/occurrences/clean/",
#                              species_list[i], "_clean.csv"), delim=";", col_names=T)
#   table1 <- rbind(table, table1)
# }
# 
# table_sp <- table1%>%
#   dplyr::filter(!is.na(bio15))
# rm(table1,table)
# 
# 
# ## join and save data from field_df the species
# 
# fwrite(table_sp, file= paste0("./niche_modelling/occurrences/clean/field_df_sp_climate.csv"),
#        sep=";", row.names=F, quote=F, dec=".")

table_sp <- read_delim("./niche_modelling/occurrences/clean/field_df_sp_climate.csv",  
                       delim=";", col_names=T)


# 2. Characterize univariate niches ---------------------------------------


perc_95_list <- list()
perc_90_list <- list()
mode_list <- list()
centroid_list <- list()
hpi_list <- list()#to save niches univariate bandwidth

sp_list <- table_sp$species%>%
  as.factor()%>%
  sort()%>%
  unique()

table_sp <- as.data.frame(table_sp)#hscv function do not work with tibble

for(i in 1:length(sp_list)){
  
  centroid_sp <- list()
  mode_sp <- list()
  perc_95_sp <- list()
  perc_90_sp <- list()
  hpi_df_sp <- list()
  
  for (j in 1:ncol(table_sp[, -c(1:3)])){#the Ecol. letters paper only include bio01, bio06 and bio12
    
    hpi_med <- hscv (x=table_sp[, j+3],  binned=T)# Warning! Hscv for matrix and hscv for vectors. In addition pilot="samse" is only applicable for matrix
    bg_med <- kde(x=table_sp[, j+3], H=hpi_med )
    bg_volume <- bg_med
    # plot(bg_volume, col="black",
    #      ylim=c(0, 3*max(bg_volume$estimate)), 
    #      main=paste0(sp_list[i],
    #                  "_", names(table_sp)[j+3]))#,  xlim=c(5950,6070), ylim=c(0.00197,0.0022)
    # 
    work_specie <- table_sp[table_sp$species== sp_list[i], j+3]# c(2:4) to plot 3d
    hpi_sp <- hscv(x=work_specie, binned=T) #, plot=TRUE
    niche_volume <- kde(x=work_specie, h=if(hpi_sp<0.5){
                                           0.5
                                          }else{
                                            hpi_sp
                                          }, gridsize = 70000)
    plot(niche_volume, col="yellowgreen",  cont=95,
         main=paste0(sp_list[i],
                     "_", names(table_sp)[j+3]))#add=T,
    niche_df0 <- data.frame(cbind(niche_volume$eval.points, niche_volume$estimate))
    names(niche_df0) <- c("x", "density")
    niche_df <- niche_df0[niche_df0$density>=niche_volume$cont[["5%"]], ]
    lines(x=niche_df$x, y=niche_df$density, col="orange")
    
    
    ## Add centroid point
    
    centroid <- wt.mean(x=niche_df$x, wt=niche_df$density^10)
    abline(v=centroid, col="orange")
    
    centroid_df <- data.frame(centroid)%>%
      mutate(species=sp_list[i],
             variable=names(table_sp)[j+3])%>%
      dplyr::select(species, variable, centroid)
    
    centroid_sp[[j]] <- centroid_df
    
    # add maximum density point
    
    mode <- niche_df%>%
      dplyr::filter(density==max(density))
    abline(v=mode$x, col="red")
    
    mode_df <- mode%>%
      # mutate(!!paste0(names(table_sp)[j+3]):= mode$x)
      mutate(species=sp_list[i],
             variable=names(table_sp)[j+3])%>%
      rename(mode=x)%>%
      dplyr::select(species, variable, mode)
    
    mode_sp[[j]] <- mode_df
    
    # add 95 percentile
    
    step1 <- niche_df0[2, "x"]-niche_df0[1, "x"]
    
    perc_95 <- niche_df%>%
      dplyr::filter(density>=niche_volume$cont[["95%"]])
    
    gap_95_df <- intra.limit.univariate(perc_95, step1)
    gap_95_df <- gap_95_df%>%
      dplyr::filter(continuity==F)%>%
      select(x, density, position)
    
    perc_95_df <- perc_95%>%# identify high density areas 
      dplyr::filter(x==max(x)|
                    x==min(x))%>%
      mutate(position= case_when(
        x==min(x)~"start",
        x==max(x)~"end"
      ))
    
    perc_95_df <- rbind(perc_95_df, gap_95_df)
      
    perc_95_df <- perc_95_df%>%
      mutate(species=sp_list[i],
             variable=names(table_sp)[j+3])%>%
      rename(perc_95=x)%>%
      dplyr::select(species, variable,
                    perc_95, position) %>% 
      dplyr::arrange(perc_95)
      
      
    abline(v=perc_95_df$perc_95, col="darkblue")
    abline(h=niche_volume$cont[["95%"]], col="red")
    
    perc_95_sp[[j]] <- perc_95_df
    
    # add 90 percentile
    
    step2 <- niche_df0[2, "x"]-niche_df0[1, "x"]
    
    perc_90 <- niche_df%>%
      dplyr::filter(density>=niche_volume$cont[["90%"]])
    
    gap_90_df <- intra.limit.univariate(perc_90, step2)
    gap_90_df <- gap_90_df%>%
      dplyr::filter(continuity==F)%>%
      select(x, density, position)
    
    perc_90_df <- perc_90%>%#identify high density areas
      dplyr::filter(x==max(x)|
                    x==min(x)) %>%
      mutate(position= case_when(
                      x==min(x)~"start",
                      x==max(x)~"end"
                    ))
    
    perc_90_df <- rbind(perc_90_df, gap_90_df)
    
    perc_90_df <- perc_90_df%>%
      mutate(species=sp_list[i],
             variable=names(table_sp)[j+3])%>%
      rename(perc_90=x)%>%
      dplyr::select(species, variable, 
                    perc_90, position) %>% 
      dplyr::arrange(perc_90)
    
    abline(v=perc_90_df$perc_90, col="darkblue")
    
    perc_90_sp[[j]] <- perc_90_df
    
    hpi_df <- mode_df%>%
      dplyr::select(species, variable)%>%
      mutate(hpi=hpi_sp)
    
    hpi_df_sp[[j]] <- hpi_df
    
    
    
    }
  
  centroid_sp_all <- bind_rows(centroid_sp)
  mode_sp_all <- bind_rows(mode_sp)
  perc_95_sp_all <- bind_rows(perc_95_sp)
  perc_90_sp_all <- bind_rows(perc_90_sp)
  hpi_df_sp_all <- bind_rows(hpi_df_sp)
  
  
  centroid_list[[i]] <- centroid_sp_all
  mode_list[[i]] <- mode_sp_all
  perc_95_list[[i]] <- perc_95_sp_all
  perc_90_list[[i]] <- perc_90_sp_all
  hpi_list[[i]] <- hpi_df_sp_all

  }


centroid_species <- bind_rows(centroid_list)
mode_species <- bind_rows(mode_list)
perc_95_species <- bind_rows(perc_95_list)
perc_90_species <- bind_rows(perc_90_list)
hpi_species <- bind_rows(hpi_list)

hpi_sp_spread <- hpi_species%>%
  spread(variable, hpi)

write.table(centroid_species, "./niche_modelling/output/centroid_niche_chelsa_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
write.table(mode_species, "./niche_modelling/output/mode_centroid_niche_chelsa_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
write.table(perc_95_species, "./niche_modelling/output/perc95_niche_chelsa_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
write.table(perc_90_species, "./niche_modelling/output/perc90_niche_chelsa_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
write.table(hpi_sp_spread, "./niche_modelling/output/hpi_chelsa_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)



# 3. Get distance to population observed climate --------------------------

# can take several hours  to run all the species 

source('niche_required_functions.R')

#** load population tables ****

# field_df <- read.table("./niche_modelling/data/field_raw_data.csv", 
#                        sep=";", dec=".", h=T)
# 
# coords <- distinct(field_df[, c("community","long", "lat")])
# 
# 
# #** load climatic maps ****
# 
# dir <- "./niche_modelling/climate/chelsa_bio_v1.2/"
# list_variables <- list.files(path=dir, full.names=TRUE)
# vars_def <- brick(stack(list_variables), progress="text")
# 
# 
# ## extract climatic values for study sites
# 
# coords_climate <- raster::extract(vars_def, coords[, c("long", "lat")])
# pop_climate <- cbind(coords, coords_climate)

#write.table(pop_climate, "./niche_modelling/data/populations_climate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


## load multivariate centroid and join 

all_pop <- read.table("./niche_modelling/data/field_raw_data.csv", 
                    sep=";", dec=".", h=T)

pop_climate <- read_delim("./niche_modelling/data/populations_climate.csv", delim=";", col_names=T)

names(pop_climate)[4:22] <- paste0("pop_", names(pop_climate)[4:22])

all_pop <- all_pop%>%
  dplyr::select(community, code, lat, long,
                recruit_species, nurse, n,
                canopy_percent, gap_percent, site,
                genus, sp, sub_sp, recruit_sp)


all_climate <- left_join(all_pop, pop_climate%>%
                           select(-c(lat, long)), 
                         by="community")


# 3.1 Distance to niche centroid and mode -----------------------------------

# centroid_uni <- read_delim("./niche_modelling/output/centroid_niche_chelsa_univariate.csv", delim=";", col_names=T)
# 
# centroid_uni <- centroid_uni%>%
#   mutate(species= forcats::fct_recode(species,
#                                       'Rosmarinus officinalis' = 'Salvia rosmarinus',
#                                       'Periploca angustifolia' = 'Periploca laevigata'))%>%
#   rename(recruit_sp=species)
# 
# all_clim <- all_climate
# names(all_clim)[15:33] <- substr(names(all_clim)[15:33], 5,9)
# all_clim <- all_clim%>%
#   gather(key="variable", value="population", 
#          bio01, bio02, bio03, bio04,
#          bio05, bio06, bio07, bio08,
#          bio09, bio10, bio11, bio12,
#          bio13, bio14, bio15, bio16,
#          bio17, bio18, bio19)
# 
# setdiff(all_clim$recruit_sp, centroid_uni$recruit_sp)
# 
# all_clim_ctr <- left_join(all_clim, centroid_uni,
#                           by= c("recruit_sp", "variable"))%>%
#   mutate(dist=population-centroid)
# 
# all_clim_ctr1 <- all_clim_ctr%>%
#   pivot_wider(id_cols = c(community,lat,long,facilitated_species,
#                           nurse, n, canopy_percent,gap_percent,
#                           site,genus,sp,sub_sp,recruit_sp), 
#               names_from = variable, 
#               values_from = c("dist", "centroid","population"))

#write.table(all_clim_ctr1, "./niche_modelling/output/centroid_distances_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)



### mode distances ****

# mode_uni <- read_delim("./niche_modelling/output/mode_centroid_niche_chelsa_univariate.csv", delim=";", col_names=T)
# 
# mode_uni <- mode_uni%>%
#   mutate(species= forcats::fct_recode(species,
#                                       'Rosmarinus officinalis' = 'Salvia rosmarinus',
#                                       'Periploca angustifolia' = 'Periploca laevigata'))%>%
#   rename(recruit_sp=species)
# 
# all_clim <- all_climate
# names(all_clim)[15:33] <- substr(names(all_clim)[15:33], 5,9)
# all_clim <- all_clim%>%
#   gather(key="variable", value="population", 
#          bio01, bio02, bio03, bio04,
#          bio05, bio06, bio07, bio08,
#          bio09, bio10, bio11, bio12,
#          bio13, bio14, bio15, bio16,
#          bio17, bio18, bio19)
# 
# setdiff(all_clim$recruit_sp, mode_uni$recruit_sp)
# 
# all_clim_ctr <- left_join(all_clim, mode_uni,
#                           by= c("recruit_sp", "variable"))%>%
#   mutate(dist=population-mode)
# 
# all_clim_ctr1 <- all_clim_ctr%>%
#   pivot_wider(id_cols = c(community,lat,long,facilitated_species,
#                           nurse, n, canopy_percent,gap_percent,
#                           site,genus,sp,sub_sp,recruit_sp), 
#               names_from = variable, 
#               values_from = c("dist", "mode","population"))

#write.table(all_clim_ctr1, "./niche_modelling/output/mode_centroid_distances_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


# 3.2 Distance to closest point from niche discarding bottom 95 density percentile ------------

perc_95 <- read_delim("./niche_modelling/output/perc95_niche_chelsa_univariate.csv", delim=";", col_names=T)

perc_95 <- perc_95%>%
  mutate(species= forcats::fct_recode(species,
                                      'Rosmarinus officinalis' = 'Salvia rosmarinus',
                                      'Periploca angustifolia' = 'Periploca laevigata'))

all_clim <- all_climate
setdiff(all_clim$recruit_sp, perc_95$species)


for(i in 1:nrow(all_clim)){
  
  specie <- all_clim[i, "recruit_sp"]
  
  print(paste(specie,i))
  
  for(k in 1:19){
    
    print(k)
    
    biovariable <- if(k<=9){
      paste0("bio0", k )}else{
        paste0("bio", k )}
    
    perc_95_spb0 <- perc_95%>%
      dplyr::filter(species %in% specie &
             variable == biovariable)
    
    if(nrow(perc_95_spb0)==3){
      perc_95_spb0 <- perc_95_spb0[c(1,3), ]
    }
    
    perc_95_spb <- perc_95_spb0%>%
      dplyr::filter(position%in%c("start","end"))# not to consider peaks in interval in_out analyses
   
    if(nrow(perc_95_spb)>1){
    perc_95_spb <- perc_95_spb%>%
      mutate(auto= rep(c("start", "end"), 
                       nrow(perc_95_spb)/2))
    vector_num <- sort(rep(seq(1:12),2))
    
    if(unique(perc_95_spb$position==perc_95_spb$auto)){
      perc_95_spb <- perc_95_spb%>%
        mutate(position_spe= paste0(position, 
                                    vector_num[1:nrow(perc_95_spb)]))
    }else{    
      stop("Error in percentil vector order")
    }
    
    }
    
    
    
    # estimate if population is inside or outside percentile 
    # I repeated the condition as many times as potential discontinuities may be
    
    if(nrow(perc_95_spb)==1){
      in_out <- "out"
    }
    
    if(nrow(perc_95_spb)==2){
      in_out <- if (all_clim[i, k+14] > perc_95_spb[perc_95_spb$position_spe=="start1", "perc_95" ] &
                    all_clim[i, k+14] < perc_95_spb[perc_95_spb$position_spe=="end1", "perc_95" ]){
        "in"
      }else{
        "out"}
    } 
    
    if(nrow(perc_95_spb)==4){
      in_out <- if (all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start1", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end1", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start2", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end2", "perc_95" ]){
        "in"
      }else{"out"}
    }
    
    if(nrow(perc_95_spb)==6){
      in_out <- if (all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start1", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end1", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start2", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end2", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start3", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end3", "perc_95" ]){
        "in"
      }else{"out"}
    } 
    
    if(nrow(perc_95_spb)==8){
      in_out <- if (all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start1", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end1", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start2", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end2", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start3", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end3", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start4", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end4", "perc_95" ]){
        "in"
      }else{"out"}
    } 
    
    if(nrow(perc_95_spb)==10){
      in_out <- if (all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start1", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end1", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start2", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end2", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start3", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end3", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start4", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end4", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start5", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end5", "perc_95" ]){
        "in"
      }else{"out"}
    } 
    
    if(nrow(perc_95_spb)==12){
      in_out <- if (all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start1", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end1", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start2", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end2", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start3", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end3", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start4", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end4", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start5", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end5", "perc_95" ]|
                    all_clim[i, k+14]> perc_95_spb[perc_95_spb$position_spe=="start6", "perc_95" ]&
                    all_clim[i, k+14]< perc_95_spb[perc_95_spb$position_spe=="end6", "perc_95" ]){
        "in"
      }else{"out"}
    }
 
    # estimate distance
    distance <- GMINdistance.univariate(perc_95_spb0$perc_95, all_clim[i, k+14])
    
    # add sign depending on position
    dist <- if(all_clim[i, k+14]>distance$closest.x){
      distance$min.dist
    }else{-distance$min.dist}
    
    # add dist and closest perc95 to the table
    
    if(in_out=="in"){
      all_clim[i, paste0("closest_p95_", biovariable )] <- all_clim[i, paste0("pop_", biovariable )]
    }else{all_clim[i, paste0("closest_p95_", biovariable )] <- distance$closest.x}
    
    
    if(in_out=="in"){
      all_clim[i, paste0("dist_", biovariable )] <- 0
      }else{all_clim[i, paste0("dist_", biovariable )] <- dist}
                                                   
    
    }
  
  
}


#write.table(all_clim, "./niche_modelling/output/perc95_distances_in0_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)




# 3.3. Distance to closest point from niche discarding bottom 90 density percentile --------

# perc_90 <- read_delim("./niche_modelling/output/perc90_niche_chelsa_univariate.csv", delim=";", col_names=T)
# 
# perc_90 <- perc_90%>%
#   mutate(species= forcats::fct_recode(species,
#                                       'Rosmarinus officinalis' = 'Salvia rosmarinus',
#                                       'Periploca angustifolia' = 'Periploca laevigata'))
# 
# all_clim <- all_climate
# setdiff(all_clim$recruit_sp, perc_90$species)
# 
# 
# for(i in 1:nrow(all_clim)){
#   
#   specie <- all_clim[i, "recruit_sp"]
#   
#   print(paste(specie,i))
#   
#   for(k in 1:19){
#     
#     print(k)
#     
#     biovariable <- if(k<=9){
#       paste0("bio0", k )}else{
#         paste0("bio", k )}
#     
#     perc_90_spb0 <- perc_90%>%
#       dplyr::filter(species %in% specie &
#                       variable == biovariable)
#     
#     if(nrow(perc_90_spb0)==3){
#       perc_90_spb0 <- perc_90_spb0[c(1,3), ]
#     }
#     
#     perc_90_spb <- perc_90_spb0%>%
#       dplyr::filter(position%in%c("start","end"))# not to consider peaks in interval in_out analyses
#     
#     if(nrow(perc_90_spb)>1){
#       perc_90_spb <- perc_90_spb%>%
#         mutate(auto= rep(c("start", "end"), 
#                          nrow(perc_90_spb)/2))
#       vector_num <- sort(rep(seq(1:12),2))
#       
#       if(unique(perc_90_spb$position==perc_90_spb$auto)){
#         perc_90_spb <- perc_90_spb%>%
#           mutate(position_spe= paste0(position, 
#                                       vector_num[1:nrow(perc_90_spb)]))
#       }else{    
#         stop("Error in percentil vector order")
#       }
#       
#     }
#     
#     
#     
#     # estimate if population is inside or outside percentile 
#     # I repeated the condition as many times as potential discontinuities may be
#     
#     if(nrow(perc_90_spb)==1){
#       in_out <- "out"
#     }
#     
#     if(nrow(perc_90_spb)==2){
#       in_out <- if (all_clim[i, k+14] > perc_90_spb[perc_90_spb$position_spe=="start1", "perc_90" ] &
#                     all_clim[i, k+14] < perc_90_spb[perc_90_spb$position_spe=="end1", "perc_90" ]){
#         "in"
#       }else{
#         "out"}
#     } 
#     
#     if(nrow(perc_90_spb)==4){
#       in_out <- if (all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start1", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end1", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start2", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end2", "perc_90" ]){
#         "in"
#       }else{"out"}
#     }
#     
#     if(nrow(perc_90_spb)==6){
#       in_out <- if (all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start1", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end1", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start2", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end2", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start3", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end3", "perc_90" ]){
#         "in"
#       }else{"out"}
#     } 
#     
#     if(nrow(perc_90_spb)==8){
#       in_out <- if (all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start1", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end1", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start2", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end2", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start3", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end3", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start4", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end4", "perc_90" ]){
#         "in"
#       }else{"out"}
#     } 
#     
#     if(nrow(perc_90_spb)==10){
#       in_out <- if (all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start1", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end1", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start2", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end2", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start3", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end3", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start4", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end4", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start5", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end5", "perc_90" ]){
#         "in"
#       }else{"out"}
#     } 
#     
#     if(nrow(perc_90_spb)==12){
#       in_out <- if (all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start1", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end1", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start2", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end2", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start3", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end3", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start4", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end4", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start5", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end5", "perc_90" ]|
#                     all_clim[i, k+14]> perc_90_spb[perc_90_spb$position_spe=="start6", "perc_90" ]&
#                     all_clim[i, k+14]< perc_90_spb[perc_90_spb$position_spe=="end6", "perc_90" ]){
#         "in"
#       }else{"out"}
#     }
#     
#     # estimate distance
#     distance <- GMINdistance.univariate(perc_90_spb0$perc_90, all_clim[i, k+14])
#     
#     # add sign depending on position
#     dist <- if(all_clim[i, k+14]>distance$closest.x){
#       distance$min.dist
#     }else{-distance$min.dist}
#     
#     # add dist and closest perc90 to the table
#     
#     if(in_out=="in"){
#       all_clim[i, paste0("closest_p90_", biovariable )] <- all_clim[i, paste0("pop_", biovariable )]
#     }else{all_clim[i, paste0("closest_p90_", biovariable )] <- distance$closest.x}
#     
#     
#     if(in_out=="in"){
#       all_clim[i, paste0("dist_", biovariable )] <- 0
#     }else{all_clim[i, paste0("dist_", biovariable )] <- dist}
#     
#     
#   }
#   
#   
# }
# 

#write.table(all_clim, "./niche_modelling/output/perc90_distances_in0_univariate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


