
#*************************************#
### Maria Angeles Perez-Navarro
### CREAF  10/2020
### multivariate niche characterization 
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


# 2. Run PCA to create a common climatic space ----------------------------

# dudi_sp <- dudi.pca(table_sp%>%
#                     dplyr::select(c(4,7:10,13:20)), 
#                     scale=T, scannf = F, nf=5)
# 
# summary(dudi_sp)
# screeplot(dudi_sp)
# dudi_sp$co
# s.corcircle(dudi_sp$co[, c(1,2)])# if you don't extract specific axis it shows the two first
# factoextra::fviz_pca_var(dudi_sp, col.var = "contrib", 
#                          gradient.cols=c("#00AFBB", "#E7B800", "#FF0000"))

#save(dudi_sp, file = "./dudi_sp.RData")
load("./niche_modelling/niche_analyses/dudi_sp.RData")

table_sp_pca <- suprow(dudi_sp, table_sp%>%
                       dplyr::select(c(4,7:10,13:20)))$li

total_sp_li <- cbind(table_sp%>%
                     dplyr::select(1),
                     table_sp_pca[, c(1,2)])#join species and predicted 2 first pca axis values
colnames(total_sp_li)[1] <- c("species")
total_sp_li <- total_sp_li[order(total_sp_li$species),]


# 3. Species niche characterization ---------------------------------------

# can take hours or days to run all the species depending on the computer

#** create common grid space where characterizing species niche ****

x <- raster(ncol=2000, nrow=2000,
            xmn=-20, xmx=20, 
            ymn=-20, ymx=20)
values(x) <- 0

sp_list <- table_sp$species%>%
  as.factor()%>%
  sort()%>%
  unique()

perc_95 <- list()
perc_90 <- list()
mode_list <- list()

species_pca <- data.frame(sp_list)%>%
  mutate(centroid_pca1=NA,
         centroid_pca2=NA,
         niche_area_suit=NA,
         niche_area_a=NA)%>%
  rename(species=sp_list)


mypal <- colorRampPalette(c("white",
                            "orange",
                            "darkorange"))(15)
myblue <- colorRampPalette(c(#"#2c7bb6",
                             #"#00a6ca",
                             "#FFFFFF",                            
                             "#00ccbc",
                             "#90eb9d",
                             "#ffff8c",
                             "#ffff8c",
                             "#f9d057",
                             "#f29e2e",
                             "#e76818",
                             "#d7191c",
                             "#990000"))(20)
scales::show_col(myblue)

dir.create(".niche_modelling/niche_analyses/output/")
dir.create(".niche_modelling/niche_analyses/output/niche_raster/")
dir.create(".niche_modelling/niche_analyses/output/niche_plot/")

for (i in 1:length(sp_list)){
  
  work_specie <- total_sp_li[total_sp_li$species== sp_list[i], c(2,3)]# c(2:4) to plot 3d
  hpi_sp <- Hscv (x=as.matrix(work_specie),  binned=T, pilot="samse")
  niche_volume <- kde(x=work_specie, H=hpi_sp, compute.cont=T )#?gridsize=250
  
  plot(niche_volume, display="filled.contour2", drawpoints=F, cont=c(5,10,25,50,75,90,95),
       col = colorRampPalette(c("white", "orange"))(8), xlab="PCA1", ylab="PCA2", 
       xlim=c(min(as.vector(total_sp_li[, "Axis1"])),
              max(as.vector(total_sp_li[, "Axis1"]))), 
       ylim=c(min(as.vector(total_sp_li[, "Axis2"])),
              max(as.vector(total_sp_li[, "Axis2"]))), 
       main=paste0(sp_list[i]))
  
  ## get niche centroid
  niche_rast <- raster(niche_volume)# do that in the complex way (niche_volume$estimate, xmx, etc) give a especular plot
  niche_rast@data@values[which(niche_rast@data@values<niche_volume$cont[95])] <- 0# here select the desired percentile
  niche_spdf <- as(niche_rast, "SpatialPixelsDataFrame")
  niche_df <- data.frame(niche_spdf@coords, niche_spdf@data)
  #niche_df[niche_df$layer< contourLevels(niche_volume, prob=0.05), "layer"] <- 0
  cograv <- COGravity(x=niche_df[, 1], y=niche_df[, 2], 
                      wt=niche_df[, 3]^50)#opposite to decimal log?
  cograv# option a centre of mass
  centroid <- data.frame(matrix(nrow=1, ncol=2))
  centroid[, 1] <- cograv[[1]]
  centroid[, 2] <- cograv[[3]]
  points(x=centroid[, c(1,2)], pch=19, cex=0.7, col="darkorange4")
  
  ## get niche mode
  mode <- mode_value(niche_df, "layer")
  points(x=mode[, c(1,2)], pch=1, cex=0.85, col="black")
  
  ## raster suitability
  niche_suit <- niche_rast/max(values(niche_rast))
  
  ## add centroid and mode coordinates
  species_pca[species_pca$species==sp_list[i], "centroid_pca1"] <- centroid[, 1]
  species_pca[species_pca$species==sp_list[i], "centroid_pca2"] <- centroid[, 2]
  mode_df <- mode%>%
    mutate(species=sp_list[i])%>%
    rename(mode_pca1= x,
           mode_pca2= y,
           max_dens = layer)%>%
    dplyr::select(species, mode_pca1, mode_pca2, max_dens)# I do that instead of repeating previous two code lines in case of some species present more than one mode
  
  mode_list[[i]] <- mode_df
  
  ## get niche limit coordinates 5 percentile
  niche_rast1 <- niche_rast
  niche_rast1@data@values[which(is.na(niche_rast1@data@values))] <- 0
  niche_rast1@data@values[niche_rast1@data@values>0] <- 1 # obtain raster with continuous values to get only one contour
  a <- rasterToContour(niche_rast1)
  b <-  as(a, "SpatialPointsDataFrame")# alternatively we can use b <- (a[1,]@lines[[1]]@Lines[[1]]) but it has 10 time less coordinates
  b <- data.frame(b@coords)
  names(b) <- c("pca1", "pca2")
  points(b$pca1, b$pca2, col="blue", pch=19, cex=0.7)
  
  ## get niche limit coordinates 95 percentile
  niche_rast2 <- niche_rast
  niche_rast2@data@values[which(niche_rast2@data@values<niche_volume$cont[5])] <- 0# here select the desired percentile
  niche_rast2@data@values[which(is.na(niche_rast2@data@values))] <- 0
  niche_rast2@data@values[niche_rast2@data@values>0] <- 1 # obtain raster with continuous values to get only one contour
  a2 <- rasterToContour(niche_rast2)
  b2 <-  as(a2, "SpatialPointsDataFrame")# alternatively we can use b <- (a[1,]@lines[[1]]@Lines[[1]]) but it has 10 time less coordinates
  b2 <- data.frame(b2@coords)
  names(b2) <- c("pca1", "pca2")
  points(b2$pca1, b2$pca2, col="red", pch=19, cex=0.7)
  b2$species <- sp_list[i]
  b2 <- b2[, c("species", "pca1", "pca2")]
  perc_95[[i]] <- b2
  
  ## get niche limit coordinates 90 percentile
  niche_rast3 <- niche_rast
  niche_rast3@data@values[which(niche_rast3@data@values<niche_volume$cont[10])] <- 0# here select the desired percentile
  niche_rast3@data@values[which(is.na(niche_rast3@data@values))] <- 0
  niche_rast3@data@values[niche_rast3@data@values>0] <- 1 # obtain raster with continuous values to get only one contour
  a3 <- rasterToContour(niche_rast3)
  b3 <-  as(a3, "SpatialPointsDataFrame")# alternatively we can use b <- (a[1,]@lines[[1]]@Lines[[1]]) but it has 10 time less coordinates
  b3 <- data.frame(b3@coords)
  names(b3) <- c("pca1", "pca2")
  points(b3$pca1, b3$pca2, col="green", pch=19, cex=0.7)
  b3$species <- sp_list[i]
  b3 <- b3[, c("species", "pca1", "pca2")]
  perc_90[[i]] <- b3
  
  ## save raster
  a <- resample(niche_suit, x)
  unique(is.na(a@data@values))
  a@data@values[which(is.na(a@data@values))] <- 0
  
  # writeRaster(a, paste0(".niche_modelling/niche_analyses/output/niche_raster/", sp_list[i], ".tif"), 
  #             format="GTiff", overwrite=TRUE)
  # writeRaster(niche_rast2, paste0(".niche_modelling/niche_analyses/output/niche_raster/", sp_list[i], "perc_95.tif"), 
  #             format="GTiff", overwrite=TRUE)
  # writeRaster(niche_rast3, paste0(".niche_modelling/niche_analyses/output/niche_raster/", sp_list[i], "perc_90.tif"), 
  #             format="GTiff", overwrite=TRUE)
   
  ## add png plot 
  niche_suit_spdf <- as(a, 'SpatialPixelsDataFrame')
  niche_suit_df <- data.frame(niche_suit_spdf)
  

  suitsn <- c(0, 0.05, seq(0.075, 1, by=0.075))
  
  (nicheplot <- ggplot() +
    # stat_contour(data=niche_suit_df, 
    #              aes(x = x, y = y, z = layer, fill=..level..), 
    #              geom="polygon", col="transparent",
    #              breaks=suitsn
    # ) +
    geom_raster(data=niche_suit_df, 
                aes(x = x, y = y, fill=layer))+
    #scale_fill_gradient(low = mypal[2], high = mypal[15])+
    scale_fill_gradientn(colours = myblue)+
    geom_contour(data=niche_suit_df,
                 aes(x = x, y = y, z = layer) ,
                 colour="grey45",
                 # breaks = suitsn,
                 breaks=c(0.1, 0.5, 0.75),#0.05
                 linetype="longdash",
                 size=0.22,
                 alpha=0.5)+
    geom_contour(data=niche_suit_df,
                 aes(x = x, y = y, z = layer) ,
                 colour="black",
                 breaks=c(0.90, 0.95),
                 linetype="longdash",
                 size=0.22,
                 alpha=0.5)+
    geom_point(aes(x=centroid$X1, y=centroid$X2), size=0.3, col="grey40")+
    geom_point(aes(x=mode$x, y=mode$y), size=0.3, col="black")+
    xlim(-10,10)+
    ylim(-10,10)+  
    ggtitle(sp_list[i])+
    theme_bw()+
    ylab ( expression("Environmental Axis 2 (26.7%)")) +
    xlab ( expression("Environmental Axis 1 (48.2%)")) +
    labs(fill="Suitability") +
    theme(#legend.position="none",
          legend.background = element_blank(),
          legend.title = element_text(size = 8), 
          legend.text  = element_text(size = 7),
          legend.key.size = unit(0.5, "lines"),
          axis.text=element_text(color="black",
                                 size=7),
          axis.title=element_text(color="black",
                                 size=9),
          plot.title = element_text(color="black",
                                    size=12, face="bold", 
                                    family="serif",
                                    hjust=0.5),
          text=element_text(family="serif"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
   )
  
  # ggsave(paste0("./niche_modelling/niche_analyses/output/niche_plot/", sp_list[i], ".png"), 
  #        plot=nicheplot, width = 12, height = 8, dpi = 300, units="cm")
  # ggsave(paste0("./niche_modelling/niche_analyses/output/niche_plot/", sp_list[i], ".pdf"),
  #        plot=nicheplot, width = 12, height = 8, units="cm")
  
  
  
}

res(a)

perc_95_df <- bind_rows(perc_95)
perc_90_df <- bind_rows(perc_90)

species_mode <- bind_rows(mode_list)
nrow(species_mode)# in case it is equivalent to nrow(species_pca)
species_pca <- left_join(species_mode, species_pca, by="species")# full_join could work also

# write.table(species_pca, "./niche_modelling/niche_analyses/output/total_niche_chelsa_noseasonality.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
# write.table(perc_95_df, "./niche_modelling/niche_analyses/output/perc95_niche_chelsa_noseasonality.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
# write.table(perc_90_df, "./niche_modelling/niche_analyses/output/perc90_niche_chelsa_noseasonality.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
# 


# 4. Get distances to population observed climate -------------------------

source('niche_required_functions.R')

#** load population tables ****

field_df <- read.table("./niche_modelling/data/field_raw_data.csv", 
                    sep=";", dec=".", h=T)

coords <- distinct(field_df[, c("community","long", "lat")])


#** load climatic maps ****

dir <- "./niche_modelling/climate/chelsa_bio_v1.2/"
list_variables <- list.files(path=dir, full.names=TRUE)
vars_def <- brick(stack(list_variables), progress="text")


## extract climatic values for study sites

coords_climate <- raster::extract(vars_def, coords[, c("long", "lat")])
pop_climate <- cbind(coords, coords_climate)

#write.table(pop_climate, "./niche_modelling/data/populations_climate.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


## translate into pca space

load("./niche_modelling/niche_analyses/dudi_sp.RData")

pop_pca <- suprow(dudi_sp, pop_climate%>%
                       dplyr::select(c(4,7:10,13:20)))$li
pop_pca_join <- cbind( pop_climate%>%
                       dplyr::select(1:3), pop_pca)

#write.table(pop_pca_join, "./niche_modelling/data/populations_pca.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


## prepare table with accurate format

field_df_gather <- field_df%>%
  dplyr::select(community, lat, long,
                recruit_species, open_obs,
                canopy_obs, open_exp,
                canopy_exp, site, code)%>%
  rename(gap_percent=open_exp,
         canopy_percent=canopy_exp)%>%
  tidyr::gather(key="facilit", value="n", 
                open_obs, canopy_obs)%>%
  mutate(nurse= case_when(
    facilit=="open_obs" ~ "open",
    facilit=="canopy_obs" ~ "facil"
  ))%>%
  dplyr::select(community, code, lat, long,
                recruit_species, nurse, n,
                canopy_percent, gap_percent, site)

field_df_gather_pca <- full_join(field_df_gather, pop_pca_join%>%
                              dplyr::select(-c(long, lat, Axis4, Axis5)), 
                            by=c("community"))%>%
  mutate(genus = word(recruit_species, 1, sep = fixed("_")),
         sp = word(recruit_species, 2, sep = fixed("_")),
         sub_sp= word(recruit_species, 3, sep = fixed("_")))%>%
  mutate(recruit_sp= paste0(genus, " ", sp))



# 4.1. Distance to niche centroid and mode ------------------------------------


# load centroid table

# centroids <- read_delim("./niche_modelling/output/total_niche_chelsa_noseasonality.csv", delim=";", col_names=T)

# check differences in species names
# setdiff (field_df_gather_pca$recruit_sp, centroids$species)

# calculate distance to niche centroids

# pop <- field_df_gather_pca
# species_list <- unique(pop$recruit_sp)
# 
# for(i in 1:length(species_list)){
#   
#   pop_sp <- pop%>%
#     dplyr::filter(recruit_sp==species_list[i])
#     
#   centroid <- centroids%>%
#     dplyr::filter(species%in%species_list[i])%>%
#     dplyr::select(centroid_pca1, centroid_pca2,
#                   mode_pca1, mode_pca2)
#   pop[pop$recruit_sp==species_list[i], "centroid_pca1"] <- centroid[, 1]
#   pop[pop$recruit_sp==species_list[i], "centroid_pca2"] <- centroid[, 2]
#   
#   pop[pop$recruit_sp==species_list[i], "mode_pca1"] <- centroid[, 3]
#   pop[pop$recruit_sp==species_list[i], "mode_pca2"] <- centroid[, 4]
#   
#   distance2d_centroid <- GMINdistance(centroid[, 1], centroid[, 2], pop_sp$Axis1, pop_sp$Axis2)
#   distance2d_mode <- GMINdistance(centroid[, 3], centroid[, 4], pop_sp$Axis1, pop_sp$Axis2)
#   
#   pop[pop$recruit_sp==species_list[i], "centroid_dist2d"] <- distance2d_centroid$min.dist
#   pop[pop$recruit_sp==species_list[i], "mode_dist2d"] <- distance2d_mode$min.dist
#   
# }

#write.table(pop, "./niche_modelling/output/centroid_distances.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)

# 4.2. Distance to closest point from niche discarding bottom 95 density percentile --------

# load 95 percentile table

perc_95 <- read_delim("./niche_modelling/output/perc95_niche_chelsa_noseasonality.csv", delim=";", col_names=T)
perc_95 <- perc_95%>%
  mutate(species= forcats::fct_recode(species,
                  'Rosmarinus officinalis' = 'Salvia rosmarinus',
                  'Sorbus aria' = 'Sorbus edulis',
                  'Periploca angustifolia' = 'Periploca laevigata'))

# check differences in species names
setdiff (field_df_gather_pca$recruit_sp, perc_95$species)

pop <- field_df_gather_pca
species_list <- unique(pop$recruit_sp)

for(i in 1:length(species_list)){
 
  pop_sp <- pop%>%
    dplyr::filter(recruit_sp==species_list[i])
  
  perc_95_sp <- perc_95%>%
    dplyr::filter(species%in%species_list[i])
  
  sp_niche <- raster(paste0("./niche_modelling/output/niche_raster/", species_list[i], "perc_95.tif"))
  plot(sp_niche, main=species_list[i])
  
  sp_95_poly <- SpatialPolygons(list(Polygons(list(Polygon(perc_95_sp[, -1])), ID=1)))
  class(sp_95_poly)
  plot(sp_95_poly, add=T)
  points(perc_95_sp[, c("pca1", "pca2")], col="red", cex=0.1)
  
  in_out <- raster::extract(sp_niche, pop_sp[, c("Axis1", "Axis2")])
  
  distance2d_perc_95 <- GMINdistance(perc_95_sp$pca1, perc_95_sp$pca2, 
                                      pop_sp$Axis1, pop_sp$Axis2)
  
  pop[pop$recruit_sp==species_list[i], "in_out"] <- in_out
  pop[pop$recruit_sp==species_list[i], "perc_95_distance2d"] <- distance2d_perc_95$min.dist
  pop[pop$recruit_sp==species_list[i], "perc_95_pca1_closest"] <- distance2d_perc_95$closest.x
  pop[pop$recruit_sp==species_list[i], "perc_95_pca2_closest"] <- distance2d_perc_95$closest.y

  }

pop <- pop%>%
  mutate(perc_95_dist2d=case_when(
    in_out == 1 ~ 0,#when population is located within the niche surface distance is 0
    in_out == 0 ~ perc_95_distance2d),
  perc_95_pca1 = case_when(
    in_out == 1 ~ Axis1,
    in_out == 0 ~ perc_95_pca1_closest),
  perc_95_pca2 = case_when(
    in_out == 1 ~ Axis2,
    in_out == 0 ~ perc_95_pca2_closest))%>%
  select(-c(perc_95_distance2d, 
            perc_95_pca1_closest,
            perc_95_pca2_closest,
            in_out))

# write.table(pop, "./niche_modelling/output/perc_95_distances_in0.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)


# 4.3. Distance to closest point from niche discarding bottom 90 density percentile --------

# load 90 percentile table

# perc_90 <- read_delim("./niche_modelling/output/perc90_niche_chelsa_noseasonality.csv", delim=";", col_names=T)
# perc_90 <- perc_90%>%
#   mutate(species= forcats::fct_recode(species,
#                                       'Rosmarinus officinalis' = 'Salvia rosmarinus',
#                                       'Sorbus aria' = 'Sorbus edulis',
#                                       'Periploca angustifolia' = 'Periploca laevigata'))
# 
# # check differences in species names
# setdiff (field_df_gather_pca$recruit_sp, perc_90$species)
# 
# pop <- field_df_gather_pca
# species_list <- unique(pop$recruit_sp)
# 
# for(i in 1:length(species_list)){
#   
#   pop_sp <- pop%>%
#     dplyr::filter(recruit_sp==species_list[i])
#   
#   perc_90_sp <- perc_90%>%
#     dplyr::filter(species%in%species_list[i])
#   
#   sp_niche <- raster(paste0("./niche_modelling/output/niche_raster/", species_list[i], "perc_90.tif"))
#   plot(sp_niche, main=species_list[i])
#   
#   sp_90_poly <- SpatialPolygons(list(Polygons(list(Polygon(perc_90_sp[, -1])), ID=1)))
#   class(sp_90_poly)
#   plot(sp_90_poly, add=T)
#   points(perc_90_sp[, c("pca1", "pca2")], col="red", cex=0.1)
#   
#   in_out <- raster::extract(sp_niche, pop_sp[, c("Axis1", "Axis2")])
#   
#   distance2d_perc_90 <- GMINdistance(perc_90_sp$pca1, perc_90_sp$pca2, 
#                                      pop_sp$Axis1, pop_sp$Axis2)
#   
#   pop[pop$recruit_sp==species_list[i], "in_out"] <- in_out
#   pop[pop$recruit_sp==species_list[i], "perc_90_distance2d"] <- distance2d_perc_90$min.dist
#   pop[pop$recruit_sp==species_list[i], "perc_90_pca1_closest"] <- distance2d_perc_90$closest.x
#   pop[pop$recruit_sp==species_list[i], "perc_90_pca2_closest"] <- distance2d_perc_90$closest.y
#   
# }
# 
# pop <- pop%>%
#   mutate(perc_90_dist2d=case_when(
#     in_out == 1 ~ 0, # populations within niche has 0 distance
#     in_out == 0 ~ perc_90_distance2d),
#     perc_90_pca1 = case_when(
#       in_out == 1 ~ Axis1,
#       in_out == 0 ~ perc_90_pca1_closest),
#     perc_90_pca2 = case_when(
#       in_out == 1 ~ Axis2,
#       in_out == 0 ~ perc_90_pca2_closest))%>%
#   select(-c(perc_90_distance2d, 
#             perc_90_pca1_closest,
#             perc_90_pca2_closest,
#             in_out))
# 
# write.table(pop, "./niche_modelling/output/perc90_distances_in0.csv", sep=";", dec=".", col.names=T, row.names=F, quote=F)
# 
