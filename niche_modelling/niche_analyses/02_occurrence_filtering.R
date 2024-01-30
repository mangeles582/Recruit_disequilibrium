


#*******************************#
### Maria Angeles Perez-Navarro
### CREAF  09/2020
### Species filtering
#*******************************#


library(raster) 
library(dismo) 
library(sp)
library(stringr)
library(ggmap)
library(maps)
library(readr)
library(tidyverse)
#devtools::install_github("matthewkling/chelsaDL")#required package to download CHELSA v.1.2 directly from R
library(chelsaDL)

source("./niche_modelling/niche_analyses/niche_required_functions.R")#extra required functions, from deprecated packages or created by Maria A. Perez-Navarro


# 1. Load gbif and climate data -------------------------------------------------


# 1.1 load unclean occurrence data ----------------------------------------

dir <- "./niche_modelling/occurrences/unclean/" 
list_spdw <- list.files(path=dir,  pattern="*.csv$", full.names=TRUE)%>%sort()
list_spdw_name <- sub('.*./niche_modelling/occurrences/unclean/gbif_', '', list_spdw)
list_spdw_name <- sub("\\.csv.*", "", list_spdw_name)
species_list <- gsub("_", " ", list_spdw_name)

# 1.2 download and save climatic variables CHELSA v.1.2 -------------------

#* package climenv allows for diretly downloading chelsa v.2.1 
#* which was not used in this paper

# dir.create("./niche_modelling/climate/chelsa_bio_v1.2/")#create folder chelsa 79-2013_v1.2 within climate folder
# bio_folder <- "./niche_modelling/climate/chelsa_bio_v1.2/"
# 
# for(i in 1:19){
#   
#   if(i<10){ 
#     path <- noquote(paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_0", i,".tif"))
#     path_dest <- paste0(bio_folder,"CHELSA_bio10_0", i, ".tif")
#     download.file(path,destfile = path_dest, mode = "wb")
#   }else{
#     path <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_", i,".tif")
#     path_dest <- paste0(bio_folder,"CHELSA_bio10_", i,".tif")
#     download.file(path,destfile = past_dest, mode = "wb")
#   }
#  
# }
# 

#same for radiation and PET

# dir.create("./niche_modelling/climate/chelsa_radiation_v1.2/")
# rad_folder <- "./niche_modelling/climate/chelsa_radiation_v1.2/"
# 
# for(i in 1:12){
#   path <- noquote(paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/exchelsa/srad/CHELSA_srad_", i,".tif"))
#   path_dest <- paste0(rad_folder,"radiation_", i, ".tif")
#   download.file(path,destfile = path_dest, mode = "wb")
# }
# 
# list_rad <- list.files(path=rad_folder, full.names=TRUE)# igual me sobra
# rad_months <- stack(list_rad)
# avg_rad <- mean(rad_months)
# 
# writeRaster(avg_rad,paste0(rad_folder,"annual_radiation_mean.tif"),
#             format="GTiff", overwrite=TRUE)
# 
# #***
# 
# dir.create("./niche_modelling/climate/chelsa_pet_v1.2/")
# pet_folder <- "./niche_modelling/climate/chelsa_pet_v1.2/"
# 
# for(i in 1:12){
#   path <- noquote(paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/exchelsa/pet/CHELSA_pet_",i,"_1979-2013.tif"))
#   path_dest <- paste0(pet_folder,"pet_", i, ".tif")
#   download.file(path,destfile = path_dest, mode = "wb")
# }
# 
# list_pet <- list.files(path=pet_folder, full.names=TRUE)# igual me sobra
# pet_months <- stack(list_pet)
# avg_pet <- mean(pet_months)
# 
# writeRaster(avg_pet,paste0(pet_folder,"annual_pet_mean.tif"),
#             format="GTiff", overwrite=TRUE)
# 
# annual_pet <- raster("./niche_modelling/climate/chelsa_bio_v1.2/bio12.tif")
# annual_ppet <- annual_pet/avg_pet
# writeRaster(annual_ppet,paste0(pet_folder,"annual_ppet_mean.tif"),
#             format="GTiff", overwrite=TRUE)


dir <- "./niche_modelling/climate/chelsa_bio_v1.2/work_extension/"
list_variables <- list.files(path=dir, full.names=TRUE)
vars_def <- stack(list_variables)


# 2.1 Filter and save species occurrences ----------------------------------

dir.create("./niche_modelling/occurrences/clean/")#create folder and save there the downloaded zip
world <- readOGR("./niche_modelling/climate/work_area_shp/world_countries",
                 "World_Countries")#also upload world shp for plots

final_data <- data.frame(matrix(ncol=3,
                                nrow=length(species_list)))
names(final_data) <- c("species", "filtering_dist", "occurrences")


#* it is possible to run it in a loop but it is better to run each species
#* one by one and carefully check record distribution to manage potential
#* bias and errors

set.seed(3)

dir.create("./niche_modelling/occurrences/clean/")
dir.create("./niche_modelling/occurrences/clean/figures/")

for (i in 1: length(species_list) ){ 
  
  print(i)
  
  skip_to_next <- FALSE
  
  # loading species tables

  df <- if(length(names(read_delim(list_spdw[i],  delim=";", col_names=T)))>8){
    read_delim(list_spdw[i],  delim=";", col_names=T)
  }else(read_delim(list_spdw[i],  delim="\t", col_names=T))
  
  table <- tryCatch(if(any(names(df)=="lon")){
            df}else if ( any(names(df)=="decimalLongitude")) {
            df%>%
            mutate(lat= as.numeric(decimalLatitude),
                   lon= as.numeric(decimalLongitude))},
            error = function(e) { skip_to_next <<- TRUE})
                      
  
  # table <- tryCatch(read_delim(list_spdw[i],  delim=";", col_names=T)%>%
  #                     mutate(lat= as.numeric(decimalLatitude),
  #                            lon= as.numeric(decimalLongitude)), 
  #                   error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
  # species id filtering
  tryCatch( unique(table$taxonKey),error = function(e) {})
  tryCatch( barplot(table(table$taxonKey)),error = function(e) {})
  sum(is.na(table$taxonKey))
  unique(table$scientificName) #if any
  barplot(table(table$scientificName))
  table$scientificName <- unique(table$scientificName)[1]
  unique(table$speciesKey)
  unique(table$species)
  
  
  # study extension filtering and loading variables
  table_onvar<-raster::extract(x=vars_def, 
                               y=table[ , c("lon","lat")])
  
  table_onvar <- data.frame(table_onvar)
  table<-data.frame(table, table_onvar)
  nrow(table)
  rm(table_onvar)
  table<-table[!(is.na(table$bio03)), ]
  unique(is.na(table$bio15))
  nrow(table)
  
  
  #plot
  plot_map_occ(world, x_ext=c(-35, 70), y_ext=c(18,60), 
               df1= table, 
               title=species_list[i]) #1black, 2 red, 3green
  
  
  ## basis of record filtering
  tryCatch(barplot(table(table$basisOfRecord)),error = function(e) {sum(is.na(table$basisOfRecord))})
  unique(table$basisOfRecord)
  nrow((table[table$basisOfRecord == "FOSSIL_SPECIMEN", ]))
  nrow((table[table$basisOfRecord == "MATERIAL_SAMPLE", ]))
  table<-table[table$basisOfRecord != "FOSSIL_SPECIMEN", ]
  table<-table[table$basisOfRecord != "MATERIAL_SAMPLE", ]
  
  ## year filtering
  table$year <- as.numeric(table$year)
  tryCatch(hist(table$year),error = function(e) {sum(is.na(table$year))})
  sum(is.na(table$year))
  
  plot_map_occ(world, x_ext=c(-35, 70), y_ext=c(18,60), 
               df1= table, 
               df2= table[table$year>= 1970, ],
               df3= table[(is.na(table$year)), ],
               title=species_list[i]) #1black, 2 red, 3green
  
  #table<-table[!is.na(table$year), ] # if we decide to remove without year data
  #table1<-table[table$year >= 1970, ] # if we want only year above 1970
  #table1 <- table[-which(table$year <= 1970),] #also includes NA values
  table1 <- table%>%
    filter(is.na(year)|
            year>=1970)
  nrow(table)
  nrow(table1)
  table <- table1
  
  ## coordinate precision filtering
  # table$coordinatePrecision<-as.numeric(table$coordinatePrecision)
  # table$coordinateUncertaintyInMeters<-as.numeric(table$coordinateUncertaintyInMeters)
  # unique(table$coordinatePrecision)
  # barplot(table(table$coordinatePrecision))
  # unique(table$coordinateUncertaintyInMeters)# till 1000
  # barplot(table(table$coordinateUncertaintyInMeters))# till 707
  # 
  # sum(is.na(table$coordinatePrecision)& is.na(table$coordinateUncertaintyInMeters))
  # 
  # coord_na <- table%>%
  #   dplyr::filter(is.na(coordinatePrecision)& 
  #                   is.na(coordinateUncertaintyInMeters))
  # 
  # nrow(coord_na)
  # 
  # coord_no_na <- table%>%
  #   dplyr::filter(!is.na(coordinatePrecision)| 
  #                 !is.na(coordinateUncertaintyInMeters))
  # 
  # nrow(coord_no_na)
  # 
  # plot_map_occ(world, x_ext=c(-35, 70), y_ext=c(18,60), 
  #              df1= table, 
  #              df2= coord_na,
  #              df3= coord_no_na,
  #              title=species_list[i]) #1black, 2 red, 3green
  # 
  # # table <- coord_no_na
  # 
  # length(table$coordinatePrecision[table$coordinatePrecision<=707])
  # length(table$coordinateUncertaintyInMeters[table$coordinateUncertaintyInMeters<= 1000])
  # 
  # 
  # plot_map_occ(world, x_ext=c(-35, 70), y_ext=c(18,60), 
  #              df1=table, 
  #              df2=table[table$coordinateUncertaintyInMeters<=1000,],
  #              df3= table[table$coordinatePrecision<=707, ],
  #              title=species_list[i]) #1black, 2 red, 3green
  # 
  # 
  # table1<- table%>%
  #   dplyr::filter(coordinatePrecision<=707|
  #                   coordinateUncertaintyInMeters <= 1000)#remove coordinate with poorer precision
  
  #table<-table[table$coordinate_precision <= 7071, ] # if we decide to remove coordinate with poorer precision
  #table <- table1 
  
  # remove duplicated
  
  dup<-duplicated(table[ , c("lat", "lon")])
  length(dup[dup==TRUE])
  length(dup[!dup==TRUE])
 
  table_no_dup<-table[!dup, ]
  table <- table_no_dup
  
  # spatial autocorrelation filtering
  lat_colnumber <- which(colnames(table)=="lat")
  lon_colnumber <- which(colnames(table)=="lon")
  bio_colnumber <- which(colnames(table) %in% names(vars_def))
  
  min_dist <- ecospat.mantel.correlogram1(dfvar=table,colxy=c(lon_colnumber,lat_colnumber),
                                          n=300, colvar=bio_colnumber[1]:bio_colnumber[19], 
                                          max=50, nclass=300, nperm=1000)#ensure that this is the correct variables numer
  no_signif_dist <- min_dist$mgram%>%
    data.frame()%>%
    dplyr::filter(pval>=0.05 &
                  lag>=1)
  
  #select_dist <- min(no_signif_dist$lag)
  select_dist <- if(min(no_signif_dist$lag)>=1){
    min(no_signif_dist$lag)
  } else{
    1
  }
  
  xres(vars_def)
  res_degrees <-xres(vars_def)
  empty_cells  <- select_dist #number of cells between presence points In case of nrow()<500 ?consider to choose dist = 1 to avoid reducing too much the number of occurrences
  min_distance <- res_degrees*empty_cells
  min_distance*111.19# to know kilometre correspondence
  
  
  table1 <- ReduceSpatialClustering(data=table,
                                    minimum.distance=min_distance)#apply function to reduce spatial density
  nrow(table1)
  nrow(table)

  plot_map_occ(world, x_ext=c(-35, 70), y_ext=c(18,60),
                         df1=table, df2=table1, title=species_list[i])

  table <- table1
  rm(table1)

  ggocc <- plot_map_occ( world, x_ext=c(-35, 70), y_ext=c(18,60),
                         df1=table, title=species_list[i],
                         point_size=0.0005, title_size=10,
                         axis_title_size = 7, axis_text_size= 4.5)
  
  # clean table
  def <- table%>%
    dplyr::select("species", "lon", "lat",
           names(vars_def))%>%
    rename(x=lon,
           y=lat)
  
  head(def)
  write.table(def, file= paste0("./niche_modelling/occurrences/clean/", species_list[i], "_clean.csv"),
              sep=";", row.names=F, quote=F, dec=".")

  # ggsave(paste0("./niche_modelling/occurrences/clean/figures/", species_list[i], ".pdf"),
  #        width = 9, height = 6,  units = "cm", plot=ggocc)
  # ggsave(paste0("./niche_modelling/occurrences/clean/figures/", species_list[i], ".png"),
  #        width = 9, height = 6,  units = "cm", dpi= 300, plot=ggocc)
  # 
  # To load the data again
  load("./niche_modelling/occurrences/clean/final_data.RData")
  
  final_data[i, "species"] <- species_list[i]
  final_data[i, "filtering_dist"] <- empty_cells
  final_data[i, "occurrences"] <- nrow(def)
  
  save(final_data, file = "./niche_modelling/occurrences/clean/final_data.RData")
  
  
}


#final_data <- read_delim("./niche_modelling/occurrences/clean/zz_occurrence_summary.csv",  delim=";", col_names=T)
write.table(final_data, file= paste0("./niche_modelling/occurrences/clean/zz_occurrence_summary.csv"), 
            sep=";", row.names=F, quote=F, dec=".")



# 2.2 remove plantations data for required tree species -------------------

#** clean plantations in O.europaea, P.nigra and P. pinaster ****

dir1 <- "./niche_modelling/occurrences/clean/" 
list_change <- list.files(path=dir1,  pattern="*.csv$", full.names=TRUE)
list_change <- list_change[c(102,120,121)]#select the specific number of these species in your folder

# Olea europaea

olea_plant <- read_delim(list_change[1],  delim=";", col_names=T)
plot_map_occ(world, x_ext=c(5, 15), y_ext=c(18,51),
             df1=olea_plant%>%
               rename(lat=y,
                      lon=x),
             df2=olea_syl%>%
               rename(lat=y,
                      lon=x), title="Olea europaea")

olea_med <- olea_plant%>%
  dplyr::filter(x>=6.5&
                y<=46&
                y>22)

olea_af <- olea_plant%>%
  dplyr::filter(y<=37.5&
                y>22)

plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
             df1=olea_med%>%
               rename(lat=y,
                      lon=x),
             df2=olea_af%>%
               rename(lat=y,
                      lon=x),
             title="Olea europaea med_af")

olea <- rbind(olea_af, olea_med)

dup<-duplicated(olea[ , c("y", "x")])
length(dup[dup==TRUE])
length(dup[!dup==TRUE])

olea_no_dup<-olea[!dup, ]
olea <- olea_no_dup

plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
             df1=olea%>%
               rename(lat=y,
                      lon=x),
               title="Olea europaea join")

lat_colnumber <- which(colnames(olea)=="y")
lon_colnumber <- which(colnames(olea)=="x")
bio_colnumber <- which(colnames(olea) %in% names(vars_def))

min_dist <- ecospat.mantel.correlogram1(dfvar=olea,colxy=c(lon_colnumber,lat_colnumber),
                                        n=300, colvar=bio_colnumber[1]:bio_colnumber[19], 
                                        max=50, nclass=300, nperm=1000)#ensure that this is the correct variables numer
no_signif_dist <- min_dist$mgram%>%
  data.frame()%>%
  dplyr::filter(pval>=0.05 &
                  lag>=1)

select_dist <- min(no_signif_dist$lag)
xres(vars_def)
res_degrees <-xres(vars_def)
empty_cells  <- select_dist #number of cells between presence points In case of nrow()<500 ?look for limit select dist = 1
min_distance <- res_degrees*empty_cells
min_distance*111.19# to know kilometre correspondence


table1 <- ReduceSpatialClustering(data=olea%>%
                                    rename(lat=y,
                                           lon=x),
                                  minimum.distance=min_distance)#apply function to reduce spatial density
nrow(table1)
head(table1)
olea <- table1%>%
  rename(x=lon,
         y=lat)

ggocc <- plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
             df1=olea%>%
               rename(lat=y,
                      lon=x),
             title="Olea europaea join")

write.table(olea, file= paste0("./niche_modelling/occurrences/clean/Olea europaea_clean.csv"),
            sep=";", row.names=F, quote=F, dec=".")
ggsave(paste0("./niche_modelling/occurrences/clean/figures/Olea europaea.pdf"),
       width = 9, height = 6,  units = "cm", plot=ggocc)
ggsave(paste0("./niche_modelling/occurrences/clean/figures/Olea europaea.png"),
       width = 9, height = 6,  units = "cm", dpi= 300, plot=ggocc)


# pinus nigra

nigra_plant <- read_delim(list_change[3],  delim=";", col_names=T)

plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
             df1=nigra_plant%>%
               rename(lat=y,
                      lon=x),
             title="Pinus nigra")

nigra_filt <- nigra_plant%>%
  dplyr::filter(x>20 &
                y<50|
                x<20 &
                y<47)%>%
  dplyr::filter(x>8)
nigra_filt2 <- nigra_plant%>%
  dplyr::filter(x<5&
                y<44)

plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
             df1=nigra_filt%>%
               rename(lat=y,
                      lon=x),
             df2=nigra_filt2%>%
               rename(lat=y,
                      lon=x),
             title="Pinus nigra filt")

nigra <- rbind(nigra_filt, nigra_filt2)

dup<-duplicated(nigra[ , c("y", "x")])
length(dup[dup==TRUE])
length(dup[!dup==TRUE])

nigra_no_dup<-nigra[!dup, ]
nigra <- nigra_no_dup


ggocc <- plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
                      df1=nigra%>%
                        rename(lat=y,
                               lon=x),
                      title="Pinus nigra")

write.table(nigra, file= paste0("./niche_modelling/occurrences/clean/Pinus nigra_clean.csv"),
            sep=";", row.names=F, quote=F, dec=".")
ggsave(paste0("./niche_modelling/occurrences/clean/figures/Pinus nigra.pdf"),
       width = 9, height = 6,  units = "cm", plot=ggocc)
ggsave(paste0("./niche_modelling/occurrences/clean/figures/Pinus nigra.png"),
       width = 9, height = 6,  units = "cm", dpi= 300, plot=ggocc)


# pinus pinaster

pinaster_plant <- read_delim(list_change[4],  delim=";", col_names=T)

plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
             df1=pinaster_plant%>%
               rename(lat=y,
                      lon=x),
             title="Pinus pinaster")

pinaster_filt <- pinaster_plant%>%
  dplyr::filter(x<0&
                y<46|
                x>=0&
                y<43.5)

plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
             df1=pinaster_filt%>%
               rename(lat=y,
                      lon=x),
             title="Pinus pinaster filt")

pinaster <- pinaster_filt

dup<-duplicated(pinaster[ , c("y", "x")])
length(dup[dup==TRUE])
length(dup[!dup==TRUE])

pinaster_no_dup<-pinaster[!dup, ]
pinaster <- pinaster_no_dup


ggocc <- plot_map_occ(world, c(-35, 80), y_ext=c(18,60),
                      df1=pinaster%>%
                        rename(lat=y,
                               lon=x),
                      title="Pinus pinaster")

write.table(pinaster, file= paste0("./niche_modelling/occurrences/clean/Pinus pinaster_clean.csv"),
            sep=";", row.names=F, quote=F, dec=".")
ggsave(paste0("./niche_modelling/occurrences/clean/figures/Pinus pinaster.pdf"),
       width = 9, height = 6,  units = "cm", plot=ggocc)
ggsave(paste0("./niche_modelling/occurrences/clean/figures/Pinus pinaster.png"),
       width = 9, height = 6,  units = "cm", dpi= 300, plot=ggocc)


## join with rest species names ####

final_data <- read_delim("./niche_modelling/occurrences/clean/zz_occurrence_summary.csv",  delim=";", col_names=T)
final_data <- final_data%>%
  mutate(taxon_name= gsub("_.*","",species))

species_names <- read.table("../data/field_spcies.csv", sep=";", dec=".", h=T)


final_join <- dplyr::full_join(final_data, species_names, by="taxon_name")

write.table(final_join, file= paste0("./niche_modelling/occurrences/clean/occurrence_required_summary1.csv"), 
            sep=";", row.names=F, quote=F, dec=".")

