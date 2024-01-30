

#**********************************#
### Maria Angeles Perez-Navarro
### CREAF  09/2020
### Species download
#**********************************#

library(dismo) 
library(stringr)
library(ggmap)
library(maps)
library(dplyr)
library(tidyr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_
library(forcats)

source("./niche_modelling/niche_analyses/niche_required_functions.R")#extra required functions, from deprecated packages or created by Maria A. Perez-Navarro

# 1. prepare species list -------------------------------------------------

# also visit https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/ 

table <- read.table("./niche_modelling/data/field_raw_data.csv", 
                    sep=";", dec=".", h=T)

unique(table$recruit_species)%>%sort()

species_names <- table%>% 
  dplyr::select(recruit_species)%>%
  dplyr::distinct()

species_names <- species_names%>%
  mutate(genus = word(species_names$recruit_species, 1, sep = fixed("_")),
         sp = word(species_names$recruit_species, 2, sep = fixed("_")),
         sub_sp= word(species_names$recruit_species, 3, sep = fixed("_")))

sub_names <- species_names%>%
  filter(!is.na(sub_sp))%>%
  mutate(sub_sp = NA)

species_names <- rbind(species_names, sub_names)

species_names <- species_names%>%
  mutate(taxon_name= case_when(
    is.na(sub_sp) ~ paste0(genus, " ", sp),
    !is.na(sub_sp) ~ paste0(genus, " ", sp, " ", sub_sp)
  ))

nrow(species_names[is.na(species_names$sub_sp), ])#taxons identified at genus level

# write.table(species_names,"./niche_modelling/data/field_species.csv", sep=";",  row.names=FALSE, quote=FALSE )


# 2. Import occurrences from GBIF --------------------------------------------

#** gbif R function do not work for datasets with more than 200k occurrences ****
#** move to line 118 code **** 

# for(i in 1:nrow(species_names)){
#   
#   skip_to_next <- FALSE
#   
#   sp_name <- species_names[i, ]
#   specie_table <- tryCatch(dismo::gbif( if(is.na(sp_name$sub_sp)){
#                                     paste0(sp_name$taxon_name, "*") 
#                                     } else {
#                                       sp_name$taxon_name
#                                     }
#                   , geo=T), error = function(e) { skip_to_next <<- TRUE})
#   
#   if(skip_to_next) { next }# to skype species that give error due to datasets larger than 200,000
# 
#   maps::map(database="world", plot=T, xlim=c(-50,65), ylim=c(0,70))
#   points(specie_table$lon, specie_table$lat, col="green", cex=0.78, pch=19)
#   # if(any( colnames(specie_table) =="datasetKey")){
#   #   names(specie_table)[which(colnames(specie_table) =="datasetKey")] <- "datasetID"
#   # }
#   specie_df <- specie_table[ , c("acceptedScientificName", "acceptedTaxonKey",
#                                  "basisOfRecord", "coordinatePrecision", 
#                                  "coordinateUncertaintyInMeters", "country",
#                                  "datasetID","datasetKey",
#                                  "eventDate", "year", 
#                                  "geodeticDatum", "scientificName", 
#                                  "species", "speciesKey", "taxonomicStatus",
#                                  "lat", "lon")]
#   
#   dup<-duplicated(specie_df[ , c("lat", "lon")])
#   length(dup[!dup==TRUE])
#   specie_df<-subset(specie_df, !dup)
#   
#   points(specie_table$lon, specie_table$lat, col="red", cex=0.6, pch=19)#only to prove the concordance
#   
#   write.table(specie_df, file=paste0("./niche_modelling/occurrences/unclean/gbif_", sp_name$complete_specie, ".csv"), sep=";",  row.names=FALSE, quote=FALSE)
#   
# }
# 
# 
# dir <- "./niche_modelling/occurrences/unclean/" 
# list_spdw <- list.files(path=dir,  pattern="*.csv$", full.names=TRUE)
# list_spdw <- list_spdw[6:(length(list_spdw)-1)]
# list_spdw <- sub('.*./niche_modelling/occurrences/unclean/gbif_', '', list_spdw)
# list_spdw <- sub("\\.csv.*", "", list_spdw)
# list_spdw <- gsub("_", " ", list_spdw)
# 
# species_required <- species_names%>%
#   mutate(required = case_when(
#     taxon_name %in% list_spdw ~ "no",
#     !taxon_name %in% list_spdw ~ "yes"
#   ))


#*** occ_download() work better for large dataset ****

# https://cran.r-project.org/web/packages/rgbif/vignettes/downloads.html

# 
# species_large <- species_required%>%
#   filter (required == "yes")

user <- "" # your gbif.org username 
pwd <- "" # your gbif.org password
email <- "" # your email 

gbif_taxon_keys <- species_names%>%
  pull("taxon_name") %>% # use fewer names if you want to just test 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" ) %>% # get only accepted and matched names #& status == "ACCEPTED"
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys


res <- rgbif::occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in('hasCoordinate', TRUE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)


occ_download_meta(res)# wait until status becomes SUCCEEDED or KILLED
occ_download_meta(res) %>% gbif_citation()

dat <- occ_download_get("0061950-200613084148143")  %>% occ_download_import() #to obtain the doi

dir.create("./niche_modelling/occurrences/unclean/")#create folder and save there the downloaded zip


#check downloaded data
occurrences <- read_delim("./niche_modelling/occurrences/unclean/0061950-200613084148143.csv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE, progress = show_progress(), skip_empty_rows = TRUE)

sort(unique(occurrences$scientificName))
sort(unique(occurrences$species))
sort(unique(species_names$taxon_name))


# species_new_names <- species_names%>%
#   mutate(need_change = case_when(
#     taxon_name %in% occurrences$species ~ "no",
#     !taxon_name %in% occurrences$species ~ "yes"
#   ))%>%
#   filter(need_change=="yes")


# as whole dataset is too heavy, separate and save each species datasets

species_list <- species_names$taxon_name%>% unique() %>% sort()


for (i in 1: length(species_list)){
  
  species_df <- occurrences%>%
    filter (species == species_list[i])
  
  write.table(species_df, file=paste0("./niche_modelling/occurrences/unclean/gbif_", gsub(" ", "_", species_list[i]), ".csv"), sep=";",  row.names=FALSE, quote=FALSE)
  
}

# prepare table with number of downloaded occurrences per species


dir <- "./niche_modelling/occurrences/unclean/" 
list_spdw <- list.files(path=dir,  pattern="*.csv$", full.names=TRUE)
list_spdw <- list_spdw[8:177]
list_spdw_name <- sub('.*./niche_modelling/occurrences/unclean/gbif_', '', list_spdw)
list_spdw_name <- sub("\\.csv.*", "", list_spdw_name)
list_spdw_name <- gsub("_", " ", list_spdw_name)


table_occ <- data.frame (matrix( ncol=2, nrow=length(list_spdw_name)))
names(table_occ) <- c("taxon_name", "n_occ")


species_out <- species_names%>%
  filter(!taxon_name %in% list_spdw_name)
species_out #species not correctly download

for (i in 1:length(list_spdw_name) ){
  
  occurrences <- read_delim(list_spdw[i], ";", escape_double = FALSE, trim_ws = TRUE)
  
  # occurrences <- occurrences%>%
  #   filter(!is.na("decimalLatitude"))%>%
  #   filter(!is.na("decimalLongitude" ))
  # 
  # dup<-duplicated(occurrences[ , c("decimalLatitude", "decimalLongitude")])
  # length(dup[!dup==TRUE])
  # occurrences<-subset(occurrences, !dup)
  
  table_occ[i, "taxon_name"] <- list_spdw_name[i]
  table_occ[i, "n_occ"] <- nrow(occurrences)
  # write.table(occurrences, file=list_spdw[i], sep=";",  row.names=FALSE, quote=FALSE)
  
}

table_occ

nrow(table_occ)
species_occ <- table_occ
#write.table(species_occ, file="./niche_modelling/occurrences/unclean/number_occurrences_species.csv" , sep=";",  row.names=FALSE, quote=FALSE)


