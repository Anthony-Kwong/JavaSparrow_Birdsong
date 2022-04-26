#Analyzing introductions

#load libraries
pacman::p_load(magrittr, coxme)

#this script will investigate the relationships between the introductions and age

#Length of intros (by note count) vs age

#read in full unit table
full_ut <- read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv")
#read in smaller unit table with intros removed
nointro_ut <- read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/UnitTable_nointro.csv")
nointro_ut$X <- NULL

#get IDs of all birds
birds <- unique(nointro_ut$song_individual)

#compute the number of notes removed for each recording
sound_files = unique(nointro_ut$sound.files)

intro <- lapply(sound_files, function(file){
  #get bird ID corresponding to the recording
  data = full_ut %>%
    dplyr::filter(sound.files==file)
  ID = data$song_individual[1]
  
  #count notes in full unit table
  full_len <- full_ut %>%
    dplyr::filter(sound.files==file) %>%
    nrow()
  
  #count notes in the reduced table with no intros
  nointro_len <- nointro_ut %>%
    dplyr::filter(sound.files==file) %>%
    nrow()
  
  intro_len = full_len - nointro_len
  tibble::tibble(Bird.ID = ID, intro_len = intro_len, sound.file = file)
})

intro_data = do.call(rbind, intro)
#check no negative lengths
any(intro_data$intro_len < 0)

#add metadata based on Bird.ID
library(birdsong.tools)

meta.data = readr::read_csv("~/Documents/GitHub/JavaSparrow_Birdsong/data/meta_data.csv")
intro_df <- add_metadata(intro_data, meta.data, cols = c(5,6))

#regress intro_length against age, without accounting for genetic/social relationships
model1 <- lme4::lmer(intro_len ~ log(Age_Rec) + (1|Bird.ID) + (1|Clutch), data = intro_df)

#add social father values

#load in kinship matrices for lmekin models
kin.trim <- readRDS("~/Documents/GitHub/JavaSparrow_Birdsong/data/kin.trim.rds")
null.kin <- readRDS("~/Documents/GitHub/JavaSparrow_Birdsong/data/null.kin.rds")


tempo.lme = lmekin(formula = intro_len ~ sf_tempo_avg + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*kin.trim, #*2 because kinship halves the correlation
                   data = tempo_df )
