#Analyzing introductions

#load libraries
pacman::p_load(magrittr, coxme)

#this script will investigate the relationships between the father's intro length and the son's intro length

#Length of intros (by note count) vs age

#read in full unit table
full_ut <- read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv")
#read in smaller unit table with intros removed
nointro_ut <- read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/UnitTable_nointro.csv")
nointro_ut$X <- NULL

#get IDs of all birds
birds <- unique(nointro_ut$song_individual)

intro <- lapply(birds, function(bird){
  #get all recordings from the bird
  bird_data = full_ut %>%
    dplyr::filter(song_individual == bird)
  
  #get truncated unit table for the bird
  nointro_bird_data = nointro_ut %>%
    dplyr::filter(song_individual == bird)
  
  #get all the songs from the bird
  songs = unique(bird_data$sound.files)
  
  #get the intro length of each song
  intro_lengths = sapply(songs, function(song){
    full_len <- bird_data %>%
      dplyr::filter(sound.files == song) %>%
      nrow()
    
    nointro_len <- nointro_bird_data %>%
      dplyr::filter(sound.files == song) %>%
      nrow()
    
    introduction_len = full_len - nointro_len
  })
  
  #compute the mean intro length of the individual and return dataframe
  tibble::tibble(Bird.ID = bird, avg_intro_len = mean(intro_lengths))
})

#bind everything together into one dataframe
intro_data = do.call(rbind, intro)

#check no negative lengths
any(intro_data$avg_intro_len < 0)

#add metadata based on Bird.ID
library(birdsong.tools)

#add metadata
meta.data = readr::read_csv("~/Documents/GitHub/JavaSparrow_Birdsong/data/meta_data.csv")
intro_df <- add_metadata(intro_data, meta.data, cols = c(5,6))

#add social father phenotype
intro_df <- get_sf_phenotype(phenotype_table = intro_df, metadata = meta.data, phenotype_index = 2)
intro_df = dplyr::rename(intro_df, sf_intro_len= sf_phenotype)


#load in kinship matrices for lmekin models
kin.trim <- readRDS("~/Documents/GitHub/JavaSparrow_Birdsong/data/kin.trim.rds")
null.kin <- readRDS("~/Documents/GitHub/JavaSparrow_Birdsong/data/null.kin.rds")


intro.lme = lmekin(formula = log(avg_intro_len+1) ~ log(sf_intro_len+1) + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*kin.trim, #*2 because kinship halves the correlation
                   data = intro_df )

#no sig result for sf, sig result for age

#no genetics model
intro.lme2 = lmekin(formula = log(avg_intro_len+1) ~ log(sf_intro_len+1) + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = null.kin,
                   data = intro_df )

#no clutch model
intro.lme3 = lmekin(formula = log(avg_intro_len+1) ~ log(sf_intro_len+1) + log(Age_Rec) + (1|Bird.ID),
                    varlist = null.kin,
                    data = intro_df )

#full model vs no genetics model
1-pchisq(2*(intro.lme$loglik - intro.lme2$loglik),1)
#full model vs no clutch model
1-pchisq(2*(intro.lme$loglik - intro.lme3$loglik),1)

#no sig effects for genetics and clutch, clutch is close with p-value=0.07

#plot relationship

library(ggplot2)
library(ggthemes)

#base plot 
ggplot(intro_df, aes(log(avg_intro_len+1), log(Age_Rec) )) +
  geom_point() +
  theme_stata()

ggplot(intro_df, aes(log(avg_intro_len+1), log(Age_Rec) )) +
  geom_point() +
  theme_economist()

intro_vs_age <- ggplot(intro_df, aes( log(avg_intro_len+1), log(Age_Rec), color = "orangered" )) +
  geom_point() +
  geom_smooth(method = lm, color = "blue") +
  theme(legend.position="none") +
  xlab("log average introduction length") + 
  ylab("log(Age) at recording")

ggsave(intro_vs_age, filename = "~/Documents/GitHub/JavaSparrow_Birdsong/tempo/figures/intro_vs_age.png")

intro_plot <- ggplot(intro_df, aes(log(avg_intro_len+1), log(sf_intro_len+1), color = "orangered")) +
  geom_point() + 
  #add regression line
  geom_smooth(method = lm, color = "blue") +
  #remove legend as it's not needed here
  theme(legend.position="none") + 
  xlab("log average introduction length (sons)") +
  ylab("log average introduction length (social father)")

ggsave(intro_plot, filename = "~/Documents/GitHub/JavaSparrow_Birdsong/tempo/figures/son_sf_plots/intro_length_plot.png")

#We know investigate the relationship between age and intro length

songs = unique(full_ut$sound.files)

intro_df2 <- lapply(songs, function(recording){
  bird_data = full_ut %>%
    #get the unit table for one recording
    dplyr::filter(sound.files == recording)
  
  ID = bird_data$song_individual[1]
  
  #note count for full runit table
  full_len = full_ut %>%
    #get the unit table for one recording
    dplyr::filter(sound.files == recording) %>%
    nrow()
  
  #note count for smaller unit table without the intro
  nointro_len = nointro_ut %>%
    dplyr::filter(sound.files == recording) %>%
    nrow()
  
  intro_length = full_len - nointro_len
  
  tibble::tibble(Bird.ID = ID, intro_len = intro_length, sound.files = recording)
})

intro_df2 <- do.call(rbind, intro_df2)

#add in age
intro_df2 <- add_metadata(intro_df2, metadata = meta.data, cols = 5)
intro_df2$Bird.ID <-as.factor(intro_df2$Bird.ID)

#regress intro_len with age
library(lme4)
model2 <- lmer(log(intro_len+1) ~ log(Age_Rec)+ (1|Bird.ID), data = intro_df2)
summary(model2)
#significant age effect
