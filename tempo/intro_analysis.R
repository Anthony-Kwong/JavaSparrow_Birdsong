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
  
  #compute the mean intro length of the individual and return dataframe, log then average
  tibble::tibble(Bird.ID = bird, avg_intro_len = mean( log(intro_lengths+1) ))
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

#full model
intro.lme = lmekin(formula = avg_intro_len ~ sf_intro_len + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*kin.trim, #*2 because kinship halves the correlation
                   data = intro_df )

#run age adjustment as well just to check

#checking models, We fitted a standard lme because the current tools don't work with lmekin

# library(performance)
# intro.lme_simp = lme4::lmer(formula = avg_intro_len ~ sf_intro_len + log(Age_Rec) + (1|Clutch),
#                             data = intro_df)
# check_overdispersion(intro.lme_simp)
# check_zeroinflation(intro.lme_simp)
#not necessary because we aren't using count data in the models. they are averages

#no genetics model
intro.lme2 = lmekin(formula = avg_intro_len ~ sf_intro_len + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = null.kin,
                   data = intro_df )

#no clutch model
intro.lme3 = lmekin(formula = avg_intro_len ~ sf_intro_len + log(Age_Rec) + (1|Bird.ID),
                    varlist = 2*kin.trim,
                    data = intro_df )

#neither model
intro.lme4 = lmekin(formula = avg_intro_len ~ sf_intro_len + log(Age_Rec) + (1|Bird.ID),
                    varlist = null.kin,
                    data = intro_df)

#use function from tempo_stats script
intro_res = extract_coxme_table(intro.lme)

tibble::tibble(response = "intro_length(num)", beta_sf = intro_res$beta[2], p_sf = intro_res$p[2],
               beta_logAge = intro_res$beta[3], p_logAge = intro_res$p[3], 
               #report chi square tests
               pedigree = 1-pchisq(2*(intro.lme$loglik - intro.lme2$loglik),1), 
               clutch = 1-pchisq(2*(intro.lme$loglik - intro.lme3$loglik),1), 
               neither = 1-pchisq(2*(intro.lme$loglik - intro.lme4$loglik),2),
               full_loglik = intro.lme$loglik, no.ped_loglike = intro.lme2$loglik, no.clutch_loglike = intro.lme3$loglik, neither_loglike = intro.lme4$loglik)



#plot relationship

library(ggplot2)
library(ggthemes)

#base plot 
ggplot(intro_df, aes(avg_intro_len, log(Age_Rec) )) +
  geom_point() +
  theme_stata()

ggplot(intro_df, aes(avg_intro_len, log(Age_Rec) )) +
  geom_point() +
  theme_economist()

intro_vs_age <- ggplot(intro_df, aes( log(Age_Rec), avg_intro_len)) +
  geom_point() +
  theme(legend.position="none") +
  xlab("average log introduction length") + 
  ylab("log age at recording") +
  theme_bw(base_size = 14)

ggsave(intro_vs_age, filename = "~/Documents/GitHub/JavaSparrow_Birdsong/tempo/figures/intro_vs_age.png")

intro_plot <- ggplot(intro_df, aes(avg_intro_len, sf_intro_len)) +
  geom_point() + 
  #add regression line
  #geom_smooth(method = lm, color = "blue") +
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
model2 <- lmer(intro_len ~ log(Age_Rec)+ (1|Bird.ID), data = intro_df2)
summary(model2)
#significant age effect
