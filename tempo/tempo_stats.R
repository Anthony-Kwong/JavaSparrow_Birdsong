#script for doing all the analysis for the tempo paper 

#load packages
library(birdsong.tools)
pacman::p_load(coxme, magrittr)

#load functions
source("~/Documents/GitHub/JavaSparrow_Birdsong/functions/age_adjust.R")

#load in all relevant datasets

#import manual note classifications
full_ut=read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv",stringsAsFactors=FALSE,fill=TRUE)
#unit table without intros
nointro_ut <- read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/UnitTable_nointro.csv")

#import metadata on individual birds
meta.data=read.csv("~/Documents/GitHub/JavaSparrow_Birdsong/data/meta_data.csv")
meta.data$X1 <- NULL
meta.data=meta.data[order(meta.data$Bird.ID),]

#load the kinship matrices 
kin.trim = readRDS(file = "~/Documents/GitHub/JavaSparrow_Birdsong/data/kin.trim.rds")
null.kin = readRDS(file = "~/Documents/GitHub/JavaSparrow_Birdsong/data/null.kin.rds")

#create function to compute all the statistics for each individual bird----

#input: unit table, meta data table
#output: dataframe with bird ID, son's statistics, sf statistics, son's BD, son's DOR, son's Age_Rec, Clutch, SF ID, Sf DOR, sf Age_rec
compute_res <- function(unit_table, meta.data){
  
  #compute tempo----
  tempo_data = get_tempo(unit_table) %>% 
    dplyr::mutate(ID = sub("-.*", "", sound.files))
  
  #compute durations----
  durations = list()
  recordings = unique(unit_table$sound.files)
  for(i in 1:length(recordings)){
    tar = recordings[i]
    rec_data = unit_table %>%
      dplyr::filter(sound.files == tar)
    durations[[i]] <- tibble::tibble(Bird.ID = rec_data$song_individual[1], 
                                     dur = get_duration(rec_data))
  }
  durations = do.call(rbind, durations)
  durations = durations %>%
    dplyr::mutate(ID = Bird.ID)
  
  #compute gap table----
  gap_data = get_gaps(unit_table) %>% 
    dplyr::mutate(ID = sub("-.*", "", sound.files))
  
  #average tempo ,gaps, duration for each bird----
  stats_bundle = list(tempo_data, gap_data, durations)
  
  #get ID of each bird
  bird_ID = unique(tempo_data$ID)
  
  avg_stats = lapply(stats_bundle, function(tab){
    avg_stat = list()
    for(i in 1:length(bird_ID)){
      tar = bird_ID[i]
      bird_data = tab %>%
        dplyr::filter(ID == tar)
      
      #compute the mean stat for bird, 2nd column of each dataframe
      stat = as.numeric(unlist(bird_data[,2]))
      avg_stat[[i]] = tibble::tibble(Bird.ID = tar, mean(stat))
    }
    #bind measurements for all the birds together
    data = do.call(rbind, avg_stat)
    return(data)
  }) #end of lapply
  
  #rename the columns by statistic name, mutate doesn't like string arguments
  stat_names = sapply(stats_bundle, function(tab){
    colnames(tab)[2]
  })
  avg_stats[[1]] = avg_stats[[1]] %>%
    dplyr::mutate(avg_tempo = `mean(stat)`)
  avg_stats[[2]] = avg_stats[[2]] %>%
    dplyr::mutate(avg_gap_dur = `mean(stat)`)
  avg_stats[[3]] = avg_stats[[3]] %>%
    dplyr::mutate(avg_song_dur = `mean(stat)`)
  
  #drop redundant second col
  avg_stats = lapply(avg_stats, function(tab){
    subset(tab, select = -c(2))
  })
  
  #gap standard deviation
  std_gap = list()
  for(i in 1:length(bird_ID)){
    tar = bird_ID[i]
    ind_data = gap_data %>%
      dplyr::filter(ID == tar)
    std_gap[[i]] <- tibble::tibble(Bird.ID = tar, std_gap = var(ind_data$gap_dur)^0.5)
  }
  std_gap= do.call(rbind,std_gap)
  
  #gap score ----
  #compute vector of mean gaps for each transitition for the whole population
  trans = unique(gap_data$transitions)
  
  pop_means = list()
  for(i in 1:length(trans)){
    tar = trans[i]
    fil_data = gap_data %>%
      dplyr::filter(transitions == tar)
    pop_means[[i]] = tibble::tibble(
      mean_gap = mean(fil_data$gap_dur), 
      transitions = tar
    )
  }
  
  pop_means = do.call(rbind, pop_means)
  
  gap_scores = list()
  for(i in 1:length(bird_ID)){
    tar = bird_ID[i]
    ind_data = gap_data %>%
      dplyr::filter(ID == tar)
    gap_scores[[i]] = tibble::tibble(Bird.ID = tar,
                                     gap_score = gap_score(ind_data, pop_means))
  }
  
  gap_score_df = do.call(rbind, gap_scores)
  
  #gap variability score
  var_scores = list()
  for(i in 1:length(bird_ID)){
    tar = bird_ID[i]
    ind_data = gap_data %>%
      dplyr::filter(ID == tar)
    var_scores[[i]] <- tibble::tibble(Bird.ID = tar, var_score = gap_var_score(ind_data))
  }
  
  var_scores = do.call(rbind, var_scores)
  
  #total datasets (remove after finishing loop)
  # avg_stats
  # std_gap
  # gap_score_df
  # var_scores
  # 
  
  #add metadata (age, sf ID, sf stat, clutch)
  stats_data = c(list(gap_score_df, std_gap, var_scores),avg_stats) 
  
  fin_stat = lapply(stats_data, function(tab){
    #need 'the original birdsong data set to fit the model with Bird.ID, phenotype, social father's phenotype, Birth.Date, sf_Age_Rec, Age_Rec and Clutch(in that order). 
    #for building
    # tab = stats_data[[1]]
    
    #add son's Birth.Date, DOR, Age,Clutch, Social Father ID
    tab = add_metadata(tab, meta.data, cols = c(3:7))
    
    #add sf phenotype
    tab <- get_sf_phenotype(phenotype_table = tab, metadata = meta.data, phenotype_index =  2)
    tab <- dplyr::rename(tab, sf_y = sf_phenotype)
    
    #add sf DOR
    tab <- get_sf_phenotype(phenotype_table = tab, metadata = meta.data, phenotype_index =  4)
    tab <- dplyr::rename(tab, sf_DOR = sf_phenotype)
    
    #add sf Age_Rec
    tab <- get_sf_phenotype(phenotype_table = tab, metadata = meta.data, phenotype_index =  5)
    tab <- dplyr::rename(tab, sf_Age_Rec = sf_phenotype)
    
    tibble::tibble(tab[,1:2], tab[,8], tab[,3:7], tab[,9:10])
    
  })
  
  return(fin_stat)
}

#https://stackoverflow.com/questions/43720260/how-to-extract-p-values-from-lmekin-objects-in-coxme-package
extract_coxme_table <- function (mod){
  beta <- mod$coefficients$fixed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

#fit 6 models for each statistic 
#full model, no clutch model, no genetics model, with and without age adjust

#input: table with Bird.ID,y, sf_y , Birth.Date , DOR , Age_Rec, Clutch, kinship matrix to represent pedigree, null kinship matrix     
#Social.Father ,sf_DOR ,sf_Age_Rec  

#if model output is false, output statistics table instead

source("./functions/age_adjust.R")
fit_birdmodel <- function(stat_table, kin.trim , null.kin , model_output = T){
  
  #write formulas
  y = colnames(stat_table[2])
  sf_y = colnames(stat_table[3])
  
  base_M = paste(y, "~", sf_y)
  ff = as.formula(paste(base_M, "+", "log(Age_Rec) + (1|Bird.ID) + (1|Clutch)"))
  nc = as.formula(paste(base_M, "+", "log(Age_Rec) + (1|Bird.ID)"))
    
  #lm(paste0(x, '~', 'mpg'), data = mtcars)
  
  #standard models with no age adjustment----
  
  #full model
  M1 = lmekin(ff, varlist = 2*kin.trim, data = stat_table)
  #no genetics model
  M2 = lmekin(ff, varlist = 2*null.kin, data = stat_table)
  #no clutch model
  M3 = lmekin(nc, varlist = 2*kin.trim, data = stat_table)

  #age adjustment models
  AM1 <- age_adjust(model = M1, tol = 10^-6, data = stat_table, kin = kin.trim)
  AM2 <- age_adjust(model = M2, tol = 10^-6, data = stat_table, kin = null.kin)
  AM3 <- age_adjust(model = M3, tol = 10^-6, data = stat_table, kin = kin.trim)
  
  if(model_output == TRUE){
    models = list(M1, M2, M3, AM1, AM2, AM3)
    return(models)
  } else {
    #report age adjusted model
    mod_tab = extract_coxme_table(AM1)
    beta_sfy = mod_tab$beta[2]
    p_sfy = mod_tab$p[2]
    beta_age = mod_tab$beta[3]
    p_age = mod_tab$p[3]
    
    #test sig. of genetics using chi square log likelihood ratio test. The test statistic is 2*difference in the loglikelihood. Distribution is approximately chi square with df being df(full model)-df(reduced model) 
    chi_g = 1-pchisq(2*(AM1$loglik - AM2$loglik),1)
    #test sif. of clutch
    chi_c = 1-pchisq(2*(AM1$loglik - AM3$loglik),1)
    
    response = colnames(stat_table)[2]
    
    res_table = tibble::tibble(response = response, beta_sf = beta_sfy, p_sf = p_sfy, 
             beta_logAge = beta_age, p_age = p_age, pedigree = chi_g, clutch = chi_c)
    return(res_table)
  }
}


#compute all temporal stats----
bird_stats = compute_res(unit_table = full_ut, meta.data = meta.data)
#compute models
fullsong_models = lapply(bird_stats, FUN =  fit_birdmodel, kin.trim = kin.trim , null.kin = null.kin)
#compute result table
fullsong_tab = lapply(bird_stats, FUN =  fit_birdmodel, kin.trim = kin.trim , null.kin = null.kin, model_output = FALSE)
fullsong_tab = do.call(rbind, fullsong_tab)

#do the same for the nointro table----
bird_stats2 = compute_res(unit_table = nointro_ut, meta.data = meta.data)
#compute models
nointros_models = lapply(bird_stats2, FUN =  fit_birdmodel, kin.trim = kin.trim , null.kin = null.kin)
#compute result table
nointrosong_tab = lapply(bird_stats2, FUN =  fit_birdmodel, kin.trim = kin.trim , null.kin = null.kin, model_output = FALSE)
nointrosong_tab = do.call(rbind, nointrosong_tab)

