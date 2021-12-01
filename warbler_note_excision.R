#note excision using warbleR

pacman::p_load(warbleR, optimx, tidyverse)
library(birdsong.tools)

#setwd("/Users/shingkwong/work/PhD/Data/note_excision/code")

#set path for sound files
warbleR_options(wav.path="~/work/PhD/Data/note_excision/song_data/")

#example starting point

#took max and min durations from becky's notes. 
note_excision = auto_detec(bp = c(0.4,22), mindur = 0.001, maxdur = 0.6,threshold = 2,
                           ssmooth = 0, hold.time = 0,
                           parallel = 2)
#threshold is an important parameter. It stops detection if the signal reaches X% of the peak spectrum

#control for bird family

#read in becky's notes
becky = read_csv("~/work/PhD/Data/Files for Anthony/Files for Anthony/JavaSparrow_UnitTable.csv")

#loss function ----

#params are taken in the log scale

note_loss <- function(params, song_list){
  
  print(params)
    
  unit_table = becky %>%
    dplyr::filter(sound.files %in% song_list)
  
  exp_params = exp(params)
    
  #this lines takes ~8mins
  res = optimize_auto_detec(X = unit_table, bp = c(0.4,22), mindur = exp_params["mindur"], maxdur = exp_params["maxdur"],
                            threshold = exp_params["threshold"], ssmooth = exp_params["ssmooth"], hold.time = exp_params["hold.time"],
                            parallel = 2, 
                            path = "~/work/PhD/Data/note_excision/song_data/")
  
  #generate computation unit table
  sampled_recordings = unit_table$sound.files %>% unique()
  
  comp_table = auto_detec(bp = c(0.4,22), mindur = exp_params["mindur"], maxdur = exp_params["maxdur"],
                          threshold = exp_params["threshold"], ssmooth = exp_params["ssmooth"], 
                          hold.time = exp_params["hold.time"], parallel = 2,
                          #ensure auto_detec works on sampled recordings only
                          flist = sampled_recordings)
  
  #overlap scores, range (0,1), 1 is best
  oscores = compute_oscores(manual_table = unit_table, comp_table = comp_table)
  
  #minimise this
  SS = -res$true.positives + res$false.positives + res$split.positives - (sum(oscores)) 
  #since we want to maximise the TP and oscores, minimise the FP and SP
  
  msg = paste("TP:", res$true.positives, "FP:", res$false.positives,
              "SP:", res$split.positives,
              "x/y:", oscores[1],
              "x/a:", oscores[2],
              "x2/yz:", oscores[3])
  print(msg)
  print(SS)
  return(SS)
}

start_values = c(mindur = 0.001, maxdur = 0.6, threshold = 2,
ssmooth = 1, hold.time = 0.05) %>% log()

#optimize on a subset of songs to get a better starting point ----

wavs <- list.files(pattern = ".wav$", ignore.case = TRUE,
                   path = "~/work/PhD/Data/note_excision/song_data/")

#account for tutor lines to ensure sample song batches represent the different lines proportionally

birds = sub("\\-.*","",wavs) 
meta_data = read_csv("~/work/PhD/Data/Files for Anthony/Files for Anthony/JavaSparrow_Metadata.csv")
tutor_lines = lapply(birds, function(ID){ meta_data$Line[which(meta_data$`Bird ID` == ID)]} ) %>% unlist()
meta_table = tibble(wavs, tutor_lines = as.factor(tutor_lines) )

candidate_points = list()
set.seed(1066)
iterations = 5

for(i in 1:iterations){
  msg = paste("Optimization Round: ",i)
  print(msg)

  sample = meta_table %>%
    dplyr::group_by(tutor_lines) %>%
    sample_n(2)
  
  recording_sample = sample$wavs
  
  candidate_points[[i]] = optimx(par = start_values, fn = note_loss, song_list = recording_sample , method = "Nelder-Mead", itnmax = 5)
}

candidate_points = do.call(rbind, candidate_points)
write_csv(candidate_points, file = "~/Documents/GitHub/JavaSparrow_Birdsong/results/candidate_excision_param.csv")

# Set a seed so we all have the same results
# set.seed(1)
# sub <- wavs[sample(1:length(wavs), 10)]
# 
# initial_res = optimx(par = start_values, fn = note_loss, song_list = sub , method = "CG", itnmax = 5)

#optimize excision on whole data set ----

point = candidate_points[1,]
new_start_point = c(mindur = point$mindur, maxdur = point$maxdur, threshold = point$threshold,
                    ssmooth = point$ssmooth, hold.time = point$hold.time)
  
#minimize by default
#best_param = optimx(par = start_values, fn = note_loss , method = "CG", itnmax = 5)
#14 hours...

best_param = optimx(par = new_start_point, fn = note_loss , song_list = wavs, method = "Nelder-Mead", itnmax = 5)

write_csv(best_param, file = "~/Documents/GitHub/JavaSparrow_Birdsong/results/best_excision_param.csv" )

opt_res = optimize_auto_detec(X = becky, bp = c(0.4,22), mindur = exp(best_param$mindur), maxdur = exp(best_param$maxdur),
                              threshold = exp(best_param$threshold), ssmooth = exp(best_param$ssmooth), hold.time = exp(best_param$hold.time),
                              parallel = 2, 
                              path = "~/work/PhD/Data/note_excision/song_data/")

#Finally generate computational unit table

#excise notes computationally, 50% overlap, no power tuning, 
note_excision = auto_detec(bp = c(0.4,22), mindur = exp(best_param$mindur), maxdur = exp(best_param$maxdur),
                           threshold = exp(best_param$threshold), ssmooth = exp(best_param$ssmooth),
                           hold.time = exp(best_param$hold.time),
                           parallel = 2)

write_csv(note_excision, file = "~/Documents/GitHub/JavaSparrow_Birdsong/results/comp_unit_table.csv" )
