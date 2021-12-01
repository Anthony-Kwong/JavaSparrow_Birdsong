#read in data
library(tidyverse)
library(birdsong.tools)

becky = read_csv("~/work/PhD/Data/Files for Anthony/Files for Anthony/JavaSparrow_UnitTable.csv")
#anthony = read_csv("~/work/PhD/Data/note_excision/datasets/unit_table.csv")
setwd("/Users/shingkwong/work/PhD/Data/note_excision/code")

#specify optimization weights
int_penalty = 0
mismatch_penalty = 1

#specify objective function to minimise
loss <- function(min_amp){
  print(min_amp)
  cmd = paste("/Applications/Praat.app/Contents/MacOS/Praat --run excision_only.praat",min_amp)
  system(cmd)
  notes = read_csv("~/work/PhD/Data/note_excision/datasets/unit_table.csv")
  
  #filter out all silences
  unit_table = notes %>%
    dplyr::filter(interval_label == "sounding") %>%
    mutate(note_label = "unclassified")
  colnames(unit_table) = c("sound.files", "interval_label","start","end","duration",
                           "min_amp","note_label")  
  
  #compare unit tables
  diff_table = compare_excision(unit_table, becky)
  
  #compute mean square difference
  sq_diff = sum(diff_table$diff)/ncol(diff_table)
  
  #compute proportion of mismatches
  mismatch = which(diff_table$matched == F) %>% length()
  
  #compute weighted difference
  diff = int_penalty*sq_diff + mismatch_penalty*mismatch
  
  print(diff)
  return(diff)
}

xmin = optimize(loss , interval = c(-35,-15), tol = 10)
#takes tol=10 takes 1.3 hours

#save final unit table
notes = read_csv("~/work/PhD/Data/note_excision/datasets/unit_table.csv")
#filter out all silences
unit_table = notes %>%
  dplyr::filter(interval_label == "sounding") %>%
  mutate(note_label = "unclassified")
colnames(unit_table) = c("sound.files", "interval_label","start","end","duration",
                         "min_amp","note_label") 
write_csv(unit_table, file = "~/work/PhD/Data/note_excision/datasets/optimised_table.csv" )