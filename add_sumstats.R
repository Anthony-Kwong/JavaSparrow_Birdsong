#add summary statistics for computational unit table

library(warbleR)
library(tidyverse)
library(birdsong.tools)

#read in data
anthony = read_csv("~/work/PhD/Data/note_excision/datasets/optimised_table.csv")

#add in sel column
anthony = add_selections(anthony)

test = anthony %>%
  dplyr::rename(selec = sel)

#compute SS
SS = spectro_analysis(test, bp = c(0.4,999) , parallel = 2, 
                      path = "~/work/PhD/Data/note_excision/song_data/")

write_csv(SS, file = "~/work/PhD/Data/note_excision/datasets/comp_data.csv")

