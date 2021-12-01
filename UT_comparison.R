#This script compares our computational unit table with becky's manual unit table

library(tidyverse)
comp_ut = read_csv("./results/comp_unit_table.csv")
becky = read_csv("~/work/PhD/Files for Anthony/Files for Anthony/JavaSparrow_UnitTable.csv")

library(birdsong.tools)

#Find the notes in becky that correspond to the computational notes
match_vec = note_match(comp_table = comp_ut, manual_table = becky)

lab_comp_ut = cbind(comp_ut, note_label = as.factor(match_vec))

#Find the notes in computational table that match to becky...we need a new function

