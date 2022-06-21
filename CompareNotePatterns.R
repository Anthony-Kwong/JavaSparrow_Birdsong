###################################################################
#	Have a go with warbleR and annotating spectrograms
#
# mrm: Whaley Range, 6 July 2021
###################################################################

rm( list=ls() )	# Wipe the slate clean
source( "AnnotateSpectrogram.R" ) 

###################################################################
#	Read the data
###################################################################

# Common elements of paths
#wav.subdir <- "~/work/PhD/Data/note_excision/song_data/"
wav.subdir <- "~/Dropbox (The University of Manchester)/FINAL FILES/RenamedKoe/"

# Read a CSV file describing all notes as found by Anthony
# Also add a song ID
anthony.note.df <- read.csv( "./results/comp_unit_table.csv", header=TRUE )
anthony.note.df$song <- sapply( anthony.note.df$sound.files,
	function( filename ) { sub( ".wav$", "", filename ) } )
head( anthony.note.df )

# Read a CSV file describing notes as found by Becky
becky.note.df <- read.csv( "~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Cleaned Files/UnitTable_20201122_Cleaned2.csv"
                           , header=TRUE )
head( becky.note.df )

###################################################################
#	Draw a spectrogram with bars indicating note boundaries
###################################################################

# Draw an example
my.song <- "JS0171-20180511-001"
draw.spectrogram( my.song, 
	anthony.note.df, becky.note.df,
	paste( wav.subdir )
)

# Blurt it to a file, creating a Figures directory if necessary
fig.subdir <- "Figures"
if( !dir.exists( fig.subdir ) ) { dir.create( fig.subdir ) }

fig.path <- paste( fig.subdir, "/", my.song, ".pdf", sep="" ) ;
dev.copy2pdf( file=fig.path ) 

#Generate plots for a representative sample (contains different colors) ----
library(tidyverse)

anthony.note.df$song <- sapply( anthony.note.df$sound.files,
                                function( filename ) { sub( ".wav$", "", filename ) } )

wavs <- list.files(pattern = ".wav$", ignore.case = TRUE,
                   path = "~/work/PhD/Data/note_excision/song_data/")
birds = sub("\\-.*","",wavs) 
meta_data = read_csv("~/work/PhD/Data/Files for Anthony/Files for Anthony/JavaSparrow_Metadata.csv")
tutor_lines = lapply(birds, function(ID){ meta_data$Line[which(meta_data$`Bird ID` == ID)]} ) %>% unlist()
meta_table = tibble(wavs, tutor_lines = as.factor(tutor_lines) )

set.seed(1800)
table_sample = meta_table %>%
  dplyr::group_by(tutor_lines) %>%
  sample_n(2)

song_list = lapply( table_sample$wavs,
        function( filename ) { sub( ".wav$", "", filename ) } ) %>% unlist()

#loop through each song and plot
for(i in 1:length(song_list)){
  current_song = song_list[i]
  draw.spectrogram( current_song,
                    anthony.note.df, becky.note.df,
                    paste( wav.subdir )
  )
  
  fig.subdir <- "Figures"
  if( !dir.exists( fig.subdir ) ) { dir.create( fig.subdir ) }

  fig.path <- paste( fig.subdir, "/", current_song, ".pdf", sep="" ) ;
  dev.copy2pdf( file=fig.path )
}

