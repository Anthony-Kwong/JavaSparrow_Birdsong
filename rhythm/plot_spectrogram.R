#This script plots spectrograms using existing functions in this project

#first written to draw spectrograms for the rhythm project

# Common elements of paths
wav.subdir <- "~/Dropbox (The University of Manchester)/FINAL FILES/RenamedKoe/"

# Read a CSV file describing all notes as found by Anthony
# Also add a song ID, old code, this is for the computational note excision
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

#plot all songs from 171 and 170 for investigation
library(magrittr)

birds = c("JS0170", "JS0171","JS0037")

lapply(birds, function(b){
  #filter for the individual
  data = becky.note.df %>%
    dplyr::filter(song_individual == b)
  
  songs = unique(data$song)
  
  #loop through all songs and plot
  lapply(songs, function(s){
    draw.spectrogram( s,
                      #don't worry about different note excisions for now
                      becky.note.df, becky.note.df,
                      paste( wav.subdir ))

    fig.path = paste("./rhythm/figures/", s, ".pdf", sep="")
    dev.copy2pdf(file=fig.path)
  })
})


#170,171 have a series of super short followed by super long notes

#37 is pretty evenly spaced out as suggested by the plot
