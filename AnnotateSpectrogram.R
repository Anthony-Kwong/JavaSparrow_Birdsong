###################################################################
#	Define a circle of functions to plot a spectrogram and add 
# coloured bars indicating note boundaries
#
# mrm: Whaley Range, 6 July 2021
###################################################################

# Load packages we'll need
library( warbleR )
library( seewave )
library( RColorBrewer ) 

# Define a function to draw a single coloured rectangle
draw.rectangle <- function( t.sec, y, my.color ) {
	# unpack the args
	x.left <- t.sec[1] ;
	x.right <- t.sec[2] ;
	y.bottom <- y[1] ;
	y.top <- y[2] ;
	
	rect( x.left, y.bottom, x.right, y.top, col=my.color, border=NA )
}

# Define a function to draw a bar to indicate intervals of silence and singing
# Here note.df should have just two columns, start and end, in that order.
draw.note.bar <- function( note.df, duration.secs, bar.midline.kHz, bar.width.kHz, silence.col, note.col ) 
{
	# Set the vertical position of the bar 
	bar.half.width <- bar.width.kHz / 2
	y.top <- bar.midline.kHz + bar.half.width
	y.bottom <- bar.midline.kHz - bar.half.width
	y.lim <- c( y.bottom, y.top )
	
	# Draw a long strip in the color for silences
	draw.rectangle( c(0, duration.secs), y.lim, silence.col ) ;

	# Add smaller rectangles representing notes
	apply( note.df, 1, 
		function( my.row ) { draw.rectangle( my.row[1:2], y.lim, note.col ) }
	) ;
}

# Define a function to draw a whole spectrogram and add the coloured bars
draw.spectrogram <- function( 
	my.song, anthony.df, becky.df, wav.dir, 
	min.freq.kHz = 0, max.freq.kHz = 6.5,
	bar.width.frac = 0.04, 
	bar.midline.kHz = c( 0.125, 0.3 )
) {
	# Construct a path to a WAV file and read that puppy
	wav.file <- paste( my.song, ".wav", sep="" ) 
	wav.path <- paste( wav.dir, wav.file, sep="/" ) 
	wav.obj <- read_sound_file( wav.path ) 
	
	# Draw the spectrogram
	spectro( wav.obj, 
		fftw = TRUE,	# Use the mighty FFTW 
		flim =c(min.freq.kHz, max.freq.kHz), 	# Limit frequency range to where there's action
		osc = TRUE,		# Plot sound pressure
		palette = reverse.topo.colors,
		main = my.song	# Use the song's name as a title
	)

	# Add coloured bars to indicate intervals of silence and singing
	pair.pal <- brewer.pal( 12, "Paired" ) 
	
	duration.secs <- length( wav.obj ) / wav.obj@samp.rate
	bar.width.kHz <- (max.freq.kHz - min.freq.kHz) * bar.width.frac/2
	draw.note.bar( 
		subset( becky.df, song==my.song, select=c(start, end) ),
		duration.secs,
		bar.midline.kHz[1], bar.width.kHz,
		pair.pal[1], pair.pal[2]
	) ;
		
	draw.note.bar( 
		subset( anthony.df, song==my.song, select=c(start, end) ),
		duration.secs,
		bar.midline.kHz[2], bar.width.kHz,
		pair.pal[7], pair.pal[8]
	) ;
}
