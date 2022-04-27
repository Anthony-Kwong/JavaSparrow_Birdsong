#script to generate plots for tempo project

setwd("./son_sf_plots/")

#run all the code chunks from tempo_project first

#datasets generated from the code chunks are..
tempo_df
gap_df
std_gap_df
varscore_df 
durations_df
gap_score_df

data = list(tempo_df,
            gap_df,
            std_gap_df,
            varscore_df,
            durations_df,
            gap_score_df)

#rename columns to son and sf phenotype 
data = lapply(data, function(df){
  temp <- df %>%
    dplyr::rename(son_phenotype = 2)
  res <- temp %>%
    dplyr::rename(sf_phenotype = 3)
  return(res)
})

#get name of each statistics, ensure this matches with the dataframes collected in list
names = c("tempo","gap_length", "gap_length_std", "varscore", "avg_song_duration", "gap_score")

#generate plots
for(i in 1:length(data)){
  #name of statistic
  statistic = names[i]
  #make plot
  plot <- ggplot(data[[i]], aes(son_phenotype, sf_phenotype)) +
    geom_point() + 
    ggtitle(statistic)
  #title plot using the statistics
  plot_name = paste(statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name)
}

