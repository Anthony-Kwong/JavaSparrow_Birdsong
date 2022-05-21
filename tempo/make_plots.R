#script to generate plots for tempo project

setwd("~/Documents/GitHub/JavaSparrow_Birdsong/tempo/figures/son_sf_plots/")

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
names = c("tempo","log gap length", "gap length standard deviation", "variance score",
          "log average song duration", "gap score")

#generate plots
for(i in 1:length(data)){
  #name of statistic
  statistic = names[i]
  
  #make plot
  plot <- ggplot(data[[i]], aes(x = sf_phenotype, y = son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    ggtitle("With Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("son's",statistic)) +
    xlab(paste("social fathers's",statistic)) +
    theme_bw()
  
  #title plot using the statistics
  plot_name = paste(statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name)
}

#plots for data without intros----

nointro_data <- list(nointro_tempo_df,
                     nointro_gap_df,
                     nointro_std_gap_df,
                     nointro_varscore_df,
                     nointro_durations_df,
                     nointro_gap_score_df)

#rename columns to son and sf phenotype 
nointro_data = lapply(nointro_data, function(df){
  temp <- df %>%
    dplyr::rename(son_phenotype = 2)
  res <- temp %>%
    dplyr::rename(sf_phenotype = 3)
  return(res)
})

#generate plots
for(i in 1:length(nointro_data)){
  #name of statistic
  statistic = names[i]
  
  #make plot
  plot <- ggplot(nointro_data[[i]], aes( x=sf_phenotype, y=son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    ggtitle("Without Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("son's",statistic)) +
    xlab(paste("social fathers's",statistic)) +
    theme_bw()
  
  #title plot using the statistics
  plot_name = paste("nointro ",statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name)
}

#plot for phenotype vs age ----

setwd("../son_age_plots/")
#full songs
for(i in 1:length(data)){
  #name of statistic
  statistic = names[i]
  
  #make plot
  plot <- ggplot(data[[i]], aes(x = Age_Rec, y = son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    ggtitle("With Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("son's",statistic)) +
    xlab(paste("son's age at recording")) +
    theme_bw()
  
  #title plot using the statistics
  plot_name = paste(statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name)
}

#nointro_songs
for(i in 1:length(nointro_data)){
  #name of statistic
  statistic = names[i]
  
  #make plot
  plot <- ggplot(nointro_data[[i]], aes( x=Age_Rec, y=son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    ggtitle("Without Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("son's",statistic)) +
    xlab(paste("son's age at recording")) +
    theme_bw()
  
  #title plot using the statistics
  plot_name = paste("nointro ",statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name)
}

#boxplots to show gap length for each transition type ----
head(gap_data)

library(ggplot2)
#plot with all transitions types, all gaps
ggplot(gap_data, aes(col = transitions, gap_dur)) +
  geom_boxplot()

#This large plot has too many classes, it's hard to read

#Lets group the transitions by the first note class for more clarity

library(stringr)

#Obtain starting note classes for each gap

gap_box_df = gap_data %>%
  dplyr::mutate(start_class = str_extract(transitions, "[^-]+" ))

#replace / in classes with * for ease of saving when it comes to file names
gap_box_df$start_class <- str_replace(gap_box_df$start_class, "[/]", "*")


#make a plot for each start class
start_classes = unique(gap_box_df$start_class) 

setwd("~/Documents/GitHub/JavaSparrow_Birdsong/tempo/figures/transition_plots/")

for(i in 1:length(start_classes)){
  class = start_classes[i]
  
  data = gap_box_df %>%
    dplyr::filter(start_class == class)
  
  plot <- ggplot(data, aes(y = log(gap_dur), col = transitions)) +
    geom_boxplot() +
    ylab("log gap durations") + 
    xlab(" transitions")
  #using class here would break code because it has /
  plot_name = paste("transition_plot_", class,".png" , sep ="")
  ggsave(plot, file = plot_name)
}

big_boxplot <- ggplot(gap_box_df, aes(y = log(gap_dur), col = transitions)) +
  geom_boxplot() +
  ylab("log gap durations") + 
  xlab(" transitions")

ggsave(big_boxplot, file = "all_transition_plot.pdf", height = 20, width = 20)

ggsave(big_boxplot, file = "all_transition_plot.png", height = 20, width = 20)


### code to help write previous loop
class = "Slope"

data = gap_box_df %>%
  dplyr::filter(start_class == class)

ggplot(data, aes(y = gap_dur, col = transitions)) +
  geom_boxplot() +
  ylab("gap durations") + 
  xlab(" transitions")

x = "Slope-Curve"
str_extract(x, "[^-]+")

