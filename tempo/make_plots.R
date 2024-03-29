#script to generate plots for tempo project
library(ggplot2)
library(magrittr)
setwd("~/Documents/GitHub/JavaSparrow_Birdsong/")

#load dataset generated by tempo_stats.R
full_data = readRDS("./tempo/data/full_bird_stats.rds")
nointro_data = readRDS("./tempo/data/nointro_bird_stats.rds")

#get name of each statistics, ensure this matches with the dataframes collected in list
names = sapply(full_data, function(df){
  colnames(df)[2]
})
#replace underscores with spaces for neater plots
names = sub("_"," ",names)

#just redo the names for neatness
names[2] = "gap length standard deviation"
names[3] = "variance score"
names[4] = "tempo"
names[5] = "gap duration"
names[6] = 'song duration'

#rename columns to son and sf phenotype to work with ggplot
full_data = lapply(full_data, function(df){
  temp <- df %>%
    dplyr::rename(son_phenotype = 2)
  res <- temp %>%
    dplyr::rename(sf_phenotype = 3)
  return(res)
})

setwd("~/Documents/GitHub/JavaSparrow_Birdsong/tempo/figures/son_sf_plots/")

#generate plots
for(i in 1:length(full_data)){
  #name of statistic
  statistic = names[i]
  
  #compute limits of each plot
  sel_data = full_data[[i]]
  max_limit = max(c(sel_data$son_phenotype,na.omit(sel_data$sf_phenotype)))
  min_limit = min(c(sel_data$son_phenotype,na.omit(sel_data$sf_phenotype)))
  
  #make plot
  plot <- ggplot(sel_data, aes(x = sf_phenotype, y = son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    #ggtitle("With Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("son's",statistic,"")) +
    xlab(paste("social father's ",statistic,sep="")) +
    theme_bw(base_size = 14) +
    xlim(c(min_limit,max_limit)) +
    ylim(c(min_limit,max_limit))
    
  #title plot using the statistics
  plot_name = paste(statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name, width = 5, height = 5)
}

#plots for data without intros----

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
  
  #compute limits of each plot
  sel_data = nointro_data[[i]]
  max_limit = max(c(sel_data$son_phenotype,na.omit(sel_data$sf_phenotype)))
  min_limit = min(c(sel_data$son_phenotype,na.omit(sel_data$sf_phenotype)))
  
  #make plot
  plot <- ggplot(nointro_data[[i]], aes( x=sf_phenotype, y=son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    #ggtitle("Without Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("son's",statistic)) +
    xlab(paste("social father's",statistic)) +
    theme_bw(base_size = 14) +
    xlim(c(min_limit,max_limit)) +
    ylim(c(min_limit,max_limit))
  
  #title plot using the statistics
  plot_name = paste("nointro ",statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name, width = 5, height = 5)
}

#plot for phenotype vs age ----

setwd("../son_age_plots/")
#full songs
for(i in 1:length(full_data)){
  #name of statistic
  statistic = names[i]
  
  #make plot
  plot <- ggplot(full_data[[i]], aes(x = log(Age_Rec), y = son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    #ggtitle("With Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("bird's",statistic)) +
    xlab(paste("bird's log age at recording")) +
    theme_bw(base_size = 14)
  
  #title plot using the statistics
  plot_name = paste(statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name,width = 5, height = 5)
}

#nointro_songs
for(i in 1:length(nointro_data)){
  #name of statistic
  statistic = names[i]
  
  #make plot
  plot <- ggplot(nointro_data[[i]], aes( x=log(Age_Rec), y=son_phenotype)) +
    geom_point() + 
    #add regression line
    #geom_smooth(method = "lm", color = "blue")+
    #ggtitle("Without Introductions") +
    #remove redundant legend because we have no classes
    theme(legend.position="none") +
    ylab(paste("bird's",statistic)) +
    xlab(paste("bird's log age at recording")) +
    theme_bw(base_size = 14)
  
  #title plot using the statistics
  plot_name = paste("nointro ",statistic,"_plot.png",sep = "")
  #save plot in tempo directory
  ggsave(plot, file = plot_name,width = 5, height = 5)
}

#boxplots to show gap length for each transition type ----
#import manual note classifications
full_ut=read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv",stringsAsFactors=FALSE,fill=TRUE)
#unit table without intros
nointro_ut <- read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/UnitTable_nointro.csv")

#filter out NLP/Other
full_ut = full_ut %>%
  dplyr::filter(note_label != "NLP_Other")

gap_data = get_gaps(full_ut) %>% 
  dplyr::mutate(ID = sub("-.*", "", sound.files))

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
    xlab(" transitions") +
    theme_bw()
  #using class here would break code because it has /
  plot_name = paste("transition_plot_", class,".png" , sep ="")
  ggsave(plot, file = plot_name)
}

big_boxplot <- ggplot(gap_box_df, aes(y = log(gap_dur), col = transitions)) +
  geom_boxplot() +
  ylab("log gap durations") + 
  xlab(" transitions") +
  theme_bw()

ggsave(big_boxplot, file = "all_transition_plot.pdf", height = 20, width = 20)

ggsave(big_boxplot, file = "all_transition_plot.png", height = 20, width = 20)

#summary statistics for gaps of each transition type

gap_data
trans = unique(gap_data$transitions)
trans_data = lapply(trans, function(k){
  d = gap_data %>% 
    dplyr::filter(transitions == k)
  tibble::tibble(mean_duration = mean(d$gap_dur), std_duration = sd(d$gap_dur), transition = k, number = nrow(d))
})
trans_data = do.call(rbind, trans_data)
write.csv(trans_data, file = "../../data/transition_data.csv")


### code to help write previous loop----
class = "Slope"

data = gap_box_df %>%
  dplyr::filter(start_class == class)

ggplot(data, aes(y = gap_dur, col = transitions)) +
  geom_boxplot() +
  ylab("gap durations") + 
  xlab(" transitions")

x = "Slope-Curve"
str_extract(x, "[^-]+")

