pacman::p_load(tidyverse)

#load in raw data sets ----
acoustic_params = read_csv("~/work/PhD/Data/20210303/20210303/Acoustic_parameter_outputs_t2_o50_bp400_notfast_20201122.csv")

manual_class = read_csv("~/work/PhD/Data/20210303/20210303/UnitTable_20201122_Cleaned2.csv")

#bind vector of responses (manually clustered notes)
labelled_data = cbind(acoustic_params, label_family = manual_class$label_family)

#check variables
str(acoustic_params)

#data cleaning ----
acoustic_params$X1 <- NULL
acoustic_params <- acoustic_params %>%
  mutate_if(is.character,factor)

manual_class$label_family = as.factor(manual_class$label_family)
#clustering ----

#take out the soundfiles and select columns
features = acoustic_params[,3:ncol(acoustic_params)]

set.seed(12)
clus_res = kmeans(features, centers = 16, iter.max = 100, algorithm="MacQueen")

factoextra::fviz_cluster(clus_res, geom = "point", data = features)
