#supervised approach for clustering notes

pacman::p_load(tidyverse, Rfast)
library(tidymodels)
library(birdsong.tools)

#load in raw data sets ----
acoustic_params = read_csv("~/work/PhD/Data/20210303/20210303/Acoustic_parameter_outputs_t2_o50_bp400_notfast_20201122.csv")

manual_class = read_csv("~/work/PhD/Data/20210303/20210303/UnitTable_20201122_Cleaned2.csv")

#data cleaning

manual_class = manual_class %>% 
  rename(note_class = label_family)

#bind vector of responses (manually clustered notes)
labelled_data = cbind(acoustic_params, note_class = manual_class$note_class)

#check variables
str(labelled_data)

#Exploratory DA

#check correlations
x = acoustic_params[,4:ncol(acoustic_params)]

res = cor(x) 
res2 = Hmisc::rcorr(as.matrix(x))

source("./functions/flattenCorrMatrix.R")
cor_res = flattenCorrMatrix(res2$r, res2$P)

corrplot::corrplot(res)

heatmap(x = res)
#strong correlation between duration and time measurements. 

#data cleaning ----

#remove non-feature/response columns
labelled_data = labelled_data[,4:ncol(labelled_data)]

labelled_data <- labelled_data %>%
  mutate_if(is.character,factor)

#Partition dataset into training and testing set -----
set.seed(1066)
note_split = initial_split(labelled_data,prop=0.8)
note_train = training (note_split)
note_test = testing (note_split)

#recipe ----
std_recipe <- recipe(note_class ~., data = note_train) %>% #set sweep as response variable. everything else is a predictor.
  step_corr(all_predictors(),threshold = 0.7) %>% #remove all highly correlated predictors
  step_normalize(all_predictors()) %>% #normalize all predictors
  prep()

#RDA ----
# note_rda <- discrim::discrim_regularized(
#   mode = 'classification', 
#   frac_common_cov = tune(), #lambda
#   frac_identity = tune() #gamma
# ) %>%
#   set_engine("klaR")
# 
# rda_grid <- grid_regular(frac_common_cov=discrim::frac_common_cov(range=c(0,1)),
#                          discrim::frac_identity(range=c(0,1)),
#                          levels= 2)
# 
# #hack for weird bug
# names(rda_grid)[1] <- "frac_identity"

# source("./functions/model_tune.R")
# 
# rda_results = model_tune(recipe = std_recipe,
#                         train_data = note_train,
#                         cv_folds = 5,
#                         model = note_rda ,
#                         tuning_params = rda_grid,
#                         seed = 162,
#                         response = "note_class")

library(doParallel)
all_cores <- parallel::detectCores(logical = FALSE)
registerDoParallel(cores = all_cores)

#Random Forest

note_rf<-rand_forest(
  mode="classification",
  mtry=tune(),
  trees= 500, 
  min_n=tune()
) %>%
  set_engine("ranger")

rf_grid<-grid_regular(mtry(range=c(3,11)),min_n(range=c(100,1000)),levels= 3)

source("./functions/model_tune.R")

rf_results = model_tune(recipe = std_recipe,
                        train_data = note_train,
                        cv_folds = 10,
                        model = note_rf ,
                        tuning_params = rf_grid,
                        seed = 162,
                        response = "note_class")

predict(rf_results$fitted_model, note_test, type = "prob") %>% 
  bind_cols(note_test$note_class) %>% 
  glimpse()

results = predict(rf_results$fitted_model, note_test) %>% 
  bind_cols(note_class = note_test$note_class)

#confusion matrix
conf_mat(data = results, truth = note_class, estimate = .pred_class)

#overall testing results
predict(rf_results$fitted_model, note_test, type = "prob") %>% 
  cbind(predict(rf_results$fitted_model, note_test) ) %>%
  cbind(note_class = note_test$note_class) %>% 
  metrics(note_class, .pred_Chip:.pred_Whine, estimate = .pred_class)

#plotting gain curves
rf_probs = predict(rf_results$fitted_model, note_test, type = "prob") %>% 
  bind_cols(note_test)

rf_probs %>%
  gain_curve(note_class, .pred_Chip:.pred_Whine) %>%
  autoplot()

#CV plots

rf_cv = rf_results$tune_tibble

library(ggplot2)

ggplot(data = rf_cv, aes(x = min_n, y = mean)) + 
  geom_point( aes(color = mtry))

#Prediction on computationally excised notes

comp_notes = read_csv("~/work/PhD/Data/note_excision/datasets/comp_data.csv")

correct_labels = note_match(comp_table = comp_notes, manual_table = manual_class)

FP = note_match(manual_class, comp_notes)

comp_notes$note_label = NULL

comp_notes = comp_notes %>%
  cbind(note_class = correct_labels) %>%
  dplyr::filter(note_class != "unclassified")

#percentage of Becky's notes detected computationally
nrow(comp_notes)/nrow(manual_class)

comp_results = predict(rf_results$fitted_model, comp_notes) %>%
  bind_cols(note_class = comp_notes$note_class)

#confusion matrix
conf_mat(data = comp_results, truth = note_class, estimate = .pred_class)
 
#overall testing results
predict(rf_results$fitted_model, comp_notes, type = "prob") %>%
  cbind(predict(rf_results$fitted_model, comp_notes) ) %>%
  cbind(note_class = comp_notes$note_class) %>%
  metrics(note_class, .pred_Chip:.pred_Whine, estimate = .pred_class)

  
#SVM

#~4 hrs to make one SVM......

library(kernlab)

note_svm<-svm_poly(
  mode="classification",
  cost=tune(),
  degree=1
) %>%
  set_engine("kernlab")

svm_grid<-grid_regular(cost(range=c(10,20)),
                       levels= 2,
                       original = T)

svm_results = model_tune(recipe = std_recipe,
                        train_data = note_train,
                        cv_folds = 5,
                        model = note_svm ,
                        tuning_params = svm_grid,
                        seed = 1,
                        response = "note_class")

#xgboost

library(xgboost)
xgboost_model <- boost_tree(
    mode = "classification",
    trees = 100,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune(),
    mtry = tune()
  ) %>%
  set_engine("xgboost", objective = 'multi:softmax')



#specify parameter set
# xgboost_params <- dials::parameters(
#   min_n(),
#   tree_depth(),
#   learn_rate(),
#   loss_reduction(),
#   mtry()
# )

set.seed(52)
xgboost_grid <- dials::grid_latin_hypercube(
  min_n(),
  tree_depth(),
  learn_rate(),
  loss_reduction(), 
  finalize(mtry(), note_train),
  size = 10
)

knitr::kable(head(xgboost_grid))

xgboost_results = model_tune(recipe = std_recipe,
                         train_data = note_train,
                         cv_folds = 5,
                         model = xgboost_model ,
                         tuning_params = xgboost_grid,
                         seed = 1,
                         response = "note_class")

predict(xgboost_results$fitted_model, note_test, type = "prob") %>% 
  cbind(predict(xgboost_results$fitted_model, note_test) ) %>%
  cbind(note_class = note_test$note_class) %>% 
  metrics(note_class, .pred_Chip:.pred_Whine, estimate = .pred_class)

rf_probs = predict(rf_results$fitted_model, note_test, type = "prob") %>% 
  bind_cols(note_test)

rf_probs %>%
  gain_curve(note_class, .pred_Chip:.pred_Whine) %>%
  autoplot()
  
#computational excised notes
comp_notes = read.csv("~/work/PhD/Data/song_features/datasets/measurements.csv")
