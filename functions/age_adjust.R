#' Age adjust function
#' 
#' Used to account for the social father's phenotype at recording being potentially different to
#' its phenotype when it tutored the son's. Uses dates from the data and beta coefficients from
#' the model. The function uses the beta coefficient for age, the son's date of birth and the social
#' father's date of recording to compute an adjusted phenotype for the social father. The adjusted 
#' values replace the one's in the original dataset and a new mixed effects model is fitted. This 
#' process iterates until the change in the beta coefficients in the next step is smaller than
#' the input tolerance level. 
#' 
#' This function is NOT standalone and is designed to age adjust for the models in the tempo project.  
#'
#' @param model : mixed effect model regressing social father phenotype with the son's phenotype
#' @param tol : Tolerance level for the convergence of the beta coefficients for age. 
#' @param data : the original birdsong data set to fit the model with Bird.ID, phenotype, social father's phenotype, Birth.Date, sf_Age_Rec, Age_Rec and Clutch(in that order). 
#' @param kin_mat: The kinship matrix associated with the data set. 
#'
#' @return: The final mixed effects model of class lmekin.
#' @export
#'
#' @examples
age_adjust <- function(model, tol, data, kin){
  #check inputs
  if(class(model)!="lmekin"){
    stop("wrong input model type")
  }
  if(is.numeric(tol)==F){
    stop("tol must be numeric scalar")
  }
  
  #change name of the phenotype columns for ease of formulas
  names(data)[2] <- c("y")
  names(data)[3] <- c("sf_y")
  
  
  #extract the beta for log(Age)
  beta = model$coefficients$fixed[3]
  #check the correct coefficient was grabbed
  if(names(beta)!="log(Age_Rec)"){
    stop("Beta coefficient is not for log(Age_Rec). Check correct column in beta was chosen.")
  }
  
  #We now compute age adjusted phenotypes----
  #first we cull the birds with no social fathers, e.g. first generation
  cull_data = na.omit(data)
  
  #the son's were raised at a different time than when the social fathers were recorded
  #since songs can change over time, we want to adjust for the sf at recording sounding different from 
  #when they tutored their sons
  
  
  #check these 4 lines
  #compute time difference in years
  time_diff = as.numeric( (as.Date(cull_data$sf_DOR) - as.Date(cull_data$Birth.Date))/365 )
  
  #compute sf age on son's birthday
  sf_adj_age = cull_data$sf_Age_Rec - time_diff
  #convert difference to log space
  log_diff = log(cull_data$sf_Age_Rec) - log(sf_adj_age)
  #compute age adjusted phenotype, negative sign because when the sf was recorded
  #after son's birthday, sf phenotype at recording would have gone up by beta*diff compared to 
  #son's birthday. If sf was recorded before son's birthday, sf phenotype would be beta*diff
  #lower than what it was at son's birthday so we would add beta*diff
  y_new = cull_data$sf_y - beta*(log_diff)
  
  #create new dataframe with adjusted phenotypes for sf
  adj_data = cull_data
  adj_data$sf_y = y_new
  #adj_data$Bird.ID <- as.factor(adj_data$Bird.ID)
  
  #make a reduced kinship matrix to suit adj_data ----
  ID = adj_data$Bird.ID
  kin.adj=kin[match(ID,rownames(kin)),match(ID,colnames(kin))]
  
  #fit model using adjusted ages as response
  
  new_model = lmekin(formula = y ~ sf_y + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                     #need to make a new kinship matrix for the reduced population
                     varlist = 2*kin.adj,
                     data = adj_data)
  
  new_beta = new_model$coefficients$fixed[3]
  #new_beta = (new_model$coefficients$fixed[3] + beta)/2

  #iterate until the betas converge
  while( abs(new_beta - beta) > tol){
    print(paste("beta", beta))
    print(paste("new_beta", new_beta))
    
    #update beta
    beta = new_beta
    #age adjusted phenotype
    y_new = cull_data$sf_y - beta*(log_diff)
    #update data with new adjusted phenotype
    adj_data$sf_y = y_new
    #refit model
    new_model = lmekin(formula = y ~ sf_y + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                       #need to make a new kinship matrix for the reduced population
                       varlist = 2*kin.adj,
                       data = adj_data)
    new_beta = new_model$coefficients$fixed[3]
    #new_beta = (new_model$coefficients$fixed[3] + beta)/2
  }
  
  #return final model
  return(new_model)
}
