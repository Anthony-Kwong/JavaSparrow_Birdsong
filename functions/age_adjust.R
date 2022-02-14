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
#' This function is NOT standalone and is designed for use in this context only. 
#'
#' @param model : mixed effect model regressing social father phenotype with the son's phenotype
#' @param tol : Tolerance level for the converge of the beta coefficients for age. 
#' @param data : the original birdsong data set to fit the model with Bird.ID, phenotype, social father's phenotype, Age_Rec and Clutch(in that order). 
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
  
  #convert time differences back into years
  sf_age_at_birth = as.numeric( (cull_data$Birth.Date - cull_data$sf_DOB)/365 )
  #compute difference and transform into log space
  adj_age_diff = log(sf_age_at_birth) - log(cull_data$sf_Age_Rec)
    
  #age adjusted phenotype
  y_new = cull_data$sf_y - beta*(adj_age_diff)
  
  #create new dataframe with adjusted ages
  adj_data = cull_data
  adj_data$y = y_new
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
  
  #iterate until the betas converge
  while( abs(new_beta - beta) > tol){
    #update beta
    beta = new_beta
    #age adjusted phenotype
    y_new = cull_data$sf_y - beta*(adj_age_diff)
    #update data with new adjusted phenotype
    adj_data$y = y_new
    #refit model
    new_model = lmekin(formula = y ~ sf_y + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                       #need to make a new kinship matrix for the reduced population
                       varlist = 2*kin.adj,
                       data = adj_data)
    new_beta = new_model$coefficients$fixed[3]
    #print(abs(new_beta - beta))
  }
  
  #return final model
  return(new_model)
}