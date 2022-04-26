#This is the script we used to for the tempo project. Here we are investigating the extent which 
#various temporal features may be socially/genetically inherited. We check for social inheritance 
#using the coefficients and p-values of the social father's phenotype in random effects models. 
#We compare the model with the pedigree vs the model without the pedigree to test for the
#genetic effect. (chi square test)

#load packages

library(birdsong.tools)
pacman::p_load(kinship2,vegan,ggfortify,lme4,lmerTest,kinship2,
               coxme,ggplot2,mclust,scatterplot3d,AICcmodavg,warbleR)
colMax <- function(data) sapply(data, max, na.rm = TRUE)
source("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/get.kaz.complexity.R")


#read in Becky's unit table of manually excised notes
# notes = read_csv("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Files for Anthony/JavaSparrow_UnitTable.csv")
# meta.data = read.csv("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Cleaned Files/Java Sparrow Pedigree Info_0911.csv", header = TRUE)

##following code taken from Lewis et. al. 2021

setwd("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Cleaned Files/Final/")

#import manual note classifications
manual.note.classes=read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv",stringsAsFactors=FALSE,fill=TRUE)
#import spectral data on individual notes
spectral.note.data=read.csv("Lewisetal2021_AcousticParameters.csv",header=TRUE)
#import metadata on individual birds
meta.data=read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_metadata.csv")
meta.data$X1 <- NULL
meta.data=meta.data[order(meta.data$Bird.ID),]

#convert data to characters
meta.data$Bird.ID=as.character(meta.data$Bird.ID)
meta.data$Genetic.Father=as.character(meta.data$Genetic.Father)
meta.data$Genetic.Mother=as.character(meta.data$Genetic.Mother)
meta.data$GGF_G=as.character(meta.data$GGF_F)
meta.data$GGM_F=as.character(meta.data$GGM_F)
meta.data$GGF_M=as.character(meta.data$GGF_M)
meta.data$GGM_M=as.character(meta.data$GGM_M)

#create pedigree to use for regressions
#clean the meta.data to replace unknown parents with NAs
cols.w.ors=grep(" or ",meta.data)
for (i in 1:length(cols.w.ors)){
  meta.data[grep(" or ",meta.data[,cols.w.ors[i]]),cols.w.ors[i]]=NA
}
cols.w.ors=grep("nknown",meta.data)
for (i in 1:length(cols.w.ors)){
  meta.data[grep("nknown",meta.data[,cols.w.ors[i]]),cols.w.ors[i]]=NA
}
recover.meta.data=meta.data

#use information from the grandparents columns to extend the pedigree
male.count=sum(!is.na(unique(c(meta.data$Bird.ID,meta.data$Genetic.Father))))
ped.id=c(meta.data$Bird.ID,unique(meta.data$Genetic.Father[which(!(meta.data$Genetic.Father %in% meta.data$Bird.ID) & !is.na(meta.data$Genetic.Father))]),unique(meta.data$Genetic.Mother[which(!(meta.data$Genetic.Mother %in% meta.data$Bird.ID) & !is.na(meta.data$Genetic.Mother))]))
ped.gf=c(meta.data$Genetic.Father,rep(NA,length(ped.id)-length(meta.data$Genetic.Father)))
ped.gm=c(meta.data$Genetic.Mother,rep(NA,length(ped.id)-length(meta.data$Genetic.Mother)))
#use information about grandparents
for (i in 1:length(ped.id)){
  if (ped.id[i] %in% ped.gf){
    if (is.na(ped.gf[i])){
      ped.gf[i]=as.character(meta.data$GGF_F[match(ped.id[i],meta.data$Genetic.Father)])
    }
    if (is.na(ped.gm[i])){
      ped.gm[i]=as.character(meta.data$GGM_F[match(ped.id[i],meta.data$Genetic.Father)])
    }
  }
  if (ped.id[i] %in% ped.gm){
    if (is.na(ped.gf[i])){
      ped.gf[i]=as.character(meta.data$GGF_M[match(ped.id[i],meta.data$Genetic.Mother)])
    }
    if (is.na(ped.gm[i])){
      ped.gm[i]=as.character(meta.data$GGM_M[match(ped.id[i],meta.data$Genetic.Mother)])
    }
  }
}
#include new additions in ped.id
par.list=unique(c(ped.gf,ped.gm))
par.list=par.list[!is.na(par.list)]
new.ids=setdiff(par.list,ped.id)
new.ids=new.ids[which(new.ids!="Unknown")]
new.sex=rep(2,length(new.ids))
new.sex[which(new.ids %in% ped.gf)]=1
ped.id=c(ped.id,new.ids)
ped.gf=c(ped.gf,rep(NA,length(new.ids)))
ped.gm=c(ped.gm,rep(NA,length(new.ids)))

#create the pedigree
full.ped=pedigree(id=ped.id,dadid=ped.gf,momid=ped.gm,sex=c(rep(1,male.count),rep(2,length(ped.id)-male.count-length(new.sex)),new.sex),affected=0+(ped.id %in% manual.note.classes$song_individual),missid="Unknown")
plot.pedigree(full.ped,id=substr(ped.id,4,6),affected=full.ped$affected,cex=.7,align=c(1.5,2))
#build the kinship matrix
kin=kinship(full.ped)
#trim kinship matrix, removing mothers and grandparents we have no recordings of
ID=as.character(meta.data$Bird.ID)
kin.trim=kin[match(ID,rownames(kin)),match(ID,colnames(kin))]
#create the null kinship matrix, this is just the identity
null.kin=diag(nrow(kin.trim))
colnames(null.kin)=colnames(kin.trim)
rownames(null.kin)=rownames(kin.trim)

#save the kinship matrices for future use
saveRDS(kin.trim, file = "~/Documents/GitHub/JavaSparrow_Birdsong/data/kin.trim.rds" )
saveRDS(null.kin, file = "~/Documents/GitHub/JavaSparrow_Birdsong/data/null.kin.rds" )

#clutch is nest they grew up in

#Consistency increases with times produced for all three predictors. We will estimate ages for the six birds with missing ages. This will allow us to estimate their residual consistency (after controlling for note features and age) to use as a response variable, and then to estimate the observed consistency and residual consistency for social fathers when they are raising sons through recursive fitting. 

#Birth dates correlate closely with bird numbers, so I will assign birth dates to match the bird with the most similar number.
meta.data=recover.meta.data
meta.data$Birth.Date=as.character(meta.data$Birth.Date)
meta.data$Birth.Date[which(meta.data$Bird.ID=="JS0002")]="15/09/2007"
meta.data$Birth.Date[which(meta.data$Bird.ID=="JS0005")]="15/09/2007"  #bird not in pedigree or social father
meta.data$Birth.Date[which(meta.data$Bird.ID=="JS0036")]="01/07/2009"
meta.data$Birth.Date[which(meta.data$Bird.ID=="JS0064")]="12/12/2010"  #was father on 20/08/2011
meta.data$Birth.Date[which(meta.data$Bird.ID=="JS0119")]="25/12/2010"
meta.data$Birth.Date[which(meta.data$Bird.ID=="JS0155")]="06/10/2010"
for (i in 1:length(meta.data$Birth.Date)){
  if (substr(meta.data$Birth.Date[i],1,2)=="x/"){
    meta.data$Birth.Date[i]=paste("15",substr(meta.data$Birth.Date[i],2,nchar(meta.data$Birth.Date[i])),sep="")
  }
  if (substr(meta.data$Birth.Date[i],3,5)=="/x/"){
    meta.data$Birth.Date[i]=paste("01/07",substr(meta.data$Birth.Date[i],5,nchar(meta.data$Birth.Date[i])),sep="")
  }
  if (nchar(meta.data$Birth.Date[i])>10){
    meta.data$Birth.Date[i]=substr(meta.data$Birth.Date[i],1,10)
  }
}
meta.data$Birth.Date=as.Date(meta.data$Birth.Date,"%d/%m/%Y")
meta.data$DOR=as.Date(meta.data$DOR,"%d/%m/%Y")

#add missing ages
meta.data$Age_Rec=as.numeric(meta.data$DOR-meta.data$Birth.Date)/365.25

#save processed metadata table
readr::write_csv(meta.data, file = "~/Documents/GitHub/JavaSparrow_Birdsong/data/meta_data.csv")
#---- code from Lewis et. al. 2021 ends here

#tempo analysis----

library(magrittr)

#add tempo
tempo_data = get_tempo(manual.note.classes) %>% 
  dplyr::mutate(ID = sub("-.*", "", sound.files))

#compute average tempo for each bird
bird_ID = unique(tempo_data$ID)

avg_tempo = list()
for(i in 1:length(bird_ID)){
  tar = bird_ID[i]
  ind_data = tempo_data %>%
    dplyr::filter(ID == tar)
  avg_tempo[[i]] = tibble::tibble(Bird.ID = tar,tempo_avg = mean(ind_data$tempo))
}

avg_tempo= do.call(rbind,avg_tempo)

tempo_df = get_sf_phenotype(phenotype_table = avg_tempo, metadata = meta.data, phenotype_index = 2)
#copy meta data across
tempo_df = add_metadata(tempo_df, meta.data, cols = c(3:7))
tempo_df = dplyr::rename(tempo_df, sf_tempo_avg = sf_phenotype)
#get social father's date of recording
tempo_df = get_sf_phenotype(phenotype_table = tempo_df, metadata = meta.data, phenotype_index = 5 )
tempo_df = dplyr::rename(tempo_df, sf_DOR = sf_phenotype)
#get social father's birth dates
tempo_df = get_sf_phenotype(phenotype_table = tempo_df, metadata = meta.data, phenotype_index = 4 )
tempo_df = dplyr::rename(tempo_df, sf_DOB = sf_phenotype)
#get social father's age
tempo_df = get_sf_phenotype(phenotype_table = tempo_df, metadata = meta.data, phenotype_index = 6 )
tempo_df = dplyr::rename(tempo_df, sf_Age_Rec = sf_phenotype)

#fit model

tempo.lme = lmekin(formula = tempo_avg ~ sf_tempo_avg + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*kin.trim, #*2 because kinship halves the correlation
                   data = tempo_df )

tempo.lme2 = lmekin(formula = tempo_avg ~ sf_tempo_avg + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*null.kin,
                   data = tempo_df )


tempo.lme3 = lmekin(formula = tempo_avg ~ sf_tempo_avg + log(Age_Rec) + (1|Bird.ID),
                    varlist = 2*kin.trim,
                    data = tempo_df )

#age adjust example

#assemble small dataset with the right columns for age_adjust()
small_tempo_df <- cbind(tempo_df[,c(1:4)], tempo_df[,c(10,11)], tempo_df[,c(6,7)])
adj_tempo.lme <- age_adjust(model = tempo.lme,tol = 0.05, data = small_tempo_df, kin = kin.trim)
adj_tempo.lme2 <- age_adjust(model = tempo.lme2, tol = 0.05, data = small_tempo_df, kin = null.kin)
adj_tempo.lme3 <- age_adjust(model = tempo.lme3, tol = 0.05, data = small_tempo_df, kin = kin.trim)

#no changes to results


#no effect of clutch, indep of genetics

#coefficients the same, residual error different because some of the variance is hidden when a kinship matrix is included

1-pchisq(2*(tempo.lme$loglik - tempo.lme2$loglik),1)
1-pchisq(2*(tempo.lme$loglik - tempo.lme3$loglik),1) #no effect of clutch, indep of genetics


#test if there's any effect at for both clutch and pedigree
tempo.lme4 = lmekin(formula = tempo_avg ~ sf_tempo_avg + log(Age_Rec) + (1|Bird.ID),
                    varlist = 2*null.kin,
                    data = tempo_df )

1-pchisq(2*(tempo.lme$loglik - tempo.lme4$loglik),1)

#avg gap length ----

gap_data = get_gaps(manual.note.classes) %>% 
  dplyr::mutate(ID = sub("-.*", "", sound.files))

#compute average gap length for each bird
bird_ID = unique(gap_data$ID)

avg_gap = list()
for(i in 1:length(bird_ID)){
  tar = bird_ID[i]
  ind_data = gap_data %>%
    dplyr::filter(ID == tar)
  avg_gap[[i]] <- tibble::tibble(Bird.ID = tar, log_avg_gap = mean( log(ind_data$gap_dur) ))
}

avg_gap= do.call(rbind,avg_gap)
gap_df = get_sf_phenotype(phenotype_table = avg_gap, metadata = meta.data, phenotype_index = 2)
gap_df = add_metadata(gap_df, meta.data, cols = c(5,6))
gap_df = dplyr::rename(gap_df, sf_log_avg_gap = sf_phenotype)

#copy meta data across, Birth.Date and DOR
gap_df = add_metadata(gap_df, meta.data, cols = c(3:4))

#get social father's date of recording
t<- get_sf_phenotype(phenotype_table = gap_df, metadata = meta.data, phenotype_index = 7 )

gap_df = get_sf_phenotype(phenotype_table = gap_df, metadata = meta.data, phenotype_index = 7 )
gap_df = dplyr::rename(gap_df, sf_DOR = sf_phenotype)
#get social father's birth dates
gap_df = get_sf_phenotype(phenotype_table = gap_df, metadata = meta.data, phenotype_index = 6 )
gap_df = dplyr::rename(gap_df, sf_DOB = sf_phenotype)

#get social father's age
gap_df = get_sf_phenotype(phenotype_table = gap_df, metadata = meta.data, phenotype_index = 4 )
gap_df = dplyr::rename(gap_df, sf_Age_Rec = sf_phenotype)

#fit model
gap.lme = lmekin(formula = log_avg_gap ~ sf_log_avg_gap + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*kin.trim, #*2 because kinship halves the correlation
                   data = gap_df )

gap.lme2 = lmekin(formula = log_avg_gap ~ sf_log_avg_gap + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                      varlist = 2*null.kin, #*2 because kinship halves the correlation
                      data = gap_df )

gap.lme3 = lmekin(formula = log_avg_gap ~ sf_log_avg_gap + log(Age_Rec) + (1|Bird.ID),
                 varlist = 2*kin.trim, #*2 because kinship halves the correlation
                 data = gap_df )

gap.lme4 = lmekin(formula = log_avg_gap ~ sf_log_avg_gap + log(Age_Rec) + (1|Bird.ID),
                  varlist = 2*null.kin, #*2 because kinship halves the correlation
                  data = gap_df )


1-pchisq(2*(gap.lme$loglik - gap.lme.null$loglik),1)
1-pchisq(2*(gap.lme$loglik - gap.lme3$loglik),1)
1-pchisq(2*(gap.lme$loglik - gap.lme4$loglik),1)

#check to see if age adjustment changes results

#assemble small dataset with the right columns for age_adjust()
small_gap_df = cbind(gap_df[,1:3], gap_df[,6], gap_df[,9:10], gap_df[,4:5])
adj_gap.lme <- age_adjust(model = tempo.lme,tol = 0.05, data = small_gap_df, kin = kin.trim)
adj_gap.lme2 <- age_adjust(model = tempo.lme2, tol = 0.05, data = small_gap_df, kin = null.kin)
adj_gap.lme3 <- age_adjust(model = tempo.lme3, tol = 0.05, data = small_gap_df, kin = kin.trim)




hist(resid(gap.lme))

#std deviation of gap length ----

std_gap = list()
for(i in 1:length(bird_ID)){
  tar = bird_ID[i]
  ind_data = gap_data %>%
    dplyr::filter(ID == tar)
  std_gap[[i]] <- tibble::tibble(Bird.ID = tar, std_gap = var(ind_data$gap_dur)^0.5)
}

std_gap= do.call(rbind,std_gap)

std_gap_df = get_sf_phenotype(phenotype_table = std_gap, metadata = meta.data, phenotype_index = 2)
std_gap_df = add_metadata(std_gap_df, meta.data, cols = c(5,6))
std_gap_df = dplyr::rename(std_gap_df, sf_std_gap = sf_phenotype)

std.gap.lme = lmekin(formula = std_gap ~ sf_std_gap + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                 varlist = 2*kin.trim, #*2 because kinship halves the correlation
                 data = std_gap_df )

std.gap.lme.null = lmekin(formula = std_gap ~ sf_std_gap + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                          varlist = 2*null.kin, #*2 because kinship halves the correlation
                          data = std_gap_df )

std.gap.lme3 = lmekin(formula = std_gap ~ sf_std_gap + log(Age_Rec) + (1|Bird.ID),
                     varlist = 2*kin.trim, #*2 because kinship halves the correlation
                     data = std_gap_df )

std.gap.lme4 = lmekin(formula = std_gap ~ sf_std_gap + log(Age_Rec) + (1|Bird.ID),
                      varlist = 2*null.kin, #*2 because kinship halves the correlation
                      data = std_gap_df )


1-pchisq(2*(std.gap.lme$loglik - std.gap.lme.null$loglik),1)
1-pchisq(2*(std.gap.lme$loglik - std.gap.lme3$loglik),1)
1-pchisq(2*(std.gap.lme$loglik - std.gap.lme4$loglik),1)


#variance score

#recall our tibble of gaps
gap_data

var_scores = list()
for(i in 1:length(bird_ID)){
  tar = bird_ID[i]
  ind_data = gap_data %>%
    dplyr::filter(ID == tar)
  var_scores[[i]] <- tibble::tibble(Bird.ID = tar, var_score = gap_var_score(ind_data))
}

var_scores = do.call(rbind, var_scores)

varscore_df = get_sf_phenotype(phenotype_table = var_scores, metadata = meta.data, phenotype_index = 2)
varscore_df = add_metadata(varscore_df, meta.data, cols = c(5,6))
varscore_df = dplyr::rename(varscore_df, sf_var_score = sf_phenotype)

varscore.lme = lmekin(formula = var_score ~ sf_var_score + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                     varlist = 2*kin.trim, #*2 because kinship halves the correlation
                     data = varscore_df )

varscore.lme.null = lmekin(formula = var_score ~ sf_var_score + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                          varlist = 2*null.kin, #*2 because kinship halves the correlation
                          data = varscore_df )


1-pchisq(2*(varscore.lme$loglik - varscore.lme.null$loglik),1)

varscore.lme3 = lmekin(formula = var_score ~ sf_var_score + log(Age_Rec) + (1|Bird.ID),
                      varlist = 2*kin.trim, #*2 because kinship halves the correlation
                      data = varscore_df )

varscore.lme4 = lmekin(formula = var_score ~ sf_var_score + log(Age_Rec) + (1|Bird.ID),
                       varlist = 2*null.kin, #*2 because kinship halves the correlation
                       data = varscore_df )

1-pchisq(2*(varscore.lme$loglik - varscore.lme3$loglik),1)
1-pchisq(2*(varscore.lme$loglik - varscore.lme4$loglik),1)


?gap_var_score()

# We want to know: If a social father has more variable gaps, compared to the average. How variable will the son be. 

#More variable, less consistent songs. (investigate performance based gaps, as opposed to deliberate, musical rests)

#score measures internal consistency for each bird

#no explicit comparison with population because that would not give any additional information given
#reasonable biological assumptions

#avg song duration

durations = list()
recordings = unique(manual.note.classes$sound.files)
for(i in 1:length(recordings)){
  tar = recordings[i]
  rec_data = manual.note.classes %>%
    dplyr::filter(sound.files == tar)
  durations[[i]] <- tibble::tibble(Bird.ID = rec_data$song_individual[1], 
                                   dur = get_duration(rec_data))
}

durations = do.call(rbind, durations)

mean_dur = lapply(bird_ID, function(ID){
  ind_data = durations %>%
    dplyr::filter(Bird.ID == ID)
  tibble::tibble(Bird.ID = ID, mean_dur = mean(ind_data$dur))
})

mean_dur = do.call(rbind, mean_dur)

durations_df = get_sf_phenotype(phenotype_table = mean_dur, metadata = meta.data, phenotype_index = 2)
durations_df = add_metadata(durations_df, meta.data, cols = c(5,6))
durations_df = dplyr::rename(durations_df, sf_mean_duration = sf_phenotype)

durations.lme = lmekin(formula = mean_dur ~ sf_mean_duration + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                      varlist = 2*kin.trim, #*2 because kinship halves the correlation
                      data = durations_df )

durations.lme.null = lmekin(formula = mean_dur ~ sf_mean_duration + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                           varlist = 2*null.kin, #*2 because kinship halves the correlation
                           data = durations_df )

durations.lme3 = lmekin(formula = mean_dur ~ sf_mean_duration + log(Age_Rec) + (1|Bird.ID),
                       varlist = 2*kin.trim, #*2 because kinship halves the correlation
                       data = durations_df )


1-pchisq(2*(durations.lme$loglik - durations.lme.null$loglik),1)
1-pchisq(2*(durations.lme$loglik - durations.lme3$loglik),1)

#check relationship between number of notes and song duration----

recordings = unique(manual.note.classes$sound.files)

dur_data = lapply(recordings, function(rec){
  rec_data = manual.note.classes %>%
    dplyr::filter(sound.files == rec)
  dur = get_duration(rec_data)
  note_num = nrow(rec_data)
  id = rec_data$song_individual
  tibble::tibble(note_num = note_num, duration = dur, Bird.ID = id)
})

dur_data = do.call(rbind, dur_data)
dur_data$duration <- as.numeric(dur_data$duration)

ggplot2::ggplot( data = dur_data, aes(note_num,duration)) +
  geom_point()

M1 = lm(duration~note_num, data = dur_data)

#mean dur with mean notes

mean_dur = lapply(bird_ID, function(ID){
  ind_data = dur_data %>%
    dplyr::filter(Bird.ID == ID)
  tibble::tibble(Bird.ID = ID, mean_dur = mean(ind_data$duration),  mean_count = mean(ind_data$note_num))
})

mean_dur = do.call(rbind, mean_dur)

ggplot2::ggplot(data = mean_dur, aes(mean_count,mean_dur)) +
  geom_point()

M2 = lm(mean_count~mean_dur, data = mean_dur)

#mean log scores

#compute vector of mean gaps for each transitition for the whole population
trans = unique(gap_data$transitions)

pop_means = list()
for(i in 1:length(trans)){
  tar = trans[i]
  fil_data = gap_data %>%
    dplyr::filter(transitions == tar)
  pop_means[[i]] = tibble::tibble(
    mean_gap = mean(fil_data$gap_dur), 
    transitions = tar
  )
}

pop_means = do.call(rbind, pop_means)

gap_scores = list()
for(i in 1:length(bird_ID)){
  tar = bird_ID[i]
  ind_data = gap_data %>%
    dplyr::filter(ID == tar)
  gap_scores[[i]] = tibble::tibble(Bird.ID = tar,
                                   gap_score = gap_score(ind_data, pop_means))
}

gap_score_df = do.call(rbind, gap_scores)

gap_score_df = get_sf_phenotype(phenotype_table = gap_score_df, metadata = meta.data, phenotype_index = 2)
gap_score_df = add_metadata(gap_score_df, meta.data, cols = c(5,6))
gap_score_df = dplyr::rename(gap_score_df, sf_gap_score = sf_phenotype)

gapscore.lme = lmekin(formula = gap_score ~ sf_gap_score + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                       varlist = 2*kin.trim, #*2 because kinship halves the correlation
                       data = gap_score_df )

gapscore.lme.null = lmekin(formula = gap_score ~ sf_gap_score + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                            varlist = 2*null.kin, #*2 because kinship halves the correlation
                            data = gap_score_df)

gapscore.lme3 = lmekin(formula = gap_score ~ sf_gap_score + log(Age_Rec) + (1|Bird.ID),
                      varlist = 2*kin.trim, #*2 because kinship halves the correlation
                      data = gap_score_df )


1-pchisq(2*(gapscore.lme$loglik - gapscore.lme.null$loglik),1)
1-pchisq(2*(gapscore.lme$loglik - gapscore.lme3$loglik),1)

## exploratory plots for investigating note duration

plot(x = durations_df$mean_dur, y = durations_df$sf_mean_duration)

durations.lme$residuals
predict(durations.lme)

plot(x = durations.lme$residuals, )

library(coxme)
predict(durations.lme, data = durations_df)

model = lmer(formula = log(mean_dur) ~ log(sf_mean_duration) + log(Age_Rec) + (1|Clutch),
             data = durations_df )


plot(predict(model),resid(model))






#D_s ~ D_f + age + clutch + pedigree 

#gap_len_s ~ clutch + pedigree + transition 
#log mean gap lengths/variances + song duration + tempo  to roughly normalise them. 

#father score, geometric mean of his gap length/population mean gap length (accounting for transition)
#transition type found in one individual.....that will just be 0

#If social father’s gap length for a transition type (e.g. chip to chip) is more
# variable compared to the population mean, his son will also be more
# variable

#father_gap_length

#For every note a father sings, divide by the mean of his A-B notes 

#son's individual gap lengths ~ transition type + sf mean + the rest + position in song

#fractional position plots gap length

