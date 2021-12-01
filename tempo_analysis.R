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

#clutch is nest they grew up in

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
tempo_df = add_metadata(tempo_df, meta.data, cols = c(5,6))
tempo_df = dplyr::rename(tempo_df, sf_tempo_avg = sf_phenotype)

#fit model

tempo.lme = lmekin(formula = tempo_avg ~ sf_tempo_avg + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*kin.trim, #*2 because kinship halves the correlation
                   data = tempo_df )

tempo.lme2 = lmekin(formula = tempo_avg ~ sf_tempo_avg + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*null.kin,
                   data = tempo_df )


1-pchisq(2*(tempo.lme$loglik - tempo.lme2$loglik),1)


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

#fit model

gap.lme = lmekin(formula = log_avg_gap ~ sf_log_avg_gap + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                   varlist = 2*kin.trim, #*2 because kinship halves the correlation
                   data = gap_df )

gap.lme.null = lmekin(formula = log_avg_gap ~ sf_log_avg_gap + log(Age_Rec) + (1|Bird.ID) + (1|Clutch),
                      varlist = 2*null.kin, #*2 because kinship halves the correlation
                      data = gap_df )

1-pchisq(2*(gap.lme$loglik - gap.lme.null$loglik),1)

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


1-pchisq(2*(std.gap.lme$loglik - std.gap.lme.null$loglik),1)


#avg song duration

#variance/expected value scores


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

