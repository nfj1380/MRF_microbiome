#Microbiome Analysis
##--------------------------------------------------------------------------

setwd("~/Desktop/Moose microbiome/microbiome")
#need to read in the BIOM file from QIIME

library(biomformat)
library(OTUtable)
library(dplyr)

#-----------------------------------------------------------#
###########------------- Import Data  ------------------###########
#-----------------------------------------------------------#
#if in biom format read the data in this way"

dat <- read_biom("otu_table_rarefied.biom")
str(dat)
#make a matrix
otu_table <- as.data.frame(as.matrix(biom_data(dat)));str(otu_table)


taxonomy <- observation_metadata(dat):str(taxonomy)
metadata <- sample_metadata(dat)
head(otu_table[,1:10])

#note to Nick/Liz this is the updated files from Liz correctly rarified. Will reomove the above to avoid confusion
#if in csv form, import data as perusual

otu_table1 <- read.csv('OTU_table_updated.csv', head=T, row.names=1)

taxonomy1 <- read.csv('OTU_taxonomy_updated.csv')

#remove duplicates at same level. Level options are from column names from taxonomic dataset. No need to do this
#combine_otus('Genus', table, taxonomy)
#save table
#write.csv(otu_table, file = "otuDataTable.csv")
#To make a table containing only phyla with at least 10% abundance
# in any one sample and were observed at any abundance in at least 10% of samples

otuFilter1 <- filter_taxa(otu_table1, abundance=0.01, persistence =5)

#compare multiple filter runs
setdiff(otuFilter, otuFilter1)

#make presence/absence
otuFilter1[otuFilter1>0] <-1
#sum rows
otuFilter1$new <- rowSums(otuFilter1[,1:55] )
#remove all common OTUs
otuFilterNoCommon<-subset(otuFilter1, new < 45)

#remove 'new' sort column
otuFilterNoCommon$new <- NULL

#remove 3 samples from dead moose where no metadata was available
otuFilterNoCommon$X052410EH <- NULL
otuFilterNoCommon$X100309MS <- NULL
otuFilterNoCommon$X011210MM <- NULL

write.csv(otuFilterNoCommon, file = "OTUtableFiltered1.csv")
#-----------------------------------------------------------#
###########-------------MDS  ------------------###########
#-----------------------------------------------------------#
#use nMDS to focus on OTUs that were driving the compostional patterns

library(vegan)
library(vegan3d)
#create a jaccard similarity matrix from the presence/absence data

jdist <- vegdist(t(otuFilterNoCommon), method="jaccard")

#isoMDS function performs non-metric MDS
library(MASS)

#run nMDS
vare.mds_k2 <- isoMDS(jdist, k=2, maxit = 999) #this gets a stress pf 0.22 which is acceptable but high.
vare.mds_k3 <- isoMDS(jdist, k=3, maxit = 999) #can try for other dimensions (k) 3 dimensions lowers stress. 3 dimensions is capturing the data variability

# Check stress using a Shepard plot. You should see a relatively tight linear relationship. 
stress_plot_values<-Shepard(jdist,vare.mds_k3$points)
plot(stress_plot_values, xlab="Original Jaccard Distance", ylab="Distance in NMDS Space", main=sprintf("k=2 stress=%.3f",vare.mds$stress))

#plot
ordiplot3d(vare.mds_k3)# 3d plot - can update this code
ordiplot(vare.mds_k2)#just 2 dimensions

#get species vector correlation coefficent p
vec.sp <-envfit(vare.mds_k3$points, t(otuFilterNoCommon), perm=999)
str(vec.sp)

#extract p values
vec.sp.pval<-as.data.frame(vec.sp$vectors$pvals)
#check r values
vec.sp.r<-as.data.frame(vec.sp$vectors$r); str(vec.sp.r)

#filter OTUs not contibution to the compostional pattern
library(tibble)

#filter based on correlation coefficents
filt0.2 <- vec.sp.r %>% rownames_to_column(var="OTU") %>% filter(vec.sp$vectors$r >=0.2)  
filt0.25 <- vec.sp.r %>% rownames_to_column(var="OTU") %>% filter(vec.sp$vectors$r >=0.25)

str(filt0.25)

OTUdata <- as.data.frame(t(otuFilterNoCommon));str(OTUdata )
OTUReduced <- OTUdata[which(names(OTUdata) %in% filt0.25$OTU)]; str(OTUReduced)

#save the file
write.csv(OTUReduced, file = "OTUReducedJune2019_0.25cutoff.csv")
write.csv(OTUReduced, file = "OTUReducedJune2019.csv")
#can remove OTUs if necessary
#OTUReduced <-OTUReduced[-c(1:3),]; dim(OTUReduced)
#OTUReduced <-OTUReduced[,-41]; dim(OTUReduced);str(OTUReduced)
#-----------------------------------------------------------#
###########-------------Import env data------------------###########
#-----------------------------------------------------------#

env <- read.csv('envAll.csv', head=T, row.names =1);str(env)
env$Sex <- as.factor(env$Sex)
env<- env[ order(row.names(env)), ]
env$Sex <- as.integer(env$Sex)
str(env)

#mege data together
mooseMB <- cbind(OTUReduced, env)
write.csv(mooseMB, file= 'mooseMB.csv')
#-----------------------------------------------------------#
###########-------------MRFcov------------------###########
#-----------------------------------------------------------#
library(MRFcov)

# Remove very rare OTUs that may flag errors
mooseMB$New.CleanUp.ReferenceOTU82 <- NULL

# Extract coordinates columns (for now, these need to be named 'Latitude' and 'Longitude')
Latitude <- mooseMB$Capture.Location.UTM.Easting
mooseMB$Capture.Location.UTM.Easting <- NULL
Longitude <- mooseMB$Capture.Location.UTM.Northing
mooseMB$Capture.Location.UTM.Northing <- NULL
coords <- (data.frame(Latitude = Latitude,
                      Longitude = Longitude))

# Prep covariates for CRF analysis
# Change categorical covariates to factor format
mooseMB$Sex <- as.factor(mooseMB$Sex)
levels(mooseMB$Sex)[1]

analysis.data = mooseMB %>%
  cbind(.,data.frame(model.matrix(~.[,'Sex'],.)[,-1])) %>%
  dplyr::select(-Sex) %>%
  dplyr::rename_all(funs(gsub("\\.|model.matrix", "", .)))

# Scale numeric covariates to unit variance and add the coordinates back in
analysis.data %>%
  dplyr::mutate_at(vars(CaptureDate:BodyCondition),
                   funs(as.vector(scale(.)))) -> analysis.data

# First, compare an nonspatial CRF to a nonspatial MRF (no covariates),
# using the raw coordinates as covariates (scaled to unit variance), to
# see if predictive capacity increases. This function builds both models,
# then splits the data into folds (5 in this case) and uses the predict_MRF function
# to predict occurrences. It then generates a comparative diagnostic plot. We would
# hope that most metrics increase when using the CRF vs the MRF model
analysis.data %>%
  dplyr::bind_cols(data.frame(Latitude = Latitude,
                              Longitude = Longitude)) %>%
  dplyr::mutate_at(vars(Latitude, Longitude),
                   funs(as.vector(scale(.)))) -> analysis.data.nonspatial

moose.mrf <- MRFcov(data = analysis.data.nonspatial[,
                                                    1:42], n_nodes = 42, family = "binomial",
                    n_cores = 3)

cv_MRF_diag_rep_spatial(data = analysis.data[,
                                             1:42], n_nodes = 42, family = "binomial",
                        coords = coords, n_cores = 3, compare_null = T,
                        plot = T, n_fold_runs = 100)

# It seems that the CRF fits the data better. Now we compare a spatial vs
# nonspatial CRF. First, the nonspatial model
moose.crf <- MRFcov(data = analysis.data.nonspatial, 
                    n_nodes = 42, family = 'binomial',
                    n_cores = 3)

# Now run the spatial CRF. Here, the coordinates are used to produce spatial 
# regression splines with the call 
# mgcv::smooth.construct2(object = mgcv::s(Latitude, Longitude,
#                        bs = "gp", k = 5), data = coords, knots = NULL)
# need the mgcv library
if(!require(mgcv)){
  install.packages('mgcv')
}

# also need the pbapply library for nice progress bars
if(!require(pbapply)){
  install.packages('pbapply')
}

#makse sure the following is in the directory
source('spatial_crf.R')
moose.crf.spatial <- spatial_crf(data = analysis.data, 
                                 n_nodes = 42, family = 'binomial', 
                                 coords = coords, n_cores = 3)

# Examine base interactions for the two models
plotMRF_hm(MRF_mod = moose.crf)
plotMRF_hm(MRF_mod = moose.crf.spatial)

# View full coefficient matrix
View(moose.crf$direct_coefs)
View(moose.crf.spatial$direct_coefs)

# View spatial effects (these tend to have large effects, suggesting
# their inclusion in the model is warranted)
moose.crf.spatial$direct_coefs %>%
  dplyr::select(contains('Spatial'))

# Inspect important predictors for each OTU
moose.crf$key_coefs$X322906

# Repeat for the spatial model. These are the important predictors
# AFTER accounting for spatial effects
moose.crf.spatial$key_coefs$X322906

moose.crf$key_coefs$X325706
moose.crf.spatial$key_coefs$X325706

moose.crf$key_coefs$X338145
moose.crf.spatial$key_coefs$X338145
# etc..

# Check how important each covariate is for predicting changing interactions
# mean of covariate absolute effect sizes
cov.imp.mean <- lapply(seq_along(moose.crf$indirect_coefs), function(x){
  graph <- moose.crf$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  round(mean(abs(coefs[which(coefs != 0 )])), 4)
})
names(cov.imp.mean) <- names(moose.crf$indirect_coefs)
cov.imp.mean

cov.imp.mean.spatial <- lapply(seq_along(moose.crf.spatial$indirect_coefs), function(x){
  graph <- moose.crf.spatial$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  round(mean(abs(coefs[which(coefs != 0 )])), 4)
})
names(cov.imp.mean.spatial) <- names(moose.crf.spatial$indirect_coefs)
cov.imp.mean.spatial

# total number of positively-affected interactions for each cov
cov.pos.changes <- lapply(seq_along(moose.crf$indirect_coefs), function(x){
  graph <- moose.crf$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  coefs.nonzero <- coefs[which(coefs != 0 )]
  length(coefs.nonzero[which(coefs.nonzero > 0)])
})
names(cov.pos.changes) <- names(moose.crf$indirect_coefs)
cov.pos.changes

# total number of negatively-affected interactions for each cov
cov.neg.changes <- lapply(seq_along(moose.crf$indirect_coefs), function(x){
  graph <- moose.crf$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  coefs.nonzero <- coefs[which(coefs != 0 )]
  length(coefs.nonzero[which(coefs.nonzero < 0)])
})
names(cov.neg.changes) <- names(moose.crf$indirect_coefs)
cov.neg.changes

#-----------------------------------------------------------#
###########-------------MRFcov bootstrap test------------------###########
#-----------------------------------------------------------#

#have to remove x339838 as its too rare for the bootstrap
analysis.data$X339838<- NULL
booted_MRF <- bootstrap_MRF(data = analysis.data.nonspatial, 
                            family= c('binomial'), 
                            n_nodes = 41, 
                            n_bootstraps = 50, 
                            n_its = 10,
                            n_cores = 3)

plotMRF_hm(MRF_mod = booted_MRF)
save(booted_MRF,file="MooseMB_bootstrapTest.RData")

#look at some OTUs
booted_MRF$mean_key_coefs$X322906
booted_MRF$mean_key_coefs$X325706
# loop function for all OTUs