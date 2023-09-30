rm(list = ls())
library(ggplot2)
library(data.table)
library(ggpubr)
library(multiGPATree)
library(doParallel)
library(foreach)
library(qqman)
library(tibble)
library(dplyr)
library(tidyr)
registerDoParallel(cores = 2, cl = 2)
# setwd("path_to_data_code") # set working directory to path

# read in the data 
datfinal <- fread("multiGPATree_application_data.csv")

# 1. models for SLE + RA using GS annotations ####
names(datfinal)
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(SLE, RA))
head(gmat)

pdat <- datfinal %>% 
  tibble::column_to_rownames(., var='rsid') %>%
  dplyr::select(., c("GS_Blood","GS_Brain","GS_Epithel","GS_GI","GS_Heart","GS_Lung","GS_Muscle")) %>%
  dplyr::rename(Blood = GS_Blood,
         Brain = GS_Brain,
         Epithelium = GS_Epithel,
         GI = GS_GI,
         Heart = GS_Heart,
         Lung = GS_Lung,
         Muscle = GS_Muscle)
head(pdat)


## joint analysis of SLE + RA using GS annotations, multiGPATREE ####
stime <- proc.time()
g1 <- multiGPATree::multiGPATree(gwasPval = gmat,
                                  annMat = pdat,
                                  initAlpha = 0.1, 
                                  cpTry = 0.0000001)
etime <- proc.time()
mgpatree_sle_ra_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
mgpatree_sle_ra_gs_comptime/60 # mins
save(g1, file = 'SLE_RA_mGPATree_GS.rdata')
# load('SLE_RA_mGPATree_GS_cpTry0.0000001.rdata')

## joint analysis of SLE + RA using GS annotations, LPM ####
library(pbivnorm)
library(LPM)
LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

stime <- proc.time()
bLPMfit <- LPM::bLPM(data = LPMData$data, X = LPMData$X)
LPMfit <- LPM::LPM(bLPMfit)
etime <- proc.time()
LPM_sle_ra_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
LPM_sle_ra_gs_comptime/60 # mins
save(LPMfit, file = 'SLE_RA_LPM_GS.rdata')

## SLE only model using GS annotations, GPATree####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(SLE))
head(gmat)
stime <- proc.time()
sle <- GPATree::GPATree(gwasPval = gmat, 
                        annMat = pdat,
                        initAlpha = 0.1, 
                        cpTry = 0.0000001)
etime <- proc.time()
gpatree_sle_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_sle_gs_comptime/60
save(sle, file = 'SLE_GPATree_GS.rdata')


## RA only model using GS annotations, GPATree####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(RA))
head(gmat)
stime <- proc.time()
ra <- GPATree::GPATree(gwasPval = gmat,
                       annMat = pdat,
                       initAlpha = 0.1, 
                       cpTry = 0.0000001)
etime <- proc.time()
gpatree_ra_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_ra_gs_comptime/60
save(ra, file = 'RA_GPATree_GS.rdata')

# 2. Models for SLE + RA using GSP annotations ####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(SLE, RA))
head(gmat)
names(datfinal)
pdat <- datfinal %>% 
  tibble::column_to_rownames(., var='rsid') %>%
  select(., c(14:23)) %>%
  rename('Helper_memory_T' = `GSP_Tcells_helperMemory`,
         'Helper_naive_T' = `GSP_Tcells_helperNaive`,
         'Effector_memory_enriched_T' =  `GSP_Tcells_effectorMemory`,
         'Regulatory_T' = `GSP_Tcells_regulatory`,
         'CD8plus_naive_T' = `GSP_Tcells_CD8plusNaive`,             
         'CD8plus_memory_T' = `GSP_Tcells_CD8plusMemory`,
         'Monocytes' = `GSP_monocytes`,                     
         'Neutrophils' = `GSP_neutrophils`,
         'Primary_B' = `GSP_Bcells`,
         'Natural_killers' = `GSP_NaturalKiller`) 
head(pdat)
names(pdat)
names(gmat)

## joint analysis of SLE + RA using GSP annotations, multiGPATREE ####
stime <- proc.time()
g1 <- multiGPATree::multiGPATree(gwasPval = gmat,
                                  annMat = pdat,
                                  initAlpha = 0.1, 
                                  cpTry = 0.0000001)
etime <- proc.time()
mgpatree_sle_ra_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
mgpatree_sle_ra_gsp_comptime/60 #mins
save(g1, file = 'SLE_RA_mGPATree_GSP.rdata')

## joint analysis of SLE + RA using GSP annotations, LPM ####
library(pbivnorm)
library(LPM)
LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

stime <- proc.time()
bLPMfit <- LPM::bLPM(data = LPMData$data, X = LPMData$X)
LPMfit <- LPM::LPM(bLPMfit)
etime <- proc.time()
LPM_sle_ra_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
LPM_sle_ra_gsp_comptime/60 # mins
save(LPMfit, file = 'SLE_RA_LPM_GSP.rdata')


## SLE only model using GSP annotations####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(SLE))
stime <- proc.time()
sle <- GPATree::GPATree(gwasPval = gmat, 
                        annMat = pdat,
                        initAlpha = 0.1, 
                        cpTry = 0.0000001)
etime <- proc.time()
gpatree_sle_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_sle_gsp_comptime/60 # mins
save(sle, file = 'SLE_GPATree_GSP.rdata')


## RA only model using GSP annotations####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(RA))
head(gmat)
stime <- proc.time()
ra <- GPATree::GPATree(gwasPval = gmat,
                       annMat = pdat,
                       initAlpha = 0.1, 
                       cpTry = 0.0000001)
etime <- proc.time()
gpatree_ra_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_ra_gsp_comptime/60 # mins
save(ra, file = 'RA_GPATree_GSP.rdata')


# 3. Models for UC + CD using GS annotations ####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(UC, CD))
head(gmat)

pdat <- datfinal %>% 
  tibble::column_to_rownames(., var='rsid') %>%
  select(., c("GS_Blood","GS_Brain","GS_Epithel","GS_GI","GS_Heart","GS_Lung","GS_Muscle")) %>%
  rename(Blood = GS_Blood,
         Brain = GS_Brain,
         Epithelium = GS_Epithel,
         GI = GS_GI,
         Heart = GS_Heart,
         Lung = GS_Lung,
         Muscle = GS_Muscle)
head(pdat)


## joint analysis of UC + CD using GS annotations, multiGPATREE ####
stime <- proc.time()
g1 <- multiGPATree::multiGPATree(gwasPval = gmat,
                                  annMat = pdat,
                                  initAlpha = 0.1, 
                                  cpTry = 0.0000001)
etime <- proc.time()
mgpatree_uc_cd_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
mgpatree_uc_cd_gs_comptime/60 # mins

save(g1, file = 'UC_CD_mGPATree_GS.rdata')

## joint analysis of UC + CD using GS annotations, LPM ####
library(pbivnorm)
library(LPM)
LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

stime <- proc.time()
bLPMfit <- LPM::bLPM(data = LPMData$data, X = LPMData$X)
LPMfit <- LPM::LPM(bLPMfit)
etime <- proc.time()
LPM_uc_cd_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
LPM_uc_cd_gs_comptime/60 # mins
save(LPMfit, file = 'UC_CD_LPM_GS.rdata')


## UC only model using GS annotations####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(UC))
stime <- proc.time()
sle <- GPATree::GPATree(gwasPval = gmat, 
                        annMat = pdat,
                        initAlpha = 0.1, 
                        cpTry = 0.0000001)
etime <- proc.time()
gpatree_uc_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_uc_gs_comptime/60 # mins

save(sle, file = 'UC_GPATree_GS.rdata')


## CD only model using GS annotations####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(CD))
head(gmat)
stime <- proc.time()
ra <- GPATree::GPATree(gwasPval = gmat,
                       annMat = pdat,
                       initAlpha = 0.1, 
                       cpTry = 0.0000001)
etime <- proc.time()
gpatree_cd_gs_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_cd_gs_comptime/60 # mins

save(ra, file = 'CD_GPATree_GS.rdata')

# 4. Models for UC + CD using GSP annotations ####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(UC, CD))
head(gmat)

pdat <- datfinal %>% 
  tibble::column_to_rownames(., var='rsid') %>%
  select(., c(14:23)) %>%
  rename('Helper_memory_T' = `GSP_Tcells_helperMemory`,
         'Helper_naive_T' = `GSP_Tcells_helperNaive`,
         'Effector_memory_enriched_T' =  `GSP_Tcells_effectorMemory`,
         'Regulatory_T' = `GSP_Tcells_regulatory`,
         'CD8plus_naive_T' = `GSP_Tcells_CD8plusNaive`,             
         'CD8plus_memory_T' = `GSP_Tcells_CD8plusMemory`,
         'Monocytes' = `GSP_monocytes`,                     
         'Neutrophils' = `GSP_neutrophils`,
         'Primary_B' = `GSP_Bcells`,
         'Natural_killers' = `GSP_NaturalKiller`)
head(pdat)
names(pdat)

## joint analysis of UC + CD using GSP annotations, multiGPATREE ####
stime <- proc.time()
g1 <- multiGPATree::multiGPATree(gwasPval = gmat,
                                  annMat = pdat,
                                  initAlpha = 0.1, 
                                  cpTry = 0.0000001)
etime <- proc.time()
mgpatree_uc_cd_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
mgpatree_uc_cd_gsp_comptime/60 # mins
save(g1, file = 'UC_CD_mGPATree_GSP.rdata')

## joint analysis of UC + CD using GSP annotations, LPM ####
library(pbivnorm)
library(LPM)
LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

stime <- proc.time()
bLPMfit <- LPM::bLPM(data = LPMData$data, X = LPMData$X)
LPMfit <- LPM::LPM(bLPMfit)
etime <- proc.time()
LPM_uc_cd_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
LPM_uc_cd_gsp_comptime/60 # mins
save(LPMfit, file = 'UC_CD_LPM_GSP.rdata')


## UC only model using GSP annotations####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(UC))
stime <- proc.time()
sle <- GPATree::GPATree(gwasPval = gmat, 
                        annMat = pdat,
                        initAlpha = 0.1, 
                        cpTry = 0.0000001)
etime <- proc.time()
gpatree_uc_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_uc_gsp_comptime/60 # mins
save(sle, file = 'UC_GPATree_GSP.rdata')


## CD only model using GSP annotations####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(CD))
head(gmat)
stime <- proc.time()
ra <- GPATree::GPATree(gwasPval = gmat,
                       annMat = pdat,
                       initAlpha = 0.1, 
                       cpTry = 0.0000001)
etime <- proc.time()
gpatree_cd_gsp_comptime <- as.numeric(etime[[3]] - stime[[3]])
gpatree_cd_gsp_comptime/60 # mins
save(ra, file = 'CD_GPATree_GSP.rdata')



