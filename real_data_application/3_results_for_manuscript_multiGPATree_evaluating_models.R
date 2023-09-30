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
library(LPM)
registerDoParallel(cores = 2, cl = 2)
# setwd("path_to_data_code") # set working directory to path

# Evaluating the models  ####
## SLE + RA + GS, multiGPATree results ####
load("SLE_RA_mGPATree_GS.rdata")
slerags <- g1
p1 <- multiGPATree::plot(slerags@fit$SLE_RA)

assoc3_sle_ra_gs_mgpatree <- multiGPATree::assoc(slerags@fit$SLE_RA, 
                                                 FDR = 0.05, 
                                                 fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')
table(assoc3_sle_ra_gs_mgpatree$SLE_RA) # joint association
table(assoc3_sle_ra_gs_mgpatree$SLE) # marginal association
table(assoc3_sle_ra_gs_mgpatree$RA) # marginal association

table(assoc3_sle_ra_gs_mgpatree$SLE_RA, assoc3_sle_ra_gs_mgpatree$leaf)
table(assoc3_sle_ra_gs_mgpatree$SLE, assoc3_sle_ra_gs_mgpatree$leaf)
table(assoc3_sle_ra_gs_mgpatree$RA, assoc3_sle_ra_gs_mgpatree$leaf)


## SLE + RA + GS, LPM results ####
load("SLE_RA_LPM_GS.rdata")
LPMfit
LPMfit$beta
LPMfit$R

### to run association mapping for LPM we need to pull in raw data ####
datfinal <- fread("multiGPATree_application_data.csv")
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
LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

# statistical inference for identification of risk SNPs
# posterior1 <- post(LPMData$data[1], X = LPMData$X, id = 1, LPMfit)
# assoc1 <- assoc(posterior1, FDRset = 0.05, fdrControl = "global")
# posterior2 <- post(LPMData$data[2], X = LPMData$X, id = 2, LPMfit)
# str(posterior2)

posterior12 <- LPM::post(LPMData$data[c(1, 2)], X = LPMData$X, id = c(1, 2), LPMfit)
str(posterior12)
assoc3_sle_ra_gs_lpm <- assoc12 <- assoc(posterior12, FDRset = 0.05, fdrControl = "global")
head(assoc12)
table(assoc12$eta.marginal1)
table(assoc12$eta.marginal2)
table(assoc12$eta.joint)

# hypothesis testing of annotation
result_test_beta1 <- test_beta(LPMData$data, X = LPMData$X, id = 1, LPMfit)
result_test_beta1
result_test_beta2 <- test_beta(LPMData$data, X = LPMData$X, id = 2, LPMfit)
result_test_beta2
result_test_beta1$p_value<0.05
result_test_beta2$p_value<0.05
result_test_beta1$p_value<0.05 & result_test_beta2$p_value<0.05

### snps common and unique between mgpatree and LPM ####
table(assoc3_sle_ra_gs_mgpatree$SLE_RA, assoc3_sle_ra_gs_lpm$eta.joint)
table(assoc3_sle_ra_gs_mgpatree$SLE, assoc3_sle_ra_gs_lpm$eta.marginal1)
table(assoc3_sle_ra_gs_mgpatree$RA, assoc3_sle_ra_gs_lpm$eta.marginal2)

## SLE + GS, GPATree results ####
load("SLE_GPATree_GS.rdata")
slegs <- sle
slegspruned <- GPATree::prune(slegs, cp = 0.05)
p2 <- GPATree::plot(slegspruned)
print(p2)
assoc3_sle_gs_gpatree <- assoc3 <- GPATree::assoc(slegspruned, FDR = 0.05, fdrControl = 'global' )%>%
  rownames_to_column(., var = 'rsid')
dim(assoc3_sle_gs_gpatree)
head(assoc3)
table(assoc3$P1)
table(assoc3$P1,assoc3$leaf)

### snps common and unique between mgpatree and gpatree ####
table(assoc3_sle_ra_gs_mgpatree$SLE, assoc3_sle_gs_gpatree$P1)

## RA + GS, GPATree results ####
load("RA_GPATree_GS.rdata")
rags <- ra

ragspruned <- GPATree::prune(rags, cp = 0.05)
p3 <- GPATree::plot(ragspruned)
print(p3)
assoc3_ra_gs_gpatree <- assoc3 <- GPATree::assoc(ragspruned, FDR = 0.05, fdrControl = 'global' )%>%
  rownames_to_column(., var = 'rsid')
dim(assoc3)
table(assoc3$P1)
table(assoc3$P1,assoc3$leaf)

### snps common and unique between mgpatree and gpatree ####
table(assoc3_sle_ra_gs_mgpatree$RA, assoc3_ra_gs_gpatree$P1)

## SLE + RA + GSP, multiGPATree results ####
load("SLE_RA_mGPATree_GSP.rdata")
sleragsp <- g1
p2 <- multiGPATree::plot(sleragsp@fit$SLE_RA)

png(filename = 'Fig6AB_SLE_RA_GS_GSP_treeplot.png',
    width=650,
    height = 240)
par(mfcol=c(1,2))
multiGPATree::plot(slerags@fit$SLE_RA)
multiGPATree::plot(sleragsp@fit$SLE_RA)
dev.off()

assoc3_sle_ra_gsp_mgpatree <- assoc3 <- multiGPATree::assoc(sleragsp@fit$SLE_RA, 
                                                            FDR = 0.05, 
                                                            fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')
dim(assoc3)
head(assoc3)
table(assoc3$SLE_RA)
table(assoc3$SLE)
table(assoc3$RA)

table(assoc3$SLE_RA, assoc3$leaf)
table(assoc3$SLE, assoc3$leaf)
table(assoc3$RA, assoc3$leaf)

## SLE + RA + GSP, LPM results ####
load("SLE_RA_LPM_GSP.rdata")
LPMfit
LPMfit$beta
LPMfit$R
### to run association mapping for LPM we need to pull in raw data ####
pdat <- datfinal %>% 
  tibble::column_to_rownames(., var='rsid') %>%
  dplyr::select(., c(14:23)) %>%
  dplyr::rename('Helper_memory_T' = `GSP_Tcells_helperMemory`,
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
LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

# statistical inference for identification of risk SNPs
# posterior1 <- post(LPMData$data[1], X = LPMData$X, id = 1, LPMfit)
# assoc1 <- assoc(posterior1, FDRset = 0.05, fdrControl = "global")
# posterior2 <- post(LPMData$data[2], X = LPMData$X, id = 2, LPMfit)
# str(posterior2)

posterior12 <- post(LPMData$data[c(1, 2)], X = LPMData$X, id = c(1, 2), LPMfit)
str(posterior12)
assoc3_sle_ra_gsp_lpm <- assoc12 <- assoc(posterior12, FDRset = 0.05, fdrControl = "global")
head(assoc12)
table(assoc12$eta.marginal1)
table(assoc12$eta.marginal2)
table(assoc12$eta.joint)

# hypothesis testing of annotation
result_test_beta1 <- test_beta(LPMData$data, X = LPMData$X, id = 1, LPMfit)
result_test_beta2 <- test_beta(LPMData$data, X = LPMData$X, id = 2, LPMfit)
result_test_beta1$p_value<0.05 & result_test_beta2$p_value<0.05

### snps common and unique between mgpatree and LPM ####
table(assoc3_sle_ra_gsp_mgpatree$SLE_RA, assoc3_sle_ra_gsp_lpm$eta.joint)
table(assoc3_sle_ra_gsp_mgpatree$SLE, assoc3_sle_ra_gsp_lpm$eta.marginal1)
table(assoc3_sle_ra_gsp_mgpatree$RA, assoc3_sle_ra_gsp_lpm$eta.marginal2)


## SLE+GSP GPATree results ####
load("SLE_GPATree_GSP.rdata")
slegsp <- sle
slegsppruned <- GPATree::prune(slegsp, cp = 0.005)
p2 <- GPATree::plot(slegsppruned)
print(p2)
assoc3_sle_gsp_gpatree <- assoc3 <- GPATree::assoc(slegsppruned, FDR = 0.05, fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')

dim(assoc3)
head(assoc3)
table(assoc3$P1)
table(assoc3$P1,assoc3$leaf)

### snps common between mgpatree and gpatree ####
table(assoc3_sle_ra_gsp_mgpatree$SLE, assoc3_sle_gsp_gpatree$P1)


## RA+GSP GPATree results ####
load("RA_GPATree_GSP.rdata")
ragsp <- ra

ragsppruned <- GPATree::prune(ragsp, cp = 0.009)
p3 <- GPATree::plot(ragsppruned)
print(p3)
assoc3_ra_gsp_gpatree <- assoc3 <- GPATree::assoc(ragsppruned, FDR = 0.05, fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')

head(assoc3)
table(assoc3$P1)
table(assoc3$P1,assoc3$leaf)

### snps common between mgpatree and gpatree ####
table(assoc3_sle_ra_gsp_mgpatree$RA, assoc3_ra_gsp_gpatree$P1)


## UC + CD + GS, multiGPATree results ####
load("UC_CD_mGPATree_GS.rdata")
uccdgs <- g1
p1 <- multiGPATree::plot(uccdgs@fit$UC_CD)

assoc3_uc_cd_gs_mgpatree <- assoc3 <- multiGPATree::assoc(uccdgs@fit$UC_CD, 
                                                          FDR = 0.05, 
                                                          fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')

head(assoc3)
table(assoc3$UC_CD)
table(assoc3$UC)
table(assoc3$CD)

table(assoc3$UC_CD, assoc3$leaf)
table(assoc3$UC, assoc3$leaf)
table(assoc3$CD, assoc3$leaf)

## UC + CD + GS, LPM results ####
load("UC_CD_LPM_GS.rdata")
LPMfit
LPMfit$beta
LPMfit$R

### to run association mapping for LPM we need to pull in raw data ####
gmat <- datfinal %>%
  tibble::column_to_rownames(., var='rsid') %>%
  select(c(UC, CD))
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

LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

# statistical inference for identification of risk SNPs
# posterior1 <- post(LPMData$data[1], X = LPMData$X, id = 1, LPMfit)
# assoc1 <- assoc(posterior1, FDRset = 0.05, fdrControl = "global")
# posterior2 <- post(LPMData$data[2], X = LPMData$X, id = 2, LPMfit)
# str(posterior2)

posterior12 <- LPM::post(LPMData$data[c(1, 2)], X = LPMData$X, id = c(1, 2), LPMfit)
str(posterior12)
assoc3_uc_cd_gs_lpm <- assoc12 <- LPM::assoc(posterior12, FDRset = 0.05, fdrControl = "global")
head(assoc12)
table(assoc12$eta.marginal1)
table(assoc12$eta.marginal2)
table(assoc12$eta.joint)

# hypothesis testing of annotation
result_test_beta1 <- test_beta(LPMData$data, X = LPMData$X, id = 1, LPMfit)
result_test_beta2 <- test_beta(LPMData$data, X = LPMData$X, id = 2, LPMfit)
result_test_beta1$p_value<0.05 & result_test_beta2$p_value<0.05

### snps common and unique between mgpatree and LPM ####
table(assoc3_uc_cd_gsp_mgpatree$UC_CD, assoc3_uc_cd_gsp_lpm$eta.joint)
table(assoc3_uc_cd_gsp_mgpatree$UC, assoc3_uc_cd_gsp_lpm$eta.marginal1)
table(assoc3_uc_cd_gsp_mgpatree$CD, assoc3_uc_cd_gsp_lpm$eta.marginal2)


## UC+GS, GPATree results ####
load("UC_GPATree_GS.rdata")
ucgs <- sle
ucgspruned <- GPATree::prune(ucgs, cp = 0.05)
p2 <- GPATree::plot(ucgspruned)
print(p2)
assoc3_uc_gs_gpatree <- assoc3 <- GPATree::assoc(ucgspruned, FDR = 0.05, fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')

### snps common between mgpatree and gpatree ####
table(assoc3_uc_cd_gs_mgpatree$UC, assoc3_uc_gs_gpatree$P1)

## CD+GS, GPATree results ####
load("CD_GPATree_GS.rdata")
cdgs <- ra

cdgspruned <- GPATree::prune(cdgs, cp = 0.05)
p3 <- GPATree::plot(cdgspruned)
print(p3)
assoc3_cd_gs_gpatree <- assoc3 <- GPATree::assoc(cdgspruned, FDR = 0.05, fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')
table(assoc3$P1)
table(assoc3$P1,assoc3$leaf)

### snps common between mgpatree and gpatree ####
table(assoc3_uc_cd_gs_mgpatree$CD, assoc3_cd_gs_gpatree$P1)


## UC + CD + GSP, multiGPATree results ####
load("UC_CD_mGPATree_GSP.rdata")
uccdgsp <- g1
p1 <- multiGPATree::plot(uccdgsp)

assoc3_uc_cd_gsp_mgpatree <- assoc3 <- multiGPATree::assoc(uccdgsp@fit$UC_CD, 
                                                           FDR = 0.05, 
                                                           fdrControl = 'global') %>%
  rownames_to_column(., var = 'rsid')

head(assoc3)
table(assoc3$UC_CD)
table(assoc3$UC)
table(assoc3$CD)

table(assoc3$UC_CD,assoc3$leaf)
table(assoc3$UC,assoc3$leaf)
table(assoc3$CD,assoc3$leaf)


## UC + CD + GSP, LPM results ####
load("UC_CD_LPM_GSP.rdata")
LPMfit
LPMfit$beta
LPMfit$R
### to run association mapping for LPM we need to pull in raw data ####
pdat <- datfinal %>% 
  tibble::column_to_rownames(., var='rsid') %>%
  dplyr::select(., c(14:23)) %>%
  dplyr::rename('Helper_memory_T' = `GSP_Tcells_helperMemory`,
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
LPMData <- list('data' = list(),
                'X' = data.frame(cbind(SNP = 1:nrow(gmat), pdat))) # LPM data setup
for (i in 1:ncol(gmat)) {
  LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gmat), p = gmat[, i]))
  names(LPMData$data)[[i]] <- colnames(gmat)[i]
}

# statistical inference for identification of risk SNPs
# posterior1 <- post(LPMData$data[1], X = LPMData$X, id = 1, LPMfit)
# assoc1 <- assoc(posterior1, FDRset = 0.05, fdrControl = "global")
# posterior2 <- post(LPMData$data[2], X = LPMData$X, id = 2, LPMfit)
# str(posterior2)

posterior12 <- LPM::post(LPMData$data[c(1, 2)], X = LPMData$X, id = c(1, 2), LPMfit)
str(posterior12)
assoc3_uc_cd_gsp_lpm <- assoc12 <- LPM::assoc(posterior12, FDRset = 0.05, fdrControl = "global")
head(assoc12)
table(assoc12$eta.marginal1)
table(assoc12$eta.marginal2)
table(assoc12$eta.joint)

# hypothesis testing of annotation
result_test_beta1 <- test_beta(LPMData$data, X = LPMData$X, id = 1, LPMfit)
result_test_beta2 <- test_beta(LPMData$data, X = LPMData$X, id = 2, LPMfit)
result_test_beta1$p_value<0.05 & result_test_beta2$p_value<0.05

### snps common and unique between mgpatree and LPM ####
table(assoc3_uc_cd_gsp_mgpatree$UC_CD, assoc3_uc_cd_gsp_lpm$eta.joint)
table(assoc3_uc_cd_gsp_mgpatree$UC, assoc3_uc_cd_gsp_lpm$eta.marginal1)
table(assoc3_uc_cd_gsp_mgpatree$CD, assoc3_uc_cd_gsp_lpm$eta.marginal2)


## UC+GSP, GPATree results ####
load("UC_GPATree_GSP.rdata")
ucgsp <- sle
ucgsppruned <- GPATree::prune(ucgsp, cp = 0.005)
p2 <- GPATree::plot(ucgsppruned)
print(p2)
assoc3_uc_gsp_gpatree <- assoc3 <- GPATree::assoc(ucgsppruned, FDR = 0.05, fdrControl = 'global' ) %>%
  rownames_to_column(., var = 'rsid')
head(assoc3)
table(assoc3$P1)
table(assoc3$P1,assoc3$leaf)

### snps common between mgpatree and gpatree ####
table(assoc3_uc_cd_gsp_mgpatree$UC, assoc3_uc_gsp_gpatree$P1)


## CD+GSP GPATree results ####
load("CD_GPATree_GSP.rdata")
cdgsp <- ra

cdgsppruned <- GPATree::prune(cdgsp, cp = 0.002)
p3 <- GPATree::plot(cdgsppruned)
print(p3)
assoc3_cd_gsp_gpatree <- assoc3 <- GPATree::assoc(cdgsppruned, FDR = 0.05, fdrControl = 'global' )%>%
  rownames_to_column(., var = 'rsid')
head(assoc3)
table(assoc3$P1)
table(assoc3$P1,assoc3$leaf)

### snps common between mgpatree and gpatree ####
table(assoc3_uc_cd_gsp_mgpatree$CD, assoc3_cd_gsp_gpatree$P1)
table(assoc3_uc_cd_gsp_mgpatree$UC_CD)

# FIGURE 6, COMBINED ####
jpeg(filename = 'Fig6_SLE_RA_UC_CD_GS_GSP_treeplot.jpeg',
     width=650,
     height = 450,
     quality = 100)
par(mfrow=c(2,2))
multiGPATree::plot(slerags@fit$SLE_RA)
multiGPATree::plot(sleragsp@fit$SLE_RA)
multiGPATree::plot(uccdgs@fit$UC_CD)
multiGPATree::plot(uccdgsp@fit$UC_CD)
dev.off()

# 1. pull in anno data for later analysis ####
# datanno <- fread("/Volumes/Elements/1. MUSC/Dissertation/download from server 1 - 10222021/khatia/dissertation/real_data/SLE_RA_UC_CD/all_chr_RA_UC_CD_SLE_anno.out")
datanno <- fread("/Users/khatiwadaa/Documents/2. MUSC/Dissertation/Aim 2/Results/real data results/all_chr_RA_UC_CD_SLE_anno.out")
names(datanno)
head(datanno$rsid)
datanno$Biotype[1:15]
datanno$SYMBOL[1:15]
datanno$Consequence[1:15]
# comparing findings between multigpatree and gpatree for UC, CD, GSP
table(assoc3_uc_cd_gsp_mgpatree$CD, assoc3_cd_gsp_gpatree$P1)
table(assoc3_uc_cd_gsp_mgpatree$UC, assoc3_uc_gsp_gpatree$P1)

head(assoc3_uc_cd_gsp_mgpatree)

# UC snps discovered by multi-gpa-tree but not by gpa-tree
uc_rsid_list <- assoc3_uc_cd_gsp_mgpatree$rsid[which(assoc3_uc_cd_gsp_mgpatree$UC == 1 & assoc3_uc_gsp_gpatree$P1 == 0)]
head(uc_rsid_list)
length(uc_rsid_list)
sub_uc <- datfinal %>%
  dplyr::filter(rsid %in% uc_rsid_list)
select_anno <- datanno %>%
  dplyr::filter(rsid %in% uc_rsid_list) %>%
  dplyr::select(., c(rsid, SYMBOL, Consequence, Biotype)) 

merged_uc <- merge(select_anno, sub_uc, by ='rsid')%>%
  dplyr::filter(., Biotype == 'protein_coding')
table(merged_uc$SYMBOL)
table(merged_uc$Biotype)

testt1 <- merged_uc %>% filter(., SYMBOL %in% c("MUC19", 'THADA', "CDKAL1", "AGBL4"))
table(testt1$chr, testt1$SYMBOL)
test.table <- as.data.frame(table(merged_uc$chr, merged_uc$SYMBOL))


# CD snps discovered by multi-gpa-tree but not by gpa-tree
uc_rsid_list <- assoc3_uc_cd_gsp_mgpatree$rsid[which(assoc3_uc_cd_gsp_mgpatree$CD == 1 & assoc3_cd_gsp_gpatree$P1 == 0)]
head(uc_rsid_list)
length(uc_rsid_list)
sub_uc <- datfinal %>%
  dplyr::filter(rsid %in% uc_rsid_list)
select_anno <- datanno %>%
  dplyr::filter(rsid %in% uc_rsid_list) %>%
  dplyr::select(., c(rsid, SYMBOL, Consequence, Biotype)) 

merged_uc <- merge(select_anno, sub_uc, by ='rsid')%>%
  dplyr::filter(., Biotype == 'protein_coding')
table(merged_uc$SYMBOL)
table(merged_uc$Biotype)

testt1 <- merged_uc %>% filter(., SYMBOL %in% c("AGBL4", 'CADM2', "USP34", "BANK1"))
table(testt1$chr, testt1$SYMBOL)
test.table <- as.data.frame(table(merged_uc$chr, merged_uc$SYMBOL))
