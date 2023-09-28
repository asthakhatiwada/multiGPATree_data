rm(list = ls())
library(multiGPATree)
library(doParallel)
library(foreach)
library(LPM)
library(GPATree)
registerDoParallel(cores = detectCores(), cl = detectCores())
setwd("/Users/khatiwadaa/Documents/2. MUSC/Dissertation/Aim 2/R Code/simulation study")
source('aim2_simulation_study_functions_tree_case.R')

options(digits = 10)

# testing parameters 
M = 100; nGWAS = 2; nAnn = 25
percent_ones_ann = 0.10
percent_overlap_ann = 0.50
percent_overlap_gwas = 1
trueAlphaVec = c(0.4, 0.4)
cpTry = 0.001; initAlpha = 0.1

# simulation study for setting 1 ####
simulated_data <- setting1_aim2_sim_data(M = M,
                                         nGWAS = nGWAS,
                                         nAnn = nAnn,
                                         percent_ones_ann = percent_ones_ann,
                                         percent_overlap_ann = percent_overlap_ann,
                                         percent_overlap_gwas = percent_overlap_gwas,
                                         trueAlphaVec = trueAlphaVec)
sim_res <- aim2_multiGPATree_sim_study(simulated_data = simulated_data, 
                                       percent_ones_ann = percent_ones_ann, 
                                       percent_overlap_ann = percent_overlap_ann, 
                                       percent_overlap_gwas = percent_overlap_gwas,
                                       trueAlphaVec = trueAlphaVec, 
                                       cpTry = cpTry, 
                                       initAlpha = initAlpha,
                                       ncore = detectCores())

write.table(matrix(names(unlist(sim_res)), nrow = 1), sep = ',', 
            'setting1_simstudy_aim2_tree_case_09152023.txt', append = FALSE, col.names = FALSE, row.names = FALSE, quote = TRUE )
t.st.s1 <- proc.time()
multiGPATRee_simulation_results_setting1 <- foreach(nrep = 11:19, .combine = 'cbind') %:%
  foreach(M = 10000, .combine = 'cbind') %:%
  foreach(nGWAS = 2, .combine = 'cbind') %:%
  foreach(nAnn = 25, .combine = 'cbind') %:%
  foreach(percent_ones_ann = c(0.10), .combine = 'cbind') %:% # percent ones = 0.10
  foreach(percent_overlap_ann = c(0.35, 0.50, 0.75), .combine = 'cbind') %:% 
  foreach(percent_overlap_gwas = c(1), .combine = 'cbind') %:% 
  foreach(alpha1 = c(0.4), .combine = 'cbind') %:% # c(0.1, 0.2, 0.3, 0.4, 0.6)
  foreach(alpha2 = c(0.4), .combine = 'cbind') %:% # c(0.6, 0.7, 0.8, 0.9)
  foreach(cpTry = 0.001, .combine = 'cbind') %dopar% {
     # if (alpha1 == alpha2) {

  set.seed(12345 + nrep)
  trueAlphaVec <- c(alpha1, alpha2)
  simulated_data <- setting1_aim2_sim_data(M = M,
                                           nGWAS = nGWAS,
                                           nAnn = nAnn,
                                           percent_ones_ann = percent_ones_ann,
                                           percent_overlap_ann = percent_overlap_ann,
                                           percent_overlap_gwas = percent_overlap_gwas,
                                           trueAlphaVec = trueAlphaVec)
  sim_res <- aim2_multiGPATree_sim_study(simulated_data = simulated_data, 
                                         percent_ones_ann = percent_ones_ann, 
                                         percent_overlap_ann = percent_overlap_ann, 
                                         percent_overlap_gwas = percent_overlap_gwas,
                                         trueAlphaVec = trueAlphaVec, 
                                         cpTry = cpTry, 
                                         initAlpha = initAlpha,
                                         ncore = detectCores())

  write.table(matrix(unlist(sim_res), nrow =1), 
              sep = ',', 
              file = 'setting1_simstudy_aim2_tree_case_09152023.txt',
              append = T, 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = TRUE);
  
  
  return(sim_res)
  
  }
t.end.s1 <- proc.time()
(t.end.s1-t.st.s1)/60


# simulation study for setting 2 ####
simulated_data <- setting2_aim2_sim_data(M = M,
                                         nGWAS = nGWAS,
                                         nAnn = nAnn,
                                         percent_ones_ann = percent_ones_ann,
                                         percent_overlap_ann = percent_overlap_ann,
                                         percent_overlap_gwas = percent_overlap_gwas,
                                         trueAlphaVec = trueAlphaVec)
sim_res <- aim2_multiGPATree_sim_study(simulated_data = simulated_data, 
                                       percent_ones_ann = percent_ones_ann, 
                                       percent_overlap_ann = percent_overlap_ann, 
                                       percent_overlap_gwas = percent_overlap_gwas,
                                       trueAlphaVec = trueAlphaVec, 
                                       cpTry = cpTry, 
                                       initAlpha = initAlpha,
                                       ncore = detectCores())

write.table(matrix(names(unlist(sim_res)), nrow = 1), sep = ',', 
            'setting2_simstudy_aim2_tree_case_09152023.txt', append = FALSE, col.names = FALSE, row.names = FALSE, quote = TRUE )
t.st.s2 <- proc.time()
multiGPATRee_simulation_results_setting2 <- foreach(nrep = 1:50, .combine = 'cbind') %:%
  foreach(M = 10000, .combine = 'cbind') %:%
  foreach(nGWAS = 2, .combine = 'cbind') %:%
  foreach(nAnn = 25, .combine = 'cbind') %:%
  foreach(percent_ones_ann = c(0.10), .combine = 'cbind') %:% # percent ones = 0.10
  foreach(percent_overlap_ann = c(0.35, 0.50, 0.75), .combine = 'cbind') %:% 
  foreach(percent_overlap_gwas = c(1), .combine = 'cbind') %:% 
  foreach(alpha1 = c(0.4), .combine = 'cbind') %:% # c(0.1, 0.2, 0.3, 0.4, 0.6)
  foreach(alpha2 = c(0.4), .combine = 'cbind') %:% # c(0.6, 0.7, 0.8, 0.9)
  foreach(cpTry = 0.001, .combine = 'cbind') %dopar% {
    # if (alpha1 == alpha2) {
    
    set.seed(12345 + nrep)
    trueAlphaVec <- c(alpha1, alpha2)
    simulated_data <- setting2_aim2_sim_data(M = M,
                                             nGWAS = nGWAS,
                                             nAnn = nAnn,
                                             percent_ones_ann = percent_ones_ann,
                                             percent_overlap_ann = percent_overlap_ann,
                                             percent_overlap_gwas = percent_overlap_gwas,
                                             trueAlphaVec = trueAlphaVec)
    sim_res <- aim2_multiGPATree_sim_study(simulated_data = simulated_data, 
                                           percent_ones_ann = percent_ones_ann, 
                                           percent_overlap_ann = percent_overlap_ann, 
                                           percent_overlap_gwas = percent_overlap_gwas,
                                           trueAlphaVec = trueAlphaVec, 
                                           cpTry = cpTry, 
                                           initAlpha = initAlpha,
                                           ncore = detectCores())
    
    write.table(matrix(unlist(sim_res), nrow =1), 
                sep = ',', 
                file = 'setting2_simstudy_aim2_tree_case_09152023.txt',
                append = T, 
                row.names = FALSE, 
                col.names = FALSE, 
                quote = TRUE);
    
    
    return(sim_res)
    
  }


t.end.s2 <- proc.time()
(t.end.s2-t.st.s2x)/60


# simulation study for setting 3 ####
simulated_data <- setting3_aim2_sim_data(M = M,
                                         nGWAS = nGWAS,
                                         nAnn = nAnn,
                                         percent_ones_ann = percent_ones_ann,
                                         percent_overlap_ann = percent_overlap_ann,
                                         percent_overlap_gwas = percent_overlap_gwas,
                                         trueAlphaVec = trueAlphaVec)
sim_res <- aim2_multiGPATree_sim_study(simulated_data = simulated_data, 
                                       percent_ones_ann = percent_ones_ann, 
                                       percent_overlap_ann = percent_overlap_ann, 
                                       percent_overlap_gwas = percent_overlap_gwas,
                                       trueAlphaVec = trueAlphaVec, 
                                       cpTry = cpTry, 
                                       initAlpha = initAlpha,
                                       ncore = detectCores())

write.table(matrix(names(unlist(sim_res)), nrow = 1), sep = ',', 
            'setting3_simstudy_aim2_tree_case_09152023.txt', append = FALSE, col.names = FALSE, row.names = FALSE, quote = TRUE )
t.st.s3 <- proc.time()
multiGPATRee_simulation_results_setting3 <- foreach(nrep = 1:50, .combine = 'cbind') %:%
  foreach(M = 10000, .combine = 'cbind') %:%
  foreach(nGWAS = 2, .combine = 'cbind') %:%
  foreach(nAnn = 25, .combine = 'cbind') %:%
  foreach(percent_ones_ann = c(0.10), .combine = 'cbind') %:% # percent ones = 0.10
  foreach(percent_overlap_ann = c(0.35, 0.50, 0.75), .combine = 'cbind') %:% 
  foreach(percent_overlap_gwas = c(1), .combine = 'cbind') %:% 
  foreach(alpha1 = c(0.4), .combine = 'cbind') %:% # c(0.1, 0.2, 0.3, 0.4, 0.6)
  foreach(alpha2 = c(0.4), .combine = 'cbind') %:% # c(0.6, 0.7, 0.8, 0.9)
  foreach(cpTry = 0.001, .combine = 'cbind') %dopar% {
    # if (alpha1 == alpha2) {
    
    set.seed(12345 + nrep)
    trueAlphaVec <- c(alpha1, alpha2)
    simulated_data <- setting3_aim2_sim_data(M = M,
                                             nGWAS = nGWAS,
                                             nAnn = nAnn,
                                             percent_ones_ann = percent_ones_ann,
                                             percent_overlap_ann = percent_overlap_ann,
                                             percent_overlap_gwas = percent_overlap_gwas,
                                             trueAlphaVec = trueAlphaVec)
    sim_res <- aim2_multiGPATree_sim_study(simulated_data = simulated_data, 
                                           percent_ones_ann = percent_ones_ann, 
                                           percent_overlap_ann = percent_overlap_ann, 
                                           percent_overlap_gwas = percent_overlap_gwas,
                                           trueAlphaVec = trueAlphaVec, 
                                           cpTry = cpTry, 
                                           initAlpha = initAlpha,
                                           ncore = detectCores())
    
    write.table(matrix(unlist(sim_res), nrow =1), 
                sep = ',', 
                file = 'setting3_simstudy_aim2_tree_case_09152023.txt',
                append = T, 
                row.names = FALSE, 
                col.names = FALSE, 
                quote = TRUE);
    
    
    return(sim_res)
    
  }


t.end.s3 <- proc.time()
(t.end.s3-t.st.s3)/60




