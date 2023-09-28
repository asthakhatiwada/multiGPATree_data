# 0. data simulation ####
## setting 1: data simulation ####
setting1_aim2_sim_data <- function(M, nGWAS, nAnn, percent_ones_ann, percent_overlap_ann, percent_overlap_gwas, trueAlphaVec){
  
  # M = 1000
  # nGWAS = 2
  # nAnn = 10
  # percent_ones_ann = 0.10 # 10, 15, 20
  # percent_overlap_ann = 0.35 # >=35 percent
  # percent_overlap_gwas = 1 # >= 25, 50, 75
  # trueAlphaVec = c(0.2, 0.2)
  
  
  # start here
  
  freq_ones <- ceiling(percent_ones_ann * M)
  freq_overlap_ann <- ceiling(percent_overlap_ann * freq_ones)
  freq_overlap_gwas <- ceiling(percent_overlap_gwas * freq_overlap_ann)
  A1 <- c( rep(1, freq_ones), 
           rep(0, M-freq_ones) )
  A2 <- c( rep(0, freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(2*freq_ones-freq_overlap_ann)) )
  A3 <- c( rep(0, 2*freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(3*freq_ones-freq_overlap_ann)) )
  A4 <- c( rep(0, 3*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(4*freq_ones-2*freq_overlap_ann)) )
  A5 <- c( rep(0, 4*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones-freq_overlap_gwas), 
           rep(0, M - (5*freq_ones-2*freq_overlap_ann-freq_overlap_gwas) - freq_overlap_gwas ), 
           rep(1, freq_overlap_gwas ) )
  A6 <- c( rep(0, 5*freq_ones-3*freq_overlap_ann ), 
           rep(1, freq_ones-freq_overlap_gwas),
           rep(0, M-(6*freq_ones-3*freq_overlap_ann-freq_overlap_gwas)-freq_overlap_gwas ),
           rep(1, freq_overlap_gwas) )
  A <- cbind(A1, A2, A3, A4, A5, A6)
  # heatmap(as.matrix(A), Colv = NA, Rowv = NA, scale = 'column')
  
  if (nAnn >6){
    for (i in 7:nAnn) {
      prop_bin <- runif(M, 0.1, 0.3)
      new_ann <- rbinom(M, 1, prop_bin)
      A <- cbind(A, new_ann)
      colnames(A)[i] <- paste('A',i,sep = '')
    }
  }
  
  # generating latent Z
  
  binaryList <- vector( "list", nGWAS )
  for ( k in 1:nGWAS ) {
    binaryList[[k]] <- c( 0, 1 )
  }
  
  binaryMat <- as.matrix(expand.grid( binaryList ))
  nComp <- nrow(binaryMat)
  combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
  Zmat <- matrix(0, M, nComp)
  colnames(Zmat) <- combVec
  
  Z <- rep(0, M)
  Z[A1 == 1 & A2 == 1] <- 1
  Z[A3 == 1 & A4 == 1] <- 2
  Z[(M-freq_overlap_gwas+1):M] <- 3
  
  for (i in 1:ncol(Zmat)) {
    Zmat[Z==i-1, i] <- 1
  }
  
  
  alpha <- trueAlphaVec
  gwasPval <- matrix(runif(nGWAS*M, 0, 1), nrow = M, ncol = nGWAS)
  colnames(gwasPval) <- paste('P', 1:nGWAS, sep = '')
  gwasPval[ which(Z==1 | Z==3), 1 ] <- rbeta( length(which(Z==1 | Z==3)), alpha[1], 1) 
  gwasPval[ which(Z==2 | Z==3), 2 ] <- rbeta( length(which(Z==2 | Z==3)), alpha[2], 1) 
  rownames(A) <- rownames(gwasPval) <- paste('SNP_', 1:M, sep = '')
  
  Ztrue <- matrix(0, nrow = M, ncol = nGWAS)
  Ztrue[Z==1 | Z==3, 1] <- 1
  Ztrue[Z==2 | Z==3, 2] <- 1
  colnames(Ztrue) <- colnames(gwasPval)
  
  return(list(annMat = as.matrix(A),
              gwasPval = gwasPval,
              Ztrue = Ztrue,
              zMat = Zmat
  ))
  
}

## setting 2: data simulation ####
setting2_aim2_sim_data <- function(M, nGWAS, nAnn, percent_ones_ann, percent_overlap_ann, percent_overlap_gwas, trueAlphaVec){
  
  # M = 1000
  # nGWAS = 2
  # nAnn = 10
  # percent_ones_ann = 0.10 # 10, 15, 20
  # percent_overlap_ann = 0.5 # >=35 percent
  # percent_overlap_gwas = 1 # >= 25, 50, 75
  # trueAlphaVec = c(0.2, 0.2)
  # M*0.1*0.5
  
  # start here
  
  freq_ones <- ceiling(percent_ones_ann * M);freq_ones
  freq_overlap_ann <- ceiling(percent_overlap_ann * freq_ones);freq_overlap_ann
  freq_overlap_gwas <- ceiling(percent_overlap_gwas * freq_overlap_ann)
  a1anda2 <- ceiling(freq_overlap_ann*2/3);a1anda2
  a1anda5 <- freq_overlap_ann-a1anda2;a1anda5
  a5anda6 <- a1anda2
  A1 <- c( rep(1, freq_ones-a1anda5),
           rep(0, M-(freq_ones)), 
           rep(1, a1anda5))
  # length(A1)
  A2 <- c( rep(0, freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(2*freq_ones-freq_overlap_ann)) 
  )
  # length(A2)
  # table(A1)
  # table(A2)
  # table(A1, A2)
  A3 <- c( rep(0, 2*freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(3*freq_ones-freq_overlap_ann)) 
  )
  # table(A3)
  A4 <- c( rep(0, 3*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(4*freq_ones-2*freq_overlap_ann)) 
  )
  # length(A4)
  # table(A3,A4)
  A5 <- c( rep(0, 4*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones-a1anda5), 
           rep(0, M - (5*freq_ones-2*freq_overlap_ann)), 
           rep(1, a1anda5) 
  )
  # table(A5)
  # table(A1,A5)
  # a1anda5
  
  A6 <- c( rep(0, 5*freq_ones-3*freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M - (6*freq_ones-3*freq_overlap_ann))
  )
  # length(A6)
  # table(A6)
  # table(A5, A6)
  # table(A3, A4)
  # table(A1,A5)
  
  A <- cbind(A1, A2, A3, A4, A5,A6)
  # heatmap(as.matrix(A), Colv = NA, Rowv = NA, scale = 'column')
  
  if (nAnn >6){ # (nAnn >6)
    for (i in 7:nAnn) {
      prop_bin <- runif(M, 0.1, 0.3)
      new_ann <- rbinom(M, 1, prop_bin)
      A <- cbind(A, new_ann)
      colnames(A)[i] <- paste('A',i,sep = '')
    }
  }
  
  # generating latent Z
  
  binaryList <- vector( "list", nGWAS )
  for ( k in 1:nGWAS ) {
    binaryList[[k]] <- c( 0, 1 )
  }
  
  binaryMat <- as.matrix(expand.grid( binaryList ))
  nComp <- nrow(binaryMat)
  combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
  Zmat <- matrix(0, M, nComp)
  colnames(Zmat) <- combVec
  
  Z <- rep(0, M)
  Z[A1 == 1 & A2 == 1] <- 1
  Z[A3 == 1 & A4 == 1] <- 2
  Z[A1 == 1 & A5 == 1 | A5 == 1 & A6 == 1] <- 3
  # Z[(M-freq_overlap_gwas+1):M] <- 3
  
  for (i in 1:ncol(Zmat)) {
    Zmat[Z==i-1, i] <- 1
  }
  # Zmat[(M-freq_overlap_gwas+1):M, 2] <- 1
  
  # testmat <- as.data.frame(cbind(A,Zmat))
  # heatmap(as.matrix(testmat), Colv = NA, Rowv = NA, scale = 'column')
  
  alpha <- trueAlphaVec
  gwasPval <- matrix(runif(nGWAS*M, 0, 1), nrow = M, ncol = nGWAS)
  colnames(gwasPval) <- paste('P', 1:nGWAS, sep = '')
  gwasPval[ which(Z==1 | Z==3), 1 ] <- rbeta( length(which(Z==1 | Z==3)), alpha[1], 1) 
  gwasPval[ which(Z==2 | Z==3), 2 ] <- rbeta( length(which(Z==2 | Z==3)), alpha[2], 1) 
  rownames(A) <- rownames(gwasPval) <- paste('SNP_', 1:M, sep = '')
  # testmat <- as.data.frame(cbind(A,Zmat, gwasPval))
  # heatmap(as.matrix(testmat), Colv = NA, Rowv = NA, scale = 'column')
  
  Ztrue <- matrix(0, nrow = M, ncol = nGWAS)
  Ztrue[Z==1 | Z==3, 1] <- 1
  Ztrue[Z==2 | Z==3, 2] <- 1
  colnames(Ztrue) <- colnames(gwasPval)
  
  return(list(annMat = as.matrix(A),
              gwasPval = gwasPval,
              Ztrue = Ztrue,
              zMat = Zmat
  ))
  
}

## setting 3: data simulation ####
setting3_aim2_sim_data <- function(M, nGWAS, nAnn, percent_ones_ann, percent_overlap_ann, percent_overlap_gwas, trueAlphaVec){
  
  # M = 1000
  # nGWAS = 2
  # nAnn = 10
  # percent_ones_ann = 0.10 # 10%
  # percent_overlap_ann = 0.5 # >=35 percent
  # percent_overlap_gwas = 1 # =100%
  # trueAlphaVec = c(0.2, 0.2)
  # M*0.1*0.5
  
  # start here
  
  freq_ones <- ceiling(percent_ones_ann * M);freq_ones
  freq_overlap_ann <- ceiling(percent_overlap_ann * freq_ones);freq_overlap_ann
  freq_overlap_gwas <- ceiling(percent_overlap_gwas * freq_overlap_ann)
  # a1anda2 <- ceiling(freq_overlap_ann*1/2);a1anda2
  # a1anda5 <- freq_overlap_ann-a1anda2;a1anda5
  A1 <- c( rep(1, freq_ones),
           rep(0, 3*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_overlap_ann),
           rep(0, M-(4*freq_ones-freq_overlap_ann))
  )
  # length(A1)
  # table(A1)
  A2 <- c( rep(0, freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(2*freq_ones-freq_overlap_ann)) 
  )
  # length(A2)
  # table(A1)
  # table(A2)
  # table(A1, A2)
  A3 <- c( rep(0, 2*freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(3*freq_ones-freq_overlap_ann)) 
  )
  # table(A3)
  A4 <- c( rep(0, 3*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(4*freq_ones-2*freq_overlap_ann)) 
  )
  # length(A4)
  # table(A3,A4)
  A5 <- c( rep(0, 4*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M - (5*freq_ones-2*freq_overlap_ann))
  )
  # table(A5)
  # table(A1,A5)
  # table(A3, A4)
  # table(A1,A5)
  # table(A1,A2)
  
  A <- cbind(A1, A2, A3, A4, A5)
  # heatmap(as.matrix(A), Colv = NA, Rowv = NA, scale = 'column')
  
  if (nAnn >=6){ # (nAnn >6)
    for (i in 6:nAnn) {
      prop_bin <- runif(M, 0.1, 0.3)
      new_ann <- rbinom(M, 1, prop_bin)
      A <- cbind(A, new_ann)
      colnames(A)[i] <- paste('A',i,sep = '')
    }
  }
  
  # generating latent Z
  
  binaryList <- vector( "list", nGWAS )
  for ( k in 1:nGWAS ) {
    binaryList[[k]] <- c( 0, 1 )
  }
  
  binaryMat <- as.matrix(expand.grid( binaryList ))
  nComp <- nrow(binaryMat)
  combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
  Zmat <- matrix(0, M, nComp)
  colnames(Zmat) <- combVec
  
  Z <- rep(0, M)
  Z[A1 == 1 & A2 == 1] <- 1
  Z[A3 == 1 & A4 == 1] <- 2
  Z[A1 == 1 & A5 == 1] <- 3
  # Z[(M-freq_overlap_gwas+1):M] <- 3
  
  for (i in 1:ncol(Zmat)) {
    Zmat[Z==i-1, i] <- 1
  }
  # Zmat[(M-freq_overlap_gwas+1):M, 2] <- 1
  
  # testmat <- as.data.frame(cbind(A,Zmat))
  # heatmap(as.matrix(testmat), Colv = NA, Rowv = NA, scale = 'column')
  
  alpha <- trueAlphaVec
  gwasPval <- matrix(runif(nGWAS*M, 0, 1), nrow = M, ncol = nGWAS)
  colnames(gwasPval) <- paste('P', 1:nGWAS, sep = '')
  gwasPval[ which(Z==1 | Z==3), 1 ] <- rbeta( length(which(Z==1 | Z==3)), alpha[1], 1) 
  gwasPval[ which(Z==2 | Z==3), 2 ] <- rbeta( length(which(Z==2 | Z==3)), alpha[2], 1) 
  rownames(A) <- rownames(gwasPval) <- paste('SNP_', 1:M, sep = '')
  # testmat <- as.data.frame(cbind(A,Zmat, gwasPval))
  # heatmap(as.matrix(testmat), Colv = NA, Rowv = NA, scale = 'column')
  
  Ztrue <- matrix(0, nrow = M, ncol = nGWAS)
  Ztrue[Z==1 | Z==3, 1] <- 1
  Ztrue[Z==2 | Z==3, 2] <- 1
  colnames(Ztrue) <- colnames(gwasPval)
  
  return(list(annMat = as.matrix(A),
              gwasPval = gwasPval,
              Ztrue = Ztrue,
              zMat = Zmat
  ))
  
}

# 1. function to get performance statistic for multiGPATree aim 2 ####
perf_summary_func <- function(Ztrue, Zpred){
  
  # test <- aim2_sim_data_setting1(M = 1000,
  #                                nGWAS = 2,
  #                                nAnn = 10,
  #                                percent_ones_ann = 0.10,
  #                                percent_overlap_ann = 0.5,
  #                                percent_overlap_gwas = 1,
  #                                trueAlphaVec = c(0.2, 0.2))
  # fit <- multiGPATree(gwasPval = test$gwasPval,
  #                     annMat = test$annMat,
  #                     initAlpha = 0.1,
  #                     cpTry = 0.001,
  #                     ncore = 1)
  # class(fit)
  # fit@fit$P1_P2
  # Zpred <- multiGPATree::assoc(fit@fit$P1_P2, FDR = 0.05, fdrControl = 'global')
  # Ztrue <- test$Ztrue
  # head(Zpred)
  # dim(Zpred)
  
  # head(Ztrue)
  # dim(Ztrue)
  # colSums(Zpred[, 1:3])
  # table(Ztrue[, 1], Ztrue[, 2])
  # table(Zpred[, 1], Ztrue[, 1])
  # table(Zpred[, 2], Ztrue[, 2])
  # table(Zpred[, 3], Ztrue[, 1] == 1 & Ztrue[, 2] == 1 )
  
  
  # start here
  
  # all phenotype 1 ####
  TP_10 <- length(which(Ztrue[, 1] == 1 & Zpred[, 1] == 1)) # true positives
  TN_10 <- length(which(Ztrue[, 1] == 0 & Zpred[, 1] == 0)) # true negatives
  FN_10 <- length(which(Ztrue[, 1] == 1 & Zpred[, 1] == 0)) # false negatives
  FP_10 <- length(which(Ztrue[, 1] == 0 & Zpred[, 1] == 1)) # false positives
  power_p1 <- sensitivity_p1 <- TP_10/(TP_10 + FN_10 + 1) 
  specificity_p1 <- TN_10/(TN_10 + FP_10 + 1)
  pred_gFDR_p1 <- FP_10/(FP_10 + TP_10 + 1)
  
  # all phenotype 2 ####
  TP_01 <- length(which(Ztrue[, 2] == 1 & Zpred[, 2] == 1)) # true positives
  TN_01 <- length(which(Ztrue[, 2] == 0 & Zpred[, 2] == 0)) # true negatives
  FN_01 <- length(which(Ztrue[, 2] == 1 & Zpred[, 2] == 0)) # false negatives
  FP_01 <- length(which(Ztrue[, 2] == 0 & Zpred[, 2] == 1)) # false positives
  power_p2 <- sensitivity_p2 <- TP_01/(TP_01 + FN_01 + 1)
  specificity_p2 <- TN_01/(TN_01 + FP_01 + 1)
  pred_gFDR_p2 <- FP_01/(FP_01 + TP_01 + 1)
  
  # phenotype 1 and 2 ####
  # Zpred11 <- rep(0, nrow(Zpred))
  # Zpred11_ind <- which(Zpred[, 1] == 1 & Zpred[, 2] == 1)
  # Zpred11[Zpred11_ind] <- 1
  Ztrue11 <- rep(0, nrow(Ztrue))
  Ztrue11_ind <- which(Ztrue[, 1] == 1 & Ztrue[, 2] == 1)
  Ztrue11[Ztrue11_ind] <- 1
  
  Zpred11 <- Zpred[, 3]
  
  TP_11 <- length(which(Ztrue11 == 1 & Zpred11 == 1)) # true positives
  TN_11 <- length(which(Ztrue11 == 0 & Zpred11 == 0)) # true negatives
  FN_11 <- length(which(Ztrue11 == 1 & Zpred11 == 0)) # false negatives
  FP_11 <- length(which(Ztrue11 == 0 & Zpred11 == 1)) # false positives
  power_p1p2 <- sensitivity_p1p2 <- TP_11/(TP_11 + FN_11 + 1)
  specificity_p1p2 <- TN_11/(TN_11 + FP_11 + 1)
  pred_gFDR_p1p2 <- FP_11/(FP_11 + TP_11 + 1)
  
  out <- c(TP_10, TN_10, FN_10, FP_10, power_p1, sensitivity_p1, specificity_p1, pred_gFDR_p1,
           TP_01, TN_01, FN_01, FP_01, power_p2, sensitivity_p2, specificity_p2, pred_gFDR_p2,
           TP_11, TN_11, FN_11, FP_11, power_p1p2, sensitivity_p1p2, specificity_p1p2, pred_gFDR_p1p2)
  names(out) <- c("TP_10", "TN_10", "FN_10", "FP_10", "power_p1", "sensitivity_p1", "specificity_p1", "pred_gFDR_p1",
                  "TP_01", "TN_01", "FN_01", "FP_01", "power_p2", "sensitivity_p2", "specificity_p2", "pred_gFDR_p2",
                  "TP_11", "TN_11", "FN_11", "FP_11", "power_p1p2", "sensitivity_p1p2", "specificity_p1p2", "pred_gFDR_p1p2")
  
  return(out)
  
  
}

# 2. function to get AUC for multiGPATree aim 2 ####
aim2_perf_summary <- function(Ztrue, multiGPATree_fit){ # multiGPATree_fit is an object of class multigpatree (not multigpatreemethod)
  
  # test <- aim2_sim_data_setting1(M = 1000,
  #                                nGWAS = 2,
  #                                nAnn = 10,
  #                                percent_ones_ann = 0.10,
  #                                percent_overlap_ann = 0.5,
  #                                percent_overlap_gwas = 1,
  #                                trueAlphaVec = c(0.2, 0.2))
  # multiGPATree_fit <- multiGPATree::multiGPATree(gwasPval = test$gwasPval,
  #                                                annMat = test$annMat,
  #                                                initAlpha = 0.1,
  #                                                cpTry = 0.001,
  #                                                ncore = 1)
  # Ztrue <- test$Ztrue
  # multiGPATree_fit <- multiGPATree_fit@fit$P1_P2 # multiGPATree_fit is an object of class multiGPATree so only storing element of that class
  
  
  # start here
  global_fdr_control <- round(seq(0, 1, 0.005), 4)
  Zpred <- multiGPATree::assoc(object = multiGPATree_fit, 
                               FDR = global_fdr_control[1],
                               fdrControl = "global")
  outall <- perf_summary_func(Ztrue = Ztrue, Zpred = Zpred)
  
  Zpred_lfdr <- multiGPATree::assoc(object = multiGPATree_fit, 
                                    FDR = global_fdr_control[1], 
                                    fdrControl = 'local')
  outall_lfdr <- perf_summary_func(Ztrue = Ztrue, Zpred = Zpred_lfdr)
  names(outall_lfdr) <- paste(names(outall_lfdr), '_lfdr', sep = '')
  outall  <- c(global_fdr_control[1], outall, outall_lfdr)
  names(outall)[1] <- 'gFDR_control'
  
  for (i in 2:length(global_fdr_control) ) {
    Zpred <- multiGPATree::assoc(object = multiGPATree_fit, 
                                 FDR = global_fdr_control[i], 
                                 fdrControl = 'global')
    out <- perf_summary_func(Ztrue = Ztrue, Zpred = Zpred)
    Zpred_lfdr <- multiGPATree::assoc(object = multiGPATree_fit, 
                                      FDR = global_fdr_control[i], 
                                      fdrControl = 'local')
    out_lfdr <- perf_summary_func(Ztrue = Ztrue, Zpred = Zpred_lfdr)
    
    out  <- c(global_fdr_control[i], out, out_lfdr)
    outall <- rbind(outall, out)
  }
  
  outall <- as.data.frame(outall)
  
  out <- list()
  
  # control at 0.01
  out$multiGPATree_sensitivity_gFDRcontrol0.01_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.01_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.01_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  out$multiGPATree_specificity_gFDRcontrol0.01_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_specificity_gFDRcontrol0.01_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_specificity_gFDRcontrol0.01_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  out$multiGPATree_power_gFDRcontrol0.01_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_power_gFDRcontrol0.01_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_power_gFDRcontrol0.01_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  out$multiGPATree_power_lFDRcontrol0.01_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_power_lFDRcontrol0.01_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_power_lFDRcontrol0.01_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  
  
  out$multiGPATree_predgFDR_gFDRcontrol0.01_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.01_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.01_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  out$multiGPATree_predlFDR_lFDRcontrol0.01_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.01_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.01_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  
  # control at 0.05
  out$multiGPATree_sensitivity_gFDRcontrol0.05_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.05_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.05_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  out$multiGPATree_specificity_gFDRcontrol0.05_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_specificity_gFDRcontrol0.05_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_specificity_gFDRcontrol0.05_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  out$multiGPATree_power_gFDRcontrol0.05_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_power_gFDRcontrol0.05_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_power_gFDRcontrol0.05_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  out$multiGPATree_power_lFDRcontrol0.05_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_power_lFDRcontrol0.05_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_power_lFDRcontrol0.05_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  
  out$multiGPATree_predgFDR_gFDRcontrol0.05_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.05_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.05_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  out$multiGPATree_predlFDR_lFDRcontrol0.05_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.05_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.05_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  
  
  # control at 0.10
  out$multiGPATree_sensitivity_gFDRcontrol0.10_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.10_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.10_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.10) ]
  
  out$multiGPATree_specificity_gFDRcontrol0.10_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_specificity_gFDRcontrol0.10_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_specificity_gFDRcontrol0.10_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.10) ]
  
  out$multiGPATree_power_gFDRcontrol0.10_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_power_gFDRcontrol0.10_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_power_gFDRcontrol0.10_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.10) ]
  
  out$multiGPATree_power_lFDRcontrol0.10_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_power_lFDRcontrol0.10_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_power_lFDRcontrol0.10_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.10) ]
  
  out$multiGPATree_predgFDR_gFDRcontrol0.10_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.10_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.10) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.10_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.10) ]
  
  out$multiGPATree_predlFDR_lFDRcontrol0.10_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.1) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.10_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.10_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  
  # control at 0.15
  out$multiGPATree_sensitivity_gFDRcontrol0.15_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.15_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.15_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  out$multiGPATree_specificity_gFDRcontrol0.15_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_specificity_gFDRcontrol0.15_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_specificity_gFDRcontrol0.15_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  out$multiGPATree_power_gFDRcontrol0.15_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_power_gFDRcontrol0.15_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_power_gFDRcontrol0.15_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  out$multiGPATree_power_lFDRcontrol0.15_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_power_lFDRcontrol0.15_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_power_lFDRcontrol0.15_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  
  out$multiGPATree_predgFDR_gFDRcontrol0.15_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.15_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.15_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  out$multiGPATree_predlFDR_lFDRcontrol0.15_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.15_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.15_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  
  # control at 0.20
  out$multiGPATree_sensitivity_gFDRcontrol0.20_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.20_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_sensitivity_gFDRcontrol0.20_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.20) ]
  
  out$multiGPATree_specificity_gFDRcontrol0.20_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_specificity_gFDRcontrol0.20_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_specificity_gFDRcontrol0.20_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.20) ]
  
  out$multiGPATree_power_gFDRcontrol0.20_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_power_gFDRcontrol0.20_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_power_gFDRcontrol0.20_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.20) ]
  
  out$multiGPATree_power_lFDRcontrol0.20_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_power_lFDRcontrol0.20_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_power_lFDRcontrol0.20_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.20) ]
  
  out$multiGPATree_predgFDR_gFDRcontrol0.20_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.20_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_predgFDR_gFDRcontrol0.20_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.20) ]
  
  out$multiGPATree_predlFDR_lFDRcontrol0.20_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.20_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.20) ]
  out$multiGPATree_predlFDR_lFDRcontrol0.20_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.20) ]
  
  out$multiGPATree_AUC_P1 <- round(pracma::trapz(1-outall$specificity_p1, outall$sensitivity_p1), 4)
  out$multiGPATree_AUC_P2 <- round(pracma::trapz(1-outall$specificity_p2, outall$sensitivity_p2), 4)
  out$multiGPATree_AUC_P1_P2 <- round(pracma::trapz(1-outall$specificity_p1p2, outall$sensitivity_p1p2), 4)
  
  
  out$multiGPATree_AUC_P1_lfdr <- round(pracma::trapz(1-outall$specificity_p1_lfdr, outall$sensitivity_p1_lfdr), 4)
  out$multiGPATree_AUC_P2_lfdr <- round(pracma::trapz(1-outall$specificity_p2_lfdr, outall$sensitivity_p2_lfdr), 4)
  out$multiGPATree_AUC_P1_P2_lfdr <- round(pracma::trapz(1-outall$specificity_p1p2_lfdr, outall$sensitivity_p1p2_lfdr), 4)
  
  return(out)
  
}

# 3. function for LPM results ####
aim2_LPM_result <- function(gwasPval, annMat, Ztrue){
  
  # # Creating the dataset
  # simulated_data <- setting1_aim2_sim_data(M = 1000,
  #                                         nGWAS = 2,
  #                                         nAnn = 10,
  #                                         percent_ones_ann = 0.1,
  #                                         percent_overlap_ann = 0.50,
  #                                         percent_overlap_gwas = 1,
  #                                         trueAlphaVec = c(0.2, 0.2))
  # gwasPval <- simulated_data$gwasPval
  # annMat <- simulated_data$annMat
  # head(gwasPval)
  # head(annMat)
  # Ztrue <- simulated_data$Ztrue
  
  # start here
  library(pbivnorm)
  library(LPM)
  LPMData <- list('data' = list(),
                  'X' = data.frame(cbind(SNP = 1:nrow(gwasPval), annMat))) # LPM data setup
  for (i in 1:ncol(gwasPval)) {
    LPMData$data[[i]] <- as.data.frame(cbind(SNP = 1:nrow(gwasPval), p = gwasPval[, i]))
    names(LPMData$data)[[i]] <- colnames(gwasPval)[i]
  }
  
  LPM_start_time <- proc.time()
  bLPMfit <- LPM::bLPM(data = LPMData$data, X = LPMData$X)
  LPMfit <- LPM::LPM(bLPMfit)
  LPM_end_time <- proc.time()
  LPM_time <- as.numeric(LPM_end_time[[3]] - LPM_start_time[[3]])
  print('okay LPM-1')
  
  # LPM power when gFDR controlled at 0.05,  predicted gFDR, AUC
  
  post12 <- LPM::post(LPMData$data[c(1,2)], X = LPMData$X, id = c(1,2), LPMfit)
  
  global_fdr_control <- round(seq(0, 1, 0.005), 4)
  outall <- as.data.frame(t(rep(NA, 25)))
  outall <- outall[-1, ]
  for (i in 1:length(global_fdr_control) ) {
    # i = 1
    assoc12 <- LPM::assoc(post12, FDRset = global_fdr_control[i], fdrControl = "global")
    
    TP_10 <- length(which(Ztrue[, 1] == 1 & assoc12$eta.marginal1 == 1)) # true positives
    TN_10 <- length(which(Ztrue[, 1] == 0 & assoc12$eta.marginal1 == 0)) # true negatives
    FN_10 <- length(which(Ztrue[, 1] == 1 & assoc12$eta.marginal1 == 0)) # false negatives
    FP_10 <- length(which(Ztrue[, 1] == 0 & assoc12$eta.marginal1 == 1)) # false positives
    power_p1 <- sensitivity_p1 <- TP_10/(TP_10 + FN_10 + 1) 
    specificity_p1 <- TN_10/(TN_10 + FP_10 + 1)
    pred_gFDR_p1 <- FP_10/(FP_10 + TP_10 + 1)
    
    
    TP_01 <- length(which(Ztrue[, 2] == 1 & assoc12$eta.marginal2 == 1)) # true positives
    TN_01 <- length(which(Ztrue[, 2] == 0 & assoc12$eta.marginal2 == 0)) # true negatives
    FN_01 <- length(which(Ztrue[, 2] == 1 & assoc12$eta.marginal2 == 0)) # false negatives
    FP_01 <- length(which(Ztrue[, 2] == 0 & assoc12$eta.marginal2 == 1)) # false positives
    power_p2 <- sensitivity_p2 <- TP_01/(TP_01 + FN_01 + 1)
    specificity_p2 <- TN_01/(TN_01 + FP_01 + 1)
    pred_gFDR_p2 <- FP_01/(FP_01 + TP_01 + 1)
    
    Ztrue11 <- rep(0, nrow(Ztrue))
    Ztrue11[Ztrue[, 1] == 1 & Ztrue[, 2] == 1] <- 1
    Zpred11 <- assoc12$eta.joint
    
    TP_11 <- length(which(Ztrue11 == 1 & Zpred11 == 1)) # true positives
    TN_11 <- length(which(Ztrue11 == 0 & Zpred11 == 0)) # true negatives
    FN_11 <- length(which(Ztrue11 == 1 & Zpred11 == 0)) # false negatives
    FP_11 <- length(which(Ztrue11 == 0 & Zpred11 == 1)) # false positives
    power_p1p2 <- sensitivity_p1p2 <- TP_11/(TP_11 + FN_11 + 1)
    specificity_p1p2 <- TN_11/(TN_11 + FP_11 + 1)
    pred_gFDR_p1p2 <- FP_11/(FP_11 + TP_11 + 1)
    
    # local fdr control
    assoc12_lfdr <- LPM::assoc(post12, FDRset = global_fdr_control[i], fdrControl = "local")
    
    TP_10_lfdr <- length(which(Ztrue[, 1] == 1 & assoc12_lfdr$eta.marginal1 == 1)) # true positives
    TN_10_lfdr <- length(which(Ztrue[, 1] == 0 & assoc12_lfdr$eta.marginal1 == 0)) # true negatives
    FN_10_lfdr <- length(which(Ztrue[, 1] == 1 & assoc12_lfdr$eta.marginal1 == 0)) # false negatives
    FP_10_lfdr <- length(which(Ztrue[, 1] == 0 & assoc12_lfdr$eta.marginal1 == 1)) # false positives
    power_p1_lfdr <- sensitivity_p1_lfdr <- TP_10_lfdr/(TP_10_lfdr + FN_10_lfdr + 1)
    specificity_p1_lfdr <- TN_10_lfdr/(TN_10_lfdr + FP_10_lfdr + 1)
    pred_gFDR_p1_lfdr <- FP_10_lfdr/(FP_10_lfdr + TP_10_lfdr + 1)
    
    TP_01_lfdr <- length(which(Ztrue[, 2] == 1 & assoc12_lfdr$eta.marginal2 == 1)) # true positives
    TN_01_lfdr <- length(which(Ztrue[, 2] == 0 & assoc12_lfdr$eta.marginal2 == 0)) # true negatives
    FN_01_lfdr <- length(which(Ztrue[, 2] == 1 & assoc12_lfdr$eta.marginal2 == 0)) # false negatives
    FP_01_lfdr <- length(which(Ztrue[, 2] == 0 & assoc12_lfdr$eta.marginal2 == 1)) # false positives
    power_p2_lfdr <- sensitivity_p2_lfdr <- TP_01_lfdr/(TP_01_lfdr + FN_01_lfdr + 1)
    specificity_p2_lfdr <- TN_01_lfdr/(TN_01_lfdr + FP_01_lfdr + 1)
    pred_gFDR_p2_lfdr <- FP_01_lfdr/(FP_01_lfdr + TP_01_lfdr + 1)
    
    # Ztrue11 <- rep(0, nrow(Ztrue))
    # Ztrue11 <- ifelse(Ztrue[, 1] == 1 & Ztrue[, 2] == 1, 1, Ztrue11)
    Zpred11 <- assoc12_lfdr$eta.joint
    TP_11_lfdr <- length(which(Ztrue11 == 1 & Zpred11 == 1)) # true positives
    TN_11_lfdr <- length(which(Ztrue11 == 0 & Zpred11 == 0)) # true negatives
    FN_11_lfdr <- length(which(Ztrue11 == 1 & Zpred11 == 0)) # false negatives
    FP_11_lfdr <- length(which(Ztrue11 == 0 & Zpred11 == 1)) # false positives
    power_p1p2_lfdr <- sensitivity_p1p2_lfdr <- TP_11_lfdr/(TP_11_lfdr + FN_11_lfdr + 1)
    specificity_p1p2_lfdr <- TN_11_lfdr/(TN_11_lfdr + FP_11_lfdr + 1)
    pred_gFDR_p1p2_lfdr <- FP_11_lfdr/(FP_11_lfdr + TP_11_lfdr + 1)
    
    out  <- c(global_fdr_control[i], sensitivity_p1, specificity_p1, sensitivity_p2, specificity_p2,
              sensitivity_p1p2, specificity_p1p2, power_p1, power_p2, power_p1p2, pred_gFDR_p1,
              pred_gFDR_p2, pred_gFDR_p1p2, sensitivity_p1_lfdr, specificity_p1_lfdr, 
              sensitivity_p2_lfdr, specificity_p2_lfdr, sensitivity_p1p2_lfdr, specificity_p1p2_lfdr, 
              power_p1_lfdr, power_p2_lfdr, power_p1p2_lfdr, pred_gFDR_p1_lfdr,
              pred_gFDR_p2_lfdr, pred_gFDR_p1p2_lfdr)
    
    outall <- rbind(outall, out)
  }
  
  outall <- as.data.frame(outall)
  names(outall) <- c('gFDR_control', 'sensitivity_p1', "specificity_p1", "sensitivity_p2", "specificity_p2",
                     "sensitivity_p1p2", "specificity_p1p2", "power_p1", "power_p2", "power_p1p2",
                     "pred_gFDR_p1", "pred_gFDR_p2", "pred_gFDR_p1p2", 'sensitivity_p1_lfdr', 
                     "specificity_p1_lfdr", "sensitivity_p2_lfdr", "specificity_p2_lfdr",
                     "sensitivity_p1p2_lfdr", "specificity_p1p2_lfdr", "power_p1_lfdr", "power_p2_lfdr", "power_p1p2_lfdr",
                     "pred_gFDR_p1_lfdr", "pred_gFDR_p2_lfdr", "pred_gFDR_p1p2_lfdr")
  
  LPM_AUC_10 <- round(pracma::trapz(1-outall$specificity_p1, outall$sensitivity_p1), 4)
  LPM_AUC_01 <- round(pracma::trapz(1-outall$specificity_p2, outall$sensitivity_p2), 4)
  LPM_AUC_11 <- round(pracma::trapz(1-outall$specificity_p1p2, outall$sensitivity_p1p2), 4)
  
  
  # output  
  
  LPM_out <- list()
  LPM_out$LPM_alpha_P1 <- LPMfit$alpha[[1]]
  LPM_out$LPM_alpha_P2 <- LPMfit$alpha[[2]]
  # LPM_out$LPM_beta <- LPMfit$beta
  # LPM_out$LPM_corr <- LPMfit$R
  LPM_out$LPM_time <- LPM_time
  
  # gFDR control = 0.01
  LPM_out$LPM_sensitivity_gFDRcontrol0.01_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.01_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.01_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  LPM_out$LPM_specificity_gFDRcontrol0.01_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_specificity_gFDRcontrol0.01_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_specificity_gFDRcontrol0.01_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  LPM_out$LPM_power_gFDRcontrol0.01_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_power_gFDRcontrol0.01_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_power_gFDRcontrol0.01_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  LPM_out$LPM_power_lFDRcontrol0.01_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_power_lFDRcontrol0.01_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_power_lFDRcontrol0.01_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  
  LPM_out$LPM_predgFDR_gFDRcontrol0.01_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.01_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.01_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.01) ]
  
  LPM_out$LPM_predlFDR_lFDRcontrol0.01_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.01_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.01_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  
  # gFDR control = 0.05
  LPM_out$LPM_sensitivity_gFDRcontrol0.05_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.05_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.05_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  LPM_out$LPM_specificity_gFDRcontrol0.05_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_specificity_gFDRcontrol0.05_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_specificity_gFDRcontrol0.05_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  LPM_out$LPM_power_gFDRcontrol0.05_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_power_gFDRcontrol0.05_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_power_gFDRcontrol0.05_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  LPM_out$LPM_power_lFDRcontrol0.05_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_power_lFDRcontrol0.05_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_power_lFDRcontrol0.05_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  
  LPM_out$LPM_predgFDR_gFDRcontrol0.05_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.05_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.05_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.05) ]
  
  LPM_out$LPM_predlFDR_lFDRcontrol0.05_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.05_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.05_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  
  # gFDR control = 0.10
  LPM_out$LPM_sensitivity_gFDRcontrol0.1_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.1_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.1_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.1) ]
  
  LPM_out$LPM_specificity_gFDRcontrol0.1_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_specificity_gFDRcontrol0.1_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_specificity_gFDRcontrol0.1_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.1) ]
  
  LPM_out$LPM_power_gFDRcontrol0.1_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_power_gFDRcontrol0.1_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_power_gFDRcontrol0.1_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.1) ]
  
  LPM_out$LPM_power_lFDRcontrol0.1_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_power_lFDRcontrol0.1_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_power_lFDRcontrol0.1_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  
  LPM_out$LPM_predgFDR_gFDRcontrol0.1_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.1_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.1_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.1) ]
  
  LPM_out$LPM_predlFDR_lFDRcontrol0.1_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.1_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.1_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  
  
  # gFDR control = 0.15
  LPM_out$LPM_sensitivity_gFDRcontrol0.15_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.15_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.15_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  LPM_out$LPM_specificity_gFDRcontrol0.15_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_specificity_gFDRcontrol0.15_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_specificity_gFDRcontrol0.15_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  LPM_out$LPM_power_gFDRcontrol0.15_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_power_gFDRcontrol0.15_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_power_gFDRcontrol0.15_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  LPM_out$LPM_power_lFDRcontrol0.15_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_power_lFDRcontrol0.15_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_power_lFDRcontrol0.15_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  
  LPM_out$LPM_predgFDR_gFDRcontrol0.15_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.15_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.15_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.15) ]
  
  LPM_out$LPM_predlFDR_lFDRcontrol0.15_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.15_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.15_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  
  # gFDR control = 0.20
  LPM_out$LPM_sensitivity_gFDRcontrol0.2_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.2_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_sensitivity_gFDRcontrol0.2_P1_P2 <- outall$sensitivity_p1p2[ which(outall$gFDR_control == 0.2) ]
  
  LPM_out$LPM_specificity_gFDRcontrol0.2_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_specificity_gFDRcontrol0.2_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_specificity_gFDRcontrol0.2_P1_P2 <- outall$specificity_p1p2[ which(outall$gFDR_control == 0.2) ]
  
  LPM_out$LPM_power_gFDRcontrol0.2_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_power_gFDRcontrol0.2_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_power_gFDRcontrol0.2_P1_P2 <- outall$power_p1p2[ which(outall$gFDR_control == 0.2) ]
  
  LPM_out$LPM_power_lFDRcontrol0.2_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_power_lFDRcontrol0.2_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_power_lFDRcontrol0.2_P1_P2 <- outall$power_p1p2_lfdr[ which(outall$gFDR_control == 0.2) ]
  
  LPM_out$LPM_predgFDR_gFDRcontrol0.2_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.2_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_predgFDR_gFDRcontrol0.2_P1_P2 <- outall$pred_gFDR_p1p2[ which(outall$gFDR_control == 0.2) ]
  
  LPM_out$LPM_predlFDR_lFDRcontrol0.2_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.2_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.2) ]
  LPM_out$LPM_predlFDR_lFDRcontrol0.2_P1_P2 <- outall$pred_gFDR_p1p2_lfdr[ which(outall$gFDR_control == 0.2) ]
  
  LPM_out$LPM_AUC_P1 <- LPM_AUC_10
  LPM_out$LPM_AUC_P2 <- LPM_AUC_01
  LPM_out$LPM_AUC_P1_P2 <- LPM_AUC_11
  
  LPM_out$LPM_AUC_P1_lfdr <- round(pracma::trapz(1-outall$specificity_p1_lfdr, outall$sensitivity_p1_lfdr), 4)
  LPM_out$LPM_AUC_P2_lfdr <- round(pracma::trapz(1-outall$specificity_p2_lfdr, outall$sensitivity_p2_lfdr), 4)
  LPM_out$LPM_AUC_P1_P2_lfdr <- round(pracma::trapz(1-outall$specificity_p1p2_lfdr, outall$sensitivity_p1p2_lfdr), 4)
  
  return(LPM_out)
  
}

# 3. function for GPATree results ####
aim2_gpatree_result <- function(gwasPval, annMat, Ztrue, initAlpha, cpTry){
  
  # # Creating the dataset
  # simulated_data <- setting1_aim2_sim_data(M = 1000,
  #                                         nGWAS = 2,
  #                                         nAnn = 10,
  #                                         percent_ones_ann = 0.1,
  #                                         percent_overlap_ann = 0.50,
  #                                         percent_overlap_gwas = 1,
  #                                         trueAlphaVec = c(0.2, 0.2))
  # gwasPval <- simulated_data$gwasPval
  # annMat <- simulated_data$annMat
  # Ztrue <- simulated_data$Ztrue
  # initAlpha <- 0.1
  # cpTry <- 0.001
  # head(gwasPval)
  # head(annMat)
  
  
  # start here
  library(GPATree)

  
  gpatree_start_time <- proc.time()
  gpatree_fit_p1 <- GPATree::GPATree(gwasPval = gwasPval[, 1], annMat = annMat, initAlpha = initAlpha, cpTry = cpTry)
  gpatree_fit_p2 <- GPATree::GPATree(gwasPval = gwasPval[, 2], annMat = annMat, initAlpha = initAlpha, cpTry = cpTry)
  gpatree_end_time <- proc.time()
  gpatree_time <- as.numeric(gpatree_end_time[[3]] - gpatree_start_time[[3]])
  print('okay LPM-1')
  
  # GPATree power when gFDR controlled at 0.05,  predicted gFDR, AUC
  
  global_fdr_control <- round(seq(0, 1, 0.005), 4)
  outall <- as.data.frame(t(rep(NA, 17)))
  dim(outall)
  outall <- outall[-1, ]
  for (i in 1:length(global_fdr_control) ) {
    # i = 1
    assoc1 <- GPATree::assoc(gpatree_fit_p1, FDR = global_fdr_control[i], fdrControl = "global")
    head(assoc1)
    TP_10 <- length(which(Ztrue[, 1] == 1 & assoc1$P1 == 1)) # true positives
    TN_10 <- length(which(Ztrue[, 1] == 0 & assoc1$P1 == 0)) # true negatives
    FN_10 <- length(which(Ztrue[, 1] == 1 & assoc1$P1 == 0)) # false negatives
    FP_10 <- length(which(Ztrue[, 1] == 0 & assoc1$P1 == 1)) # false positives
    power_p1 <- sensitivity_p1 <- TP_10/(TP_10 + FN_10 + 1) 
    specificity_p1 <- TN_10/(TN_10 + FP_10 + 1)
    pred_gFDR_p1 <- FP_10/(FP_10 + TP_10 + 1)
    
    
    assoc2 <- GPATree::assoc(gpatree_fit_p2, FDR = global_fdr_control[i], fdrControl = "global")
    head(assoc2)
    TP_01 <- length(which(Ztrue[, 2] == 1 & assoc2$P1 == 1)) # true positives
    TN_01 <- length(which(Ztrue[, 2] == 0 & assoc2$P1 == 0)) # true negatives
    FN_01 <- length(which(Ztrue[, 2] == 1 & assoc2$P1 == 0)) # false negatives
    FP_01 <- length(which(Ztrue[, 2] == 0 & assoc2$P1 == 1)) # false positives
    power_p2 <- sensitivity_p2 <- TP_01/(TP_01 + FN_01 + 1)
    specificity_p2 <- TN_01/(TN_01 + FP_01 + 1)
    pred_gFDR_p2 <- FP_01/(FP_01 + TP_01 + 1)
    
    # local fdr control
    assoc1_lfdr <- GPATree::assoc(gpatree_fit_p1, FDR = global_fdr_control[i], fdrControl = "local")
    head(assoc1_lfdr)
    TP_10_lfdr <- length(which(Ztrue[, 1] == 1 & assoc1_lfdr$P1 == 1)) # true positives
    TN_10_lfdr <- length(which(Ztrue[, 1] == 0 & assoc1_lfdr$P1 == 0)) # true negatives
    FN_10_lfdr <- length(which(Ztrue[, 1] == 1 & assoc1_lfdr$P1 == 0)) # false negatives
    FP_10_lfdr <- length(which(Ztrue[, 1] == 0 & assoc1_lfdr$P1 == 1)) # false positives
    power_p1_lfdr <- sensitivity_p1_lfdr <- TP_10_lfdr/(TP_10_lfdr + FN_10_lfdr + 1)
    specificity_p1_lfdr <- TN_10_lfdr/(TN_10_lfdr + FP_10_lfdr + 1)
    pred_gFDR_p1_lfdr <- FP_10_lfdr/(FP_10_lfdr + TP_10_lfdr + 1)
    
    assoc2_lfdr <- GPATree::assoc(gpatree_fit_p2, FDR = global_fdr_control[i], fdrControl = "local")
    # head(assoc2_lfdr)
    TP_01_lfdr <- length(which(Ztrue[, 2] == 1 & assoc2_lfdr$P1 == 1)) # true positives
    TN_01_lfdr <- length(which(Ztrue[, 2] == 0 & assoc2_lfdr$P1 == 0)) # true negatives
    FN_01_lfdr <- length(which(Ztrue[, 2] == 1 & assoc2_lfdr$P1 == 0)) # false negatives
    FP_01_lfdr <- length(which(Ztrue[, 2] == 0 & assoc2_lfdr$P1 == 1)) # false positives
    power_p2_lfdr <- sensitivity_p2_lfdr <- TP_01_lfdr/(TP_01_lfdr + FN_01_lfdr + 1)
    specificity_p2_lfdr <- TN_01_lfdr/(TN_01_lfdr + FP_01_lfdr + 1)
    pred_gFDR_p2_lfdr <- FP_01_lfdr/(FP_01_lfdr + TP_01_lfdr + 1)
    
    out  <- c(global_fdr_control[i], 
              sensitivity_p1, specificity_p1, sensitivity_p2, specificity_p2,
              power_p1, power_p2, 
              pred_gFDR_p1, pred_gFDR_p2, 
              sensitivity_p1_lfdr, specificity_p1_lfdr, sensitivity_p2_lfdr, specificity_p2_lfdr,
              power_p1_lfdr, power_p2_lfdr, 
              pred_gFDR_p1_lfdr, pred_gFDR_p2_lfdr)
    # length(out)
    
    outall <- rbind(outall, out)
  }
  
  outall <- as.data.frame(outall)
  names(outall) <- c('gFDR_control', 'sensitivity_p1', "specificity_p1", "sensitivity_p2", "specificity_p2",
                     "power_p1", "power_p2", 
                     "pred_gFDR_p1", "pred_gFDR_p2", 
                     'sensitivity_p1_lfdr', "specificity_p1_lfdr", "sensitivity_p2_lfdr", "specificity_p2_lfdr",
                     "power_p1_lfdr", "power_p2_lfdr",
                     "pred_gFDR_p1_lfdr", "pred_gFDR_p2_lfdr")
  
  gpatree_AUC_10 <- round(pracma::trapz(1-outall$specificity_p1, outall$sensitivity_p1), 4)
  gpatree_AUC_01 <- round(pracma::trapz(1-outall$specificity_p2, outall$sensitivity_p2), 4)
  
  
  # output  
  # gpatree_fit_p1@fit$alpha
  gpatree_out <- list()
  gpatree_out$gpatree_alpha_P1 <- gpatree_fit_p1@fit$alpha
  gpatree_out$gpatree_alpha_P2 <- gpatree_fit_p2@fit$alpha
  # LPM_out$LPM_beta <- LPMfit$beta
  # LPM_out$LPM_corr <- LPMfit$R
  gpatree_out$gpatree_time <- gpatree_time
  
  # gFDR control = 0.01
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.01_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.01) ]
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.01_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.01) ]
  
  gpatree_out$gpatree_specificity_gFDRcontrol0.01_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.01) ]
  gpatree_out$gpatree_specificity_gFDRcontrol0.01_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.01) ]
  
  gpatree_out$gpatree_power_gFDRcontrol0.01_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.01) ]
  gpatree_out$gpatree_power_gFDRcontrol0.01_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.01) ]
  
  gpatree_out$gpatree_power_lFDRcontrol0.01_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.01) ]
  gpatree_out$gpatree_power_lFDRcontrol0.01_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.01_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.01) ]
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.01_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.01) ]
  
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.01_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.01) ]
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.01_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.01) ]
  
  # gFDR control = 0.05
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.05_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.05) ]
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.05_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.05) ]
  
  gpatree_out$gpatree_specificity_gFDRcontrol0.05_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.05) ]
  gpatree_out$gpatree_specificity_gFDRcontrol0.05_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.05) ]
  
  gpatree_out$gpatree_power_gFDRcontrol0.05_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.05) ]
  gpatree_out$gpatree_power_gFDRcontrol0.05_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.05) ]
  
  gpatree_out$gpatree_power_lFDRcontrol0.05_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.05) ]
  gpatree_out$gpatree_power_lFDRcontrol0.05_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.05_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.05) ]
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.05_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.05) ]
  
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.05_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.05) ]
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.05_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.05) ]
  
  # gFDR control = 0.10
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.1_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.1) ]
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.1_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.1) ]
  
  gpatree_out$gpatree_specificity_gFDRcontrol0.1_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.1) ]
  gpatree_out$gpatree_specificity_gFDRcontrol0.1_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.1) ]
  
  gpatree_out$gpatree_power_gFDRcontrol0.1_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.1) ]
  gpatree_out$gpatree_power_gFDRcontrol0.1_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.1) ]
  
  gpatree_out$gpatree_power_lFDRcontrol0.1_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.1) ]
  gpatree_out$gpatree_power_lFDRcontrol0.1_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.1_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.1) ]
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.1_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.1) ]
  
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.1_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.1) ]
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.1_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.1) ]
  
  
  # gFDR control = 0.15
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.15_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.15) ]
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.15_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.15) ]
  
  gpatree_out$gpatree_specificity_gFDRcontrol0.15_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.15) ]
  gpatree_out$gpatree_specificity_gFDRcontrol0.15_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.15) ]
  
  gpatree_out$gpatree_power_gFDRcontrol0.15_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.15) ]
  gpatree_out$gpatree_power_gFDRcontrol0.15_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.15) ]
  
  gpatree_out$gpatree_power_lFDRcontrol0.15_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.15) ]
  gpatree_out$gpatree_power_lFDRcontrol0.15_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.15_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.15) ]
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.15_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.15) ]
  
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.15_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.15) ]
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.15_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.15) ]
  
  # gFDR control = 0.20
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.2_P1 <- outall$sensitivity_p1[ which(outall$gFDR_control == 0.2) ]
  gpatree_out$gpatree_sensitivity_gFDRcontrol0.2_P2 <- outall$sensitivity_p2[ which(outall$gFDR_control == 0.2) ]
  
  gpatree_out$gpatree_specificity_gFDRcontrol0.2_P1 <- outall$specificity_p1[ which(outall$gFDR_control == 0.2) ]
  gpatree_out$gpatree_specificity_gFDRcontrol0.2_P2 <- outall$specificity_p2[ which(outall$gFDR_control == 0.2) ]
  
  gpatree_out$gpatree_power_gFDRcontrol0.2_P1 <- outall$power_p1[ which(outall$gFDR_control == 0.2) ]
  gpatree_out$gpatree_power_gFDRcontrol0.2_P2 <- outall$power_p2[ which(outall$gFDR_control == 0.2) ]
  
  gpatree_out$gpatree_power_lFDRcontrol0.2_P1 <- outall$power_p1_lfdr[ which(outall$gFDR_control == 0.2) ]
  gpatree_out$gpatree_power_lFDRcontrol0.2_P2 <- outall$power_p2_lfdr[ which(outall$gFDR_control == 0.2) ]
  
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.2_P1 <- outall$pred_gFDR_p1[ which(outall$gFDR_control == 0.2) ]
  gpatree_out$gpatree_predgFDR_gFDRcontrol0.2_P2 <- outall$pred_gFDR_p2[ which(outall$gFDR_control == 0.2) ]
  
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.2_P1 <- outall$pred_gFDR_p1_lfdr[ which(outall$gFDR_control == 0.2) ]
  gpatree_out$gpatree_predlFDR_lFDRcontrol0.2_P2 <- outall$pred_gFDR_p2_lfdr[ which(outall$gFDR_control == 0.2) ]
  
  gpatree_out$gpatree_AUC_P1 <- gpatree_AUC_10
  gpatree_out$gpatree_AUC_P2 <- gpatree_AUC_01
  
  gpatree_out$gpatree_AUC_P1_lfdr <- round(pracma::trapz(1-outall$specificity_p1_lfdr, outall$sensitivity_p1_lfdr), 4)
  gpatree_out$gpatree_AUC_P2_lfdr <- round(pracma::trapz(1-outall$specificity_p2_lfdr, outall$sensitivity_p2_lfdr), 4)
  
  return(gpatree_out)
  
}

# 5. combining multiGPATree, LPM and GPATree models for simulation study of aim 2 ####
aim2_multiGPATree_sim_study <- function(simulated_data, percent_ones_ann, percent_overlap_ann, percent_overlap_gwas,
                                        trueAlphaVec, cpTry, initAlpha, ncore) {
  
  # percent_ones_ann = 0.10
  # percent_overlap_ann = 0.50
  # percent_overlap_gwas = 1
  # trueAlphaVec = c(0.4, 0.4)
  # initAlpha <- 0.1
  # cpTry <- 0.001
  # ncore <- 1
  # simulated_data <- setting1_aim2_sim_data(M = 100,
  #                                         nGWAS = 2,
  #                                         nAnn = 25,
  #                                         percent_ones_ann = percent_ones_ann,
  #                                         percent_overlap_ann = percent_overlap_ann,
  #                                         percent_overlap_gwas = percent_overlap_gwas,
  #                                         trueAlphaVec = trueAlphaVec)
  # testmat <- as.data.frame(cbind(simulated_data$annMat,simulated_data$zMat))
  # heatmap(as.matrix(testmat), Colv = NA, Rowv = NA, scale = 'column')
  

  # start here
  # multiGPATree implementation ####
  
  print('start: multiGPATree Method')
  multiGPATree_time_start <- proc.time()
  multiGPATree_out <- multiGPATree:::multiGPATree(gwasPval = simulated_data$gwasPval,
                                                  annMat = simulated_data$annMat, 
                                                  initAlpha = initAlpha, 
                                                  cpTry = cpTry,
                                                  ncore = ncore)
  multiGPATree_time_end <- proc.time()
  multiGPATree_time <- as.numeric(multiGPATree_time_end[[3]] - multiGPATree_time_start[[3]])
  print('multiGPATree done in main')
  
  multigpatree_class_obj <- multiGPATree_out@fit$P1_P2
  multiGPATree_perf_data <- aim2_perf_summary(Ztrue = simulated_data$Ztrue, multiGPATree_fit = multigpatree_class_obj)
  print('multiGPATree performance plots done')
  
  
  # multiGPATree_selcombCART <- multiGPATree::leaf(object = multigpatree_class_obj)
  # cartmod_frame_ind <- sort(unique(multigpatree_class_obj@fit$fit$where))
  # snp_which_leaf <- rep(NA, nrow(simulated_data$gwasPval))
  # 
  # for (i in 1:length(cartmod_frame_ind)) {
  #   
  #   snp_which_leaf[which(multigpatree_class_obj@fit$fit$where == cartmod_frame_ind[i])] <- paste('LEAF', i)
  #   
  # }
  # 
  # Z_ann <- cbind(as.data.frame(simulated_data$zMat), 'Leaf' = snp_which_leaf)
  # leaf_ann <- multiGPATree::leaf(multigpatree_class_obj)
  # leaf_ann <- cbind('Leaf' = rownames(leaf_ann), as.data.frame(leaf_ann[, 3:ncol(leaf_ann)]))
  # leaf_ann$var_is_1 <- rep(NA, nrow(leaf_ann))
  # 
  # for (i in 1:nrow(leaf_ann)) {
  #   # i = 2
  #   var_list <- colnames(leaf_ann)[which(leaf_ann[i, ] == 1)]
  #   
  #   if (length(var_list) > 0){ 
  #     var_list <- paste("A", sort(readr::parse_number(var_list)), sep = '')
  #     leaf_ann$var_is_1[i] <- paste(var_list, collapse = ", ")    
  #   } else {
  #     leaf_ann$var_is_1[i] <- 'none'
  #   }
  # }
  # 
  # Z_ann_all <- merge(Z_ann, leaf_ann, by = 'Leaf')
  # head(Z_ann_all)
  
  # multiGPATree_perf_data$ann_P1_TP <- length(which(Z_ann_all$var_is_1 == 'A1, A2' & Z_ann_all$`10` == 1))
  # multiGPATree_perf_data$ann_P1_TN <- length(which(Z_ann_all$var_is_1 != 'A1, A2' & Z_ann_all$`10` == 0))
  # multiGPATree_perf_data$ann_P1_FP <- length(which(Z_ann_all$var_is_1 == 'A1, A2' & Z_ann_all$`10` == 0))
  # multiGPATree_perf_data$ann_P1_FN <- length(which(Z_ann_all$var_is_1 != 'A1, A2' & Z_ann_all$`10` == 1))
  # multiGPATree_perf_data$ann_P1_sensitivity <- multiGPATree_perf_data$ann_P1_TP/(multiGPATree_perf_data$ann_P1_TP + multiGPATree_perf_data$ann_P1_FN)
  # multiGPATree_perf_data$ann_P1_specificity <- multiGPATree_perf_data$ann_P1_TN/(multiGPATree_perf_data$ann_P1_TN + multiGPATree_perf_data$ann_P1_FP)
  # 
  # 
  # multiGPATree_perf_data$ann_P2_TP <- length(which(Z_ann_all$var_is_1 == 'A3, A4' & Z_ann_all$`01` == 1))
  # multiGPATree_perf_data$ann_P2_TN <- length(which(Z_ann_all$var_is_1 != 'A3, A4' & Z_ann_all$`01` == 0))
  # multiGPATree_perf_data$ann_P2_FP <- length(which(Z_ann_all$var_is_1 == 'A3, A4' & Z_ann_all$`01` == 0))
  # multiGPATree_perf_data$ann_P2_FN <- length(which(Z_ann_all$var_is_1 != 'A3, A4' & Z_ann_all$`01` == 1))
  # multiGPATree_perf_data$ann_P2_sensitivity <- multiGPATree_perf_data$ann_P2_TP/(multiGPATree_perf_data$ann_P2_TP + multiGPATree_perf_data$ann_P2_FN)
  # multiGPATree_perf_data$ann_P2_specificity <- multiGPATree_perf_data$ann_P2_TN/(multiGPATree_perf_data$ann_P2_TN + multiGPATree_perf_data$ann_P2_FP)
  # 
  # multiGPATree_perf_data$ann_P1_P2_TP <- length(which(Z_ann_all$var_is_1 == 'A5, A6' & Z_ann_all$`11` == 1)) 
  # multiGPATree_perf_data$ann_P1_P2_TN <- length(which(Z_ann_all$var_is_1 != 'A5, A6' & Z_ann_all$`11` == 0)) 
  # multiGPATree_perf_data$ann_P1_P2_FP <- length(which(Z_ann_all$var_is_1 == 'A5, A6' & Z_ann_all$`11` == 0))
  # multiGPATree_perf_data$ann_P1_P2_FN <- length(which(Z_ann_all$var_is_1 != 'A5, A6' & Z_ann_all$`11` == 1)) 
  # multiGPATree_perf_data$ann_P1_P2_sensitivity <- multiGPATree_perf_data$ann_P1_P2_TP/(multiGPATree_perf_data$ann_P1_P2_TP + multiGPATree_perf_data$ann_P1_P2_FN)
  # multiGPATree_perf_data$ann_P1_P2_specificity <- multiGPATree_perf_data$ann_P1_P2_TN/(multiGPATree_perf_data$ann_P1_P2_TN + multiGPATree_perf_data$ann_P1_P2_FP)
  
  # multiGPATree_perf_data
  print('end: multiGPATree method')
  print('start: performance plot for LPM')
  
  
  # LPM implementation ####
  print('start: LPM')
  LPM_out <- aim2_LPM_result(gwasPval = simulated_data$gwasPval, 
                             annMat = simulated_data$annMat, 
                             Ztrue = simulated_data$Ztrue)
  print('end: LPM')
  
  # GPATree implementation ####
  print('start: GPATree')
  gpatree_out <- aim2_gpatree_result(gwasPval = simulated_data$gwasPval, 
                                     annMat = simulated_data$annMat, 
                                     Ztrue = simulated_data$Ztrue,
                                     initAlpha = initAlpha, cpTry = cpTry)
  print('end: GPATree')
  # things to return ####
  
  simout <- list()
  simout$M <- nrow(simulated_data$gwasPval)
  simout$nGWAS <- ncol(simulated_data$gwasPval)
  simout$nAnn <- ncol(simulated_data$annMat)
  simout$percent_ones_ann <- percent_ones_ann
  simout$percent_overlap_ann <- percent_overlap_ann
  simout$percent_overlap_gwas <- percent_overlap_gwas
  simout$P1_alpha_true <- trueAlphaVec[1]
  simout$P2_alpha_true <- trueAlphaVec[2]
  simout$multiGPATree_numIterConvg <- multigpatree_class_obj@fit$numIterConvergence
  simout$multiGPATree_alpha_P1 <- multigpatree_class_obj@fit$alpha[[1]]
  simout$multiGPATree_alpha_P2 <- multigpatree_class_obj@fit$alpha[[2]]
  simout$multiGPATree_fitSelectVar <- multigpatree_class_obj@fit$fitSelectVar
  simout$cpTry <- cpTry
  simout$multiGPATree_time <- multiGPATree_time
  simout[(length(simout)+1):(length(simout)+length(multiGPATree_perf_data))] <- unlist(multiGPATree_perf_data)
  names(simout)[(length(simout)-length(multiGPATree_perf_data)+1):length(simout)] <- names(multiGPATree_perf_data)
  
  simout[(length(simout)+1):(length(simout)+length(LPM_out))] <- unlist(LPM_out)
  names(simout)[(length(simout)-length(LPM_out)+1):length(simout)] <- names(LPM_out)
  
  simout[(length(simout)+1):(length(simout)+length(gpatree_out))] <- unlist(gpatree_out)
  names(simout)[(length(simout)-length(gpatree_out)+1):length(simout)] <- names(gpatree_out)
  
  return(simout)
  
}



