rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
setwd("/Users/khatiwadaa/Documents/2. MUSC/Dissertation/Aim 2/Results/simulation study")
t1 <- fread('setting3_simstudy_aim2_tree_case_09152023.txt')
table(t1$percent_overlap_ann, t1$P2_alpha_true) # only parameters that are changing
table(t1$multiGPATree_fitSelectVar)
table(t1$M)
table(t1$nGWAS)
table(t1$percent_overlap_ann)
# store curated simulation results data for manuscript submission ####
unique(t1$nAnn)
unique(t1$percent_ones_ann)
unique(t1$percent_overlap_ann)
unique(t1$percent_overlap_gwas)
names(t1)
sim_results <- t1 %>%
  dplyr::select(., c(P1_alpha_true, P2_alpha_true, percent_overlap_ann, 
                     multiGPATree_alpha_P1, multiGPATree_alpha_P2, multiGPATree_fitSelectVar,
                     multiGPATree_power_lFDRcontrol0.20_P1, multiGPATree_power_lFDRcontrol0.20_P2,
                     multiGPATree_power_lFDRcontrol0.20_P1_P2,multiGPATree_predlFDR_lFDRcontrol0.20_P1,
                     multiGPATree_predlFDR_lFDRcontrol0.20_P2, multiGPATree_predlFDR_lFDRcontrol0.20_P1_P2,
                     multiGPATree_AUC_P1_lfdr, multiGPATree_AUC_P2_lfdr, multiGPATree_AUC_P1_P2_lfdr,
                     LPM_alpha_P1, LPM_alpha_P2, LPM_power_lFDRcontrol0.2_P1, LPM_power_lFDRcontrol0.2_P2,
                     LPM_power_lFDRcontrol0.2_P1_P2, LPM_predlFDR_lFDRcontrol0.2_P1,
                     LPM_predlFDR_lFDRcontrol0.2_P2, LPM_predlFDR_lFDRcontrol0.2_P1_P2,
                     LPM_AUC_P1_lfdr, LPM_AUC_P2_lfdr, LPM_AUC_P1_P2_lfdr,
                     gpatree_alpha_P1, gpatree_alpha_P2,
                     gpatree_power_lFDRcontrol0.2_P1, gpatree_power_lFDRcontrol0.2_P2,
                     gpatree_predlFDR_lFDRcontrol0.2_P1, gpatree_predlFDR_lFDRcontrol0.2_P2, 
                     gpatree_AUC_P1_lfdr, gpatree_AUC_P2_lfdr)
  )
save(sim_results, file = "setting3_multiGPATree_simulation_results_data_09172023.RData")


datfile <- sim_results
# 1.1 A1-A6 simultaneously selected ####
d <- datfile[,c('P2_alpha_true', "percent_overlap_ann", "multiGPATree_fitSelectVar")]
d$all6ann_sel <- ifelse(d$multiGPATree_fitSelectVar == 'A1, A2, A3, A4, A5', 1, 0)

propt <- d %>% count(P2_alpha_true, percent_overlap_ann,all6ann_sel) %>% 
group_by(P2_alpha_true,percent_overlap_ann) %>%
mutate(prop_all6ann = n / sum(n))
names(propt)
tt <- as.data.frame(table(propt$P2_alpha_true, 
                        propt$percent_overlap_ann, 
                        propt$all6ann_sel))
head(tt)
tt <- tt[ ,1:3]
names(tt)[1:3] <- c('P2_alpha_true', 'percent_overlap_ann', 'all6ann_sel')
tt[, 1:3] <- lapply(lapply(tt[, 1:3], as.character), as.numeric) # change to character class and then to numeric

proptnew <- tt %>% full_join(propt, by = c('P2_alpha_true', 'percent_overlap_ann', 'all6ann_sel'))
proptnew[is.na(proptnew)==TRUE] <- 0
# propt <- subset(propt, prop == 1 )
propt2 <- subset(proptnew, all6ann_sel == 1)
propt1 <- propt2[, c('P2_alpha_true', 'percent_overlap_ann', 'prop_all6ann')]

# 1.2. % of noise annotations among selected annotations, % of true annotations (A1-A6) that are selected ####

dnew1 <- datfile[, c("P2_alpha_true", "percent_overlap_ann", "multiGPATree_fitSelectVar")]
dnew1$total_num_ann_select <- dnew1$wrong_num_ann_select <- dnew1$FDR_ann <- dnew1$TDR_ann <-rep(NA, nrow(dnew1))

for (i in 1:nrow(dnew1)) {
# i = 1
ann_select <- as.numeric(stringr::str_extract_all(dnew1$multiGPATree_fitSelectVar[i],"\\(?[0-9]+\\)?")[[1]])
ann_select <- ann_select[order(ann_select, decreasing = FALSE)]
# ann_select
dnew1$total_num_ann_select[i] <- length(ann_select)
dnew1$wrong_num_ann_select[i] <- length(which(ann_select > 5))
dnew1$FDR_ann[i] <- dnew1$wrong_num_ann_select[i]/dnew1$total_num_ann_select[i] 
dnew1$TDR_ann[i] <- (dnew1$total_num_ann_select[i] - dnew1$wrong_num_ann_select[i])/5
}
names(dnew1)
propt <- dnew1 %>% group_by(P2_alpha_true, percent_overlap_ann) %>% 
summarize(mean_FDR = mean(FDR_ann, na.rm=TRUE),mean_TDR = mean(TDR_ann, na.rm=TRUE))


if ( nrow(propt1) > 0 ) {
merge_data <- merge(propt1, propt, by = c("P2_alpha_true", "percent_overlap_ann"))
} else {
merge_data <- cbind(propt[, 1:2], 'prop_all6ann' = 0, propt[, 3:4])
}  

# change data from wide to long format for plotting the three lines together

# 1.3. 3 annotation selection plots combined into one ####
merge_data_l <- merge_data %>% pivot_longer(cols=c('prop_all6ann', 'mean_FDR', 'mean_TDR'),
                                          names_to='criteria',
                                          values_to='prop') %>%
mutate(criteria = recode(criteria, 
                         mean_FDR = "Mean % of noise annotations (A6-A25) \namong selected annotations", 
                         prop_all6ann = "% of simulated data for which A1-A5 \nare simultaneously selected ", 
                         mean_TDR = "Mean % of true annotations (A1-A5) \nthat are selected"))
table(merge_data_l$criteria)

p410 <- ggplot(merge_data_l, aes(x = as.factor(percent_overlap_ann), 
                            y = prop, 
                            group = as.factor(criteria),
                            color = as.factor(criteria)
)) + 
ylim(0, 1) +
geom_line(size = 0.30)+
geom_point(size = 0.35) +
xlab("% overlap in annotated SNPs (v)") +
ggtitle("") + 
ylab ('Proportion') +
theme_bw() +
theme(axis.text.x=element_text(angle=50, size=8, vjust=0.5),
      panel.spacing.x = unit(0.5, "lines"),
      axis.text.y=element_text(size = 8),
      strip.text.x = element_text(size = 8),
      # plot.margin = unit(c(1,1,1,1), "cm"),
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.position = 'bottom') + 
guides(color = guide_legend(nrow = 3,
                            title=""))
p410
  
# 2. alpha 1 and alpha 2 combined plot ####
alpha1 <- alpha2 <- 0.4
data_sub_alpha <- datfile[,c("percent_overlap_ann", 
                           "multiGPATree_alpha_P1", "LPM_alpha_P1", "gpatree_alpha_P1",
                           "multiGPATree_alpha_P2", "LPM_alpha_P2", "gpatree_alpha_P2")] 
names(data_sub_alpha)
data_sub_alpha$multiGPATree_alpha_P1_P2 <- 2.1
data_sub_alpha$LPM_alpha_P1_P2 <- 2.1
data_sub_alpha$gpatree_alpha_P1_P2 <- 2.1
head(data_sub_alpha)
data_sub_alpha_long <- melt(data_sub_alpha, 
                          id.vars = c("percent_overlap_ann"), 
                          measure.vars = c("multiGPATree_alpha_P1", "LPM_alpha_P1", "gpatree_alpha_P1",
                                           "multiGPATree_alpha_P2", "LPM_alpha_P2", "gpatree_alpha_P2",
                                           "multiGPATree_alpha_P1_P2", "LPM_alpha_P1_P2", "gpatree_alpha_P1_P2"), 
                          variable.name = 'Method', value.name = 'est_alpha')
levels(data_sub_alpha_long$Method)
data_sub_alpha_long <- data_sub_alpha_long %>%
mutate(Method = recode(Method,
                       multiGPATree_alpha_P1 = "multiGPATree_P1", 
                       LPM_alpha_P1 = "LPM_P1",
                       gpatree_alpha_P1 = "GPATree_P1",
                       multiGPATree_alpha_P2 = "multiGPATree_P2", 
                       LPM_alpha_P2 = "LPM_P2",
                       gpatree_alpha_P2 = "GPATree_P2",
                       GPATree_alpha_P1_P2 = 'multiGPATree_P1_P2',
                       LPM_alpha_P1_P2 = "LPM_P1_P2",
                       gpatree_alpha_P1_P2 = "GPATree_P1_P2",))


levels(data_sub_alpha_long$Method)

malpha <- ggplot(data_sub_alpha_long, aes(x = as.factor(percent_overlap_ann), 
                                    y = est_alpha, 
                                    fill = Method)) +
geom_boxplot(size=0.15, outlier.size = 0.001, key_glyph = "boxplot") + 
ylim(min(data_sub_alpha_long$est_alpha, alpha1)+0.05, 1) +
geom_hline(yintercept = alpha1, linetype = 'dashed', color = 'Chartreuse', size = 0.5) +
theme_bw() + 
theme(axis.text.x=element_text(angle=50, size=8, vjust=0.5),
      legend.title = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      axis.text.y=element_text(size = 8),
      plot.margin = unit(c(0,1,0,0), "cm"),
      legend.position = 'bottom') +
scale_fill_manual(values=c("pink", "lightblue", 'green', "red", "blue", 'darkgreen'),
                  name = "") +
guides(fill=guide_legend(nrow=3)) +
xlab("") +
ylab(expression(hat(alpha)))
# labs(fill='') +
malpha  

  
  
# 3. AUC for P1 and P2 combined plot ####
names(datfile)
data_sub_alpha <- datfile[,c("percent_overlap_ann",
                           "multiGPATree_AUC_P1_lfdr", "LPM_AUC_P1_lfdr","gpatree_AUC_P1_lfdr",
                           "multiGPATree_AUC_P2_lfdr", "LPM_AUC_P2_lfdr","gpatree_AUC_P2_lfdr",
                           "multiGPATree_AUC_P1_P2_lfdr", "LPM_AUC_P1_P2_lfdr")]

data_sub_alpha_long <- melt(data_sub_alpha, 
                          id.vars = c("percent_overlap_ann"), 
                          measure.vars = c("multiGPATree_AUC_P1_lfdr", "LPM_AUC_P1_lfdr","gpatree_AUC_P1_lfdr",
                                           "multiGPATree_AUC_P2_lfdr", "LPM_AUC_P2_lfdr","gpatree_AUC_P2_lfdr",
                                           "multiGPATree_AUC_P1_P2_lfdr", "LPM_AUC_P1_P2_lfdr"), 
                          variable.name = 'Method', value.name = 'est_alpha')
data_sub_alpha_long <- data_sub_alpha_long %>%
mutate(Method = recode(Method,
                       multiGPATree_AUC_P1_lfdr = "multiGPATree_P1",
                       multiGPATree_AUC_P2_lfdr = "multiGPATree_P2",
                       multiGPATree_AUC_P1_P2_lfdr = "multiGPATree_P1_P2",
                       LPM_AUC_P1_lfdr = "LPM_P1",
                       LPM_AUC_P2_lfdr = "LPM_P2",
                       LPM_AUC_P1_P2_lfdr = "LPM_P1_P2",
                       gpatree_AUC_P1_lfdr = "GPATree_P1",
                       gpatree_AUC_P2_lfdr = "GPATree_P2"
))
table(data_sub_alpha_long$Method)

pauc <- ggplot(data_sub_alpha_long, aes(x = as.factor(percent_overlap_ann), 
                                      y = est_alpha, 
                                      fill = Method)) +
geom_boxplot(size=0.15, outlier.size = 0.001, key_glyph = "boxplot") + 
ylim(min(data_sub_alpha_long$est_alpha)-0.005, 1) +
# xlim(0.83, 1) +
theme_bw() + 
scale_fill_manual(values=c("pink", "lightblue", "green", "red", "blue", "darkgreen", 'brown4', 'navyblue'),
                  name = "") +
  guides(fill=guide_legend(nrow=3)) +
theme(axis.text.x=element_text(angle=50, size=8, vjust=0.5),
      legend.title = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      axis.text.y=element_text(size = 8),
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.position = 'bottom') +
xlab("") +
ylab('AUC')+ 
labs(fill='')
pauc  
  
  
# 4. Power for P1 and P2 when local FDR controlled at 0.20 combined plot ####
names(datfile)
data_sub_alpha <- datfile[,c("percent_overlap_ann",
                           "multiGPATree_power_lFDRcontrol0.20_P1", 
                           "LPM_power_lFDRcontrol0.2_P1","gpatree_power_lFDRcontrol0.2_P1",
                           "multiGPATree_power_lFDRcontrol0.20_P2", "LPM_power_lFDRcontrol0.2_P2", "gpatree_power_lFDRcontrol0.2_P2",
                           "multiGPATree_power_lFDRcontrol0.20_P1_P2", "LPM_power_lFDRcontrol0.2_P1_P2")]

data_sub_alpha_long <- melt(data_sub_alpha, 
                          id.vars = c("percent_overlap_ann"), 
                          measure.vars = c("multiGPATree_power_lFDRcontrol0.20_P1", 
                                           "LPM_power_lFDRcontrol0.2_P1",
                                           "gpatree_power_lFDRcontrol0.2_P1",
                                           "multiGPATree_power_lFDRcontrol0.20_P2", 
                                           "LPM_power_lFDRcontrol0.2_P2",
                                           "gpatree_power_lFDRcontrol0.2_P2",
                                           "multiGPATree_power_lFDRcontrol0.20_P1_P2", 
                                           "LPM_power_lFDRcontrol0.2_P1_P2"), 
                          variable.name = 'Method', value.name = 'est_alpha')
table(data_sub_alpha_long$Method)
data_sub_alpha_long <- data_sub_alpha_long %>%
mutate(Method = recode(Method,
                       multiGPATree_power_lFDRcontrol0.20_P1 = "multiGPATree_P1",
                       multiGPATree_power_lFDRcontrol0.20_P2 = "multiGPATree_P2",
                       multiGPATree_power_lFDRcontrol0.20_P1_P2 = "multiGPATree_P1_P2",
                       LPM_power_lFDRcontrol0.2_P1 = "LPM_P1",
                       LPM_power_lFDRcontrol0.2_P2 = "LPM_P2",
                       LPM_power_lFDRcontrol0.2_P1_P2 = "LPM_P1_P2",
                       gpatree_power_lFDRcontrol0.2_P1 = "GPATree_P1",
                       gpatree_power_lFDRcontrol0.2_P2 = "GPATree_P2"))
table(data_sub_alpha_long$Method)
names(data_sub_alpha_long)
ppower <- ggplot(data_sub_alpha_long, aes(x = as.factor(percent_overlap_ann), 
                                        y = est_alpha, 
                                        fill = Method)) +
geom_boxplot(size=0.15, outlier.size = 0.001, key_glyph = "boxplot") + 
ylim(min(data_sub_alpha_long$est_alpha)-0.1, 1) +
theme_bw() + 
theme(axis.text.x=element_text(angle=50, size=8, vjust=0.5),
      legend.title = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      axis.text.y=element_text(size = 8),
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.position = 'bottom') +
scale_fill_manual(values=c("pink", "lightblue", "green", "red", "blue", "darkgreen",  'brown4', 'navyblue'),
                  name = "") +
  guides(fill=guide_legend(nrow=3)) +
xlab("% overlap in annotated SNPs (v)") +
ylab(expression(Power[paste('lfdr=',0.20)]))+ 
labs(fill='')
ppower  
  
# 5. predicted localFDR for P1 and P2 when local FDR controlled at 0.20 combined plot ####
names(datfile)
data_sub_alpha <- datfile[,c("percent_overlap_ann",
                           "multiGPATree_predlFDR_lFDRcontrol0.20_P1",
                           "LPM_predlFDR_lFDRcontrol0.2_P1",
                           "gpatree_predlFDR_lFDRcontrol0.2_P1",
                           "multiGPATree_predlFDR_lFDRcontrol0.20_P2",
                           "LPM_predlFDR_lFDRcontrol0.2_P2",
                           "gpatree_predlFDR_lFDRcontrol0.2_P2",
                           "multiGPATree_predlFDR_lFDRcontrol0.20_P1_P2",
                           "LPM_predlFDR_lFDRcontrol0.2_P1_P2")]

data_sub_alpha_long <- melt(data_sub_alpha,
                          id.vars = c("percent_overlap_ann"),
                          measure.vars = c("multiGPATree_predlFDR_lFDRcontrol0.20_P1",
                                           "LPM_predlFDR_lFDRcontrol0.2_P1",
                                           "gpatree_predlFDR_lFDRcontrol0.2_P1",
                                           "multiGPATree_predlFDR_lFDRcontrol0.20_P2",
                                           "LPM_predlFDR_lFDRcontrol0.2_P2",
                                           "gpatree_predlFDR_lFDRcontrol0.2_P2",
                                           "multiGPATree_predlFDR_lFDRcontrol0.20_P1_P2",
                                           "LPM_predlFDR_lFDRcontrol0.2_P1_P2"),
                          variable.name = 'Method', value.name = 'est_alpha')
table(data_sub_alpha_long$Method)
data_sub_alpha_long <- data_sub_alpha_long %>%
mutate(Method = recode(Method,
                       multiGPATree_predlFDR_lFDRcontrol0.20_P1 = "multiGPATree_P1",
                       multiGPATree_predlFDR_lFDRcontrol0.20_P2 = "multiGPATree_P2",
                       multiGPATree_predlFDR_lFDRcontrol0.20_P1_P2 = "multiGPATree_P1_P2",
                       LPM_predlFDR_lFDRcontrol0.2_P1 = "LPM_P1",
                       LPM_predlFDR_lFDRcontrol0.2_P2 = "LPM_P2",
                       LPM_predlFDR_lFDRcontrol0.2_P1_P2 = "LPM_P1_P2",
                       gpatree_predlFDR_lFDRcontrol0.2_P1 = "GPATree_P1",
                       gpatree_predlFDR_lFDRcontrol0.2_P2 = "GPATree_P2"))
table(data_sub_alpha_long$Method)
names(data_sub_alpha_long)
plfdr <- ggplot(data_sub_alpha_long, aes(x = as.factor(percent_overlap_ann),
                                       y = est_alpha,
                                       fill = Method)) +
geom_boxplot(size=0.15, outlier.size = 0.001, key_glyph = "boxplot") +
geom_hline(yintercept = 0.20, linetype = 'dashed', color = 'Chartreuse', size = 0.5) +
ylim(min(data_sub_alpha_long$est_alpha), 0.25) +
theme_bw() +
theme(axis.text.x=element_text(angle=50, size=8, vjust=0.5),
      legend.title = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      axis.text.y=element_text(size = 8),
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.position = 'bottom') +
scale_fill_manual(values=c("pink", "lightblue", "green", "red", "blue", "darkgreen", 'brown4', 'navyblue'),
                  # guide_legend(nol = 3),
                  name = "") +
  guides(fill=guide_legend(nrow=3)) +
xlab("") +
ylab(expression(Predicted ~lfdr[paste('lfdr=',0.20)]))+
labs(fill='')
plfdr
  

# combine the 5 plots ####
ggall1 <-  {pauc + ppower + plfdr + 
    plot_annotation(tag_levels = 'A')+
    plot_layout(guides = 'collect',nrow = 1, ncol = 3) &
    # labs(x = '% overlap in annotated SNPs (v)') &
    theme_bw() +
    theme(legend.position = 'bottom',
          axis.text.x=element_text(angle=50, size=10, vjust=0.5),
          axis.text.y=element_text(size=10),
          legend.text = element_text(size = 10),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10))} 
ggall1  
ggall2 <- {malpha + p410 +
    plot_annotation(tag_levels = list(c('D', 'E', 'F'))) +
    plot_layout(ncol = 3, nrow=1) &
    theme_bw() +
    theme(legend.position = 'bottom',
          axis.text.x=element_text(angle=50, size=10, vjust=0.5),
          axis.text.y=element_text(size=10),
          legend.text = element_text(size = 10),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10)) }
ggall2
saveplot <- cowplot::plot_grid(ggall1, ggall2, nrow = 2)
saveplot

ggsave(# file ='setting1_simulationResult_09052023.png',
       file ='setting3_simulationResult_09172023.png',
       saveplot,
       width = 10.5, height = 7)

 

