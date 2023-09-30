rm(list = ls())
# install.packages("/home/khatia/dissertation/multiGPATree_0.0.0.9000.tar.gz", repos = NULL, type = "source")
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

datfinal <- fread("multiGPATree_application_data.csv")
dim(datfinal)
head(datfinal)
# Fig 3: histogram  and manhattan plots ####
jpeg(filename = 'Fig3_hist_manhattan_plot.jpeg', width=800, height = 540, quality = 100)
par(mfrow=c(2,4))
## histogram plots ####
hist(datfinal$SLE, breaks = 100, xlab = 'p-values',  main = 'SLE')
hist(datfinal$RA, breaks = 100, xlab = 'p-values',  main = 'RA')
hist(datfinal$UC, breaks = 100, xlab = 'p-values',  main = 'UC')
hist(datfinal$CD, breaks = 100, xlab = 'p-values',  main = 'CD')
## manhattan plot ####
man1 <- manhattan(datfinal, chr="chr_num", bp="position", snp="rsid", p="SLE" ,cex=0.4, cex.axis=0.8,suggestiveline = FALSE)
man2 <- manhattan(datfinal, chr="chr_num", bp="position", snp="rsid", p="RA" ,cex=0.4, cex.axis=0.8,suggestiveline = FALSE )
man3 <- manhattan(datfinal, chr="chr_num", bp="position", snp="rsid", p="UC"  ,cex=0.4, cex.axis=0.8,suggestiveline = FALSE)
man4 <- manhattan(datfinal, chr="chr_num", bp="position", snp="rsid", p="CD"  ,cex=0.4, cex.axis=0.8,suggestiveline = FALSE)
dev.off()

#Fig 4: 3 plots to show characteristics GS annotations ####
names(datfinal)
pdat <- datfinal %>% 
  select(., c(8:14)) %>%
  dplyr::rename(Blood = GS_Blood,
         Brain = GS_Brain,
         Epithelium = GS_Epithel,
         GI = GS_GI,
         Heart = GS_Heart,
         Lung = GS_Lung,
         Muscle = GS_Muscle) %>%
  mutate(total_tiss = rowSums(.))
  #492,557 SNPs included 

### plot1: proportion of SNPs annotated for specific number of annotations ####
ptab <- as.data.frame(prop.table(table(pdat$total_tiss)))
ptab$percent <- round(ptab$Freq*100, 2)

ggpropsnps <- ggplot(ptab, aes(x= Freq, y=Var1, fill='orange')) +
  geom_bar(stat="identity", width = 0.75) + 
  geom_text(aes(label= paste(percent, "%", sep='')), 
            size = 1.1, 
            position=position_dodge(width=0.5), 
            hjust=-0.05) + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, vjust=0.5, colour = 'black',angle=45),
        axis.text.y=element_text(colour = 'black'),
        text = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(size = .15), 
        axis.ticks.y = element_line(size = .15),
        axis.line = element_line(colour = "black", size = .15)) +
  labs(x = "Proportion of SNP", y = "Number of GS annotations") +
  xlim(c(0, 0.81))
ggpropsnps


### plot2: proportion of functional SNP based on annotations ####

bardat <- as.data.frame(cbind('tissue' = colnames(pdat)[1:(ncol(pdat)-1)], 
                              'prop' = round(colMeans(pdat[,1:(ncol(pdat)-1)]), 2), 
                              'percent' = round(colMeans(pdat[, 1:(ncol(pdat)-1)])*100, 2) ))
bardat$prop <- round(as.numeric(bardat$prop), 2)
bardat$percent <- as.numeric(bardat$percent)
bardat$tissue <- as.factor(bardat$tissue)
bardat$tissue <- factor(bardat$tissue, levels = colnames(pdat)[1:(ncol(pdat)-1)])
options(digits = 2)
ggpropfunc <- ggplot(bardat, aes(x= prop, 
                                 y = tissue,
                                 # y=stringr::str_wrap(tissue, 45), 
                                 fill=tissue)) +
  geom_bar(stat="identity", width = 0.85) + 
  geom_text(aes(label= paste(percent, "%", sep='')), 
            size = 1.5, 
            position=position_dodge(width=0.1), 
            hjust=-0.05) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, vjust=0.5, colour = 'black',angle = 45),
        axis.text.y=element_text(colour = 'black'),
        text = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_line(size = .15), 
        axis.ticks.y = element_line(size = .15),
        axis.line = element_line(colour = "black", size = 0.15)) +
  labs(x = "Proportion of SNP", y = "") +
  xlim(c(0, 0.18))
ggpropfunc

### plot3: overlap between annotations ####

logORtab <- REtab <- matrix(0, nrow = ncol(pdat)-1, ncol = ncol(pdat)-1)
rownames(logORtab) <- colnames(logORtab) <- rownames(REtab) <- colnames(REtab) <- colnames(pdat)[1:ncol(pdat)-1]
class(pdat)
pdat <- as.data.frame(pdat)
for (i in 1:7) {
  for (j in 1:7) {
    
    # i = 1
    # j = 1
    tab <- table(pdat[, i], pdat[, j])
    
    # relative enrichment
    print(tab)
    REtab[i, j] <- (tab[2,2]/(tab[2,2] + tab[2,1]))/(tab[1,2]/(tab[1,2] + tab[1,1])) 
    print(REtab[i, j])
    
    # log OR
    logORtab[i, j] <- log((tab[2,2]/(tab[2,1]))/(tab[1,2]/(tab[1,1])) )
    print(logORtab[i, j])
    
  }
}

logORtab[is.infinite(logORtab)] <- NA
co <- melt(logORtab)
gglogOR <- ggplot(co, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value)) + 
  geom_text(aes(fill = co$value, 
                label = round(co$value, 2)),
            size=1.5) +
  scale_fill_gradient2(low = "blue", 
                       high = "red",
                       mid = "white",
                       midpoint=mean(co$value, na.rm = TRUE)) +
  xlab('') +
  ylab('') +
  labs(fill = 'log(OR)') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, colour = 'black',angle=45,hjust=0.75),
        axis.text.y=element_text(colour = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_line(size = .15), 
        axis.ticks.y = element_line(size = .15),
        legend.key.width = unit(0.25, "cm"),
        axis.line = element_line(colour = "black", size = 0.15))
gglogOR
ggtosave <- ggpropsnps + ggpropfunc+ gglogOR +  
  # plot_layout(ncol = 3, widths = c(.75, .75, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7),
        axis.title.x = element_text(size=6)) 
ggtosave
ggsave(filename = 'Fig4_GS_EDA.png', plot = ggtosave, dpi = 300, 
       scale = 1, limitsize = T, height = 2, width = 5)

#Fig 5: 3 plots to show characteristics GSP annotations ####
names(datfinal)
pdat <- datfinal %>% 
  dplyr::select(., c(15:24)) %>%
  dplyr::rename('Helper memory T' = `GSP_Tcells_helperMemory`,
         'Helper naive T' = `GSP_Tcells_helperNaive`,
         'Effector/memory enriched T' =  `GSP_Tcells_effectorMemory`,
         'Regulatory T' = `GSP_Tcells_regulatory`,
         'CD8+ naive T' = `GSP_Tcells_CD8plusNaive`,             
         'CD8+ memory T' = `GSP_Tcells_CD8plusMemory`,
         'Monocytes' = `GSP_monocytes`,                     
         'Neutrophils' = `GSP_neutrophils`,
         'Primary B' = `GSP_Bcells`,
         'Natural killers' = `GSP_NaturalKiller`) %>%
  mutate(total_tiss = rowSums(.))
dim(pdat) #492,557 SNPs included 


### plot1: proportion of SNPs annotated for specific number of annotations ####
ptab <- as.data.frame(prop.table(table(pdat$total_tiss)))
ptab$percent <- round(ptab$Freq*100, 2)

ggpropsnps <- ggplot(ptab, aes(x= Freq, y=Var1, fill='orange')) +
  geom_bar(stat="identity", width = 0.75) + 
  geom_text(aes(label= paste(percent, "%", sep='')), 
            size = 1.1, 
            position=position_dodge(width=0.5), 
            hjust=-0.05) + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, vjust=0.5, colour = 'black',angle=45),
        axis.text.y=element_text(colour = 'black'),
        text = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(size = .15), 
        axis.ticks.y = element_line(size = .15),
        axis.line = element_line(colour = "black", size = .15)) +
  labs(#x = "Proportion of SNP", 
        x = '',
       y = "Number of GSP annotations") +
  xlim(c(0, 1))
ggpropsnps


### plot2: proportion of functional SNP based on annotations ####

bardat <- as.data.frame(cbind('tissue' = colnames(pdat)[1:(ncol(pdat)-1)], 
                              'prop' = round(colMeans(pdat[,1:(ncol(pdat)-1)]), 2), 
                              'percent' = round(colMeans(pdat[, 1:(ncol(pdat)-1)])*100, 2) ))
bardat$prop <- round(as.numeric(bardat$prop), 2)
bardat$percent <- as.numeric(bardat$percent)
bardat$tissue <- as.factor(bardat$tissue)
bardat$tissue <- factor(bardat$tissue, levels = colnames(pdat)[1:(ncol(pdat)-1)])
bardat
options(digits = 2)
ggpropfunc <- ggplot(bardat, aes(x= prop, 
                                 y = tissue,
                                 # y=stringr::str_wrap(tissue, 45), 
                                 fill=tissue)) +
  geom_bar(stat="identity", width = 0.85) + 
  geom_text(aes(label= paste(percent, "%", sep='')), 
            size = 1.1, 
            position=position_dodge(width=0.1), 
            hjust=-0.05) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, vjust=0.5, colour = 'black',angle = 45),
        axis.text.y=element_text(colour = 'black'),
        text = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_line(size = .15), 
        axis.ticks.y = element_line(size = .15),
        axis.line = element_line(colour = "black", size = 0.15)) +
  labs(#x = "Proportion of SNP", 
    x = '', 
    y = "") +
  xlim(c(0, 0.1))
ggpropfunc

### plot3: overlap between annotations ####

logORtab <- REtab <- matrix(0, nrow = ncol(pdat)-1, ncol = ncol(pdat)-1)
rownames(logORtab) <- colnames(logORtab) <- rownames(REtab) <- colnames(REtab) <- colnames(pdat)[1:ncol(pdat)-1]
class(pdat)
pdat <- as.data.frame(pdat)
for (i in 1:10) {
  for (j in 1:10) {
    
    # i = 1
    # j = 1
    tab <- table(pdat[, i], pdat[, j])
    
    # relative enrichment
    print(tab)
    REtab[i, j] <- (tab[2,2]/(tab[2,2] + tab[2,1]))/(tab[1,2]/(tab[1,2] + tab[1,1])) 
    print(REtab[i, j])
    
    # log OR
    logORtab[i, j] <- log((tab[2,2]/(tab[2,1]))/(tab[1,2]/(tab[1,1])) )
    print(logORtab[i, j])
    
  }
}

logORtab[is.infinite(logORtab)] <- NA
co <- melt(logORtab)
gglogOR <- ggplot(co, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value)) + 
  geom_text(aes(fill = co$value, 
                label = round(co$value, 1)),
            size=1.1) +
  scale_fill_gradient2(low = "blue", 
                       high = "red",
                       mid = "white",
                       midpoint=mean(co$value, na.rm = TRUE)) +
  xlab('') +
  ylab('') +
  labs(fill = 'log(OR)') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, colour = 'black',angle=45,vjust=1, hjust = 1),
        axis.text.y=element_text(colour = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_line(size = .15), 
        axis.ticks.y = element_line(size = .15),
        legend.key.width = unit(0.25, "cm"),
        axis.line = element_line(colour = "black", size = 0.15))
gglogOR
ggtosave <- ggpropsnps + ggpropfunc+ gglogOR +  
  # plot_layout(ncol = 3, widths = c(.75, .75, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        text = element_text(size = 4),
        axis.title.x = element_text(size=3)) 
ggtosave
ggsave(filename = 'Fig5_GSP_EDA.png', plot = ggtosave, dpi = 300, 
       scale = 1, limitsize = T, height = 2, width = 5)

