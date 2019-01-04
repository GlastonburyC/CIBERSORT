library(ggplot2)
bmi=read.table('All.Genes.CPM.txt',head=F)
macro=read.table('All.Genes.CPM.MacroAdj.txt',head=F)

all=merge(bmi,macro,by=1)

bon=0.05/dim(bmi)[1]

BMI_TWAS <- ggplot(all,aes(V3.x,V3.y)) +
      geom_point(colour="steelblue",size=2,alpha=0.6) +
      xlim(c(0,100)) +
      ylim(c(0,100)) +
      labs(x="Unadjusted BMI TWAS",size=10) +
      labs(y="Macrophage Adjusted BMI TWAS",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('BMI TWAS adjustment') + 
      geom_hline(yintercept=-log10(bon),linetype="dotted") +
      geom_vline(xintercept=-log10(bon),linetype="dotted")

BMI_TWAS
######################################################################################################

library(stats);library(ggplot2)

ests=read.table('Adipose.CPM.CellEsts.noMSC.txt',head=T)
adipose=read.table('TwinsUK/CPM.Adipose.90per.txt',head=T)

adipo_PC = prcomp(t(adipose[,2:ncol(adipose)]),scale=T)

qc=read.table('TwinsUK/qc_F_freezev2.txt',head=T)
qc$EurobatsID=gsub("EB","",qc$EurobatsID)

dxa=read.table('2006_09_eurobats_DXA_freeze.txt',head=T,sep="\t")

all = data.frame(adipo=ests$Adipocytes_SRR1296133, macro=ests$M1_SRR2910670+ests$M2_SRR2939149,BMI=qc$BMI,PC2=adipo_PC$x[,2])

# Filter cell types based on distance from mean.
M2_sd2=all[-which(all$macro > (mean(all$macro)+sd(all$macro)*2)),]

M2_sd2$sample=gsub("EB","",row.names(M2_sd2))

M2_sd2=M2_sd2[,c(5,1,2,3,4)]

M2_sd2=merge(M2_sd2,qc,by=1)

M2_sd2=M2_sd2[mixedorder(M2_sd2$sample),]


adipo_sd=all[-which(all$adipo < (mean(all$adipo)-sd(all$adipo)*2)),]

adipo_sd$sample=gsub("EB","",row.names(adipo_sd))

adipo_sd=adipo_sd[,c(5,1,2,3,4)]

adipo_sd=merge(adipo_sd,qc,by=1)

adipo_sd=adipo_sd[mixedorder(adipo_sd$sample),]


# merge macrophage dataframe with DXA variables
macro_dxa=merge(M2_sd2,dxa,by=1)
macro_dxa$AG_ratio=macro_dxa$ANDROID_FAT/macro_dxa$GYNOID_FAT

adipo_dxa=merge(adipo_sd,dxa,by=1)
adipo_dxa$AG_ratio=adipo_dxa$ANDROID_FAT/adipo_dxa$GYNOID_FAT

PC_celltype <- ggplot(M2_sd2,aes(macro,PC2)) +
      geom_point(colour="green4",size=2,alpha=0.6) +
      #xlim(c(0,100)) +
      #ylim(c(0,100)) +
      labs(x="Macrophage proportion",size=10) +
      labs(y="TwinsUK Adipose PC2",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('Macrophage proportion vs PC2')

pdf('Macro.vs.adiposePC2.pdf')
PC_celltype


##############################################################################################################


BMI_vs_Macro <- ggplot(M2_sd2,aes(BMI,macro)) +
      geom_point(colour="red4",size=2,alpha=0.6) +
      #xlim(c(0,100)) +
      #ylim(c(0,100)) +
      labs(x="BMI",size=10) +
      labs(y="Macrophage proportion",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('Macrophage proportion vs BMI')

pdf('BMI.vs.Macro.pdf')
BMI_vs_Macro


AG_ratio_vs_Macro <- ggplot(macro_dxa,aes(AG_ratio,macro)) +
      geom_point(colour="red4",size=2,alpha=0.6) +
      #xlim(c(0,100)) +
      #ylim(c(0,100)) +
      labs(x="Android/Gynoid Ratio",size=10) +
      labs(y="Macrophage proportion",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('Macrophage proportion vs AG_ratio')

pdf('BMI.vs.AG_ratio.pdf')
AG_ratio_vs_Macro
dev.off()
