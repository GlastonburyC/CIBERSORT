library(ggplot2)
library(grid)
library(gridExtra)

ests1=read.csv('inSilico.CellEsts.Contam.1per.txt',head=T,sep="\t")
ests10=read.csv('inSilico.CellEsts.Contam.10per.txt',head=T,sep="\t")
ests50=read.csv('inSilico.CellEsts.Contam.50per.txt',head=T,sep="\t")
truth1=read.table('1000.MixingProps.contamination.1per.txt',head=T)
truth10=read.table('1000.MixingProps.contamination.10per.txt',head=T)
truth50=read.table('1000.MixingProps.contamination.50per.txt',head=T)


Adipocytes1 <- ggplot(ests1,aes(ests1$Adipocytes_SRR1296133,truth1$adipo1)) +
      geom_point(colour="steelblue",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated Adipocytes",size=10) +
      labs(y="Ground Truth Adipocytes",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('1% contamination')



CD41 <- ggplot(ests1,aes(ests1$CD4_SRR1422906,truth1$tcell3)) +
      geom_point(colour="green4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated CD4+",size=10) +
      labs(y="Truth CD4+",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('1% contamination')



Macrophages1 <- ggplot(ests1,aes(ests1$M1_SRR2910670+ests1$M2_SRR2939149,truth1$macro3)) +
      geom_point(colour="turquoise4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated Macrophages",size=10) +
      labs(y="Truth Macrophages",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('1% contamination')



HMVEC1 <- ggplot(ests1,aes(ests1$HMVEC_SRR2776477,truth1$dermal_MVEC)) +
      geom_point(colour="red4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated MVEC",size=10) +
      labs(y="Truth MVEC",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('1% contamination')


#########

Adipocytes10 <- ggplot(ests10,aes(ests10$Adipocytes_SRR1296133,truth10$adipo1)) +
      geom_point(colour="steelblue",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated Adipocytes",size=10) +
      labs(y="Ground Truth Adipocytes",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('10% contamination')


CD410 <- ggplot(ests10,aes(ests10$CD4_SRR1422906,truth10$tcell3)) +
      geom_point(colour="green4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated CD4+",size=10) +
      labs(y="Truth CD4+",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('10% contamination')


Macrophages10 <- ggplot(ests10,aes(ests10$M1_SRR2910670+ests10$M2_SRR2939149,truth10$macro3)) +
      geom_point(colour="turquoise4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated Macrophages",size=10) +
      labs(y="Truth Macrophages",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('10% contamination')


HMVEC10 <- ggplot(ests10,aes(ests10$HMVEC_SRR2776477,truth10$dermal_MVEC)) +
      geom_point(colour="red4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated MVEC",size=10) +
      labs(y="Truth MVEC",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('10% contamination')


########

Adipocytes50 <- ggplot(ests50,aes(ests50$Adipocytes_SRR1296133,truth50$adipo1)) +
      geom_point(colour="steelblue",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated Adipocytes",size=10) +
      labs(y="Ground Truth Adipocytes",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('50% contamination')


CD450 <- ggplot(ests50,aes(ests50$CD4_SRR1422906,truth50$tcell3)) +
      geom_point(colour="green4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated CD4+",size=10) +
      labs(y="Truth CD4+",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('50% contamination')


Macrophages50 <- ggplot(ests50,aes(ests50$M1_SRR2910670+ests50$M2_SRR2939149,truth50$macro3)) +
      geom_point(colour="turquoise4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated Macrophages",size=10) +
      labs(y="Truth Macrophages",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('50% contamination')


HMVEC50 <- ggplot(ests50,aes(ests50$HMVEC_SRR2776477,truth50$dermal_MVEC)) +
      geom_point(colour="red4",size=2,alpha=0.6) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x="Estimated MVEC",size=10) +
      labs(y="Truth MVEC",size=10)+
      geom_smooth(method = "lm", se = TRUE,colour="black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")) + 
      ggtitle('50% contamination')

pdf('Simulations.withContaminates.all.pdf',width=8,height=8)


grid.arrange(Adipocytes1,Adipocytes10,Adipocytes50,HMVEC1,HMVEC10,HMVEC50,Macrophages1,Macrophages10,Macrophages50,CD41,CD410,CD450, ncol = 3)

dev.off()








pdf('Simulations.withContaminates.10per.pdf')
grid.arrange(Adipocytes, HMVEC,Macrophages,CD4, ncol = 2, top="10% contaminating cells")

dev.off()

mean(abs(ests$HMVEC_SRR2776477-truth$dermal_MVEC))
cor.test(ests$HMVEC_SRR2776477,truth$dermal_MVEC)

mean(abs((ests$M1_SRR2910670+ests$M2_SRR2939149)-truth$macro3))
cor.test(ests$M1_SRR2910670+ests$M2_SRR2939149,truth$macro3)

mean(abs(ests$CD4_SRR1422906-truth$tcell3))
cor.test(ests$CD4_SRR1422906,truth$tcell3)

mean(abs(ests$Adipocytes_SRR1296133-truth$adipo1))
cor.test(ests$Adipocytes_SRR1296133,truth$adipo1)
