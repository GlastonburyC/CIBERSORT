p = "BMI" #i.e. 'BMI'
tis <- "F" # i.e. 'F'
reads_dir <- "/home/DTR/Expression/EUROBATS/Counts"
out_dir <- "/home/glastonc/new_bams/counts/TWAS/"

inversenormal <- function(x) {

 # inverse normal if you have missing data

 return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))

}

calc_models <- function(tis, reads_dir, out_dir){
 #I split the read files so I can massively parralise the job. - split each tissue expression file into 100 files. Genome wide expression analysis for 120,000 in 700 individuals ~ 20mins.
file.p<-Sys.getenv("SGE_TASK_ID") # iterate through file names using the SUN GRID ENGINE ID   qsub -t ....
#file.p="X"
reads = read.table(file=paste("TWAS/split_counts/CPM.Adipose.90per.",file.p,".txt",sep=""),sep="\t", header=T)
covs = read.table(file=paste("TwinsUK/qc_",tis,"_freezev2.txt",sep=""),header=T,sep="",stringsAsFactors=F)

insert_size=read.table('TWAS/QC_metrics/all.insert.distance.txt',stringsAsFactors=F,head=T)
gc_mean=read.table('TWAS/QC_metrics/all.GC_mean.txt',stringsAsFactors=F,head=F)
pct_mismatch=read.table('TWAS/QC_metrics/all_pct.mismatch.txt',head=T,stringsAsFactors=F)
#ncb=read.table('/home/glastonc/new_bams/QC_metrics/NCB.estimates.txt',stringsAsFactors=F,head=T)

cells=read.csv('TwinsUK/Macrophage.ests.filtered.txt',head=T,sep="\t")

insert_size=insert_size[substring(insert_size[,1],3,10) %in% cells[,1],]
gc_mean=gc_mean[substring(gc_mean[,1],3,10) %in% cells[,1],]
pct_mismatch=pct_mismatch[substring(pct_mismatch[,1],3,10) %in% cells[,1],]
covs=covs[substring(covs[,1],3,10) %in% cells[,1],]


gene=reads[,1]
reads=reads[,substring(colnames(reads),3,10) %in% cells[,1]]

tmp=cbind(gene,reads)
reads=tmp

  ##########################################################################
  
  # Set all variables for model
  family <- as.factor(as.matrix(covs['Family'])) # This is coded with the family ID.
  zygosity <- as.factor(as.matrix(covs['Zygosity'])) #  Same number for MZ twins (family ID),  for DZ twins (twin ID + 00000)
  INSERT_SIZE_MODE <- as.numeric(as.matrix(insert_size[,2]))
  GC_mean <- as.numeric(as.matrix(gc_mean[,2]))
  pct_mismatch<- as.numeric(as.matrix(pct_mismatch[,2]))
 # ncb <- as.numeric(as.matrix(ncb[,2]))
  PrimerIndex <- as.factor(as.matrix(covs['PrimerIndex']))
  date <- as.factor(as.matrix(covs['date']))
  age <- as.numeric(as.matrix(covs['AGE']))
# batch <- as.factor(as.matrix(covs['Set'])) # ONLY BLOOD HAS BATCH EFFECT
  BMI <- as.numeric(as.matrix(covs['BMI']))
	
  M2=as.numeric(as.matrix(cells$macro))    
  if(tis == "F"){n_ids = dim(covs)[1]}

  ###################################

    sink(file=paste(out_dir,"/Genes.CPM_chr",file.p,"_",p,"_",tis,".tab", sep=""))
        #######################################################################
        run_models(reads, family, zygosity, INSERT_SIZE_MODE, GC_mean, PrimerIndex, date, age, n_ids, BMI, age_sq, Phenotype,ncb, pct_mismatch,M2)
        #######################################################################
    sink()
}

run_models <- function(reads, family, zygosity, INSERT_SIZE_MODE, GC_mean, PrimerIndex, date , age, n_ids, BMI, age_sq, Phenotype,ncb, pct_mismatch,M2){
  n_reads <- dim(reads)[1]
  reads_id <- NULL
     require(lme4)
     require(Matrix)
    i_reads <- NULL
    for(i_reads in 2:n_reads){
       reads_id <- reads[i_reads,c(1)]

       y_in = as.matrix(as.double(reads[i_reads,2:ncol(reads)]))
       y_in = inversenormal(y_in)
       Fm1 = Fm2 = NULL
          #######################################################################
          Fm1 <- lmer(y_in ~ 1 + scale(BMI) + scale(age) + scale(GC_mean) + scale(INSERT_SIZE_MODE) + (1 | date) + (1|family) + (1 | PrimerIndex) + (1|zygosity), REML=FALSE)
          Fm2 <- lmer(y_in ~ 1 + scale(age) + scale(GC_mean) + scale(INSERT_SIZE_MODE)  + (1|family) + (1 | date) + (1 | PrimerIndex) + (1|zygosity), REML=FALSE)
          #######################################################################
          
                                                                                # Output the exon name, Chisq value, p-value and all the fixed effect (betas) from model 1
if(i_reads ==1){
  cat(paste("Gene","ChiSq","pvalue",p,"age_beta","GC_beta","insert_beta","NCB_beta","pctMis_beta",sep="\t"))
}else{
cat(paste(as.character(reads[i_reads,1]),anova(Fm1,Fm2)$Chis[2],anova(Fm1,Fm2)$Pr[2], fixef(Fm1)[2],fixef(Fm1)[3],fixef(Fm1)[4]),fixef(Fm1)[5],fixef(Fm1)[6],fixef(Fm1)[7],"\n")
    }}
}
################################################################################################################################

calc_models(tis, reads_dir, out_dir)

