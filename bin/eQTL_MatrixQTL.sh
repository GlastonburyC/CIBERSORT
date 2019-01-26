chr <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) 
tis="F"
pheno="adipocytes"
library("MatrixEQTL");

useModel = modelLINEAR #  This model includes an interaction term for genotype vs covar

SNP_file_name = paste("/home/ubuntu/dosages/macrophages.chr",chr,"_1000g.F.dosage.Final.maf5",sep="");
snps_location_file_name = paste("/home/ubuntu/dosages/info/chr",chr,"_1000g.dosage.Final.maf5.snppos",sep="");

expression_file_name = paste("/home/ubuntu/cibersort_bams/TwinsUK/TMM.adipose.PEER.residuals.k30.invNorm.withGeno.Macro.txt",sep="");
gene_location_file_name = paste("/home/ubuntu/cibersort_bams/TwinsUK/gene.info.txt",sep="");

covariates_file_name = paste("/home/ubuntu/cibersort_bams/TwinsUK/macrophage.withGeno.MatrixQTL.txt",sep="");

output_file_name = paste("/home/ubuntu/cibersort_bams/TwinsUK/baseline_eqtl/cis_eqtl_",tis,"_chr",chr,".txt",sep="");
output_file_name.cis = paste("/home/ubuntu/cibersort_bams/TwinsUK/baseline_eqtl/cis_eqtl_",tis,"_chr",chr,".txt",sep="");
errorCovariance = numeric();
cisDist = 1e6;
pvOutputThreshold_cis = 1; 
pvOutputThreshold = 0

snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1 # no column of row labels
snps$fileSliceSize = 2000 # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name);

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;a
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

#cvrt = SlicedData$new();
#cvrt$fileDelimiter = "\t";      # the TAB character
#cvrt$fileOmitCharacters = "NA"; # denote missing values;
#cvrt$fileSkipRows = 1;          # one row of column labels
#cvrt$fileSkipColumns = 0;       # one column of row labels
#if(length(covariates_file_name)>0) {
#cvrt$LoadFile(covariates_file_name);
#}

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos$chr=paste0("chr",genepos$chr)

me = Matrix_eQTL_main(
                snps = snps, 
                #cvrt = cvrt,
                gene = gene, 
                useModel = useModel, 
                errorCovariance = errorCovariance, 
                verbose = TRUE, 
                output_file_name = output_file_name,
                output_file_name.cis = output_file_name.cis,
		pvOutputThreshold.cis = pvOutputThreshold_cis,
		pvOutputThreshold = pvOutputThreshold,
                snpspos = snpspos,      
                genepos = genepos,
                cisDist = cisDist,
                pvalue.hist = FALSE, #To make a hist of all p-values set
                noFDRsaveMemory = TRUE
    );
