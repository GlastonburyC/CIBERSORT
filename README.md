# CIBERSORT
## Using CIBERSORT to estimate cell type composition.

## TwinsUK

adapter/polyA trimming steps:
```
trim_galore -stringency 5 -q 1 -o '$line' --paired '$line'/'$line'.f1.fq '$line'/'$line'.f2.fq
perl prinseq-lite.pl -fastq '$line'/'$line'.f1_val_1.fq -fastq2 '$line'/'$line'.f2_val_2.fq -out_good '$line'/'$line' -trim_tail_left 5 -trim_tail_right 5 -min_len 20
```
STAR alignment:
```
./STAR --runThreadN '$THREAD_NO' --runMode alignReads --readFilesIn '$line'/'$line'_1.fastq '$line'/'$line'_2.fastq --genomeDir hg19 --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang 48 --outFilterMismatchNmax 6 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:'$line'_maternal PU:Illumina PL:Illumina LB:'$line'_maternal SM:'$line'_maternal CN:Seq_centre --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '$line'_ref.
```
read filtering:
```
samtools view -@ '$THREAD_NO' -b -F4 -q 30 '$line'/reference/'$line'_ref.Aligned.sortedByCoord.out.bam -o '$line'/reference/'$line'.filtered.bam
```
gene counts with Gencode v19
```
featureCounts -p -T '$THREAD_NO' -a gencode.v19.annotation.gtf -o '$line'/reference/'$line'.GeneCount_Ref.txt '$line'/reference/'$line'.filtered.bam
```
Gene counts then converted to CPM:
```
expr[,2:ncol(expr)]=apply(expr,2,function(x) (x/sum(x)) *1e6)
write.table(expr,"TwinsUK.adiposeSamples.txt",col.names=T,row.names=F,sep="\t",quote=F)
```

CIBERSORT used to estimate fractions:

``` 
R –-no-restore
library(Rserve)
Rserve()
q()
java -jar CIBERSORT.jar -M DeCODE.adiposeSamples.txt -B Adipose.sigMatrix.txt –n 1000 >> Decode.Adipose.CellEsts.txt
```


# Other cohorts
Thanks for using our signature matrix to estimate the cellular composition of your adipose tissue RNA-seq samples.

The adipose tissue signature matrix was constructed from datasets that are independent from TwinsUK and consist of purified cell type RNA-seq data obtained from macrophages, adipocytes, endothelial cells and CD4+ t-cells – These are thought to be some of the most dominant cell types present in adipose tissue biopsies.

The signature matrix was constructed using Counts Per Million (CPM) gene expression levels that have not been filtered for lowly expressed genes. It is important that your adipose tissue RNA-seq data is therefore in the same space. Each deconvolution happens on a sample by sample basis, therefore it is not necessary to filter genes that are expressed in at least 10% of samples. Additionally, rare cell types will be represented by lowly expressed genes relative to other genes expressed in more common cell types (linearity assumption of all cell estimation methods) so filtering lowly expressed genes could bias cell type estimates upwards. CIBERSORT implements support vector regression, so a gene that is not expressed in a mixture will not be selected as a maximally discriminating feature.

It is important you produce a raw gene count matrix obtained from your BAMs (we used featureCount for example) and transform them into CPM using the following R code:

```
expr[,2:ncol(expr)]=apply(expr,2,function(x) (x/sum(x)) *1e6)
write.table(expr,"DeCODE.adiposeSamples.txt",col.names=T,row.names=F,sep="\t",quote=F)
```

where: 

- expr is an NxM matrix with the first column being gene (ensembl ID) and all other columns being adipose RNA-seq samples (raw gene counts).

Cell type deconvolution is very sensitive to the method/normalization used and therefore if we observe vast differences in cell estimates I may reproduce the signature matrix using a pipeline closer to your RNA-seq study (i.e. the RNA-seq aligner / gencode annotation) and ask you to re-run.
 
## CIBERSORT to estimate cell fractions:
``` 
R –-no-restore
library(Rserve)
Rserve()
q()

java -jar CIBERSORT.jar -M DeCODE.adiposeSamples.txt -B Adipose.sigMatrix.txt –n 1000 >> Decode.Adipose.CellEsts.txt
```

Where:	

- DeCODE.adiposeSamples.txt - NxM matrix with the first column being gene (ensembl ID) and all other columns being samples.

- Adipose.sigMatrix.txt – signature matrix containing median CPM expression levels of genes used to estimate relative cell type proportions.

This will output the cell type estimates, with some additional information in the header (This can be discarded).

Please send these cell estimates to us so we can see if they are similar to our predictions. If that is the case, we ask you to continue this pipeline.
 
## Correlation of cell type estimates with obesity metrics:

Please test the following correlations with BMI both with and without removing outlier cell estimates (+- 3 S.D from the mean) for each cell type.

- All macrophages (M1+M2 macrophages).
- Adipocytes vs BMI.
- You can also test the other cell types present in the signature matrix if you observe they are present in your samples.

If you have any other obesity metric, such as central adiposity (We have used DEXA scans for this) we would also be grateful to see a correlation with cell type proportion.

