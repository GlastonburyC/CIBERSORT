# Using CIBERSORT to estimate cell type composition.

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
java -jar CIBERSORT.jar -M TwinsUK.adiposeSamples.txt -B Adipose.sigMatrix.txt –n 1000 >> TwinsUK.Adipose.CellEsts.txt
```
