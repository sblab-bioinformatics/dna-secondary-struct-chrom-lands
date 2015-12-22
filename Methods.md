# Data analysis

For ease of presentation, paths to programs and to data files are omitted in the following scripts.
Input and output file names are often replaced by placeholder variables (strings starting with `$` _as per_ bash syntax).


## Software tools and public data files

Data processing was performed in Linux environment with GNU coreutils tools.
The following software and data files are required for the analysis described here:

* [cutadapt](http://cutadapt.readthedocs.org/en/stable/guide.html) version 1.8

* [bwa](https://github.com/lh3/bwa) version 0.7

* [tophat2](http://ccb.jhu.edu/software/tophat/index.shtml) version 2.0

* [samtools](http://www.htslib.org/) version 1.1

* [Picard](https://broadinstitute.github.io/picard/command-line-overview.html) version 1.140

* [bedtools](http://bedtools.readthedocs.org/en/latest/) version 2.24

* Human reference genome version hg19 in fasta format, indexes for bwa and bowtie2, and gene annotation file (`genes.gtf`)
were downloaded from [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)

* [F-Seq](https://github.com/aboyle/F-seq) version 1.84

* [macs2](https://github.com/taoliu/MACS/) version 2.1.0.20150731

* [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) version 0.6

* [R](https://cran.r-project.org/) with packages:
[edgeR](https://www.bioconductor.org/packages/3.3/bioc/html/edgeR.html),
[data.table](https://cran.r-project.org/web/packages/data.table/index.html),
[ggplot2](http://ggplot2.org/),

* Convenience custom scripts: [mergePeaks.sh](https://github.com/dariober/bioinformatics-cafe/blob/master/mergePeaks.sh),
[sortBedAsBam.py](https://github.com/dariober/bioinformatics-cafe/blob/master/sortBedAsBam.py), [tableCat.py](https://github.com/dariober/bioinformatics-cafe/tree/master/tableCat)


* List of regions excluded from peak calling:
[hg19.wgEncodeDukeMapabilityRegionsExcludable.bed.gz](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz).

The intervals in this file were subtracted from the reference genome to obtain a file of regions to include (`hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed`):

``` 
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -A -e \
    "select chrom, size from hg19.chromInfo" \
| awk -v OFS="\t" '{print $1, 0, $2}' \
| subtractBed -a - -b wgEncodeDukeMapabilityRegionsExcludable.bed.gz \
| sort -k1,1 -k2,2n > hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed
```

## Processing FAIRE-Seq data

Raw fastq files were trimmed to remove adapter contamination and were aligned to the reference genome.
Aligned files were filtered to remove secondary and not primary alignments as well as reads with mapping quality below 10:

```
cutadapt -f fastq -e 0.1 -q 20 -O 3 -a CTGTCTCTTATACACATCT $fq > /dev/stdout 2> $bname.cutadapt.txt \
| bwa mem -M genome.fa /dev/stdin \
| samtools view -S -u -F2304 -q 10 -L hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed - \
| samtools sort - ${bname}.hg19.tmp &&
java -Xmx5g -jar ~/bin/picard.jar MarkDuplicates I=${bname}.hg19.bam O=$bamclean/${bname}.hg19.clean.bam M=$logdir/$bname.md.txt
```

Technical replicates were then merged in a single file. The header in the merged bam files edited with [addRGtoSAMHeader.py](https://github.com/dariober/bioinformatics-cafe/blob/f524a45ea0d85b9be7cf24508b32c9488b74ae95/addRGtoSAMHeader.py)

```
java -Xmx1g -jar picard.jar MergeSamFiles VALIDATION_STRINGENCY=SILENT AS=true I=$toMerge1 I=$toMerge2 O=$bam &&
addRGtoSAMHeader.py -i $bam -o $outbam &&
samtools index ${bam%.bam}.rg.bam
```

### Mapping FAIRE-Seq peaks

Before peak calling, alignment bam files were first split by chromosome and the alignments mapped as duplicate removed. Finally, bam files were
converted to bed format. Peak calling was perfomed with F-Seq as follows for each FAIRE-Seq bam file:

```
bff='crg_100bp_hg19' ## Directory obtained from http://fureylab.web.unc.edu/software/fseq/
fseq -of bed -t 5 -l 800 -v -b $bff -d $bdir -o $outputdir
```

Peak files from individual chromosomes were then concateated.

## Processing ATAC-Seq data

Raw reads were trimmed, aligned, filtered and duplicates marked using the same procedure as for FAIRE-Seq. 
Prior to mapping of open chromatin, reads mapping to chrM were removed and technical replicates merged in a single bam file.

Peaks of read enrichment were mapped using macs2 as follows:

```
samtools view -h -F1024 $bam | grep -v -P '\tchrM\t' | samtools view -b - > $tmpBam
macs2 callpeak --keep-dup all -t $tmpBam -n ${bname} 
```

### Differential chromatin state between entinostat and control cells

Differences in ATAC-Seq signal and hence chromatin state between entinostat treated and control cells were mapped by comparing the read counts in the two sets
of libraries. First, consensus peaks within entinostat and control libraries were defined as those shared between at least two libraries:

```
## Consensus peaks in entinostat
mergePeaks.sh rhh_ATAC_entinostat_bio1_23102015_peaks.narrowPeak.gz \
    rhh_ATAC_entinostat_bio2_23102015_peaks.narrowPeak.gz | awk '$5 > 1' > entino.merge.narrowPeak

## Consensus peaks control
mergePeaks.sh rhh145-146_K3_K4_atac_hacat_peaks.narrowPeak \
    rhh147-148_untreat_14102014_atac_hacat_peaks.narrowPeak \
    rhh149-150_untreat_27052014_atact_hacat_peaks.narrowPeak | awk '$5 > 1' > ctrl.merge.narrowPeak

```

The union of the two consensus peak sets was tested for differential chromatin state:

```
mergePeaks.sh ctrl.merge.narrowPeak entino.merge.narrowPeak \
| sortBedAsBam.py -i - -b /nas/sblab_data1/berald01/repository/bam_clean/rhh_ATAC_entinostat_bio1_23102015.bam > union.merge.bed

# Count reads in peaks
for bam in rhh_ATAC_entinostat_bio1_23102015.bam \
           rhh_ATAC_entinostat_bio2_23102015.bam \
           rhh145-146_K3_K4_atac_hacat.bam \
           rhh147-148_untreat_14102014_atac_hacat.bam \
           rhh149-150_untreat_27052014_atact_hacat.bam
do
coverageBed -g genome.txt -counts -sorted -a union.merge.bed -b $bam > ${bam%%.bam}.union.bed
done

echo 'chrom start end peakFiles nfiles nreads sample_id' | tr ' ' '\t' > union.merge.counts.bed
tableCat.py -i rhh*.union.bed -r '.union.bed' >> union.merge.counts.bed
rm rhh*union.bed
rm union.merge.bed
```

`genome.txt` is tab delimited file giving the of chromosomes in the input files. See the documentation of `coverageBed` for further details.

Read counts were adjusted by libraries size defined as number of reads mapped in all chromosomes but excluding chrM:

```
for bam in rhh_ATAC_entinostat_bio1_23102015.bam \
           rhh_ATAC_entinostat_bio2_23102015.bam \
           rhh145-146_K3_K4_atac_hacat.bam \
           rhh147-148_untreat_14102014_atac_hacat.bam \
           rhh149-150_untreat_27052014_atact_hacat.bam
do
samtools idxstats $bam | awk -v bam=$bam '$1 != "chrM" {sum+=$3}END{print bam, sum}'
done
```

Testing for differential chromatin state:

```
R
library(data.table)
library(edgeR)
library(reshape2)
cnt<- fread('union.merge.counts.bed')
cnt[, locus := paste(chrom, start, end, sep= '_')]
cntct<- dcast.data.table(dat= cnt, locus ~ sample_id, value.var= 'nreads') 
y<- data.frame(cntct[, 2:ncol(cntct), with= FALSE])
row.names(y)<- cntct$locus

# Lib size *without* chrM from above:
libSize<- c(
    rhh145.146_K3_K4_atac_hacat= 81791133,
    rhh147.148_untreat_14102014_atac_hacat= 65675634,
    rhh149.150_untreat_27052014_atact_hacat= 141758614,
    rhh_ATAC_entinostat_bio1_23102015= 17546666,
    rhh_ATAC_entinostat_bio2_23102015= 20667324)

group <- factor(c('ctrl', 'ctrl', 'ctrl', 'ent', 'ent'))

y<- DGEList(counts=y, group=group)
stopifnot(rownames(y$samples) == names(libSize))
y$samples$lib.size<- libSize
y<- calcNormFactors(y, method= 'none')
y<- estimateDisp(y)
y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)

et<- exactTest(y, pair= levels(y$samples$group))

detable<- data.frame(topTags(et, n= Inf)$table)
detable$locus<- rownames(detable)
detable<- data.table(detable)
detable<- merge(detable, unique(cnt[, list(chrom, start, end, locus)]), by= 'locus')

pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
pdf('maplot.atac-entinostat.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
    main= "Differential chromatin (ATAC)\n[entinostat - ctrl]", colramp= pal, col= 'blue')
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= 0, col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.table(detable, "diff.atac-entinostat.txt", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```

### Differential chromatin state between HaCaT and NHEK cells

Differences in ATAC signal between HaCaT and NHEK cells were detected following the same approach as above for entinostat vs control.
First, consensus peak sets were produced for the two cell lines, then a union peak set was generated for read count and testing:

```
## Consensus peaks HEK
mergePeaks.sh rhh_HEKnp_ATAC_24022015_peaks.narrowPeak.gz \
    rhh_HEKnp_ATAC_27032015_peaks.narrowPeak.gz | awk '$5 > 1' > hek_atac.merge.narrowPeak

## Consensus peaks HACAT
mergePeaks.sh rhh145-146_K3_K4_atac_hacat_peaks.narrowPeak \
    rhh147-148_untreat_14102014_atac_hacat_peaks.narrowPeak \
    rhh149-150_untreat_27052014_atact_hacat_peaks.narrowPeak | awk '$5 > 1' > hacat_atac.merge.narrowPeak

## Union set
mergePeaks.sh hacat_atac.merge.narrowPeak  hek_atac.merge.narrowPeak \
| sortBedAsBam.py -i - -b /nas/sblab_data1/berald01/repository/bam_clean/rhh145-146_K3_K4_atac_hacat.bam > union_atac.merge.bed

# Count reads in peaks
for bam in rhh_HEKnp_ATAC_24022015.mrg.bam \
           rhh_HEKnp_ATAC_27032015.mrg.bam \
           rhh145-146_K3_K4_atac_hacat.bam \
           rhh147-148_untreat_14102014_atac_hacat.bam \
           rhh149-150_untreat_27052014_atact_hacat.bam
do
    coverageBed -g genome.txt -counts -sorted -a union_atac.merge.bed -b $bam > ${bam%%.bam}.union.bed
done

echo 'chrom start end peakFiles nfiles nreads sample_id' | tr ' ' '\t' > union_atac.merge.counts.bed
tableCat.py -i rhh*.union.bed -r '.union.bed' >> union_atac.merge.counts.bed
rm rhh*union.bed
rm union_atac.merge.bed
```

Differential chromatin state:

```
R
library(data.table)
library(edgeR)
library(reshape2)
cnt<- fread('union_atac.merge.counts.bed')
cnt[, locus := paste(chrom, start, end, sep= '_')]
cntct<- dcast.data.table(dat= cnt, locus ~ sample_id, value.var= 'nreads') 
y<- data.frame(cntct[, 2:ncol(cntct), with= FALSE])
row.names(y)<- cntct$locus

# Lib size *without* chrM
libSize<- c(
    rhh145.146_K3_K4_atac_hacat=             81791133,
    rhh147.148_untreat_14102014_atac_hacat=  65675634,
    rhh149.150_untreat_27052014_atact_hacat= 141758614,
    rhh_HEKnp_ATAC_24022015.mrg=             153476303,
    rhh_HEKnp_ATAC_27032015.mrg=             91762466 
)

group <- factor(c('hacat', 'hacat', 'hacat', 'hek', 'hek'))

y<- DGEList(counts=y, group=group)
stopifnot(rownames(y$samples) == names(libSize))
y$samples$lib.size<- libSize
y<- calcNormFactors(y, method= 'none')
y<- estimateDisp(y)
y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)

plotMDS(y, main= 'MDS for differential chromatin state (ATAC-Seq)')

et<- exactTest(y, pair= levels(y$samples$group))

detable<- data.frame(topTags(et, n= Inf)$table)
detable$locus<- rownames(detable)
detable<- data.table(detable)
detable<- merge(detable, unique(cnt[, list(chrom, start, end, locus)]), by= 'locus')

pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
pdf('maplot.atac-hek_vs_hacat.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
    main= "Differential chromatin (ATAC)\n[HEK - HACAT]", colramp= pal, col= 'blue')
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= 0, col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.table(detable, "diff.atac-hek_vs_hacat.txt", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```

## Processing G4-ChIP data

Raw reads were trimmed, aligned, filtered and duplicates marked using the same procedure as for FAIRE-Seq. 
Prior to mapping of open chromatin technical replicates were merged in a single bam file.

Peaks were identified with macs2 with the appropriate input control for each pull-down library:

```
macs2 callpeak --keep-dup all \
    -t $bamdir/rhh_25cyc_BG4_12082015.bam \
    -c $bamdir/rhh155_25cyc_input_703_503_12082015.bam -n rhh_25cyc_BG4_12082015

macs2 callpeak --keep-dup all \
    -t $bamdir/rhh175_ChIPwthacat_704_502_entst_26082015.bam \
    -c $bamdir/rhh178_inputwthacat_703_517_entst_26082015 .bam -n rhh175_ChIPwthacat_704_502_entst_26082015

for bam in HEKnp_Lonza_1472015_BG4.md.bam HEKnp_Lonza_1572015_BG4.md.bam
do
macs2 callpeak --keep-dup all -p 0.0001 \
    -t $bam \
    -c merged_14_and_15072015_input_heknplonza.md.bam -n ${bam%.md.bam}.1e4
done
```

### Differential BG4 binding between entinostat and control cells

Differences in BG4 binding between entinostat treated and control cells have been tested for the union of intervals between the four
G4-ChIP libraries. The union set of intervals and read count in each interval has been produced as follows:

```
mergePeaks.sh rhh_25cyc_BG4_12082015_peaks.narrowPeak \
    rhh175_ChIPwthacat_704_502_entst_26082015_peaks.narrowPeak \
    rhh_ChIP_entst_17082015_peaks.narrowPeak \
    rhh_ChIP_entst_26082015_peaks.narrowPeak \
| sortBedAsBam.py -i - -b rhh_ChIP_entst_17082015.bam > union.bed

for bam in rhh_25cyc_BG4_12082015.bam \
           rhh175_ChIPwthacat_704_502_entst_26082015.bam \
           rhh_ChIP_entst_17082015.bam \
           rhh_ChIP_entst_26082015.bam
do
    coverageBed -sorted -a union.bed -b $bam > ${bam%%.bam}.union.bed
done

echo 'chrom start end peakFiles nfiles nreads nonzero len frac sample_id' | tr ' ' '\t' > union.counts.bed
tableCat.py -i rhh*.union.bed -r '.union.bed' >> union.counts.bed
rm rhh*union.bed
rm union.bed
```

The resulting matrix of counts having G4-ChIP sites as rows and libraries as columns was analyzed for differential binding.
Note that the number of reads in each alignment file was used as library size.

```
R
library(data.table)
library(edgeR)
library(reshape2)
cnt<- fread('union.counts.bed')
cnt[, locus := paste(chrom, start, end, sep= '_')]
cntct<- dcast.data.table(dat= cnt, locus ~ sample_id, value.var= 'nreads') 
y<- data.frame(cntct[, 2:ncol(cntct), with= FALSE])
row.names(y)<- cntct$locus
group <- factor(c('ctrl', 'ctrl', 'ent', 'ent'))

## Library sizes as n. reads in each bam file:
libSize<- c(rhh175_ChIPwthacat_704_502_entst_26082015= 243965455,
    rhh_25cyc_BG4_12082015= 135996072,
    rhh_ChIP_entst_17082015= 99481706,
    rhh_ChIP_entst_26082015= 152700155)

y<- DGEList(counts=y, group=group)
stopifnot(rownames(y$samples) == names(libSize))
y$samples$lib.size<- libSize
y<- calcNormFactors(y, method= 'none')
y<- estimateDisp(y)

y<- estimateDisp(y)
y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)
et<- exactTest(y, pair= levels(y$samples$group))
detable<- data.frame(topTags(et, n= Inf)$table)
detable$locus<- rownames(detable)
detable<- data.table(detable)
detable<- merge(detable, unique(cnt[, list(chrom, start, end, locus)]), by= 'locus')

pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
pdf('maplot.entinostat.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC', main= "BG4 differential binding [entinostat - ctrl]", colramp= pal, col= 'blue')
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= 0, col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.table(detable, "expr_diff.entinostat.txt", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```

### Differential BG4 binding between HaCaT and NHEK cells

Union set of testable site was produced merging the four libraries from HaCaT and NHEK cells:

```
mergePeaks.sh \
    rhh_25cyc_BG4_12082015_peaks.narrowPeak \
    rhh175_ChIPwthacat_704_502_entst_26082015_peaks.narrowPeak \
    HEKnp_Lonza_1472015_BG4.1e4_peaks.narrowPeak \
    HEKnp_Lonza_1572015_BG4.1e4_peaks.narrowPeak \
| sortBedAsBam.py -i - -b rhh_25cyc_BG4_12082015.bam > union_bg4.merge.bed

# Count reads in union peak set
for bam in rhh_25cyc_BG4_12082015.bam \
           rhh175_ChIPwthacat_704_502_entst_26082015.bam \
           HEKnp_Lonza_1472015_BG4.md.bam \
           HEKnp_Lonza_1572015_BG4.md.bam
do
    coverageBed -g genome.txt -sorted -a union_bg4.merge.bed -b $bam > ${bam%%.bam}.union.bed
done

echo 'chrom start end peakFiles nfiles nreads nonzero len frac sample_id' | tr ' ' '\t' > union_bg4.merge.counts.bed
tableCat.py -i rhh*.union.bed HEKnp*.union.bed -r '\..*' >> union_bg4.merge.counts.bed
rm rhh*union.bed
rm union_bg4.merge.bed
```

Testing for differential binding:

```
R
library(data.table)
library(edgeR)
library(reshape2)
cnt<- fread('union_bg4.merge.counts.bed')
cnt[, locus := paste(chrom, start, end, sep= '_')]
cntct<- dcast.data.table(dat= cnt, locus ~ sample_id, value.var= 'nreads') 
y<- data.frame(cntct[, 2:ncol(cntct), with= FALSE])
row.names(y)<- cntct$locus
group <- factor(c('hek', 'hek', 'hacat', 'hacat'))

## Libray sizes:
libSize<- c(
    HEKnp_Lonza_1472015_BG4= 260449974,
    HEKnp_Lonza_1572015_BG4= 261673243,
    rhh175_ChIPwthacat_704_502_entst_26082015 = 243965455,
    rhh_25cyc_BG4_12082015= 135996072
)

y<- DGEList(counts=y, group=group)
stopifnot(rownames(y$samples) == names(libSize))
y$samples$lib.size<- libSize
y<- calcNormFactors(y, method= 'none')
y<- estimateDisp(y)

y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)
et<- exactTest(y, pair= levels(y$samples$group))
detable<- data.frame(topTags(et, n= Inf)$table)
detable$locus<- rownames(detable)
detable<- data.table(detable)
detable<- merge(detable, unique(cnt[, list(chrom, start, end, locus)]), by= 'locus')

pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
pdf('maplot.BG4-hek_vs_hacat.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
    main= "BG4 differential binding [HEK - HACAT]", colramp= pal, col= 'blue')
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= 0, col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.table(detable, "diff.BG4-hek_vs_hacat.txt", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```

## Processing RNA-Seq data

Raw reads were trimmed, aligned to the reference genome and assigned to genes as follows:

```
cutadapt -O 3 -a AGATCGGAAGAGC -o cutadapt/$fq $fq" 
tophat2 -o ${outdir}/${fq%%.*}_hg19 -p 8 --library-type fr-unstranded -G genes.gtf genome ${fq}

# Count reads in genes
samtools view $bam | htseq-count -m intersection-strict -s no -t exon -i gene_id - genes.gtf > ${htseqcount}"
```

A matrix of counts of reads in genes, `all.htseq`, was generated by binding side by side the counts from the individual files. The first row of `all.htseq`
is header, gene names are in the first column.

### Differential gene expression: Effect of entinostat in hacat cells

```
R
library(data.table)
library(edgeR)
library(reshape2)
cnt<- read.table('all.htseq', row.names= 1, header= TRUE)

## Entinostat treated libraries are Entinostat_[1-4], controls are HaCaT_[1-4]
cnt<- cnt[!rownames(cnt) %in% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"),
    c("Entinostat_1", "Entinostat_2", "Entinostat_3", "Entinostat_4", "HaCaT_1", "HaCaT_2", "HaCaT_3", "HaCaT_4")]
keep<- rowSums(cpm(cnt) > 1) >= 2
cnt<- cnt[keep, ]
group <- factor(c('ent', 'ent', 'ent', 'ent', 'ctrl', 'ctrl', 'ctrl', 'ctrl'))

y<- DGEList(counts= cnt, group= group)
y<- calcNormFactors(y)
y<- estimateDisp(y)
y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)

et<- exactTest(y, pair= levels(y$samples$group))
detable<- data.frame(topTags(et, n= Inf)$table)
detable$gene_id<- rownames(detable)
detable<- data.table(detable)

## MA-plot of gene expression
## --------------------------
pal<- colorRampPalette("transparent", space = "Lab") # Do not colour NS genes
pdf('maplot.entinostat.rnaseq.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
    main= "Differential gene [entinostat - ctrl]", colramp= pal, col= 'blue', nrpoints= 0)
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= c(-2, 0, 2), col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.') # Mark DE genes
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.table(detable, "expr_diff.entinostat.rnaseq.txt", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```

The following R function was used to normalize raw read counts in genes to transcripts per million (TPM):

```
tpm<- function(counts, gene_length){
    stopifnot(length(counts) == length(gene_length))
    stopifnot(is.numeric(counts))
    stopifnot(gene_length > 0)
    rate<- counts / gene_length
    denom<- sum(rate)          
    xtpm<- rate / denom * 1e6  
    return(xtpm)
}
```
Gene lengths were extracted from annotation GTF file `genes.gtf` using custom script [geneLengthFromGTF.py](https://github.com/dariober/bioinformatics-cafe/blob/master/geneLengthFromGTF.py).

### Differential gene expression between HaCaT and NHEK

Differential gene expression between cell lines was tested with the same procedure as for the difference between entinostat and control in HaCaT:

```
R
library(data.table)
library(edgeR)
library(reshape2)
cnt<- read.table('all.htseq', row.names= 1, header= TRUE)
cnt<- cnt[!rownames(cnt) %in% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"), 
    c("HEK_Gibco_1", "HEK_Gibco_2", "HEK_Gibco_3", "HEK_Gibco_4", "HaCaT_1", "HaCaT_2", "HaCaT_3", "HaCaT_4")]
keep<- rowSums(cpm(cnt) > 1) >= 2
cnt<- cnt[keep, ]
group <- factor(c('hek', 'hek', 'hek', 'hek', 'hacat', 'hacat', 'hacat', 'hacat'))

y<- DGEList(counts= cnt, group= group)
y<- calcNormFactors(y)
y<- estimateDisp(y)
y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)

et<- exactTest(y, pair= levels(y$samples$group))
detable<- data.frame(topTags(et, n= Inf)$table)
detable$gene_id<- rownames(detable)
detable<- data.table(detable)

pal<- colorRampPalette("transparent", space = "Lab")
pdf('maplot.hek_vs_hacat.rnaseq.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
    main= "Differential gene [HEK - HACAT]", colramp= pal, col= 'blue', nrpoints= 0)
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= c(-2, 0, 2), col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')
dev.off()

write.table(detable, "expr_diff.hek_vs_hacat.rnaseq.txt", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
```

## Gene expression in relation to ATAC and G4-ChIP peaks in promoters

Gene promoters were annotated with the presence ATAC-Seq and G4-ChIP peaks in order to assess the effect of ATAC-Seq and G4-ChIP
peaks on gene expression.

Promoters were defined as the regions spanning the transcription start sites by 1000 bp up- and down-stream. From the gene annotation
file `genes.gtf`, promoters were extracted as follows:


```
# Check column 9 is always gene_id and column 11 transcript_id
nrec=`wc -l genes.gtf | cut -d ' ' -f 1`
ngene=`awk '$9 == "gene_id"' genes.gtf | wc -l | cut -d ' ' -f 1`
ntx=`awk '$11 == "transcript_id"' genes.gtf | wc -l | cut -d ' ' -f 1`
echo $nrec $ngene $ntx

#  + Strand: Get start of first exons on each transcript
awk -v OFS="\t" '$7 == "+" {print $1, $4, $5, $7, $10, $12}' genes.gtf \
| sed 's/"//g' | sed 's/ //g' | sed 's/;//g' \
| sort -k 6,6 -k2,2n \
| groupBy -g 6 -c 1,2,5 -o first,first,first \
| awk -v OFS="\t" '{print $2, $3-1000, $3+1000, $4, $1, "+"}' > tss.plus.bed

# - strand: Get end of last exons on each transcript
awk -v OFS="\t" '$7 == "-" {print $1, $4, $5, $7, $10, $12}' genes.gtf \
| sed 's/"//g' | sed 's/ //g' | sed 's/;//g' \
| sort -k 6,6 -k3,3nr \
| groupBy -g 6 -c 1,3,5 -o first,first,first \
| awk -v OFS="\t" '{print $2, $3-1000, $3+1000, $4, $1, "-"}' > tss.minus.bed

# Merge promoters within genes
cat tss.minus.bed tss.plus.bed \
| cut -f1,2,3,4 \
| sort \
| uniq \
| awk -v OFS="\t" '{print $1 "_" $4, $2, $3, $4}' \
| sortBed \
| mergeBed \
| sed 's/_/\t/' \
| awk -v OFS="\t" '{print $1, $3, $4, $2}' \
| sortBed > hg19.gene_name.promoters.bed

## Get only genes with one merged promoter
R
library(data.table)
bed<- fread('hg19.gene_name.promoters.bed')
setnames(bed, names(bed), c('chrom', 'start', 'end', 'gene_name'))
cntgene<- bed[, .N, by= gene_name]
bed<- bed[gene_name %in% cntgene[N == 1, gene_name] ]
write.table(bed, 'hg19.gene_name.promoters.bed', sep= '\t', row.names= FALSE, col.names= FALSE, quote= FALSE)
quit(save= 'no')
##
rm tss.minus.bed tss.plus.bed
```

The obtained promoters were then annotated with ATAC and G4-ChIP peaks and observed quadruplex sequences. 

```
## ATAC sites
mergePeaks.sh \
   rhh145-146_K3_K4_atac_hacat_peaks.narrowPeak \
   rhh147-148_untreat_14102014_atac_hacat_peaks.narrowPeak \
   rhh149-150_untreat_27052014_atact_hacat_peaks.narrowPeak \
| awk '$5 >= 2' > atac_hacat.narrowPeak 

annotateBed -counts -i hg19.gene_name.promoters.bed -files \
    atac_hacat.narrowPeak \
    ../OQs/Na_PDS_hits_intersect.bed.gz \
    BG4_hacat.1rep.narrowPeak \
    BG4_hacat.2rep.narrowPeak \
| sortBed > hg19.promoters.ant.bed
rm hg19.gene_name.promoters.bed
```

## Correlation in read density between replicates

The consistency between replicates was assessed by comparing the tag density in windows of 1000 bp.

```
# Analyze only chr19:
grep 'chr19' genome.bed \
| windowMaker -b - -w 1000 > chr19.windows.bed 

for bam in  rhh175_ChIPwthacat_704_502_entst_26082015.bam \
            rhh_25cyc_BG4_12082015.bam \
            rhh145-146_K3_K4_atac_hacat.bam \
            rhh147-148_untreat_14102014_atac_hacat.bam \
            rhh149-150_untreat_27052014_atact_hacat.bam \
            rhh_hacat_05082014_FAIRE.bam \
            rhh_hacat_14062014_FAIRE.bam \
            rhh_ChIP_entst_17082015.bam \
            rhh_ChIP_entst_26082015.bam \
            rhh_ATAC_entinostat_bio1_23102015.bam \
            rhh_ATAC_entinostat_bio2_23102015.bam \
            HEKnp_Lonza_1472015_BG4.md.bam \
            HEKnp_Lonza_1572015_BG4.md.bam \
            rhh_HEKnp_ATAC_24022015.mrg.bam \
            rhh_HEKnp_ATAC_27032015.mrg.bam \
            rhh_hek_09112014_FAIRE.bam \
            rhh_hek_18092014_FAIRE.bam
do
echo "coverageBed -g genome.txt -sorted \
    -a chr19.windows.bed -b $bamdir/$bam -counts > ${bam%%.bam}.cnt.bed" > ${bam%%.bam}.tmp.sh
done
ls *.tmp.sh | xargs -P 0 -n 1 bash && rm *tmp.sh

tableCat.py -i *.cnt.bed -r '\.cnt\.bed' > windows.cnt.cat.bed

R
library(data.table)
library(ggplot2)
reps<- fread('windows.cnt.cat.bed')
design<- fread('design.txt')
setnames(reps, names(reps), c('chrom', 'start', 'end', 'cnt', 'library_id'))
design[, library_id := sub('.bam', '', file_name)]

stopifnot(design$library_id == unique(design$library_id))
stopifnot(unique(reps$library_id) %in% design$library_id)

reps<- merge(reps, design, by= 'library_id')
nrm<- reps[, list(chrom, start, end, cnt_nrm= (cnt+1) / sum(.SD$cnt+1)), by= list(library_id, cell_line, chip, rept)]
nrm[, list(sum= sum(cnt_nrm)), by= library_id]$sum

rep13<- nrm[cell_line == 'HACAT' & chip == 'ATAC' & rept == 1, list(cell_line, chip, chrom, start, end, rep_1= cnt_nrm)]
rep13$rep_2<- nrm[cell_line == 'HACAT' & chip == 'ATAC' & rept == 3, cnt_nrm]
rep13[, cmp := 'rep 1 vs 3']

rep23<- nrm[cell_line == 'HACAT' & chip == 'ATAC' & rept == 2, list(cell_line, chip, chrom, start, end, rep_1= cnt_nrm)]
rep23$rep_2<- nrm[cell_line == 'HACAT' & chip == 'ATAC' & rept == 3, cnt_nrm]
rep23[, cmp := 'rep 2 vs 3']

nrm<- dcast.data.table(data= nrm[rept %in% c(1, 2)], cell_line + chip + chrom + start + end ~ rept, value.var= 'cnt_nrm')
setnames(nrm, c('1', '2'), c('rep_1', 'rep_2'))
nrm[, cmp := 'rep 1 vs 2']

nrm<- rbindlist(list(nrm, rep13, rep23))
nrm[, xtitle := paste(cell_line, chip, cmp)]
gg<- ggplot(nrm[seq(1, nrow(nrm), length.out= 100000)], aes(x= log10(rep_1), y= log10(rep_2))) +
    geom_point(size= 0.5, alpha= 0.20) +
    geom_abline(yintercept= 0,  slope= 1, colour= 'red', linetype= 'dashed') +
    facet_wrap(~xtitle) +
    xlab('log10(count)') + 
    ylab('log10(count)') +
    ggtitle('Correlation in normalized read count in 1 kb bins between replicates') +
    geom_text(data= nrm[, list(r= cor(rep_1, rep_2)), by= xtitle],
        aes(x= -6, y= -2, label= sprintf('r= %s', round(r, 2))), size= 5)
ggsave('bam_correl.pdf', w= 24, h= 24, units= 'cm')
```

## Visualization in IGV

For the purpose of visualization onn the Integrative Genome Browser ([IGV](https://www.broadinstitute.org/igv/)), bam files were converted to TDF format:

```
igvtools count $bam $tdf $hg19ChromSizes
```
