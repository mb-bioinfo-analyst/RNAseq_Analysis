---
title: RNA-Seq Analysis Workflow
type: docs
---

## Prepare Data

### Download read data

For this analysis, Mouse mammary data (fastq files) can be downloaded from [figshare](https://figshare.com/s/f5d63d8c265a05618137) or [ncbi (GSE60450)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450)

I have used the package Rsubread in this analysis. You can find more information about this package [here](https://bioconductor.org/packages/release/bioc/html/Rsubread.html).

The raw reads downloaded from SRA are in .sra format. The sra toolkit from NCBI is used to convert the .sra files to .fastq files using the fastq-dump command. You can find more details on sra-toolkit [here](https://hpc.nih.gov/apps/sratoolkit.html). 

### Downloading genome files

For this analysis, the index files for chromosome 1 of mouse genome build mm10 is provided in the data files from [figshare](https://figshare.com/s/f5d63d8c265a05618137). However, you can obtain full genome fasta files for various genomes from the [UCSC genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu) or [here](https://hgdownload.soe.ucsc.edu/downloads.html); or from [NCBI](https://www.ncbi.nlm.nih.gov/genome); or from [ENSEMBL]
(https://www.ensembl.org/info/data/ftp/index.html).

### Import Data for Alignment

Mapping reads to the genome is a very important task, and many different aligners are available, such as STAR [@Dobin2013], topHat [@trapnell2009tophat], bowtie [@Langmead2012] and Rsubread [@liao2013subread]. **Rsubread** is the only aligner that can run in R. The majority of alignment programs run on Linux and need a lot of processing power. The majority of mapping activities demand more powerful processors than a typical laptop, hence read mapping is often carried out on a server that has linux-like environment. We will only map 1,000 reads from each sample of our mouse lactation dataset from [@Fu2015](https://www.nature.com/articles/ncb3117), but because of time and limitation of computational power, I will only be mapping to chromosome 1.

First, let's load the Rsubread package into R.

```r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")
library(Rsubread)

```

Make sure that all of the sequencing read data (.fastq.gz files) ar placed in in the data directory. `list.files` command can be used to search for all .fastq.gz files in the data directory.
In the code below `$` symbol represents the end of the string, so that we will get the files that end in ".fastq.gz" in their names.

```r

fastq.files <- list.files(path = "./data", pattern = ".fastq.gz$", full.names = TRUE)
# fastq.files

```

## Alignment

### Build the index

Before performing alignment, an index of the reference genome is to be build. To build an index you need need the fasta file (.fa), which can be downloaded from the links given above. For the reasons mentioned earlier, here I am building the index only for chromosome 1. It usually takes around 30 minutes to 1 hour for the whole mouse genome index if run on a dedicated server. The command below builds the index of the chr1 genome for mm10 using the "chr1.fa" file.

```r

buildindex(basename="chr1_mm10",reference="chr1.fa")

```

Once the index is generated, you can see the additional files generated in the same directory. 

```r

dir()

```

### Aligning reads to chromosome 1 of reference genome

Now we can align our reads using the `align` command. We can tweak with the parameters in accordance to our data, but the default parameters for the `align` function are usually fine. This input data is comprised of 100 base pairs single end reads, but if we have paired end data, we need to specify the second read file/s using the `readfile2` argument.

`Rsubread` sets the output file names by default by using the filename with addition of ".subread.BAM" at the end.

we can align our 12 fastq.gz files using the `align` command.


```r

align(index="data/chr1_mm10",readfile1=fastq.files)

```

The 12 samples will be aligned to the reference genome one after the other and The BAM files are saved in the working directory.

The default setting for `align` allows us to only keeps reads that uniquely map to the reference genome. **Differential expression of genes, this is important**, because the reads are assigned only to one place in the genome, which is easier to interpret. But if the circumstances are that your data needs different parameteric settings, you can get more information from looking at the align documentation.

```r

?align

```

To get a summary of the proportion of reads that are mapped to the reference genome, `propmapped` function can be used.

```r

# enlist the bam files
bam.files <- list.files(path = "./data", pattern = ".BAM$", full.names = TRUE)
# bam.files

# summarize the proportion of mapped reads
props <- propmapped(files=bam.files)
props

```


## Quality control

Let's extract quality scores for the file "SRR1552450.fastq.gz" using 100 reads. A quality score of 30 corresponds to a 1 in 1000 chance of an incorrect base call. (A quality score of 10 is a 1 in 10 chance of an incorrect base call.) The overall distribution of quality scores can be seen using a boxplot.

```r

# Extract quality scores
qs <- qualityScores(filename="data/SRR1552450.fastq.gz",nreads=100)
# Check dimension of qs
dim(qs)
# Check first few elements of qs with head
head(qs)

boxplot(qs)

```


## Read Counts

The alignment step provides us with a set of BAM files, where every file contains read alignments for each library. BAM file contains a chromosome location for all of the uniquely mapped reads. These mapped reads can be accounted for mouse genes with the help of `featureCounts` function. `featureCounts` provides a  built-in annotation for mouse and human genome assemblies (NCBI refseq annotation).

Number of reads that are map to exons of genes are summed to obtain the count for each gene. `featureCounts` takes a list of BAM files as input, and returns an object that includes the count matrix. This matrix is made up of samples in columns and genes in rows.

```r

fc <- featureCounts(bam.files, annot.inbuilt="mm10")

```

fc$stats reports the numbers and reasons of unassigned reads (eg. ambiguity, multi-mapping, secondary alignment, mapping quality, fragment length, chimera, read duplicate, non-junction and so on), along with the number of successful assignment of reads for each library. See [subread documentation](http://bioconductor.org/packages/release/bioc/html/Rsubread.html)

```r

## Take a look at the featurecounts stats
fc$stat

```

The counts for the samples are stored in fc$counts. Take a look at that.

```r

## Take a look at the dimensions to see the number of genes
dim(fc$counts)

# head(fc$counts)

```

Here, row names of the matrix represent the Entrez gene ids and the column names are the samples. The `annotation` slot shows the annotation information that `featureCounts` used to summarize reads over genes.

```r

head(fc$annotation)

```

## Tips and Tricks

{{< columns >}}

The Phred offset determines the encoding for the base-calling quality string in the
FASTQ file. For the Illumina 1.8 format onwards, this encoding is set at +33.
However, older formats may use a +64 encoding. Users should ensure that the
correct encoding is specified during alignment. If unsure, one can examine the
first several quality strings in the FASTQ file. A good rule of thumb is to check
whether lower-case letters are present (+64 encoding) or absent (+33).

<--->


featureCounts` requires gene annotation specifying the genomic start and end
position of each exon of each gene. *Rsubread* contains built-in gene annotation
for mouse and human. For other species, users will need to read in a data frame
in GTF format to define the genes and exons. Users can also specify a custom annotation file in SAF format. 
{{< /columns >}}

