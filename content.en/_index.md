---
title: "RNA-Seq Analysis Workflow"
type: docs
date: "June 22, 2023"
output: pdf_document
---


# RNA-Seq Analysis: A Beginner's Guide to Understanding Gene Expression Using R
**Author: Bilal Mustafa PhD**

**Date: June 22, 2023**

***  

![](https://upload.wikimedia.org/wikipedia/commons/6/68/MIcroarray_vs_RNA-Seq.png)



***

## Introduction

Welcome to the lesson on the RNA-Seq analysis procedure! A potent method for analyzing transcriptome patterns and measuring gene expression levels is RNA-Seq. In this lesson, we'll walk you through each crucial step of evaluating RNA-Seq data, including data preparation and assessing differential gene expression. This lesson will provide you a strong basis to build on whether you are new to bioinformatics or want to improve your abilities in RNA-Seq analysis.

In this guide, you will be able to understand the process that includes obtaining raw data (fastq files), preprocessing, alignment of reads to a reference genome, counting the expression levels of genes, normalization techniques, Differential expression analysis, functional annotation and interpretation through visualization techniques. 

You will find code, step-by-step directions, and helpful hints to aid you in completing the analysis process throughout this guide book. By the time you finish this exercise, you will be well-versed in the fundamental ideas and resources needed for R-based RNA-Seq data processing.

So let's get started and begin this fascinating adventure into the world of RNA-Seq analysis!

***

## Prepare Data

### Download read data

For this analysis, Mouse mammary data (fastq files) can be downloaded from [figshare](https://figshare.com/s/f5d63d8c265a05618137) or [ncbi (GSE60450)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450)

I have used the package Rsubread in this analysis. You can find more information about this package [here](https://bioconductor.org/packages/release/bioc/html/Rsubread.html).

The raw reads downloaded from SRA are in .sra format. The sra toolkit from NCBI is used to convert the .sra files to .fastq files using the fastq-dump command. You can find more details on sra-toolkit [here](https://hpc.nih.gov/apps/sratoolkit.html). 

### Downloading genome files

For this analysis, the index files for chromosome 1 of mouse genome build mm10 is provided in the data files from [figshare](https://figshare.com/s/f5d63d8c265a05618137). However, you can obtain full genome fasta files for various genomes from the [UCSC genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu) or [here](https://hgdownload.soe.ucsc.edu/downloads.html); or from [NCBI](https://www.ncbi.nlm.nih.gov/genome); or from [ENSEMBL](https://www.ensembl.org/info/data/ftp/index.html).

### Import Data for Alignment

Mapping reads to the genome is a very important task, and many different aligners are available, such as [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537), [topHat](https://academic.oup.com/bioinformatics/article/25/9/1105/203994), [bowtie](https://www.nature.com/articles/nmeth.1923) and [Rsubread](https://academic.oup.com/nar/article/41/10/e108/1075719). **Rsubread** is the only aligner that can run in R. The majority of alignment programs run on Linux and need a lot of processing power. The majority of mapping activities demand more powerful processors than a typical laptop, hence read mapping is often carried out on a server that has linux-like environment. We will only map 1,000 reads from each sample of our mouse lactation dataset from [@Fu2015](https://www.nature.com/articles/ncb3117), but because of time and limitation of computational power, I will only be mapping to chromosome 1.

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

***

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

***

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

***

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

***

## Tips and Tricks

{{< columns >}}

The Phred offset determines the encoding for the base-calling quality string in the
FASTQ file. For the Illumina 1.8 format onwards, this encoding is set at +33.
However, older formats may use a +64 encoding. Users should ensure that the
correct encoding is specified during alignment. If unsure, one can examine the
first several quality strings in the FASTQ file. A good rule of thumb is to check
whether lower-case letters are present (+64 encoding) or absent (+33).

<--->


`featureCounts` requires gene annotation specifying the genomic start and end position of each exon of each gene.   
`Rsubread` contains built-in gene annotation for mouse and human. For other species, users will need to read in a data frame in GTF format to define the genes and exons. Users can also specify a custom annotation file in SAF format. 
{{< /columns >}}

***

## Preprocessing of Count Data

You may want to download the Data files from [Figshare](https://figshare.com/s/1d788fd384d33e913a2a) and place them in your `/data` directory.

The data for this tutorial comes from a Nature Cell Biology paper, [*EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival*](http://www.ncbi.nlm.nih.gov/pubmed/25730472). This study is focused on the expression profiles of **B**asal stem-cell enriched cells and committed **L**uminal cells in the mammary gland of virgin, pregnant and lactating mice. Each group comprise of two biological replicates.

Let's start of with installing and loading all the packages required to analyse the data.

```r

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

```

Basic information about the samples can be found in the **sampleinfo** file contains.

```r

# Read the sample information into R
sampleinfo <- read.delim("data/SampleInfo.txt")
View(sampleinfo)
sampleinfo

```

Let's take a look at the count data. In RStudio we can used either `head` or `View` command to have a glimpse of the dataframe. The `dim` command shows dimensions of the matrix (how many rows and columns the data frame has).

```r

# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
head(seqdata)
View(seqdata)
dim(seqdata)

```

Upto this point we have `seqdata` variable that contains count information about genes (one gene per row), the first column is the Entrez gene id, the second column is gene length and the remaining columns counts of reads aligning to the gene in each sample. Detailsabout the sample can be found in **GSE60450_series_matrix.txt** from the [GEO website](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450)). 


### Data Formatting

This count matrix needs to be manipulated and formatted so that it will be suitable for downstream analysis.   
1. Make a new matrix that only contains the counts   
2. Store the gene identifiers as rownames

Let's create a new data object, `countdata`, that contains only the counts for the 12 samples.  

```r
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
```

Here I will also shorten the column names for easy understanding about each sample.  

```r

# Extract the first 7 characters and use these as the colnames
# substr() does the job
colnames(countdata) <- substr(colnames(countdata), 1, 7)
colnames(countdata)

```

We should always keep the column names in count matrix same as `SampleName` in the `sampleinfo` file and in the same order. This will help in making the contrasts for the downstream analysis.


```r
# Check the column names and order in count data and sampleinfo
colnames(countdata)==sampleinfo$SampleName

all(colnames(countdata)==sampleinfo$SampleName)

```

***

## Data Filtering

A few of the statistical approximations that are utilized in the downstream analysis are affected by genes with extremely low counts in all libraries, which holds no value/significance for differential expression, rather serves as noise data. And they increase the burden of multiple testing when evaluating false discovery rates, decreasing the ability to identify genes with differential expression. Hence, these genes should be **eliminated**.

There are numerous approaches typically used in RNA-Seq analysis to eliminate genes with low expression. These strategies seek to exclude genes with extremely low read counts or expression levels, as they may not contribute much to the total biological signal and may bring noise into downstream analysis. Here are some of the frequently used approaches for eliminating low-expression genes:

1. Count-based filtering:
    * Thresholding: Set a minimum threshold for read counts below which genes are considered as lowly expressed and filtered out.
    * Proportional filtering: Remove genes that have low expression relative to the total library size or the expression distribution across samples. This can be done using quantile-based methods or by applying filters based on a certain percentile cutoff.   

2. Normalization-based filtering:
    * Trimmed Mean of M-values (TMM): Use TMM normalization to estimate the scaling factors for each sample, which can then be used to filter out genes with low expression values relative to the library size.
    * Upper Quartile (UQ) or Median of Ratios (MOR): Normalize the expression values using the UQ or MOR method, and filter out genes with expression values below a certain threshold.  

3. Variance-based filtering:
    * Coefficient of Variation (CV): Calculate the CV for each gene across samples and remove genes with low expression variability, as they may represent noise or technical artifacts.
    * Dispersion estimation: Estimate the dispersion of gene expression using methods like the negative binomial or Poisson models, and filter out genes with low dispersion.  

4. Machine learning-based filtering:
    * Use machine learning algorithms, such as random forests or support vector machines, to predict the significance or importance of genes based on their expression patterns. Genes with low predicted significance can be filtered out.

It's important to note that the choice of method may depend on the specific characteristics of your data, experimental design, and the goals of your analysis. It's often recommended to combine multiple filtering methods or perform exploratory data analysis to determine the appropriate cutoffs for gene filtering.

Additionally, it's crucial to consider the impact of gene filtering on downstream analyses and interpretation of results. Filtering out lowly expressed genes can affect differential expression analysis, pathway enrichment analysis, and other downstream analyses. Therefore, it's essential to carefully evaluate the trade-off between retaining biological signal and removing noise when applying gene expression filtering methods.


For this problem, I prefer filtering on a minimum counts per million threshold present in at least two samples where there are biological replicates in each group, which is the case in this case because each group has two samples. Here, the sample size for each group is two. I opted to preserve genes in this dataset if their counts-per-million (CPM) expression levels in at least two samples is more than 0.5. (Normal practice is `remove if CPM < 1`)  

We'll use the `cpm` function from the [**edgeR** library](https://academic.oup.com/bioinformatics/article/26/1/139/182458) to generate the CPM values and then filter. 

*It is important to note that by converting to CPMs, we are normalizing for the various sequencing depths for each sample.*

```r

# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

col1sum <- sum(countdata[,1])/1000000
countdata[1,1]/col1sum

```

Now we want to Check how many genes are left after we filter out the low expressed genes based on CPM value less than 0.5. 
For this, calculate a threshold for gene expression based on counts per million (CPM), where genes with CPM values greater than 0.5 are considered to have sufficient expression.

```r

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
rowSums(head(thresh))
table(rowSums(thresh))

```

Then determine which genes meet this criterion in at least two samples and creates a logical vector indicating which genes to keep. Finally, subset the original count data to retain only the rows (genes) that meet the expression threshold and have at least two instances of meeting the threshold across all samples.

```r
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# summary(keep)

# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
# dim(countdata)
# dim(counts.keep)
```

If the count falls below a certain threshold, it is regarded very low, suggesting that the related gene is not expressed in that sample. As there are two replicas in each group, there is a demand for expression in two or more libraries. This assures that if a gene is exclusively expressed in one group, it will be maintained. For bigger libraries, lower CPM criteria are generally adequate. A decent threshold may be set as a general rule by determining the CPM that corresponds to a count of 10, which in this case is roughly 0.5. CPMs should be used instead of counts since counts do not account for changes in library sizes across samples.


```r
# Check if CPM of 0.5 correspond to count of about 10-15
# for the first sample
plot(myCPM[,1], countdata[,1], xlab="CPM", ylab="Raw Count", main=colnames(myCPM)[1])
# Add a vertical line at 0.5 CPM
abline(v=0.5)


# limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1], xlab="CPM", ylab="Raw Count", main=colnames(myCPM)[1], 
     ylim=c(0,50), xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

```

***

## Differential expression

Prepare the data for Differential expression analysis using edgeR package.  
```r

# Convert to an edgeR object
dgeObj <- DGEList(counts.keep)

# Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)

# sample information
sampleinfo <- read.delim("data/SampleInfo_Corrected.txt")
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group

```


### Design matrix

We start with making a design matrix for the groupings. We have two variables: status and cell type. On these basis, We will fit two models under two assumptions: no interaction and interaction of these two components. 

1. A model that has only major effects and no interaction. 
    * The fundamental assumption here is that the effect of the status is the same in both types of cells.

```r

# Create the two variables
group <- as.character(group)
type <- sapply(strsplit(group, ".", fixed=T), function(x) x[1])
status <- sapply(strsplit(group, ".", fixed=T), function(x) x[2])

# Specify a design matrix with an intercept term
design <- model.matrix(~ type + status)
design

```

### Data exploration

The MDS plot provides a visualization of the distances between samples based on their biological coefficient of variation (BCV), which reflects the variability of gene expression across samples. By observing the MDS plot, we can assess the quality of the data in terms of sample similarity and clustering patterns. If the samples cluster tightly together and show distinct groups, it suggests good data quality with clear differences between conditions. On the other hand, if the samples are scattered and lack clear grouping, it may indicate higher variability or confounding factors affecting the gene expression patterns.

Regarding the anticipation of the importance of the interaction term, it would depend on the specific research question and design of the experiment. If the interaction term represents a biologically relevant interaction between two factors, and if the MDS plot reveals distinct groupings of samples based on those factors, it is likely that the interaction term will be important in explaining the variation in gene expression. However, further statistical analysis and hypothesis testing would be necessary to determine the significance and impact of the interaction term in the data.

```r

plotMDS(dgeObj, labels=group, cex=0.75, xlim=c(-4, 5))

```

### Estimating the dispersion

The common dispersion is a measure that estimates the average biological coefficient of variation (BCV) across all genes in the dataset. It provides an indication of the overall variability or dispersion of gene expression levels within the dataset. By averaging the BCV values of all genes, the common dispersion provides a summary measure of the typical level of variation observed in gene expression across samples. It helps researchers understand the overall level of variability in the dataset and can be useful for quality assessment and subsequent statistical analysis of differential gene expression.

```r

dgeObj <- estimateCommonDisp(dgeObj)

```

Estimating gene-wise dispersion involves calculating the variability or dispersion of gene expression levels for each individual gene in the dataset. However, in addition to considering the raw expression values, this estimation takes into account a possible trend with the average count size. This means that it considers the relationship between the dispersion and the average expression level of each gene. By accounting for this trend, the gene-wise dispersion estimates provide a more accurate assessment of the variability specific to each gene, considering its expression level. This information is valuable for downstream statistical analyses, such as identifying differentially expressed genes, as it helps to account for the varying dispersion patterns across genes and provides a more precise modeling of the data.

```r

dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)

```

Plot the estimated dispersions:

```r

plotBCV(dgeObj)

```


### Differential Expression

At this point we need to fit the model using *glmFit()* functionction. 

```r

# Fit the linear model
fit <- glmFit(dgeObj, design)
# names(fit)
# head(coef(fit))

```

Likelihood ratio tests are conducted to compare the goodness-of-fit of two competing models: luminal and basal. These models represent different gene expression patterns associated with different biological characteristics. By comparing the likelihoods of these models, the test determines which model better explains the observed gene expression data. The top genes are the ones that exhibit the most significant differences in expression between the luminal and basal subtypes, based on the likelihood ratio test. These genes are considered the most relevant and informative in distinguishing between the two subtypes and can potentially provide insights into the underlying biological processes and characteristics associated with each subtype.

```r
lrt.B_L <- glmLRT(fit, coef=2)
topTags(lrt.B_L)

```

### Contrasts

When comparing the gene expression between two groups, such as pregnant and virgin, where we don't have a specific parameter for testing the hypothesis, we can construct a contrast. A contrast is a linear combination of the coefficients in a statistical model that allows us to specifically test the difference between the two groups of interest. In this case, we would define a contrast that captures the difference in gene expression between pregnant and virgin samples. By applying this contrast, we can perform statistical tests, such as t-tests or likelihood ratio tests, to identify genes that are differentially expressed between the two groups. The contrast helps us focus on the specific comparison of interest and extract relevant information from the data for hypothesis testing and inference.

```r

P_V <- makeContrasts(statuspregnant-statusvirgin, levels=design)
lrt.P_V <- glmLRT(fit, contrast=P_V)
topTags(lrt.P_V)

```

***


```r

results <- as.data.frame(topTags(lrt.BvsL,n = Inf))
# results
# dim(results)

```

The plotSmear function in the edgeR package is a visualization tool used to display the results of a differential expression analysis. It generates a plot that represents the log-fold change (logFC) on the x-axis and the log-counts per million (logCPM) on the y-axis. The plot helps in identifying genes that are differentially expressed by highlighting them. By examining the position and distribution of points on the plot, researchers can gain insights into the magnitude and direction of gene expression changes between different conditions or groups. This is an equivalent of [*MA-plot* for microarray data](https://en.wikipedia.org/wiki/MA_plot) which is commonly used in microarray data analysis, which visualizes the log-fold change against the average intensity of gene expression.


```r

summary(de <- decideTestsDGE(lrt.BvsL))
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.BvsL, de.tags=detags)

```

To add labels to highlight specific genes on the plotSmear plot, you would need to map the identifiers from the edgeR output to more familiar gene names or symbols. This mapping can be done using external databases or annotations that provide the correspondence between the identifiers and gene names.


## Adding annotation to the edgeR results
***
  

There are other techniques for adding annotation to your analysis, and in this case, we will use the org.Mm.eg.db package. You may also use [BioMart](http://www.biomart.org/), which acts as an interface to the BioMart database resource. Although BioMart provides a greater range of annotations, organism-level packages such as org.Mm.eg.db are more aligned with the Bioconductor methodology.  When deciding between the two approaches, you may want to consider the level of annotation detail you require and choose the option that best fits your needs.


```r

ann <- select(org.Mm.eg.db,keys=rownames(results),columns=c("ENTREZID","SYMBOL","GENENAME"))

```

To incorporate the annotation information into the results data frame, you can merge or bind the annotation data based on matching identifiers. However, it is important to be cautious when dealing with a one-to-many mapping scenario, as simply appending the annotation to the fit object may not be appropriate in such cases.


```r

results.annotated <- cbind(results, ann)
# results.annotated

# save the results
# write.csv(results.annotated,file="B.PregVsLacResults.csv",row.names=FALSE)
```


> When determining which genes are considered significant in a differential expression analysis, it is common practice to use an adjusted p-value cutoff of 0.05 rather than the raw p-value. This is necessary because when testing a large number of genes, such as in the case of over 15,000 genes, the probability of identifying differentially expressed genes by chance alone increases significantly. By controlling the false discovery rate (FDR), indicated by the adjusted p-value column in the results table, we aim to limit the number of false positives. The decideTests function helps identify significant genes at a 5% FDR, providing a threshold for determining differential expression.



```r

# Annotate the plotsmear figure
plotSmear(lrt.B_L, de.tags=detags)
N <- 200
text(results.annotated$logCPM[1:N],results.annotated$logFC[1:N],labels = results.annotated$SYMBOL[1:N],col="blue")

```

We can visualize the results of a differential expression analysis using a volcano plot, which provides a graphical representation of the significance of gene expression changes. In a volcano plot, the y-axis represents a measure of significance (such as the negative logarithm of the adjusted p-value), while the x-axis represents the fold-change in gene expression between experimental conditions. This plot allows us to easily identify genes that exhibit both significant changes in expression and large fold-changes, as they appear farther away from the center of the plot.

```r

signif <- -log10(results.annotated$FDR)
plot(results.annotated$logFC,signif,pch=16)
points(results.annotated[detags,"logFC"],-log10(results.annotated[detags,"FDR"]),pch=16,col="red")

```

***

## Retrieving Genomic Locations

Install the required package.

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)

```

Create an object for the taxonomy database.

```r

tx <- TxDb.Mmusculus.UCSC.mm10.knownGene
# columns(tx)

```

The txdb packages' easy integration with GenomicRanges, a commonly used object type in Bioconductor for manipulating genomic intervals, is one of their advantages. We can use GenomicRanges to conduct operations on intervals such as identifying overlaps between areas and counting the number of intervals that intersect with a specific region. This object type depicts genomic intervals by identifying the chromosome, start location, and end position of each section, making working with genomic data in R more straightforward.

```r

library(GenomicRanges)

```

***

## Retrieving Gene Coordinates as GenomicRanges


The output of `exonsBy` is a list, where each item in the list is the exon co-ordinates of a particular gene. 

```r

exo <- exonsBy(tx,"gene")
exo

```

***

## Exporting tracks

The outcomes of a Bioconductor investigation may be saved in a browser, enabling interactive analysis, integration with other data types, and sharing with collaborators. Making a browser track that shows the locations of genes with differential expression is one approach to do this. The bed format, which enables us to exhibit genomic ranges together with extra data from our study, such as fold-change and significance, may be used to do this. With this method, you may study the genomic context of differentially expressed genes in a visually appealing and interactive fashion. It also makes it easier to combine various data types for thorough research.

Create a data frame for differentially expressed genes.

```r

sigGenes <- results.annotated[detags,]
# sigGenes

```

To start, we use the 'range' function to get a single range for each gene, which we then transform using 'unlist' into a more useful object.

```r

exoRanges <- unlist(range(exo))
sigRegions <- exoRanges[na.omit(match(sigGenes$ENTREZID, names(exoRanges)))]
# sigRegions

```

We can also color each range using the .bed format based on particular property of the analysis such as its fold chang or its direction of dysregulation. we can also utilize the metadata with the help of  GenomicRanges.

```r

mcols(sigRegions) <- sigGenes[match(names(sigRegions), rownames(sigGenes)),]
# sigRegions

```


```r
sigRegions[order(sigRegions$LR,decreasing = TRUE)]
```


```r
seqlevels(sigRegions)
sigRegions <- keepSeqlevels(sigRegions, paste0("chr", c(1:19,"X","Y")))
```

To enhance the visualization of the browser track, we can create a score based on the p-values. One common approach is to use the negative logarithm (base 10) of the adjusted p-values as the score. This transformation helps emphasize significant results by amplifying their magnitude. Additionally, we can assign a color scheme to the regions based on the fold-change values. This allows for a visual representation of the direction and magnitude of the fold-change, providing further insights into the differential expression patterns.


```r

Score <- -log10(sigRegions$FDR)

```

Adjust the color scheme by, creating a color ramp palette rbPal using the colorRampPalette function, specifying the colors "red" and "blue". The fold-change values (logFC) are then adjusted to be within the range of -5 to 5 using the pmax and pmin functions. Next, the adjusted fold-change values are divided into 10 equally spaced intervals using the cut function, and the corresponding colors from the color ramp palette are assigned to each interval. The resulting colors are stored in the variable Col. This process allows for the generation of a color scheme based on the fold-change values for the regions in the browser track.

```r
rbPal <-colorRampPalette(c("red", "blue"))
logfc <- pmax(sigRegions$logFC, -5)
logfc <- pmin(logfc , 5)

Col <- rbPal(10)[as.numeric(cut(logfc, breaks = 10))]
```

save the generated colors and score as the GRanges object as `score` and `itemRgb` columns respectively, which later can be used in the browser track. Export the signifcant results from the DE analysis as a `.bed` track and use the file as input for [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/).


```r

mcols(sigRegions)$score <- Score
mcols(sigRegions)$itemRgb <- Col
sigRegions
library(rtracklayer)
export(sigRegions , con = "topHits.bed")

```

## Extracting Reads

In Bioconductor, the Rsamtools package provides a low-level interface for working with bam/sam files. Additionally, the GenomicAlignments package allows for efficient retrieval of reads that map to specific genomic regions. These packages enable users to access and manipulate the aligned reads in their analysis, providing greater flexibility and control over the data exploration process.

```r

library(GenomicAlignments)

```

Make sure that the all bam files are named properly, along with their corresponding index file (`.bai`) is present in the same location.


```r

list.files("bam/")

```

The readGAlignments function effectively gets just the reads that align to a specific genomic region of interest by using the index file associated with the bam file. The output of readGAlignments comprises the CIGAR string, which contains information on matches, insertions, and deletions in the alignment, as well as details about the genomic location of each aligned read.

```r

generegion <- exo[["100421"]]

my.reads <- readGAlignments(file="bam/MCL1.DG.bam",
                       param=ScanBamParam(which=generegion))
# my.reads

```
And you can also use the QC information in the same manner;

```r

my.reads <- readGAlignments(file="bam/MCL1.DG.bam",
                       param=ScanBamParam(which=generegion,
                                          what=c("seq","mapq","flag")))
# my.reads

```

You can find more information on flags [here](https://broadinstitute.github.io/picard/explain-flags.html).


***

## Gene Set Testing

It might be difficult and time-consuming to evaluate the findings of gene set testing when working with a huge number of differentially expressed genes. Researchers frequently use gene set testing to solve this issue and learn more about the pathways or gene networks that these differentially expressed genes are connected to.

Gene set testing may be divided into two categories:   
1. competitive tests   
2. Self-contained tests  

Competitive gene set tests look at whether the differentially expressed genes are overrepresented within a certain gene set in comparison to all other genes in the experiment. For this we prefer tools like GOseq and camera. This kind of test aids in determining if the pathway's or gene set's genes are enriched among the genes with differential expression.

On the other hand, self-contained assays, such as the ROAST technique, seek to determine if a group of genes or a pathway as a whole exhibits differential expression. These tests evaluate the overall pattern of gene expression, revealing information on the behavior of the pathway's genes as a whole.

In pathway analysis, both competitive and self-contained gene set tests are useful, and which test to apply depends on the particular research topic and the kind of analysis.  

***


## Competitive Gene Set Analysis

### GOseq Method

GOseq is a specialized method used for conducting Gene Ontology (GO) analysis specifically tailored for RNA-seq data. It is designed to address the potential bias introduced by differences in gene length when detecting over-representation of GO terms.

In RNA-seq experiments, gene length can affect the number of reads mapped to a gene, leading to a potential bias in the analysis of GO term enrichment. GOseq takes into account this bias by incorporating gene length information as a correction factor in the statistical analysis. This correction helps to mitigate the impact of gene length on the detection of enriched GO terms, resulting in more accurate and reliable results.

By using GOseq, researchers can perform GO analysis on RNA-seq data while considering the inherent bias introduced by gene length, leading to more robust and biologically meaningful interpretations of the enriched GO terms.


```r

# enlist all DE genes
results <- as.data.frame(topTags(lrt.B_L, n = Inf))
# print(head(results))

# Filter out significant DEGs on FDR:
genes <- results$FDR < 0.01

# Add gene names to that list:
names(genes) <- rownames(results)

# print(head(genes))

```

The Probability Weighting Function (PWF) is a crucial step in the GOseq package for conducting Gene Ontology (GO) analysis on RNA-seq data. It estimates the probability distribution of gene lengths to account for the bias introduced by gene length variation. By fitting the PWF, the GOseq package adjusts for this bias and allows for accurate assessment of GO term enrichment by considering the expected distribution of gene lengths within the analysis. This step ensures more reliable and biologically relevant results in GO analysis of RNA-seq data.

```r
library(goseq)

#print(supportedGeneIDs())
#print(supportedGenomes())

pwf <- nullp(genes, "mm10","knownGene")

```

Conduct gene set enrichment analysis:

```r
# goseq
go.results <- goseq(pwf, "mm10","knownGene")
# go.results
```

### FGSEA Method

The fgsea method is a gene set enrichment analysis approach that assesses the enrichment of gene sets in high-throughput genomic data. It uses a pre-ranked gene list based on a specified statistic, such as differential expression scores or fold changes. The fgsea algorithm calculates an enrichment score for each gene set by comparing the positions of genes from the gene set in the ranked list. This method accounts for the gene set size and gene set-specific statistics to identify gene sets that are significantly enriched or depleted in the data, providing valuable insights into biological processes and pathways associated with the analyzed genes.

```r

library(fgsea)

# Create ranks
results.ord <- results[ order(-results[,"logFC"]), ]
# head(results.ord)
ranks <- results.ord$logFC
names(ranks) <- rownames(results.ord)
# head(ranks)


#plot(ranks)
barplot(ranks)
```

Load pathways.

```r

load("data/mouse_H_v5.rdata")
pathways <- Mm.H

```

Conduct analysis:

```r

# fgsea
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nperm=1000)
# class(fgseaRes)
# dim(fgseaRes)
#head(fgseaRes)
# head(fgseaRes[order(padj), ])
```

Plot the results for specific pathways (like HALLMARK_MYOGENESIS). Find rank of the 'HALLMARK_MYOGENESIS' pathway genes in the sorted genes.

```r

# We will create a barplot of logFC for the sorted genes and add one vertical red bar for each gene in the 'HALLMARK_MYOGENESIS' pathway

#pathways[["HALLMARK_MYOGENESIS"]]

tmpInd <- match(pathways[["HALLMARK_MYOGENESIS"]],names(ranks))
tmpInd <- tmpInd[!is.na(tmpInd)]

#tmpInd

ranks2 <- rep(0,length(ranks))
ranks2[tmpInd] <- ranks[tmpInd]

barplot(ranks2)
```

Create enrichment score plot:

```r

plotEnrichment(pathways[["HALLMARK_MYOGENESIS"]],
               ranks)
               
```

You can find more information through [GSEA article](http://www.pnas.org/content/102/43/15545.full).

Select top pathways and plot outcome for all these:

```r
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)
# plotGseaTable
```

### CAMERA Method

CAMERA (Competitive Gene Set Analysis) is a gene set testing method that evaluates the enrichment of gene sets using pre-defined gene sets curated by the Broad Institute ([MSigDB](http://software.broadinstitute.org/gsea/msigdb)). It addresses the question of whether a given gene set is differentially expressed as a whole. The analysis begins by fitting a linear model for each gene set, incorporating the experimental design and other relevant factors. CAMERA then computes the empirical Bayes shrinkage estimators of the coefficients, which capture the differential expression information. Next, it calculates a score for each gene set based on the contrast between the coefficients of the genes within the set and those outside the set. This score reflects the overall differential expression pattern of the gene set. To assess the significance of the score, CAMERA employs a permutation-based procedure to generate null distributions and computes p-values accordingly. Additionally, CAMERA allows for multiple testing correction using methods such as the Benjamini-Hochberg procedure. The result is a list of gene sets ranked by their significance, enabling researchers to identify the gene sets that are most likely to be differentially expressed. By using curated gene sets from the Broad Institute, CAMERA provides a valuable tool for investigating the functional implications and biological processes associated with gene expression changes in high-throughput genomic data.

```r

# camera.DGEList

# Load in the mouse c2 gene sets
# The R object is called Mm.c2
load("data/mouse_c2_v5.rdata")

```

The ids2indices function allows us to map the Entrez gene IDs between the list of gene sets and our DGEList object. By providing the Entrez gene IDs from the gene sets and the row names of our DGEList object as inputs, the function returns the corresponding indices or positions of the genes in the DGEList object. This mapping enables us to match the genes between the gene sets and our differential gene expression data, facilitating downstream analysis and interpretation of the results.

```r

c2.ind <- ids2indices(Mm.c2, rownames(dgeObj$counts))

```

In the CAMERA gene set testing method, the input includes the DGEList object dgeObj, the indexed list of gene sets c2.ind, the design matrix, and the contrast being tested, among other arguments. By default, CAMERA estimates the correlation for each gene set separately. However, it is often beneficial to set a small inter-gene correlation, typically around 0.05, using the inter.gene.cor argument. This helps to account for potential correlations among genes within the same gene set and improves the accuracy of the testing procedure.

```r

# Conduct analysis for the luminal-vs-basal contrast:
group <- as.character(group)
type <- sapply(strsplit(group, ".", fixed=T), function(x) x[1])
status <- sapply(strsplit(group, ".", fixed=T), function(x) x[2])
# Specify a design matrix without an intercept term
design <- model.matrix(~ type + status)

# Check contrasts:
# print(colnames(design))

# Run analysis:
gst.camera <- camera.DGEList(dgeObj,index=c2.ind,design=design,contrast=2,inter.gene.cor=0.05)

```

Save camera results to a csv file to open in excel.

```r

write.csv(gst.camera,file="gst_LumVsBas.csv")

```

   
***



## Self-Contained Gene Set Analysis

### ROAST Method

[ROAST  (Rotation Gene Set Testing)](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq401) is a gene set testing method that aims to assess whether a gene set is differentially expressed as a whole. It is particularly useful for small sample sizes or when there is heterogeneity in the expression patterns within a gene set. ROAST uses a rotation-based approach to generate a null distribution by randomly rotating the sample labels. This accounts for the dependence structure among genes within the gene set. The method estimates the proportion of rotations that result in a more extreme test statistic than the observed data, allowing for the calculation of p-values. ROAST can handle both competitive and self-contained gene sets and offers flexibility in defining the gene sets of interest. The method provides robust statistical inference, taking into account the complex correlation structure among genes. It is implemented in the limma package in R and can be applied to various types of high-dimensional genomic data beyond RNA-seq. Overall, ROAST is a valuable tool for gene set analysis, aiding in the identification of functionally relevant gene sets and providing insights into the biological processes underlying the observed differential expression patterns.


Let's see if there are any MYC signalling pathways in MsigDB C2 collection. We can do this with the `grep` command on the names of the gene sets.

```r

grep("MYC_",names(c2.ind))

# Let's save these so that we can subset c2.ind to test all gene sets with MYC in the name
myc <- grep("MYC_",names(c2.ind))

# What are these pathways called?
names(c2.ind)[myc]

```

We can use ROAST to see if these MYC related gene sets tend to be differentially expressed. 


```r

myc.rst <- roast(dgeObj,index=c2.ind[myc],design=design,contrast=3,nrot=999)
myc.rst[1:15,]

```

The output of ROAST analysis consists of a table with several columns providing information about each gene set. Each row in the table represents a single gene set. The "NGenes" column indicates the number of genes included in each set. The "PropDown" and "PropUp" columns indicate the proportions of genes in the set that are down-regulated and up-regulated, respectively, with fold changes greater than 2. The "Direction" column shows the net direction of change in the gene set, determined from the significance of changes in each direction.

The "PValue" column provides evidence for whether the majority of genes in the set are differentially expressed (DE) in the specified direction (either up or down). A lower p-value indicates stronger evidence of DE in the given direction. The "PValue.Mixed" column provides evidence for whether the majority of genes in the set are DE in any direction, regardless of up or down regulation.

False discovery rates (FDRs) are computed from the corresponding p-values across all gene sets. The FDR represents the expected proportion of false positives among the significant results. It is a measure of the statistical significance, taking into account the multiple hypothesis testing conducted in gene set analysis.

Overall, the output of ROAST analysis provides valuable information about each gene set, including the number of genes, proportions of up- and down-regulated genes, net direction of change, and statistical evidence for differential expression. It helps researchers identify gene sets that are significantly enriched for differentially expressed genes, providing insights into the underlying biological processes related to the observed expression changes.


In the context of ROAST analysis, a common scenario is to utilize a set of differentially expressed (DE) genes that have been identified from a separate analysis conducted on an independent dataset. This set of DE genes serves as a reference or benchmark for evaluating the similarity of gene expression changes in the current dataset under a specific contrast of interest.

By applying ROAST to the current dataset and the chosen contrast, it becomes possible to assess whether similar patterns of gene expression changes, as observed in the reference set of DE genes, are also present in the contrast of interest in the current dataset. In other words, ROAST evaluates whether there is a significant overlap or concordance between the gene expression changes identified in the reference set and the genes showing differential expression in the current dataset under the specified contrast.

This approach allows researchers to validate and compare the consistency of gene expression changes across different datasets or experimental conditions. It provides a means to assess whether the identified DE genes in the current dataset exhibit similar expression patterns to those observed in a previously defined set of DE genes, thereby supporting the reproducibility or generalizability of the gene expression changes across different contexts.

***

