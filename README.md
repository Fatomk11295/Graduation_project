# A comparative analyses between GTEX healthy samples obtained from Bone marrow and Whole blood tissues in comparison to TCGA Acute Myeloid Leukemia database.


## Project description:


- Biological questions

  Which samples are more reliable to compare to TCGA AML, samples collected from Bone marrow or whole blood?

  Several articles used samples collected from Bone marrow as control when studying Acute Myeloid Leukemia (AML) gene expression profile, while others used samples obtained from Whole blood. Therefore, the first question we asked ourselves was “In a world full of data, which data set is the best to use ?”. If we want to study the AML gene expression profile to find a new targeting pathway to treat AML we have to use samples as a case. Selecting the optimum control sample will form the foundation of our study. 


-----------

### **Step 1: Collect the best available data**

To initiate the analysis we downloaded 3 databases :

1. [GTEX phenotype database](https://xenabrowser.net/datapages/?dataset=GTEX_phenotype&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) 
2. [GTEX raw counts database](https://xenabrowser.net/datapages/?dataset=gtex_gene_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
3. [TCGA AML raw counts database](https://xenabrowser.net/datapages/?dataset=TCGA.LAML.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

GTEX raw counts database were subsetted based on tissue type (Bone marrow or Whole blood). The following steps were the same for each sub-dataset.

### **Step 2: Preprocessing the dataframe**

The finalized dataset (which contains cases and controls samples) underwent a few preprocessing steps to fit the data into the suggested analysis pipeline. The preprocessing steps were :

1. *Unit conversion* - from log2(x+1) to raw counts

```python
def convert_unit(df):
  '''This function convert normalized counts log2(x+1) to raw counts'''
  #The equation is as following:
  # 2 ** value (right) = x + 1 (left)
  '''Step 1 : Find x value (Reads value before normalization)'''
  right_array = pow(df.values,2)
  '''Step 2 : substract all values in right array by 1'''
  raw_array = right_array - 1
  '''Step 3 : return the dataframe with raw counts'''
  rawCounts_df = pd.DataFrame(raw_array, index=df.index,
                        columns=df.columns)
  return rawCounts_df
```
2. *Remove duplicated genes*

```python
rawCounts_df_WithoutDuplicates = rawCounts_df.drop_duplicates()
```

3. *Label the samples* - based on it type into 'Case' or 'Control'

```python
#Create a list to label sample as case or control 
kind= []
for x in df_Transposed.index:
  if x[0] == 'K': # K if you are using samples from Whole blood / G if you are using samples from Bone Marrow
    kind.append('Control')
  if x[0] == 'T':
    kind.append('Case')
#Define the new column Type based on kind list
df_Transposed['Type'] = kind
```

**The python code used for the preprocessing steps can be found [here](https://github.com/Fatomk11295/Graduation_project/blob/main/prepare_the_data.py).**

----

The final dataset contains **20148 genes** with total samples **221** when control is BM-dataset, and **488** samples for WB-dataset.

-----------------------
### **Step 3: Differential gene expression analysis**

For this analysis we used [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) package. 

```r
# Load essential libraries
library(DESeq2)
library(readr)
library(SummarizedExperiment)
library(ggplot2)

# Step 1: prepare data for DESeq
# link to working directory and read data files

colData=read.csv("Sample_list.csv", row.names = 1)
counts_data=read.csv("counts_data_Ready.csv",header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors =  F)

# making sure the row names in colData matches to column names in counts_data

all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?

all(colnames(counts_data) == rownames(colData))

# Step 2: construct a DESeqDataSet object

#Round gene expression values and convert them into intgers instead of float
counts_data[1:4,1:5]
range(counts_data)
counts_data<- round(counts_data, digits = 0)
range(counts_data)

# construct the dataset object for the DESeq analysis
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design =~Type)

# Set the factor level to compare between control and case samples, means, telling DESeq which level to compare against , here, we used "control" as our reference level

dds$Type <- relevel(dds$Type, ref = "Control")

## Check which fitType is the best 

ddsPlot <- estimateSizeFactors(dds)
par <- estimateDispersions(ddsPlot, fitType = "parametric")
meanfit <- estimateDispersions(ddsPlot, fitType = "mean")
loc <- estimateDispersions(ddsPlot, fitType = "local")
plotDispEsts(par, main= "dispEst: parametric")
plotDispEsts(meanfit, main= "dispEst: mean")
plotDispEsts(loc, main= "dispEst: local")

# Step 3: Run DESeq 
# For our data we chose the fitType local

dds <- DESeq(dds, fitType = "local" )
final_results <- results(dds)

# Explore Results

summary(final_results)

# Subset the data based on pvalue -only take below threshold 0.05-
final_results0.05 <- results(dds, alpha = 0.05)
summary(final_results0.05)

# contrasts
resultsNames(dds)

# MA plot
plotMA(final_results0.05, colNonSig = "salmon", colSig = "orchid1")

# Saving adjusted p_value for each gene to csv file
write.csv(final_results0.05, "DESeqResults_Jul12.2022_localFitType_BoneMarrow.csv")

# Volcano Plot

df = read.csv("DESeqResults_Jul12.2022_localFitType_BoneMarrow.csv")

# The basic scatter plot: x is "log2FoldChange", y is "padj"
ggplot(data=df, aes(x=log2FoldChange, y=padj)) + geom_point()

# Convert directly in the aes()
p <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point()

# Add more simple "theme"
p <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() + theme_minimal()

# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-1, 1), col="salmon") + geom_hline(yintercept=-log10(0.05), col="turquoise ")

# The significantly differential expressed genes are the ones found in the upper-left and upper-right corners.

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
# add a column of NOs
df$diffexpressed <- "NO"
# if log2Foldchange > 1and padj < 0.05, set as "UP" 
df$diffexpressed[df$log2FoldChange > 1 & df$padj < 0.05] <- "UP"
# if log2Foldchange < -1 and padj < 0.05, set as "DOWN"
df$diffexpressed[df$log2FoldChange < -1 & df$padj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"

ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=dflabel)) + 
  geom_point() + 
  theme_minimal() + scale_color_manual(values=c("salmon", "slategray2", "turquoise"))
```

![My Remote Image]([https://www.dropbox.com/s/.../my-remote-image.jpg?dl=0](https://www.dropbox.com/s/whst55bkk9hscoj/volcano_WB_18Jul22.png?dl=0))

--------------
### **Step 4: Pathway enrichment analysis**

To find the enriched pathways in each group of samples we used Gene Set Enrichment analysis with [Cluster profiler](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html) package in r.

```R
#Import packages 

library(clusterProfiler)
library(stats, warn.conflicts = FALSE)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# STEP 1 : set orgnaism for annotation

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE, warn.conflicts = FALSE)

#Check which options are available with the keytypes command
keytypes(org.Hs.eg.db)

# STEP 2 : Prepare the dataset 

# reading in data from deseq2
df = read.csv("DESeqResults_Jul12.2022_localFitType_WholeBlood.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

#STEP 3 : Convert gene IDs for gseKEGG function 

# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", 
          toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$ENTREZID = dedup_ids$ENTREZID


#STEP 4 : prepare the data for KEEG enrichment analysis 

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$ENTREZID

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


# STEP 5 : Create a KEEG object 

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


# STEP 6 : Visualize the data 

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)


# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, categorySize="padj", foldChange=gene_list)


ridgeplot(kk2) + labs(x = "enrichment distribution") + theme(axis.text = element_text(size = 8))
```
Write description/interpretation of the analysis output.


