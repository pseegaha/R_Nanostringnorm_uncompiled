# Installation
## clone from git
```
git clone https://github.com/pseegaha/R_Nanostringnorm_uncompiled.git
```
## R install 
```
install.packages("R_Nanostringnorm_uncompiled/",repos=NULL,type="source")
```
## Load the required R packages
library(NanoStringNorm)

## Handling and analyzation of NanoString data
This package contains functions for parsing, handling, quality control and analyzation of data created on the
NanoString nCounter system in the form of RCC files.

## Reading File 
It is best to store all the RCC files (obtained from the NanoString's nCounter system) of interest in a separate directory, or at least with in a directory with no RCC files
from other experiments. The read.markup.RCC() function will then allow to transfer your files into a df object:
```
df <-read.markup.RCC(rcc.path = "path/to/RCC/directory/",rcc.pattern = "*.RCC|*.rcc",exclude = NULL,include = NULL,nprobes = -1)
```

The df object is a list containing the raw counts and the header from each RCC file. After normalization, other 
objects will be attached to it, e.g. normalization factors.

### Normalising Nanostring data
This is a quick R guide to learn about Nanostring technology (nCounter) and how to pre-process the data profiled on this platform.

## Description
The nCounter system from Nanostring Technologies is a direct, reliable and highly sensitive multiplexed measurement of nucleic acids (DNA and RNA) based on a novel digital barcode technology. It involves Custom Codeset of genes or off-the-shelf preassembled panels and on single cell (more details on NanoString website).

Each mRNA Expression CodeSet contains probes designed against fourteen ERCC transcript sequences.

– Six of these sequences are used as positive hybridization controls and eight are designed as negative controls.

– These positive controls and negative controls are present in each CodeSet independent of the sample. These help in normalising for any technical/systemic variability.

– In addition, the codesets can contain some housekeeping genes which can be used for normalising sample variability (biological normalisation) i.e. to correct for differences in sample input between assays. It is based on the assumption that the target sequences of the house keeping genes are consistent in their expression levels.

Note:  Read the nCounter guide available in the the link for more
details: (https://www.nanostring.com/application/files/1214/8942/4642/MAN-C0011-03_nCounter_Gene_Expression_Data_Analysis_Guidelines.pdf)

## Load the dataset
The data produced by the nCounter Digital Analyzer (nanostring) are exported as a Reporter Code Count (RCC) file which is a comma-separated text (.csv) file that contains the counts for each gene in a sample. Each cartridge has 12 lanes  i.e. 12 samples can be profiled on one nanostring cartridge.

For processing the data one can apply the normalization steps recommended by the company (using NanoStringNorm R package). Alternatively, the data can be treated as regular digital counts (RNA-seq) and can be analysed using edgeR TMM normalisation approach. However, in our experience former works better then the latter as it accounts for cross-hybridization related biases by allowing user to do background correction.

You can read the RCC files in two different ways i.e.use the excel import function read.xls.RCC to read directly from nCounter output files if provided in .xls format by the facility. However, do ensure that you are using the worksheet with the raw counts and not something that has been processed. An example dataset can be downloaded from GEO (GSE51488).

## Read the raw counts from the RCC excel spreadsheet output by the nCounter platform
```df <-read.xls.RCC("GSE51488_GAMA_Nanostring_RAW_Spleen_1.xls", sheet = 1)```
or,

you can use the following to process single sample markup RCC files (example:GSE95100) and merge the individual .RCC files together in one variable.

## Read the raw counts from individual RCC files from the directory (path of .RCC files )
```df <-read.markup.RCC(rcc.path = ".",rcc.pattern = "*.RCC|*.rcc",exclude = NULL,include = NULL,nprobes = -1)```
Pre-processing
Firstly, remove systemic biases by using geometric mean.

## Use geometric mean for technical normalisation
```all_samples_gm <- NanoStringNorm(x = df,anno = NA,CodeCount = 'geo.mean',Background = 'none',SampleContent = 'none', round.values = FALSE, take.log =FALSE,return.matrix.of.endogenous.probes =FALSE)```
Then, correct for cross-hybridization and normalise for sample variability by using background correction and house keeping genes respectively.

## Use housekeeping genes along with background correction(mean+2SD) for biological normalisation---#
```normalised_df <- NanoStringNorm(x = all_samples_gm,anno = NA,CodeCount = 'none',Background = 'mean.2sd',SampleContent = 'housekeeping.geo.mean', round.values = FALSE,is.log = FALSE, take.log = TRUE, return.matrix.of.endogenous.probes = TRUE )```
This returns the normalised values in log2 scale. If you want the data to be on linear scale then change take.log = FALSE

## Save the normalised data in a file---#
```write.table(normalised_df,"Normalised_data_nanostring.csv",sep=",",quote=F,row.names = T,col.names = T)```
The information about the R packages can be found below.

## Print the package versions used ---#
```sessionInfo()```

## Background Correction and Positive Control Normalization
The first step should be the calculation and subtraction of the background and positive control normalization, i.e. normalizing using the given spike-ins.
By default, the background is calculated using the mean of all background probes, and two standard deviations are added to it. You can adjust this using the "bm" and "sd.factor"
parameters. Similarly, you can choose which method should be used to calculate the positive control factors using the "pcm" parameter in the nsPositiveControlNormalization() function.
By default, the geometric mean is used.

```
df <- nsBackgroundCorrect(df)
df <- nsPositiveControlNormalization(df)
```


## Quality Control
The three basic plots that shold always be looked at are visualizing the FOV (Field Of View) counts, binding density and positive scaling factors.

```
plotBindingDensities(df)
plotFOV(df)
plotPositiveScalingFactors(df)
```

# Content Normalization

## Decision over normalization method
Currently, three different methods are implemented for content normalization: top100, housekeeping and total (aka global).
The primary choice of the normalizaton method should be dependent on the experimental design. Generally, the following should be considered:

- Housekeeping gene normalization is a good choice if the trust in the housekeepers is high
- Global normalization can be used if only a relatviely small fraction of genes is expected to be differentially expressed
- Top100 normalization is usefull if only few genes of the panel are expressed above threshold (miRNA data)

Some functions have been included to compare different normalization methods. To asses correlation of the housekeepers, you can use:

```
plotHousekeepingGenes(df)
housekeepingCorrelation(df)
```

The first function will generate a line chart, which allows to see whether all houesekeeping genes follow the same pattern over all genes. Additionally,
the houseekepingCorrelation function will show a correlation matrix for the housekeeping genes. Housekeeping genes should show high correlation. To compare the
effect of all three normalization methods, you can use the following function:

```
plotNormDistanceRatio(df, groups)
```

Where groups is a vector that assigns every sample to a specific group. The function then calculates the ratio of the mean inner-group distance to the mean global-distance
for every normalization and for the un-normalized counts. A smaller ratio after normalization implies better clustering of the groups. If groups are expected to be very different,
then this might be of interest.

## Normalization
After a content normalization method has been chosen, the counts can be normalized using the following function (houseekping normalization):
```
df <- nsNormalize(df, method="housekeeping")
```


## Downstream Analysis Functionality
After normalization, there exist several functions to visualize and explore the data. The following functions can be used for data visualization:

```
# Define a group vectors (examplary)
groups <- c("Group1", "Group1", "Group1", "Group1", "Group2", "Group2", "Group2", "Group2")
# Define a mean count that genes should have for consideration (default is 3 or 5 for most functions)
my_cutoff = 3

# Plot group-wise boxplot (log count values)
plotGroupBoxplot(df, groups, title="my_title", countCutoff=my_cutoff)

# Plot scaled Heatmap of expression values
plotHeatmapExpression(df, groups, countCutoff=my_cutoff)

# Plot a PCA
plotPCA(df, groups, countCutoff=my_cutoff)

```

To make differential expression analysis, a limma based analysis and the t-test are both implemented.
To be continued ...
