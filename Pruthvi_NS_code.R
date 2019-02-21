#Please set it to your directory
#reading the directory
read.csv("/Users/pruthvi/Downloads/case_study_datascientist_data/case_study_annotations.csv")

install.packages("NanoStringNorm")
install.packages("heatmaply")
#installing the missing vsn
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("vsn", version = "3.8")

#revoking heatmaply plot
library(heatmaply)

# Load the required R packages
library(NanoStringNorm)

#please set it to your directory
# read the raw counts from individual RCC files from the directory (path of .RCC files )
df <-read.markup.RCC(rcc.path = "/Users/pruthvi/Downloads/case_study_datascientist_data",rcc.pattern = "*.RCC|*.rcc",exclude = NULL,include = NULL,nprobes = -1)

# use geometric mean for technical normalisation
all_samples_gm <- NanoStringNorm(x = df,anno = NA,CodeCount = 'geo.mean',Background = 'none',SampleContent = 'none', round.values = FALSE, take.log =FALSE,return.matrix.of.endogenous.probes =FALSE)

# use housekeeping genes along with background correction(mean+2SD) for biological normalisation---#
normalised_df <- NanoStringNorm(x = all_samples_gm,anno = NA,CodeCount = 'none',Background = 'mean.2sd',SampleContent = 'housekeeping.geo.mean', round.values = FALSE,is.log = FALSE, take.log = TRUE, return.matrix.of.endogenous.probes = TRUE )

# save the normalised data in a file, using this improvised command—# write.table with proper column number in the header
#please set it to your directory
write.table(normalised_df,"Normalised_data_nanostring.csv",sep="\t",quote=F,row.names = T,col.names = NA)


# print the package versions used ---#
sessionInfo()

#we’ll load the data using read.csv().
#please set it to your directory
df<- read.csv("/Users/pruthvi/Downloads/case_study_datascientist_data/Normalised_data_nanostring.csv", sep="\t")

#when we type just df and enter, it shows the table on R console because the df represents the csv file now. We are now transposing the columns to rows and vice versa because the task assigned by the assignment here is to see the heat map with positive and negative control genes in columns and samples in rows. For this purpose we are performing the transpose and saving the data as dv now#
dv=t(df)
dv

#edit the (0,0) string in the header of the column which is missing to a character before proceeding with the following steps
rownames(dv)[rownames(dv)=="X"]="Genes"

#checking the rows
rownames(dv)

# we are now saving the transposed table data i.e, dv to a different csv file name #
#Please set it to your directory
write.table(dv,"transposed.csv",sep=",",quote=F,row.names = T,col.names = F)

#we will load the transposed file now to the R console
#please set it to your directory
df<- read.csv("/Users/pruthvi/Downloads/case_study_datascientist_data/transposed.csv")

#susequent commands to checkout the string(0)
df
lapply(df, class)
table(df$class)
class(0)

#r command to check if all the elements in table are same or not
sapply(df,class)

#to check if the factor is true or not
sapply(df,is.factor)

#Is it NA?
is.na(df)

#which one is NA?
which (is.na(df))

#the above code will neglect the character of the string(0) and plot the map with the x and y legend names
rownames(df) <- df[, 1]
df <- df[, -1]

#install.packages("heatmaply")
#revoking heatmaply plot
#library(heatmaply)

#command for heat map
heatmaply(df,cexRow=0.5, cexCol=0.5,Rowv=FALSE,Colv=FALSE,scale="none",xlab = "Genes", ylab = "Sample_id", main="NanoString Normalised Data", margings=c(60,100,40,20))

#making a new data frame for the requested two genes
input <- df[,c('MCL1','CXCL1')]
trial=as.data.frame(input)

#boxplot command
boxplot(trial,xlab="MCL1",ylab="CXCL1",col=c("lightblue","lavender"))

#abline command for minimum i.e, 0 percentile
A=abline(h=quantile(trial$MCL1,0),col="violet",lty=2)
B=abline(h=quantile(trial$CXCL1,0),col="green",lty=2)

#abline command for 25 percentile
C=abline(h=quantile(trial$MCL1,0.25),col="Navy",lty=2)
D=abline(h=quantile(trial$CXCL1,0.25),col="red",lty=2)

#abline command for 50 percentile
E=abline(h=quantile(trial$MCL1,0.5),col="blue",lty=2)
F=abline(h=quantile(trial$CXCL1,0.5),col="brown",lty=2)

#abline command for 75 percentile
G=abline(h=quantile(trial$MCL1,0.75),col="yellow",lty=2)
H=abline(h=quantile(trial$CXCL1,0.75),col="pink",lty=2)

#abline command for maximum.i.e 100 percentile
I=abline(h=quantile(trial$MCL1,1),col="orange",lty=2)
J=abline(h=quantile(trial$CXCL1,1),col="maroon",lty=2)

#CXCL1 quantile 
a=quantile(trial$CXCL1,0)
b=quantile(trial$CXCL1,0.25)
c=quantile(trial$CXCL1,0.5)
d=quantile(trial$CXCL1,0.75)
j=quantile(trial$CXCL1,1)

 #CXCL1 row
CXCL1row<-data.frame(col1=c(a),col2=c(b),col3=c(c),col4=c(d),col5=c(j))

#MCL1 quantile
e=quantile(trial$MCL1,0)
f=quantile(trial$MCL1,0.25)
g=quantile(trial$MCL1,0.5)
h=quantile(trial$MCL1,0.75)
i=quantile(trial$MCL1,1)
#MCL1 row
MCL1row<-data.frame(col1=c(e),col2=c(f),col3=c(g),col4=c(h),col5=c(i))

#Creating Dataframe
Data_analysis_summary=rbind(MCL1row,CXCL1row)

#naming the columns
colnames(Data_analysis_summary)[1]=c("min")
colnames(Data_analysis_summary)[2]=c("25th_percentile")
colnames(Data_analysis_summary)[3]=c("50th_percentile")
colnames(Data_analysis_summary)[4]=c("75th_percentile")
colnames(Data_analysis_summary)[5]=c("max")
colnames(Data_analysis_summary)
Data_analysis_summary

#naming the rows
rownames(Data_analysis_summary)[1]=c("MCL1")
rownames(Data_analysis_summary)[2]=c("CXCL1")
Data_analysis_summary

#saving the data_analysis_summary##please set it to your system directory
write.table(Data_analysis_summary,"Summary_Stat.csv",sep="\t",row.name=T,col.name=NA)

