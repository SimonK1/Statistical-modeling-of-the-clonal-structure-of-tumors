# Loading required libraries
library(sciClone)
library(tidyverse)


# Read vaf data
v = read.table("./inputFiles/SciClone/ScicloneVafFile.tsv", header=T);
v1 = v[1:100,]

# Read copy number variants data
cn1 = read.table("./inputFiles/SciClone/ScicloneCopyFile.tsv")
names = c("Sample1")
sc = sciClone(vafs=v,
           copyNumberCalls=cn1,
           sampleNames=names[1])

# Create output
writeClusterTable(sc, "./results/SciClone/ResultSciClone")
sc.plot1d(sc,"./results/SciClone/ResultSciClone.1d.pdf")