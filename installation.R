install.packages("remotes")

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("TitanCNA")

BiocManager::install("IRanges")
remotes::install_github("genome/bmm")
remotes::install_github("genome/sciClone")