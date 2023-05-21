install.packages("doMC")
install.packages("ragg")
install.packages("remotes")
remotes::install_github("trevorld/r-optparse")
remotes::install_github('r-lib/ragg')

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
    
BiocManager::install("IRanges")
BiocManager::install("limma")
BiocManager::install("TitanCNA")


remotes::install_github("genome/bmm")
remotes::install_github("genome/sciClone")




