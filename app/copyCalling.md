#  CNVpytor / CNVnator

To use the functions of the main pipeline, it is necessary to prepare a copy calling file, through a tool that finds all CNVs. We recommend using CNVpytor or CNVnator. Both libraries are suitable for preparing variants as they work identically and return the same formats. CNVpytor serves as an extension of the CNVnator library for compatibility with the Python programming language.



Before using this tool, it is necessary to have an indexed bam file in the same folder as regular bam file. If you do not have such file available, one must be created using the instructions in the file indexBam.md



Below we find all the commands needed to run and calculate CNVs. The process can be time consuming. For an 80GB BAM file, the process took about an hour (Calculated on a MacBook pro with 16GB RAM). The link to the official documentation for this tool is below:



[CNVpytor ](https://github.com/abyzovlab/CNVpytor)- Official documentation

[CNVpytor - Installer](https://anaconda.org/bioconda/cnvpytor) - Installation through conda

[CNVnator](https://github.com/abyzovlab/CNVnator) - Official documentation



```
cnvpytor -root file.pytor -rd Input.bam
cnvpytor -root file.pytor -his 1000 10000 100000
cnvpytor -root file.pytor -partition 1000 10000 100000
cnvpytor -root file.pytor -call 1000 10000 100000
```



After entering the view command, the interface of the view function appears. It is possible to enter different values to calculate different granularity of CNVs. We describe the specific benefits and disadvantages of granularity in more detail in the attached work. We include functions to process three different types of granularities:



## **Granularity 100000**

```
!cnvpytor -root file.pytor -view 100000
```



## Granularity 10000

```
cnvpytor -root file.pytor -view 10000
```



## Granularity 1000

```
cnvpytor -root file.pytor -view 1000
```



After entering one of the previous commands, the following commands must be entered in the view interface:

```
print calls
set Q0_range 0 0.5
set size_range 100000 inf
print calls
set p_range 0 0.00001
set print_filename output.xls
print calls
set print_filename output.vcf
print calls
```

