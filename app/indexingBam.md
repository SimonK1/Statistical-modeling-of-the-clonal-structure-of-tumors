# **Indexing Bam File**

For indexing, it is necessary to use the Samtools tool, which can index even large BAM files. Below is the command where you have to replace the paths to the BAM files needed for indexing. The official documentation for this feature can be found here:

[Index function](*http://www.htslib.org/doc/samtools-index.html*)

```
samtools index ./BamSource/Input.bam  ./BamSource/Input.bam.bai
```



(Optional) Below is a command that serves to divide the bam file according to chromosomes. This command is useful only if the computer on which we perform the indexing cannot process the entire file (all chromosomes at once). This command specifically selects the first chromosome and saves it to the file chr1bam.bam. The official documentation for this feature can be found here:

[View Function](*http://www.htslib.org/doc/samtools-view.html*)

```
samtools view -bo chr1.bam -s 123.4 ./BamSource/Input.bam 1
```

