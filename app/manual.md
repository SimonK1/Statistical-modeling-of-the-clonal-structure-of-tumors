# Accompanying manual

This manual belongs to the main part of this bachelor thesis and corresponds to the main.ipynb file.

This manual is preceded by copy calling steps, for which indexing of bam file is necessary. These two steps must be completed for the pipeline to work. Instructions for their preparation can be found in the files copyCalling.md and indexBam.md.

#### **However, for running the pipeline with sample data, these steps are not necessary, as all the necessary data is supplied!!!**

### User interface

The entire solution is concentrated within this jupyter notebook. To start, just start the whole jupyter notebook or start cell by cell.

### Importing libraries

Import all necessary libraries to run the pipeline. The rest of the imports are located in the main file, which is responsible for storing all the necessary transformers.

```python
import main
from sklearn.pipeline import Pipeline 
import pandas as pd
import cnvpytor
```



### Reading VCF file

The VCF file is loaded using a custom function that distributes the data to the assigned columns. The function is based on standard mandatory VCF columns with single sample. It creates a dataframe from these columns.

```python
main_vcf = main.load_vcf_file('./DO52567.vcf')
main_vcf.head()
```



### Reading Copy Caller results

Reading copy caller results follows the CNVpytor library columns, all of which are required. It creates a dataframe from these columns.

```python
cnvnator_df = pd.read_csv('./copyCaller/results/output.tsv',  names=['file_name', 'method', 'CNV Type', 'chr', 'CNV Region Start', 'CNV Region End', 'CNV size', 'CNV level', 'e-val1', 'e-val2', 'e-val3', 'e-val4', 'q0', 'pN', 'dG'], delimiter=r"\s+")
cnvnator_df.head()
```



### Pipeline

The entire main pipeline consists of 12 transformers that together prepare and run data for all tools except for running the basic PyClone. The operation of the transformers is explained below and with the help of detailed comments directly in the code.

- VcfDataExtractionTransformer()
  - Transformer that extracts all the necessary information from the supplied VCF file. Splits the last column according to the FORMAT column standard and converts the values to integer. It will also format the correct GENOTYPE according to the PyClone sample file.
- FilterQualityTransformer(*percentage*=90)
  - This transformer is responsible for filtering the percentage of the highest quality samples. By default, its input parameter is set to 90%, which filters out samples with a quality higher than 90%. This parameter can be changed as needed.
- CopyCallsMergeTransformer(cnvnator_df)
  - A transformer that combines the supplied vcf file with the file that is the result of the copy calling process. This transformer searches the CNV file and looks for whether the given read from the VCF file fits into one of the segments found by the copy caller. If it fits into any segment, it indicates the corresponding copy number. If the copy number does not match, it will be set to the default 2.
- PyCloneTransformer(*samples*=300)
  - This transformer is responsible for preparing data for the PyClone and PyClone-VI tools. The transformer receives the samples parameter, which represents the number of reads that it will randomly select from the VCF file. This parameter can be changed as needed, but by default it is set to 300, when we get correct results and the time requirement is not very high.
  - Transformer also prepares two files with the corresponding columns. One for PyClone and the other for PyClone-VI. These files are exported dataframes that contain the necessary columns extracted from the VCF file. Transfomer also changes the form of the GENOTYPE column and calculates the variant frequency with which PyClone gives more accurate results.
- SciCloneTransformer(cnvnator_df)
  - This transformer is responsible for preparing data for the SciClone tool. It prepares specifically two different files. The first one is a file containing calculated variant frequencies and the second one is a file containing copy number values. This transformer extracts this data directly from the VCF file.
- CoverageFileTransformer(cnvnator_df)
  - This transformer is responsible for creating the coverage file for the TitanCNA tool. The coverage file must be calculated based on the calculated segments from the CNV file. This transformer inserts segments that have a default copy number among the found segments and thus creates a coverage file. It also calculates logR values that TitanCNA requires instead of copy number values.
- TitanCNATransformer()
  - This transformer is responsible for creating another necessary file for the TitanCNA tool. Transformer extracts the necessary data directly from the VCF file and converts it into the necessary formats for TitanCNA. This file contains information about individual reads without copy number values. TitanCNA combines these two files and calculates the segments itself.
- FastCloneTransformer()
  - This transformer prepares similar data as the PyClone transformer, but with minor changes. This data is suitable for FastClone.

All transformers listed below run the respective tools using command line commands using the OS library. The TitanCNA and Sciclone tools specifically run R scripts that execute the necessary commands. These scripts are stored in the scripts folder.

- PyCloneVI()

- TitanCNA()

- SciClone()

- FastClone()

  

### Results

**PyClone**

PyClone brings a lot of results, but the one that is the best for the purposes of this bachelor's thesis is located in the plots folder and then in the cluster subfolder. Inside there is a file density.pdf which represents the necessary cellular prevalence. It shows all found clusters and their cellular prevalence value. We can also find in the picture which cluster is the largest and therefore the most important.

**PyClone-VI**

PyClone-VI does not provide any display method by itself, but the main pipeline will display the necessary density graph, which is in the corresponding folder. This image represents individual clusters, their cellular prevalence value and also the number of reads they contain. The higher the column, the more significant the cluster. Simply, similar representation as the older version of PyClone.

**SciClone**

The SciClone library for single sample produces only one main file type, which contains a density graph that indicates cellular prevalence based on variant allele frequency.

**TitanCNA**

The TitanCNA library provides a large number of display results. However, the folder generated by TitanCNA according to the name of the experiment contains all the clusters that TItan found according to the chromosomes. These charts show prevalence values for individual chromosomes. Each horizontal line represents one cluster. They are marked as Z1, Z2...

**FastClone**

Unfortunately, the FastClone library does not produce the desired results from the supplied data. For this reason, we cannot visualize or analyze the results in any way.