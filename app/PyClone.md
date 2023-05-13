# PyClone installation

To install the pyclone library, it is necessary to have Python version 2.7 available, therefore the best idea is to create the necessary virtual environment using the conda tool. It is then possible to start the installation in the following way:

[Pyclone](https://github.com/Roth-Lab/pyclone) - Official documentation

```
conda install pyclone -c bioconda -c conda-forge
```



After installation and activation of the environment, it is necessary to run following command:

```
PyClone run_analysis_pipeline --in_files ./inputFiles/PyClone/PyCloneInput.tsv  --working_dir ./results/PyClone/
```



This operation is quite time consuming. For 300 samples, this operation takes 22 minutes on a MacBook with 16GBs of RAM.