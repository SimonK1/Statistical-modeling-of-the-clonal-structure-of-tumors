# Statistical modeling of the clonal structure of tumors

This repository serves as an accompanying repository for submitting a bachelor's thesis with the topic "Statistical modeling of the clonal structure of tumors". Due to the large volume of input genomic data, it was necessary to create a repository from which they can be easily accessed. The bachelor's thesis and the entire repository were developed by a FIIT STU (Faculty of Informatics and Information Technologies) student named Simon Kokavec.

## Launch 

As a first step, it is necessary to make sure that docker is running.

To successfully launch the container, it is necessary to clone the repository and open it using the terminal. To be specific, it is necessary to open the terminal and navigate to the directory in which the clone of this repository is located. The following command must be run inside the repository:

```
docker build -t jupyter . && docker run --rm -ti -p 8888:8888 jupyter
```

This command prepares the entire necessary installation and subsequently starts jupyter lab. All the installation commands are in the **dockerfile**.

This process can take a considerable amount of time. More than 100 R dependencies must be installed to run the R libraries used by this container. These libraries must be downloaded and built. The entire process of initially starting the container takes about half an hour. After successful building of the container, the terminal displays a message that must be copied into the browser window. Subsequently,  jupyter lab will open.

The working directory is set to a directory called **app** where all the necessary files are located. Jupyter lab will open in this directory.



### Main pipeline

The main pipeline is located in the **main.ipynb** file and the corresponding documentation is located in the **manual.md** file. 

Pipeline prepares all the files needed to run all libraries except the PyClone library. This library is written for a very old version of python 2.7, so neither the container nor the code can be used with this library. Therefore, the pipeline only prepares input files for this tool. To start the PyClone library, it is necessary to follow the instructions in the **PyClone.md** file.



### Indexing bam files

If the supplied data has not to been indexed, it is necessary to follow the instructions in the file **indexBam.md**. To start the main pipeline, this step is not necessary, as the files are already supplied.



### Copy Calling 

If the data does not contain CNVs, it is necessary to follow the instructions in the file **copyCalling.md**. To start the main pipeline, this step is not necessary, as the files are already supplied.



### File structure

##### Directories

- app - Directory in which the entire functioning of the program is located. There are both input and output files. Docker container sets this directory as the **working directory**.

##### Files

- Dockerfile - File needed to initialize the docker container
- reqs.txt - This file serves as a requirements file for the installation of the docker container and all python dependencies.
- installation.R - This file server as an installation file for all R packages and dependencies.

### App directory structure

App directory structure is described below:

##### Directories

- copyCaller - The directory where the input and output files for the copy caller are stored. Inside there is subdirectory called **inputFile** where you need to insert the input file for CNVpytor/CNVnator derived from BAM file (We do not attach this file as it reaches large sizes). After successful run of copy caller, the resulting .tsv file must be saved in a subdirectory called **results**. 
- inputFiles - Directory in which all input files for used tools are stored. This folder and all its subfolders are used by the main pipeline for storing data preparation results. The internal subdirectories are divided nominally according to the names of the tools.
- results - The directory used by the main pipeline to store the final results. The internal subdirectories are divided nominally according to the names of the tools.
- scripts - All R scripts used by the pipeline are stored in this directory. Do not change this directory.

##### Files

- main.ipynb  - It is the main file from which the main pipeline is launched.
- manual.md - This file serves as a manual for the main pipeline.
- DO52567.vcf - Input VCF file. (Added through GIT LFS. Visible as symbolic link until cloned)
- copyCalling.md - This file serves as a manual for the optional use of copy caller.
- indexBam.md - This file server as a manual for the optional indexing of bam file.
- preprocessing.py - This file contains all the important transformers that the main pipeline uses.
