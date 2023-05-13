# Statistical modeling of the clonal structure of tumors





## Launch 

As a first step, it is necessary to make sure that docker is running.

To successfully launch the container, it is necessary to clone the repository and open it using the terminal. It is necessary to open the terminal and navigate to the directory in which the clone of this repository is located. The following command must be run inside the repository:

```
docker build -t jupyter . && docker run --rm -ti -p 8888:8888 jupyter
```

This command prepares the entire necessary installation and subsequently starts jupyter lab. All the installation commands are in the **dockerfile**.

This process can take a considerable amount of time. More than 100 R dependencies must be installed to run the R libraries used by this container. These libraries must be downloaded and built. The entire process of initially starting the container takes about half an hour. After successful building of the container, the terminal displays a message that must be copied into the browser window. Subsequently,  jupyter lab will open.



### Main pipeline

The main pipeline is located in the **main.ipynb** file and the corresponding documentation is located in the **manual.md** file.



### Indexing bam files

If the supplied data has not to been indexed, it is necessary to follow the instructions in the file **IndexBam.md**.



### Copy Calling 

If the data does not contain CNVs, it is necessary to follow the instructions in the file **CopyCalling.md**

