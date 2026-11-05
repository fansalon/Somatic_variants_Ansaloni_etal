# Long-read sequencing reveals recurrent and tissue-specific somatic variants in adult mouse tissues
Ansaloni F. and Simi A. *et al.*

This repository contains the code to reproduce the analyses performed in this manuscript - currently under revision.


## Step 1 -> get the data from Zenodo
Due to the GitHub file size limit, some of the files needed as input were deposited on Zenodo (https://zenodo.org/records/23989158). Get the data and place the files in ```zenodo/input_data/```.

```
# install this zenodo utility, in case not already installed
pip install zenodo-get

# get all the data associated to this repo
zenodo_get 23989158
```


## Step 2 -> run the code the reproduce the figures
In the ```figures``` folder, the code to re-generate each figure is deposited. Running the code via ```Rscript FigX.R``` should be enough to reproduce all the figures.



## Step 3 -> run the code the reproduce all the analyses
In the ```analyses``` folder, each analysis is contained in a numbered folder. Simply follow the order and run the code, if you'd like to reproduce it yourself. Some of the analysis are meant to be run, and were actually run, on a PC, some others were run on a cluster. There might be some reproducibility issues with some of the software versioning, but, at least in our hands, this does not change much.
