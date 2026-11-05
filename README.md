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
