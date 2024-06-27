# Nextflow pipeline to perform AlphaFold-derived template modelling

This repository contains the scripts to reproduce the results published in ...

## Graphical diagram

![Graphical scheme](figures/figureS2.png)

## Dependencies

The pipeline is designed to run on an high perfromance computing cluster through [LSF](https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=management-job-submission) via [NEXTFLOW](https://www.nextflow.io/) (version 23.04.1).

### ColabFold database
Download ColabFold databases as explained [here](https://github.com/sokrypton/ColabFold) and place them in [assets](assets/)

### PfamxAlphaFold
The alphafold-pfam database should be created using [these instructions](https://github.com/matthiasblum/pfam-alphafold) and placed in [assets](assets/). 

### Pfam
The Pfam 35 [MySQL](https://www.mysql.com/) database can be downloaded [here](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/) and and each file of the database imported as follows:
```
zcat file.sql.gz | mysql -u <user> -p pfam_35
```
Once the pfam database has been created, pass the hostname in "prams.pfamhost", the username in "params.pfamuser",  the password in "params.pfampassword" and the port in "params.pfamport" for the database in the [pipeline.nf](./pipeline.nf) file header.

### Change script permissions

Allow execution of scripts by nextflow with:
```
chmod +x bin/*
```

## Run the pipeline

To reproduce the results obtained in the paper, run the following:

```
bsub < sh/run.sh
```

To run the pipeline on a customised set of proteins, use [this script](sh/run_custom_proteins.sh) as example. Note that you will need to specify the proteins and the domain to be used for each protein.

## Other scripts 

The scripts contained in [scripts](scripts/) were also adopted:

- estimate_co2.sh was used to estimate the amount of CO2 produced with the computation;
- get_distribution.py was used to extract the whole pLDDT pfam distributions;
- get_domain_info.sh was used to extract the information about domains considered in the paper;

The [images](images) can be reproduced using [this notebook](generate_figures.ipynb) after downloading the results ...
Image figureS3 was obtained with [draw.io](https://www.drawio.com/).
Dependencies needed: python 3.8, seaborn (version 12.2), pandas (version 1.5.3), matplotlib (version 3.6.2), numpy (version 1.24.3), Biopython (version 1.81), scikit-learn (version 1.2.2), 


[The diagram](figures/AF2Fix_diagram.png) was generated with https://app.diagrams.net
