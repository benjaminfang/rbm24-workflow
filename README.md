# RNA-seq and scRNA-seq workflow involed in rbm24 study

This repository contains codes and information to re-produce the RNA-seq and scRNA-seq analyses in our published paper.

Any inquiry or assitance needs please contact me without hizatation at benjaminfang.ol@outlook.com

## Requirements

- Python >= 3.11

- biobrary == 0.1.4

    `pip install biobrary==0.1.4`

- annoread == 0.1

    `pip install annoread == 0.1`

- matplotlib

- numpy

- cutadapt == 4.6

- bowtie2 == 2.5.2

- samtools == 1.18

- DESeq2 == 1.44.0

- cellranger == 7.2.0

- trimmomatic == 0.39

- scanpy == 1.9.6

- DecontX

    Included in singleCellTK(the version 2.12.1 was used in our workflow).

- Scrubleti == 0.2.3

## Installation

All requirements need be installed. The version number list in **Requirements** is 
the software version in our workflow, other version may works too.

After the installation of requirements, The bash script can run directly.

## Usage

### RNA-seq

1. Download the RNA-seq data from NCBI(accession ID: XXX), and put it under RNA-seq/raw-data/.

2. Change directory to RNA-seq/scripts, and edit bash scripts in the directory, changing *working_dire* to
your actural working path.

3. The bash scripts should be executed in numerical order.


### scRNA-seq

1. Download the scRNA-seq data from NCBI(accession ID: XXX), and put it under scRNA-seq/raw-data/6hpf.

2. Change directory to scRNA-seq/scripts/6hpf, and changing *working_dire* to your actural working path.

3. Run the bash scripts in numerical order.

## Citations

Yizhuang Zhang. et al. Rbm24 is an organizer component of germ plasm to determine germ cell fate. (under review)
