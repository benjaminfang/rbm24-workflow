# The RNA-seq and scRNA-seq workflow involved in the rbm24 study

This repository contains code and information to reproduce
the RNA-seq and scRNA-seq analyses in our *paper* currently
under submission.

For any inquiries or assistance needed, please don't hesitate
to contact us at shaoming@sdu.edu.cn.

## Requirements

- Python >= 3.11

- biobrary == 0.1.5

    `pip install biobrary==0.1.5`

- annoread == 0.1.5

    `pip install annoread==0.1.5`

- matplotlib

- numpy

- cellranger == 7.2.0

- trimmomatic == 0.39

- scanpy == 1.9.6

- DecontX

    Included in singleCellTK(the version 2.12.1 was used in our workflow).

- Scrublet == 0.2.3

## Installation

All requirements need to be installed. The version number listed in
**Requirements** is the software version used in our workflow,
but other versions may work as well.

After the installation of the requirements, the bash script
can be executable directly.

## Usage

### scRNA-seq

1. Download the scRNA-seq data from NCBI(accession ID: XXX), and place it under scRNA-seq/raw-data/6hpf.

2. Change directory to scRNA-seq/scripts/6hpf. Edit bash scripts and change *working_dire* to your actual working path.

3. Run the bash scripts in numerical order.

## Citations

Yizhuang Zhang. et al. Rbm24 is an organizer component of germ plasm to determine germ cell fate. (under submission)
