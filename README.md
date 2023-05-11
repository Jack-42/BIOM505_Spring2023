# BIOM505 Spring 2023
### Code and final data for BIOM505 Spring 2023. 
### By Jack Ringer

This work investigates how predictions from ESMFold compare to AlphaFold in predicting protein stability change for mutant proteins using ACDC-NN and DDGun. For additional information please see the project report (docs/BIOM505_FinalReport.pdf). For questions please post an issue.

# Code
The **code/** directory contains the all of the code used to arrive at the findings described in the final report.
* The notebook file (Protein_Stability_Prediction.ipynb) contains code for data cleaning, visualization, analyzing results, and miscellaneous intermediary tasks. Note that not all of the cells will be runnable using only the data provided in this repo.
* The subdirectory **code/carc_scripts** contains scripts used to generate protein structures using ESMFold, protein profiles using HH-Suite/HH-Blits, and ddG predictions with DDGun and ACDC-NN. These scripts were used on the Xena machine at CARC and are not intended for use on local machines.
* Note: The script esmfold_inference.py was not written by me and is only slightly modified from code from the ESM repo: https://github.com/facebookresearch/esm

# Data
The **data/** directory contains the final ddG prediction files as well as a few intermediary files. Note that this repo **does not** include all of the data used in this project due to size limitations (e.g., structure pdbs and protein profiles are not included). However, these files can be regenerated using the code contained in this repo. 

Note that the ddG experiment files can be downloaded from the ThermoMut (https://biosig.lab.uq.edu.au/thermomutdb/) and FireProt (https://loschmidt.chemi.muni.cz/fireprotdb/) websites.

# Setup
Code contained in the notebook ought to be self-contained (assuming you have the relevant data). 

The sub-sections below describe how each method can be setup for use at CARC. It is recommended to use a separate conda environment for each method when installing packages.
### ESMFold
The document docs\esmfold_setup_carc.pdf describes how one can setup ESMFold on CARC to generate structure predictions. Then the scripts contained in the **esmfold_scripts/** subdirectory can be used to generate predictions using ESMFold.

### HH-Suite
The user guide for HH-Suite provides documentation on how HH-Suite can be installed: https://github.com/soedinglab/hh-suite/wiki

### DDGun
Clone the DDGun repository (https://github.com/biofold/ddgun) and make sure the DDGUN_SCRIPT_PATH variable is set appropriately. 

### ACDC-NN
Follow instructions from the ACDC-NN repo: https://github.com/compbiomed-unito/acdc-nn

# Acknowledgement 
Much thanks to CARC (http://carc.unm.edu/) for allowing us access to computational resources throughout the course of this project. As well, thanks to the authors of the repositories and websites listed above for their work and making their projects publicly-available!
