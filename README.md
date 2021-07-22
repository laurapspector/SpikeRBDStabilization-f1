SpikeRBDStabilization README File

Last Updated: July 21 2021\
Current Version: 1.0

## CONTENTS: ##
1. About
2. Detailed Instructions
3. Running the Code
4. Output
5. References


### 1. ABOUT: ###
- - - -
This software identifies potential escape mutants for a given neutralizing antibody (nAb) that binds to the Spike receptor binding domain (S RBD) on SARS-CoV-2. The code was initially used in Francino-Urdaniz et al. (2021). To accompany this publication, a protocol paper was written with explicit instructions for completing experiments and running this software correctly **(link to be added upon publication)**.

The software includes two separate modules: deep mutational scanning (abbreviated 'dms') and 'analysis'. They are run separately, first the dms module and then the analysis module. The dms module is responsible for compiling all data from deep sequencing FASTQ files and calculating basic metrics. The analysis section calculates statistics, determines escape mutants and generates the final CSV files and Microsoft Excel heatmaps.


### 2. DETAILED INSTRUCTIONS: ###
- - - -
For a comprehensive description and instructions on the setup, inputs, and optional arguments for this software, please refer to the protocol **(link to be added upon publication)**.


### 3. RUNNING THE CODE: ###
- - - -
Once the set up and configuration file has been completed, the software can be run easily from the command line:
1. python3 -m dms --config configuration_file_name [other arguments]
2. python3 -m analysis --config configuration_file_name [other arguments]

Note: For each new antibody, the [Analysis] section of the configuration file will need to be modified.


### 4. OUTPUT: ###
- - - -
This program will create two folders in the output directory: 'Output' (from dms module) and 'Processed' (from analysis module). The default output directory is the root directory of the software. The quantity of output files from this software depends entirely on the number of experiments, samples, and proteins defined in the configuration file.

The Output folder contains:
1. A CSV file for each nAb that displays each variant, the reference and selected counts and enrichment ratios.
2. A CSV file for the control population that displays reference and selected counts as well as enrichment ratios.

The Processed folder contains:
1. A CSV file for each nAb that now includes a critical enrichment ratio, false discovery rate and a p-value for the significance of each enrichment ratio.
2. A Microsoft Excel file for each nAb that shows a heatmap of the variants with any mutants classified as escape mutant hits are highlighted in blue.

More details on output files and a complete explanation of the analysis can be found in the protocol **(link to be added upon publication)**.


### 5. REFERENCES: ###
- - - -
Francino-Urdaniz, I. et al. (2021) ‘One-shot identification of SARS-CoV-2 S RBD escape mutants using yeast screening’, bioRxiv, 1(303), pp. 1–12.
