SpikeRBDStabilization README File

Last Updated: October 14, 2021

## CONTENTS: ##
1. About
2. Detailed Instructions
3. Running the Code
4. Output
5. References


### 1. ABOUT: ###
- - - -
This software identifies potential escape mutants for a given neutralizing antibody (nAb) that binds to the Spike receptor binding domain (S RBD) on SARS-CoV-2. The code was initially used in [Francino-Urdaniz et al. (2021)](https://doi.org/10.1016/j.celrep.2021.109627). To accompany this publication, a protocol paper was written with explicit instructions for completing experiments and running this software correctly ([Haas et al., 2021](https://doi.org/10.1016/j.xpro.2021.100869)).

The software includes two separate modules: deep mutational scanning (abbreviated 'dms') and 'analysis'. They are run separately, first the dms module and then the analysis module. The dms module is responsible for compiling all data from deep sequencing FASTQ files and calculating basic metrics. The analysis section calculates statistics, determines escape mutants and generates the final CSV files and Microsoft Excel heatmaps.


### 2. DETAILED INSTRUCTIONS: ###
- - - -
For a comprehensive description and instructions on the setup, inputs, and optional arguments for this software, please refer to the protocol ([Haas et al., 2021](https://doi.org/10.1016/j.xpro.2021.100869)).


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

More details on output files and a complete explanation of the analysis can be found in [Haas et al., 2021](https://doi.org/10.1016/j.xpro.2021.100869).


### 5. REFERENCES: ###
- - - -
Francino-Urdaniz, I. M., Steiner, P. J., Kirby, M. B., Zhao, F., Haas, C. M., Barman, S., Rhodes, E. R., Leonard, A. C., Peng, L., Sprenger, K. G., Jardine, J. G., & Whitehead, T. A. (2021). One-shot identification of SARS-CoV-2 S RBD escape mutants using yeast screening. *Cell Reports (Cambridge)*, 36(9), 109627-109627. https://doi.org/10.1016/j.celrep.2021.109627

Haas, C. M., Francino-Urdaniz, I. M., Steiner, P. J., & Whitehead, T. A. (2021). Identification of SARS-CoV-2 S RBD escape mutants using yeast screening and deep mutational scanning. *STAR Protocols*, 2(4), 100869-100869. https://doi.org/10.1016/j.xpro.2021.100869
