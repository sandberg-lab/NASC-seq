# NASC-seq analysis pipeline

## System Requirements

No custom hardware was used, but the use of a computing cluster is recommended.
Running time for the entire pipeline for a single cell on a desktop computer is estimated at 24 hours.
To facilitate easy implementation and rapid processing we suggest to move the entire analysis pipeline to an EC2 instance from Amazon Web Services (AWS). We have preloaded an Amazon Machine Image (AMI) with the code, dependencies, hg38, and sample data (AMI-XXXXXXXX).

## Software Dependencies

The analysis pipeline has been tested with the following supporting software.

```
R 3.5.1
STAR (2.5.4b)
Samtools (1.9)
Python2.7 2.7.14
Python3 3.6.5
Java 1.8
Picard (MarkDuplicates 2.18.17-SNAPSHOT)
Cutadapt 1.18
fastQC 0.11.8
python3.6 devtools
PyStan 2.18.0.0 (python module)
Pandas 0.23.4 (python module)
Pysam 0.15.1 (python module)
joblib 0.13.0 (python module)
scipy xx.xx.xx (python module)
ggplot2 3.1.0 (R-package)
RsubRead 1.32.1 (R-package)
tidyr 0.8.2 (R-package)
dplyr 0.7.8 (R-package)
cowplot 0.9.3 (R-package)
```

## Usage

  -h, --help                        show this help message and exit
  
  -e, --experimentdir [dir]         experiment directory
  
  -p, --numCPU [integer]            number of CPUs to use
  
  -f, --flag [flag]                 flag specifying the step in the NASC-seq pipeline that will be performed:

    trim                            Trim fastq files using trimgalore
    
    align                           Align fastq files using STAR to hg38
    
    removegenome                    Removes genome from shared memory
    
    removeduplicates                Removes duplicates from aligned bam files using picard
    
    annotate                        Annotate features in bam files using Rsubread and index files
    
    conversiontag                   Tag conversions in bam file headers
    
    vcfFilter                       Select possible SNPs from shared mismatches between cells
    
    tagFilter                       Remove SNPs from conversion tags and index files
    
    cellQC                          Perform basic QC visualization to decide on QC cutoffs
    
    cellFilter                      Filter bam files based on QC cutoffs in the config.py file
    
    calculatePE                     Calculate error probability
    
    prepareData                     Prepares pickles from data for scalable processing using AWS or similar

## Installation

The git repository can be cloned directly (<1 minute). 

## Directory Structure

Raw fastq files (paired-end, 2x 150 cycles) should be moved to the following directories:
    
```
    Experimentdir
        fastqFiles
            rawdata
                P1_A01_S1 (Plate1, Row A, column 01)
                    ...R1.fastq.gz
                    ...R2.fastq.gz
                P1_A02_S2
                    ...R1.fastq.gz
                    ...R2.fastq.gz
                ...
```
## Performing NASC-seq analysis

A configuration file should be prepared in the root of the experiment directory. For the layout of the configuration file see the example config file (/NASC-seq/data/config_example.py). In addition to the locations of some of the dependencies, users will have to refer to a genome using STAR (gnv), a gtf file (gtf) and an SJBD file (sjdb) to facilitate memory sharing while aligning. The config file furthermore includes links to a file with strand information for all features in the genome (strandFile), a stan model file (stanFile).


## Output

While running the pipeline a number of fastq files as well as bam files will be produced and saved in the fastqFiles and bamFiles folders respectively.

The folders will contain the following partially processed datafiles:
    
```
Experimentdir
    bamFiles
        aligned_bam             STAR output
        duplRemoved_bam         PICARD duplicate removal output
        annotated_bam           Rsubread annotated output
        annotated_sorted_bam    Sorted Rsubread annotated output
        tagged_bam              Conversion-tagged output
        filteredTagged_bam      SNP-filtered conversion-tagged output
```
    
## Example Data

### Scalable data processing

To facilitate rapid processing and easy implementation we suggest to run the entire pipeline on Amazon Web Services
