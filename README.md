# ViralRecall
ViralRecall is a flexible command-line tool for predicting prophage and other virus-like regions in genomic data.

# Dependencies
ViralRecall is written in Python 3.5.2 and requires numpy, pandas, and biopython. Additionally, matplotlib is required if you wish to use the plotting option. 
ViralRecall uses Prodigal and HMMER3 for protein prediction and HMM searches, respectively. Please ensure these tools are installed in your PATH before using. 
One a Unix system you should be able to install these tools with: 
> sudo apt install prodigal
> sudo apt install hmmer

or if you don't have sudo priveleges, you can try with conda:
>conda install prodigal -c bioconda
>conda install hmmer3 -c bioconda

## Installation
Please ensure you are using > Python 3.5.2 and have the appropriate python modules installed. If this is an issue please create a Python environment using conda (see here https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) 
ViralRecall was tested on Ubuntu 16.04 and should work on most Unix-based systems. To see the help menu use:
> python viralrecall.py -h

## Databases
Viralrecall uses two main HMM databases to analyze viral signatures in genomic data. The first is Pfam, which is a broad-specificity database that detects many protein families that are common in the genomes of cellular organisms. The Pfam database here has been modified to remove common viral protein families. The second HMM database is VOGDB, which is a comprehensive database of known viral protein families. The full VOG HMM database is quite large, and in some situations users may wish to use smaller sets to speed up runtime. We have created three VOG database that vary in size (all, large, and small). Users wishing to perform an exhaustive search should use the "all" database, while those preferring a faster search may wish to use the "small" database (the "large" database is an intermediate size option). 

The database files are available for download from the Virginia Tech library system. To download and unpack, navigate to the folder that contains the viralrecall.py script and type:
> wget XXX
> tar -xvfz XXX

## Usage
To test if ViralRecall will run properly type:
> python viralrecall.py -i examples/test_seq.fna -p test_outdir -t 2 -db small -f

Results should be located in the test_outdir folder. 
The output files contain:
*.faa:                   The proteins predicted from the input file using Prodigal
*.full_annot.tsv:        A full annotation table of the predicted ORFs. 
*.prophage.annot.tsv:    An annotation of only the viral regions (only present if some viral regions found)
*summary.tsv             Summary statistics for the predicted viral regions. 
*.pfamout                Raw output of the Pfam HMMER3 search
*.vogout                 Raw output of the VOGDB HMMER3 search

Additionally, for each viral region viralrecall will print out .faa and .fna files for the proteins and nucleotide sequences for the regions found. 


## Other Examples

There are several parameters you can change in viralrecall depending on your preferences 



## Citation


