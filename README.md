
# ViralRecall

*Note bene*, dear reader, a new and better version of this tool (v3.0) has been written by the intrepid Abdeali Jivaji and can be found here: https://github.com/abdealijivaji/ViralRecall_3.0

This new version is faster, easier to install, and has some new features that may be of interest. The 2.0 code will remain here but will no longer be maintained. 



ViralRecall is a flexible command-line tool for detecting signatures of
giant viruses (NCLDV) in genomic data. Version 2 has been updated to
focus more on NCLDV compared to version 1, but the original options are
still available.

### Dependencies

ViralRecall is written in Python 3.5.6 and requires biopython,
matplotlib, numpy, and pandas. ViralRecall uses Prodigal and HMMER3 for
protein prediction and HMM searches, respectively. Please ensure these
tools are installed in your PATH before using.

A requirements.yml file is provided in this repo to solve some of the
issues related to outdated pandas version. The requirements.yml file
specifies conda environment dependencies so you don't have to install
each separately. 
After cloning repository, please follow these steps:

> cd viralrecall

and 

> conda env create -f requirements.yml

If you wish to proceed without the requirements.yml file, simply create a conda 
environment by typing `conda create -n viralrecall`. In that case you might have
to install some dependencies yourself. On a Unix system you should be
able to install these tools with:

> sudo apt install prodigal

and

> sudo apt install hmmer

or if you don't have sudo privileges, you can try with conda:

> conda install prodigal -c bioconda

and

> conda install hmmer -c bioconda

### Installation and Database Download

Please ensure you are using \> Python 3.5.2 and have the appropriate
python modules installed. If this is an issue please create a Python
environment using conda (see here
<https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>)

To start: \> git clone <https://github.com/faylward/viralrecall>

> cd viralrecall

ViralRecall was tested on Ubuntu 16.04 and should work on most
Unix-based systems. To see the help menu use: \> python viralrecall.py
-h

Viralrecall can be run using either of two viral HMM databases: 1)
GVOGs, a custom set of Giant Virus Orthologous Groups that are fairly
specific to Nucleo-Cytoplasmic Large DNA Viruses (NCLDV), or 2) VOGDB,
which contains a wide collection of viral orthologous groups and is
useful for broad characterization of viral signatures.

In addition to either the GVOG or VOGDB searches, ViralRecall matches
proteins against the Pfam database (Pfam v. 32), which is a
broad-specificity database that detects many protein families that are
common in the genomes of cellular organisms.

The database files are available for download from Zenodo. To download
and unpack, navigate to the folder that contains the viralrecall.py
script and type:

> wget -O hmm.tar.gz
> <https://zenodo.org/records/12666277/files/hmm.tar.gz?download=1>

and then

> tar -xvzf hmm.tar.gz

This should create a hmm/ directory with the appropriate HMM files,
including the gvog.hmm database and the vogdb.hmm database (downloaded
from the vogdb.org website on 12/14/2020). This directory should be
located in the same directory as the acc/ directory and the
viralrecall.py script.

Note: The GVOG database was updated in May 2021, and it is recommended
to use this latest version. The new database is substantially smaller,
among other things, and this helps with runtime. But if you'd like to
access the original GVOG database you can do so here:
<https://zenodo.org/record/4710691/files/hmm.tar.gz?download=1>

After you have downloaded and unpacked the database files you should be
able to run viralrecall.py. To see the help menu you can run: \> python
viralrecall.py -h

### Basic Usage

To test if ViralRecall will run properly type: \> python viralrecall.py
-i examples/arm29B.fna -p test_outdir -t 2 -f

Results should be located in the test_outdir folder. The output folder
will contain:

\*.faa: The proteins predicted from the input file using Prodigal

\*.fasta: The DNA coding sequence of the predicted proteins from
Prodigal

\*.full_annot.tsv: A full annotation table of the predicted ORFs. This
includes descriptions of the GVOG and Pfam annotations, so it can be
useful if you want to look at certain annoatations in more depth.

\*.vregion.annot.tsv: An annotation of only the regions with some viral
signatures

\*summary.tsv: Summary statistics for the predicted viral regions (or
contig-level stats if the -c flag was used). This also includes the
NCLDV marker output (marker hit: bit score)

\*.pfamout: Raw output of the Pfam HMMER3 search

\*.vogout: Raw output of the GVOG or VOGDB HMMER3 search

\*.markerout: Raw output of the NCLDV marker gene HMMER3 search

Additionally, for each viral region viralrecall will print out .faa and
.fna files for the proteins and nucleotide sequences for the regions
found. Please be sure to use only .fna files as input.

### Options

There are several parameters you can change in viralrecall depending on
your preferences and the data you're analyzing. The important parameters
that will influence the results are:

**-db, --database** This is the database usef for viral detection. The
default is the GVOG database (also can be specifified with "GVOG").
"VOG" can be specified for the vogdb database, and "marker" to only
search against 10 NCLDV marker genes. GVOGs are more useful for
NCLDV-specific searches. The "marker" option is much faster and may be
useful for quickly screening large datasets.

**-s, --minscore** This is the mean score that a genomic regions needs
to have in order to pass the filter and get reported as a viral region.
The score is calculated from the HMMER3 scores, with higher scores
indicating more and better matches to the GVOG database, and lower
scores indicating more and higher matches to the Pfam database. The
default is 10.

**-w, --window** Size of the sliding window to use for calculating
moving averages. A smaller window may help predict short viral regions,
but may split large viral regions into several pieces.

**-m, --minsize** Minimum size, in kilobases, of the viral regions to
report.

**-g, --minvog** Minimum number of hits against the GVOG database that
must be recorded in a region in order for it to be reported (larger
values == higher confidence).

**-c, --contiglevel** If this option is used, mean ViralRecall scores
will be provided for the input contigs rather than viral regions. This
is useful for screening contigs for viral signatures.

**-r, --redo** If you have already run ViralRecall and you want to
re-run it with different parameters, you can use the -r flag to avoid
re-running Prodigal and HMMER, which are the most time-consuming steps.

**-b, --batch** Use this flag if the input is a folder of .fna files to
search, rather than a single .fna file.

For example, if we wanted to recover regions of a eukaryotic contig with
signatures of NCLDV, we could use the following command:

> python viralrecall.py -i examples/arm29B.fna -p test_outdir -s 15 -m
> 30 -g 10

Here we are asking for only regions that have a mean score \>= 15, are
at least 30 kilobases long, and have at least 10 GVOG hits.

If we want to quickly re-do the above analysis with different
parameters, but without re-doing gene predictions and HMMER3 searches,
we can use the -r flag:

> python viralrecall.py -i examples/arm29B.fna -p test_outdir -s 15 -m
> 15 -g 15 -w 20 -r

Maybe we want to re-do the analysis using a different e-value. The
default is 1e-10, which is fairly stringent, so we can relax it a bit:

python viralrecall.py -i examples/arm29B.fna -p test_outdir -s 15 -m 15
-g 15 -w 20 -r -e 1e-5

So once you finished the hmmer searches you can easily re-calculate
things with the -r flag.

### Batch mode

If you have many sequences you wish to test you can put them all in a
folder and use batch mode. Here the input (-i) should point to a folder
with .fna files in it. Basic usage is:

> python viralrecall.py -i examples/testfolder -p folderout -b

All of the output files should have their own folder in the folderout
directory. You can also use the -b flag with the -r flag for quick
re-calculations.

### Miscellaneous

Prodigal is intended to predict genes on prokaryotic genomes, and it
therefore will draw an error if used on very long eukaryotic contigs (\>
32 Mbp in length). I've re-compiled a binary of Prodigal that will run
on longer contigs, and this is available in the bin/ folder of this
GitHub repo. This is NOT the default version of prodigal that is used -
if you wish to use this binary you will need to make sure it is located
in your PATH (and not any other version of prodigal you may have
installed).

Some users have noticed errors or warnings involving Pandas, which uses
slightly different syntax in different Python versions. These can can
usually be resolved by changing the Python version used to 3.5.4.

### References

ViralRecall: A Flexible Command-Line Tool for the Detection of Giant
Virus Signatures in Omic Data, FO Aylward and M Moniruzzaman, Viruses,
2021; 13(2):150. <https://doi.org/10.3390/v13020150>.

For questions or comments feel free to email Frank Aylward at faylward
*at* vt dot edu

This tool requires Prodigal and HMMER3. Their citations are:

Hyatt et al. “Prodigal: prokaryotic gene recognition and translation
initiation site identification”. BMC bioinformatics, 2010.

Eddy, "A new generation of homology search tools based on probabilistic
inference". Genome Informatics, 2009.
