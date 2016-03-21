Phylomark
=========
Phylomark - a tool to find phylogenetic markers from whole genome alignments***

contact: jasonsahl@gmail.com

#Changelog:

v_1_2: fixed potential problem with read orientation prior to fragment alignment
v_1_3: fixed the script structure, added additional scripts and README details
v_1_4: Huge overhaul. Blastall changed to blast+. Took away the need for using MUGSY
       alignments. Now scales better to hundreds to thousands of genomes. Command
       line arguments are completely different

-to install Phylomark, enter the directory and type:

(sudo) python setup.py install

-If you don't have sudo privalages, try:

python setup.py install --user

Then set your PYTHONPATH to include Phylomark

export PYTHOPATH=/Users/jsahl/Phylomark:$PYTHONPATH

To make this permanent, add this to your .bashrc or .profile

Dependencies

1. Biopython (www.biopython.org) #in .bashrc, point PYTHONPATH variable to Bio location
 (e.g. PYTHONPATH=/home/jsahl/biopython-1.53:$PYTHONPATH; export PYTHONPATH)
<<<<<<< HEAD
2. bx-python-tools (https://bitbucket.org/james_taylor/bx-python/wiki/Home) #add to PYTHONPATH
 (e.g. PYTHONPATH=/home/jsahl/bx-python-central/lib:$PYTHONPATH; export PYTHONPATH)
3. HashRF (http://code.google.com/p/hashrf/)
=======
 -blast+, can be obtained from NCBI
>>>>>>> development

the following scripts are included with Phylomark.  If you have an architecture different
than i86linux64, then you may need to re-compile on your system
<<<<<<< HEAD
1. FastTree (http://www.microbesonline.org/fasttree/)
2. mothur (http://www.mothur.org)
3. muscle (http://www.drive5.com/muscle/)

Phylomark requires 3 files to run correctly:

1. Reference genome in FASTA format, or set of genes in FASTA format
2. whole genome phylogeny in Newick format
3. PATH to directory of genome FASTA files

See the manual for detailed instructions
=======
-FastTree (http://www.microbesonline.org/fasttree/)
-mothur (http://www.mothur.org)
-muscle (http://www.drive5.com/muscle/)
-Dendropy (see the end of this document for license information)

Phylomark requires 3 arguments to run correctly:

1. whole genome phylogeny (can be generated with multiple methods)
2. directory of genomes that went into your phylogeny
3. Reference genome from one isolate from the whole genome alignment

-Now you want to alter the file, phylomark_env.sh, to set the Phylomark_DIR environment variable.
Then you can set the environment by:

source phylomark_env.sh

Once the files are generated and your environment is correct, Phylomark can be run by:

>phylomark.py -r <reference genome> -d <genome directory> -t <wga.tree>

Other parameters that can be changed include:
-s : step_size (integer).  The sliding window will move this many bases
-l : frag_length (integer).  Length of genomic fragments to include
--parallel_workers= (integer) : number of processors to use

Known issues:

-if the genome name is too long, muscle will truncate it.  Dendropy will then throw an error because
the names don't match compared to the whole genome phylogeny.

-If there are dashes in your genome names, they get changed to underscores somewhere in the pipeline.
While I'm trying to track down the source of the error, the current fix is to either remove dashes
or replace them with underscores

*An additional script Phylomark_v1_1_R.py is included to provide more detailed analysis about
nucleotide frequencies in each genomic fragment

-Two new dependencies are required for this script:

R (tested version = 2.14.1 - R v3 appears to not be supported at this time)
bioStrings (http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)

-The snps.r script must be in the same directory as your other input files
(script is modified from http://manuals.bioinformatics.ucr.edu/home/ht-seq)

-A new directory is created (R_output).  For each fragment, two files are created.  One file
is a table showing the base frequencies at each position in the alignment.  The second file is
a .pdf showing a cluster dendrogram, and a plot showing the nucleotide conservation across
the length of the fragment.  As this directory can fill up rapidly, I recommend that you use
a larger step size (e.g. 10) and a larger fragment size (e.g. 800).  Then I would only look
at the plots and tables for my best performing fragments. [Currently broken]

-Output
The two files that are of interest include:

seqs_shredded.txt (all of your potential markers)
results.txt (a list of your markers and the RF values)

Now you might be interested in combinations of markers that give you the lowest RF value.
For this you can do:

python find_marker_sets.py -d genomes/ -s tmp.fasta -t wga.tree

The "-d" flag points to a directory of your genomes in FASTA (*.fasta) format.  The "-t"
flag points to the WGA tree used in Phylomark.  You can set the number of markers to test
("-m" flag, default is 3), and the number of random iterations ("-i" flag, default is 20).
The "-s" flag points to a set of sequences. These should be filtered based on the Phylomark output.

The resulting file ("marker_results.txt") has the IDs of your markers and the resulting
RF value.  For example if you have three markers, your output would look like:

2018    2057    2523    32
864     867     2056    31

In this case, the markers 864, 867, 2056 produce the lowest RF values.  The sequence
for these markers are in your all_reads.fasta file
>>>>>>> development
