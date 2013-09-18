Phylomark
=========
***Phylomark - a tool to find phylogenetic markers from whole genome alignments***
*contact: jasonsahl@gmail.com

#Changelog:

v_1_2: fixed potential problem with read orientation prior to fragment alignment
v_1_3: fixed the script structure, added additional scripts and README details

#to install Phylomark, enter the directory and type:

python setup.py install

Then set your PYTHONPATH to include Phylomark

export PYTHOPATH=/Users/jsahl/Phylomark:$PYTHONPATH

To make this permanent, add this to your .bashrc or .profile

#Dependencies

-Biopython (www.biopython.org) #in .bashrc, point PYTHONPATH variable to Bio location
 (e.g. PYTHONPATH=/home/jsahl/biopython-1.53:$PYTHONPATH; export PYTHONPATH)
-bx-python-tools (https://bitbucket.org/james_taylor/bx-python/wiki/Home) #add to PYTHONPATH
 (e.g. PYTHONPATH=/home/jsahl/bx-python-central/lib:$PYTHONPATH; export PYTHONPATH)
-HashRF (http://code.google.com/p/hashrf/)

#the following scripts are included with Phylomark.  If you have an architecture different
than i86linux64, then you may need to re-compile on your system
-FastTree (http://www.microbesonline.org/fasttree/)
-mothur (http://www.mothur.org)
-muscle (http://www.drive5.com/muscle/)

Phylomark requires 5 files to run correctly:

1. concatenated alignment from the whole genome maf file
2. whole genome phylogeny
3. input mask from mothur showing polymorphic positions
4. combined multi-fastA of all genomes that went into the whole genome alignment
5. Reference genome from one isolate from the whole genome alignment

Files 1-4 can be created with the Phylomark_prey.py script included.  All you need to have is
the input MAF file, and a directory of genomes that went into the alignment.  Examples
of these files for E. coli are included on SourceForge

tools/Phylomark_prep.py --input-maf=your.maf --fasta-dir=fasta_dir

#Now you want to alter the file, phylomark_env.sh, to set the Phylomark_DIR environment variable.
Then you can set the environment by:

source phylomark_env.sh

Once the files are generated and your environment is correct, Phylomark can be run by:

Phylomark.py -a <concatenated_alignment> -m <mothur_mask> -t <wga.tree> -r <reference_genome>
-c <combined_multi_fasta> 

Other parameters that can be changed include:
-s : step_size (integer).  The sliding window will move this many bases
-l : frag_length (integer).  Length of genomic fragments to include
-k : keep_length (integer).  Keep fragments if they contain this many polymorphisms
--parallel_workers= (integer) : number of processors to use

Known issues:

-if the genome name is too long, muscle will truncate it.  HashRF will then throw an error because
the names don't match compared to the whole genome phylogeny. 

***An additional script Phylomark_v1_1_R.py is included to provide more detailed analysis about
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
at the plots and tables for my best performing fragments.

###Output
The two files that are of interest include:

all_reads.fasta (all of your potential markers)
results.txt (a list of your markers and the RF values)

To select markers that fall below your chosen RF threshold (30 in this case), you can to:

python tools/filter_phylomark_output.py -i all_reads.fasta -d results.txt -o tmp.fasta -r 30

Now you might be interested in combinations of markers that give you the lowest RF value.
For this you can do:

python find_marker_sets.py -d genomes/ -s tmp.fasta -t wga.tree 

The "-d" flag points to a directory of your genomes in FASTA (*.fas) format.  The "-t"
flag points to the WGA tree used in Phylomark.  You can set the number of markers to test
("-m" flag, default is 3), and the number of random iterations ("-i" flag, default is 20).

The resulting file ("marker_results.txt") has the IDs of your markers and the resulting
RF value.  For example if you have three markers, your output would look like:

2018    2057    2523    32
864     867     2056    31

In this case, the markers 864, 867, 2056 produce the lowest RF values.  The sequence
for these markers are in your all_reads.fasta file
