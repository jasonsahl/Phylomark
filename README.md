### Phylomark - a tool to find phylogenetic markers from whole genome alignments

contact: jasonsahl@gmail.com

### Changelog:

v_1_2: fixed potential problem with read orientation prior to fragment alignmenti

v_1_3: fixed the script structure, added additional scripts and README details

v_1_4: Huge overhaul. Blastall changed to blast+. Took away the need for using MUGSY alignments. Now scales better to hundreds to thousands of genomes. Command line arguments are completely different

v_1_5: Code rewritten for consistency with new programs as well as working with Python 3.7+

### Installation:

-The easiest way to install Phylomark is through Conda:

`conda create -n phylomark python=3.7`  
`conda activate phylomark`  
`conda install -c bioconda muscle fasttree blast biopython dendropy mothur`  

-Clone the repo from github:

`git clone https://github.com/jasonsahl/Phylomark.git`  
`cd Phylomark`  
`python setup.py install --user`  

Then set your PYTHONPATH to include Phylomark

`export PYTHOPATH=/Users/jsahl/Phylomark:$PYTHONPATH`

To make this permanent, add this to your .bashrc or .profile

### Running Phylomark:

Phylomark requires 3 arguments to run correctly:

1. whole genome phylogeny (can be generated with multiple methods)  
2. directory of genomes that went into your phylogeny  
3. Reference genome from one isolate from the whole genome phylogeny  

Phylomark can be run by:  

`phylomark.py -r reference_genome -d genome_directory -t wga.tree`  

Other parameters that can be changed include:
-s : step_size (integer).  The sliding window will move this many bases
-l : frag_length (integer).  Length of genomic fragments to include
--parallel_workers= (integer) : number of processors to use

### Known issues:

-if the genome name is too long, muscle will truncate it. Dendropy will then throw an error because
the names don't match compared to the whole genome phylogeny.

-If there are dashes in your genome names, they get changed to underscores somewhere in the pipeline.
While I'm trying to track down the source of the error, the current fix is to either remove dashes
or replace them with underscores

-Output
The two files that are of interest include:

1. query_sequences.fasta (all of your potential markers)  
2. results.txt (a list of your markers, the RF values, the eulcidian distance between the two trees, and the number of SNPs)  

One approach would be to find ~10 markers that have the lowest RF values and put them into a new FASTA (tmp.fasta in this example)

Now you might be interested in combinations of markers that give you the lowest RF value.  For this you can do:

`python find_marker_sets.py -d genomes/ -s tmp.fasta -t wga.tree`

The "-d" flag points to a directory of your genomes in FASTA (*.fasta) format. The "-t"
flag points to the WGA tree used in Phylomark. You can set the number of markers to test
("-m" flag, default is 3), and the number of random iterations ("-i" flag, default is 20).

The resulting file ("marker_results.txt") has the IDs of your markers and the resulting
RF value.  For example if you have three markers, your output would look like:

2018    2057    2523    32  
864     867     2056    31  

In this case, the markers 864, 867, 2056 produce the lowest RF values.  The sequence
for these markers are in your all_reads.fasta file
