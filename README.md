Phylomark
=========
Phylomark - a tool to find phylogenetic markers from whole genome alignments***

contact: jasonsahl@gmail.com

#Changelog:

v_1_2: fixed potential problem with read orientation prior to fragment alignment
v_1_3: fixed the script structure, added additional scripts and README details

-to install Phylomark, enter the directory and type:

python setup.py install

Then set your PYTHONPATH to include Phylomark

export PYTHOPATH=/Users/jsahl/Phylomark:$PYTHONPATH

To make this permanent, add this to your .bashrc or .profile

Dependencies

1. Biopython (www.biopython.org) #in .bashrc, point PYTHONPATH variable to Bio location
 (e.g. PYTHONPATH=/home/jsahl/biopython-1.53:$PYTHONPATH; export PYTHONPATH)
2. bx-python-tools (https://bitbucket.org/james_taylor/bx-python/wiki/Home) #add to PYTHONPATH
 (e.g. PYTHONPATH=/home/jsahl/bx-python-central/lib:$PYTHONPATH; export PYTHONPATH)
3. HashRF (http://code.google.com/p/hashrf/)

the following scripts are included with Phylomark.  If you have an architecture different
than i86linux64, then you may need to re-compile on your system
1. FastTree (http://www.microbesonline.org/fasttree/)
2. mothur (http://www.mothur.org)
3. muscle (http://www.drive5.com/muscle/)

Phylomark requires 3 files to run correctly:

1. Reference genome in FASTA format, or set of genes in FASTA format
2. whole genome phylogeny in Newick format
3. PATH to directory of genome FASTA files

See the manual for detailed instructions
