#!/usr/bin/env python

import sys
import optparse
import errno
import subprocess
from optparse import OptionParser
try:
    from phylomark.util import *
except:
    print("your phylomark environment is not set. Add the phylomark directory to your PYTHONPATH")
    sys.exit()

def is_tool(name):
    from distutils.spawn import find_executable
    return find_executable(name) is not None

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s cannot be opened' % option)
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("directory of fastas cannot be found")
        sys.exit()

def test_options(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("option not supported.  Only select from T and F")
        sys.exit()

def main(ref,genes,genomes,tree,step_size,frag_length,parallel_workers):
    dependencies = ['mothur','muscle','FastTree','blastn','makeblastdb']
    if "NULL" not in ref and "NULL" not in genes:
        logging.logPrint("You can't select both genes and a reference sequence..exiting")
        sys.exit()
    elif "NULL" in ref and "NULL" in genes:
        logging.logPrint("You must choose either genes or reference..exiting")
        sys.exit()
    logging.logPrint("Checking the path of dependencies")
    for dependency in dependencies:
        result = is_tool(dependency)
        if result is False:
            print("%s is not in your path, but needs to be!" % dependency)
            sys.exit()
    genome_path = os.path.abspath("%s" % genomes)
    #Need to split reads from reference and also generate combined file
    if "NULL" not in ref:
        logging.logPrint("Prepping sequences")
        reads = split_sequence_by_window(ref, step_size, frag_length)
        #This is a dictionary of ID:read for the reference sequence
        fasta_dict = read_sequences(reads)
        outfile = open("query_sequences.fasta", "w")
        for k,v in fasta_dict.items():
            outfile.write(">%s\n%s\n" % (k,v))
        outfile.close()
    elif "NULL" not in genes:
        """Need to rename these sequences before I process"""
        outfile = open("query_sequences.fasta", "w")
        with open(genes) as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                outfile.write(">%d\n" % record_count_1.next())
                outfile.write(str(record.seq)+"\n")
        outfile.close()
    if os.path.isfile("combined.seqs"):
        pass
    else:
        process_fastas(genome_path, "combined.seqs")
    check_tree_and_reads("combined.seqs", tree)
    os.system("makeblastdb -in combined.seqs -dbtype nucl > /dev/null 2>&1")
    num_refs = get_ref_numbers("combined.seqs")
    """This function will likely need to be also changed if gene sequences are to be used"""
    if "NULL" not in genes:
        fasta_dict = {}
        with open("query_sequences.fasta") as my_fasta:
            for record in SeqIO.parse(my_fasta, "fasta"):
                fasta_dict.update({record.id:record.seq})
    logging.logPrint("Number of sequences to process = %s" % len(fasta_dict))
    logging.logPrint("Starting the loop")
    tree_loop(fasta_dict, "combined.seqs", tree, parallel_workers, num_refs)
    logging.logPrint("Loop finished")
    outfile = open("tmp.txt", "w")
    outfile.write("sequence\tRF\tEUC\t#polymorphisms\tcontig_length\n")
    outfile.close()
    os.system("cat tmp.txt all_distances.txt > results.txt")
    process_tmp_trees()
    logging.logPrint("Cleaning up")
    try:
        subprocess.check_call("rm length.txt distance.txt name.txt polys.txt tmp.txt *euc_dist.txt all_distances.txt combined.seqs*", shell=True)
    except:
        pass
    if "NULL" not in genes:
        os.system("rm query_sequences.fasta")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage="usage: %prog [options]",version="%prog 1.8")
    parser.add_option("-r", "--ref_file", dest="ref",
                      help="/path/to/reference_genome [optional], use if genes not used",
                      action="callback", callback=test_file, type="string", default="NULL")
    parser.add_option("-g", "--genes", dest="genes",
                      help="/path/to/genes file, if desired [optional], use if reference not used",
                      action="callback", callback=test_file, type="string", default="NULL")
    parser.add_option("-d", "--genome_directory", dest="genomes",
                      help="/path/to/genome_directory, files ending in .fasta [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-t", "--tree", dest="tree",
                      help="/path/to/reference_tree [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-s", "--step", dest="step_size",
                      help="step size for shredding sequences, defaults to 5",
                      default="5", type="int")
    parser.add_option("-l", "--frag_length", dest="frag_length",
                      help="shred sequences into given length, defaults to 500",
                      default="500", type="int")
    parser.add_option("-p", "--parallel_workers", dest="parallel_workers",
                      help="How much work to do in parallel, defaults to 2, should number of CPUs your machine has",
                      default="2", type="int")
    parser.add_option("", "--debug", dest="debug",
                      help="Turn debug statements on",
                      action="store_true", default=False)

    options, args = parser.parse_args()

    mandatories = ["genomes","tree"]
    for m in mandatories:
        if not getattr(options, m, None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    logging.DEBUG = options.debug

    main(options.ref, options.genes, options.genomes, options.tree, options.step_size,
         options.frag_length, options.parallel_workers)
