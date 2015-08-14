#!/usr/bin/env python

import sys
import optparse
import errno
try:
    from phylomark.util import *
except:
    print "your phylomark environment is not set.  Add this directory to your PYTHONPATH"
    sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s cannot be opened' % option
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastas cannot be found"
        sys.exit()

def test_options(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported.  Only select from T and F"
        sys.exit()

def main(ref, genomes, tree, step_size, frag_length, parallel_workers, run_r):
    if "T" in run_r:
        try:
            os.makedirs('./R_output')
        except OSError, e:
            if e.errno != errno.EEXIST:
        	raise
    else:
        pass
    if "T" in run_r:
        rb = suprocess.call(['which', 'R'])
        if rb == 0:
            pass
        else:
            print "R is not in your path, but needs to be!"
            sys.exit()
    dependencies = ['mothur','muscle','FastTree','blastn','makeblastdb','mothur']
    logging.logPrint("Checking the path of dependencies")
    for dependency in dependencies:
        ra = subprocess.call(['which', '%s' % dependency])
        if ra == 0:
            pass
        else:
            print "%s is not in your path, but needs to be!" % dependency
            sys.exit()
    genome_path = os.path.abspath("%s" % genomes)
    logging.logPrint("Prepping sequences")
    #Need to split reads from reference and also generate combined file
    reads = split_sequence_by_window(ref, step_size, frag_length)
    write_sequences(reads)
    process_fastas(genome_path, "combined.seqs")
    check_tree_and_reads("combined.seqs", tree)
    os.system("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % combined)
    num_refs = get_ref_numbers(combined)
    fastadir = split_seqs("seqs_shredded.txt")
    logging.logPrint("Starting the loop")
    tree_loop(fastadir, "combined.seqs", tree, parallel_workers, run_r, num_refs)
    logging.logPrint("Loop finished")
    outfile = open("tmp.txt", "w")
    print >> outfile, "sequence\tRF\t#polymorphisms\tcontig_length"
    outfile.close()
    os.system("cat tmp.txt all_distances.txt > results.txt")
    logging.logPrint("Cleaning up")
    try:
        subprocess.check_call("rm length.txt distance.txt name.txt polys.txt tmp.txt all_distances.txt", shell=True)
    except:
        sys.exc_clear()
    cleanup_tmpdirs(fastadir)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-r", "--ref_file", dest="ref",
                      help="/path/to/reference_genome [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-d", "--genome_directory", dest="genomes",
                      help="/path/to/genome_directory [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-t", "--tree", dest="tree",
                      help="/path/to/reference_tree [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-s", "--step", dest="step_size",
                      help="step size for shredding sequences",
                      default="5", type="int")
    parser.add_option("-l", "--frag_length", dest="frag_length",
                      help="shred sequences into given length",
                      default="500", type="int")
    parser.add_option("-p", "--parallel_workers", dest="parallel_workers",
                      help="How much work to do in parallel, defaults to 2, should number of CPUs your machine has",
                      default="2", type="int")
    parser.add_option("-u", "--run_r", dest="run_r",
                      help="Run R implementation?  Default is F, modify to T to run",
                      default="F", action="callback", callback=test_options, type="string")
    parser.add_option("", "--debug", dest="debug",
                      help="Turn debug statements on",
                      action="store_true", default=False)

    options, args = parser.parse_args()

    mandatories = ["ref", "genomes", "tree"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    logging.DEBUG = options.debug

    main(options.ref, options.genomes, options.tree, options.step_size,
         options.frag_length, options.parallel_workers, options.run_r)
