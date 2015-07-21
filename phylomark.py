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
        print 'genes file cannot be opened'
        sys.exit()

def test_options(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported.  Only select from T and F"
        sys.exit()

def main(alignment, mask, ref, combined, tree, step_size, frag_length, keep_length, parallel_workers, run_r):
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
    dependencies = ['mothur','muscle','FastTree','blastn','makeblastdb']
    logging.logPrint("Checking the path of dependencies")
    for dependency in dependencies:
        ra = subprocess.call(['which', '%s' % dependency])
        if ra == 0:
            pass
        else:
            print "%s is not in your path, but needs to be!" % dependency
            sys.exit()
    logging.logPrint("Prepping sequences")
    check_tree_and_reads(combined, tree)
    #reads = split_sequence_by_window(alignment, step_size, frag_length)
    #write_sequences(reads)
    #qual_reads = split_quality_values(mask, step_size, frag_length)
    #write_qualities(qual_reads)
    #split_read("quals_shredded.txt")
    #sum_qual_reads("padded_quals.txt")
    #filter_lines_by_value("summed_qualities.txt", keep_length)
    #get_seqs_by_id("seqs_shredded.txt", "seq_names_over_value.txt", "query_sequences.fas")
    #os.system("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % ref)
    #logging.logPrint("Blasting to find contiguous sequences")
    #blast_against_single("query_sequences.fas", ref, "6")
    #filter_blast_report("blast_one.out", frag_length)
    os.system("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % combined)
    #os.system("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % alignment)
    num_refs = get_ref_numbers(combined)
    #fastadir = get_reduced_seqs_by_id("query_sequences.fas", "continuous_seq_names.txt")
    fastadir = split_seqs(alignment)
    logging.logPrint("Starting the loop")
    tree_loop(fastadir, combined, tree, parallel_workers, run_r, num_refs)
    logging.logPrint("Loop finished")
    subprocess.check_call("awk '{print $1}' all_distances.txt > names.txt", shell=True) #should I sort?
    pull_line("names.txt", "summed_qualities.txt", "reduced_quals.txt")
    merge_files_by_column(0, "all_distances.txt", "summed_qualities.txt", "results.txt")
    logging.logPrint("Cleaning up")
    try:
        subprocess.check_call("rm reduced_quals.txt summed_qualities.txt quals_shredded.txt padded_quals.txt blast* continuous* distance.txt name.txt seq_names_over_value.txt", shell=True)
    except:
        sys.exc_clear()
    cleanup_tmpdirs(fastadir)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-a", "--alignment", dest="alignment",
                      help="/path/to/alignment [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-m", "--mask", dest="mask",
                      help="/path/to/filter_mask [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-r", "--ref_file", dest="ref",
                      help="/path/to/reference_genome [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-c", "--combined_seqs", dest="combined",
                      help="/path/to/multifasta_references [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-t", "--tree", dest="tree",
                      help="/path/to/reference_tree [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-s", "--step", dest="step_size",
                      help="step size for shredding sequences",
                      default="5", type="int")
    parser.add_option("-l", "--frag_length", dest="frag_length",
                      help="shred sequences into given length",
                      default="500", type="int")
    parser.add_option("-k", "--keep_length", dest="keep_length",
                      help="keep if polymorphisms greater than value",
                      default="50", type="int")
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

    mandatories = ["alignment", "combined", "tree"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    logging.DEBUG = options.debug

    main(options.alignment, options.mask, options.ref, options.combined, options.tree, options.step_size,
         options.frag_length, options.keep_length, options.parallel_workers, options.run_r)
