#!/usr/bin/env python

"""pull sequences out of a set of genomes,
concatenate them for phylogeny"""

import optparse
import glob
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import os
import subprocess
import random

def get_seq_name(fasta_in):
    name = os.path.basename(fasta_in)
    return name

def split_multifasta(fasta_in):

def combine_seqs(dir_path):
    handle = open("combined.seqs", "w")
    for infile in glob.glob(os.path.join(dir_path, '*.fas')):
        names = get_seq_name(infile)
        reduced = names.replace('.fas','')
        print >> handle, ">"+str(reduced)
        for record in SeqIO.parse(open(infile), "fasta"):
            print >> handle, record.seq
    handle.close()
     
def run_blast(start_dir):
    for infile in glob.glob(os.path.join(start_dir, 'subsample_*.txt')):
        names = get_seq_name(infile)
        reduced = names.replace('.fasta','')
        cmd = ["blastall",
               "-p", "blastn",
               "-i", infile,
               "-d", "combined.seqs",
               "-o", "%s.blast.out" % reduced,
               "-m", "7",
               "-q", "-4",
               "-r", "5",
               "-a", "2",
               "-F", "F"]
        subprocess.check_call(cmd)

def parse_blast_xml_report(start_dir):
    """uses biopython to split the output
    from a blast file with xml output"""
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(start_dir, '*.blast.out')):
        names = get_seq_name(infile)
        reduced = names.replace('.blast.out','')
        result_handle = open(infile, "U")
        blast_records = NCBIXML.parse(result_handle)
        blast_record = blast_records.next()    
        handle = open("%s.blast.parsed" % reduced, "w") 
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                test = Seq(hsp.sbjct)
                if int(hsp.query_start)<int(hsp.query_end):
                    print >> handle, ">", alignment.title, test
                if int(hsp.query_start)>int(hsp.query_end):
                    print >> handle, ">", alignment.title, test.reverse_complement()
        handle.close()
        result_handle.close()
        os.system("sort -u -k 3,3 %s.blast.parsed > %s.blast.unique" % (reduced, reduced))
        
def parsed_blast_to_seqs(start_dir):
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(start_dir, '*.blast.unique')):
        names = get_seq_name(infile)
        reduced = names.replace('.blast.unique','')
        outfile = open("%s.extracted.seqs" % reduced, "w")
        for line in open(infile, "U"):
            fields = line.split(" ")
            print >> outfile, fields[0] + fields[2], "\n", fields[3],
        outfile.close()

def get_names(fasta_file):
   names = [ ]
   for record in SeqIO.parse(open(fasta_file), "fasta"):
      names.append(record.id)
   return names

def split_files():
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, '*.extracted.seqs')):
        names = get_seq_name(infile)
        reduced = names.replace('.extracted.seqs','')
        names = get_names(infile)
        for name in names:
            handle = open("%s.seqs.fasta" % name, "a") 
            for record in SeqIO.parse(open(infile), "fasta"):
                if name == record.id:
                    print >> handle, ">"+str(record.id)
                    print >> handle, record.seq
            handle.close()

def process_fastas():
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, '*seqs.fasta')):
        handle = open("%s.concat" % infile, "w") 
        names = get_seq_name(infile)
        print >> handle, ">"+str(names),"\n",
        for record in SeqIO.parse(open(infile), "fasta"):
            seqs = [ ]
            seqs.append(record.seq)
            print >> handle, "".join([str(x) for x in seqs]),"\n",

def select_random_seqs(seq_path, markers):
    all_seqs = [ ]
    for record in SeqIO.parse(open(seq_path, "U"), "fasta"):
        all_seqs.append(record.id)
    for i in range(1, 20):
        seqrecords=[ ]
        outfile=open("subsample_%s.txt" % i, "w")
        outseqs=random.sample(set(all_seqs), int(markers))
        for record in SeqIO.parse(open(seq_path, "U"), "fasta"):
            if record.id in outseqs:
                seqrecords.append(record)
        SeqIO.write(seqrecords, outfile, "fasta")
        outfile.close()
        
def main(directory, seqs, markers, tree):
    dir_path=os.path.abspath("%s" % directory)
    seq_path=os.path.abspath("%s" % seqs)
    start_dir = os.getcwd()
    combine_seqs(dir_path)
    os.system("formatdb -i combined.seqs -p F")
    select_random_seqs(seq_path, markers)
    run_blast(start_dir)
    parse_blast_xml_report(start_dir)
    parsed_blast_to_seqs(start_dir)
    #os.system("rm *.blast.out *.blast.parsed *.blast.unique")
    #split_files()
    #process_fastas()
    #os.system("cat *.seqs.fasta.concat > tmp_concatenated")
    #os.system('sed "s/ //g" tmp_concatenated > all_concatenated')
    #os.system("muscle -in all_concatenated -out all_concatenated_aligned.fasta")
    #os.system("rm *seqs.fasta* tmp_concatenated all_concatenated")
    
def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastas cannot be found"
        sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--directory", dest="directory",
                      help="path to genome files in fasta format [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-s", "--seqs", dest="seqs",
                      help="path to all_reads.fasta [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-m", "--markers", dest="markers",
                      help="number of markers to keep, defaults to 3",
                      type="int", action="store", default="3")
    parser.add_option("-t", "--tree", dest="tree",
                      help="path to WGA tree file [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    options, args = parser.parse_args()
    
    mandatories = ["directory", "seqs", "tree"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory, options.seqs, options.markers, options.tree)
