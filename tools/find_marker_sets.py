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

def get_seq_name(fasta_in):
    name = os.path.basename(fasta_in)
    return name

def combine_seqs(dir_path):
    handle = open("combined.seqs", "w")
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        names = get_seq_name(infile)
        reduced = names.replace('.fasta','')
        print >> handle, ">"+str(reduced)
        for record in SeqIO.parse(open(infile), "fasta"):
            print >> handle, record.seq
    handle.close()
     
def run_blast(seq_path):
    for infile in glob.glob(os.path.join(seq_path, '*.fasta')):
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

def parse_blast_xml_report():
    """uses biopython to split the output
    from a blast file with xml output"""
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, '*.blast.out')):
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
        
def parsed_blast_to_seqs():
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, '*.blast.unique')):
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

def main(directory, seqs):
    dir_path=os.path.abspath("%s" % directory)
    seq_path=os.path.abspath("%s" % seqs)
    combine_seqs(dir_path)
    os.system("formatdb -i combined.seqs -p F")
    run_blast(seq_path)
    parse_blast_xml_report()
    parsed_blast_to_seqs()
    os.system("rm *.blast.out *.blast.parsed *.blast.unique")
    split_files()
    process_fastas()
    os.system("cat *.seqs.fasta.concat > tmp_concatenated")
    os.system('sed "s/ //g" tmp_concatenated > all_concatenated')
    os.system("muscle -in all_concatenated -out all_concatenated_aligned.fasta")
    os.system("rm *seqs.fasta* tmp_concatenated all_concatenated")
    
def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastas cannot be found"
        sys.exit()

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--directory", dest="directory",
                      help="path to genome files in fasta format [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-s", "--seqs", dest="seqs",
                      help="path to sequence directory [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-m", "--markers", dest="markers",
                      help="path to genome files in fasta format [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    options, args = parser.parse_args()
    
    mandatories = ["directory", "seqs"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory, options.seqs)
