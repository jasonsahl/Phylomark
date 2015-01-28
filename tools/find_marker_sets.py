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
import errno

def get_seq_name(fasta_in):
    name = os.path.basename(fasta_in)
    return name

def split_multi_fasta(fasta_in):
    curr_dir=os.getcwd()
    for record in SeqIO.parse(open(fasta_in, "U"), "fasta"):
        f_out = os.path.join(curr_dir,record.id+'.fasta')
        SeqIO.write([record],open(f_out,'w'),"fasta")
    
def combine_seqs(dir_path):
    handle = open("combined.seqs", "w")
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        names = get_seq_name(infile)
        reduced = names.replace('.fasta','')
        print >> handle, ">"+str(reduced)
        for record in SeqIO.parse(open(infile), "fasta"):
            print >> handle, record.seq
    handle.close()
     
def run_blast(infile):
    names = get_seq_name(infile)
    reduced = names.replace('.fasta','')
    cmd = ["blastall",
            "-p", "blastn",
            "-i", infile,
            "-d", "combined.seqs",
            "-o", "%s.blast.out" % reduced,
            "-m", "7",
            "-b", "2000",
            "-v", "2000",
            "-q", "-4",
            "-r", "5",
            "-a", "8",
            "-F", "F"]
    subprocess.check_call(cmd)

def parse_blast_xml_report(blast_output):
    """uses biopython to split the output
    from a blast file with xml output"""
    names = get_seq_name(blast_output)
    reduced = names.replace('.blast.out','')
    result_handle = open(blast_output, "U")
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
        
def parsed_blast_to_seqs(parsed):
    names = get_seq_name(parsed)
    reduced = names.replace('.blast.unique','')
    outfile = open("%s.extracted.seqs" % reduced, "w")
    for line in open(parsed, "U"):
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
        
def parse_hashrf_file(infile):
    for line in open(infile, "U"):
        if "<0,1>" in line:
            fields = line.split(" ")
            rf = fields[1]
    return rf
    handle.close()

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

def run_loop(seq_path, markers, start_dir, tree_path, iterations):
    all_seqs = [ ]
    new_dir=os.getcwd()
    out_results = open("marker_results.txt", "a")
    for record in SeqIO.parse(open(seq_path, "U"), "fasta"):
        all_seqs.append(record.id)
    for i in range(1, iterations):
        seqrecords=[ ]
        outfile=open("subsample_%s.txt" % i, "w")
        outseqs=random.sample(set(all_seqs), int(markers))
        for record in SeqIO.parse(open(seq_path, "U"), "fasta"):
            if record.id in outseqs:
                seqrecords.append(record)
        SeqIO.write(seqrecords, outfile, "fasta")
        outfile.close()
    for infile in glob.glob(os.path.join(new_dir, 'subsample_*')):
        name=get_seq_name(infile)
        headers=[ ]
        for record in SeqIO.parse(infile, "fasta"):
            headers.append(record.id)
        split_multi_fasta(infile)
        for my_file in glob.glob(os.path.join(new_dir, '*.fasta')):
            run_blast(my_file)
        for blast_output in glob.glob(os.path.join(new_dir, '*blast.out')):
            parse_blast_xml_report(blast_output)
        for parsed in glob.glob(os.path.join(new_dir, '*blast.unique')):
            parsed_blast_to_seqs(parsed)
        os.system("rm *.blast.out *.blast.parsed *.blast.unique")
        split_files()
        process_fastas()
        os.system("cat *.seqs.fasta.concat > tmp_concatenated")
        os.system('sed "s/ //g" tmp_concatenated > all_concatenated')
        os.system('sed "s/ //g" tmp_concatenated > all_concatenated')
        os.system('sed "s/.seqs.fasta//g" all_concatenated > reduced_concatenated')
        os.system("muscle -in reduced_concatenated -out all_concatenated_aligned.fasta > /dev/null 2>&1")
        os.system("FastTree -nt -noboot all_concatenated_aligned.fasta > tmp.tree 2> /dev/null")
        os.system("cat %s %s > combined.tree" % (tree_path, "tmp.tree"))
        os.system("hashrf %s 2 -p list -o %s > /dev/null 2>&1" % ("combined.tree", "result.rf"))
        rf=parse_hashrf_file("result.rf")
        print >> out_results,"\t".join(headers),"\t",rf,
        print "%s processed" % name
        os.system("rm *fasta* tmp_concatenated all_concatenated")
        
def main(directory, seqs, markers, tree, iterations):
    dir_path=os.path.abspath("%s" % directory)
    seq_path=os.path.abspath("%s" % seqs)
    tree_path=os.path.abspath("%s" % tree)
    start_dir = os.getcwd()
    ap=os.path.abspath("%s" % start_dir)
    combine_seqs(dir_path)
    try:
        os.makedirs('%s/scratch' % ap)
    except OSError, e:
     	if e.errno != errno.EEXIST:
            raise
    os.system("mv combined.seqs %s/scratch" % ap)
    os.chdir("%s/scratch" % ap)
    os.system("formatdb -i combined.seqs -p F")
    run_loop(seq_path, markers, ap, tree_path, iterations)
    os.system("cp marker_results.txt %s" % ap)
    os.chdir(ap)
    os.system('rm -rf %s/scratch' % ap)
    
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
    parser.add_option("-i", "--iterations", dest="iterations",
                      help="number of iterations to process, defaults to 20",
                      type="int", action="store", default="20")
    options, args = parser.parse_args()
    
    mandatories = ["directory", "seqs", "tree"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory, options.seqs, options.markers, options.tree, options.iterations)
