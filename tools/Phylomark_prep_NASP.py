#!/usr/bin/env python

""""Creates Phylomark input
from a NASP matrix and a
directory of genomes"""


import os
import sys
import optparse
import subprocess

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print 'genes file cannot be opened'
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastas cannot be found"
        sys.exit()

def filter_matrix(in_matrix):
    """filter NASP matrix for a user provided list of genome
    coordinates"""
    infile=iter(fileinput.input([in_matrix]))
    out_file=open("nasp_matrix_dashes.txt", "w")
    firstLine = infile.readline()
    print >> out_file, firstLine,
    lines = [ ]
    for line in infile:
        new_fields = [ ]
        fields=line.split("\t")
        for field in fields:
            if field == "X":
                new_fields.append("-")
            elif field == "N":
                new_fields.append("-")
            elif field == "n":
                new_fields.append("-")
            elif field == ".":
                new_fields.append("-")
            else:
                new_fields.append(field)
        print >> out_file,"\t".join(new_fields),

def get_field_index(matrix_in):
    """index to find where the SNP calls end"""
    matrix=open(matrix_in, "rU")
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    matrix.close()
    return last
        
def matrix_to_fasta(matrix_in, last):
    """converts an ISG matrix to fasta format"""
    reduced = [ ]
    out_fasta = open("nasp_in.fasta", "w")
    for line in open(matrix_in):
        fields = line.split("\t")
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        out_fasta.write(">"+str(x[0])+"\n")
        out_fasta.write("".join(x[1:])+"\n")
    out_fasta.close()
                
def main(nasp_matrix, directory, reference):
    last = get_field_index(nasp_matrix)
    filter_matrix(nasp_matrix)
    matrix_to_fasta("nasp_matrix_dashes.txt", last)
    subprocess.check_call(['mothur',
                           '#filter.seqs(fasta="nasp_in.fasta", vertical=F, trump=-)' % final_fasta])
    subprocess.check_call(['mothur',
                           '#filter.seqs(fasta=nasp_in.filtered.fasta, soft=100, vertical=F)' % final_filter_fasta])
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-m", "--nasp_matrix", dest="nasp_matrix",
                      help="/path/to/NASP matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-d", "--directory", dest="directory",
                      help="/path/to/genomes in FASTA format [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-r", "--reference", dest="reference",
                      help="name of reference in FASTA directory (exclude the .fasta) [REQUIRED]",
                      action="store", type="string")                  
    options, args = parser.parse_args()
    
    mandatories = ["nasp_matrix", "directory", "reference"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    logging.DEBUG = options.debug
            
    main(options.nasp_matrix, directory, reference)
