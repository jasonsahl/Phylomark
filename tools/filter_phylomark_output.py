#!/usr/bin/python

"""takes a list of record.ids and returns to you the sequences
from a fasta list that are part of the list"""

from Bio import SeqIO
from optparse import OptionParser
import sys

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def get_list(ids,rf):
    infile = open(ids)
    to_keep = []
    with open(infile) as my_file:
        for line in my_file:
            fields = line.split()
            if int(fields[1])>=int(rf):
                to_keep.append(fields[0])
    return to_keep

def filter_reads(in_fasta, to_keep, out_fasta):
    output_handle = open(out_fasta, "w")
    seqrecords=[]
    with open(in_fasta) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in to_keep:
                seqrecords.append(record)
    SeqIO.write(seqrecords, output_handle, "fasta")
    output_handle.close()

def main(in_fasta,ids,out_fasta,rf):
    to_keep=get_list(ids,rf)
    filter_reads(in_fasta,to_keep,out_fasta)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input_fasta", dest="in_fasta",
                    help="/path/to/input fasta [REQUIRED]",
                    action="callback", callback=test_file, type="string")
    parser.add_option("-d", "--headers", dest="ids",
                    help="/path/to/results.txt file [REQUIRED]",
                    action="callback", callback=test_file, type="string")
    parser.add_option("-o", "--output_fasta", dest="out_fasta",
                    help="/path/to/output fasta [REQUIRED]",
                    action="store", type="string")
    parser.add_option("-r", "--rf", dest="rf",
                    help="upper bound for RF value to keep, defaults to 30",
                    default="30", action="store", type="int")

    options, args = parser.parse_args()

    mandatories = ["in_fasta","ids","out_fasta"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.in_fasta,options.ids,options.out_fasta,options.rf)
