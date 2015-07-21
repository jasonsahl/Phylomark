#!/usr/bin/env python

import sys
import os
import subprocess
import tempfile
import string
import itertools
import threading
import optparse
import dendropy

try:
    from igs.utils import functional as func
    from igs.utils import logging
    from igs.threading import functional as p_func
except:
    print "you need to add Phylomark to your PYTHONPATH!"
    sys.exit()

try:
    from Bio import SeqIO
    from Bio.Blast import NCBIXML
    from Bio.Seq import reverse_complement
    from Bio.Seq import Seq
    from Bio import Phylo
except:
    print "BioPython isn't in your PATH, but needs to be"
    sys.exit()

class Increments:
    def __init__(self, start, increment):
        self.state = start
        self.p_start = start
        self.p_increment = increment

    def next(self):
        self.state += self.p_increment
        return self.state

    def reset(self):
        self.state = self.p_start

record_count_1 = Increments(1, 1)
record_count_2 = Increments(1, 1)

def paste_files(name_file, distance_file, all_distance_file):
    handle = open(all_distance_file, "w")
    output = []
    distance_file_lines = open(distance_file, "rU").readlines()
    for lines in zip(open(name_file, "rU"), open(distance_file, "rU")):
        handle.write("\t".join([s.strip() for s in lines]) + "\n")
    handle.close()

def split_sequence_by_window(input_file, step_size, frag_length):
    """cuts up fasta sequences into given chunks"""
    infile = open(input_file, "rU")
    first_record = list(itertools.islice(SeqIO.parse(infile,"fasta"), 1))[0]
    return sliding_window(first_record.seq, frag_length, step_size)

def write_sequences(reads):
    """write shredded fasta sequences to disk"""
    handle = open("seqs_shredded.txt", "w")
    for read in reads:
        print >> handle, ">%d\n%s" % (record_count_1.next(), read)
    handle.close()

def split_quality_values(qual_file, step_size, frag_length):
    infile = open(qual_file, "rU")
    instring = infile.readlines()
    infile.close()
    qual_reads = sliding_window(','.join(instring), frag_length, step_size)
    return qual_reads

def sliding_window(sequence, frag_length, step_size=5):
    """cuts up sequence into a given length"""
    numOfChunks = (len(sequence) - frag_length) + 1
    for i in range(0, numOfChunks, step_size):
        yield sequence[i:i + frag_length]

def write_qualities(qual_reads):
    """write shredded quality files to disk"""
    handle = open("quals_shredded.txt", "w")
    for read in qual_reads:
        print >> handle, read
    handle.close()

def split_read(input_file, output_file):
    """insert gaps into quality files - needed to put into array"""
    handle = open(output_file, "w")
    qual_lines = open(input_file, "rW")
    for line in qual_lines:
        line = line.strip()
        print >> handle, ' '.join(line)
    handle.close()

def sum_qual_reads(input_file, output_file):
    """adds up quality values for each sequence"""
    handle = open(output_file, "w")
    padded_qual_lines = open(input_file, "rU")
    for line in padded_qual_lines.xreadlines():
        sum_values = sum([int(s) for s in line.split()])
        print >> handle, sum_values
    handle.close()

def filter_lines_by_value(filter_in, keep_length):
    handle = open("seq_names_over_value.txt", "w")
    for line in open(filter_in):
        fields = line.split(" ")
        if int(fields[1]) >= keep_length:
            print >> handle, line,
    handle.close()

def get_seqs_by_id(fasta_file, names_file, out_file):
    """retrieves sequences from a large fasta file
    with matching fasta header"""
    names = [">" + l.strip().split()[0] for l in open(names_file)]
    fout = open(out_file, "w")
    print_line = False
    for line in open(fasta_file):
        line = line.strip()
        if line[0] == ">":
            if line in names:
                print_line = True
            else:
                print_line = False

        if print_line:
            fout.write(line)
            fout.write("\n")
    fout.close()

def format_blast_database(ref_file):
    cmd = ["formatdb",
           "-i", ref_file,
           "-p", "F"]
    subprocess.check_call(cmd)

def blast_against_reference(blast_in, combined, outfile):
    try:
        #print 'blastn -query %s -db %s -out %s -dust no -num_alignments 2000 -outfmt "7 std sseq"' % (blast_in,combined,outfile)
        os.system('blastn -query %s -db %s -out %s -dust no -num_alignments 2000 -outfmt "7 std sseq"' % (blast_in,combined,outfile))
    except:
        print "blast problem!"

def blast_against_single(blast_in, ref, blast_type):
    cmd = ["blastn",
           "-query", blast_in,
           "-db", ref,
           "-dust", "no",
           "-evalue", "0.1",
           "-out", "blast_one.out",
           "-outfmt", str(blast_type),
           "-num_alignments", "2000",
           "-num_threads", "2"]
    subprocess.check_call(cmd)

def check_tree_and_reads(combined, tree):
    combined_ids = []
    tree_ids = []
    for record in SeqIO.parse(open(combined), "fasta"):
        combined_ids.append(record.id)
    mytree = Phylo.read(tree, 'newick')
    for clade in mytree.find_clades():
        if clade.name:
            tree_ids.append(clade.name)
    combined_length = len(combined_ids)
    if len(set(combined_ids).intersection(tree_ids)) == int(combined_length):
        pass
    else:
        print "tree names and combined file do not match! exiting"
        sys.exit()

def filter_blast_report(blast_file, frag_length):
    """only return sequences that show a complete
    blast alignment with reference sequence
    will only accept 100% of the frag_length"""
    min_frag_length = int(1 * frag_length)
    handle = open("continuous_seq_names.txt", "w")
    for line in open(blast_file):
        fields = line.split("\t")
        if int(fields[3]) >= min_frag_length and float(fields[2]) >= 0.99:
            print >> handle, fields[0]
    handle.close()

def split_seqs(fasta_file):
    fastadir = tempfile.mkdtemp()
    for record in SeqIO.parse(open(fasta_file), "fasta"):
        f_out = os.path.join(fastadir, record.id + '.fasta')
        SeqIO.write([record], open(f_out, "w"), "fasta")
    return fastadir

def get_reduced_seqs_by_id(fasta_file, names_file):
    """retrieves sequences based on fasta header
    then splits the sequences into a temporary folder"""
    fastadir = tempfile.mkdtemp()
    get_seqs_by_id(fasta_file, names_file, "all_reads.fasta")
    for record in SeqIO.parse(open("all_reads.fasta"), "fasta"):
            f_out = os.path.join(fastadir, record.id + '.fasta')
            SeqIO.write([record], open(f_out, "w"), "fasta")
    return fastadir

def get_ref_numbers(combined):
    records = []
    for record in SeqIO.parse(open(combined), "fasta"):
        records.append(record.id)
    return len(records)


def run_dendropy(tmp_tree, wga_tree, outfile):
    out = open(outfile, "w")
    tree_one = dendropy.Tree.get_from_path(wga_tree,schema="newick",preserve_underscores=True)
    tree_two = dendropy.Tree.get_from_path(tmp_tree,schema="newick",preserve_underscores=True, taxon_set=tree_one.taxon_set)
    #RFs = dendropy.treecalc.robinson_foulds_distance(tree_one, tree_two)
    RFs = tree_one.symmetric_difference(tree_two)
    print >> out, RFs

def tree_loop(fastadir, combined, tree, parallel_workers, run_r, num_refs):
    def _temp_name(t, f):
        return t + '_' + f

    def _perform_workflow(data):
        tn, f = data
        logging.debugPrint(lambda : "Processing file: %s" % f)
        #print tn, f
        blast_against_reference(f, combined, _temp_name(tn, "blast_parsed.txt"))
        subprocess.check_call("sort -u -k 2,2 %s > %s" % (_temp_name(tn, "blast_parsed.txt"),
                                                          _temp_name(tn, "blast_unique.parsed.txt")),
                              shell=True)
        parsed_blast_to_seqs(_temp_name(tn, "blast_unique.parsed.txt"), _temp_name(tn, "seqs_in.fas"))
        subprocess.check_call("muscle -in %s -out %s > /dev/null 2>&1" % (_temp_name(tn, "seqs_in.fas"),
                                                                          _temp_name(tn, "seqs_aligned.fas")),
                              shell=True)
        subprocess.check_call(['mothur',
                               '#filter.seqs(fasta=%s, soft=100, vertical=F)' % _temp_name(tn, "seqs_aligned.fas"), '>', '/dev/null 2>&1'])
        subprocess.check_call('sed "s/[^1]/0/g" %s | sed "s/0/2/g" | sed "s/1/0/g" | sed "s/2/1/g" > %s' % (_temp_name(tn, "seqs_aligned.filter"),
                                                                                                            _temp_name(tn, "mask.txt")), shell=True)
        split_read(_temp_name(tn, "mask.txt"),_temp_name(tn, "padded.txt"))
        sum_qual_reads(_temp_name(tn, "padded.txt"), _temp_name(tn,"polys.txt"))
        if "T" == run_r:
            name = get_seq_name(f)
            subprocess.check_call("cat snps.r | R --slave --args %s %s.table %s.pdf 2> /dev/null" % (_temp_name(tn, "seqs_aligned.fas"), name, name),
        					  shell=True)
            os.system("mv %s.table ./R_output/%s.table.txt" % (name, name))
            os.system("mv %s.pdf ./R_output/%s.plots.pdf" % (name, name))
        else:
            pass
        subprocess.check_call("FastTree -nt -noboot %s > %s 2> /dev/null" % (_temp_name(tn, "seqs_aligned.fas"),
                                                                             _temp_name(tn, "tmp.tree")),
                              shell=True)
        run_dendropy("%s" % (_temp_name(tn, "tmp.tree")), tree, "%s" % (_temp_name(tn, "tmp.RF")))
        thread_id = id(threading.current_thread())
        thread_distance_file = str(thread_id) + '_distance.txt'
        parse_rf_file(_temp_name(tn, "tmp.RF"), thread_distance_file)
        thread_name_file = str(thread_id) + '_name.txt'
        write_strip_name(f, thread_name_file)
        polys_name_file = str(thread_id) + '_polys.txt'
        parse_poly_file(_temp_name(tn, "polys.txt"), polys_name_file)
        subprocess.check_call(["rm",
                               _temp_name(tn, "blast_parsed.txt"),
                               _temp_name(tn, "blast_unique.parsed.txt"),
                               _temp_name(tn, "seqs_in.fas"),
                               _temp_name(tn, "seqs_aligned.fas"),
                               _temp_name(tn, "tmp.tree"),
                               _temp_name(tn, "tmp.RF"),
                               _temp_name(tn, "mask.txt"),
                               _temp_name(tn, "padded.txt"),
                               _temp_name(tn, "polys.txt")])
        return (thread_distance_file, thread_name_file, polys_name_file)

    files = os.listdir(fastadir)
    files_and_temp_names = [(str(idx), os.path.join(fastadir, f))
                            for idx, f in enumerate(files)]

    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=parallel_workers))

    #I do this to make sure and remove any old files that are setting around
    subprocess.call("rm distance.txt name.txt", shell=True, stderr=open(os.devnull, 'w'))

    for files in func.chunk(5, results):
        print files
        distances = [d for d, _ in files]
        names = [n for _, n in files]
        polys = [p for _, p in files]
        subprocess.check_call("cat %s >> distance.txt" % " ".join(distances), shell=True)
        subprocess.check_call("cat %s >> name.txt" % " ".join(names), shell=True)
        subprocess.check_call("cat %s >> polys.txt" % " ".join(polys), shell=True)
        subprocess.check_call("rm %s" % " ".join(distances), shell=True)
        subprocess.check_call("rm %s" % " ".join(names), shell=True)
        subprocess.check_call("rm %s" % " ".join(polys), shell=True)
    paste_files("name.txt", "distance.txt", "all_distances.txt")

def pull_line(names_in, quality_in, out_file):
    handle = open(out_file, "w")
    counts = open(quality_in)
    names = [l.rstrip("\n") for l in open(names_in)]
    for line in counts:
        fields = line.strip().split("\t")
        if fields[0] in names:
            print >> handle, line,
    handle.close()

def merge_files_by_column(column, file_1, file_2, out_file):
    """Takes 2 file and merge their columns based on the column. It is assumed
    that line ordering in the files do not match, so we read both files into memory
    and join them"""
    join_map = {}
    for line in open(file_1):
        row = line.split()
        column_value = row.pop(column)
        join_map[column_value] = row

    for line in open(file_2):
        row = line.split()
        column_value = row.pop(column)
        if column_value in join_map:
            join_map[column_value].extend(row)

    fout = open(out_file, 'w')
    for k, v in join_map.iteritems():
        fout.write('\t'.join([k] + v) + '\n')

    fout.close()

def cleanup_tmpdirs(fastadir):
    subprocess.check_call(["rm", "-rf", fastadir])

def parse_blast_xml_report(blast_file, outfile):
    """uses biopython to split the output
    from a blast file with xml output"""
    result_handle = open(blast_file)
    blast_records = NCBIXML.parse(result_handle)
    blast_record = blast_records.next()
    handle = open(outfile, "w")
    for alignment in blast_record.alignments:
         for hsp in alignment.hsps:
             test = Seq(hsp.sbjct)
             if int(hsp.query_start)<int(hsp.query_end):
                 print >> handle, ">", alignment.title, test
             if int(hsp.query_start)>int(hsp.query_end):
                 print >> handle, ">", alignment.title, test.reverse_complement()
    handle.close()
    result_handle.close()

def parsed_blast_to_seqs(parsed_file, outfile):
    infile = open(parsed_file, "rU")
    handle = open(outfile, "w")
    for line in infile:
        if line.startswith("#"):
            pass
        else:
            fields = line.split()
            print >> handle, ">"+str(fields[1])
            print >> handle, fields[12]
    handle.close()

def parse_rf_file(infile, outfile):
    handle = open(outfile, "a")
    for line in open(infile):
        print >> handle, line,
    handle.close()

def parse_poly_file(infile, outfile):
    handle = open(outfile, "a")
    for line in open(infile):
        print >> handle, line,
    handle.close()

def write_strip_name(filename, outfile):
    handle = open(outfile, "a")
    filename = os.path.splitext(os.path.basename(filename))[0]
    print >> handle, filename
    handle.close()
