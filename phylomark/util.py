#!/usr/bin/env python

import sys
import os
import subprocess
import tempfile
import string
import itertools
import threading
import optparse
try:
    import dendropy
    from dendropy.calculate import treecompare
    #import dendropy.treecalc
except:
    print("dendropy is not installed or needs to be updated!")
    sys.exit()
import glob

"""This dependency will need to be removed"""
try:
    from igs.utils import functional as func
    from igs.utils import logging
    from igs.threading import functional as p_func
except:
    print("IGS not found in your PYTHONPATH!")
    sys.exit()
try:
    from Bio import SeqIO
    from Bio.Seq import reverse_complement
    from Bio.Seq import Seq
    from Bio import Phylo
except:
    print("BioPython isn't in your PATH, but needs to be")
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

def paste_files(name_file, distance_file, euc_file, poly_file, length_file, all_distance_file):
    handle = open(all_distance_file, "w")
    output = []
    distance_file_lines = open(distance_file).readlines()
    for lines in zip(open(name_file), open(distance_file), open(euc_file), open(poly_file), open(length_file)):
        handle.write("\t".join([s.strip() for s in lines]) + "\n")
    handle.close()

def split_read(input_file, output_file):
    """insert gaps into quality files - needed to put into array"""
    handle = open(output_file, "w")
    with open(input_file) as qual_lines:
        for line in qual_lines:
            line = line.strip()
            handle.write(' '.join(line))
            handle.write("\n")
    handle.close()

def sum_qual_reads(input_file, output_file):
    """adds up quality values for each sequence"""
    handle = open(output_file, "w")
    with open(input_file) as padded_qual_lines:
        for line in padded_qual_lines.readlines():
            sum_values = sum([int(s) for s in line.split()])
            handle.write(str(sum_values))
            handle.write("\n")
    handle.close()

def get_seqs_by_id(fasta_file, names_file, out_file):
    """retrieves sequences from a large fasta file
    with matching fasta header"""
    names = [">" + l.strip().split()[0] for l in open(names_file)]
    fout = open(out_file, "w")
    print_line = False
    with open(fasta_file) as my_file:
        for line in my_file:
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
        os.system('blastn -task blastn -query %s -db %s -out %s -dust no -num_alignments 2000 -outfmt "7 std sseq"' % (blast_in,combined,outfile))
    except:
        print("blast problem!")
        sys.exit()

def blast_against_single(blast_in, ref, blast_type):
    cmd = ["blastn",
           "-task", "blastn",
           "-query", blast_in,
           "-db", ref,
           "-dust", "no",
           "-evalue", "0.1",
           "-out", "blast_one.out",
           "-outfmt", str(blast_type),
           "-num_alignments", "2000",
           "-num_threads", "2"]
    subprocess.check_call(cmd)

def check_tree_and_reads(combined,tree):
    combined_ids = []
    tree_ids = []
    with open(combined) as my_combined:
        for record in SeqIO.parse(my_combined, "fasta"):
            combined_ids.append(record.id)
    mytree = Phylo.read(tree, 'newick')
    for clade in mytree.find_clades():
        if clade.name:
            tree_ids.append(clade.name)
    combined_length = len(combined_ids)
    if len(set(combined_ids).intersection(tree_ids)) == int(combined_length):
        pass
    else:
        print("tree names and combined file do not match! exiting")
        sys.exit()

def filter_blast_report(blast_file, frag_length):
    """only return sequences that show a complete
    blast alignment with reference sequence
    will only accept 100% of the frag_length"""
    min_frag_length = int(1 * frag_length)
    handle = open("continuous_seq_names.txt", "w")
    with open(blast_file) as my_blast:
        for line in my_blast:
            fields = line.split("\t")
            if int(fields[3]) >= min_frag_length and float(fields[2]) >= 0.99:
                handle.write(fields[0]+"\n")
    handle.close()

def split_seqs(fasta_file):
    fastadir = tempfile.mkdtemp()
    with open(fasta_file) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            f_out = os.path.join(fastadir, record.id + '.fasta')
            SeqIO.write([record], open(f_out, "w"), "fasta")
    return fastadir

def get_reduced_seqs_by_id(fasta_file, names_file):
    """retrieves sequences based on fasta header
    then splits the sequences into a temporary folder"""
    fastadir = tempfile.mkdtemp()
    get_seqs_by_id(fasta_file, names_file, "all_reads.fasta")
    with open("all_reads.fasta") as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            f_out = os.path.join(fastadir, record.id + '.fasta')
            SeqIO.write([record], open(f_out, "w"), "fasta")
    return fastadir

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def get_ref_numbers(combined):
    records = []
    with open(combined) as my_combined:
        for record in SeqIO.parse(my_combined, "fasta"):
            records.append(record.id)
    return len(records)

def process_fastas(directory, out_fasta):
    """make the combined fasta file"""
    fout = open(out_fasta, "w")
    for infile in glob.glob(os.path.join(directory, '*.fasta')):
        names = get_seq_name(infile)
        reduced = names.rstrip('.fasta')
        fout.write('>' + str(reduced) + '\n')
        with open(infile) as my_file:
            for record in SeqIO.parse(my_file, "fasta"):
                fout.write(str(record.seq) + '\n')
        fout.write('\n')
    fout.close()

def get_contig_length(in_fasta, outfile):
    my_out = open(outfile, "w")
    length = []
    with open(in_fasta) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
           length.append(str(len(record.seq)))
    my_out.write("".join(length)+"\n")

def split_sequence_by_window(input_file, step_size, frag_length):
    """cuts up fasta sequences into given chunks"""
    infile = open(input_file)
    first_record = list(itertools.islice(SeqIO.parse(infile,"fasta"), 1))[0]
    return sliding_window(first_record.seq, frag_length, step_size)

def sliding_window(sequence, frag_length, step_size):
    """cuts up sequence into a given length"""
    numOfChunks = (len(sequence) - frag_length) + 1
    for i in range(0, numOfChunks, step_size):
        yield sequence[i:i + frag_length]

def write_sequences(reads):
    """write shredded fasta sequences to disk"""
    handle = open("seqs_shredded.txt", "w")
    for read in reads:
        handle.write(">%d\n%s" % (record_count_1.next(), read))
    handle.close()

def read_sequences(reads):
    read_dict = {}
    for read in reads:
        read_dict.update({"%d" % record_count_1.next():"%s" % read})
    return read_dict

def run_dendropy(tmp_tree,wga_tree,outfile):
    out = open(outfile, "w")
    tree_one = dendropy.Tree.get_from_path(wga_tree,schema="newick",preserve_underscores=True)
    tree_two = dendropy.Tree.get_from_path(tmp_tree,schema="newick",preserve_underscores=True, taxon_namespace=tree_one.taxon_namespace)
    #RFs = tree_one.symmetric_difference(tree_two)
    RFs = treecompare.symmetric_difference(tree_one,tree_two)
    out.write(str(RFs))
    out.write("\n")
    out.close()

def run_dendropy_euclidian(tmp_tree,wga_tree,outfile):
    out = open(outfile, "w")
    tns = dendropy.TaxonNamespace()
    #tree_one = dendropy.Tree.get_from_path(wga_tree,schema="newick",preserve_underscores=True)
    tree_one = dendropy.Tree.get_from_path(wga_tree,"newick",taxon_namespace=tns)
    #tree_two = dendropy.Tree.get_from_path(tmp_tree,schema="newick",preserve_underscores=True, taxon_namespace=tree_one.taxon_namespace)
    tree_two = dendropy.Tree.get_from_path(tmp_tree,"newick",taxon_namespace=tree_one.taxon_namespace)
    eu_dist = dendropy.calculate.treecompare.euclidean_distance(tree_one,tree_two)
    out.write(str(eu_dist))
    out.write("\n")
    out.close()

def check_and_align_seqs(infile, num_genomes, outfile):
    lengths = []
    with open(infile) as my_file:
        for record in SeqIO.parse(my_file, "fasta"):
            lengths.append(len(record.seq))
    lengths.sort(key=int)
    if (lengths[0]/lengths[-1]) > 0.75 and len(lengths) == num_genomes:
        os.system("muscle -in %s -out %s > /dev/null 2>&1" % (infile,outfile))

def tree_loop(fasta_dict, combined, tree, parallel_workers, num_refs):
    def _temp_name(t, f):
        return t + '_' + f

    def _perform_workflow(data):
        tn, f = data
        #tn is the sequence ID, f is the actual sequence
        outfile = open("%s.fasta" % tn, "w")
        outfile.write(">%s\n%s\n" % (tn,f))
        outfile.close()
        logging.logPrint("Processing sequence: %s" % tn)
        blast_against_reference("%s.fasta" % tn, combined, _temp_name(tn, "blast_parsed.txt"))
        subprocess.check_call("sort -u -k 2,2 %s > %s" % (_temp_name(tn, "blast_parsed.txt"),
                                                          _temp_name(tn, "blast_unique.parsed.txt")),shell=True)
        parsed_blast_to_seqs(_temp_name(tn, "blast_unique.parsed.txt"), _temp_name(tn, "seqs_in.fas"))
        check_and_align_seqs(_temp_name(tn, "seqs_in.fas"), num_refs, _temp_name(tn, "seqs_aligned.fas"))
        if os.path.isfile(_temp_name(tn, "seqs_aligned.fas")):
            """What if there are NO SNPs in a given region"""
            subprocess.call(['mothur',
                                   '#filter.seqs(fasta=%s, soft=100, vertical=F)' % _temp_name(tn, "seqs_aligned.fas")], stdout=subprocess.PIPE)
            subprocess.check_call('sed "s/[^1]/0/g" %s | sed "s/0/2/g" | sed "s/1/0/g" | sed "s/2/1/g" > %s' % (_temp_name(tn, "seqs_aligned.filter"),
                                                                                                                _temp_name(tn, "mask.txt")), shell=True)
            split_read(_temp_name(tn, "mask.txt"),_temp_name(tn, "padded.txt"))
            sum_qual_reads(_temp_name(tn, "padded.txt"), _temp_name(tn,"polys.txt"))
            subprocess.check_call("FastTree -nt -noboot %s > %s 2> /dev/null" % (_temp_name(tn, "seqs_aligned.fas"),
                                                                                 _temp_name(tn, "tmp.tree")),shell=True)
            run_dendropy("%s" % (_temp_name(tn, "tmp.tree")), tree, "%s" % (_temp_name(tn, "tmp.RF")))
            run_dendropy_euclidian("%s" % (_temp_name(tn, "tmp.tree")), tree, "%s" % (_temp_name(tn, "tmp.EU")))
            get_contig_length("%s.fasta" % tn, _temp_name(tn, "length.txt"))
            thread_id = id(threading.current_thread())
            thread_distance_file = str(thread_id) + '_distance.txt'
            parse_rf_file(_temp_name(tn, "tmp.RF"), thread_distance_file)
            thread_euclidian_file = str(thread_id) + "_euc_dist.txt"
            parse_rf_file(_temp_name(tn, "tmp.EU"), thread_euclidian_file)
            thread_name_file = str(thread_id) + '_name.txt'
            write_strip_name("%s.fasta" % tn, thread_name_file)
            polys_name_file = str(thread_id) + '_polys.txt'
            parse_poly_file(_temp_name(tn, "polys.txt"), polys_name_file)
            length_name_file = str(thread_id) + '_length.txt'
            parse_poly_file(_temp_name(tn, "length.txt"), length_name_file)
            try:
                subprocess.check_call("rm mothur*", shell=True, stderr=open(os.devnull, 'w'))
            except:
                pass
            subprocess.check_call(["rm",
                                   _temp_name(tn, "blast_parsed.txt"),
                                   "%s.fasta" % tn,
                                   _temp_name(tn, "blast_unique.parsed.txt"),
                                   _temp_name(tn, "seqs_in.fas"),
                                   _temp_name(tn, "seqs_aligned.fas"),
                                   _temp_name(tn, "tmp.tree"),
                                   _temp_name(tn, "tmp.RF"),
                                   _temp_name(tn, "tmp.EU"),
                                   _temp_name(tn, "mask.txt"),
                                   _temp_name(tn, "padded.txt"),
                                   _temp_name(tn, "polys.txt"),
                                   _temp_name(tn, "seqs_aligned.filter"),
                                   _temp_name(tn, "length.txt"),
                                   _temp_name(tn, "seqs_aligned.filter.fasta")])
            return (thread_distance_file, thread_name_file, polys_name_file, length_name_file,
                    thread_euclidian_file)
        else:
            subprocess.check_call(["rm",
                                   _temp_name(tn, "blast_parsed.txt"),
                                   "%s.fasta" % tn,
        #                           _temp_name(tn, "euc_dist.txt"),
                                   _temp_name(tn, "blast_unique.parsed.txt"),
                                   _temp_name(tn, "seqs_in.fas")])

    files_and_temp_names = [(str(idx), f)
                             for idx, f in fasta_dict.items()]
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=parallel_workers))

    #I do this to make sure and remove any old files that are setting around
    subprocess.call("rm distance.txt name.txt polys.txt length.txt euc_dist.txt", shell=True, stderr=open(os.devnull, 'w'))

    for files in func.chunk(5, results):
        distances = []
        names = []
        polys = []
        lengths = []
        euc_dist = []
        for value in files:
            if value:
                distances.append(value[0])
                names.append(value[1])
                polys.append(value[2])
                lengths.append(value[3])
                euc_dist.append(value[4])
        if distances:
            subprocess.check_call("cat %s >> distance.txt" % " ".join(distances), shell=True)
            subprocess.check_call("cat %s >> name.txt" % " ".join(names), shell=True)
            subprocess.check_call("cat %s >> polys.txt" % " ".join(polys), shell=True)
            subprocess.check_call("cat %s >> length.txt" % " ".join(lengths), shell=True)
            subprocess.check_call("cat %s >> euc_dist.txt" % " ".join(euc_dist), shell=True)
            subprocess.check_call("rm %s" % " ".join(distances), shell=True)
            subprocess.check_call("rm %s" % " ".join(names), shell=True)
            subprocess.check_call("rm %s" % " ".join(polys), shell=True)
            subprocess.check_call("rm %s" % " ".join(lengths), shell=True)
    paste_files("name.txt","distance.txt","euc_dist.txt","polys.txt","length.txt","all_distances.txt")

def pull_line(names_in, quality_in, out_file):
    handle = open(out_file, "w")
    names = [l.rstrip("\n") for l in open(names_in)]
    with open(quality_in) as counts:
        for line in counts:
            fields = line.strip().split("\t")
            if fields[0] in names:
                handle.write(line)
                handle.write("\n")
    handle.close()

def merge_files_by_column(column, file_1, file_2, out_file):
    """Takes 2 file and merge their columns based on the column. It is assumed
    that line ordering in the files do not match, so we read both files into memory
    and join them"""
    join_map = {}
    with open(file_1) as f1:
        for line in f1:
            row = line.split()
            column_value = row.pop(column)
            join_map[column_value] = row

    with open(file_2) as f2:
        for line in f2:
            row = line.split()
            column_value = row.pop(column)
            if column_value in join_map:
                join_map[column_value].extend(row)

    fout = open(out_file, 'w')
    for k, v in join_map.items():
        fout.write('\t'.join([k] + v) + '\n')
    fout.close()

def parsed_blast_to_seqs(parsed_file, outfile):
    infile = open(parsed_file, "rU")
    handle = open(outfile, "w")
    with open(parsed_file) as infile:
        for line in infile:
            if line.startswith("#"):
                pass
            else:
                fields = line.split()
                handle.write(">"+str(fields[1])+"\n")
                handle.write(fields[12])
                handle.write("\n")
    handle.close()

def parse_rf_file(infile, outfile):
    handle = open(outfile, "a")
    for line in open(infile):
        handle.write(line)
    handle.close()

def parse_poly_file(infile, outfile):
    handle = open(outfile, "a")
    with open(infile) as my_file:
        for line in my_file:
            handle.write(line)
    handle.close()

def process_tmp_trees():
    import os
    import glob
    cwd = os.getcwd()
    outfile = open("all_trees.txt", "w")
    for infile in glob.glob(os.path.join(cwd, '*_tmp.tree')):
        path_split = infile.split("/")
        name_fields = path_split[-1].split("_")
        name = name_fields[0]
        with open(infile) as my_file:
            first_line = my_file.readline()
            outfile.write("%s" % name + "\t" + first_line)
    outfile.close()
    try:
        os.system("rm %s/*_tmp.tree" % cwd)
    except:
        pass

def write_strip_name(filename, outfile):
    handle = open(outfile, "a")
    filename = os.path.splitext(os.path.basename(filename))[0]
    handle.write(filename)
    handle.write("\n")
    handle.close()
