#!/usr/bin/env python

"""generate all files needed to input
into Phylomark"""
import sys
import os
import re
import glob
import subprocess

try:
    from bx.align.tools.thread import get_components_for_species
    from bx.align import maf
except:
    print "bx-python needs to be installed!"
    sys.exit()
try:  
    from Bio import SeqIO
except:
    print "biopython not installed, but needs to be!"
    sys.exit()
try:
    from igs.utils import cli
    from igs.utils import functional as func
except:
    print "You need to add Phylomark to your PYTHONPATH!"
    sys.exit()

OPTIONS = [
    ('input_maf', '', '--input-maf', 'Input MAF to process', cli.notNone),
    ('fasta_dir', '', '--fasta-dir', 'Directory FASTA files are located, MUST HAVE .fas ENDING', cli.defaultIfNone(os.getcwd())),
    ('tmp_dir', '', '--tmp-dir', 'Temporary directory to put all out files', cli.defaultIfNone('/tmp')),
    ('keep_files', '', '--keep-files', 'Keep any temporary files', cli.defaultIfNone(False), cli.BINARY)
    ]

def maf_to_fasta(in_maf, out_fasta):
    """converts a maf format file to fasta.  Script is taken from
    bx-python-tools found at https://bitbucket.org/james_taylor/bx-python/wiki/Home"""
    fout = open(out_fasta, "w")
    comps = None
    maf_reader = maf.Reader(open(in_maf))
    for i, m in enumerate(maf_reader):
        if comps: l = [m.components[i] for i in comps]
        else: l = m.components
        for c in l:
            fout.write('>%s:%d-%d\n' % (c.src, c.start, c.end))
            fout.write(c.text + '\n')

    fout.close()
            
def get_record_ids(in_fasta, out_names):
    """parse the record names from the MAF-converted FASTA file"""
    fout = open(out_names, "w")
    for record in SeqIO.parse(open(in_fasta), "fasta"):
        R=record.id
        M=re.search('(?<!.)\w+', R)
        fout.write(M.group(0) + '\n')
    fout.close()

def thread_for_species(in_maf, in_names, output):
    """pulls out MAF blocks that contain sequence from all desired genomes.
     Script is taken from bx-python-tools found at 
     https://bitbucket.org/james_taylor/bx-python/wiki/Home"""
    species = open(in_names, "rU").read().splitlines()
    fout = open(output, "w")
    maf_reader = maf.Reader(open(in_maf))
    maf_writer = maf.Writer(fout)
    for m in maf_reader:
        new_components = get_components_for_species(m, species)
        if new_components:
            m.components = new_components
            m.score = 0.0
            maf_writer.write(m)
    maf_reader.close()
    maf_writer.close()
    fout.close()

def maf_to_concat_fasta(in_maf, in_name, output):
    """converts a MAF file to a concatenated FASTA.
    Script is taken from bx-python-tools found at 
     https://bitbucket.org/james_taylor/bx-python/wiki/Home"""
    species = open(in_name, "rU").read().splitlines()
    fout = open(output, "w")
    wrap = 50
    fill = ""
    texts = {}
    for s in species: texts[s] = []
    maf_reader = maf.Reader(open(in_maf))
    for m in maf_reader:
        for s in species:
            c = m.get_component_by_src_start(s)
            if c: texts[s].append(c.text)
            else: texts[s].append('-' * m.text_size)
    for s in species:
        fout.write('>' + s + '\n')
        print_n(fill.join(texts[s]), wrap, fout)
        
def print_n(s, n, fout):
    if n <= 0:
        fout.write(s + '\n')
    else:
        for chunk in func.chunk(n, s):
            out_s = ''.join(chunk)
            fout.write(out_s + '\n')

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)
                
def process_fastas(directory, out_fasta):
    """make the combined fasta file"""
    fout = open(out_fasta, "w")
    for infile in glob.glob(os.path.join(directory, '*.fas')):
        names = get_seq_name(infile)
        reduced = names.rstrip('.fas')
        fout.write('>' + str(reduced) + '\n')
        for record in SeqIO.parse(open(infile), "fasta"):
            fout.write(str(record.seq) + '\n')
        fout.write('\n')
    fout.close()
           

def main(options, _args):
    rc=subprocess.call(['which', 'mothur'])
    if rc == 0:
        pass
    else:
        print "mothur needs to be in your path"
        sys.exit()
    all_fasta = os.path.join(options('general.tmp_dir'), 'all.fasta')
    maf_to_fasta(options('general.input_maf'), all_fasta)
    record_names = os.path.join(options('general.tmp_dir'), 'record_names')
    get_record_ids(all_fasta, record_names)
    unique_lines = set(open(record_names).readlines())
    new_records = open(record_names, "w").writelines(unique_lines)
    reduce_maf = os.path.join(options('general.tmp_dir'), 'reduced.maf')
    thread_for_species(options('general.input_maf'), record_names, reduce_maf)
    final_fasta = os.path.join(options('general.tmp_dir'), 'final.fasta')
    maf_to_concat_fasta(reduce_maf, record_names, final_fasta)
    subprocess.check_call(['mothur',
                           '#filter.seqs(fasta=%s, vertical=F, trump=-)' % final_fasta])
    final_filter_fasta = os.path.join(options('general.tmp_dir'), 'final.filter.fasta')
    subprocess.check_call(['mothur',
                           '#filter.seqs(fasta=%s, vertical=F, trump=.)' % final_filter_fasta])
    final_filter_filter_fasta = os.path.join(options('general.tmp_dir'),
                                             'final.filter.filter.fasta')
    concatenated_alignment_fasta = 'concatenated_alignment.fasta'
    subprocess.check_call(['mv',
                           final_filter_filter_fasta,
                           concatenated_alignment_fasta])
    subprocess.check_call(['mothur',
                           '#filter.seqs(fasta=%s, soft=100, vertical=F)' % concatenated_alignment_fasta])
    concatenated_alignment_filter = 'concatenated_alignment.filter'
    subprocess.check_call('sed "s/[^1]/0/g" %s | sed "s/0/2/g" | sed "s/1/0/g" | sed "s/2/1/g" > mask_in.txt' % concatenated_alignment_filter, shell=True)
    subprocess.check_call('FastTree -nt -noboot concatenated_alignment.fasta > wga.tree', shell=True)
    process_fastas(options('general.fasta_dir'), 'combined.fasta')
    if not options('general.keep_files'):
        subprocess.check_call(['rm',
                               record_names,
                               reduce_maf,
                               final_fasta,
                               final_filter_fasta,
                               concatenated_alignment_filter])

if __name__ == "__main__":
    main(*cli.buildConfigN(OPTIONS))

