#!/usr/bin/env python3

import sys

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord

from intervaltree import IntervalTree

import pysam
import mappy as mp

TE_CUT_SITES = {}
TE_SEQS = {}
TE_ALIGNERS = {}

CHROMOSOMES = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
               "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
               "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

def print_alignment(hit, has_signal):
    print("        [{}] [{}/{} bases]\tread:{}-{}\t{}:{}-{}\tQ{}\t{}".format(
        "!!!" if has_signal else "   ",
        hit.mlen, hit.q_en - hit.q_st,
        hit.q_st, hit.q_en,
        hit.ctg, hit.r_st, hit.r_en,
        hit.mapq,
        hit.cigar_str)
    )

def te_sequence(seq, cut_site, slop=500):
    lo = max(0, cut_site - slop)
    hi = min(len(seq), cut_site + slop)
    return seq[lo:hi], cut_site - lo

def make_aligner(consensus, cut_site):
    r = next(SeqIO.parse(consensus, "fasta"))
    id = r.id
    seq, cs = te_sequence(r.seq, cut_site)

    cut_site_ref = f'cut-site_{consensus}'
    SeqIO.write([SeqRecord(Seq(seq), id=id)], cut_site_ref, "fasta")

    aligner = mp.Aligner(cut_site_ref,
                         k=9, w=15, best_n=10,
                         n_threads=1,
                         min_cnt=2)

    if not aligner:
        raise Exception("ERROR: failed to load/build index")

    TE_CUT_SITES[id] = cs
    TE_SEQS[id] = seq
    TE_ALIGNERS[id] = aligner

def make_aligners():
    make_aligner('l1.fa', 5932)
    make_aligner('aluya5.fa', 242)
    make_aligner('aluyb8.fa', 254)

def check_alignment_quality(alignment):
    # TODO figure out decent values for these
    # print((alignment.mapq, alignment.mlen))
    return (alignment.mapq > 1 and alignment.mlen > 40)

def check_alignment_read_location(alignment, query_len, is_5):
    # The alignment should be near the end of the read, allowing some slop for e.g. adapters.
    distance = alignment.q_st if is_5 else (query_len - alignment.q_en)

    # TODO configure threshold
    # print((is_5, alignment.q_st, alignment.q_en, query_len, distance))
    return distance < 80

def check_alignment_te_location(alignment, expected_cut_site):
    # The alignment should align near the cut site, but could be going the long OR the
    # short way through:
    #
    # L1Hs: -------------------------------X----------
    #                   ||||||||||||||||||||
    #                                      |||||||||||
    distance = min(abs(alignment.r_st - expected_cut_site),
                   abs(alignment.r_en - expected_cut_site))
    # print(distance)
    # c411b5e9-78da-44a8-8e53-759757439390
    return distance < 20

def check_alignment(alignment, expected_cut_site, query_len, is_5):
    return (check_alignment_quality(alignment) and
            check_alignment_read_location(alignment, query_len, is_5) and
            check_alignment_te_location(alignment, expected_cut_site))

class Signal:
    def __init__(self, read, is_5, te_id):
        self.is_5 = is_5
        self.read_id = read.query_name
        self.contig = read.reference_name
        self.te_id = te_id

        if not read.is_reverse:
            # ------------------------------------------------- ref
            #         rs               re
            #
            #         s5               s3
            # |||||||||[_______________]||||||||||>
            #         qs               qe
            if is_5:
                self.pos = read.reference_start
                self.clip_length = read.query_alignment_start
            else:
                self.clip_length = read.query_length - read.query_alignment_end
                self.pos = read.reference_end
        else:
            # ------------------------------------------------- ref
            #          rs               re
            #
            #          s3               s5
            # <|||||||||[_______________]|||||
            #          qs               qe
            if is_5:
                self.clip_length = read.query_length - read.query_alignment_end
                self.pos = read.reference_end
            else:
                self.clip_length = read.query_alignment_start
                self.pos = read.reference_start

        self.is_reference = self.clip_length < 60

    def __str__(self):
        return f'<{self.te_id} signal ({"5p" if self.is_5 else "3p"}) at {self.contig}:{self.pos} length {self.clip_length} ({self.read_id})>'

    def __repr__(self):
        return str(self)

def check_read(read, te_id, min_length=500):
    if read.query_length < min_length:
        return None, None

    l = min(read.query_length, 250)

    max_alignment_length = min(l, len(TE_SEQS[te_id]))
    # print()
    # print(f'    Checking {read.query_name} length {read.query_length} rev {read.is_reverse} {te_id} (cut site at {TE_CUT_SITES[te_id]} in TE, max alignment length {max_alignment_length})...')

    r_5 = read.query_sequence[:l]
    r_3 = read.query_sequence[-l:]

    if read.is_reverse:
        r_5, r_3 = reverse_complement(r_3), reverse_complement(r_5)

    aligner = TE_ALIGNERS[te_id]

    signal_5 = []
    hits_5 = list(aligner.map(r_5))
    if hits_5:
        # print(f"        5' end: {r_5}")
        for hit in hits_5:
            has_signal = check_alignment(hit, TE_CUT_SITES[te_id], l, True)
            # print_alignment(hit, has_signal)
            if has_signal:
                signal_5.append(hit)

    signal_3 = []
    hits_3 = list(aligner.map(r_3))
    if hits_3:
        # print(f"        3' end: {r_3}")
        for hit in hits_3:
            has_signal = check_alignment(hit, TE_CUT_SITES[te_id], l, False)
            # print_alignment(hit, has_signal)
            if has_signal:
                signal_3.append(hit)

    s5 = Signal(read, True, te_id)  if signal_5 else None
    s3 = Signal(read, False, te_id) if signal_3 else None

    return s5, s3

def cluster_alignment(db, signal, slop=15):
    if signal.contig not in db:
        db[signal.contig] = IntervalTree()

    tree = db[signal.contig]

    lo = max(0, signal.pos - slop)
    hi = (signal.pos + slop) # technically we should bound at chromosome end…

    existing = tree[lo:hi]
    if not existing:
        tree[lo:hi] = [signal]
    else:
        merged = [signal]
        for i in existing:
            merged.extend(i.data)
            tree.remove(i)
        tree[lo:hi] = merged

def is_representative_alignment(r):
    # TODO make this smarter in the future
    is_primary = not (r.is_supplementary or r.is_secondary)
    return (r.is_mapped
            and is_primary
            and r.reference_name in CHROMOSOMES
            and r.mapping_quality >= 20
            and r.get_tag('qs') > 10)

def run(alignment_path):
    te_ids = ["LINE1", "AluYa5", "AluYb8"]
    make_aligners()
    clustering_dbs = {te_id: {} for te_id in te_ids}

    i = 0
    with pysam.AlignmentFile(alignment_path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            i += 1
            if i % 1000 == 0:
                print('.', end='', flush=True)

            if is_representative_alignment(r):
                for te_id in te_ids:
                    s5, s3 = check_read(r, te_id)
                    if s5: cluster_alignment(clustering_dbs[te_id], s5)
                    if s3: cluster_alignment(clustering_dbs[te_id], s3)
    print()

    for te_id, db in clustering_dbs.items():
        with \
            open(f'nanopal-nonref-{te_id}.all.bed', 'w') as f_ins_all, \
            open(f'nanopal-nonref-{te_id}-multiple.bed', 'w') as f_ins_multiple, \
            open(f'nanopal-ref-{te_id}.all.bed', 'w') as f_ref_all, \
            open(f'nanopal-ref-{te_id}.multiple.bed', 'w') as f_ref_multiple:
            for contig, tree in db.items():
                    for i in tree:
                        lo = i.begin
                        hi = i.end
                        signals = i.data

                        ref_signals = [s for s in signals if s.is_reference]
                        ins_signals = [s for s in signals if not s.is_reference]

                        ref_support = len(ref_signals)
                        ins_support = len(ins_signals)

                        def _ids(sigs):
                            return ' '.join(signal.read_id for signal in sigs)

                        bed_line = f'{contig}\t{lo}\t{hi}\t{ref_support} {_ids(ref_signals)}'
                        if ref_support > 0:
                            print(bed_line, file=f_ref_all)
                        if ref_support > 1:
                            print(bed_line, file=f_ref_multiple)

                        bed_line = f'{contig}\t{lo}\t{hi}\t{ins_support} {_ids(ins_signals)}'
                        if ins_support > 0:
                            print(bed_line, file=f_ins_all)
                        if ins_support > 1:
                            print(bed_line, file=f_ins_multiple)

if __name__ == '__main__':
    run(sys.argv[1])
