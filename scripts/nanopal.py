#!/usr/bin/env python3

import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord

from intervaltree import IntervalTree

import pysam
import mappy as mp

import numpy as np

TE_CUT_SITES = {}
TE_SEQS = {}
TE_ALIGNERS = {}

CHROMOSOMES = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
               "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
               "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

SW_UP = 1
SW_LEFT = 2
SW_DIAG = 3

SW_MATCH = -1
SW_GAP = -1

def smith_waterman(seq1, seq2):
    rows = len(seq2) + 1
    cols = len(seq1) + 1

    grid = np.zeros((rows, cols), dtype=np.int32)
    prev = np.zeros((rows, cols), dtype=np.uint8)

    best_score = 0
    best_row, best_col = -1, -1

    for row in range(1, rows):
        ch2 = seq2[row-1]
        for col in range(1, cols):
            ch1 = seq1[col-1]

            score_match = grid[row-1][col-1] + (1 if ch1 == ch2 else -1)
            score_gap1 = grid[row][col-1] + SW_GAP
            score_gap2 = grid[row-1][col] + SW_GAP
            best = max(score_match, score_gap1, score_gap2)

            grid[row][col] = best

            if best == score_match:
                prev[row][col] = SW_DIAG
            elif best == score_gap1:
                prev[row][col] = SW_LEFT
            else:
                prev[row][col] = SW_UP

            if best > best_score:
                best_score = best
                best_row = row
                best_col = col

    row, col = best_row, best_col
    e2, e1 = row, col
    s2, s1 = None, None
    # path = []

    while grid[row][col] > 0:
        s2, s1 = row, col
        # path.append((row, col))
        if prev[row][col] == SW_UP:
            row -= 1
        elif prev[row][col] == SW_LEFT:
            col -= 1
        else:
            row -= 1
            col -= 1

    # if not path:
    #     return [], None, None
    # else:
    #     return path[::-1], (s1-1, e1-1), (s2-1, e2-1)

    if s2:
        return (s1-1, e1-1), (s2-1, e2-1)
    else:
        return None, None


def print_smith_waterman(seq1, seq2):
    path, coords1, coords2 = smith_waterman(seq1, seq2)
    print(f'{seq1=}')
    print(f'{seq2=}')
    print()
    print(f'seq1[{coords1}], seq2[{coords2}]')
    print()
    chars1 = []
    charsM = []
    chars2 = []
    p1, p2 = -1, -1
    for i2, i1 in path:
        # print(((i1, i2), (p1, p2)))
        if i1 == p1 or i2 == p2:
            charsM.append(' ')
        else:
            charsM.append('|' if seq1[i1-1] == seq2[i2-1] else 'x')

        if i1 > p1:
            p1 = i1
            chars1.append(seq1[i1-1])
        else:
            chars1.append('-')

        if i2 > p2:
            p2 = i2
            chars2.append(seq2[i2-1])
        else:
            chars2.append('-')

    # print(chars1)
    # print(charsM)
    # print(chars2)
    print (
        ''.join(chars1) + '\n' +
        ''.join(charsM) + '\n' +
        ''.join(chars2)
    )


def print_alignment(hit, has_signal):
    print("        [{}] [{}/{} bases]\tread:{}-{}\t{}:{}-{}\tQ{}\t{}".format(
        "!!!" if has_signal else "   ",
        hit.mlen, hit.q_en - hit.q_st,
        hit.q_st, hit.q_en,
        hit.ctg, hit.r_st, hit.r_en,
        hit.mapq,
        hit.cigar_str)
    )

def te_sequence(seq, cut_site, slop=100):
    # Split the whole TE reference sequence seq at `cut_site` and return the reference
    # sequence for the upstream and downstream way, respectively.  Both returned
    # reference sequences will have the cut site at index 0, so good alignments should
    # start there.

    # Using AluYa5 as an example:
    #
    #                                                                               cut site
    #                                                                                  |
    #                                                                                  V   downstream -->
    # 5'-…………CTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCA-3' 5'-CTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-3'
    # 3'-…………GAACCCTCCGACTCCGTCCTCTTACCGCACTTGGGCCCTCCGCCTCGAACGTCACTCGGCTCTAGGGCGGT-5' 3'-GACGTGAGGTCGGACCCGCTGTCTCGCTCTGAGGCAGAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-5'
    #                                                                   <-- upstream

    lo = max(0, cut_site - slop)
    hi = min(len(seq), cut_site + slop)

    us = reverse_complement(seq[lo:cut_site])
    ds = seq[cut_site:hi]

    return us, ds

class Aligner:
    def __init__(self, consensus_seq):
        self.consensus_seq = consensus_seq



def make_te_aligners(consensus, cut_site):
    r = next(SeqIO.parse(consensus, "fasta"))
    id = r.id

    seq_us, seq_ds = te_sequence(r.seq, cut_site)

    us_aligner = Aligner(seq_us)
    ds_aligner = Aligner(seq_ds)

    # us_ref  = f'us_{consensus}'
    # ds_ref = f'ds_{consensus}'
    # SeqIO.write([SeqRecord(Seq(seq_us),  id=id)], us_ref,  "fasta")
    # SeqIO.write([SeqRecord(Seq(seq_ds), id=id)], ds_ref, "fasta")

    # us_aligner  = mp.Aligner(us_ref,  preset="map-ont", k=7, w=8, best_n=10, n_threads=1, min_cnt=2)
    # ds_aligner = mp.Aligner(ds_ref, preset="map-ont", k=7, w=8, best_n=10, n_threads=1, min_cnt=2)

    # if not us_aligner:
    #     raise Exception("ERROR: failed to load/build index")

    # if not ds_aligner:
    #     raise Exception("ERROR: failed to load/build index")

    TE_ALIGNERS[id] = (us_aligner, ds_aligner)

def make_aligners():
    make_te_aligners('l1.fa', 5932)
    make_te_aligners('aluya5.fa', 242)
    make_te_aligners('aluyb8.fa', 254)

def check_alignment_quality(read_id, alignment):
    # TODO figure out decent values for these
    ok = alignment.mapq > 1 and alignment.mlen > 30

    if not ok:
        print(f'excluding {read_id} for alignment quality', file=sys.stderr)

    return ok

def check_alignment_read_location(read_id, alignment, query_len, is_5):
    # The alignment should be near the appropriate end of the read, allowing some slop
    # for e.g. adapters.
    distance = alignment.q_st if is_5 else (query_len - alignment.q_en)

    # TODO configure slop
    ok = distance < 80

    if not ok:
        print(f'excluding {read_id} for read location', file=sys.stderr)

    return ok

def check_alignment_te_location(read_id, alignment):
    # The reference sequence was generated such that the cut site is always at the
    # start, so we should expect a qualifying alignment to start there.
    ok = alignment.r_st < 20

    if not ok:
        print(f'excluding {read_id} for te location', file=sys.stderr)

    return ok


def check_alignment(read_id, alignment, query_len, is_5):
    return (check_alignment_quality(read_id, alignment) and
            check_alignment_read_location(read_id, alignment, query_len, is_5) and
            check_alignment_te_location(read_id, alignment))

class Signal:
    def __init__(self, read, te_id, is_5, is_us):
        self.is_5 = is_5
        self.read_id = read.query_name
        self.contig = read.reference_name
        self.te_id = te_id
        self.is_us = is_us

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
        return f'<{self.te_id} {"upstream" if self.is_us else "downstream"} signal ({"5p" if self.is_5 else "3p"}) at {self.contig}:{self.pos} length {self.clip_length} ({self.read_id})>'

    def __repr__(self):
        return str(self)

def alignments(read_id, te_id, aligner, sequence, query_len, is_5):
    hits = list(aligner.map(sequence))
    # hits = list(aligner.map(sequence))

    # if read_id == 'b5a85522-300e-4816-bc6e-ddb557e040e7':
    #     print((read_id, te_id, is_5))
    #     print(hits)
    #     if hits:
    #         print(hits[0])
    #         print(hits[0].mapq)

    # if read_id == '533eb45e-c548-4826-8f06-273597e52bb1':
    #     print((read_id, te_id, is_5))
    #     print(sequence)
    #     print(hits)
    #     print()
    if not hits:
        print(f"no {te_id} hits for {read_id} {is_5}", file=sys.stderr)

    return [hit for hit in hits if check_alignment(read_id, hit, query_len, True)]

def check_read(read, te_id, min_length=150):
    read_id = read.query_name

    if read.query_length < min_length:
        print(f"excluding {read_id} for min read length", file=sys.stderr)
        return []

    l = min(read.query_length, 120)

    r_5 = read.query_sequence[:l]
    r_3 = read.query_sequence[-l:]

    if read.is_reverse:
        r_5, r_3 = reverse_complement(r_3), reverse_complement(r_5)

    # r_5 = r_5[35:-30]

    us_aligner, ds_aligner = TE_ALIGNERS[te_id]

    us_signal_5 = alignments(read_id, te_id, us_aligner, r_5, l, True)
    ds_signal_5 = alignments(read_id, te_id, ds_aligner, r_5, l, True)
    us_signal_3 = alignments(read_id, te_id, us_aligner, r_3, l, False)
    ds_signal_3 = alignments(read_id, te_id, ds_aligner, r_3, l, False)

    signals = []
    if us_signal_5: signals.append(Signal(read, te_id, is_5=True,  is_us=True))
    if us_signal_3: signals.append(Signal(read, te_id, is_5=False, is_us=True))
    if ds_signal_5: signals.append(Signal(read, te_id, is_5=True,  is_us=False))
    if ds_signal_3: signals.append(Signal(read, te_id, is_5=False, is_us=False))

    return signals

def cluster_alignment(db, signal, slop=15):
    if signal.contig not in db:
        db[signal.contig] = IntervalTree()

    tree = db[signal.contig]

    lo = max(0, signal.pos - slop)
    hi = (signal.pos + slop) # TODO technically we should bound at chromosome end…

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
            and r.mapping_quality >= 10
            and r.get_tag('qs') > 10)

def run(alignment_path, output_dir):
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
                    for signal in check_read(r, te_id):
                        cluster_alignment(clustering_dbs[te_id], signal)
    print()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_dir = output_dir.rstrip('/')
    for te_id, db in clustering_dbs.items():
        with \
            open(f'{output_dir}/nanopal-{te_id}-all.bed', 'w') as f_all, \
            open(f'{output_dir}/nanopal-{te_id}-multiple.bed', 'w') as f_multiple:
            for contig, tree in db.items():
                    for i in tree:
                        lo = i.begin
                        hi = i.end
                        signals = i.data
                        support = len(signals)

                        us_signals  = [s for s in signals if s.is_us]
                        ds_signals = [s for s in signals if not s.is_us]

                        us_support  = len(us_signals)
                        ds_support = len(ds_signals)

                        def _ids(sigs):
                            return ' '.join(signal.read_id for signal in sigs)

                        desc = f'us {us_support} / ds {ds_support} / upstream reads / {_ids(long_signals)} / downstream reads / {_ids(ds_signals)}'

                        bed_line = f'{contig}\t{lo}\t{hi}\t{desc}'

                        print(bed_line, file=f_all)
                        if support >= 2:
                            print(bed_line, file=f_multiple)

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])
