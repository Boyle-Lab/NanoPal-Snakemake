#!/usr/bin/env python3

import csv
import math
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement, complement
from Bio.SeqRecord import SeqRecord

from ssw import AlignmentMgr

from intervaltree import IntervalTree

import pysam

import numpy as np

TE_CUT_SITES = {}
TE_SEQS = {}
TE_ALIGNERS = {}

CHROMOSOMES = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
               "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
               "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}


def smith_waterman(seq1, seq2, print_alignment=False):
    SW_NA = 0
    SW_UP = 1
    SW_LEFT = 2
    SW_DIAG = 3

    SW_MATCH = 1
    SW_MISMATCH = -2
    SW_GAP = -2

    rows = len(seq2) + 1
    cols = len(seq1) + 1

    grid = np.zeros((rows, cols), dtype=np.int32)
    prev = np.zeros((rows, cols), dtype=np.uint8)

    best_score = 0
    best_row, best_col = -1, -1

    for row, ch2 in enumerate(seq2, start=1):
        for col, ch1 in enumerate(seq1, start=1):
            score_match = grid[row-1][col-1] + (SW_MATCH if ch1 == ch2 else SW_MISMATCH)
            score_gap1 = grid[row][col-1] + SW_GAP
            score_gap2 = grid[row-1][col] + SW_GAP

            # I kid you not, this max() call was taking 30% of our processing time when
            # I profiled it.
            # best = max(score_match, score_gap1, score_gap2, 0)

            if score_match <= 0 and score_gap1 <= 0 and score_gap2 <= 0:
                best = 0
            elif score_match >= score_gap1 and score_match >= score_gap2:
                best = score_match
                prev[row][col] = SW_DIAG
            elif score_gap1 >= score_gap2:
                best = score_gap1
                prev[row][col] = SW_LEFT
            else:
                best = score_gap2
                prev[row][col] = SW_UP

            grid[row][col] = best

            if best > best_score:
                best_score = best
                best_row = row
                best_col = col

    row, col = best_row, best_col
    e2, e1 = row, col
    s2, s1 = None, None
    if print_alignment:
        path = []

    while grid[row][col] > 0:
        s2, s1 = row, col

        if print_alignment:
            path.append((row, col))

        if prev[row][col] == SW_UP:
            row -= 1
        elif prev[row][col] == SW_LEFT:
            col -= 1
        else:
            row -= 1
            col -= 1

    if print_alignment:
        if not path:
            print("No alignment found.")
        else:
            print_smith_waterman(seq1, seq2, path[::-1], (s1-1, e1-1), (s2-1, e2-1))

    if s2:
        return best_score, (s1-1, e1-1), (s2-1, e2-1)
    else:
        return None


def print_smith_waterman(seq1, seq2, path, coords1, coords2):
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

    print (
        ''.join(chars1) + '\n' +
        ''.join(charsM) + '\n' +
        ''.join(chars2)
    )


def te_sequence(seq, cut_site, slop=100):
    # Split the whole TE reference sequence seq at `cut_site` and return the reference
    # sequence for the long and short way, respectively.  Both returned reference
    # sequences will have the cut site at index 0, so good alignments should start
    # there.
    #
    # Using AluYa5 as an example:
    #
    #                                                                               cut site
    #                                                                                  |
    #                                                                                  V   short -->
    # 5'-…………CTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCA-3' 5'-CTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-3'
    # 3'-…………GAACCCTCCGACTCCGTCCTCTTACCGCACTTGGGCCCTCCGCCTCGAACGTCACTCGGCTCTAGGGCGGT-5' 3'-GACGTGAGGTCGGACCCGCTGTCTCGCTCTGAGGCAGAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-5'
    #                                                                       <-- long

    lo = max(0, cut_site - slop)
    hi = min(len(seq), cut_site + slop)

    long = reverse_complement(seq[lo:cut_site])
    short = seq[cut_site:hi]

    # Need to str() here to remove the Seq wrapper, because using Seq objects in
    # Smith-Waterman dominates the CPU time with __getitem__() calls.
    return str(long), str(short)

class Aligner:
    def __init__(self, te_id, direction, consensus_seq):
        self.te_id = te_id
        self.direction = direction
        self.consensus_seq = consensus_seq

        self.align_mgr_5 = AlignmentMgr(match_score=1, mismatch_penalty=2)
        self.align_mgr_5.set_reference(consensus_seq)

        self.align_mgr_3 = AlignmentMgr(match_score=1, mismatch_penalty=2)
        self.align_mgr_3.set_reference(complement(consensus_seq))

    def align(self, query_seq, is_5, print_alignment=False):
        m = self.align_mgr_5 if is_5 else self.align_mgr_3
        m.set_read(query_seq)

        a = m.align(gap_open=2, gap_extension=1)

        if print_alignment:
            print(f'Alignment of {"5" if is_5 else "3"}\' end of read for {self.te_id} ({self.direction}-way):')
            print(f'    {a}')
            print()

        if a.reference_start == -1:
            return None

        return a.optimal_score, (a.reference_start, a.reference_end), (a.read_start, a.read_end)

def make_te_aligners(consensus, cut_site):
    r = next(SeqIO.parse(consensus, "fasta"))
    id = r.id

    seq_lway, seq_sway = te_sequence(r.seq, cut_site)

    lway_aligner = Aligner(id, 'long',  seq_lway)
    sway_aligner = Aligner(id, 'short', seq_sway)

    TE_ALIGNERS[id] = (lway_aligner, sway_aligner)

def make_aligners():
    make_te_aligners('l1.fa', 5918)
    make_te_aligners('aluya5.fa', 242)
    make_te_aligners('aluyb8.fa', 254)

def check_alignment_quality(read_id, alignment, score_threshold=15):
    # ref (te) query
    score, (rs, re), (qs, qe) = alignment

    consensus_len = re - rs
    query_len = qe - qs

    # TODO figure out decent values for these
    ok = score >= score_threshold and consensus_len > 20 and query_len > 20

    if not ok:
        print(f'excluding {read_id} for alignment quality', file=sys.stderr)

    return ok

def check_alignment_read_location(read_id, alignment, query_len, is_5):
    _, (_, _), (qs, _) = alignment

    # The alignment should be near the appropriate end of the read, allowing some slop
    # for e.g. adapters.
    ok = qs < 60

    if not ok:
        print(f'excluding {read_id} for read location', file=sys.stderr)

    return ok

def check_alignment_te_location(read_id, alignment):
    _, (rs, _), (_, _) = alignment

    # The reference sequence was generated such that the cut site is always at the
    # start, so we should expect a qualifying alignment to start there.
    # TODO choose threshold
    ok = rs < 20

    if not ok:
        print(f'excluding {read_id} for te location', file=sys.stderr)

    return ok

def check_alignment(read_id, alignment, query_len, is_5):
    return (check_alignment_quality(read_id, alignment) and
            check_alignment_read_location(read_id, alignment, query_len, is_5) and
            check_alignment_te_location(read_id, alignment))

class Signal:
    def __init__(self, read, te_id, is_5, is_lway):
        self.is_5 = is_5
        self.read_id = read.query_name
        self.contig = read.reference_name
        self.te_id = te_id
        self.is_lway = is_lway

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
        return f'<{self.te_id} {"long" if self.is_lway else "short"}-way signal ({"5p" if self.is_5 else "3p"}) at {self.contig}:{self.pos} clip-length {self.clip_length} ({self.read_id})>'

    def __repr__(self):
        return str(self)

def alignment(read_id, te_id, aligner, sequence, query_len, is_5):
    debug = read_id == 'ef44d580-e37b-460d-b2f0-566ae9e9e8aa'
    # debug = False
    alignment = aligner.align(sequence, is_5, print_alignment=debug)

    if not alignment:
        print(f"no {te_id} alignment for {read_id} {is_5}", file=sys.stderr)
        return

    if not check_alignment(read_id, alignment, query_len, is_5):
        # print(f"failed quality check for {te_id} {read_id} {is_5}", file=sys.stderr)
        return

    return alignment

def check_read(read, te_id, end_length=140):
    read_id = read.query_name
    l = end_length

    if read.query_length < l:
        print(f"excluding {read_id} for min read length", file=sys.stderr)
        return []

    r_5 = read.query_sequence[:l]
    r_3 = read.query_sequence[-l:]

    if read.is_reverse:
        r_5, r_3 = reverse_complement(r_3), reverse_complement(r_5)

    # We've structured the TE references to start with the cut site and proceed
    # left-to-right, so we need to reverse the 3' end of the read to match this
    # convention.
    r_3 = r_3[::-1]

    lway_aligner, sway_aligner = TE_ALIGNERS[te_id]

    lway_signal_5 = alignment(read_id, te_id, lway_aligner, r_5, l, True)
    sway_signal_5 = alignment(read_id, te_id, sway_aligner, r_5, l, True)
    lway_signal_3 = alignment(read_id, te_id, lway_aligner, r_3, l, False)
    sway_signal_3 = alignment(read_id, te_id, sway_aligner, r_3, l, False)

    signals = []
    if lway_signal_5: signals.append(Signal(read, te_id, is_5=True,  is_lway=True))
    if lway_signal_3: signals.append(Signal(read, te_id, is_5=False, is_lway=True))
    if sway_signal_5: signals.append(Signal(read, te_id, is_5=True,  is_lway=False))
    if sway_signal_3: signals.append(Signal(read, te_id, is_5=False, is_lway=False))

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
            and r.mapping_quality >= 20
            and r.get_tag('qs') >= 9.0)

def mean_sd(data):
    if not data:
        return 0, 0

    n = len(data)
    mean = sum(data) / n
    sd = math.sqrt(sum(math.pow(x-mean, 2) for x in data) / n)
    return mean, sd

def output_results(output_dir, clustering_dbs):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_dir = output_dir.rstrip('/')
    for te_id in clustering_dbs:
        with open(f'{output_dir}/nanopal-{te_id}-all.bed', 'w') as f_all, \
             open(f'{output_dir}/nanopal-{te_id}-multiple.bed', 'w') as f_multiple, \
             open(f'{output_dir}/nanopal-{te_id}.csv', 'w') as f_csv:
            w_csv = csv.writer(f_csv)
            w_csv.writerow(['te_id', 'contig', 'start', 'end',
                            'total_support',
                            'long_way_support', 'long_way_mean', 'long_way_sd',
                            'short_way_support', 'short_way_mean', 'short_way_sd',
                            'long_way_read_ids', 'short_way_read_ids'])
            for contig, tree in clustering_dbs[te_id].items():
                for i in tree:
                    lo = i.begin
                    hi = i.end
                    signals = i.data
                    support = len(signals)

                    lway_signals = [s for s in signals if s.is_lway]
                    sway_signals = [s for s in signals if not s.is_lway]

                    lway_support = len(lway_signals)
                    sway_support = len(sway_signals)

                    def _ids(sigs):
                        return ' '.join(signal.read_id for signal in sigs)

                    def _clips(sigs):
                        return [signal.clip_length for signal in sigs]

                    lway_read_ids = _ids(lway_signals)
                    sway_read_ids = _ids(sway_signals)
                    lway_mean, lway_sd = mean_sd(_clips(lway_signals))
                    sway_mean, sway_sd = mean_sd(_clips(sway_signals))

                    lway_desc = f'long {lway_support} m {lway_mean:.1f} sd {lway_sd:.1f}'
                    sway_desc = f'short {sway_support} m {sway_mean:.1f} sd {sway_sd:.1f}'

                    desc = f'{lway_desc} / {sway_desc} / {lway_read_ids} {sway_read_ids}'

                    mid = int((lo+hi)/2)
                    thick_start = mid-1
                    thick_end = mid+1

                    bed_line = f'{contig}\t{lo}\t{hi}\t{desc}\t{support}\t.\t{thick_start}\t{thick_end}'

                    print(bed_line, file=f_all)
                    if support >= 2:
                        print(bed_line, file=f_multiple)

                    w_csv.writerow([te_id, contig, lo, hi, support,
                                    lway_support, lway_mean, lway_sd,
                                    sway_support, sway_mean, sway_sd,
                                    lway_read_ids, sway_read_ids])


def run(te_id, alignment_path, output_dir):
    # te_ids = ["LINE1", "AluYa5", "AluYb8"]
    te_ids = [te_id]
    make_aligners()
    clustering_dbs = {te_id: {} for te_id in te_ids}

    i = 0
    with pysam.AlignmentFile(alignment_path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            i += 1
            if i % 10000 == 0:
                print(i, flush=True)

            # if i == 50000:
            #     break

            if is_representative_alignment(r):
                for te_id in te_ids:
                    for signal in check_read(r, te_id):
                        cluster_alignment(clustering_dbs[te_id], signal)
    print()

    output_results(output_dir, clustering_dbs)

if __name__ == '__main__':
    te_id = sys.argv[1]
    bam_path = sys.argv[2]
    results_path = sys.argv[3]
    run(te_id, bam_path, results_path)
