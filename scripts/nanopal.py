#!/usr/bin/env python3

import collections
import csv
import math
import os
import sys

from Bio import SeqIO
from Bio.Seq import reverse_complement, complement

from ssw import AlignmentMgr

from intervaltree import IntervalTree

import pysam

# Maximum amount of softclip to allow for signals inside a reference element before we
# consider them to be non-reference insertions inside a reference element.
REF_SOFTCLIP = 40

# Maximum amount of softclip to allow for signals OUTSIDE a reference element before we
# consider them to be non-reference insertions.  If a signal has LESS than this as the
# median softclip we consider it an unannotated reference element.
NO_SOFTCLIP = 10

TE_ALIGNERS = {}

CHROMOSOMES = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
               "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
               "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}


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
    return str(long).upper(), str(short).upper()

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


CUT_SITES = {
    'LINE': 5918,
    'AluYa': 242,
    'AluYb': 254,
}

def make_te_aligners(te_id, consensus_path):
    r = next(SeqIO.parse(consensus_path, "fasta"))

    seq_lway, seq_sway = te_sequence(r.seq, CUT_SITES[te_id])

    lway_aligner = Aligner(te_id, 'long',  seq_lway)
    sway_aligner = Aligner(te_id, 'short', seq_sway)

    TE_ALIGNERS[te_id] = (lway_aligner, sway_aligner)

def check_alignment_quality(read_id, alignment, score_threshold=15):
    # ref (te) query
    score, (rs, re), (qs, qe) = alignment

    consensus_len = re - rs
    query_len = qe - qs

    # TODO figure out decent values for these
    ok = score >= score_threshold and consensus_len > 20 and query_len > 20

    # if not ok:
    #     print(f'excluding {read_id} for alignment quality', file=sys.stderr)

    return ok

def check_alignment_read_location(read_id, alignment, query_len, is_5):
    _, (_, _), (qs, _) = alignment

    # The alignment should be near the appropriate end of the read, allowing some slop
    # for e.g. adapters.
    ok = qs < 60

    # if not ok:
    #     print(f'excluding {read_id} for read location', file=sys.stderr)

    return ok

def check_alignment_te_location(read_id, alignment):
    _, (rs, _), (_, _) = alignment

    # The reference sequence was generated such that the cut site is always at the
    # start, so we should expect a qualifying alignment to start there.
    # TODO choose threshold
    ok = rs < 20

    # if not ok:
    #     print(f'excluding {read_id} for te location', file=sys.stderr)

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
        # print(f"no {te_id} alignment for {read_id} {is_5}", file=sys.stderr)
        return

    if not check_alignment(read_id, alignment, query_len, is_5):
        # print(f"failed quality check for {te_id} {read_id} {is_5}", file=sys.stderr)
        return

    return alignment

def check_read(read, te_id, end_length=140):
    read_id = read.query_name
    l = end_length

    if read.query_length < l:
        # print(f"excluding {read_id} for min read length", file=sys.stderr)
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

def is_valid_read(r, excluded_reads):
    # TODO make this smarter in the future
    return (r.query_name not in excluded_reads
            and r.get_tag('qs') >= 9.0)

def is_representative_alignment(r):
    # TODO make this smarter in the future
    is_primary = not (r.is_supplementary or r.is_secondary)
    return (and r.is_mapped
            and is_primary
            and r.reference_name in CHROMOSOMES
            and r.mapping_quality >= 20)

def median_sd(data):
    if not data:
        return 0, 0

    data = sorted(data)

    n = len(data)

    if n % 2 == 0:
        a = data[n // 2 - 1]
        b = data[n // 2]
        median = (a + b) / 2
    else:
        median = data[n // 2]

    mean = sum(data) / n
    sd = math.sqrt(sum(math.pow(x-mean, 2) for x in data) / n)

    return median, sd

def build_reference_db(reference_te_path):
    db = collections.defaultdict(IntervalTree)
    with open(reference_te_path) as f:
        # Lines look like:
        # chr1	100000000	100000637	L1M2	ref	-	6514	7153
        for line in f:
            contig, start, end, te_id, _, _, _, _ = line.split('\t')
            if contig in CHROMOSOMES:
                db[contig][int(start):int(end)] = te_id
    return db

def find_reference_elements(reference_db, contig, start, end):
    return sorted(x.data for x in reference_db[contig][start:end])

def output_results(te_id, output_dir, clustering_db, reference_db):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    count_ref, count_all, count_multi = 0, 0, 0

    output_dir = output_dir.rstrip('/')
    with open(f'{output_dir}/nanopal-{te_id}-all.bed', 'w') as f_all, \
         open(f'{output_dir}/nanopal-{te_id}-multiple.bed', 'w') as f_multiple, \
         open(f'{output_dir}/nanopal-{te_id}-reference.bed', 'w') as f_ref, \
         open(f'{output_dir}/nanopal-{te_id}.csv', 'w') as f_csv:
        w_csv = csv.writer(f_csv)
        w_csv.writerow(['te_id', 'contig', 'start', 'end', 'total_support',
                        'is_reference', 'reference_elements',
                        'long_way_support', 'long_way_median', 'long_way_sd',
                        'short_way_support', 'short_way_median', 'short_way_sd',
                        'long_way_read_ids', 'short_way_read_ids'])
        for contig, tree in clustering_db.items():
            for i in tree:
                lo = i.begin
                hi = i.end
                signals = i.data
                support = len(signals)

                refs = find_reference_elements(reference_db, contig, lo, hi)

                lway_signals = [s for s in signals if s.is_lway]
                sway_signals = [s for s in signals if not s.is_lway]

                lway_support = len(lway_signals)
                sway_support = len(sway_signals)

                def _ids(sigs):
                    return ' '.join(signal.read_id for signal in sigs)

                def _clips(sigs):
                    return [signal.clip_length for signal in sigs]

                def _refs(refs):
                    return ' '.join(refs)

                lway_read_ids = _ids(lway_signals)
                sway_read_ids = _ids(sway_signals)
                lway_median, lway_sd = median_sd(_clips(lway_signals))
                sway_median, sway_sd = median_sd(_clips(sway_signals))

                is_reference = refs and (lway_median <= REF_SOFTCLIP) and (sway_median <= REF_SOFTCLIP)
                no_soft_clip = lway_median <= NO_SOFTCLIP and sway_median <= NO_SOFTCLIP

                lway_desc = f'long {lway_support} m {lway_median:.1f} sd {lway_sd:.1f}'
                sway_desc = f'short {sway_support} m {sway_median:.1f} sd {sway_sd:.1f}'

                desc = f'{lway_desc} / {sway_desc} / {lway_read_ids} {sway_read_ids}'

                if is_reference:
                    desc = f'{_refs(refs)} / ' + desc
                elif no_soft_clip:
                    desc = f'unannotated / ' + desc

                mid = int((lo+hi)/2)
                thick_start = mid-1
                thick_end = mid+1

                bed_line = f'{contig}\t{lo}\t{hi}\t{desc}\t{support}\t.\t{thick_start}\t{thick_end}'

                if is_reference or no_soft_clip:
                    print(bed_line, file=f_ref)
                    count_ref += 1
                else:
                    print(bed_line, file=f_all)
                    count_all += 1

                    if support >= 2:
                        count_multi += 1
                        print(bed_line, file=f_multiple)

                w_csv.writerow([te_id, contig, lo, hi, support,
                                is_reference, _refs(refs),
                                lway_support, lway_median, lway_sd,
                                sway_support, sway_median, sway_sd,
                                lway_read_ids, sway_read_ids])

    with open(f'{output_dir}/summary.txt', 'w') as f:
        print(f'{count_ref} reference events', file=f)
        print(f'{count_all} non-reference events, {count_multi} with multiple read support', file=f)

def parse_minimera(minimera_path):
    excluded_reads = set()

    with open(minimera_path) as f:
        r = csv.reader(f)
        _ = next(r) # header
        # read-id,read-length,classification,monotony,foldback-point,mean-qscore,llqr,llqr-start,llqr-end,lq-total,processing-time-microsec

        for row in r:
            id = row[0]
            classification = row[2]

            if classification == "foldback":
                excluded_reads.add(id)

    return excluded_reads


def run(te_id, te_fasta_path, reference_te_path, alignment_path, minimera_path, output_dir):
    make_te_aligners(te_id, te_fasta_path)
    clustering_db = collections.defaultdict(IntervalTree)
    reference_db = build_reference_db(reference_te_path)
    excluded_reads = parse_minimera(minimera_path)

    i = 0
    with pysam.AlignmentFile(alignment_path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            i += 1
            if i % 10000 == 0:
                print(f'{i: 14,}', flush=True)

            # if i == 50000:
            #     break

            if is_valid_read(r, excluded_reads) and is_representative_alignment(r):
                for signal in check_read(r, te_id):
                    cluster_alignment(clustering_db, signal)
    print()

    output_results(te_id, output_dir, clustering_db, reference_db)

if __name__ == '__main__':
    run(*sys.argv[1:])
