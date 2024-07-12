#!/usr/bin/env bash

set -eu
set -x

palmer_reads="$1" # read.all.palmer.final.txt
ref_mei=$(readlink -f "$2") # hg38.RM.L1.ref
orig_mei=$(readlink -f "$3") # PALMER.NA12878.L1.txt

out_dir="$4"
out_summary="$5" # summary.final.txt

mkdir -p "$out_dir"
cd "$out_dir"

# We start with the results from Palmer and Nanopal's own BLAST queries:
#
# palmer_reads:
#                            1           2            3        4         5        6         7        8         9        10        11     12         13
#                            read_name   read_length  n_plus5  n_minus5  n_plus3  n_minus3  p_plus5  p_minus5  p_plus3  p_minus3  chr    start_pos  end_pos
# f5edeb8c-f225-4b47-b7ee-1c0df4dea038   172594       0        1         0        0         0        0         0        0         chr18  74826199   74999245
#
# I *think* that the chr/start/end will be NON 0 0 if the end of the read does
# not align to a position in the genome, which would be the case for
# a non-reference insertion.

# Extract the 5' position of the alignment of any MAPPED reads (because we're REMOVING the NON mapped).
#
#          chr start start+1 read_id
awk '{print $11, $12, $12+1, $1}' "$palmer_reads"  | grep -v "NON" > read.5.Q.txt

# Extract the 3' position of the alignment of any MAPPED reads (because we're REMOVING the NON mapped).
#
#           chr end   end+1  read_id
awk '{print $11, $13, $13+1, $1}' "$palmer_reads"  | grep -v "NON" > read.3.Q.txt

# 5' events -------------------------------------------------------------------

# Intersect the mapped reads' 5' positions with the reference MEIs, producing read.5.Q.ref.txt.
cp read.5.Q.txt Q.txt

# ref_mei:
# chr     start           end             type    ?       strand  ?       ?
# chr1    100000000       100000637       L1M2    ref     -       6514    7153
cp "${ref_mei}"  S.txt

intersect Q.txt S.txt > inter.txt
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt

collapse input_rm_cluster.txt > output_rm_cluster.txt
mv output_rm_cluster.txt read.5.Q.ref.txt

# Resulting read.5.Q.ref.txt
# chr3    121088939       121089040       1b684087-3b2a-492f-9dfb-e04180da03dc    L1HS
# chr12   37769514        37769615        1b685572-7ac1-42c9-866b-4ab091cf7c2e    NON
# chr11   50504827        50504928        1b685a87-e92a-443c-a998-596fe8e59490    NON
# chr4    111488479       111488580       1b685b54-06d3-43ad-9084-c598629c94dc    L1PA8A/L1PA13
# chr6    33425878        33425979        1b6864bd-8009-412b-9e66-0156a95fa391    NON
# chr6    54579950        54580051        1b6864dc-eeba-40c6-a70b-e9a792371bb6    L1ME2

# Now intersect the original input (reads' 5' positions) with P&P results, which look like this:
#
# 1      2          3          4                                                                                                       5                           6     7  8     9     10    11    12 13     14    15
# chr10  108693690  108693716  25.21.0.409836.0/1.9.9.53.53.21.0.0.cluster0_chr10_108693709_108693713_108693709_108693713_NA12878      TTAAAACAGCTCTCTTAAAGTCTCCA  898   0  5214  5229  6039  6053  +  L1HS   5214  6053
# chr10  122695701  122695713  26.8.0.454148.0/1.13.13.69.69.12.0.0.cluster0_chr10_122695701_122695716_122695701_122695716_NA12878     GAAAAACTTAAAGGGGAGGAGCCAAG  6104  0  2     5     6036  6054  +  L1HS   2     6054
# chr10  25418754   25418778   34.7.0.604201.0/1.9.9.0.0.20.0.0.cluster0_chr10_25418758_25418779_25418758_25418779_NA12878             GGCTCCTCCTCCCGAAAATAATCTTT  6045  0  1     11    6040  6051  -  L1HS   1     6051
cp "${orig_mei}" S.txt
intersect Q.txt S.txt > inter.txt

# In rare cases, two or more P&P targets might be within 50bp one of our reads'
# ends, which `intersect` will output as separate lines.  We need to deduplicate
# these before we pass them to `paste` below, otherwise the lines will get out
# of sync.  Note that the other file (the one with the reference events, as
# opposed to the P&P events) is already deduplicated by `collapse`.

mv inter.txt prededupe.read.5.Q.P.txt

# This is deduping on field FIVE, beacause -f 4 means SKIP the first four
# fields.  So it's deduping based on read ID.
cat prededupe.read.5.Q.P.txt \
    | awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' \
    | awk '{print $1,$2,$3,$5,$4}' \
    | sort -n | uniq -f 4 \
    | awk '{print $1,$2,$3,$5,$4}' \
    > read.5.Q.P.txt

# So now we have:
#
#     read.5.Q.txt:     The 5' positions of all alignments of mapped reads.
#     read.5.Q.ref.txt: The 5' positions of all alignments of mapped reads, with any reference events noted.
#     read.5.Q.P.txt:   The 5' positions of all alignments of mapped reads, with any P&P events noted.

# 3' events -------------------------------------------------------------------

# Intersect the mapped reads' 3' positions with the reference MEIs, producing read.3.Q.ref.txt.
cp read.3.Q.txt Q.txt

cp "${ref_mei}"  S.txt

intersect Q.txt S.txt > inter.txt
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt

collapse input_rm_cluster.txt > output_rm_cluster.txt
mv output_rm_cluster.txt read.3.Q.ref.txt

# Result format is the same as read.5.Q.txt above

# Now intersect the original input (reads' 3' positions) with P&P results.
cp "${orig_mei}"  S.txt
intersect Q.txt S.txt > inter.txt

# And dedupe as before.
mv inter.txt prededupe.read.3.Q.P.txt

cat prededupe.read.3.Q.P.txt \
    | awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' \
    | awk '{print $1,$2,$3,$5,$4}' \
    | sort -n | uniq -f 4 \
    | awk '{print $1,$2,$3,$5,$4}' \
    > read.3.Q.P.txt

# So now we have:
#
#     read.3.Q.txt:     The 3' positions of all alignments of mapped reads.
#     read.3.Q.ref.txt: The 3' positions of all alignments of mapped reads, with any reference events noted.
#     read.3.Q.P.txt:   The 3' positions of all alignments of mapped reads, with any P&P events noted.

# Combine ---------------------------------------------------------------------

# Sort all the files by read ID.
sort -k 4 read.5.Q.P.txt   > 5.P.inter.txt
sort -k 4 read.3.Q.P.txt   > 3.P.inter.txt
sort -k 4 read.5.Q.ref.txt > 5.ref.inter.txt
sort -k 4 read.3.Q.ref.txt > 3.ref.inter.txt

paste 5.P.inter.txt 3.P.inter.txt 5.ref.inter.txt 3.ref.inter.txt | awk '{print $4,$5,$10,$15,$20}' > read.RM.txt

join -1 1 -2 1 \
     -a 1 \
     -o auto -e NA \
     --check-order \
     <(cat "$palmer_reads" | sort) \
     <(cat read.RM.txt | sort) \
     > "$out_summary"

