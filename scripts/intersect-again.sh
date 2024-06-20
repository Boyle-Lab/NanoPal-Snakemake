#!/usr/bin/env bash

set -eu
set -x

bam="$1" # Nanopore.sorted.bam
ref_mei=$(readlink -f "$2") # hg38.RM.L1.ref
pp_mei=$(readlink -f "$3") # L1.inter.fi
in_summary="$4" # summary.final.txt
valid_read_ids="$5" # RC.all.list
mei="$6" # LINE

out_dir="$7"
out_summary="$8" # summary.final.2.txt
out_result_log="$9" # formerly stdout

mkdir -p "$out_dir"
cd "$out_dir"

# $ref_mei contains the locations of reference MEI of various subtypes, e.g.:
#
#     chr1    100005056       100006760       L1M2    ref     -       1502    3236
#     chr1    100012520       100012872       L1ME3G  ref     +       165     481
#     chr1    100013062       100013131       L1ME3G  ref     +       519     592
#     chr1    100031348       100031482       L1MC4   ref     -       2451    2588
#     chr1    100031717       100032555       L1MC4   ref     -       6950    7495
#
# Here we split apart the file into N+1 subsets: the N different subtypes we're
# specifically interested, and one more for all the others.
#
# TODO this is going to be janky for Alus, because it's not as simple to grep
# them out (e.g. we want AluYa, AluYb, AluY[^ab], and other).  But could we
# split this step entirely into a separate offline generation step?  There's
# nothing in the reference that's changing at runtime, right?
#
# For now I'm just going to hardcode the stuff in here to move forward.
if test "$mei" = "LINE"; then
    cat "$ref_mei" | grep    L1PA > ref.L1PA
    cat "$ref_mei" | grep    L1HS > ref.L1HS
    cat "$ref_mei" | grep -v L1HS| grep -v L1PA | grep L1 > ref.L1
elif test \( "$mei" = "AluYa" \) -o \( "$mei" = "AluYb" \); then
    cat "$ref_mei" | grep AluYa > ref.AluYa5
    cat "$ref_mei" | grep AluYb > ref.AluYb8
    cat "$ref_mei" | grep -v AluYb | grep -v AluYa | grep AluY > ref.AluY
    cat "$ref_mei" | grep -v AluYb | grep -v AluYa | grep -v AluY | grep Alu > ref.Alu
else
    echo Unknown MEI type "$mei", expected LINE, AluYa, or AluYb.
    exit 1
fi

# Tack on an extra column to the end of the summary file, with 1 if the read ID
# is in the list of valid mapped reads, 0 if it's not.
join -1 1 -2 1 \
     -a 1 \
     -o auto -e "0" \
     --check-order \
     <(sort -k1 "$in_summary") \
     <(sort -k1 "$valid_read_ids" | sed -e 's/$/ 1/') \
     > "$out_summary"

# TODO Refactor the rest of this into separate chunks.

# Checking for signal on single end only
SIGNAL1=$(awk '(($3+$4)!=0 && ($5+$6)==0) || (($3+$4)==0 && ($5+$6)!=0)' "$out_summary" | wc -l)

# Checking for signal on both ends
SIGNAL2=$(awk '($3+$4)!=0 && ($5+$6)!=0' "$out_summary" | wc -l)

# Checking for no signal
FAIL=$(awk '($3+$4+$5+$6)==0' "$out_summary" | wc -l)

# non-reference events
cat "$out_summary" | awk '$11=="NON"' > summary.final.unmap.read.txt

# reference event, nanopal got no signal, palmer did get signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)==0&&($7+$8+$9+$10)!=0' > summary.final.odd.read.txt

# reference event, nanopal got no signal, palmer got no signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6+$7+$8+$9+$10)==0' > summary.final.no.L1.read.txt

# reference event, nanopal got signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0' > summary.final.all.L1.read.txt

# Both palmer and nanopal have signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0&&($7+$8+$9+$10)!=0' > summary.final.PALMER.read.txt

# only nanopal has signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0&&($7+$8+$9+$10)==0' > summary.final.ref.L1.read.txt


# Getting chr, start, clusteringinfo, validreadbit, 5/3, ± for reads where palmer and nanopal both have signal
cat summary.final.PALMER.read.txt | awk '{
                                                        # chr, start, clusteringinfo, valid_mapped_read_bit, 5 ±
    if      (($7+$8)  > 0 && ($9+$10) == 0 && $7>0) {print $11,$12,$14,$18,"5","+"}
    else if (($7+$8)  > 0 && ($9+$10) == 0 && $8>0) {print $11,$12,$14,$18,"5","-"}
                                                        # chr, end, clusteringinfo, valid_mapped_read_bit, 3 ±
    else if (($7+$8) == 0 && ($9+$10)  > 0 && $9>0) {print $11,$13,$15,$18,"3","+"}
    else if (($7+$8) == 0 && ($9+$10)  > 0 && $10>0) {print $11,$13,$15,$18,"3","-"}

    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $7>0&&$9>0) {print  $11,$12,$14,$18,"5","+""\n"$11,$13,$15,$18,"3","+"}
    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $7>0&&$10>0) {print $11,$12,$14,$18,"5","+""\n"$11,$13,$15,$18,"3","-"}
    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $8>0&&$9>0) {print  $11,$12,$14,$18,"5","-""\n"$11,$13,$15,$18,"3","+"}
    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $8>0&&$10>0) {print $11,$12,$14,$18,"5","-""\n"$11,$13,$15,$18,"3","-"}
}' | sort -k 2 -n | sort -k 1 > capture.loci.palmer
# Note: this appears to be trying to sort by two keys, but since it doesn't pass
# --stable sort will do a last-resort comparison, so this is actually only
# sorted by field 1 to end of line (but NOT numeric for field 2).  The proper
# way to do what I think this is tryign to do is: sort -k 1,1 -k 2,2n

# Getting chr, start, clusteringinfo, validreadbit, 5/3, for a mix of palmer/nanopal signals
cat summary.final.PALMER.read.txt | awk '{
    if      (($7+$8)>0&&($9+$10)==0&&$18==0&&($5+$6)>0) {print $11,$13,$17,$18,"3"}
    else if (($7+$8)>0&&($9+$10)==0&&$18==1&&($3+$4)>0) {print $11,$13,$17,$18,"3"}
    else if (($7+$8)==0&&($9+$10)>0&&$18==0&&($3+$4)>0) {print $11,$12,$16,$18,"5"}
    else if (($7+$8)==0&&($9+$10)>0&&$18==1&&($5+$6)>0) {print $11,$12,$16,$18,"5"}
}' | sort -k 2 -n | sort -k 1 > capture.loci.ref.add

# Getting ... for when only nanopal has signal
cat summary.final.ref.L1.read.txt | awk '{
    if      (($3+$4)  > 0 && ($5+$6) == 0 && $18 == 0) {print $11,$12,$16,$18,"5"}
    else if (($3+$4)  > 0 && ($5+$6) == 0 && $18 == 1) {print $11,$13,$17,$18,"3"}
    else if (($3+$4) == 0 && ($5+$6)  > 0 && $18 == 0) {print $11,$13,$17,$18,"3"}
    else if (($3+$4) == 0 && ($5+$6)  > 0 && $18 == 1) {print $11,$12,$16,$18,"5"}
    else if (($3+$4)  > 0 && ($5+$6)  > 0) {print $11,$12,$16,$18,"5""\n"$11,$13,$17,$18,"3"}
}' | sort -k 2 -n | sort -k 1 > capture.loci.ref


cat capture.loci.palmer | awk ' /cluster/ {print $1,$2,$2+1,$3,$4,$5,"Nanopore"}'    > capture.loci.palmer.process
cat capture.loci.palmer | awk '!/cluster/ {print $1,$2,$2+1,$3,$4,$5,$6,"Nanopore"}' > capture.loci.potential.process

cat capture.loci.ref capture.loci.ref.add > capture.loci.ref.all

if test "$mei" = "LINE"; then
    cat capture.loci.ref.all | grep L1HS |                                awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1hs
    cat capture.loci.ref.all | grep -v L1HS | grep L1PA |                 awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1pa
    cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep L1 |    awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1other
    cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep -v L1 | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1non
elif test \( "$mei" = "AluYa" \) -o \( "$mei" = "AluYb" \); then
    cat capture.loci.ref.all | grep AluYa |                                   awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYa5
    cat capture.loci.ref.all | grep AluYb |                                   awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYb8
    cat capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep AluY    | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYother
    cat capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep -v AluY | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYnon
else
    echo Unknown MEI type "$mei", expected LINE, AluYa, or AluYb.
    exit 1
fi

#######MEI
cp "${pp_mei}" Q.txt
cp capture.loci.palmer.process S.txt
# inter
intersect Q.txt S.txt > inter.txt
mv inter.txt inter.palmer.txt

cp capture.loci.potential.process S.txt
# inter
intersect Q.txt S.txt > inter.txt
mv inter.txt inter.palmer.add.txt

if test "$mei" = "LINE"; then
    cp ref.L1HS Q.txt
    cp capture.loci.r.l1hs S.txt
    # inter
    intersect Q.txt S.txt > inter.txt
    mv inter.txt inter.r.l1hs.txt

    cp ref.L1PA Q.txt
    cp capture.loci.r.l1pa S.txt
    # inter
    intersect Q.txt S.txt > inter.txt
    mv inter.txt inter.r.l1pa.txt

    cp ref.L1 Q.txt
    cp capture.loci.r.l1other S.txt
    # inter
    intersect Q.txt S.txt > inter.txt
    mv inter.txt inter.r.l1.txt
elif test \( "$mei" = "AluYa" \) -o \( "$mei" = "AluYb" \); then
    cp ref.AluYa5 Q.txt
    cp capture.loci.r.AluYa5 S.txt
    # inter
    intersect Q.txt S.txt > inter.txt
    mv inter.txt inter.r.AluYa5.txt

    cp ref.AluYb8 Q.txt
    cp capture.loci.r.AluYb8 S.txt
    # inter
    intersect Q.txt S.txt > inter.txt
    mv inter.txt inter.r.AluYb8.txt

    cp ref.AluY Q.txt
    cp capture.loci.r.AluYother S.txt
    # inter
    intersect Q.txt S.txt > inter.txt
    mv inter.txt inter.r.AluY.txt

    cp ref.Alu Q.txt
    cp capture.loci.r.AluYnon S.txt
    # inter
    intersect Q.txt S.txt > inter.txt
    mv inter.txt inter.r.Alu.txt
else
    echo Unknown MEI type "$mei", expected LINE, AluYa, or AluYb.
    exit 1
fi

# From Nanopal script: report the numbers of captured events
cat inter.palmer.txt inter.palmer.add.txt | grep Nano | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort | uniq -c > p.txt.fi

if test "$mei" = "LINE"; then
    cat inter.r.l1hs.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1hs.txt.fi
    cat inter.r.l1pa.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1pa.txt.fi
    cat inter.r.l1.txt   | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1.txt.fi
elif test \( "$mei" = "AluYa" \) -o \( "$mei" = "AluYb" \); then
    cat inter.r.AluYa5.txt | grep Nano | awk '{print $1,$2,$3}'    | sort | uniq -c > r.AluYa5.txt.fi
    cat inter.r.AluYb8.txt | grep Nano | awk '{print $1,$2,$3}'    | sort | uniq -c > r.AluYb8.txt.fi
    cat inter.r.AluY.txt   | grep Nano | awk '{print $1,$2,$3,$4}' | sort | uniq -c > r.AluY.txt.fi
    cat inter.r.Alu.txt    | grep Nano | awk '{print $1,$2,$3,$4}' | sort | uniq -c > r.Alu.txt.fi
else
    echo Unknown MEI type "$mei", expected LINE, AluYa, or AluYb.
    exit 1
fi

cp capture.loci.potential.process Q.txt

cp "${pp_mei}" S.txt
# inter
intersect Q.txt S.txt > inter.txt
awk '!/cluster/' inter.txt > inter.potential.txt

cp inter.potential.txt input_cluster.txt
cluster
cp clustered.txt potential.clustered.txt.fi

rm -f "$out_result_log"

echo There are ${SIGNAL1} reads caputuring putative L1 signals on one end. >> "$out_result_log"
echo There are ${SIGNAL2} reads caputuring putative L1 signals on both ends. >> "$out_result_log"
echo There are ${FAIL} reads having no putative L1 signals. >> "$out_result_log"

if test "$mei" = "LINE"; then
    cat p.txt.fi | awk '{sum+=$1} END {print "Non-reference L1Hs sum = ", sum}' >> "$out_result_log"
    cat r.l1hs.txt.fi | awk '{sum+=$1} END {print "Reference L1Hs sum = ", sum}' >> "$out_result_log"
    cat r.l1pa.txt.fi | awk '{sum+=$1} END {print "Reference L1PA sum = ", sum}' >> "$out_result_log"
    cat r.l1.txt.fi | awk '{sum+=$1} END {print "Other reference L1 sum = ", sum}' >> "$out_result_log"
    cat potential.clustered.txt.fi | awk '{sum+=$4} END {print "Potential non-reference L1Hs sum = ", sum}' >> "$out_result_log"
    echo "Number of non_LINE reads" >> "$out_result_log"
    cat capture.loci.r.l1non | wc -l | awk '{print $1}' >> "$out_result_log"

    echo "Number of non-reference L1Hs and the file (number of supporting reads + coordinate + overlap information)" >> "$out_result_log"
    wc -l p.txt.fi >> "$out_result_log"
    echo "Number of reference L1Hs and the file (number of supporting reads + coordinate)" >> "$out_result_log"
    wc -l r.l1hs.txt.fi >> "$out_result_log"
    echo "Number of reference L1PA and the file (number of supporting reads + coordinate)" >> "$out_result_log"
    wc -l r.l1pa.txt.fi >> "$out_result_log"
    echo "Number of other reference L1 and the file (number of supporting reads + coordinate)" >> "$out_result_log"
    wc -l r.l1.txt.fi >> "$out_result_log"
    echo "Number of potential specific non-reference L1Hs and the file (coordinate + number of supporting reads + number of left supporting reads + number of right supporting reads + strand)" >> "$out_result_log"
    wc -l potential.clustered.txt.fi >> "$out_result_log"
elif test \( "$mei" = "AluYa" \) -o \( "$mei" = "AluYb" \); then
    cat p.txt.fi | awk '{sum+=$1} END {print "Non-reference AluY reads sum = ", sum}' >> "$out_result_log"
    cat r.AluYa5.txt.fi | awk '{sum+=$1} END {print "Reference AluYa reads sum = ", sum}' >> "$out_result_log"
    cat r.AluYb8.txt.fi | awk '{sum+=$1} END {print "Reference AluYb reads sum = ", sum}' >> "$out_result_log"
    cat r.AluY.txt.fi | awk '{sum+=$1} END {print "Other reference AluY reads sum = ", sum}' >> "$out_result_log"
    cat capture.loci.r.AluYnon | grep Alu | wc -l | awk '{print "Other reference Alu reads sum = "$1}' >> "$out_result_log"
    cat potential.clustered.txt.fi | awk '{sum+=$4} END {print "Potential specific non-reference AluY reads sum = ", sum}' >> "$out_result_log"
    echo "Number of non_Alu reads" >> "$out_result_log"
    cat capture.loci.r.AluYnon | grep -v Alu | wc -l | awk '{print $1}' >> "$out_result_log"

    echo "Number of non-reference AluY and the file (number of supporting reads + coordinate + overlap information)" >> "$out_result_log"
    wc -l p.txt.fi >> "$out_result_log"
    echo "Number of reference AluYa and the file (number of supporting reads + coordinate)" >> "$out_result_log"
    wc -l r.AluYa5.txt.fi >> "$out_result_log"
    echo "Number of reference AluYb and the file (number of supporting reads + coordinate)" >> "$out_result_log"
    wc -l r.AluYb8.txt.fi >> "$out_result_log"
    echo "Number of other reference AluY and the file (number of supporting reads + coordinate)" >> "$out_result_log"
    wc -l r.AluY.txt.fi >> "$out_result_log"
    echo "Number of other reference Alu and the file (number of supporting reads + coordinate)" >> "$out_result_log"
    wc -l r.Alu.txt.fi >> "$out_result_log"
    echo "Number of potential specific non-reference AluY and the file (coordinate + number of supporting reads + number of left supporting reads + number of right supporting reads + strand)" >> "$out_result_log"
    wc -l potential.clustered.txt.fi >> "$out_result_log"
else
    echo Unknown MEI type "$mei", expected LINE, AluYa, or AluYb.
    exit 1
fi



