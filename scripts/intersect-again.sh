#!/usr/bin/env bash

set -eu
set -x

ref_mei=$(readlink -f "$1") # hg38.RM.L1.ref
pp_mei=$(readlink -f "$2") # L1.inter.fi
in_summary="$3" # summary.final.txt
valid_read_ids="$4" # RC.all.list
mei="$5" # LINE

out_dir="$6"
out_summary="$7" # summary.final.2.txt
out_result_log="$8" # formerly stdout

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

case "$mei" in
    LINE)
        cat "$ref_mei" | grep    L1PA > ref.L1PA
        cat "$ref_mei" | grep    L1HS > ref.L1HS
        cat "$ref_mei" | grep -v L1HS| grep -v L1PA | grep L1 > ref.L1
        ;;
    AluYa | AluYb)
        cat "$ref_mei" | grep AluYa > ref.AluYa5
        cat "$ref_mei" | grep AluYb > ref.AluYb8
        cat "$ref_mei" | grep -v AluYb | grep -v AluYa | grep AluY > ref.AluY
        cat "$ref_mei" | grep -v AluYb | grep -v AluYa | grep -v AluY | grep Alu > ref.Alu
        ;;
    SVA_E | SVA_F)
        cat "$ref_mei" | grep SVA_E > ref.SVA_E
        cat "$ref_mei" | grep SVA_F > ref.SVA_F
        cat "$ref_mei" | grep -v SVA_E | grep -v SVA_F | grep SVA > ref.SVA
        ;;
    *)
        echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
        exit 1
        ;;
esac

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

case "$mei" in
    LINE)
        cat capture.loci.ref.all | grep L1HS |                                awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1hs
        cat capture.loci.ref.all | grep -v L1HS | grep L1PA |                 awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1pa
        cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep L1 |    awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1other
        cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep -v L1 | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1non
        ;;
    AluYa | AluYb)
        cat capture.loci.ref.all | grep AluYa |                                   awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYa5
        cat capture.loci.ref.all | grep AluYb |                                   awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYb8
        cat capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep AluY    | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYother
        cat capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep -v AluY | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYnon
        ;;
    SVA_E | SVA_F)
        cat capture.loci.ref.all | grep SVA_E |                                  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.SVA_E
        cat capture.loci.ref.all | grep SVA_F |                                  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.SVA_F
        cat capture.loci.ref.all | grep -v SVA_E | grep -v SVA_F | grep SVA    | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.SVAother
        cat capture.loci.ref.all | grep -v SVA_E | grep -v SVA_F | grep -v SVA | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.SVAnon
        ;;
    *)
        echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
        exit 1
        ;;
esac

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

case "$mei" in
    LINE)
        intersect ref.L1HS capture.loci.r.l1hs    > inter.r.l1hs.txt
        intersect ref.L1PA capture.loci.r.l1pa    > inter.r.l1pa.txt
        intersect ref.L1   capture.loci.r.l1other > inter.r.l1.txt
        ;;
    AluYa | AluYb)
        intersect ref.AluYa5 capture.loci.r.AluYa5    > inter.r.AluYa5.txt
        intersect ref.AluYb8 capture.loci.r.AluYb8    > inter.r.AluYb8.txt
        intersect ref.AluY   capture.loci.r.AluYother > inter.r.AluY.txt
        intersect ref.Alu    capture.loci.r.AluYnon   > inter.r.Alu.txt
        ;;
    SVA_E | SVA_F)
        intersect ref.SVA_E capture.loci.r.SVA_E    > inter.r.SVA_E.txt
        intersect ref.SVA_F capture.loci.r.SVA_F    > inter.r.SVA_F.txt
        intersect ref.SVA   capture.loci.r.SVAother > inter.r.SVA.txt
        ;;
    *)
        echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
        exit 1
        ;;
esac

# From Nanopal script: report the numbers of captured events
cat inter.palmer.txt inter.palmer.add.txt | grep Nano | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort | uniq -c > p.txt.fi

case "$mei" in
    LINE)
        cat inter.r.l1hs.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1hs.txt.fi
        cat inter.r.l1pa.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1pa.txt.fi
        cat inter.r.l1.txt   | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1.txt.fi
        ;;
    AluYa | AluYb)
        cat inter.r.AluYa5.txt | grep Nano | awk '{print $1,$2,$3}'    | sort | uniq -c > r.AluYa5.txt.fi
        cat inter.r.AluYb8.txt | grep Nano | awk '{print $1,$2,$3}'    | sort | uniq -c > r.AluYb8.txt.fi
        cat inter.r.AluY.txt   | grep Nano | awk '{print $1,$2,$3,$4}' | sort | uniq -c > r.AluY.txt.fi
        cat inter.r.Alu.txt    | grep Nano | awk '{print $1,$2,$3,$4}' | sort | uniq -c > r.Alu.txt.fi
        ;;
    SVA_E | SVA_F)
        cat inter.r.SVA_E.txt | grep Nano | awk '{print $1,$2,$3}'     | sort | uniq -c > r.SVA_E.txt.fi
        cat inter.r.SVA_F.txt | grep Nano | awk '{print $1,$2,$3}'     | sort | uniq -c > r.SVA_F.txt.fi
        cat inter.r.SVA.txt   | grep Nano | awk '{print $1,$2,$3, $4}' | sort | uniq -c > r.SVA.txt.fi

        cat r.SVA.txt.fi | grep SVA_D    > r.SVA_D.txt.fi
        cat r.SVA.txt.fi | grep -v SVA_D > r.SVAother.txt.fi
        ;;
    *)
        echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
        exit 1
        ;;
esac

# awk ! instead of grep -v because grep will exit nonzero if nothing is found.
intersect capture.loci.potential.process "${pp_mei}" \
    | awk '!/cluster/' \
    > inter.potential.txt

cp inter.potential.txt input_cluster.txt
cluster
cp clustered.txt potential.clustered.txt.fi

{
    echo There are "${SIGNAL1}" reads caputuring putative L1 signals on one end.
    echo There are "${SIGNAL2}" reads caputuring putative L1 signals on both ends.
    echo There are "${FAIL}" reads having no putative L1 signals.

    case "$mei" in
        LINE)
            cat p.txt.fi | awk '{sum+=$1} END {print "Non-reference L1Hs sum = ", sum}'
            cat r.l1hs.txt.fi | awk '{sum+=$1} END {print "Reference L1Hs sum = ", sum}'
            cat r.l1pa.txt.fi | awk '{sum+=$1} END {print "Reference L1PA sum = ", sum}'
            cat r.l1.txt.fi | awk '{sum+=$1} END {print "Other reference L1 sum = ", sum}'
            cat potential.clustered.txt.fi | awk '{sum+=$4} END {print "Potential non-reference L1Hs sum = ", sum}'
            echo "Number of non_LINE reads"
            cat capture.loci.r.l1non | wc -l | awk '{print $1}'

            echo "Number of non-reference L1Hs and the file (number of supporting reads + coordinate + overlap information)"
            wc -l p.txt.fi
            echo "Number of reference L1Hs and the file (number of supporting reads + coordinate)"
            wc -l r.l1hs.txt.fi
            echo "Number of reference L1PA and the file (number of supporting reads + coordinate)"
            wc -l r.l1pa.txt.fi
            echo "Number of other reference L1 and the file (number of supporting reads + coordinate)"
            wc -l r.l1.txt.fi
            echo "Number of potential specific non-reference L1Hs and the file (coordinate + number of supporting reads + number of left supporting reads + number of right supporting reads + strand)"
            wc -l potential.clustered.txt.fi
            ;;
        AluYa | AluYb)
            cat p.txt.fi | awk '{sum+=$1} END {print "Non-reference AluY reads sum = ", sum}'
            cat r.AluYa5.txt.fi | awk '{sum+=$1} END {print "Reference AluYa reads sum = ", sum}'
            cat r.AluYb8.txt.fi | awk '{sum+=$1} END {print "Reference AluYb reads sum = ", sum}'
            cat r.AluY.txt.fi | awk '{sum+=$1} END {print "Other reference AluY reads sum = ", sum}'
            cat capture.loci.r.AluYnon | grep -c Alu | awk '{print "Other reference Alu reads sum = "$1}'
            cat potential.clustered.txt.fi | awk '{sum+=$4} END {print "Potential specific non-reference AluY reads sum = ", sum}'
            echo "Number of non_Alu reads"
            cat capture.loci.r.AluYnon | grep -v Alu | awk '{print $1}'

            echo "Number of non-reference AluY and the file (number of supporting reads + coordinate + overlap information)"
            wc -l p.txt.fi
            echo "Number of reference AluYa and the file (number of supporting reads + coordinate)"
            wc -l r.AluYa5.txt.fi
            echo "Number of reference AluYb and the file (number of supporting reads + coordinate)"
            wc -l r.AluYb8.txt.fi
            echo "Number of other reference AluY and the file (number of supporting reads + coordinate)"
            wc -l r.AluY.txt.fi
            echo "Number of other reference Alu and the file (number of supporting reads + coordinate)"
            wc -l r.Alu.txt.fi
            echo "Number of potential specific non-reference AluY and the file (coordinate + number of supporting reads + number of left supporting reads + number of right supporting reads + strand)"
            wc -l potential.clustered.txt.fi
            ;;
        SVA_E | SVA_F)
            cat p.txt.fi | awk '{sum+=$1} END {print "Non-reference SVA reads sum = ", sum}'
            cat r.SVA_E.txt.fi | awk '{sum+=$1} END {print "Reference SVA_E reads sum = ", sum}'
            cat r.SVA_F.txt.fi | awk '{sum+=$1} END {print "Reference SVA_F reads sum = ", sum}'
            cat r.SVA_D.txt.fi | awk '{sum+=$1} END {print "Reference SVA_D reads sum = ", sum}'
            cat r.SVAother.txt.fi | awk '{sum+=$1} END {print "Other reference SVA reads sum = ", sum}'
            cat potential.clustered.txt.fi | awk '{sum+=$4} END {print "Potential specific non-reference SVA reads sum = ", sum}'
            echo "Number of non_SVA reads"
            cat capture.loci.r.SVAnon | wc -l | awk '{print $1}'

            echo "Number of non-reference SVA and the file (number of supporting reads + coordinate + overlap information)"
            wc -l p.txt.fi
            echo "Number of reference SVA_E and the file (number of supporting reads + coordinate)"
            wc -l r.SVA_E.txt.fi
            echo "Number of reference SVA_F and the file (number of supporting reads + coordinate)"
            wc -l r.SVA_F.txt.fi
            echo "Number of reference SVA_D and the file (number of supporting reads + coordinate)"
            wc -l r.SVA_D.txt.fi
            echo "Number of other reference SVA and the file (number of supporting reads + coordinate)"
            wc -l r.SVAother.txt.fi
            echo "Number of potential specific non-reference SVA and the file (coordinate + number of supporting reads + number of left supporting reads + number of right supporting reads + strand)"
            wc -l potential.clustered.txt.fi
            ;;
        *)
            echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
            exit 1
            ;;
    esac
} > "$out_result_log"

