#!/usr/bin/env bash

set -eu
set -x

dataset_id="$1"
shift
ref_mei=$(readlink -f "$1") # hg38.RM.L1.ref
shift
pp_mei=$(readlink -f "$1") # L1.inter.fi
shift
in_summary="$1" # summary.final.txt
shift
revcomp_read_ids="$1" # RC.all.list
shift
mei="$1" # LINE
shift

out_dir="$1"
shift
out_summary="$1" # summary.final.2.txt
shift
out_result_log="$1" # formerly stdout
shift

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
# is in the list of reverse complemented reads, 0 if it's not.
join -1 1 -2 1 \
     -a 1 \
     -o auto -e "0" \
     --check-order \
     <(sort -k1 "$in_summary") \
     <(sort -k1 "$revcomp_read_ids" | sed -e 's/$/ 1/') \
     > "$out_summary"

# TODO Refactor the rest of this into separate chunks.

# Checking for signal on single end only
SIGNAL1=$(awk '(($3+$4)!=0 && ($5+$6)==0) || (($3+$4)==0 && ($5+$6)!=0)' "$out_summary" | wc -l)

# Checking for signal on both ends
SIGNAL2=$(awk '($3+$4)!=0 && ($5+$6)!=0' "$out_summary" | wc -l)

# Checking for no signal
FAIL=$(awk '($3+$4+$5+$6)==0' "$out_summary" | wc -l)

ENRICHMENT=$(echo | awk -v s1="$SIGNAL1" -v s2="$SIGNAL2" -v fail="$FAIL" '
    END {
        total = s1 + s2 + fail
        print ((s1 + s2) * 100.0) / total
    }
')

# The summary file now looks like this:                                                                (of the ALIGNMENT)
#                                                 nanopal --- palmer ----                              P&P ----  Ref ---------
#     read id ----------------------------        5+ 5- 3+ 3- 5+ 5- 3+ 3- position ------------------  5'   3'   5'     3'      RevComp
#     0000090b-e5db-457a-ac96-f111a076a4c8  1177  1  0  0  0  0  0  0  0  chr6   124810530  124811705  NON  NON  L1PA3  L1MA10  1
#     00000b31-b008-4cc5-b475-0e9f34753c63  9610  0  1  0  0  0  0  0  0  chr7   24100312   24109938   NON  NON  L1PA3  NON     1
#     0000146a-ced5-4e84-809b-862bb1aa087e  1885  0  0  0  0  0  0  0  0  NON    0          0          NA   NA   NA     NA      1
#     00001f2f-420b-4f1c-86e6-c9e6115968fa  484   0  0  0  0  0  0  0  0  chr11  2132414    2132857    NON  NON  NON    NON     0
#
#     col 1                                 2     3  4  5  6  7  8  9 10  11     12         13         14   15   16     17      18
#
# Note that the values of the P&P and Ref columns are both marking "an
# already-known event is here", but the formatting will look different.
#
# The Ref columns will have values of the type of event (e.g. L1HS), optionally
# clustered if there's more than one close by (e.g. L1HS/L1PA3).
#
# The P&P columns will only have a single value, and it'll be the long Palmer ID
# format like: 22.20.0.427778.0/1.12.12.0.0.31.0.0.cluster0_chr15_27930643_27930658_27930645_27930658_NA12878

# unmapped reads
cat "$out_summary" | awk '$11=="NON"' > summary.final.unmap.read.txt

# mapped read, nanopal got no signal, palmer did get signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)==0&&($7+$8+$9+$10)!=0' > summary.final.odd.read.txt

# mapped read, nanopal got no signal, palmer got no signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6+$7+$8+$9+$10)==0' > summary.final.no.L1.read.txt

# mapped read, nanopal got signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0' > summary.final.all.L1.read.txt

# mapped read, both palmer and nanopal have signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0&&($7+$8+$9+$10)!=0' > summary.final.PALMER.read.txt

# mapped read, only nanopal has signal
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0&&($7+$8+$9+$10)==0' > summary.final.ref.L1.read.txt

# Start with the set of reads where both Nanopal and Palmer have signal (from above), and extract out
# chr, start, already-palmer-known-event, revcomp, 5/3, ±, and id.
#
# TODO Do we need to check `revcomp` here as well, like in the next two `awk` calls?  If not, why not?
cat summary.final.PALMER.read.txt | awk '{
    #       palmer 5         palmer 3
    if      (($7+$8)  > 0 && ($9+$10) == 0 &&  $7 > 0) {print $11,$12,$14,$18,"5","+",$1} # chr, start, clusteringinfo, revcomp, 5, ±, id
    else if (($7+$8)  > 0 && ($9+$10) == 0 &&  $8 > 0) {print $11,$12,$14,$18,"5","-",$1}
    else if (($7+$8) == 0 && ($9+$10)  > 0 &&  $9 > 0) {print $11,$13,$15,$18,"3","+",$1} # chr, end,   clusteringinfo, revcomp, 3, ±, id
    else if (($7+$8) == 0 && ($9+$10)  > 0 && $10 > 0) {print $11,$13,$15,$18,"3","-",$1}

    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $7 > 0 &&  $9 > 0) {print $11,$12,$14,$18,"5","+",$1"\n"$11,$13,$15,$18,"3","+",$1}
    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $7 > 0 && $10 > 0) {print $11,$12,$14,$18,"5","+",$1"\n"$11,$13,$15,$18,"3","-",$1}
    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $8 > 0 &&  $9 > 0) {print $11,$12,$14,$18,"5","-",$1"\n"$11,$13,$15,$18,"3","+",$1}
    else if (($7+$8)  > 0 && ($9+$10)  > 0 && $8 > 0 && $10 > 0) {print $11,$12,$14,$18,"5","-",$1"\n"$11,$13,$15,$18,"3","-",$1}
}' | sort -k 2 -n | sort -k 1 > capture.loci.palmer

# Note: this appears to be trying to sort by two keys, but since it doesn't pass
# --stable sort will do a last-resort comparison, so this is actually only
# sorted by field 1 to end of line (but NOT numeric for field 2).  The proper
# way to do what I think this is trying to do is: sort -k 1,1 -k 2,2n

# Getting chr, start, clusteringinfo, revcomp, 5/3, for a mix of palmer/nanopal signals
cat summary.final.PALMER.read.txt | awk '{
    #       palmer 5         palmer 3         revcomp       nano 3/5
    if      (($7+$8)  > 0 && ($9+$10) == 0 && $18 == 0 && ($5+$6) > 0) {print $11,$13,$17,$18,"3",".",$1} # palmer 5  true, palmer 3 false, revcomp false, nano 3 true, id
    else if (($7+$8)  > 0 && ($9+$10) == 0 && $18 == 1 && ($3+$4) > 0) {print $11,$13,$17,$18,"3",".",$1} # palmer 5  true, palmer 3 false, revcomp  true, nano 5 true, id
    else if (($7+$8) == 0 && ($9+$10)  > 0 && $18 == 0 && ($3+$4) > 0) {print $11,$12,$16,$18,"5",".",$1} # palmer 5 false, palmer 3  true, revcomp false, nano 5 true, id
    else if (($7+$8) == 0 && ($9+$10)  > 0 && $18 == 1 && ($5+$6) > 0) {print $11,$12,$16,$18,"5",".",$1} # palmer 5 false, palmer 3  true, revcomp  true, nano 3 true, id
}' | sort -k 2 -n | sort -k 1 > capture.loci.ref.add

# Getting ... for when only nanopal has signal
cat summary.final.ref.L1.read.txt | awk '{
    #        nano 5           nano 3          revcomp
    if      (($3+$4)  > 0 && ($5+$6) == 0 && $18 == 0) {print $11,$12,$16,$18,"5",".",$1} # nano 5  true, nano 3 false, revcomp false, id
    else if (($3+$4)  > 0 && ($5+$6) == 0 && $18 == 1) {print $11,$13,$17,$18,"3",".",$1} # nano 5  true, nano 3 false, revcomp  true, id

    else if (($3+$4) == 0 && ($5+$6)  > 0 && $18 == 0) {print $11,$13,$17,$18,"3",".",$1} # nano 5 false, nano 3  true, revcomp false, id
    else if (($3+$4) == 0 && ($5+$6)  > 0 && $18 == 1) {print $11,$12,$16,$18,"5",".",$1} # nano 5 false, nano 3  true, revcomp  true, id

    else if (($3+$4)  > 0 && ($5+$6)  > 0            ) {print $11,$12,$16,$18,"5",".",$1"\n"$11,$13,$17,$18,"3",".",$1} # nano 5 true, nano 3 true, revcomp whatever, id
}' | sort -k 2 -n | sort -k 1 > capture.loci.ref


cat capture.loci.palmer | awk ' /cluster/ {print $1,$2,$2+1,$3,$4,$5,$6,$7,"Nanopore"}' > capture.loci.palmer.process
cat capture.loci.palmer | awk '!/cluster/ {print $1,$2,$2+1,$3,$4,$5,$6,$7}' > capture.loci.potential.process

cat capture.loci.ref capture.loci.ref.add > capture.loci.ref.all

case "$mei" in
    LINE)
        cat capture.loci.ref.all | grep L1HS |                                awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.l1hs
        cat capture.loci.ref.all | grep -v L1HS | grep L1PA |                 awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.l1pa
        cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep L1 |    awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.l1other
        cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep -v L1 | awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.l1non
        ;;
    AluYa | AluYb)
        cat capture.loci.ref.all | grep AluYa |                                   awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.AluYa5
        cat capture.loci.ref.all | grep AluYb |                                   awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.AluYb8
        cat capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep AluY    | awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.AluYother
        cat capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep -v AluY | awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.AluYnon
        ;;
    SVA_E | SVA_F)
        cat capture.loci.ref.all | grep SVA_E |                                  awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.SVA_E
        cat capture.loci.ref.all | grep SVA_F |                                  awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.SVA_F
        cat capture.loci.ref.all | grep -v SVA_E | grep -v SVA_F | grep SVA    | awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.SVAother
        cat capture.loci.ref.all | grep -v SVA_E | grep -v SVA_F | grep -v SVA | awk '{print $1,$2,$2+1,$3,$4,$5,$7,"Nanopore"}' > capture.loci.r.SVAnon
        ;;
    *)
        echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
        exit 1
        ;;
esac

intersect "${pp_mei}" capture.loci.palmer.process    > inter.palmer.txt
intersect "${pp_mei}" capture.loci.potential.process > inter.palmer.add.txt

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
    echo "There are ${SIGNAL1} reads capturing putative signals on one end."
    echo "There are ${SIGNAL2} reads capturing putative signals on both ends."
    echo "There are ${FAIL} reads having no putative signals."
    echo "Enrichment: ${ENRICHMENT}% of reads have signal(s)."

    case "$mei" in
        LINE)
            cat p.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Non-reference L1Hs sum = ", sum}'
            cat r.l1hs.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Reference L1Hs sum = ", sum}'
            cat r.l1pa.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Reference L1PA sum = ", sum}'
            cat r.l1.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Other reference L1 sum = ", sum}'
            cat potential.clustered.txt.fi | awk 'BEGIN {sum=0} {sum+=$4} END {print "Potential non-reference L1Hs sum = ", sum}'
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
            echo "Number of potential Nanopore-specific non-reference L1Hs and the file (coordinate + number of supporting reads + number of left supporting reads + number of right supporting reads + strand)"
            wc -l potential.clustered.txt.fi
            ;;
        AluYa | AluYb)
            cat p.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Non-reference AluY reads sum = ", sum}'
            cat r.AluYa5.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Reference AluYa reads sum = ", sum}'
            cat r.AluYb8.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Reference AluYb reads sum = ", sum}'
            cat r.AluY.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Other reference AluY reads sum = ", sum}'
            cat capture.loci.r.AluYnon | grep -c Alu | awk '{print "Other reference Alu reads sum = "$1}'
            cat potential.clustered.txt.fi | awk 'BEGIN {sum=0} {sum+=$4} END {print "Potential specific non-reference AluY reads sum = ", sum}'
            echo "Number of non_Alu reads"
            cat capture.loci.r.AluYnon | grep -cv Alu | awk '{print $1}'

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
            cat p.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Non-reference SVA reads sum = ", sum}'
            cat r.SVA_E.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Reference SVA_E reads sum = ", sum}'
            cat r.SVA_F.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Reference SVA_F reads sum = ", sum}'
            cat r.SVA_D.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Reference SVA_D reads sum = ", sum}'
            cat r.SVAother.txt.fi | awk 'BEGIN {sum=0} {sum+=$1} END {print "Other reference SVA reads sum = ", sum}'
            cat potential.clustered.txt.fi | awk 'BEGIN {sum=0} {sum+=$4} END {print "Potential specific non-reference SVA reads sum = ", sum}'
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
