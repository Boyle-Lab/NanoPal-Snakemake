#!/usr/bin/env bash

set -eu

bam="$1" # Nanopore.sorted.bam
ref_mei=$(readlink -f "$2") # hg38.RM.L1.ref
pp_mei=$(readlink -f "$3") # L1.inter.fi
in_summary="$4" # summary.final.txt

out_dir="$5"
out_summary="$6" # summary.final.2.txt
out_result_log="$7" # formerly stdout

mkdir -p "$out_dir"
cd "$out_dir"

##different types of LINE
#######MEI
grep    L1PA "$ref_mei" > ref.L1PA
grep    L1HS "$ref_mei" > ref.L1HS
grep -v L1HS "$ref_mei"| grep -v L1PA | grep L1 > ref.L1

##Valid mapping reads
samtools view "$bam" -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 -f 0x10 \
    | awk '{print $1}' \
    >  RC.all.list

rm -f RC.tag
while read -r aa b c d e f g h i j  k l m n o p q r
do
awk -v aa="$aa" '
    $1 == aa { print 1; exit 0 }
    END { print 0 }
' RC.all.list >> RC.tag

wait
done < "$in_summary"

paste "$in_summary" RC.tag > "$out_summary"

SIGNAL1=$(awk '(($3+$4)!=0 && ($5+$6)==0) || (($3+$4)==0 && ($5+$6)!=0)' "$out_summary" | wc -l)
SIGNAL2=$(awk '($3+$4)!=0 && ($5+$6)!=0' "$out_summary" | wc -l)
FAIL=$(awk '($3+$4+$5+$6)==0' "$out_summary" | wc -l)

cat "$out_summary" | awk '$11=="NON"' > summary.final.unmap.read.txt

cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)==0&&($7+$8+$9+$10)!=0' > summary.final.odd.read.txt

cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6+$7+$8+$9+$10)==0' > summary.final.no.L1.read.txt

cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0' > summary.final.all.L1.read.txt

cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0&&($7+$8+$9+$10)!=0' > summary.final.PALMER.L1.read.txt
cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0&&($7+$8+$9+$10)!=0' > summary.final.PALMER.read.txt

cat "$out_summary" | awk '$11!="NON"' | awk '($3+$4+$5+$6)!=0&&($7+$8+$9+$10)==0' > summary.final.ref.L1.read.txt

cat summary.final.PALMER.read.txt | awk '{
    if(($7+$8)>0&&($9+$10)==0&&$7>0) {print $11,$12,$14,$18,"5","+"}
    else if (($7+$8) > 0  && ($9+$10) == 0 &&$8>0) {print $11,$12,$14,$18,"5","-"}
    else if (($7+$8) == 0 && ($9+$10) > 0 && $9>0) {print $11,$13,$15,$18,"3","+"}
    else if (($7+$8) == 0 && ($9+$10) > 0 && $10>0) {print $11,$13,$15,$18,"3","-"}
    else if (($7+$8) > 0  && ($9+$10) > 0 && $7>0&&$9>0) {print $11,$12,$14,$18,"5","+""\n"$11,$13,$15,$18,"3","+"}
    else if (($7+$8) > 0  && ($9+$10) > 0 && $7>0&&$10>0) {print $11,$12,$14,$18,"5","+""\n"$11,$13,$15,$18,"3","-"}
    else if (($7+$8) > 0  && ($9+$10) > 0 && $8>0&&$9>0) {print $11,$12,$14,$18,"5","-""\n"$11,$13,$15,$18,"3","+"}
    else if (($7+$8) > 0  && ($9+$10) > 0 && $8>0&&$10>0) {print $11,$12,$14,$18,"5","-""\n"$11,$13,$15,$18,"3","-"}
}' | sort -k 2 -n | sort -k 1 > capture.loci.palmer

cat summary.final.PALMER.L1.read.txt | awk '{
    if(($7+$8)>0&&($9+$10)==0&&$18==0&&($5+$6)>0) {print $11,$13,$17,$18,"3"}
    else if (($7+$8)>0&&($9+$10)==0&&$18==1&&($3+$4)>0) {print $11,$13,$17,$18,"3"}
    else if (($7+$8)==0&&($9+$10)>0&&$18==0&&($3+$4)>0) {print $11,$12,$16,$18,"5"}
    else if (($7+$8)==0&&($9+$10)>0&&$18==1&&($5+$6)>0) {print $11,$12,$16,$18,"5"}
}' | sort -k 2 -n | sort -k 1 > capture.loci.ref.add

cat summary.final.ref.L1.read.txt | awk '{
    if(($3+$4)>0&&($5+$6)==0&&$18==0) {print $11,$12,$16,$18,"5"}
    else if(($3+$4)>0&&($5+$6)==0&&$18==1) {print $11,$13,$17,$18,"3"}
    else if (($3+$4)==0&&($5+$6)>0&&$18==0) {print $11,$13,$17,$18,"3"}
    else if (($3+$4)==0&&($5+$6)>0&&$18==1) {print $11,$12,$16,$18,"5"}
    else if(($3+$4)>0&&($5+$6)>0) {print $11,$12,$16,$18,"5""\n"$11,$13,$17,$18,"3"}
}' | sort -k 2 -n | sort -k 1 > capture.loci.ref


cat capture.loci.palmer | grep cluster    | awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}'    > capture.loci.palmer.process
cat capture.loci.palmer | grep -v cluster | awk '{print $1,$2,$2+1,$3,$4,$5,$6,"Nanopore"}' > capture.loci.potential.process


cat capture.loci.ref capture.loci.ref.add > capture.loci.ref.all
cat capture.loci.ref.all | grep L1HS |                                 awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1hs
cat capture.loci.ref.all | grep -v L1HS | grep L1PA |                  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1pa
cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep L1 |     awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1other
cat capture.loci.ref.all | grep -v L1HS | grep -v L1PA | grep -v L1 |  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.l1non

# TODO these inter calls were originally inter.0118.o, not sure if that's the
# same as the inter binary I have or something slightly different?

#######MEI
cp "${pp_mei}" Q.txt
cp capture.loci.palmer.process S.txt
inter
mv inter.txt inter.palmer.txt

cp capture.loci.potential.process S.txt
inter
mv inter.txt inter.palmer.add.txt

#######MEI
cp ref.L1HS Q.txt
cp capture.loci.r.l1hs S.txt
inter
mv inter.txt inter.r.l1hs.txt

#######MEI
cp ref.L1PA Q.txt
cp capture.loci.r.l1pa S.txt
inter
mv inter.txt inter.r.l1pa.txt

#######MEI
cp ref.L1 Q.txt
cp capture.loci.r.l1other S.txt
inter
mv inter.txt inter.r.l1.txt

# From Nanopal script: report the numbers of captured events
cat inter.palmer.txt inter.palmer.add.txt | grep Nano | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort | uniq -c > p.txt.fi

cat inter.r.l1hs.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1hs.txt.fi

cat inter.r.l1pa.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1pa.txt.fi

cat inter.r.l1.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.l1.txt.fi

cp capture.loci.potential.process Q.txt

#######MEI
cp "${pp_mei}" S.txt
inter
cat inter.txt | grep -v cluster > inter.potential.txt

#######MEI
cp inter.potential.txt input_cluster.txt
cluster
cp clustered.txt potential.clustered.txt.fi

rm -f "$out_result_log"

echo There are ${SIGNAL1} reads caputuring putative L1 signals on one end. >> "$out_result_log"
echo There are ${SIGNAL2} reads caputuring putative L1 signals on both ends. >> "$out_result_log"
echo There are ${FAIL} reads having no putative L1 signals. >> "$out_result_log"

cat p.txt.fi | awk '{sum+=$1} END {print "Non-reference L1Hs sum = ", sum}' >> "$out_result_log"
cat r.l1hs.txt.fi | awk '{sum+=$1} END {print "Reference L1Hs sum = ", sum}' >> "$out_result_log"
cat r.l1pa.txt.fi | awk '{sum+=$1} END {print "Reference L1PA sum = ", sum}' >> "$out_result_log"
cat r.l1.txt.fi | awk '{sum+=$1} END {print "Other reference L1 sum = ", sum}' >> "$out_result_log"
cat potential.clustered.txt.fi | awk '{sum+=$4} END {print "Potential non-reference L1Hs sum = ", sum}' >> "$out_result_log"
echo "Number of non_SVA reads" >> "$out_result_log"
cat capture.loci.r.l1non | wc -l | awk '{print $1}' >> "$out_result_log"


echo "Number of non-reference SVA and the file (number of suppoting reads + coordinate + overlap information)" >> "$out_result_log"
wc -l p.txt.fi >> "$out_result_log"
echo "Number of reference L1Hs and the file (number of suppoting reads + coordinate)" >> "$out_result_log"
wc -l r.l1hs.txt.fi >> "$out_result_log"
echo "Number of reference L1PA and the file (number of suppoting reads + coordinate)" >> "$out_result_log"
wc -l r.l1pa.txt.fi >> "$out_result_log"
echo "Number of other reference L1 and the file (number of suppoting reads + coordinate)" >> "$out_result_log"
wc -l r.l1.txt.fi >> "$out_result_log"
echo "Number of potential specific non-reference L1Hs and the file (coordinate + number of suppoting reads + number of left suppoting reads + number of right suppoting reads + strand)" >> "$out_result_log"
wc -l potential.clustered.txt.fi >> "$out_result_log"
