#!/usr/bin/env bash

set -eu

palmer_reads="$1" # read.all.palmer.final.txt
ref_mei=$(readlink -f "$2") # hg38.RM.L1.ref
orig_mei=$(readlink -f "$3") # PALMER.NA12878.L1.txt

out_dir="$4"
out_summary="$5" # summary.final.txt

mkdir -p "$out_dir"
cd "$out_dir"

awk '{print $11, $12, $12+1, $1}' "$palmer_reads"  | grep -v "NON" > read.5.Q.txt
awk '{print $11, $13, $13+1, $1}' "$palmer_reads"  | grep -v "NON" > read.3.Q.txt

cp read.5.Q.txt Q.txt

#######MEI
cp ${ref_mei}  S.txt
inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt
RM_collapse
mv output_rm_cluster.txt read.5.Q.ref.txt

#######MEI
cp ${orig_mei}  S.txt
inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > read.5.Q.P.txt
cp read.3.Q.txt Q.txt

#######MEI
cp ${ref_mei}  S.txt
inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt
RM_collapse
mv output_rm_cluster.txt read.3.Q.ref.txt

#######MEI
cp ${orig_mei}  S.txt
inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > read.3.Q.P.txt

sort -k 4 read.5.Q.P.txt > 5.P.inter.txt
sort -k 4 read.3.Q.P.txt > 3.P.inter.txt
sort -k 4 read.5.Q.ref.txt > 5.ref.inter.txt
sort -k 4 read.3.Q.ref.txt > 3.ref.inter.txt

paste 5.P.inter.txt 3.P.inter.txt 5.ref.inter.txt 3.ref.inter.txt | awk '{print $4,$5,$10,$15,$20}' > read.RM.txt

rm -f read.RM.2.txt
while read -r a b c d e f g h i j k l m; do
    awk '
        BEGIN { score = 0 }
        {
            if($1=="'$a'") {
                score += 1
                print $2,$3,$4,$5
            }
        }
        END {
            if (score == 0) {
                print "NA","NA","NA","NA"
            }
        }' read.RM.txt >> read.RM.2.txt
done < "$palmer_reads"

# Summary table of MEIs for each read
paste "$palmer_reads" read.RM.2.txt > summary.final.txt
