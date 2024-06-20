#!/usr/bin/env bash

set -eu
set -x

palmer_reads="$1" # read.all.palmer.final.txt
ref_mei=$(readlink -f "$2") # hg38.RM.L1.ref
orig_mei=$(readlink -f "$3") # PALMER.NA12878.L1.txt
mei="$4"

out_dir="$5"
out_summary="$6" # summary.final.txt

mkdir -p "$out_dir"
cd "$out_dir"

awk '{print $11, $12, $12+1, $1}' "$palmer_reads"  | grep -v "NON" > read.5.Q.txt
awk '{print $11, $13, $13+1, $1}' "$palmer_reads"  | grep -v "NON" > read.3.Q.txt

cp read.5.Q.txt Q.txt

#######MEI
cp ${ref_mei}  S.txt
# inter
intersect Q.txt S.txt > inter.txt
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt
# RM_collapse
collapse input_rm_cluster.txt > output_rm_cluster.txt
mv output_rm_cluster.txt read.5.Q.ref.txt

#######MEI
cp ${orig_mei}  S.txt
# inter
intersect Q.txt S.txt > inter.txt

if test "$mei" = "LINE"; then
    awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' inter.txt > read.5.Q.P.txt
elif test \( "$mei" = "AluYa" \) -o \( "$mei" = "AluYb" \); then
    cp inter.txt inter-pre-dedupe-1.txt

    # This is deduping on field FIVE, beacause -f 4 means SKIP the first four
    # fields.  So it's deduping based on read ID.
    cat inter.txt \
        | awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' \
        | awk '{print $1,$2,$3,$5,$4}' \
        | sort -n | uniq -f 4 \
        | awk '{print $1,$2,$3,$5,$4}' \
        > read.5.Q.P.txt
else
    echo Unknown MEI type "$mei", expected LINE, AluYa, or AluYb.
    exit 1
fi

cp read.3.Q.txt Q.txt

#######MEI
cp ${ref_mei}  S.txt
# inter
intersect Q.txt S.txt > inter.txt
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt
# RM_collapse
collapse input_rm_cluster.txt > output_rm_cluster.txt
mv output_rm_cluster.txt read.3.Q.ref.txt

#######MEI
cp ${orig_mei}  S.txt
# inter
intersect Q.txt S.txt > inter.txt

if test "$mei" = "LINE"; then
    awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > read.3.Q.P.txt
elif test \( "$mei" = "AluYa" \) -o \( "$mei" = "AluYb" \); then
    cp inter.txt inter-pre-dedupe-2.txt

    cat inter.txt \
        | awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' \
        | awk '{print $1,$2,$3,$5,$4}' \
        | sort -n | uniq -f 4 \
        | awk '{print $1,$2,$3,$5,$4}' \
        > read.3.Q.P.txt
else
    echo Unknown MEI type "$mei", expected LINE, AluYa, or AluYb.
    exit 1
fi

sort -k 4 read.5.Q.P.txt > 5.P.inter.txt
sort -k 4 read.3.Q.P.txt > 3.P.inter.txt
sort -k 4 read.5.Q.ref.txt > 5.ref.inter.txt
sort -k 4 read.3.Q.ref.txt > 3.ref.inter.txt

paste 5.P.inter.txt 3.P.inter.txt 5.ref.inter.txt 3.ref.inter.txt | awk '{print $4,$5,$10,$15,$20}' > read.RM.txt

# rm -f read.RM.2.txt
# while read -r a b c d e f g h i j k l m; do
#     awk '
#         BEGIN { score = 0 }
#         {
#             if($1=="'$a'") {
#                 score += 1
#                 print $2,$3,$4,$5
#             }
#         }
#         END {
#             if (score == 0) {
#                 print "NA","NA","NA","NA"
#             }
#         }' read.RM.txt >> read.RM.2.txt
# done < "$palmer_reads"

join -1 1 -2 1 \
     -a 1 \
     -o auto -e NA \
     --check-order \
     <(cat "$palmer_reads" | sort) \
     <(cat read.RM.txt | sort) \
     > "$out_summary"

     # > read.RM.2.txt

# Summary table of MEIs for each read
# paste "$palmer_reads" read.RM.2.txt > "$out_summary"
