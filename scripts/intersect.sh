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

# We start with the results from Palmer and Nanopal's own BLAST queries:
#
# palmer_reads:
#                            1           2            3        4         5        6         7        8         9        10        11     12         13
#                            read_name   read_length  n_plus5  n_minus5  n_plus3  n_minus3  p_plus5  p_minus5  p_plus3  p_minus3  chr    start_pos  end_pos
# f5edeb8c-f225-4b47-b7ee-1c0df4dea038   172594       0        1         0        0         0        0         0        0         chr18  74826199   74999245

# Extract the 5' position of any (non-reference?  unmapped? what is that grep?) reads
# chr start start+1 read_id
awk '{print $11, $12, $12+1, $1}' "$palmer_reads"  | grep -v "NON" > read.5.Q.txt

# Extract the 3' position of any NON-reference reads
# chr end   end+1   read_id
awk '{print $11, $13, $13+1, $1}' "$palmer_reads"  | grep -v "NON" > read.3.Q.txt

cp read.5.Q.txt Q.txt

# ref_mei:
# chr     start           end             type    ?       strand  ?       ?
# chr1    100000000       100000637       L1M2    ref     -       6514    7153
cp "${ref_mei}"  S.txt

intersect Q.txt S.txt > inter.txt
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt

collapse input_rm_cluster.txt > output_rm_cluster.txt
mv output_rm_cluster.txt read.5.Q.ref.txt

#######MEI
cp "${orig_mei}"  S.txt
# inter
intersect Q.txt S.txt > inter.txt

case "$mei" in
    LINE | SVA_E | SVA_F)
        awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' \
            inter.txt > read.5.Q.P.txt
        ;;
    AluYa | AluYb)
        cp inter.txt inter-pre-dedupe-1.txt

        # This is deduping on field FIVE, beacause -f 4 means SKIP the first four
        # fields.  So it's deduping based on read ID.
        cat inter.txt \
            | awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' \
            | awk '{print $1,$2,$3,$5,$4}' \
            | sort -n | uniq -f 4 \
            | awk '{print $1,$2,$3,$5,$4}' \
            > read.5.Q.P.txt
        ;;
    *)
        echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
        exit 1
        ;;
esac

cp read.3.Q.txt Q.txt

#######MEI
cp "${ref_mei}"  S.txt
# inter
intersect Q.txt S.txt > inter.txt
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt
# RM_collapse
collapse input_rm_cluster.txt > output_rm_cluster.txt
mv output_rm_cluster.txt read.3.Q.ref.txt

#######MEI
cp "${orig_mei}"  S.txt
# inter
intersect Q.txt S.txt > inter.txt

case "$mei" in
    LINE | SVA_E | SVA_F)
        awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' \
            inter.txt > read.3.Q.P.txt
        ;;
    AluYa | AluYb)
        cp inter.txt inter-pre-dedupe-2.txt

        cat inter.txt \
            | awk '{ if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8 }' \
            | awk '{print $1,$2,$3,$5,$4}' \
            | sort -n | uniq -f 4 \
            | awk '{print $1,$2,$3,$5,$4}' \
            > read.3.Q.P.txt
        ;;
    *)
        echo Unknown MEI type "$mei", expected one of: LINE AluYa AluYb SVA_E SVA_F
        exit 1
        ;;
esac

sort -k 4 read.5.Q.P.txt > 5.P.inter.txt
sort -k 4 read.3.Q.P.txt > 3.P.inter.txt
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

