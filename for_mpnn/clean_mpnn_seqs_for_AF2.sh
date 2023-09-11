#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N clean_seqs
#$ -l h_rt=1:00:00
#$ -o log_files/


# This cleans the sequences from pMPNN to be better read for AF2 

PREFIX=="3C62_chainD_rog_40_5-8_60"
FOLD="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/$PREFIX/seqs"

# Remove first two lines (input sequence ID and seq)
sed -i -e 1,2d $FOLD/*

NEW_ID=">"
# Replace ID of sequence with filename + sequence number
for FILE in $FOLD/*; do
    if [ -f "$FILE" ]; then        
        NEW_ID+=$(basename "$FILE")
        NEW_ID=${NEW_ID:0:-3}
        sed -i "/>/s/.*/$NEW_ID/" $FILE
        awk  -i inplace 'NR % 2 == 1 {printf "%s_%.1f\n", $0, ((NR - 1) / 2)} NR % 2 == 0' $FILE
        awk  -i inplace 'NR % 2 == 1 {print substr($0, 1, length($0)-2)} NR % 2 == 0' $FILE
    fi
done



