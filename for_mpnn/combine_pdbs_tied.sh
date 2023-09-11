#!/bin/bash

#$ -S /bin/bash
#$ -o log_files/
#$ -cwd
#$ -j y
#$ -r y


#PDB1="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/3C62_chainD_40_5-8_40/3C62_chainD_40_5-8_40_0_1.pdb"
#PDB2="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/3C62_chainD_40_5-8_40/3C62_chainD_40_5-8_40_0_2.pdb"
#OUT_FN="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/3C62_chainD_40_5-8_40/3C62_chainD_40_5-8_40_0_tied.pdb"
PREFIX="3C62_chainD_rog_40_5-8_60"
FOLDER_WITH_PDBS="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/3C62_chainD_40_5-8_40/*pdb"

for PDB_FN in $FOLDER_WITH_PDBS; do
    if [[ $PDB_FN == *"1.pdb"* ]]; then
        PDB1=$PDB_FN
        PRE=${PDB_FN:0:-5}
        PDB2="${PRE}2.pdb"
        OUT="${PRE}tied.pdb"
        python combine_pdbs_tied.py --pdb_fn1=$PDB1 --pdb_fn2=$PDB2 --out_fn=$OUT_FN

    fi
done
# python combine_pdbs_tied.py --pdb_fn1=$PDB1 --pdb_fn2=$PDB2 --out_fn=$OUT_FN
