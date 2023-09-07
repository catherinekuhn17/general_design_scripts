#!/bin/bash

#$ -S /bin/bash
#$ -o log_files/
#$ -cwd
#$ -j y
#$ -r y

PDB1="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/3C62_chainD_40_5-8_40/3C62_chainD_40_5-8_40_0_1.pdb"
PDB2="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/3C62_chainD_40_5-8_40/3C62_chainD_40_5-8_40_0_2.pdb"
OUT_FN="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/3C62_chainD_40_5-8_40/3C62_chainD_40_5-8_40_0_comb.pdb"

python comb_pdbs.py --pdb_fn1=$PDB1 --pdb_fn2=$PDB2 --out_fn=$OUT_FN