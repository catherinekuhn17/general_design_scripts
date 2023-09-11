#!/bin/bash

#$ -S /bin/bash
#$ -o log_files/
#$ -cwd
#$ -j y
#$ -r y

FOLDER='/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/random_placements/testing/beta_TET_11_AP_tests/beta_TET_11_AP_place_20'
TEMP_FN='/wynton/home/kortemme/ckuhn/Desktop/kortemme/metal_binding/inputs/switch_temps/beta_TET_11_AP.pdb'
OUT_FN='/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/random_placements/testing/beta_TET_11_AP_tests/rmsd_out.csv'

python get_motif_rmsd.py --folder $FOLDER --temp_fn $TEMP_FN --out_fn $OUT_FN