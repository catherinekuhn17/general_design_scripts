#!/bin/bash
#$ -S /bin/bash
#$ -o log_files/
#$ -cwd
#$ -j y
#$ -r y 
#conda activate mpnn_env
conda activate SE3nv
cd /wynton/home/kortemme/ckuhn/Desktop/programs/ProteinMPNN
PREFIX="3C62_chainD_rog_40_5-8_60"
FOLDER_WITH_PDBS="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/$PREFIX/*pdb"
TMP_OUTPUT_DIR="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/$PREFIX"
path_for_parsed_chains=$TMP_OUTPUT_DIR"/parsed_pdbs.jsonl"
path_for_assigned_chains=$TMP_OUTPUT_DIR"/assigned_pdbs.jsonl"
path_for_fixed_positions=$TMP_OUTPUT_DIR"/fixed_pdbs.jsonl"
chains_to_design="A"
#OUT_DIR="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/random_placements/initial_diff/$PREFIX"
OUT_DIR="/wynton/home/kortemme/ckuhn/Desktop/thesis_proj/diffusion_models/switch_diff/$PREFIX"
for PDB_FN in $FOLDER_WITH_PDBS; do
	i=0
        # determine where histidines are so they are kept the same!
	for line in  `awk '/N   HIS/ {print $6}' $PDB_FN`; do
		arr[$i]=$line
		i=`expr $i + 1`
	done
	fixed_positions=${arr[*]}
	echo $fixed_positions 
	echo $PDB_FN
	echo $path_for_parsed_chains
        python helper_scripts/parse_multiple_chains_pdb.py --input_path=$PDB_FN --output_path=$path_for_parsed_chains
        python helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
	python helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"
 
	python protein_mpnn_run.py \
		--omit_AAs 'HC' \
		--jsonl_path $path_for_parsed_chains \
		--chain_id_jsonl $path_for_assigned_chains \
		--fixed_positions_jsonl $path_for_fixed_positions \
		--out_folder $OUT_DIR \
		--num_seq_per_target 1 \
		--sampling_temp "0.2" \
		--seed 37 \
		--use_soluble_model  \
		--batch_size 1 \
		--use_soluble_model=True

done

