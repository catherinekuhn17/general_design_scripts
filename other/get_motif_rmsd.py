import sys, glob, os, shutil
import numpy as np
import pandas as pd
import pickle
from argparse import ArgumentParser
from distutils.util import strtobool

def get_args(argv=None):
    p = ArgumentParser(description=__doc__)
    p.add_argument("--folder", type=str, help="path to folder you want to run analysis on (with pdb and trb files)")
    p.add_argument("--temp_fn", type=str, help="fn of template pdb with motif")
    p.add_argument("--out_fn", type=str, help="fn to write output csv to")
    args = p.parse_args()
    if argv is not None:
        args = p.parse_args(argv) # for use when testing
    else:
        args = p.parse_args()
    return args

num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    ]
aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE2"," CE3"," NE1"," CZ2"," CZ3"," CH2"), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None), # val
]

aa2num= {x:i for i,x in enumerate(num2aa)}
def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res = [(l[22:26],l[17:20]) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
    seq = [aa2num[r[1]] if r[1] in aa2num.keys() else 20 for r in res]
    pdb_idx = [( l[21:22].strip(), int(l[22:26].strip()) ) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]  # chain letter, res num
    plddt_val = [(float(l[61:66].strip()) ) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]  # chain letter, res num

    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(res), 14, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' '+l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain,resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm is not None and tgtatm.strip() == atom.strip(): # ignore whitespace
                xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break

    # save atom mask
    mask = np.logical_not(np.isnan(xyz[...,0]))
    xyz[np.isnan(xyz[...,0])] = 0.0

    # remove duplicated (chain, resi)
    new_idx = []
    i_unique = []
    for i,idx in enumerate(pdb_idx):
        if idx not in new_idx:
            new_idx.append(idx)
            i_unique.append(i)
            
    pdb_idx = new_idx
    xyz = xyz[i_unique]
    mask = mask[i_unique]
    seq = np.array(seq)[i_unique]

    out = {'xyz':xyz, # cartesian coordinates, [Lx14]
            'mask':mask, # mask showing which atoms are present in the PDB file, [Lx14]
            'idx':np.array([i[1] for i in pdb_idx]), # residue numbers in the PDB file, [L]
            'seq':np.array(seq), # amino acid sequence, [L]
            'pdb_idx': pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
            'plddt_val':plddt_val
           }

    # heteroatoms (ligands, etc)
    if parse_hetatom:
        xyz_het, info_het = [], []
        for l in lines:
            if l[:6]=='HETATM' and not (ignore_het_h and l[77]=='H'):
                info_het.append(dict(
                    idx=int(l[7:11]),
                    atom_id=l[12:16],
                    atom_type=l[77],
                    name=l[16:20]
                ))
                xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

        out['xyz_het'] = np.array(xyz_het)
        out['info_het'] = info_het

    return out
def parse_pdb(filename, **kwargs):
    '''extract xyz coords for all heavy atoms'''
    lines = open(filename,'r').readlines()
    return parse_pdb_lines(lines, **kwargs)

def calc_rmsd(xyz1, xyz2, eps=1e-6):

    # center to CA centroid
    xyz1 = xyz1 - xyz1.mean(0)
    xyz2 = xyz2 - xyz2.mean(0)

    # Computation of the covariance matrix
    C = xyz2.T @ xyz1

    # Compute optimal rotation matrix using SVD
    V, S, W = np.linalg.svd(C)

    # get sign to ensure right-handedness
    d = np.ones([3,3])
    d[:,-1] = np.sign(np.linalg.det(V)*np.linalg.det(W))

    # Rotation matrix U
    U = (d*V) @ W

    # Rotate xyz2
    xyz2_ = xyz2 @ U

    L = xyz2_.shape[0]
    rmsd = np.sqrt(np.sum((xyz2_-xyz1)*(xyz2_-xyz1), axis=(0,1)) / L + eps)

    return rmsd    

def main():
    '''
    where everything is loaded in and calculated
    '''
    # load in arguments
    args = get_args()
    folder = args.folder
    temp_fn =args.temp_fn
    out_fn = args.out_fn
    
    # load in fns
    pdb_fns = glob.glob(f'{folder}/*pdb')
    temp_pdb = parse_pdb(temp_fn)
    
    # build RMSD dict
    out_dict={}
    for fn in pdb_fns:
        prefix =fn.split('/')[-1].split('.pdb')[0]
        itter_num = int(prefix.split('_')[-1])
        out_dict[itter_num]={}
        out_dict[itter_num]['file_name'] = prefix
        #load in pdb
        model_pdb = parse_pdb(fn)
        # load in trb file
        trb_fn = f'{fn[0:-4]}.trb'
        trb_pkl = pd.read_pickle(trb_fn)
        # get motif residue info
        temp_resi = trb_pkl['con_ref_pdb_idx']
        model_resi = trb_pkl['con_hal_pdb_idx']
        # find idx of motif residues in template
        temp_motif_resi_idx=[]
        for idx, resi in enumerate(temp_pdb['pdb_idx']):
            if resi in temp_resi:
                temp_motif_resi_idx.append(idx)
        # find idx of motif residues in diffused model
        model_motif_resi_idx=[]
        for idx, resi in enumerate(model_pdb['pdb_idx']):
            if resi in model_resi:
                model_motif_resi_idx.append(idx)
        #find RMSD 
        xyz_temp = temp_pdb['xyz'][temp_motif_resi_idx]
        xyz_model = model_pdb['xyz'][model_motif_resi_idx]

        #bb only 
        bb_rmsd = calc_rmsd(xyz_temp[:,:4].reshape(len(xyz_temp)*4,3), 
                            xyz_model[:,:4].reshape(len(xyz_temp)*4,3))
        out_dict[itter_num]['bb_rmsd'] = bb_rmsd

    # turn dictionary into a datafram to be saved
    keys = list(out_dict.keys())
    keys.sort()
    sorted_dict = {i: out_dict[i] for i in keys}
    df_rmsd=pd.DataFrame()
    df_rmsd['itter'] = sorted_dict.keys()
    df_rmsd['file_name'] = [d['file_name'] for d in sorted_dict.values()]
    df_rmsd['bb_rmsd'] = [d['bb_rmsd'] for d in sorted_dict.values()]
    df_rmsd.to_csv(out_fn)

if __name__ == "__main__":
    main()