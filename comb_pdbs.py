import numpy as np
from argparse import ArgumentParser


def get_args(argv=None):
    p = ArgumentParser(description=__doc__)
    p.add_argument("--pdb_fn1", type=str, help="path to pdb 1")
    p.add_argument("--pdb_fn2", type=str, help="path to pdb 1")
    p.add_argument("--out_fn", type=str, help="fn of combined pdb")
    args = p.parse_args()
    if argv is not None:
        args = p.parse_args(argv) # for use when testing
    else:
        args = p.parse_args()
    return args


'''
DEFINE GLOBAL VARIABLES
'''

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

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_N = {a:n for n,a in enumerate(alpha_1)}

aa123 = {aa1: aa3 for aa1, aa3 in zip(alpha_1, num2aa)}
aa321 = {aa3: aa1 for aa1, aa3 in zip(alpha_1, num2aa)}


'''
FUNCTIONS FOR READING PDB FILE
'''

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

'''
FUNCTIONS FOR WRITING PDB FILE
'''

def write_pdb_string(xyz, res, pdb_idx=None):
    # xyz: (L, 3, 3)
    if pdb_idx is None:
        pdb_idx = [('A', i+1) for i in range(len(res))]

    idx = np.arange(len(res))

    wrt = ""
    atmNo = 0
    for i_res, (ch,i_pdb), aa in zip(idx, pdb_idx, res):
        for i_atm, atm in enumerate(["N", "CA", "C", "O"]):
            atmNo += 1
            wrt += "ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(
                int(atmNo), atm, aa123[aa], ch, int(i_pdb), xyz[i_res,i_atm,0], 
                xyz[i_res,i_atm,1], xyz[i_res,i_atm,2], 1.0, 0.0)

    return wrt
    

def write_pdb(xyz, prefix, res, pdb_idx=None, comments=None, sc=False, Bfacts=None, h=None):

    if sc:
        wrt = write_pdb_str_sc(res, xyz, Bfacts=Bfacts, pdb_idx=pdb_idx,h=h)
    else:
        wrt = write_pdb_string(xyz, res, pdb_idx)

    with open(prefix, 'wt') as fp:
        fp.write(wrt)
        fp.write("\n\n\n\n")
        
        if comments is not None:
            for line in comments:
                fp.write(line+'\n')
    
    return

def main():
    '''
    where everything is loaded in and calculated
    '''
    
    # load in arguments
    args = get_args()
    
    # load in pdbs
    pdb1=parse_pdb(args.pdb_fn1)
    pdb2=parse_pdb(args.pdb_fn2)
    
    # make pdb2 into chain B:
    pdb2['pdb_idx'] = [('B', e[1]) for e in pdb2['pdb_idx']]
    # translate pdb2 by (300,0,0)
    translate = np.ones((14, 3))*[300,0,0]
    pdb2['xyz'] =  pdb2['xyz'] + translate

    # combine pdbs into one pdb
    pdb_comb = {}
    xyz_list = [np.array(pdb1['xyz']), np.array(pdb2['xyz'])]
    pdb_idx_list = [np.array(pdb1['pdb_idx']), np.array(pdb2['pdb_idx'])]
    seq_list = [np.array(pdb1['seq']), np.array(pdb2['seq'])]

    pdb_comb['xyz'] = np.concatenate(xyz_list, axis=0)
    pdb_comb['pdb_idx'] = np.concatenate(pdb_idx_list, axis=0)
    pdb_comb['seq'] = np.concatenate(seq_list, axis=0)
    
    # write out pdb
    write_pdb(pdb_comb['xyz'], 
              args.out_fn, 
              [alpha_1[e] for e in pdb_comb['seq']], 
              pdb_idx = pdb_comb['pdb_idx']) 

if __name__ == "__main__":
    main()