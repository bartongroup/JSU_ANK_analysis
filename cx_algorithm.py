import prointvar
from prointvar.pdbx import PDBXreader
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import *
import os

class RightChain(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
        
    def accept_chain(self, chain):
        if chain.get_id() == self.chain_id:
            return 1
        else:
            return 0

def split_pdb_files(df, pdb_dir, pdb_split_dir):
    strucs = df.PDB_dbAccessionId_A.unique().tolist()
    if not os.path.isdir(pdb_split_dir):
        os.mkdir(pdb_split_dir)
    for struc in strucs:
        pdb_in = os.path.join(pdb_dir, "{}.pdb".format(struc))
        if os.path.isfile(pdb_in):
            parser = PDBParser()
            structure = parser.get_structure(struc, pdb_in)
            for chain in structure.get_chains():
                pdb_out = os.path.join(pdb_split_dir, "pdb{}_{}.ent".format(structure.id, chain.id))
                if os.path.isfile(pdb_out):
                    print("{} already available!".format(pdb_out))
                else:
                    io = PDBIO()
                    io.set_structure(structure)
                    io.save(pdb_out, RightChain(chain.id))
        else:
            print("{} not found!".format(pdb_in))
    print("All PDBs have been split!")        
        
def get_cx_score(natom, vatom = 20.1, r = 20):
    v_int = natom*vatom
    v_sphere = (4/3)*math.pi*r**3
    v_ext = v_sphere - v_int
    cx = (v_ext/v_int)
    return cx

def get_cx_scores(input_struct, r = 20):
    filename = input_struct.split('/')[-1]
    pdb_code = filename[3:7]
    df = PDBXreader(inputfile=input_struct).atoms(format_type="pdb", add_contacts=True, residue_agg=False, dist = r)
    df_atom = df[(df.group_PDB == 'ATOM')&(df.label_atom_id != "H")]
    if len(df_atom) == 0:
        print("No ATOM records found in {}".format(pdb_code))
        return df_atom
    df_atom['contact_indexes_list'] = df_atom['contact_indexes'].str.split(',')
    df_atom['contact_indexes_list'] = df_atom['contact_indexes_list'].apply(lambda x: map(int, x)).apply(list)
    df_atom['neighbors'] = df_atom['contact_indexes_list'].apply(len)
    df_atom['cx_atom'] = df_atom.neighbors.apply(get_cx_score, args=(20.1,r))
    df_atom.label_seq_id = df_atom.label_seq_id.astype(int)
    df_atom = df_atom.sort_values(by=['label_seq_id'])
    df_atom['res_cx_min'] = df_atom.groupby('label_seq_id')['cx_atom'].transform('min')
    df_atom['res_cx_max'] = df_atom.groupby('label_seq_id')['cx_atom'].transform('max')
    df_atom['res_cx_median'] = df_atom.groupby('label_seq_id')['cx_atom'].transform('median')
    df_atom['res_cx_mean'] = df_atom.groupby('label_seq_id')['cx_atom'].transform('mean')
    df_atom = df_atom.drop_duplicates(['label_asym_id','label_seq_id'])
    df_atom['pdb_id'] = pdb_code
    df_atom_cropped = df_atom[['pdb_id','label_asym_id','label_seq_id','res_cx_min','res_cx_mean','res_cx_median','res_cx_max']]
    df_atom_cropped.index 
    return df_atom_cropped
    
def get_cx_df(df, csv_dir, files, out_df, r = 20):
    if not os.path.isdir(csv_dir):
        os.mkdir(csv_dir)
    cx_dfs = []
    files_left = len(files)
    for file in files:
        struc_id = file.split("/")[-1].split(".")[0][3:]
        out_csv = os.path.join(csv_dir, "{}.csv".format(struc_id))
        if os.path.isfile(out_csv):
            struc_df = pd.read_csv(out_csv)
            print("{} already available!".format(out_csv))
        else:
            struc_df = get_cx_scores(file, r = r)
            if len(struc_df) == 0:
                files_left -= 1
                continue
            struc_df.to_csv(out_csv, index = False)
        cx_dfs.append(struc_df)
        files_left -= 1
        print("{} has been processed successfully! {} to go!".format(file, files_left))
    cx_all_strucs = pd.concat(cx_dfs)
    cx_all_strucs_sorted = cx_all_strucs.sort_values(by=['pdb_id','label_asym_id','label_seq_id'])
    cx_all_strucs_sorted.index = range(0,len(cx_all_strucs_sorted))
    cx_all_strucs_sorted.to_csv(out_df, index = False)
    print('{} was correctly saved'.format(out_df))
    return cx_all_strucs_sorted

def merge_with_cx(structure_df, cx_df):
    r_keys = ['pdb_id', 'label_seq_id', 'label_asym_id']
    structure_df = structure_df.astype(str)
    cx_df = cx_df.drop_duplicates(subset = r_keys).astype(str)
    merged = pd.merge(structure_df, cx_df, left_on = ["PDB_dbAccessionId_A", "PDB_dbResNum_A", "PDB_dbChainId_A"], right_on = r_keys, how = 'left')
    return merged