import proteofav
from proteofav import validation
import numpy as np
import pandas as pd
import os
import re

def download_validation(df, validation_dir):
    structures = df.PDB_dbAccessionId_A.unique().tolist()
    for struc in structures:
        out_file = os.path.join(validation_dir, "{}.csv".format(struc))
        if os.path.isfile(out_file):
            print("Validation for {} already exists!".format(struc))
            continue
        else:
            validation_table = validation.select_validation(struc)
            validation_table.to_csv(out_file, index = False)

def get_rsrz_dict(df, validation_dir):
    structures = df.PDB_dbAccessionId_A.unique().tolist()
    rsrz_dict = {}
    for struc in structures:
        table = pd.read_csv(os.path.join(validation_dir, "{}.csv".format(struc)))
        if "validation_rsrz" in list(table.columns):
            rsrz_dict[struc] = True
        else:
            rsrz_dict[struc] = False
    return rsrz_dict

def get_res_exp_from_df(df, pdb_dir):
    structures = df.PDB_dbAccessionId_A.unique().tolist()
    resolutions = {}
    p_res = re.compile("^REMARK\s*\d*\s*RESOLUTION\.\s*(\d*\.\d*)\s*ANGSTROMS\.\s*$")
    p_exp = re.compile("^EXPDTA\s*(\w-?\w*\s?\w*)\s*$")
    for structure in structures:
        structure_path = os.path.join(pdb_dir, "{}.pdb".format(structure))
        if os.path.isfile(structure_path):
            resolutions[structure] = {}
            resolutions[structure]['experiment'] = None
            resolutions[structure]['resolution'] = None
            pdb_fh = open(structure_path, 'r')
            for line in pdb_fh:
                m_res = p_res.match(line)
                m_exp = p_exp.match(line)
                if m_exp:
                    e = m_exp.groups()[0]
                    resolutions[structure]['experiment'] = e
                if m_res:
                    r = m_res.groups()[0]
                    resolutions[structure]['resolution'] = float(r)
                    break
            pdb_fh.close()
            print("{} processed successfully!".format(structure))
        else:
            print("{} file not found!".format(structure_path))
    struc_df = pd.DataFrame.from_dict(resolutions, orient = "index")
    struc_df["pdb_id"] = struc_df.index
    struc_df.index = range(0, len(struc_df))
    return struc_df

def get_validation_df(df, validation_dir):
    structures = df.PDB_dbAccessionId_A.unique().tolist()
    validation_dfs = []
    for struc in structures:
        table = pd.read_csv(os.path.join(validation_dir, "{}.csv".format(struc)))
        table["pdb_id"] = struc
        validation_dfs.append(table)
        print("{} validation appended successfully!".format(struc))
    validation_master_df = pd.concat(validation_dfs)
    return validation_master_df

def merge_with_validation(structure_df, validation_df):
    r_keys = ['pdb_id', 'validation_resnum', 'validation_chain']
    structure_df = structure_df.astype(str)
    validation_df = validation_df.drop_duplicates(subset = r_keys).astype(str)
    merged = pd.merge(structure_df, validation_df, left_on = ["PDB_dbAccessionId_A", "PDB_dbResNum_A", "PDB_dbChainId_A"], right_on = r_keys, how = 'left')
    merged = pd.merge(merged, validation_df, left_on = ["PDB_dbAccessionId_B", "PDB_dbResNum_B", "PDB_dbChainId_B"], right_on = r_keys, suffixes = ('_A', '_B'), how = 'left')
    return merged

def get_prot_prot_rows(df):
    return df[df.interaction_type == "Protein-Protein"]

def filter_rsrz_rscc(df, rsrz_t = 2, rscc_t = 0.85):
    df.validation_rsrz_A = df.validation_rsrz_A.astype(float)
    df.validation_rsrz_B = df.validation_rsrz_B.astype(float)
    df.validation_rsr_A = df.validation_rsr_A.astype(float)
    df.validation_rsr_B = df.validation_rsr_B.astype(float)
    df.validation_rscc_A = df.validation_rscc_A.astype(float)
    df.validation_rscc_B = df.validation_rscc_B.astype(float)
    return df[(df.validation_rscc_A > rscc_t) & (df.validation_rsrz_A < rsrz_t) & (df.validation_rscc_B > rscc_t) & (df.validation_rsrz_B < rsrz_t)]