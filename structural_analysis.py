import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import Bio.SeqIO
import Bio.AlignIO
import statistics
import importlib
from Bio.PDB import *
import os
import pandas as pd
import matplotlib.patches as mpatches

def get_struc_info(df):
    n_res = len(df.drop_duplicates(['UniProt_dbAccessionId_A', 'PDB_dbAccessionId_A', 'PDB_dbChainId_A', 'PDB_dbResNum_A']))
    n_res_un = len(df.drop_duplicates(["SOURCE_ID_A", "Alignment_column_A"]))
    n_strucs =  len(df.PDB_dbAccessionId_A.unique().tolist())
    n_prots = len(df.UniProt_dbAccessionId_A.unique().tolist())
    n_reps = len(df.SOURCE_ID_A.unique().tolist())
    print("The dataframe contains information of:\n{} residues, {} of which are unique\n{} PDB structures\n{} different proteins\n{} unique repeats".format(n_res, n_res_un, n_strucs, n_prots, n_reps))
    
def get_rsa_class_consensus(df, aln_cols):
    """
    """
    df_unique = df.drop_duplicates(['UniProt_dbAccessionId_A', 'PDB_dbAccessionId_A', 'PDB_dbChainId_A', 'PDB_dbResNum_A'])
    df_unique = df_unique.dropna(subset = ["RSA_CLASS_UNB_A", "RSA_UNB_A"])
    df_unique.Alignment_column_A = df_unique.Alignment_column_A.astype(int)
    df_cons_cols = df_unique[df_unique.Alignment_column_A.isin(aln_cols)]
    get_struc_info(df_cons_cols)
    grouped_col_rep = df_cons_cols.groupby(["Alignment_column_A", "SOURCE_ID_A"])
    rsa_class_dict = {}
    rsa_dict = {}
    i = 1
    for k, v in grouped_col_rep:
        col = k[0]
        if col not in rsa_class_dict:
            print("Processing consensus column {}".format(i))
            i += 1
            rsa_class_dict[col] = []
        if col not in rsa_dict:
            rsa_dict[col] = []
        rsa_list = list(map(float, v.RSA_UNB_A.tolist()))
        rsa_class_list = v.RSA_CLASS_UNB_A.tolist()
        try:
            rsa_class_cons = statistics.mode(rsa_class_list)
        except:
            rsa_class_cons = rsa_class_list[0]
        rsa_median = statistics.median(rsa_list)
        rsa_class_dict[col].append(rsa_class_cons)
        rsa_dict[col].append(rsa_median)
    return rsa_class_dict, rsa_dict

def get_rsa_class_df(rsa_class_dict):
    n = len(rsa_class_dict.keys())
    core = []
    surf = []
    part = []
    for v in rsa_class_dict.values():
        core.append(v.count("Core"))
        surf.append(v.count("Surface"))
        part.append(v.count("Part. Exposed")) 
    df_dssp = pd.DataFrame(list(zip(core, surf, part)), columns =  ["core", "surf", "part"])
    df_dssp["tot"] = df_dssp.core + df_dssp.surf + df_dssp.part
    df_dssp["p_core"] = df_dssp.core / df_dssp.tot
    df_dssp["p_surf"] = df_dssp.surf / df_dssp.tot
    df_dssp["p_part"] = df_dssp.part / df_dssp.tot
    df_dssp["p_tot"] = df_dssp.p_core + df_dssp.p_surf + df_dssp.p_part
    df_dssp.index = range(1, n + 1)
    return df_dssp

def plot_rsa_consensus(df_dssp, palette = ["darkgreen", "red", "darkred"], bwidth = 1.2, out = None):
    core = list(df_dssp.p_core)
    surf = list(df_dssp.p_surf)
    part = list(df_dssp.p_part)
    n = len(df_dssp)
    n_cols = len(df_dssp)
    bottom = []
    for i in range(0, len(core)):
        bottom.append(core[i] + part[i])

    r = list(np.arange(1, n_cols*1.5, 1.5))
    
    plt.figure(figsize=(180,80))
    plt.rcParams.update({"axes.linewidth": 10})
    plt.bar(r, core,  color=palette[0], edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Core')
    plt.bar(r, part, color=palette[1], bottom = core, edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Part')
    plt.bar(r, surf,  bottom = bottom, color=palette[2], edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Surf')
    #plt.title('DSSP RSA classification', pad = 100, fontsize = 140)
    plt.xlabel('Consensus residue position', labelpad = 100, fontsize = 160)
    plt.ylabel('p', labelpad = 100, fontsize = 160)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 80)
    legend = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 120)
    legend.get_frame().set_linewidth(10)
    legend.get_frame().set_edgecolor("black")
    for legobj in legend.legendHandles:
        legobj.set_linewidth(7.5)
    plt.tick_params(axis= 'both' , which = 'major', pad = 60, width = 15, length = 50, labelsize = 140)
    plt.xticks(np.arange(1, n_cols*1.5, 1.5), list(range(1, n + 1)))
    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.axhline(y=0.5, linewidth = 5, linestyle = "--", color = 'black')
    plt.xlim(0, (n_cols*1.5)+0.5)
    plt.ylim(0, 1.025)
    if out != None:
        plt.savefig(out)
    plt.show()

def plot_rsa_median(rsa_dict, out = None):
    n = len(rsa_dict.keys())
    ks = list(range(1, n + 1))
    vs = [statistics.median(v) for v in rsa_dict.values()]
    
    palette = {}
    c = 1
    for v in vs:
        if v <= 5:
            palette[c] = "royalblue"
        elif v > 25:
            palette[c] = "firebrick"
        else:
            palette[c] = "tomato"
        c +=1
        
    legend_dict = { 'Core' : 'royalblue', 'Part' : 'tomato', 'Surf' : 'firebrick' }
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(facecolor=legend_dict[key], label=key, edgecolor = "black", linewidth = 3)
        patchList.append(data_key)
        
    i = 1
    lowerbound = []
    upperbound = []
    for v in rsa_dict.values():
        med = statistics.median(v)
        data = pd.Series(v)
        low, up = medianCI(data, 0.95, 0.5)
        lowerbound.append(med-low)
        upperbound.append(up-med)
        i = i+1
    plt.figure(figsize=(180,80))
    plt.rcParams.update({"axes.linewidth": 10})
    ax=sns.barplot(ks, vs, edgecolor = 'black', palette = palette,linewidth = 7.5)#,
                   #ci=[lowerbound, upperbound], errcolor = "black", errwidth = 7.5, capsize = 35)
    plt.errorbar([k-1 for k in ks], vs, yerr=[lowerbound, upperbound], c = "black", linewidth = 10, linestyle="None", capsize = 35.0, capthick = 7.5)
    #ax.set_title('Median RSA per consensus position', pad =100, fontsize = 140)
    ax.tick_params(axis= 'both' , which = 'major', pad = 60, width = 15, length = 50, labelsize = 140)
    ax.set_xlabel("Domain position", labelpad = 100, fontsize = 160)
    ax.set_ylabel('Median RSA', labelpad = 100, fontsize = 160)
    ax.axhline(25, linewidth = 7.5, color = 'black', linestyle = '--')
    ax.axhline(5, linewidth = 7.5, color = 'black', linestyle = '--')
    plt.xlim(-0.7,32.7)
    plt.yticks(np.arange(0, 75, 5))
    legend = plt.legend(handles=patchList, loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 120)
    legend.get_frame().set_linewidth(10)
    legend.get_frame().set_edgecolor("black")
    for legobj in legend.legendHandles:
        legobj.set_linewidth(7.5)
    if out != None:
        plt.savefig(out)
    plt.show()
    return vs

def medianCI(data, ci, p):
    '''
    data: pandas datafame/series or numpy array
    ci: confidence level
    p: percentile' percent, for median it is 0.5
    output: a list with two elements, [lowerBound, upperBound]
    '''
    if type(data) is pd.Series or type(data) is pd.DataFrame:
        #transfer data into np.array
        data = data.values

    #flat to one dimension array
    data = data.reshape(-1)
    data = np.sort(data)
    N = data.shape[0]
    
    lowCount, upCount = stats.binom.interval(ci, N, p, loc=0)
    #given this: https://onlinecourses.science.psu.edu/stat414/node/316
    #lowCount and upCount both refers to  W's value, W follows binomial Dis.
    #lowCount need to change to lowCount-1, upCount no need to change in python indexing
    lowCount -= 1
    # print lowCount, upCount
    return data[int(lowCount)], data[int(upCount)]

def get_rsa_cons_dict(rsa_class_dict):
    c = 1
    rsa_cons_dict1 = {}
    rsa_cons_dict2 = {}
    for v in rsa_class_dict.values():
        mode = statistics.mode(v)
        if mode == 'Part. Exposed':
            rsa_cons_dict1[c] = mode
            rsa_cons_dict2[c] = "Surface"
        else:
            rsa_cons_dict1[c] = mode
            rsa_cons_dict2[c] = mode
        c += 1
    return rsa_cons_dict1, rsa_cons_dict2

def add_rsa_class(df, rsa_cons_dict):
    df['rsa_cons'] = df.index.map(rsa_cons_dict)
    colors_rsa = {}
    for k, v in rsa_cons_dict.items():
        if v == 'Core':
            colors_rsa[k] = 'mediumblue'
        elif v == 'Surface':
            colors_rsa[k] = 'darkred'
        elif v == 'Part. Exposed':
            colors_rsa[k] = 'red'
    df['color_rsa'] = df.index.map(colors_rsa)
    return df

def plot_rsa_cx(df, out = None):
    plt.figure(figsize=(140,60))
    plt.rcParams.update({"axes.linewidth": 5})
    colours = ["royalblue", "tomato", "firebrick"]
    palette = {}
    for i in df.index:
        rsa = df.loc[i,"rsa"]
        if rsa <= 5:
            palette[i] = colours[0]
        elif rsa >= 25:
            palette[i] = colours[2]
        else:
            palette[i] = colours[1]
    ax = sns.barplot(df.index, df.rsa, edgecolor = 'black',linewidth = 5, palette = palette, label = "RSA")
    ax.set_title('Median RSA and CX score', pad =100, fontsize = 140)
    ax.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 100)
    ax.set_xlabel("Consensus residue position", labelpad = 100, fontsize = 120)
    ax.set_ylabel("Median RSA", labelpad = 100, fontsize = 120)
    ax.plot([], [], color='black', lw=15, label='CX value')
    ax.axhline(y = 5, color = "black", linewidth = 5, linestyle = '--')
    ax.axhline(y = 25, color = "black", linewidth = 5, linestyle = '--')
    #plt.xticks(np.arange(1,33*1.5,1.5),list(range(1, 33 + 1)))
    ax2 = ax.twinx()
    sns.lineplot(df.index-1, df.cx, color = 'black',linewidth = 15)
    sns.scatterplot(df.index-1, df.cx, color = 'black', s = 3500, edgecolor = 'black', linewidth = 5)
    ax2.set_ylabel("CX value", fontsize = 100, labelpad = 100)
    #ax2.set_yticks(np.arange(0,125,25))
    ax2.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 80)
    #plt.xlim(0,(33*1.5)+0.5)
    plt.xlim(-0.7,32.7)
    ax.legend(loc='upper left', fontsize = 80)
    if out != None:
        plt.savefig(out)
    plt.show()
    
def get_ss_class_consensus(df, aln_cols):
    """
    """
    df_unique = df.drop_duplicates(['UniProt_dbAccessionId_A', 'PDB_dbAccessionId_A', 'PDB_dbChainId_A', 'PDB_dbResNum_A'])
    df_unique = df_unique.dropna(subset = ["SS_CLASS_A"])
    df_unique.Alignment_column_A = df_unique.Alignment_column_A.astype(int)
    df_cons_cols = df_unique[df_unique.Alignment_column_A.isin(aln_cols)]
    get_struc_info(df_cons_cols)
    grouped_col_rep = df_cons_cols.groupby(["Alignment_column_A", "SOURCE_ID_A"])
    ss_class_dict = {}
    i = 1
    for k, v in grouped_col_rep:
        col = k[0]
        if col not in ss_class_dict:
            print("Processing consensus column {}".format(i))
            i += 1
            ss_class_dict[col] = []
        ss_class_list = v.SS_CLASS_A.tolist()
        try:
            ss_class_cons = statistics.mode(ss_class_list)
        except:
            ss_class_cons = ss_class_list[0]
        ss_class_dict[col].append(ss_class_cons)
    return ss_class_dict

def get_ss_class_df(rsa_class_dict):
    n = len(rsa_class_dict.keys())
    h = []
    c = []
    e = []
    for v in rsa_class_dict.values():
        h.append(v.count("H"))
        c.append(v.count("C"))
        e.append(v.count("E")) 
    df_dssp = pd.DataFrame(list(zip(h, c, e)), columns =  ["helix", "coil", "strand"])
    df_dssp["tot"] = df_dssp.helix + df_dssp.coil + df_dssp.strand
    df_dssp["p_helix"] = df_dssp.helix / df_dssp.tot
    df_dssp["p_coil"] = df_dssp.coil / df_dssp.tot
    df_dssp["p_strand"] = df_dssp.strand / df_dssp.tot
    df_dssp["p_tot"] = df_dssp.p_helix + df_dssp.p_coil + df_dssp.p_strand
    df_dssp.index = range(1, n + 1)
    return df_dssp

def plot_ss_class_consensus(df_dssp, palette = ["mediumseagreen", "tomato", "cornflowerblue"], bwidth = 1.2, out = None):
    helix = list(df_dssp.p_helix)
    coil = list(df_dssp.p_coil)
    strand = list(df_dssp.p_strand)

    n_cols = len(df_dssp)
    bottom = []
    for i in range(0, len(helix)):
        bottom.append(strand[i]+helix[i])

    r = list(np.arange(1,n_cols*1.5,1.5))
    plt.figure(figsize=(180,80))
    plt.rcParams.update({"axes.linewidth": 10})
    plt.bar(r, strand,  color=palette[0], edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Strand')
    plt.bar(r, helix, color=palette[1], bottom = strand, edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Helix')
    plt.bar(r, coil,  bottom = bottom, color=palette[2], edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Loop')
    #plt.title('DSSP Secondary structure classification', pad = 100, fontsize = 140)
    plt.xlabel('Domain position', labelpad = 100, fontsize = 160)
    plt.ylabel('p', labelpad = 100, fontsize = 160)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 80)
    legend = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 120)
    legend.get_frame().set_linewidth(10)
    legend.get_frame().set_edgecolor("black")
    for legobj in legend.legendHandles:
        legobj.set_linewidth(7.5)
    plt.tick_params(axis= 'both' , which = 'major', pad = 60, width = 15, length = 50, labelsize = 140)
    plt.axhline(y=0.5, linewidth = 5, linestyle = "--", color = 'black')
    plt.xticks(np.arange(1,n_cols*1.5,1.5),list(range(1, n_cols + 1)))
    plt.yticks(np.arange(0,1.1,0.1))
    plt.xlim(0,(n_cols*1.5)+0.5)
    plt.ylim(0, 1.025)
    if out != None:
        plt.savefig(out)
    plt.show()

def get_ss_cons_dict(ss_dict):
    c = 1
    ss_cons_dict1 = {}
    ss_cons_dict2 = {}
    h = 0
    for v in ss_dict.values():
        mode = statistics.mode(v)
        if mode == 'H':
            if ss_cons_dict1[c-1] == "C":
                h += 1
            ss_cons_dict2[c] = mode+str(h)
            ss_cons_dict1[c] = mode
        else:
            ss_cons_dict1[c] = mode
            ss_cons_dict2[c] = mode
        c += 1
    return ss_cons_dict1, ss_cons_dict2

def add_ss_class(df, ss_cons_dict, palette = ["mediumseagreen", "firebrick", "cornflowerblue"]):
    df['ss_cons'] = df.index.map(ss_cons_dict)
    colors_ss = {}
    for k, v in ss_cons_dict.items():
        if v == 'H' or v == 'H1':
            colors_ss[k] = palette[1]
        elif v == 'H2':
            colors_ss[k] = 'orange'
        elif v == 'C':
            colors_ss[k] = palette[2]
        elif v == 'E':
            colors_ss[k] = palette[0]
    df['color_ss'] = df.index.map(colors_ss)
    return df    

def get_ss_consensus(df, aln_cols):
    """
    """
    df_unique = df.drop_duplicates(['UniProt_dbAccessionId_A', 'PDB_dbAccessionId_A', 'PDB_dbChainId_A', 'PDB_dbResNum_A'])
    df_unique = df_unique.dropna(subset = ["SS_A"])
    df_unique.Alignment_column_A = df_unique.Alignment_column_A.astype(int)
    df_cons_cols = df_unique[df_unique.Alignment_column_A.isin(aln_cols)]
    get_struc_info(df_cons_cols)
    grouped_col_rep = df_cons_cols.groupby(["Alignment_column_A", "SOURCE_ID_A"])
    ss_class_dict = {}
    i = 1
    for k, v in grouped_col_rep:
        col = k[0]
        if col not in ss_class_dict:
            print("Processing consensus column {}".format(i))
            i += 1
            ss_class_dict[col] = []
        ss_class_list = v.SS_A.tolist()
        try:
            ss_class_cons = statistics.mode(ss_class_list)
        except:
            ss_class_cons = ss_class_list[0]
        ss_class_dict[col].append(ss_class_cons)
    return ss_class_dict

def get_ss_df(ss_class_dict):
    n = len(ss_class_dict.keys())
    h = []
    b = []
    e = []
    g = []
    i = []
    t = []
    s = []
    c = []
    for v in ss_class_dict.values():
        h.append(v.count("H"))
        b.append(v.count("B"))
        e.append(v.count("E"))
        g.append(v.count("G"))
        i.append(v.count("I"))
        t.append(v.count("T")) 
        s.append(v.count("S"))
        c.append(v.count(""))
    df_dssp = pd.DataFrame(list(zip(h,b,e,g,i,t,s,c)), columns =  ["a_helix","b_bridge","strand","helix_3_10", "pi_helix","turn","bend","coil"])
    df_dssp["tot"] = df_dssp.a_helix + df_dssp.b_bridge + df_dssp.strand+df_dssp.helix_3_10 + df_dssp.pi_helix + df_dssp.turn+ df_dssp.bend + df_dssp.coil
    df_dssp["p_a_helix"] = df_dssp.a_helix / df_dssp.tot
    df_dssp["p_b_bridge"] = df_dssp.b_bridge / df_dssp.tot
    df_dssp["p_strand"] = df_dssp.strand / df_dssp.tot
    df_dssp["p_3_10_helix"] = df_dssp.helix_3_10 / df_dssp.tot
    df_dssp["p_pi_helix"] = df_dssp.pi_helix / df_dssp.tot
    df_dssp["p_turn"] = df_dssp.turn / df_dssp.tot
    df_dssp["p_bend"] = df_dssp.bend / df_dssp.tot
    df_dssp["p_coil"] = df_dssp.coil / df_dssp.tot
    df_dssp["p_tot"] = df_dssp.p_a_helix + df_dssp.p_b_bridge + df_dssp.p_strand+df_dssp.p_3_10_helix + df_dssp.p_pi_helix +df_dssp.p_turn + df_dssp.p_bend + df_dssp.p_coil
    df_dssp.index = range(1, n + 1)
    return df_dssp

def plot_ss_consensus(df_dssp, palette = ["firebrick", "orangered", "darkorange", "lawngreen","mediumseagreen","aquamarine","teal","cornflowerblue"], bwidth = 1.2, out = None):
    a_helix = list(df_dssp.p_a_helix)
    b_bridge = list(df_dssp.p_b_bridge)
    strand = list(df_dssp.p_strand)
    helix_3_10 = list(df_dssp.p_3_10_helix)
    pi_helix = list(df_dssp.p_pi_helix)
    turn = list(df_dssp.p_turn)
    bend = list(df_dssp.p_bend)
    coil = list(df_dssp.p_coil)
    n_cols = len(df_dssp)
    bottom = [0 for i in range(0, n_cols)]
    #print(bottom)
    #print(len(bottom))
    r = list(np.arange(1,n_cols*1.5,1.5))
    plt.figure(figsize=(180,80))
    plt.rcParams.update({"axes.linewidth": 10})
    plt.bar(r, a_helix,  color=palette[0], bottom = bottom, edgecolor='black', linewidth = 7.5, width=bwidth, label = r'$\alpha$-helix')
    for i in range(0, n_cols):
        bottom[i] += a_helix[i]
    plt.bar(r, helix_3_10, color=palette[1], bottom = bottom,  edgecolor='black', linewidth = 7.5, width=bwidth, label = r'$3_{10}$-helix')
    for i in range(0, n_cols):
        bottom[i] += helix_3_10[i]
    plt.bar(r, pi_helix,  color=palette[2], bottom = bottom,  edgecolor='black', linewidth = 7.5, width=bwidth, label = r'$\pi$-helix')
    for i in range(0, n_cols):
        bottom[i] += pi_helix[i]
    plt.bar(r, b_bridge,  color=palette[3], bottom = bottom, edgecolor='black', linewidth = 7.5, width=bwidth, label = r'$\beta$-bridge')
    for i in range(0, n_cols):
        bottom[i] += b_bridge[i]
    plt.bar(r, strand, color=palette[4], bottom = bottom, edgecolor='black', linewidth = 7.5, width=bwidth, label = r'$\beta$-strand')
    for i in range(0, n_cols):
        bottom[i] += strand[i]
    plt.bar(r, turn,  color=palette[5], bottom = bottom, edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Turn')
    for i in range(0, n_cols):
        bottom[i] += turn[i]
    plt.bar(r, bend,  color=palette[6], bottom = bottom,edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Bend')
    for i in range(0, n_cols):
        bottom[i] += bend[i]
    plt.bar(r, coil, color=palette[7], bottom = bottom, edgecolor='black', linewidth = 7.5, width=bwidth, label = 'Coil')
    
    #plt.title('DSSP Secondary structure assignment', pad = 100, fontsize = 140)
    plt.xlabel('Domain position', labelpad = 100, fontsize = 160)
    plt.ylabel('p', labelpad = 100, fontsize = 160)
    legend = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 120)
    #legend = ax.legend()
    legend.get_frame().set_linewidth(10)
    legend.get_frame().set_edgecolor("black")
    for legobj in legend.legendHandles:
        legobj.set_linewidth(7.5)
    plt.tick_params(axis= 'both' , which = 'major', pad = 60, width = 15, length = 50, labelsize = 140)
    plt.axhline(y=0.5, linewidth = 7.5, linestyle = "--", color = 'black')
    plt.xticks(np.arange(1,n_cols*1.5,1.5),list(range(1, n_cols + 1)))
    plt.yticks(np.arange(0,1.1,0.1))
    plt.xlim(0,(n_cols*1.5)+0.5)
    plt.ylim(0, 1.025)
    if out != None:
        plt.savefig(out)
    plt.show()
    
def get_cx_score_consensus(df, aln_cols):
    """
    """
    df_unique = df.drop_duplicates(['UniProt_dbAccessionId_A', 'PDB_dbAccessionId_A', 'PDB_dbChainId_A', 'PDB_dbResNum_A'])
    df_unique = df_unique.dropna(subset = ["res_cx_min","res_cx_max", "res_cx_mean", "res_cx_median"])
    df_unique.Alignment_column_A = df_unique.Alignment_column_A.astype(int)
    df_cons_cols = df_unique[df_unique.Alignment_column_A.isin(aln_cols)]
    grouped_cols_rep = df_cons_cols.groupby(["Alignment_column_A", "SOURCE_ID_A"])
    df_cons_cols["res_cx_min_mean"] = grouped_cols_rep.res_cx_min.transform('mean')
    df_cons_cols["res_cx_max_mean"] = grouped_cols_rep.res_cx_max.transform('mean')
    df_cons_cols["res_cx_mean_mean"] = grouped_cols_rep.res_cx_mean.transform('mean')
    df_cons_cols["res_cx_median_mean"] = grouped_cols_rep.res_cx_median.transform('mean')
    df_cons_cols_un_res = df_cons_cols.drop_duplicates(["Alignment_column_A", "SOURCE_ID_A"])
    get_struc_info(df_cons_cols)
    grouped_cols = df_cons_cols_un_res.groupby("Alignment_column_A")
    cols_cx_means = pd.concat([
        pd.DataFrame(grouped_cols.res_cx_min_mean.mean()), 
        pd.DataFrame(grouped_cols.res_cx_max_mean.mean()),
        pd.DataFrame(grouped_cols.res_cx_mean_mean.mean()),
        pd.DataFrame(grouped_cols.res_cx_median_mean.mean()),
    ], axis = 1)
    n = len(cols_cx_means)
    cols_cx_means.index = range(1, n + 1)
    return cols_cx_means

def plot_cx_consensus(df, buried_res, out = None):
    df_not_buried = df[~df.index.isin(buried_res)]
    n = len(df)
    plt.figure(figsize=(140,60))
    #plt.rcParams['axes.linewidth'] = 10
    for x in range(1,n+1):
        plt.axvline(x, linewidth = 10, color = 'grey',linestyle = '--')   
    plt.axvspan(3.5, 7.5, facecolor='grey', alpha=0.5)
    #plt.axvspan(7.5, 8.5, facecolor='red', alpha=0.5)
    plt.axvspan(8.5, 10.5, facecolor='grey', alpha=0.5)
    #plt.axvspan(10.5, 15.5, facecolor='red', alpha=0.5)
    #plt.axvspan(15.5, 16.5, facecolor='white', alpha=0.5)
    plt.axvspan(16.5, 18.5, facecolor='grey', alpha=0.5)
    #plt.axvspan(18.5, 20.5, facecolor='orange', alpha=0.5)
    plt.axvspan(20.5, 21.5, facecolor='grey', alpha=0.5)
    #plt.axvspan(21.5, 22.5, facecolor='orange', alpha=0.5)
    #plt.axvspan(22.5, 23.5, facecolor='white', alpha=0.5)
    #plt.axvspan(23.5, 32.5, facecolor='darkgreen', alpha=0.5)
    #plt.axvspan(32.5, 33, facecolor='white', alpha=0.5)
    sns.lineplot(df_not_buried.index, df_not_buried.res_cx_min_mean, color = 'darkblue',linewidth = 15, label = 'Min')
    sns.scatterplot(df.index, df.res_cx_min_mean, color = 'darkblue', s = 3500, edgecolor = 'black', linewidth = 5)
    sns.lineplot(df_not_buried.index, df_not_buried.res_cx_mean_mean, color = 'darkred',linewidth = 15, label = 'Mean')
    sns.scatterplot(df.index, df.res_cx_mean_mean, color = 'darkred', s = 3500, edgecolor = 'black', linewidth = 5)
    sns.lineplot(df_not_buried.index, df_not_buried.res_cx_median_mean, color = 'darkgreen',linewidth = 15, label = 'Median')
    sns.scatterplot(df.index, df.res_cx_median_mean, color = 'darkgreen', s = 3500, edgecolor = 'black', linewidth = 5)
    sns.lineplot(df_not_buried.index, df_not_buried.res_cx_max_mean, color = 'darkorange',linewidth = 15, label = 'Max')
    sns.scatterplot(df.index, df.res_cx_max_mean, color = 'darkorange', s = 3500, edgecolor = 'black', linewidth = 5)
    plt.title('Average CX value per consensus position', pad = 120, fontsize = 160)
    plt.xlabel('Consensus residue position', labelpad = 120, fontsize = 140)
    plt.ylabel('Average CX value', labelpad = 120, fontsize = 140)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 100, title = 'Statistic', title_fontsize = 120)
    plt.tick_params(axis= 'both' , which = 'major', pad = 60, labelsize = 120,width = 5, length = 30)
    plt.xticks(np.arange(1,n+1,1))
    plt.xlim(1,n)
    if out != None:
        plt.savefig(out)
    plt.show()
    
def get_occ_cols(aln_in, fmt_in): #returns a dictionary of the occupoancy of every column in the alignment
    aln = Bio.SeqIO.parse(aln_in, fmt_in)
    occ = {}
    for rec in aln:
        seq = str(rec.seq)
        for i in range(0, len(seq)):
            if i + 1 not in occ:
                occ[i+1] = 0
            if seq[i] != '-':
                occ[i+1] += 1
    return occ

def get_n_seq(seq_file, seq_format = "fasta"):
    """
    """
    seqs = Bio.SeqIO.parse(seq_file,seq_format)
    n = 0
    for seq in seqs:
        n += 1
    return n

def get_aln_len(aln_in, fmt_in):
    aln = Bio.AlignIO.read(aln_in, fmt_in)
    return aln.get_alignment_length()

def get_cons_cols(aln_in, fmt_in, t = 0.5): #returns the columns with a relative occupancy greater than a threshold
    """
    """
    n_seq = get_n_seq(aln_in, fmt_in)
    aln_occ = get_occ_cols(aln_in, fmt_in)
    cons_cols = [col for col, occ in aln_occ.items() if occ/n_seq > t]
    return cons_cols

def get_struc_res_occ(structure_table, col_mask, interaction_mask = None): #returns a dataframe with the occupancy of every alignment column in structure
    if interaction_mask != None:
        structure_table_interaction = structure_table[(structure_table.UniProt_dbAccessionId_A != structure_table.UniProt_dbAccessionId_B)&(structure_table.interaction_type == interaction_mask)]
        structures_mask = structure_table_interaction.PDB_dbAccessionId_A.unique().tolist()
        structure_table = structure_table[structure_table.PDB_dbAccessionId_A.isin(structures_mask)]
    table_a = structure_table[['UniProt_dbAccessionId_A','UniProt_dbResNum_A','Alignment_column_A']]
    table_a = table_a.drop_duplicates(['UniProt_dbAccessionId_A','UniProt_dbResNum_A']).dropna()
    table_a.Alignment_column_A = table_a.Alignment_column_A.astype(int)
    table_a = table_a.rename(columns={'UniProt_dbAccessionId_A':'UniProt_dbAccessionId','UniProt_dbResNum_A':'UniProt_dbResNum','Alignment_column_A':'Alignment_column'})
    
    table_b = structure_table[['UniProt_dbAccessionId_B','UniProt_dbResNum_B','Alignment_column_B']]
    table_b = table_b.drop_duplicates(['UniProt_dbAccessionId_B','UniProt_dbResNum_B']).dropna()
    table_b.Alignment_column_B = table_b.Alignment_column_B.astype(int)
    table_b = table_b.rename(columns={'UniProt_dbAccessionId_B':'UniProt_dbAccessionId','UniProt_dbResNum_B':'UniProt_dbResNum','Alignment_column_B':'Alignment_column'})

    table_ab = pd.concat([table_a,table_b])
    table_ab = table_ab.drop_duplicates(['UniProt_dbAccessionId','UniProt_dbResNum'])
    res_occ = pd.DataFrame(table_ab.Alignment_column.value_counts())
    res_occ.reset_index(inplace = True)
    res_occ = res_occ.rename(columns={'index':'aln_col','Alignment_column':'occ'})
    res_occ = res_occ[res_occ.aln_col.isin(col_mask)]
    res_occ = res_occ.sort_values(by=['aln_col'])
    res_occ.index = range(1,len(col_mask)+1)
    
    return res_occ

def get_ppis(df, aln_len, col_mask):
    df_filt = df.dropna(subset=['UniProt_dbAccessionId_A', 'UniProt_dbAccessionId_B'])
    df_filt = df_filt[(df_filt.UniProt_dbAccessionId_A != df_filt.UniProt_dbAccessionId_B)&(df_filt.interaction_type == 'Protein-Protein')]
    get_struc_info(df_filt)
    df_reps = df_filt.groupby("SOURCE_ID_A")
    cons_norm = pd.DataFrame(0, index = range(1, aln_len+1), columns = ["contacts"])
    for repeat, row in df_reps: # REPEAT
        df_strucs = row.groupby("PDB_dbAccessionId_A")
        rep_cons = pd.DataFrame(0, index = range(1, aln_len+1), columns = ["cons"])
        for struc, row in df_strucs: # STRUCTURE
            df_chains = row.groupby("PDB_dbChainId_A")
            struc_cons = pd.DataFrame(0, index = range(1, aln_len+1), columns = ["cons"])
            for chain, row in df_chains: # CHAIN
                chain_cons = pd.DataFrame(row.groupby("Alignment_column_A")["interaction_type"].count()).reindex(range(1, aln_len+1)).fillna(0)
                struc_cons.cons = struc_cons.cons + chain_cons.interaction_type
            struc_cons.cons = struc_cons.cons.map(int).map(lambda x: 1 if x > 0 else 0)
            rep_cons.cons = rep_cons.cons + struc_cons.cons
        rep_cons.cons = rep_cons.cons.map(int).map(lambda x: 1 if x > 0 else 0)
        cons_norm.contacts = cons_norm.contacts + rep_cons.cons
    cons_norm = cons_norm[cons_norm.index.isin(col_mask)]
    cons_norm.index = range(1,len(col_mask)+1)
    return cons_norm

def get_OR(df):
    df.contacts = df.contacts + 1 #pseudocount
    tot_occ = sum(df.occ)
    tot_cons = sum(df.contacts)
    for i in range(1, len(df)+1):
        i_occ = df.loc[i,"occ"]
        i_cons = df.loc[i,'contacts']
        rest_occ = tot_occ - i_occ
        rest_cons = tot_cons - i_cons
        oddsr, pval = stats.fisher_exact([[i_cons, rest_cons], [i_occ, rest_occ]])
        vals = [i_cons,rest_cons,i_occ, rest_occ]
        se_logor = 1.96*(math.sqrt(sum(list(map((lambda x: 1/x),vals)))))
        logor = math.log(oddsr)
        df.loc[i,'oddsratio'] = oddsr
        df.loc[i,'log_oddsratio'] = logor
        df.loc[i,'pvalue'] = pval
        df.loc[i,'ci_dist'] = se_logor
    return df

def get_surf_dict():
    surf = {}
    core = [4, 5, 6, 7, 9, 10, 17, 18, 21]
    convex = [13, 14, 15, 16, 19, 20, 22, 23, 24]
    concave = [1, 3, 8, 11, 12, 32, 33]
    basal = [25, 26, 27, 28, 29, 30, 31, 2]
    for i in range(1,34):
        if i in core:
            surf[i] = 'Core'
        elif i in convex:
            surf[i] = 'Convex'
        elif i in concave:
            surf[i] = 'Concave'
        elif i in basal:
            surf[i] = 'Basal'
    return surf

def add_surf_class(df):
    surf = get_surf_dict()
    df['surf'] = df.index.map(surf)
    colors_surf = {}
    for k, v in surf.items():
        if v == 'Core':
            colors_surf[k] = 'royalblue'
        elif v == 'Concave':
            colors_surf[k] = 'darkred'
        elif v == 'Convex':
            colors_surf[k] = 'orange'
        elif v == 'Basal':
            colors_surf[k] = 'darkgreen'
    df['color_surf'] = df.index.map(colors_surf)
    return df

def add_miss_class(df):
    UMD = [1,3,8,33]
    UME = [11,12,15,23,24,30,31]
    CMD = [6,9,13, 21,22]
    CME = [4,5]
    for i in range(1,34):
        if i in UMD:
            df.loc[i,'class'] = 'UMD'
            df.loc[i,'color_class'] = 'firebrick'
        elif i in UME:
            df.loc[i,'class'] = 'UME'
            df.loc[i,'color_class'] = 'orange'
        elif i in CMD:
            df.loc[i,'class'] = 'CMD'
            df.loc[i,'color_class'] = 'royalblue'
        elif i in CME:
            df.loc[i,'class'] = 'CME'
            df.loc[i,'color_class'] = 'green'
        else:
            df.loc[i,'class'] = 'None'
            df.loc[i,'color_class'] = 'grey'
    return df


def plot_ppi_enrichment(df, legend_title = None, class_col = None, color_col = None, out = None):
    plt.figure(figsize=(180,80))
    #plt.title("Enrichment in PPIs per residue position", fontsize = 140, pad = 100)
    plt.xlabel("Domain position", fontsize = 160, labelpad = 100)
    plt.ylabel("PPIES", fontsize = 160, labelpad = 100)
    plt.tick_params(axis= 'both' , which = 'major', pad = 60, width = 15, length = 50, labelsize = 140)
    plt.xticks(np.arange(1,len(df)+1,1))
    plt.xlim(0.5,len(df)+0.5)
    plt.axhline(y = 0, color = "black", linewidth = 7.5)
    lines_range = list(np.arange(5.5,len(df),5))
    for point in lines_range:
        plt.axvline(x= point, color = "grey", linestyle = '--', linewidth = 7.5)
    plt.xticks(np.arange(1,len(df)+1,1))
    if class_col != None and color_col != None:
        classes = list(df[class_col].unique())
        for c in classes:
            df_c = df[df[class_col] == c]
            color = df_c[color_col].unique().tolist()[0]
            plt.scatter(df_c.index, df_c.log_oddsratio, c = color, label = c, s = 10000, edgecolor = 'black', linewidth = 10)
            plt.errorbar(df_c.index,df_c.log_oddsratio, yerr=df_c.ci_dist, c = color, linewidth = 10, linestyle="None", capsize = 35.0, capthick = 7.5)
        legend = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 120)#, title = legend_title, title_fontsize = 120)
        legend.get_frame().set_linewidth(10)
        legend.get_frame().set_edgecolor("black")
        for legobj in legend.legendHandles:
            legobj.set_linewidth(7.5)
    else:
        plt.scatter(df.index, df.log_oddsratio, s = 7500, edgecolor = 'black', linewidth = 10)
        plt.errorbar(df.index,df.log_oddsratio, yerr=df.ci_dist, linewidth = 10, linestyle="None", capsize = 35.0, capthick = 7.5)
    if out != None:
        plt.savefig(out)
    plt.show()
    
def get_cons_cols_eq(cons_cols):
    cons_cols_eq = {}
    for i in range(1, len(cons_cols) + 1):
        cons_cols_eq[i] = cons_cols[i-1]
    return cons_cols_eq

def get_color_command(df, cons_cols, pdb_id, class_list, cols_list):
    cons_cols_eq = get_cons_cols_eq(cons_cols)
    df.Alignment_column_A = df.Alignment_column_A.astype(int)
    df_cons = df[df.Alignment_column_A.isin(cons_cols)]
    df_struc = df_cons[df_cons.PDB_dbAccessionId_A == pdb_id].drop_duplicates(["PDB_dbResNum_A"])
    command = ""
    for i, l in enumerate(class_list):
        l_cols = [cons_cols_eq[res] for res in l]
        command+="col {} : ".format(cols_list[i])+", ".join(df_struc[df_struc.Alignment_column_A.isin(l_cols)].PDB_dbResNum_A.unique().tolist())+"\n"
    return command