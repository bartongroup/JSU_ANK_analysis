import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import structural_analysis
from scipy import stats
import math
import variant_analysis

def get_cons_cols_eq(cons_cols):
    cons_cols_eq = {}
    for i, col in enumerate(cons_cols):
        cons_cols_eq[col] = i+1
    return cons_cols_eq

def format_df(df, cons_cols):
	cons_cols_eq = get_cons_cols_eq(cons_cols)
	df2 = df.copy(deep = True)
	df2_filt = df2.dropna(subset = ["Alignment_column_A", "Alignment_column_B"])
	df2_filt.Alignment_column_A = df2_filt.Alignment_column_A.astype(int)
	df2_filt.Alignment_column_B = df2_filt.Alignment_column_B.astype(int)
	df2_filt.UniProt_dbResNum_A = df2_filt.UniProt_dbResNum_A.astype(int)
	df2_filt.UniProt_dbResNum_B = df2_filt.UniProt_dbResNum_B.astype(int)
	df2_cons = df2_filt[(df2_filt.Alignment_column_A.isin(cons_cols)) & (df2_filt.Alignment_column_B.isin(cons_cols))]
	df2_cons["Alignment_column_cons_A"] = df2_cons.Alignment_column_A.map(cons_cols_eq)
	df2_cons["Alignment_column_cons_B"] = df2_cons.Alignment_column_B.map(cons_cols_eq)
	return df2_cons

def get_res_in_rep(df, rep_id, cons_eq):
    df_rep = df[(df.SOURCE_ID_A == rep_id)&(df.SOURCE_ID_B == rep_id)]
    res_a = df_rep.Alignment_column_A.unique().tolist()
    res_b = df_rep.Alignment_column_B.unique().tolist()
    res = sorted(list(set(res_a + res_b)))
    cons_res = [cons_eq[r] for r in res]
    return cons_res

def get_intra_cons_occ(df, cons_eq):
    structural_analysis.get_struc_info(df)
    reps = df.SOURCE_ID_A.unique().tolist()
    res_occ = {}
    for rep in reps:
        rep_res = get_res_in_rep(df, rep, cons_eq)
        for i in rep_res:
            for j in rep_res:
                idx = tuple([i, j])
                if idx not in res_occ:
                    res_occ[idx] = 0
                res_occ[idx] += 1
    return res_occ

def get_inter_cons_occ(df, cons_eq):
    structural_analysis.get_struc_info(df)
    df_filt = df.sort_values(by = ["UniProt_dbAccessionId_A", "UniProt_dbResNum_A"]) #sorts them by repeat order
    prots = df_filt.groupby("UniProt_dbAccessionId_A")
    res_occ = {}
    pairs = 0
    for prot, rows in prots:
        reps = rows.SOURCE_ID_A.unique().tolist()
        for i, rep in enumerate(reps):
            try:
                rep_res = get_res_in_rep(rows, rep, cons_eq)
                next_rep_res = get_res_in_rep(rows, reps[i + 1], cons_eq)
                pairs += 1
                for i in rep_res:
                    for j in next_rep_res:
                        idx = tuple([i, j])
                        if idx not in res_occ:
                            res_occ[idx] = 0
                        res_occ[idx] += 1
            except:
                continue
    print("{} pairs of ARs were used".format(pairs))
    return res_occ

def get_df_cols_eq(df):
    cols = df.columns.tolist()
    un_cols = list(set([col[:-2] for col in cols]))
    cols_eq = {}
    for col in cols:
        if col[-2:] == "_A":
            cols_eq[col] = col[:-2] + "_B"
        elif col[-2:] == "_B":
            cols_eq[col] = col[:-2] + "_A"
        else:
            cols_eq[col] = col
    return cols_eq

def fix_direction(df, rep1, rep2, cols_eq):
    df_fixed = df.copy(deep = True)
    df_r = df_fixed[df_fixed.SOURCE_ID_A == rep1]
    df_w = df_fixed[df_fixed.SOURCE_ID_A == rep2]
    df_r = df_r.reindex(sorted(df_r.columns), axis = 1)
    df_w = df_w.rename(columns = cols_eq)
    df_w = df_w.reindex(sorted(df_w.columns), axis = 1)
    df_fixed = pd.concat([df_r, df_w])
    return df_fixed

def symmetric_contact_matrix(df):
    df2 = pd.DataFrame.copy(df, deep = True)
    rows = list(df.index)
    cols = list(df.columns)
    for i in rows:
        for j in cols[i-1:]:
            df2.loc[j,i] = df2.loc[i,j]
    return df2

def get_inter_cons(df, cons_cols_eq, df_cols_eq, int_mask = None):
    df2 = df.copy(deep = True)
    if int_mask != None:
        df2 = df2[df2.Int_Types.str.contains(int_mask)]
    df_filt = df2.sort_values(by = ["UniProt_dbAccessionId_A", "UniProt_dbResNum_A"]) #sorts them by repeat order
    prots = df_filt.groupby("UniProt_dbAccessionId_A")
    real_cols = list(cons_cols_eq.keys())
    cons_cols = list(cons_cols_eq.values())
    contact_matrix_inter = pd.DataFrame(np.nan, index = cons_cols, columns = cons_cols).fillna(0)
    for prot, rows in prots:
        reps = rows.SOURCE_ID_A.unique().tolist()
        for i, rep in enumerate(reps):
            try:
                rep1 = rep
                rep2 = reps[i+1]
                reps_pair = [rep1, rep2]
                pair_df = rows[(rows.SOURCE_ID_A != rows.SOURCE_ID_B)&((rows.SOURCE_ID_A.isin(reps_pair)) & (rows.SOURCE_ID_B.isin(reps_pair)))]
                pair_df_fixed = fix_direction(pair_df, rep1, rep2, df_cols_eq)
                pair_df_fixed = pair_df_fixed.drop_duplicates(subset = ["Alignment_column_cons_A", "Alignment_column_cons_B"])
                repeat_crosstab = pd.crosstab(pair_df_fixed.Alignment_column_cons_A, pair_df_fixed.Alignment_column_cons_B)
                repeat_crosstab = repeat_crosstab.reindex(cons_cols).fillna(0) # changes index (rows) so it includes all positions (1-30)
                repeat_crosstab = repeat_crosstab.reindex(columns = cons_cols).fillna(0) # changes column index so it includes all positions (1-30)
                repeat_crosstab[repeat_crosstab >= 1] = 1 # flattening all interatomic interactions and overrepresented interactions in multiple structure to value of one
                contact_matrix_inter += repeat_crosstab # adding each crosstab for each repeat to the contact matrix dataframe to study all interaction
            except:
                continue
    return contact_matrix_inter

def get_intra_cons(df, cons_cols_eq, t = 0, int_mask = None):
    real_cols = list(cons_cols_eq.keys())
    cons_cols = list(cons_cols_eq.values())
    df2 = df.copy(deep = True)
    if int_mask != None:
        df2 = df2[df2.Int_Types.str.contains(int_mask)]
    contact_matrix_intra = pd.DataFrame(np.nan, index = cons_cols, columns = cons_cols).fillna(0)
    if t == 0:
        df_filt = df2[df2.SOURCE_ID_A == df2.SOURCE_ID_B].groupby("SOURCE_ID_A")
    else:
        df_filt = df2[(df.SOURCE_ID_A == df2.SOURCE_ID_B) & (abs(df2.Alignment_column_cons_A - df2.Alignment_column_cons_B) > t)].groupby("SOURCE_ID_A")
    for repeat, row in df_filt:
        repeat_crosstab = pd.crosstab(row.Alignment_column_cons_A,
                                      row.Alignment_column_cons_B)    # creates a cross table of how many contacts there are between any two given residues
        repeat_crosstab = repeat_crosstab.reindex(cons_cols).fillna(0)         # changes index (rows) so it includes all positions (1-33)
        repeat_crosstab = repeat_crosstab.reindex(columns = cons_cols).fillna(0) # changes column index so it includes all positions (1-33)
        repeat_crosstab[repeat_crosstab >= 1] = 1 # flattening all interatomic interactions and overrepresented interactions in multiple structure to value of one
        contact_matrix_intra += repeat_crosstab  # adding each crosstab for each repeat to the contact matrix dataframe to study all interactions
    contact_matrix_intra_symm = symmetric_contact_matrix(contact_matrix_intra)
    return contact_matrix_intra_symm

def normalize_contacts(df, res_occ_dict):
    df2 = df.copy(deep = True)
    for k, v in res_occ_dict.items():
        df2.loc[k[0],k[1]] = df2.loc[k[0],k[1]]/v
    return df2

def plot_inter_cons(df, cmap = 'Greens', out = None):
    fig, ax = plt.subplots(1, 1, figsize = (200, 160))
    ax = sns.heatmap(df, cmap = cmap,  square = True, linewidths = 5, linecolor = 'black', cbar_kws = dict(ticks = np.arange(0,1.1,0.1)))
    ax.xaxis.tick_top()
    ax.yaxis.tick_left()
    ax.tick_params(axis = 'both' , labelrotation = 'auto', which = 'major', pad = 60, width = 5, length = 30, labelsize = 180)
    ax.figure.axes[-1].yaxis.label.set_size(180)
    ax.set_ylim(len(df),0)
    cbar = ax.collections[0].colorbar
    cbar.set_label('p', labelpad = 190)
    cbar.ax.tick_params(pad = 60, width = 5, length = 30, labelsize = 180)
    plt.xlabel(r'Consensus residue position $(AR_{n+1})$', fontsize = 180, labelpad = 160)
    plt.ylabel(r'Consensus residue position $(AR_n)$', fontsize = 180, labelpad = 160)
    ax.xaxis.set_label_position('top')
    plt.title("Inter - ANK repeat contact map", fontsize = 200, pad = 160)
    if out != None:
        plt.savefig(out)
    plt.show()

def plot_intra_cons(df, cmap = 'Greens', out = None):
    fig, ax = plt.subplots(1, 1, figsize=(200, 160))
    ax = sns.heatmap(df, cmap = cmap, square = True, linewidths = 5, linecolor = 'black', cbar_kws = dict(ticks = np.arange(0,1.1,0.1)))
    ax.xaxis.tick_top()
    ax.yaxis.tick_left()
    ax.tick_params(axis= 'both' , labelrotation = 'auto', which = 'major', pad = 60, width = 5, length = 30, labelsize = 180)
    ax.figure.axes[-1].yaxis.label.set_size(180)
    ax.set_ylim(len(df),0)
    cbar = ax.collections[0].colorbar
    cbar.set_label('p', labelpad=190)
    cbar.ax.tick_params(pad = 60, width = 5, length = 30, labelsize = 180)
    plt.xlabel("Consensus residue position", fontsize = 180, labelpad = 160)
    plt.ylabel("Consensus residue position", fontsize = 180, labelpad = 160)
    ax.xaxis.set_label_position('top')
    plt.title("Intra - ANK repeat contact map", fontsize = 200, pad = 160) #450
    if out != None:
        plt.savefig(out)
    plt.show()

def get_tot_occ_res(occ_dict, n_cols, contact_map = "intra", t = 0):
    res_occ = {}
    idx = range(1, n_cols + 1)
    for i in idx:
        if i not in res_occ:
            res_occ[i] = 0
        for j in idx:
            if contact_map == "intra":
                if abs(j - i) > t:
                    res_occ[i] += occ_dict[(i,j)]
            elif contact_map == "inter":
                if i == j:
                    res_occ[i] += occ_dict[(i,j)]
                else:
                    res_occ[i] += (occ_dict[(i,j)] + occ_dict[(j,i)])         
    return res_occ

def get_tot_cons_res(cons_df, contact_map = "intra", t = 0):
    res_cons = {}
    idx = cons_df.index.tolist()
    for i in idx:
        if i not in res_cons:
            res_cons[i] = 0
        for j in idx:
            if contact_map == "intra":
                if abs(j - i) > t:
                    res_cons[i] += cons_df.loc[i,j]
            elif contact_map == "inter":
                if i == j:
                    res_cons[i] += cons_df.loc[i,j]
                else:
                    res_cons[i] += (cons_df.loc[i,j] + cons_df.loc[j,i])         
    return res_cons

def get_OR_from_cons(res_cons, res_occ):
    cons = list(res_cons.values())
    occ = list(res_occ.values())
    df = pd.DataFrame(list(zip(occ, cons)), columns =['occ', 'contacts'])
    df.index = range(1,len(cons)+1)
    df = get_OR(df)
    return df

def get_OR(df):
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

def plot_cons_enrichment(df, contact_map = "Intra", legend_title = None, class_col = None, color_col = None, out = None):
    plt.figure(figsize=(140,60))
    plt.rcParams.update({"axes.linewidth": 5})
    plt.title("Enrichment in {}-ANK contacts per residue position".format(contact_map.lower()), fontsize = 140, pad = 100)
    plt.xlabel("Consensus residue position", fontsize = 120, labelpad = 100)
    plt.ylabel("{}-ANK contacts enrichment score".format(contact_map), fontsize = 120, labelpad = 100)
    plt.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 100)
    plt.xticks(np.arange(1,len(df)+1,1))
    plt.xlim(0.5,len(df)+0.5)
    plt.axhline(y = 0, color = "black", linewidth = 5)
    lines_range = list(np.arange(5.5,len(df),5))
    for point in lines_range:
        plt.axvline(x= point, color = "grey", linestyle = '--', linewidth = 5)
    plt.xticks(np.arange(1,len(df)+1,1))
    if class_col != None and color_col != None:
        classes = list(df[class_col].unique())
        for c in classes:
            df_c = df[df[class_col] == c]
            color = df_c[color_col].unique().tolist()[0]
            plt.scatter(df_c.index, df_c.log_oddsratio, c = color, label = c, s = 5000, edgecolor = 'black', linewidth = 10)
            plt.errorbar(df_c.index,df_c.log_oddsratio, yerr=df_c.ci_dist, c = color, linewidth = 7.5, linestyle="None", capsize = 35.0, capthick = 7.5)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 80)#, title = legend_title, title_fontsize = 120)
    else:
        plt.scatter(df.index, df.log_oddsratio, s = 5000, edgecolor = 'black', linewidth = 7.5)
        plt.errorbar(df.index,df.log_oddsratio, yerr=df.ci_dist, linewidth = 7.5, linestyle="None", capsize = 35.0, capthick = 7.5)
    if out != None:
        plt.savefig(out)
    plt.show()
