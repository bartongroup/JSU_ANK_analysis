import pandas as pd
import math
import Bio
from Bio import AlignIO
import pandas as pd
import numpy as np
import seaborn as sns
import scipy as sp
import scipy.stats as stats
import matplotlib.pyplot as plt
import structural_analysis
import matplotlib.patches as mpatches

def get_freqs(col):
    abs_freqs = {
        'A':0, 'R':0, 'N':0, 'D':0,'C':0,'Q':0, 'E':0,'G':0,'H':0, 'I':0,
        'L':0, 'K':0,'M':0, 'F':0, 'P':0,'S':0,'T':0, 'W':0, 'Y':0, 'V':0, '-':0
    }
    for aa in col:
        aa = aa.upper()
        if col.count('-') == len(col):
            abs_freqs['-'] = 1
            return abs_freqs
        if aa not in ['-', 'X', 'B']:
            abs_freqs[aa] += 1
    rel_freqs = {k: v/(len(col)-(col.count('-')+col.count('X')+col.count('B'))) for k, v in abs_freqs.items()}
    return rel_freqs

def get_entropy(freqs):
    S = 0
    for f in freqs.values():
        if f != 0:
            S += f*math.log2(f)
    return -S

def get_shenkin(col):
    S = get_entropy(get_freqs(col))
    return (2**S)*6

def in_columns(aln_in, infmt):
    aln = Bio.AlignIO.read(aln_in, infmt)
    n_cols = len(aln[0])
    cols = {}
    for col in range(1,n_cols+1):
        cols[col] = []
    for row in aln:
        seq = str(row.seq)
        for i in range(0,len(seq)):
            cols[i+1].append(seq[i])
    return cols

def get_stats(col):
    n_seqs = len(col)
    gaps = col.count('-')
    occ = n_seqs - gaps
    occ_pct = occ/n_seqs
    gaps_pct = 1 - occ_pct
    return occ, gaps, occ_pct, gaps_pct

def calculate_shenkin(aln_in, infmt, t = 0.5):
    cols = in_columns(aln_in, infmt)
    scores = []
    occ = []
    gaps = []
    occ_pct = []
    gaps_pct = []
    for k, v in cols.items():
        scores.append(get_shenkin(v))
        stats = (get_stats(v))
        occ.append(stats[0])
        gaps.append(stats[1])
        occ_pct.append(stats[2])
        gaps_pct.append(stats[3])
    df = pd.DataFrame(list(zip(list(range(1,len(scores)+1)),scores, occ,gaps, occ_pct, gaps_pct)), columns = ['col','shenkin','occ','gaps','occ_pct','gaps_pct'])
    #df = df[df.occ_pct >= t]
    max_shenkin = max(list(df.shenkin))
    #df.index = range(1, len(df)+1)
    return df

def plot_cons_occ(df, qs, colors = ["darkorchid", "darkorange"], cons_score_col = 'shenkin', cons_score_label = "Shenkin score", out = None):
    plt.figure(figsize=(140,60))
    plt.rcParams.update({"axes.linewidth": 5})
    ax = sns.barplot(df.index, df[cons_score_col], edgecolor = 'black',linewidth = 5, color = colors[0], label = "Shenkin")
    ax.set_title('Shenkin divergence score', pad =100, fontsize = 140)
    ax.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 100)
    ax.set_xlabel("Consensus residue position", labelpad = 100, fontsize = 120)
    ax.set_ylabel(cons_score_label, labelpad = 100, fontsize = 120)
    ax.plot([], [], color='orange', lw=15, label='% Occupancy')
    plt.axhline(y = qs[0], color = "black", linewidth = 5, linestyle = '--')
    plt.axhline(y = qs[1], color = "black", linewidth = 5, linestyle = '--')
    ax2 = ax.twinx()
    sns.lineplot(df.index-1, df.occ_pct.map(lambda x: x*100), color = colors[1],linewidth = 15)
    sns.scatterplot(df.index-1, df.occ_pct.map(lambda x: x*100), color = colors[1], s = 3500, edgecolor = 'black', linewidth = 5)
    ax2.set_ylabel("% Occupancy", fontsize = 100, labelpad = 100)
    ax2.set_yticks(np.arange(0,125,25))
    ax2.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 80)
    plt.xlim(-0.7,len(df)-0.3)
    ax.legend(loc='best', fontsize = 80)
    if out != None:
        plt.savefig(out)
    plt.show()
    
def plot_conservation(df, colors = ["darkorchid", "darkorange"], cons_score_col = 'shenkin', cons_score_label = "Shenkin score", shenkin_threshold = 60, out = None):
    plt.figure(figsize=(140,60))
    plt.rcParams.update({"axes.linewidth": 5})
    palette = {}
    for i in df.index:
        if df.loc[i, "shenkin"] < 60:
            palette[i] = "steelblue"
        else:
            palette[i] = "orangered"
    ax = sns.barplot(df.index, df[cons_score_col], edgecolor = 'black',linewidth = 5, palette = palette, label = "Shenkin")
    ax.set_title('Shenkin divergence score', pad =100, fontsize = 140)
    ax.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 100)
    ax.set_xlabel("Consensus residue position", labelpad = 100, fontsize = 120)
    ax.set_ylabel(cons_score_label, labelpad = 100, fontsize = 120)
    #ax.plot([], [], color='orange', lw=15, label='% Occupancy')
    plt.axhline(y = shenkin_threshold, color = "black", linewidth = 5, linestyle = '--')
    #plt.axhline(y = qs[1], color = "black", linewidth = 5, linestyle = '--')
    #ax2 = ax.twinx()
    #sns.lineplot(df.index-1, df.occ_pct.map(lambda x: x*100), color = colors[1],linewidth = 15)
    #sns.scatterplot(df.index-1, df.occ_pct.map(lambda x: x*100), color = colors[1], s = 3500, edgecolor = 'black', linewidth = 5)
    #ax2.set_ylabel("% Occupancy", fontsize = 100, labelpad = 100)
    #ax2.set_yticks(np.arange(0,125,25))
    #ax2.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 80)
    plt.xlim(-0.7,len(df)-0.3)
    if out != None:
        plt.savefig(out)
    plt.show()

def plot_conservation_mod(df, colors = ["steelblue", "orange", "tomato", "orchid"], cons_score_col = 'shenkin', cons_score_label = "Shenkin score", qs = [60,], out = None):
    plt.figure(figsize=(180,80))
    plt.rcParams.update({"axes.linewidth": 10})
    palette = {}
    for i in df.index:
        if df.loc[i, cons_score_col] <= qs[0]:
            palette[i] = colors[0]
        elif df.loc[i, cons_score_col] > qs[0] and df.loc[i, cons_score_col] <= qs[1]:
            palette[i] = colors[1]
        elif df.loc[i, cons_score_col] > qs[1] and df.loc[i, cons_score_col] <= qs[2]:
            palette[i] = colors[2]
        else:
            palette[i] = colors[3]
            #r'$3_{10}$-helix'
    legend_dict = { r'$N_{s} \leq 25$' : colors[0], r'$25 < N_{s} \leq 50$' : colors[1], r'$50 < N_{s} \leq 75$' : colors[2], r'$N_{s} > 75$' : colors[3]}
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(facecolor=legend_dict[key], label=key, edgecolor = "black", linewidth = 3)
        patchList.append(data_key)
    ax = sns.barplot(df.index, df[cons_score_col], edgecolor = 'black',linewidth = 7.5, palette = palette, label = "Shenkin")
    #ax.set_title('Shenkin divergence score', pad =100, fontsize = 140)
    ax.tick_params(axis= 'both' , which = 'major', pad = 60, width = 15, length = 50, labelsize = 140)
    ax.set_xlabel("Domain position", labelpad = 100, fontsize = 160)
    ax.set_ylabel(cons_score_label, labelpad = 100, fontsize = 160)
    #ax.plot([], [], color='orange', lw=15, label='% Occupancy')
    plt.yticks(np.arange(0, 125, 25))
    for i, q in enumerate(qs):
        plt.axhline(y = q, color = "black", linewidth = 7.5, linestyle = '--') #colors[i]
    #plt.axhline(y = qs[1], color = "black", linewidth = 5, linestyle = '--')
    #ax2 = ax.twinx()
    #sns.lineplot(df.index-1, df.occ_pct.map(lambda x: x*100), color = colors[1],linewidth = 15)
    #sns.scatterplot(df.index-1, df.occ_pct.map(lambda x: x*100), color = colors[1], s = 3500, edgecolor = 'black', linewidth = 5)
    #ax2.set_ylabel("% Occupancy", fontsize = 100, labelpad = 100)
    #ax2.set_yticks(np.arange(0,125,25))
    #ax2.tick_params(axis= 'both' , which = 'major', pad = 60, width = 5, length = 30, labelsize = 80)
    plt.xlim(-0.7,len(df)-0.3)
    legend = plt.legend(handles=patchList, loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 120)
    #legend = ax.legend()
    legend.get_frame().set_linewidth(10)
    legend.get_frame().set_edgecolor("black")
    for legobj in legend.legendHandles:
        legobj.set_linewidth(7.5)
    if out != None:
        plt.savefig(out)
    plt.show()
    
def format_variant_table(df, col_mask, vep_mask = ["missense_variant"], tab_format = "gnomad"):
    df_filt = df.copy(deep = True)
    df_filt.reset_index(inplace=True)
    if tab_format == "gnomad":
        df_filt.columns = [' '.join(col).strip() for col in df_filt.columns.tolist()]
    df_filt.columns = [col.lower().replace(" ", "_") for col in df_filt.columns.tolist()]
    df_filt = df_filt[df_filt.source_id.str.contains("HUMAN")]
    df_filt = df_filt.dropna(subset=["vep_consequence"])
    #print(len(df_filt))
    #print(len(df_filt.source_id.unique().tolist()))
    #print(df_filt.vep_consequence.value_counts())
    #print(df_filt.vep_consequence.unique().tolist())
    df_filt = df_filt[df_filt.vep_consequence.isin(vep_mask) ]
    df_filt = df_filt[df_filt.alignment_column.isin(col_mask)]
    return df_filt

def generate_subset_aln(aln_in, aln_fmt, df, aln_out = None):
    seqs_ids = df.source_id.unique().tolist()
    aln = Bio.SeqIO.parse(aln_in, aln_fmt)
    variant_seqs = [rec for rec in aln if rec.id in seqs_ids]
    if aln_out == None:
        pref, fmt = aln_in.split(".")
        aln_out =  pref + "_variant_seqs."+fmt
    Bio.SeqIO.write(variant_seqs, aln_out, aln_fmt)
    return aln_out

def get_OR(df, variant_col = "variants"):
    tot_occ = sum(df.occ)
    tot_vars = sum(df[variant_col])
    idx = df.index.tolist()
    for i in idx:
        i_occ = df.loc[i,"occ"]
        i_vars = df.loc[i,variant_col]
        rest_occ = tot_occ - i_occ
        rest_vars = tot_vars - i_vars
        oddsr, pval = stats.fisher_exact([[i_vars, rest_vars], [i_occ, rest_occ]])
        #print(i_vars, rest_vars, i_occ, rest_occ)
        #print(i, oddsr, pval)
        vals = [i_vars,rest_vars,i_occ, rest_occ]
        se_logor = 1.96*(math.sqrt(sum(list(map((lambda x: 1/x),vals)))))
        logor = math.log(oddsr)
        #print(i, oddsr, logor, pval, se_logor)
        df.loc[i,'oddsratio'] = oddsr
        df.loc[i,'log_oddsratio'] = logor
        df.loc[i,'pvalue'] = pval
        df.loc[i,'ci_dist'] = se_logor
    return df


def add_miss_class(df, cons_col = "shenkin", thresholds = [30, 70], colours = ["royalblue", "green", "grey", "firebrick", "orange"]):
    for i in df.index:
        if df.loc[i, cons_col] <= thresholds[0] and df.loc[i, "log_oddsratio"] < 0:
            df.loc[i,"miss_class"] = "CMD"
        elif df.loc[i, cons_col] <= thresholds[0] and df.loc[i, "log_oddsratio"] > 0:
            df.loc[i,"miss_class"] = "CME"
        elif df.loc[i, cons_col] >= thresholds[1] and df.loc[i, "log_oddsratio"] < 0:
            df.loc[i,"miss_class"] = "UMD"
        elif df.loc[i, cons_col] >= thresholds[1] and df.loc[i, "log_oddsratio"] > 0:
            df.loc[i,"miss_class"] = "UME"
        else:
            df.loc[i,"miss_class"] = "None"
           
    coloring = {
        "CMD": colours[0],
        "CME": colours[1],
        "UMD": colours[3],
        "UME": colours[4],
        "None": colours[2]
    }
    df["miss_color"] =  df.miss_class.map(coloring) 
    return df

def get_missense_df(aln_in, aln_fmt, variants_df, cons_cols, shenkin_aln, aln_out = None, t = 0.5, get_or = True):
    variants_aln = generate_subset_aln(aln_in, aln_fmt, variants_df, aln_out)
    variants_aln_info = calculate_shenkin(variants_aln, aln_fmt, t)
    variants_aln_info = variants_aln_info[variants_aln_info.occ_pct >= t]
    cons_cols = structural_analysis.get_cons_cols(variants_aln, aln_fmt, t = t)
    vars_df = pd.DataFrame(variants_df.alignment_column.value_counts().reindex(cons_cols, fill_value=0).sort_index()).reset_index()
    vars_df.index = range(1, len(cons_cols)+1)
    vars_df.columns = ["col", "variants"]
    merged = pd.merge(variants_aln_info, vars_df, on = "col", how = 'left')
    merged.index = range(1, len(vars_df)+1)
    merged.variants = merged.variants + 1
    merged["shenkin"] = shenkin_aln["shenkin"]
    merged["norm_shenkin_rel"] = shenkin_aln["norm_shenkin_rel"] # REPLACES ONLY HUMAN DIVERGENCE SCORES BY ALL SPECIES, FULL ALIGNMENT PROVIDES A BETTER MEASURE FOR THE OVERALL DIVERGENCE
    merged["norm_shenkin_abs"] = shenkin_aln["norm_shenkin_abs"]
    if get_or == True:
        merged_or = get_OR(merged)
        return merged_or
    else:
        return merged

def plot_variants_enrichment(df, legend_title = None, class_col = None, color_col = None, out = None):
    plt.figure(figsize=(140,60))
    plt.title("Splicing region variants per residue position", fontsize = 140, pad = 100)
    plt.xlabel("Consensus residue position", fontsize = 120, labelpad = 100)
    plt.ylabel("Splicing region variants enrichment score", fontsize = 120, labelpad = 100)
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
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 80, title = legend_title, title_fontsize = 120)
    else:
        plt.scatter(df.index, df.log_oddsratio, s = 5000, edgecolor = 'black', linewidth = 7.5)
        plt.errorbar(df.index,df.log_oddsratio, yerr=df.ci_dist, linewidth = 7.5, linestyle="None", capsize = 35.0, capthick = 7.5)
    if out != None:
        plt.savefig(out)
    plt.show()
    
def plot_shenkin_logOR(df, strat_column, color_column, cons_col = "shenkin", thresholds = [30,70], out = None, markers = ["h", ".", "s", "*", "^"], colors = ["darkred", "tan", "green", "darkblue", "yellow"], label = True, error = True):
    plt.figure(figsize=(140,100))
    classes_df = sorted(df[strat_column].unique().tolist())
    #print(classes_df)
    i = 0
    for c in classes_df:
        #print(c)
        df_c = df[df[strat_column] == c]
        color = df_c[color_column].unique().tolist()[0]
        plt.scatter(df_c[cons_col], df_c.log_oddsratio,  marker = markers[i], c = colors[i], label = c, s = 15000, edgecolor = 'black', linewidth = 15)
        if error == True:
            plt.errorbar(df_c[cons_col],df_c.log_oddsratio, yerr=df_c.ci_dist, c = colors[i], linewidth = 15, linestyle="None", capsize = 35.0, capthick = 15)
        i += 1
    if label == True:
        for x, y, z in zip(df[cons_col], df.log_oddsratio, df.index):
            label = "{:d}".format(z)
            plt.annotate(label, # this is the text
                            (x,y), # this is the point to label
                            textcoords="offset points", # how to position the text
                            xytext=(-75,75), # distance from text to points (x,y)
                            ha='right',
                            fontsize = 140) # horizontal alignment can be left, right or center
    #plt.title("Missense enrichment score relative to sequence divergence", fontsize = 200, pad = 140)
    plt.xlabel("Normalised Shenkin divergence score", fontsize = 180, labelpad = 120)
    plt.ylabel("MES", fontsize = 180, labelpad = 120)
    plt.tick_params(axis= 'both' , which = 'major', pad = 60, width = 15, length = 50, labelsize = 160)
    #plt.yticks(np.arange(-0.6,0.3,0.1))
    plt.xticks(np.arange(0,105,5))
    plt.xlim(-5,105)
    plt.axhline(color = "black", linewidth = 10, linestyle = '--')
    plt.axvline(x= thresholds[0], color = "black", linewidth = 10, linestyle = '--')
    plt.axvline(x= thresholds[1], color = "black", linewidth = 10, linestyle = '--')
    legend = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 140)
    legend.get_frame().set_linewidth(10)
    legend.get_frame().set_edgecolor("black")
    for legobj in legend.legendHandles:
        legobj.set_linewidth(7.5)
    if out != None:
        plt.savefig(out)
    plt.show()
    