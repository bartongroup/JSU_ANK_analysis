import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
import csv
import Bio
import Bio.SeqIO
import pandas as pd
import gzip
import os
import urllib
import requests
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict

def get_tax_ids_dict(wd, download = False):
    """
    :param wd: working directory where the files are going to be stored
    :param download: whether to download url content
    :type wd: str
    :type download: bool
    :returns: taxon mnemonic-taxon ID relationship
    :rtype: dict
    """
    taxids_cropped_path = os.path.join(wd, "tax_ids_rev_cropped.tab.gz")
    if download == True:
        taxids_path = os.path.join(wd, "tax_ids_rev.tab.gz")
        baseURL = "http://www.uniprot.org/taxonomy/"
        payload = {"query": "reviewed:yes", "format": "tab", "compress": "yes"}
        result = requests.get(baseURL, params = payload)
        if result.ok:
            print("Server response was OK, code {}".format(result.status_code))
            open(taxids_path, "wb").write(result.content)
            taxids = pd.read_table(taxids_path)
            taxids = taxids[["Taxon","Mnemonic"]]
            taxids.to_csv(taxids_cropped_path, index = False, compression = "gzip")
        else:
            print("Something went wrong: ", result.status_code)
            return
    taxids = pd.read_csv(taxids_cropped_path)
    taxids_dict = dict(zip(taxids.Taxon,taxids.Mnemonic))
    return taxids_dict

def get_accs_tax_ids_dict(wd, download = False):
    """
    :param wd: working directory where the files are going to be stored
    :param download: whether to download url content
    :type wd: str
    :type download: bool
    :returns: UniProt accession ID-taxon ID relationship
    :rtype: dict
    """
    accsids_path = os.path.join(wd, "accs_ids_rev.tab.gz")
    if download == True:
        baseURL = "http://www.uniprot.org/uniprot/"
        payload = {"query": "reviewed:yes", "format": "tab", "compress": "yes", "columns": "id,organism-id"}
        result = requests.get(baseURL, params = payload)
        if result.ok:
            print("Server response was OK, code {}".format(result.status_code))
            open(accsids_path, "wb").write(result.content)
            accsids = pd.read_table(accsids_path)
        else:
            print("Something went wrong: ", result.status_code)
            return
    else:
        accsids = pd.read_table(accsids_path)
    accsids_dict = dict(zip(accsids.Entry, accsids["Organism ID"]))
    return accsids_dict

def retrieve_uniprot(wd, key, download = False):
    """
    Retrieves annotations from UniProt records file in xml format.

    :param wd: working directory where the files are going to be stored
    :param key: keyword for repeat family
    :param download: whether to download url content
    :param species: species of interest
    :type wd: str
    :type key: str
    :type download: bool
    :returns: dataframe containing all UniProt annotations for the keyword
    :rtype: pandas.DataFrame

    """
    rec_file = os.path.join(wd,"{}_rev_allsp.xml.gz".format(key))
    if download == True:
        baseURL = "http://www.uniprot.org/uniprot/"
        payload = {"query": "keyword:{} AND reviewed:yes".format(key), "format": "xml", "compress": "yes"}
        result = requests.get(baseURL, params=payload)
        if result.ok:
            print("Server response was OK, code {}".format(result.status_code))
            open(rec_file, "wb").write(result.content)
        else:
            print("Something went wrong: ", result.status_code)
            return
    handle = gzip.open(rec_file, "rt")
    records = Bio.SeqIO.parse(handle, "uniprot-xml")
    uniprot_dict = {
    "accession": [], "source": [], "start": [], "end": [], "length": []
    }
    for rec in records:
        for feat in rec.features:
            if feat.type == "repeat" and feat.qualifiers["description"].split()[0] == key:
                uniprot_dict["accession"].append(rec.id)
                uniprot_dict["source"].append("Uniprot")
                uniprot_dict["length"].append(len(rec.seq))
                uniprot_dict["start"].append(feat.location.start + 1)
                uniprot_dict["end"].append(feat.location.end)
    uniprot_df = pd.DataFrame(uniprot_dict)
    uniprot_df["repeat_id"] = uniprot_df.accession + "/" + uniprot_df.start.apply(str) + "-" + uniprot_df.end.apply(str)
    return uniprot_df

def output_list(sig):
    """
    Returns accession numbers of proteins having annotations for a given signature.

    :param sig: signature accession of interest
    :type sig: str
    :returns: accession of proteins with signature annotations
    :rtype: list
    """
    if "PF" in sig:
        db = "pfam"
    elif "SM" in sig:
        db = "smart"
    elif "PS" in sig:
        db = "profile"
    elif "PR" in sig:
        db = "prints"
        
    BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/{}/{}/".format(db, sig)
    accs = []
    context = ssl._create_unverified_context() #disable SSL verification to avoid config issues
    next = BASE_URL
    last_page = False
    while next:
        try:
            req = request.Request(next, headers = {"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            if res.status == 408:  # If the API times out due a long running query
                sleep(61) # wait just over a minute
                continue # then continue this loop with the same URL
            elif res.status == 204:
                break #no data so leave loop
            payload = json.loads(res.read().decode())
            next = payload["next"]
            if not next:
                last_page = True
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                raise e

        for i, item in enumerate(payload["results"]):
            accs.append(item["metadata"]["accession"])
        
        if next: # Don't overload the server, give it time before asking for more
            sleep(1)
    return accs

def get_sig_accs_dict(sigs):
    """
    Returns dictionary with all accessions (all species) for a given signature
    
    :param sigs: signatures of interest
    :type sigs: list
    
    :returns: dictionary containing all accessions containing annotations of given signature
    :rtype: dict
    """
    sig_accs = {}
    for sig in sigs:
        sig_accs[sig] = output_list(sig)
        print(sig, len(sig_accs[sig]))
    return sig_accs

def retrieve_annotations(accs, sigs, download_dir, download = False):
    """
    Retrieves annotations from different databases included in InterPro

    :param accs: all the accessions of the proteins of interest
    :param sigs: all the protein signatures of interest
    :param download: indicates whether tsv files need to be downloaded
    :type accs: list
    :type sigs: list
    :returns: dataframe containing annotations for all proteians and all signatures
    :rtype: pandas.DataFrame
    """
    if download == True:
        download_links = []
        for acc in accs:
            download_links.append('https://www.ebi.ac.uk/interpro/legacy/protein/{}?export=tsv'.format(acc))
        counter = 0
        for url in download_links:
            filename = accs[counter] + '.text'
            filepath = os.path.join(download_dir, filename)
            if not os.path.isfile(filepath): # making sure not to overwrite files with identical names
                urllib.request.urlretrieve(url, filepath) # retreiving TSV files and saving them
            counter += 1
        print(counter)
    tables = []
    for file in os.listdir(download_dir):
        with open (os.path.join(download_dir, file)) as tsv:
            row = []
            for line in csv.reader(tsv, delimiter = "\t"):
                row.append(line[:15])
            tables.append(pd.DataFrame.from_records(row[1:], columns = row[0]))
    interpro_df = pd.concat(tables)
    # removing all letter cases and gaps replaced with underscores to simplify downstream analysis when calling columns
    interpro_df.columns = interpro_df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    ann_df = interpro_df[interpro_df.signature_accession.isin(sigs)]
    ann_df.reset_index(inplace = True)
    ann_df = ann_df[["protein_accession","start_position", "stop_position", "sequence_length", "signature_accession"]]
    ann_df.rename(columns = {"protein_accession": "accession", "sequence_length": "length", "signature_accession": "source","start_position": "start","stop_position": "end"}, inplace = True)
    ann_df.start = ann_df.start.astype(int)
    ann_df.length = ann_df.length.astype(int)
    ann_df.end = ann_df.end.astype(int)
    ann_df["repeat_id"] = ann_df.accession + "/" + ann_df.start.apply(str) + "-" + ann_df.end.apply(str)
    return ann_df

def swissprot_dict(db = "/cluster/gjb_lab/2394007/db/swissprot_rev_human.fasta"):
    """
    Retrieves sequences and their data from Swiss-Prot

    :param db: absolute path to a fasta file containing sequences, Swiss-Prot database by default
    :type db: str
    :returns: dictionary containing the sequence id, description and sequence for all proteins in Swiss-Prot
    :rtpe: dict
    """
    swissprot = Bio.SeqIO.parse(db, "fasta")
    proteins = {}
    for protein in swissprot:
        acc = protein.id.split("|")[1]
        proteins[acc] = {}
        proteins[acc]["id"] = protein.id
        proteins[acc]["desc"] = protein.description
        proteins[acc]["seq"] = protein.seq
    return proteins

def df_to_dict(df): ### USE TO STRATIFY BY DIFFERENT SOURCES
    """
    Transforms a dataframe containint repeat annotations from different databases
    into a dictionary with the desired format

    Parameters:

    df (DataFrame): dataframe containing information about repeat annotations

    Returns:

    def_dict (dictionary): dictionary containing the same information as df

    """
    def_dict = {}
    species = df.groupby("species")
    for spec, accs in species:
        def_dict[spec] = {}
        accs_df = accs.groupby("accession")
        for acc, reps in accs_df:
            srcs = reps.groupby("source")
            def_dict[spec][acc] = {}
            for src, anns in srcs:
                anns.index = anns.repeat_id
                def_dict[spec][acc][src] = pd.DataFrame.to_dict(anns, orient = "index")
    return def_dict

def df_to_dict2(df): ### USE WHEN STRATIFICATION BY DIFFERENT SOURCES IS NOT DESIRED
    """
    Transforms a dataframe containint repeat annotations from different databases
    into a dictionary with the desired format

    :param df: (DataFrame): information about repeat annotations
    :type df: pandas.DataFrame
    :returns: contains the same information as df
    :rtype: dict

    """
    def_dict = {}
    specs = df.groupby("species")
    for spec, reps in specs:
        def_dict[spec] = {}
        accs = reps.groupby("accession")
        for acc, reps in accs:
            reps.index = reps["repeat_id"]
            def_dict[spec][acc] = pd.DataFrame.to_dict(reps, orient = "index")
    return def_dict

def get_repeat_seqs3(prots_dict, seq_df, out_file, out_fmt = "fasta"):
    """
    Writes sequences to a file

    Parameters:

    :param prots_dict: contains basic information for a set of sequences
    :param seq_df: contains repeat annotations 
    :param out_file: absolute path to the output file
    :param out_fmt: output format (fasta by default)
    :type prots_dict: dict
    :type seq_df: pandas.DataFrame
    :type out_file: str
    :type out_fmt: str
    """
    rep_recs = []
    s = 0
    seq_dict = df_to_dict2(seq_df)
    cols = list(seq_df.columns)
    for spec, spec_accs in seq_dict.items():
        for acc, acc_reps in spec_accs.items():
            for rep, info in acc_reps.items():
                start = info["start"]
                end = info["end"]
                rep_seq = prots_dict[acc]["seq"][start - 1 : end]
                if len(rep_seq) == 0:
                    print(rep)
                    continue
                rep_obj = Bio.SeqRecord.SeqRecord(rep_seq)
                if "origin" in cols:
                    rep_obj.id = prots_dict[acc]["id"] + "_{}".format(info["origin"]) + "/" + str(start) + "-" + str(end)
                else:
                    rep_obj.id = prots_dict[acc]["id"] + "_{}".format(info["source"]) + "/" + str(start) + "-" + str(end)
                rep_obj.description = prots_dict[acc]["desc"]
                rep_obj.annotations.update({"accession": acc})
                rep_recs.append(rep_obj)
                s += 1
    print(s)
    Bio.SeqIO.write(rep_recs, out_file, out_fmt)

def get_df_from_dict(dic):
    """
    Transforms a dictionary into a dataframe with the desired format

    Parameters:

    :param dic: contains information about repeat annotations/hits
    :type dic: dict
    :returns: dataframe containing the same information as the dict
    :rtype: pandas.DataFrame
    """
    dfs = []
    specs = list(dic.keys())
    for spec in specs:
        accs = list(dic[spec].keys())
        for acc in accs:
            df = pd.DataFrame.from_dict(dic[spec][acc], orient = "index")
            #df['repeat_id'] = df.index
            dfs.append(df)
    max_df = pd.concat(dfs,sort = False)
    max_df.index = (range(0,len(max_df)))
    max_df.start = max_df.start.astype(int)
    max_df.end = max_df.end.astype(int)
    if "length" in list(max_df.columns):
        max_df.length = max_df.length.astype(int)
    max_df.repeat_length = max_df.repeat_length.astype(int)
    return max_df

def get_n_reps(dic):
    """
    Calculates number of annotations/hits in a dictionary

    :param dic: information about repeat annotations
    :type dic: dict
    :returns: number of annotations/hits in a dictionary
    :rtype: int
    """
    h = 0
    specs = list(dic.keys())
    for spec in specs:
        accs = list(dic[spec].keys())
        for acc in accs:
            h += len(dic[spec][acc])
    return h

def jaccard(e_i, e_best):
    """
    Calculates the Jaccard index between two entries in a dictionary

    :param e_i: dictionary containing information about a repeat annotation
    :param e_best: dictionary containing information about the best
    annotation so far for a given real repeat
    :type e_i: dict
    :type e_best: dict
    :returns: Jaccard index between the two dictionary entries (repeat annotation/hit)
    :rtype: float
    """
    l1 = list(range(e_best["start"], e_best["end"] + 1))
    l2 = list(range(e_i["start"], e_i["end"] + 1))
    intersection = len(list(set(l1).intersection(l2)))
    union = (len(l1) + len(l2)) - intersection
    j = float(intersection) / union
    return j

def metric(e_i, e_best):
    """
    Calculates an overlapping metric between two entries in a dictionary

    :param e_i: dictionary containing information about a repeat annotation
    :param e_best: dictionary containing information about the best
    annotation so far for a given real repeat
    :type e_i: dict
    :type e_best: dict
    :returns: overlapping proportion according to m between the two dictionary entries
    """
    l1 = list(range(e_best["start"], e_best["end"]))
    l2 = list(range(e_i["start"], e_i["end"]))
    intersection = len(list(set(l1).intersection(l2)))
    m = float(intersection) / len(l2)
    return m

def get_overlapping_p(dic):
    """
    Calculates number of absolute and relative frequency of overlapping annotaions/hits
    on dic

    :param dic: dictionary containing information about repeat annotations
    :type dic: dict
    :returns: the absolute and relative frequency of annotations that overlap with another one,
    respectively (if A and B overlap, 1 is added, not 2)
    :rtype: (int, float)
    """
    overlap = 0
    specs = list(dic.keys())
    for spec in specs:
        accs = list(dic[spec].keys())
        for acc in accs:
            reps = list(dic[spec][acc].keys())
            if len(reps) > 1:
                for i in range(1,len(reps)):
                    if dic[spec][acc][reps[i]]['start'] <= dic[spec][acc][reps[i-1]]['end']:
                        overlap += 1
    return overlap, overlap/get_n_reps(dic)

def get_annotation_overlaps(anns_dict, accs, sigs, met):
    """
    Calculates number of absolute and relative frequency of overlapping annotaions/hits
    in anns_dict

    :param anns_dict: dictionary containing information about repeat annotations
    :param accs: all the accessions of the proteins in anns_dict
    :param sigs: all protein signatures in anns_dict
    :param met: metric m or j to be used to calculate overlapping between annotations
    :type anns_dict: dict
    :type accs: list
    :type sigs: list
    :type met: str
    :returns: list containig all the overlapping metrics, for every comparison considered,
    dictionary where the key is a tuple of the annotations compared and the value is the overlapping metric
    :rtype: (list, dict)

    """
    overlaps = []
    overlaps_dict = {}
    species = list(anns_dict.keys())
    for spec in species:
        anns_spec = anns_dict[spec]
        accs = list(anns_spec.keys())
        for acc in accs:
            anns_acc = anns_spec[acc]
            sigs = list(anns_acc.keys())
            for sig in enumerate(sigs):
                anns_sig = anns_acc[sig[1]]
                for rep in anns_sig:
                    for other_sig in sigs[sig[0]+1:]:
                        for other_rep in anns_acc[other_sig]:
                            if met == 'j':
                                m = jaccard(anns_acc[sig[1]][rep],anns_acc[other_sig][other_rep])
                            elif met == 'm':
                                m = metric(anns_acc[sig[1]][rep],anns_acc[other_sig][other_rep])
                            overlaps.append(m)
                            overlaps_dict[((spec + '_' + sig[1] + '_' + rep),(spec + '_' + other_sig + '_' + other_rep))] = m
    return overlaps, overlaps_dict

def get_overlapping_proteins(df):
    """
    Calculates the absolute and relative frequency of overlapping annotations/hits
    for every protein in df

    :param df: dataframe containing information about repeat annotations
    :type df: pandas.DataFrame
    :returns: dictionary with key accession number and value a tuple containing
    absolute and relative frequency of overlapping annotations/hits
    :rtype: dict
    """
    accs = list(df.accession.unique())
    over = {}
    for acc in accs:
        df_acc = df[df.accession == acc]
        acc_dict = df_to_dict2(df_acc)
        ov, pov = get_overlapping_p(acc_dict)
        over[acc] = (ov, pov)
    return over

def merge_dfs(dic_r, dic_s, sig_r, sig_s, t, met = 'j'):
    """
    Merges annotations from two different databases with different signatures 
    within two different signatures using an overlapping metric and a threshold

    :param dic_r: reference annotations
    :param dic_s: sample annotations
    :param sig_r: reference signature
    :param sig_s: sample signature
    :param t: overlapping metric threshold
    :param met: overlapping metric to be used, j by default
    :type dic_r: dict
    :type dic_s: dict
    :type sig_r: str
    :type sig_s: str
    :type t: float
    :type met: str
    :returns: consensus annotations resulting from the merge of dic_r and dic_s,
    list containing all different accession numbers in rows
    :rtype; (dict, list)

    """
    rnp = 0 # repeats in new protein
    bas = 0 # better annotated repeats in sample
    arp = 0 # absent repeats in known protein
    bar = 0 # better annotated repeats in reference
    rows = {}
    species = list(set(list(dic_r.keys())+list(dic_s.keys())))
    for spec in species:
        rows[spec] = {}
        if spec in dic_r and spec in dic_s:
            accs = list(set(list(dic_r[spec].keys())+list(dic_s[spec].keys())))
            for acc in accs:
                rows[spec][acc] = {}
                if acc in dic_r[spec] and acc in dic_s[spec]:
                    reps_r = list(dic_r[spec][acc][sig_r].keys())
                    reps_s = list(dic_s[spec][acc][sig_s].keys())
                    for i in range(0, len(reps_s)):
                        best_rep = ''
                        best_m = -1
                        for j in range(0,len(reps_r)):
                            if met == 'm':
                                m = metric(dic_s[spec][acc][sig_s][reps_s[i]], dic_r[spec][acc][sig_r][reps_r[j]])
                            elif met == 'j':
                                m = jaccard(dic_s[spec][acc][sig_s][reps_s[i]], dic_r[spec][acc][sig_r][reps_r[j]])
                            #if acc in ["Q10EL1", "Q86AT8", "Q9J5H8"]:
                                #print(m, best_m, dic_s[spec][acc][sig_s][reps_s[i]]["repeat_id"], dic_r[spec][acc][sig_r][reps_r[j]]["repeat_id"])
                            if m > best_m:
                                best_rep = reps_r[j]
                                best_rep_j = j
                                best_m = m
                        if best_m > t:
                            if dic_s[spec][acc][sig_s][reps_s[i]]['repeat_length'] > dic_r[spec][acc][sig_r][best_rep]['repeat_length']:
                                if sig_s == 'Uniprot' and dic_s[spec][acc][sig_s][reps_s[i]]['start'] != dic_r[spec][acc][sig_r][best_rep]['start']:
                                    continue
                                else:
                                    overlap = False
                                    for rep in dic_r[spec][acc][sig_r].keys():
                                        j = jaccard(dic_s[spec][acc][sig_s][reps_s[i]], dic_r[spec][acc][sig_r][rep])
                                        if j > 0:
                                            overlap = True
                                            break
                                    if overlap == False:
                                        bas += 1
                                        rows[spec][acc][reps_s[i]] = dic_s[spec][acc][sig_s][reps_s[i]]
                                        reps_r.pop(best_rep_j)
                                    else:
                                        continue
                            else:
                                rows[spec][acc][best_rep] = dic_r[spec][acc][sig_r][best_rep]
                                reps_r.pop(best_rep_j)
                                bar += 1
                        else:
                            #print(best_m, t, dic_s[spec][acc][sig_s][reps_s[i]]) ###test
                            rows[spec][acc][reps_s[i]] = dic_s[spec][acc][sig_s][reps_s[i]]
                            arp += 1
                    if len(reps_r) != 0:
                        for rep in reps_r:
                            rows[spec][acc][rep] = dic_r[spec][acc][sig_r][rep]
                            bar += 1
                elif acc not in dic_r[spec]:
                    reps = list(dic_s[spec][acc][sig_s].keys())
                    for i in range(0,len(reps)):
                        rnp += 1
                        rows[spec][acc][reps[i]] = dic_s[spec][acc][sig_s][reps[i]]
                elif acc not in dic_s[spec]:
                    reps = list(dic_r[spec][acc][sig_r].keys())
                    for rep in reps:
                        rows[spec][acc][rep] = dic_r[spec][acc][sig_r][rep]
                        bar += 1
        elif spec not in dic_r:
            accs = list(dic_s[spec].keys())
            for acc in accs:
                rows[spec][acc] = {}
                sigs = list(dic_s[spec][acc].keys())
                for sig in sigs:
                    reps = list(dic_s[spec][acc][sig].keys())
                    for rep in reps:
                        rnp += 1
                        rows[spec][acc][rep] = dic_s[spec][acc][sig][rep]
        elif spec not in dic_s:
            accs = list(dic_r[spec].keys())
            for acc in accs:
                rows[spec][acc] = {}
                sigs = list(dic_r[spec][acc].keys())
                for sig in sigs:
                    reps = list(dic_r[spec][acc][sig].keys())
                    for rep in reps:
                        bar += 1
                        rows[spec][acc][rep] = dic_r[spec][acc][sig][rep]          
    anns_t = get_n_reps(rows)
    anns_sum = rnp+bas+arp+bar
    print('Total number of annotations: {}\nRepeats in new protein: {}\nBetter annotated repeats in sample: {}\nAbsent repeats: {}\nBetter annotated repeats in reference: {}\nThe addition is: {}'.format(anns_t,rnp, bas, arp, bar,anns_sum))
    return rows, accs

def j_overlap(df1, df2, ref, sample, key = "ANK", maxy = 800, out = None, plot = True):
    """
    Merges annotations from two different databases with different signatures 
    within two different signatures using an overlapping metric and a threshold

    :param df1: DatFrame containing annotations
    :param df2: DatFrame containing annotations
    :param ref: reference signature
    :param sample: sample signature
    :param key: repeat keyword
    :param maxy: maximum for the y axis
    :param out: absolute path to the output file
    :param plot: whether to plot results
    :type df1: pandas.DataFrame
    :type df2: pandas.DataFrame
    :type ref: str
    :type sample: str
    :type key: str
    :type maxy: int
    :type out: str
    :type plot: bool
    :returns: list containig all the overlapping metrics, for every comparison considered,
    dictionary where the key is a tuple of the annotations compared and the value is the overlapping metric
    :rtype: (list, dict)
    """
    max_df = pd.concat([df1, df2], sort = False)
    all_accs = list(max_df.accession.unique())
    max_dict = df_to_dict(max_df)
    all_sigs = list(max_df.source.unique())
    ovs, ovs_dict = get_annotation_overlaps(max_dict, all_accs, all_sigs, "j")
    if plot == True:
        plt.figure(figsize = (140, 60))
        
        bins = np.arange(0, 1, 0.01)
        plt.rcParams.update({"axes.linewidth": 5})
        plt.xticks(np.arange(0, 1.04, 0.05))
        plt.yticks(np.arange(0, 1600, 100))
        plt.tick_params(axis = "both" , which = "major", pad = 60, width = 5, length = 30, labelsize = 100)
        plt.xticks(fontsize = 100)
        plt.yticks(fontsize = 100)
        plt.hist(ovs_dict.values(), bins = bins, edgecolor = "black", linewidth = 5)
        plt.title("Frequency of J index between {} and {} {} annotations".format(ref, sample, key), pad = 100, fontsize = 140)
        plt.xlabel("J index", labelpad = 120, fontsize = 100)
        plt.ylabel("Frequency", labelpad = 120, fontsize = 100)
        plt.xlim(-0.01, 1.01)
        plt.ylim(0, maxy)
        if out != None:
            plt.savefig(out)
        plt.show()
    return ovs, ovs_dict

def fix_overlaps(dic, common_length = 33):
    """
    Gets rid of overlap between annotations by removing trailing residues on repeats overlapping with their adjacent

    :param dic: information about repeat annotations
    :type dic: dict
    :returns: sorted dictionary (within accession) with fixed annotations
    :rtype: OrderedDict
    """
    specs = list(dic.keys())
    wrong_reps = []
    for spec in specs:
        dic_spec = dic[spec]
        accs = list(dic_spec.keys())
        for acc in accs:
            reps = list(dic[spec][acc].keys())
            if len(reps) > 1:
                for i in range(1, len(reps)):
                    if dic[spec][acc][reps[i]]["start"] <= dic[spec][acc][reps[i - 1]]["end"]:
                        if metric(dic[spec][acc][reps[i]], dic[spec][acc][reps[i-1]]) >= 1:
                            print("{} within {}".format(dic[spec][acc][reps[i]]["repeat_id"], dic[spec][acc][reps[i-1]]["repeat_id"]))
                            diff_i = abs(common_length - dic[spec][acc][reps[i]]["repeat_length"])
                            diff_0 = abs(common_length - dic[spec][acc][reps[i-1]]["repeat_length"])
                            print(diff_i,dic[spec][acc][reps[i]]["repeat_id"], diff_0, dic[spec][acc][reps[i-1]]["repeat_id"])
                            if diff_i > diff_0:
                                print("keeping {}".format(dic[spec][acc][reps[i-1]]["repeat_id"]))
                                wrong_reps.append([spec, acc, reps[i]])
                            else:
                                print("keeping {}".format(dic[spec][acc][reps[i]]["repeat_id"]))
                                wrong_reps.append([spec, acc, reps[i - 1]])
                        elif metric(dic[spec][acc][reps[i-1]], dic[spec][acc][reps[i]]) >= 1:
                            print("{} within {}".format(dic[spec][acc][reps[i-1]]["repeat_id"], dic[spec][acc][reps[i]]["repeat_id"]))
                            diff_i = abs(common_length - dic[spec][acc][reps[i]]["repeat_length"])
                            diff_0 = abs(common_length - dic[spec][acc][reps[i-1]]["repeat_length"])
                            print(diff_i,dic[spec][acc][reps[i]]["repeat_id"], diff_0, dic[spec][acc][reps[i-1]]["repeat_id"])
                            if diff_i > diff_0:
                                print("keeping {}".format(dic[spec][acc][reps[i-1]]["repeat_id"]))
                                wrong_reps.append([spec, acc, reps[i]])
                            else:
                                print("keeping {}".format(dic[spec][acc][reps[i]]["repeat_id"]))
                                wrong_reps.append([spec, acc, reps[i - 1]])
                            #print("keeping {}".format(dic[spec][acc][reps[i-1]]["repeat_id"]))
                            #wrong_reps.append([spec, acc, reps[i]])
                        else:
                            #print(dic[spec][acc][reps[i]]["repeat_id"], dic[spec][acc][reps[i-1]]["repeat_id"])
                            diff = dic[spec][acc][reps[i]]["start"] - dic[spec][acc][reps[i - 1]]["end"]
                            #print(dic[spec][acc][reps[i]]["start"], dic[spec][acc][reps[i-1]]["end"])
                            rep_id, end = reps[i - 1].split("-")
                            #print(rep_id, end)
                            new_rep_id = rep_id + "-" + str(int(end) + diff -1)
                            #print(diff)
                            #print(int(end) + diff -1)
                            dic[spec][acc][new_rep_id] = dic[spec][acc][reps[i-1]]
                            dic[spec][acc][new_rep_id]["end"] =  dic[spec][acc][reps[i-1]]["end"] + diff - 1
                            dic[spec][acc][new_rep_id]["repeat_length"] =  dic[spec][acc][new_rep_id]["end"] + 1 - dic[spec][acc][new_rep_id]["start"]
                            dic[spec][acc][new_rep_id]["repeat_id"] =  new_rep_id
                            #print(dic[spec][acc][reps[i-1]]["repeat_id"], dic[spec][acc][reps[i]]["repeat_id"])
                            #print(new_rep_id)
                            #print(dic[spec][acc][new_rep_id])

                            del dic[spec][acc][reps[i - 1]]
    print(wrong_reps)
    for wrong_rep in wrong_reps:
        if wrong_rep[2] in dic[wrong_rep[0]][wrong_rep[1]]:
            del dic[wrong_rep[0]][wrong_rep[1]][wrong_rep[2]]
    sorted_dic = {}
    for spec, accs in dic.items():
        sorted_dic[spec] = {}
        for acc, reps in dic[spec].items():
            sorted_reps = dict(OrderedDict(sorted(dic[spec][acc].items(), key = lambda x: x[1]["end"], reverse = False)))
            sorted_dic[spec][acc] = sorted_reps
    return sorted_dic

