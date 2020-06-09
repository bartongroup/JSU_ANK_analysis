import retrieve_data_allsp
import re
import os
import subprocess
import time
import Bio.SeqIO
import pandas as pd
from Bio.Align.Applications import ClustalOmegaCommandline
import random
import importlib
importlib.reload(retrieve_data_allsp)

def get_lens_dict(df):
    """
    """
    return df.repeat_length.value_counts().to_dict()

def order_seq_files(df, suf, seq_dir, conf_db, common_length = None, extra = "", seq_fmt = "fasta"):
    """
    """
    lens_dict = get_lens_dict(df)
    best_db = conf_db[0]
    lens = list(lens_dict.keys())
    if common_length == None:
        common_length = lens[0]
    shorter = sorted([l for l in lens if l < common_length], reverse = True)
    longer = sorted([l for l in lens if l > common_length])
    all_lens =  list([common_length,]) + shorter + longer
    ordered_files = []
    for seq_len in all_lens:
        for db in conf_db:
            file = os.path.join(seq_dir,"{}_{}_{}{}.{}".format(suf, db, seq_len, extra, seq_fmt))
            if os.path.isfile(file):
                ordered_files.append(file)
    return ordered_files

def get_occ_cols(aln_in, fmt_in): #returns a dictionary of the occupoancy of every column in the alignment
    """
    """
    aln = Bio.SeqIO.parse(aln_in, fmt_in)
    occ = {}
    for rec in aln:
        seq = str(rec.seq)
        for i in range(0, len(seq)):
            if i + 1 not in occ:
                occ[i + 1] = 0
            if seq[i] != "-":
                occ[i + 1] += 1
    return occ

def get_cons_cols(aln_in, fmt_in, t_occ = 0.5): #returns the columns with a relative occupancy greater than a threshold
    """
    """
    n_seq = get_n_seq(aln_in, fmt_in)
    aln_occ = get_occ_cols(aln_in, fmt_in)
    cons_cols = [col for col, occ in aln_occ.items() if occ/n_seq > t_occ]
    return cons_cols

def get_n_seq(seq_file, seq_format = "fasta"):
    """
    """
    seqs = Bio.SeqIO.parse(seq_file,seq_format)
    n = 0
    for seq in seqs:
        n += 1
    return n

def get_sequences(df, suf, seq_dir, seq_db):
    """
    """
    swissprot = retrieve_data_allsp.swissprot_dict(seq_db)
    if "origin" in list(df.columns):
        orig_col = 'origin'
        origins = list(df.origin.unique())
    else:
        orig_col = 'source'
        origins = list(df.source.unique())
    for origin in origins:
        df_origin = df[df[orig_col] == origin]
        lens = list(df_origin.repeat_length.unique())
        for seq_len in lens:
            df_len = df_origin[df_origin.repeat_length == seq_len]
            retrieve_data_allsp.get_repeat_seqs3(swissprot, df_len, os.path.join(seq_dir,'{}_{}_{}.fasta'.format(suf, origin, seq_len)),'fasta')
            print("Extracted {} sequences from {} with a length of {}".format(suf,origin,seq_len))

def is_terminal_gap(cons_seq):
    """
    """
    p = re.compile("^-*\w+-*$")
    m = p.match(cons_seq)
    if m:
        return True
    else:
        return False

def wait4output(outfile_path):
    """
    """
    while not os.path.isfile(outfile_path):
        #print('not yet')
        time.sleep(1)
    print('{} has been created!'.format(outfile_path.split('/')[-1]))

def qsub_aln_job(cline, profile_path, t):
    """
    """
    script_header = "#!/bin/bash\n#$ -V\n#$ -cwd\n#$ -b n\n#$ -N clustalo\n#$ -M 2394007@dundee.ac.uk\n#$ -m a\n#$ -pe smp {}\n".format(t)
    script = script_header + str(cline) + "\n"
    aln_dir = os.path.dirname(profile_path)
    sh_path = os.path.join(aln_dir,'clustalo.sh')
    sh_fh = open(sh_path,'w')
    sh_fh.write(script)
    sh_fh.close()
    proc = subprocess.Popen(["/opt/uge/bin/lx-amd64/qsub", sh_path],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = proc.communicate()
    print("job submitted")
    return o,e

def get_gapless_seqs2(aln_in, fmt_in, cons_cols, seq_len, seq_sig, common_length):
    """
    """
    n_cons = len(cons_cols)
    gapless = []
    gappy = []
    aln = Bio.SeqIO.parse(aln_in, fmt_in)
    for rec in aln:
        sig = rec.id.split('|')[-1].split('/')[0].split('_')[-1]
        seq = str(rec.seq)
        gapless_seq = seq.replace('-','')
        cons_seq = ''.join([seq[i] for i in range(0, len(seq)) if i+1 in cons_cols])
        gaps = cons_seq.count('-')
        #print(len(gapless_seq), n_cons)
        #print(len(gapless_seq), common_length)
        #print("--")
        if len(gapless_seq) == seq_len and sig == seq_sig :
            if len(gapless_seq) == common_length:
                if gaps == 0:
                    gapless.append(rec)
                else:
                    gappy.append(rec)
            elif len(gapless_seq) < common_length:
                if gaps == 0:
                    gapless.append(rec)
                elif gaps != 0:
                    if gaps  > (n_cons - len(gapless_seq)):
                        if is_terminal_gap(cons_seq) == False:
                            gappy.append(rec)
                        elif is_terminal_gap(cons_seq) == True:
                            gappy.append(rec)
                    elif gaps  <= (n_cons - len(gapless_seq)):
                        if is_terminal_gap(cons_seq) == False:
                            gappy.append(rec)
                        elif is_terminal_gap(cons_seq) == True:
                            gapless.append(rec)
            elif len(gapless_seq) > common_length:
                if gaps  > 0:
                    #print(is_terminal_gap(cons_seq))
                    if is_terminal_gap(cons_seq) == False:
                        #print(cons_seq)
                        gappy.append(rec)
                    elif is_terminal_gap(cons_seq) == True:
                        cons_seq_gapless = cons_seq.replace("-", "")
                        if len(cons_seq_gapless)/len(gapless_seq) < 0.5:
                            gappy.append(rec)
                        else:
                            gapless.append(rec)
                else:
                    gapless.append(rec)
        else:
            gapless.append(rec)
    print(len(gappy), len(gapless))
    if len(gappy) != 0 and len(gapless) != 0:
        Bio.SeqIO.write(gapless, aln_in[:-4]+'_gapless.sto', 'stockholm')
        Bio.SeqIO.write(gappy, aln_in[:-4]+'_gappy.fasta', 'fasta')
        return [aln_in[:-4]+'_gapless.sto',aln_in[:-4]+'_gappy.fasta']
    elif len(gappy) != 0 and len(gapless) == 0:
        #Bio.SeqIO.write(gapless, aln_in[:-4]+'_gapless.sto', 'stockholm')
        Bio.SeqIO.write(gappy, aln_in[:-4]+'_gappy.sto', 'stockholm')
        return [aln_in[:-4]+'_gappy.fasta',]
    elif len(gappy) == 0 and len(gapless) != 0:
        Bio.SeqIO.write(gapless, aln_in[:-4]+'_gapless.sto', 'stockholm')
        print('file created')
        return [aln_in[:-4]+'_gapless.sto',]
    else:
        return None

def refine_alignment2(aln_file, seq_len, seq_sig, best_db, common_length, cons_cols = None, it = 3, t_occ = 0.5, t = 1, aln_fmt_in = 'stockholm', aln_fmt_out = 'st', start = True):
    """
    """
    for i in range(0, it):
        print("Iteration " + str(i))
        if start == True:
            aln_seqs = get_n_seq(aln_file, aln_fmt_in)
            cons_cols = get_cons_cols(aln_file, aln_fmt_in, t_occ = 0.5)
        #else:
            #cons_cols = get_cons_cols_from_seq_list('/cluster/gjb_lab/2394007/all_sp_aln/alns2/anks_{}_{}_gapless.sto'.format(best_db, common_length), 'stockholm',
                                        #aln_file, 'stockholm')
        gapless_info = get_gapless_seqs2(aln_file, aln_fmt_in, cons_cols, seq_len, seq_sig,common_length)
        print(len(gapless_info))
        if len(gapless_info) == 1 and gapless_info[0].split("/")[-1].split(".")[0].split("_")[-1] == "gapless":
            print("There are no gappy sequences in {}".format(aln_file))
            return gapless_info[0]
        elif len(gapless_info) == 1 and gapless_info[0].split("/")[-1].split(".")[0].split("_")[-1] == "gappy":
            aln_file = gapless_info[0]
            #print("Refining alignment: {} {} {}".format(gapless_info[0], gapless_info[1], aln_file))
            return aln_file
        #elif gapless_info == True and it == 0 and start == True:
        #    cons_cols 
        else:
            aln_file = gapless_info[0][:-11]+'2gapless.sto'
            #print("Refining alignment: {} {} {}".format(gapless_info[0], gapless_info[1], aln_file))
            align_seqs2profile(gapless_info[0], gapless_info[1], aln_file, outfmt = aln_fmt_out, t = t)
    return aln_file

def align_seqs2profile(profile_path, infile_path, outfile_path, outfmt = 'st', t = 1):
    """
    """
    n_seq = get_n_seq(infile_path)
    if n_seq == 1:
        clustalomega_cline = ClustalOmegaCommandline(
            profile1 = profile_path, profile2 = infile_path,
            isprofile = True, seqtype = 'protein',
            outfile = outfile_path, outfmt = outfmt,
            verbose = False, auto = False, force = True, threads = t)
    else:
        clustalomega_cline = ClustalOmegaCommandline(
            profile1 = profile_path, isprofile = True,
            infile = infile_path, seqtype = 'protein',
            outfile = outfile_path, outfmt = outfmt,
            verbose = False, auto = False, force = True, threads = t)
    print(clustalomega_cline)
    o,e = qsub_aln_job(clustalomega_cline, profile_path, t)
    wait4output(outfile_path)
    time.sleep(3)
    #clustalomega_cline()
    #return outfile_path

def align_profile2profile(profile1_path, profile2_path, outfile_path, outfmt = 'st', t = 1):
    """
    """
    clustalomega_cline = ClustalOmegaCommandline(
        profile1 = profile1_path, profile2 = profile2_path,
        isprofile = True, seqtype = 'protein',
        outfile = outfile_path, outfmt = outfmt,
        verbose = False, auto = False, force = True, threads = t)
    print(clustalomega_cline)
    o,e = qsub_aln_job(clustalomega_cline, profile1_path, t)
    wait4output(outfile_path)
    time.sleep(3)
    return outfile_path
    #clustalomega_cline()

def align_seqs2hmm(hmm_path, infile_path, outfile_path, outfmt = 'st', t = 1):
    """
    """
    clustalomega_cline = ClustalOmegaCommandline(
        hmm_input = hmm_path, infile = infile_path,
        seqtype = 'protein',
        outfile = outfile_path, outfmt = outfmt,
        verbose = False, auto = False, force = True, threads = t)
    print(clustalomega_cline)
    o,e = qsub_aln_job(clustalomega_cline, hmm_path, t)
    wait4output(outfile_path)
    time.sleep(3)
    #clustalomega_cline()

def align_seqs(in_file, aln_dir, t = 1):
    """
    """
    out_file = os.path.join(aln_dir, in_file.split('/')[-1].replace(".fasta",".sto"))
    clustalomega_cline = ClustalOmegaCommandline(
        infile = in_file, seqtype = 'protein',
        outfile = out_file, outfmt = 'st',
        verbose = False, auto = False, threads = t,
        force = True)
    print(clustalomega_cline)
    #clustalomega_cline()
    o,e = qsub_aln_job(clustalomega_cline, aln_dir, t)
    wait4output(out_file)
    time.sleep(3)
    return out_file

def generate_random_ids(aln_in, aln_fmt_in, pref = "guide"):
    aln =  Bio.SeqIO.parse(aln_in, aln_fmt_in)
    clean_recs = []
    for rec in aln:
        new_id = random_with_N_digits(15)
        rec.id = "{}_{}".format(pref,str(new_id))
        rec.name = "ank_{}_{}".format(pref,str(new_id))
        rec.description = ""
        rec.annotations["start"] = 1
        rec.annotations["end"] = 33
        clean_recs.append(rec)
    aln_path, aln_fmt = aln_in.split(".")
    Bio.SeqIO.write(clean_recs, aln_path + "_uuids.{}".format(aln_fmt), aln_fmt_in)
    return aln_path + "_uuids.{}".format(aln_fmt)
#generate_random_ids("/cluster/gjb_lab/2394007/all_sp_aln/trial/anks33_guide.sto", "stockholm")

def random_with_N_digits(n):
    range_start = 10**(n-1)
    range_end = (10**n)-1
    return random.randint(range_start, range_end)

def remove_empty_cols(aln_in, fmt_in):
    aln =  Bio.SeqIO.parse(aln_in, fmt_in)
    occ = get_occ_cols(aln_in, fmt_in)
    not_empty = [k for k, v in occ.items() if v != 0]
    clean_recs = []
    for rec in aln:
        new_seq = ''.join([str(rec.seq)[i-1] for i in not_empty])
        rec.seq = Bio.Seq.Seq(new_seq)
        clean_recs.append(rec)
    aln_path, aln_fmt = aln_in.split(".")
    out_file = aln_path + "_clean.{}".format(aln_fmt)
    Bio.SeqIO. write(clean_recs, out_file, fmt_in)
    return out_file

def remove_guides(aln_in, fmt_in, pref = "guide"):
    aln =  Bio.SeqIO.parse(aln_in, fmt_in)
    clean_recs = []
    for rec in aln:
        if rec.id.startswith(pref):
            continue
        else:
            clean_recs.append(rec)
    aln_path, aln_fmt = aln_in.split(".")
    out_file = aln_path + "_noguide.{}".format(aln_fmt)
    Bio.SeqIO.write(clean_recs, out_file, fmt_in)    
    return out_file

def guided_alignment3(df, wd, suf, conf_db, common_length, it = 3, t = 1, keep_gappy = False, get_seqs = False, seq_db = "/cluster/gjb_lab/2394007/db/swissprot_rev_human.fasta"):
    """
    """
    gappy_seqs = 0
    best_db = conf_db[0]
    seq_dir = os.path.join(wd,'seqs')
    if not os.path.isdir(seq_dir):
        os.mkdir(seq_dir)
    if get_seqs == True:
        get_sequences(df, suf, seq_dir, seq_db)
    aln_dir = os.path.join(wd,'alns')
    if not os.path.isdir(aln_dir):
        os.mkdir(aln_dir)
    ordered_files = order_seq_files(df, suf, seq_dir, conf_db)
    for file in ordered_files:
        seq_sig = file.split('/')[-1].split('.')[0].split('_')[1]
        seq_len = int(file.split('/')[-1].split('_')[-1].split('.')[0])
        if seq_len == common_length:
            if seq_sig == best_db:
                out_file = align_seqs(file, aln_dir, t = t)
                out_file = refine_alignment2(out_file, seq_len, seq_sig, best_db, common_length, it = it, t = t)
            else:
                profile = clean_outfile
                out_file = os.path.join(aln_dir, file.split('/')[-1].replace(".fasta",".sto"))
                align_seqs2profile(profile, file, out_file, outfmt = 'st', t = t) 
                out_file = refine_alignment2(out_file, seq_len, seq_sig, best_db, common_length, it = it, t = t)
        elif seq_len == 0:
            continue
        else:
            pref, suf = file.split('/')[-1].split(".")
            out_file = os.path.join(aln_dir, pref + "_2rest.sto")
            #guide_file = generate_random_ids(guide_template, "stockholm")
            align_seqs2profile(clean_outfile, file, out_file, outfmt = 'st', t = t)
            out_file = refine_alignment2(out_file, seq_len, seq_sig, best_db, common_length, it, t = t)
            #outfile_path = os.path.join(aln_dir, pref + "_2rest.sto")
            #out_file = align_profile2profile(guide_template, out_file, outfile_path, outfmt = 'st', t = t)
        aln_seqs = get_n_seq(out_file,  'stockholm')
        cons_cols = get_cons_cols(out_file,  'stockholm', t_occ = 0.5)
        gapless_info = get_gapless_seqs2(out_file, 'stockholm', cons_cols, seq_len, seq_sig, common_length)
        if keep_gappy == False:
        	out_file = gapless_info[0]
        gappy_seqs += (aln_seqs -  get_n_seq(gapless_info[0],  'stockholm'))
        print("Left seqs: {}".format(gappy_seqs))
        #clean_outfile = remove_guides(out_file, "stockholm")
        #print("GUIDE SEQUENCES REMOVED")
        clean_outfile = remove_empty_cols(out_file, "stockholm")
        print("EMPTY COLUMNS REMOVED")
        #guide_template = clean_outfile













