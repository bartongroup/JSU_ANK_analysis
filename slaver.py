import os
import Bio
import Bio.SeqIO
import Bio.AlignIO
import new_aligner
import importlib
import pandas as pd

def get_cons_cols_eq(ref_file, sample_file, aln_fmt = "stockholm", t_occ = 0.4):
    """
    returns dictionary with relationship between the consensus columns of two alignments. These
    alignments must have the same number of consensus columns and in the same order.
    """
    ref_cons_cols = new_aligner.get_cons_cols(ref_file, aln_fmt, t_occ)
    sample_cons_cols = new_aligner.get_cons_cols(sample_file, aln_fmt, t_occ)
    n_cols = len(ref_cons_cols)
    cons_cols_eq = {}
    for i in range(0, n_cols):
        cons_cols_eq[sample_cons_cols[i]] = ref_cons_cols[i]
    return cons_cols_eq

def get_aln_len(aln_file, aln_fmt):
    return Bio.AlignIO.read(aln_file, aln_fmt).get_alignment_length()

def insert_gaps(seq, insertions, sample_aln_len):
    seq_list = list(str(seq))
    for i in insertions:
        seq_list.insert(i, "-")
    len_diff = sample_aln_len - len(seq_list)
    if len_diff != 0:
        seq_list.extend(["-"]*len_diff)
    return Bio.Seq.Seq("".join(seq_list))

def get_insertion_idx(ref_file, sample_file, aln_fmt = "stockholm"):
    cons_cols_eq = get_cons_cols_eq(ref_file, sample_file, aln_fmt)
    ins_ref = 0
    ins_sample = 0
    insertions_ref = []
    insertions_sample = []
    for k, v in cons_cols_eq.items():
        offset_n = (k + ins_sample) - (v + ins_ref)
        if offset_n == 0:
            continue
        elif offset_n < 0:
            for i in range(k + ins_sample - 1, v + ins_ref - 1):
                insertions_sample.append(i)
            ins_sample += abs(offset_n)
        elif offset_n > 0:
            for i in range(v + ins_ref - 1 , k + ins_sample - 1):
                insertions_ref.append(i)
            ins_ref += abs(offset_n)
    return insertions_ref, insertions_sample

def chain_to_ref(ref_file, sample_file, out_file, aln_fmt = "stockholm"):
    insertions_ref, insertions_sample = get_insertion_idx(ref_file, sample_file, aln_fmt)
    ref_aln = Bio.SeqIO.parse(ref_file, aln_fmt)
    ref_aln_len = get_aln_len(ref_file, aln_fmt) + len(insertions_ref)
    sample_aln = Bio.SeqIO.parse(sample_file, aln_fmt)
    sample_aln_len = get_aln_len(sample_file, aln_fmt) + len(insertions_sample)
    max_len = max([ref_aln_len, sample_aln_len])
    all_recs = []
    for rec_ref in ref_aln:
        rec_ref.seq = insert_gaps(rec_ref.seq, insertions_ref, max_len)
        all_recs.append(rec_ref)
    for rec_sample in sample_aln:
        if "ref" not in rec_sample.id:
            rec_sample.seq = insert_gaps(rec_sample.seq, insertions_sample, max_len)
            all_recs.append(rec_sample)
    Bio.SeqIO.write(all_recs, out_file, aln_fmt)