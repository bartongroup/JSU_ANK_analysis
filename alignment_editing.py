### EXTRACTS ALL SEQUENCES FROM DATAFRAME
import Bio
import Bio.SeqIO

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
    for acc, acc_reps in seq_dict.items():
        for rep, info in acc_reps.items():
            start = info["start"]
            end = info["end"]
            rep_seq = prots_d
            ict[acc]["seq"][start - 1 : end]
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
    
### creates dictionary to correct ids
def get_ids_dict(seqs_file, seqs_fmt):
    ids_dict = {}
    seqs_recs = list(Bio.SeqIO.parse(seqs_file,seqs_fmt))
    seqs_ids = [rec.id for rec in seqs_recs if len(rec.seq) > 0]
    for ID in seqs_ids:
        ids_dict[ID.split('|')[1]+'/'+ID.split('/')[1]] = ID
    return ids_dict

### changes the MSA sequences IDs so we can obtain missing sequences from MSA easily and realign
def change_ids(aln_file,aln_fmt,ids_dict, out_file, out_fmt):
    aln_recs = list(Bio.SeqIO.parse(aln_file,aln_fmt))
    for rec in aln_recs:
        id_key = rec.id.split('|')[1]+'/'+rec.id.split('/')[1]
        rec.id = ids_dict[id_key]
    Bio.SeqIO.write(aln_recs, out_file, out_fmt)
    
### creates a dictionary where the key is a cropped id and the value is a sequence record object
def get_seqs_dict(seqs_file, seqs_fmt):
    seqs_dict = {}
    seqs_recs = list(Bio.SeqIO.parse(seqs_file,seqs_fmt))
    for rec in seqs_recs:
        if len(rec.seq) > 0:
            seqs_dict[rec.id] = rec
    return seqs_dict

### GIVEN AN ALIGNMENT AND A SEQUENCES FILE, RETURNS THE SEQUENCES THAT ARE NOT IN ALIGNMENT
def get_missing_seqs(aln_file, seqs_file, aln_fmt = "stockholm", seqs_fmt = "fasta", write_seqs = False, out_file = None, out_fmt = "fasta"):
    aln_recs = list(Bio.SeqIO.parse(aln_file, aln_fmt))
    seqs_recs = list(Bio.SeqIO.parse(seqs_file, seqs_fmt))
    n_aln = len(aln_recs)
    n_seqs = len(seqs_recs)
    seqs_dict = get_seqs_dict(seqs_file, seqs_fmt)
    print('Alignment has {} sequences'.format(n_aln))
    print('There are {} sequencesin total'.format(n_seqs))
    print('There should be {} sequences in the output file'.format(n_seqs - n_aln))
    aln_ids = [rec.id for rec in aln_recs]
    seqs_ids = [rec.id for rec in seqs_recs]
    print(len(aln_ids), len(seqs_ids))
    missing_ids = [seq_id for seq_id in seqs_ids if seq_id not in aln_ids]
    print(len(missing_ids))
    if write_seqs:
        missing_recs = [seqs_dict[missing_id] for missing_id in missing_ids]
        Bio.SeqIO.write(missing_recs, out_file, out_fmt)
    return missing_ids