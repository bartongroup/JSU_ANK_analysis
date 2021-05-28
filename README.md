# Ankyrin repeats in context with human genetic variation
This repository contais a set of notebooks and libraries used to analyse the distribution of missense variants within the ankyrin repeat motif and explain the observed patterns with structural features such as intra-domain contacts, residue solvent accessibility (RSA) or protein-protein interactions. Included in this repository, in the _/files_ directory, are the main input files needed to run the notebook. These files are the multiple sequence alignment containing the 7,407 reviewed repeat sequences used in this analysis as well as the tables resulting from the packages _VarAlign_ and _ProIntVar_.

The data_extraction.ipynb notebook contains the functions necessary to download all the ankyrin repeat annotation records found in reviewed proteins from UniProt and InterPro. These include manually curated annotations from UniProt as well as annotations from the ProSitem SMART, PFAM and PRINTS databases. Once downloaded, the records are merged in a dataframe and saved.

In the database_integration notebook, the sets of annotations belonging to different database signatures are merged sequentially in a specific order. This order was established according to the number and quality of the annotations for each database. From higher to lower confidence: Prosite (PS50088), SMART (SM00248), UniProt, PRINTS (PR01415), PFAM (PF00023) and PFAM (PF13606). This procedure, results in a non-redundant set of 7,407 ankyrin repeat sequences.


