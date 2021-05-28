# Ankyrin repeats in context with human genetic variation
This repository contais a set of notebooks and libraries used to analyse the distribution of missense variants within the ankyrin repeat motif and explain the observed patterns with structural features such as intra-domain contacts, residue solvent accessibility (RSA) or protein-protein interactions. Included in this repository, in the _/files_ directory, are the main input files needed to run the notebook. These files are the multiple sequence alignment containing the 7,407 reviewed repeat sequences used in this analysis as well as the tables resulting from the packages _VarAlign_ and _ProIntVar_.

The data_extraction.ipynb notebook contains the functions necessary to download all the ankyrin repeat annotation records found in reviewed proteins from UniProt and InterPro. These include manually curated annotations from UniProt as well as annotations from the ProSitem SMART, PFAM and PRINTS databases. Once downloaded, the records are merged in a dataframe and saved.

In the database_integration notebook, the sets of annotations belonging to different database signatures are merged sequentially in a specific order. This order was established according to the number and quality of the annotations for each database. From higher to lower confidence: Prosite (PS50088), SMART (SM00248), UniProt, PRINTS (PR01415), PFAM (PF00023) and PFAM (PF13606). This procedure, results in a non-redundant set of 7,407 ankyrin repeat sequences.

In the notebook named upset, we plot the intersection betweeen the different sets of ankyrin repeat annotations coming from the different databases.

The methodology followed to align the 7,407 unique ankyrin repeat sequences is contained within the alignment notebook. The main functions used to align the sequences are found in the new_aligner, alignment_editing and slaver python libraries.

In variant_analysis notebook we analyse the conservation of the family and the distribution of missense variants along the ankyrin repeat motif. Specifically, we calcualte enrichment in missense variation per position in the motif.

We download the validation data for all the X-ray structures that mapped to our alignment, process them, and merge them in a single dataframe on the structure_validation.ipynb notebook.

Finally, the last two notebooks: contact_maps and structural_analysis. The former includes the creation of the contact maps, the enrichment in intra- and inter-repeat contacts per position. The latter analyses other structural features such as the relative solvent accessibility (RSA), secondary structure (SS) and enrichment in protein-substrate interactions, both on a position and surface basis.
