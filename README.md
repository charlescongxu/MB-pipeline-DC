###MB-pipeline-DC/README.md

######Linux pipeline for metabarcoding.  
######available at https://sourceforge.net/projects/metabarcode-pipe/

required software:  
* raxmlHPC (newest version, also you will need version 7.2.7 if you want to use pplacer for taxonomic placement)
* mafft and clustalw for multiple sequence alignments
* pynast and mothur for aligning query sequences to reference alignment
* raxmlEPA (included in raxmlHPC package)
* pplacer

required scripts:  
* download_mamDB16S_Entrez.pl
* create_fasta_database_from_genbank_flatfiles_ccyx.pl
* create_fasta_database_from_genbank_flatfiles_Hsapiens.pl
* format_conversion.pl
* multiple_sequence_splitter.pl
* parse_taxa_overlapping_newick_and_fasta.pl
* preprocess_fasta_database.pl
* species_filter.pl

required files:  
* BE_mammal_supertree.nwk (a prevously published supertree phylogeny of mammals by Bininda-Emonds)
