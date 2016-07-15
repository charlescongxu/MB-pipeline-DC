###MB-pipeline-DC/README.md

######Linux pipeline for metabarcoding.  
######available at https://sourceforge.net/projects/metabarcode-pipe/

required software:
* usearch4.2.66_i86linux32
(http://www.drive5.com/usearch/download.html)
* raxmlHPC (newest version, also you will need version 7.2.7 if you want to use pplacer for taxonomic placement)
(https://github.com/stamatak/standard-RAxML)
* mafft and clustalw for multiple sequence alignments
(http://mafft.cbrc.jp/alignment/software/)
(http://www.clustal.org/clustal2/#Download)
* pynast and mothur for aligning query sequences to reference alignment
(http://biocore.github.io/pynast/)
(http://www.mothur.org/wiki/Download_mothur)
* raxmlEPA (included in raxmlHPC package)
* pplacer
(http://matsen.fhcrc.org/pplacer/)

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
