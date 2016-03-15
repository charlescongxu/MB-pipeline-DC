MB-pipeline-DC/README.md

Linux pipeline for metabarcoding.  
available at https://sourceforge.net/projects/metabarcode-pipe/

required software:  
* blast+ (executables makeblastdb, blastdbcmd)
* usearch (http://www.drive5.com/usearch/download.html), old version 4.2.66_i86linux32 is recommended
* raxmlHPC (version > 7.2.8-ish)
* multiple alignment software such as clustalo, clustalw, prank, tcoffee, mafft, fsa or muscle
* pynast, mothur or sina for aligning query sequences to reference alignment

required scripts:  
* create_fasta_database_from_genbank_flatfiles.pl
* format_conversion.pl
* multiple_sequence_splitter.pl
* parse_hits.pl
* parse_ncbi_tax_database.pl
* parse_taxa_overlapping_newick_and_fasta.pl
* parse_taxon_from_fastafile.pl
* preprocess_fasta_database.pl
* species_filter.pl

required files:  
* BE_mammal_supertree.nwk (a prevously published supertree phylogeny of mammals by Bininda-emonds)
