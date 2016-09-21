#######################################################################################################################################
### SECTION ONE: Create a reference database ##########################################################################################
#######################################################################################################################################
# Note: If you already have a reference multiple sequence alignment, skip to SECTION TWO.

# use 'download_mamDB16S_Entrez.pl' to entrez download all NCBI sequences in .gb format using the following query:
# 16S[All Fields]+AND+"Mammalia"[Organism]+AND+(biomol_genomic[PROP]+AND+mitochondrion[filter])
# http://www.ncbi.nlm.nih.gov/nucleotide/ for manual download
# go to "Send to" button and download as file, genbank format
perl download_mamDB16S_Entrez.pl

# output:
#  mam16S.gb

# use 'create_fasta_database_from_genbank_flatfiles_ccyx.pl' to parse reference gb files to create a single fasta format database
# -outformat 1 for accession_ncbiTaxnumber and -outformat 2 for Taxstring_GInumber
perl create_fasta_database_from_genbank_flatfiles_ccyx.pl -in mam16S.gb.gz -out mam16S.fas -outformat 2
# gunzips a gb file, extracts a fasta file from the gb, gzips the mam16S.gb file, and also omits seqs < 150 bp (this is the whole gene)

# output:
#  mam16S.fas

# manually download the Homo_sapiens refseq in .gb format
# we have to do this because the entrez download does not contain Homo_sapiens
http://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1
# go to "Send to" button and download as file, genbank format, change name and gzip to produce
#  Homo_sapiens_refseq.gb.gz

# output:
#  Homo_sapiens_refseq.gb.gz

# use 'create_fasta_database_from_genbank_flatfiles_Hsapiens.pl' to parse Homo_sapiens gb file to create a single fasta format database
# -outformat 1 for accession_ncbiTaxnumber and -outformat 2 for Taxstring_GInumber
perl create_fasta_database_from_genbank_flatfiles_Hsapiens.pl -in Homo_sapiens_refseq.gb.gz -out Homo_sapiens_refseq.fas -outformat 2

# output:
#  Homo_sapiens_refseq.fas

# use 'multiple_sequence_splitter_ccyx.pl' to parse reference fasta database and split into many files, one for each gene
# we want multiple fasta files, one for each gene
perl multiple_sequence_splitter_ccyx.pl mam16S.fas mam16S.gb

# outputs:
#  mss_log                 (empty, not sure what its for...)
#  nameparsed.<GENE_NAME>  (one for each gene)
#  resultsgenecounts       (counts how many sequences exist for each gene)

# rename nameparsed.16S as 'mamDB16S' and delete all other nameparsed files
mv nameparsed.16S mamDB16S
rm nameparsed.*

# output:
#  mamDB16S

# use 'multiple_sequence_splitter.pl' to parse Homo_sapien fasta database and split into many files, one for each gene
perl multiple_sequence_splitter_ccyx.pl Homo_sapiens_refseq.fas Homo_sapiens_refseq.gb.gz

# output:
#  mss_log                 (empty, not sure what its for...)
#  nameparsed.<GENE_NAME>  (one for each gene)
#  resultsgenecounts       (counts how many sequences exist for each gene)

# rename nameparsed.16S as 'Homo_sapiens_refseq16S' and delete all other nameparsed files
mv nameparsed.16S Homo_sapiens_refseq16S
rm nameparsed.*

# output:
#  Homo_sapiens_refseq16S

# add the 16S gene from the Homo_sapiens refseq to mamDB16S to create a full mammal 16S reference database
cat Homo_sapiens_refseq16S mamDB16S > mamDB16S_full

# remove ambiguous sequences from reference database
cat mamDB16S_full | awk '{if (index($0,">")==1) print $0; else { if (length($0)==0){printf("\n");} else { split($0, s, "") notN == 0; for (i=1; i <= length(s); i++) { if ((s[i]=="n")||(s[i]=="x")||(s[i]=="b")||(s[i]=="d")||(s[i]=="h")||(s[i]=="v")||(s[i]=="k")||(s[i]=="y")||(s[i]=="s")||(s[i]=="w")||(s[i]=="r")||(s[i]=="m")) {continue;} notN = 1; printf("%s", s[i]) } if (notN==1) { printf("\n");} } } }' > mamDB16S_full_unambig

# output:
#  mamDB16S_full_unambig

# use 'preprocess_fasta_database.pl' to remove sequences that are too short or too long
perl preprocess_fasta_database.pl -in mamDB16S_full_unambig -filter_seq_length_outliers -lower_length_limit 300 -upper_length_limit 1600

# output:
#  mamDB16S_full_unambig.genome  (this contains seqs that are too long)
#  mamDB16S_full_unambig.ng

# use 'preprocess_fasta_database.pl' to remove from reference database:
#   duplicate entries  (same sequence and acccession)
#   identical sequences if coming from the same species  (retained if from different species)
# Note: this step requires software: usearch v5.2.32
perl preprocess_fasta_database.pl -in mamDB16S_full_unambig.ng -reduce_redundency -usearch_command usearch
# gives error "Died at preprocess_fasta_database.pl line 343", which we can ignore

# output:
#  mamDB16S_full_unambig.ng.rr  # the output we want:  no duplicate entries nor identical sequences from same species
#  mamDB16S_full_unambig.ng.sorted      (sorted by decreasing length)
#  mamDB16S_full_unambig.ng.sorted.nr   (no species name, just cluster numbers)
#  mamDB16S_full_unambig.ng.uclust_log  (log file containing info of each sequence)

# use 'species_filter.pl' to obtain only the longest sequence per species
# 1 = species_filtering and 2 = species_filtering_tobycoded
# open script and make sure $filter_by_MRS	= 0        
# the other options are:  (0 = take longest sequence for each species) 
#                         (1 = take the most representative sequence for each species)
#					                (2 = take the sequence with least ambiguous data (N's etc), these make alignment more difficult)
perl species_filter.pl mamDB16S_full_unambig.ng.rr 1

# output:
#  mamDB16S_full_unambig.ng.rr.ID_filtered
#  species_filter_LOG

# remove sequences that are erraneous or may be mislabeled
# these can be identified if they are in totally incorrect positions in the final tree
# the general command is:
#  cat whatever.fasta | awk '{if (substr($0,1,6) == ">chrUn") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' > fixed.fasta
cat mamDB16S_full_unambig.ng.rr.ID_filtered | awk '{if (substr($0,1,30) == ">Myomyscus_brockmani_389617679") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' > mamDB16S_full_unambig.ng.rr.ID_filtered.fixed
cat mamDB16S_full_unambig.ng.rr.ID_filtered.fixed | awk '{if (substr($0,1,28) == ">Canis_himalayensis_34305140") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' > tmp
cat tmp | awk '{if (substr($0,1,22) == ">Canis_indica_34305148") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' > mamDB16S_full_unambig.ng.rr.ID_filtered.fixed
rm tmp

# output:
#  mamDB16S_full_unambig.ng.rr.ID_filtered.fixed

# use mafft to create a multiple sequence alignment of the reference database
mafft mamDB16S_full_unambig.ng.rr.ID_filtered.fixed > mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft

# output:
#  mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft

# Finally, we need a reference phylogeny
# this is a quick method to get a nice-ish tree using a previously published tree as a backbone
# BE_mammal_supertree.nwk is a supertree of mammals from Bininda-Emonds (ask DC for source)
# outgroup that might be used in the RAxML search: -o Tachyglossus_aculeatus
# Note: may need to first remove previous files - rm BE_mammal_supertree.nwk.overlapping *.overlapping *.overlapping.phy R*ref_tree

# use 'parse_taxa_overlapping_newick_and_fasta.pl' to obtain a subset of the supertree that contains only overlapping taxa
# Note: make sure $include_nonoverlapping_taxa_in_fasta_file = 1 or else reference alignment will be constrained by the supertree

perl parse_taxa_overlapping_newick_and_fasta.pl BE_mammal_supertree.nwk mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft

# output files:
#  BE_mammal_supertree.nwk.overlapping
#  mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping

# use 'format_conversion.pl' to convert the reference sequences from fasta to phylip format

perl format_conversion.pl mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy fasta phylip

# output files:
#  mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy

# use RAxML to create the reference phylogeny using the supertree as a backbone
# version 7.2.7 is for placement using pplacer unless you use taxit to create a reference package first, otherwise latest version will work too

# this creates the reference tree for epa, only includes the reference sequences, not the queries
# orig raxml was: raxmlHPC-SSE3 -s mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy -n ref_tree_epa -m GTRCAT -c 25 -g BE_mammal_supertree.nwk.overlapping -p 12345
# -T 3 # 3 threads, which is the max i want to use on mac. default = 2
# -s alignment
# -n name of alignment file
# -m model
# -c number of distinct rate categories
# -g name of constraint tree

raxmlHPC-PTHREADS-SSE3 -T 2 -s mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy -n ref_tree_epa -m GTRCAT -c 25 -g BE_mammal_supertree.nwk.overlapping -p 12345
# required 37 minutes on my macbook pro

# output files:
#  RAxML_log.ref_tree_epa
#  RAxML_info.ref_tree_epa
#  RAxML_bestTree.ref_tree_epa

# this creates the reference tree for pplacer (requires raxml 7.2.7)
# only necessary if you do not use taxit to create a reference package, otherwise the latest version will work too
# orig raxml: raxmlHPC-SSE3-7.2.7 -s mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy -n ref_tree_pplacer -m GTRCAT -c 25 -g BE_mammal_supertree.nwk.overlapping -p 12345
raxmlHPC-SSE3 -s mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy -n ref_tree_pplacer -m GTRCAT -c 25 -g BE_mammal_supertree.nwk.overlapping -p 12345

# output files:
#  RAxML_log.ref_tree_pplacer
#  RAxML_info.ref_tree_pplacer
#  RAxML_bestTree.ref_tree_pplacer


