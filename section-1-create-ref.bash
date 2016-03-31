#######################################################################################################################################
### SECTION ONE: Create a reference framework from mined NCBI data ####################################################################
#######################################################################################################################################
# Note: If you already have a reference multiple sequence alignment, skip to SECTION TWO.

# change your working directory to where you want to download references

cd ~/home/wangxy/charles/data/references

# download the NCBI taxonomy database

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# unzip the NCBI taxonomy database

tar -xvzf taxdump.tar.gz

# output files:
#  citations.dmp
#  delnodes.dmp
#  division.dmp
#  gc.prt
#  gencode.dmp
#  merged.dmp
#  names.dmp
#  nodes.dmp
#  readme.txt

# download the current NCBI release for the taxonomic group you want to target
# the following is written for MAMMALS
# for other taxa you need to replace the appropriate division abbreviations (e.g. gbpri*)
# the following is a list of possible abbreviations:
#  PRI - primate sequences
#  ROD - rodent sequences
#  MAM - other mammalian sequences
#  VRT - other vertebrate sequences
#  INV - invertebrate sequences
#  PLN - plant, fungal, and algal sequences
#  BCT - bacterial sequences
#  VRL - viral sequences
#  PHG - bacteriophage sequences
#  SYN - synthetic sequences
#  UNA - unannotated sequences
#  EST - EST sequences (expressed sequence tags)
#  PAT - patent sequences
#  STS - STS sequences (sequence tagged sites)
#  GSS - GSS sequences (genome survey sequences)
#  HTG - HTG sequences (high-throughput genomic sequences)
#  HTC - unfinished high-throughput cDNA sequencing
#  ENV - environmental sampling sequences

# there are 2 ways to do this
# 1. download all parts of the NCBI release within each division

wget ftp://ftp.ncbi.nih.gov/genbank/gbmam*.seq.gz
wget ftp://ftp.ncbi.nih.gov/genbank/gbpri*.seq.gz
wget ftp://ftp.ncbi.nih.gov/genbank/gbrod*.seq.gz
wget ftp://ftp.ncbi.nih.gov/genbank/gbvrt*.seq.gz

# 2. download the NCBI release of each division in parts
# this way is easier to resume downloading if process is interupted
# however, for this you need to look on the website (ftp://ftp.ncbi.nih.gov/genbank/)
# and write in here how many files in the current release (if in doubt just write a high number)

for i in {1..37};do echo "Mammal            division $i";wget ftp://ftp.ncbi.nih.gov/genbank/gbmam$i.seq.gz;done
for i in {1..54};do echo "Primate           division $i";wget ftp://ftp.ncbi.nih.gov/genbank/gbpri$i.seq.gz;done
for i in {1..31};do echo "Rodent            division $i";wget ftp://ftp.ncbi.nih.gov/genbank/gbrod$i.seq.gz;done
for i in {1..61};do echo "Other Vertebrates division $i";wget ftp://ftp.ncbi.nih.gov/genbank/gbrod$i.seq.gz;done

# output files:
#  many gzipped sequence files of the target taxonomic groups (e.g. gbmam1.seq.gz)

# use perl script 'parse_ncbi_tax_database.pl' to parse NCBI database files (names.dmp and nodes.dmp)
# for taxonomic information of specified taxonomic group
# 40674 = ncbi tax number for mammals
# script options: $parse_species_only= 1; ignores partially labeled IDs

perl parse_ncbi_tax_database.pl 40674

# output files:
#  key_<MONTH+YEAR>_Mammalia
#  parse_ncbi_tax_database_LOG

# use perl script 'create_fasta_database_from_genbank_flatfiles.pl' to parse NCBI flatfiles for creating a single fasta format database
# -outformat 1 for accession_ncbiTaxnumber and -outformat 2 for Taxstring_GInumber
# open script to modify the following settings:
#  my @genbank_divisions      = "<all GenBank database divisions used (e.g. gbmam)>"
#  my $gene_specific_analysis = 0         (depending on if you are interested in 1 or many genes)
#  my @product_of_interests   = "<gene name (e.g. 16S)>" (?only works if $gene_specific_analysis = 1?)
#  my $print_genome_taxa      = 0         (if you want to print full genomes)
#  $upper_entry_length 		    = 10000000  (?just put a big number here?)
#  $limit_taxon               = 0         (make sure this is 0, obselete I think)
#  $limit_taxon_name          = "<name of taxon you want to limit (e.g. Mammalia)>" (obselete I think)
#  $parse_protein             = 0         (if you want to also parse proteins)
#  my $verbose                = 0         (1 for debugging purposes)
#  my $accession_keyfile_ID   = "<accession abbreviation of taxonomic group (e.g. mam)>"

perl create_fasta_database_from_genbank_flatfiles.pl -out mamDB -outformat 1

# output files:
#  accession_key.<MONTH+YEAR>.mam
#  DodgyDnaSeqsFound (should hopefully contain nothing)
#  mamDB

# use perl script 'parse_taxon_from_fastafile.pl' to consolidate files containing taxonomic information and fasta database
# to create a file with just species name and accession as the fasta ID
# open script to modify the following settings:
#  $memory_efficient = 1              (1 = read only one entry into memory at a time, allows large fasta databases to be used)
#  $id_format  = 3                    (1 = tobycode, underscore, accession, 2 = ncbi_taxon_number, underscore, accession, 3 = full species name, underscore, accession, 4 = family, underscore, full species name, underscore, accession)
#  $parse_binomial_labelled_only = 1  (ignores non-standard bionomial labeles (e.g. Genus sp. something; Genus nr. something; Genus (somename) someothername; Genus sp. code12345))

perl parse_taxon_from_fastafile.pl mamDB key_Mar2016_Mammalia

# output files:
#  mamDB.parsed
#  missingtaxids (list of ommited sequences from mamDB that did not have NCBI taxon numbers listed in key_<MONTH+YEAR>_Mammalia)
#  parse_order_from_endop_fastafile_LOG

#######################################################################################################################################

# now we need to obtain a mammal 16S reference sequence

# search NCBI Nucleotide for the following:
#  mammalia [organism] AND "complete mitochondrial genome" AND refseq
# click on 'Send:' -> 'File' -> Format: 'GenBank (full)'
# save as 'refseqs.gb'

# use perl script 'create_fasta_database_from_genbank_flatfiles.pl' to parse reference sequence for creating a single fasta format database
# -outformat 1 for accession_ncbiTaxnumber and -outformat 2 for Taxstring_GInumber
# open script to modify the following settings:
#  my @genbank_divisions      = "<all GenBank database divisions used (e.g. gbmam)>"
#  my $gene_specific_analysis = 0/1       (depending on if you are interested in 1 or many genes)
#  my @product_of_interests   = ("<gene name (e.g. 16S)>")
#  my $print_genome_taxa      = 0/1       (depending on if you want to print full genomes (e.g. D. melanogaster))
#  $upper_entry_length 		    = 10000000  (?just put a big number here?)
#  $limit_taxon               = 1/0       (?what does this do?)
#  $limit_taxon_name          = "<name of taxon you want to limit (e.g. Insecta)>"
#  $parse_protein             = 1/0       (depending on if you want to also parse proteins)
#  my $verbose                = 1/0       (1 for debugging purposes)
#  my $accession_keyfile_ID   = "<accession abbreviation of taxonomic group (e.g. inv)>"

perl create_fasta_database_from_genbank_flatfiles.pl -in refseqs.gb -out refseqs.fas -outformat 2

# output files:
#  accession_key.<MONTH+YEAR>.mam (NOTE: this file overwrites the previous version created from extracting GenBank divisions)
#  DodgyDnaSeqsFound (should hopefully contain nothing)
#  refseqs.fas

# use 'multiple_sequence_splitter.pl' to parse reference fasta file and split into many files, one for each gene
# Note: this is optimized on insect genomes so names currently used may be appropriate for eukaryote mitochondrial genomes only
#       (e.g. 'SMALLSUBUNITRIBOSOMALRNA' => '12S')
# open script to modify the following settings:
#  excised_seq_length_limit = 60  (???)
#  verbose			            = 0   (1 for debugging purposes)
#  accession_or_gi		      = 2   (1 = accession, 2 = GI number)
#  parse_protein			      = 0   (default = 0, 1 = extract translated sequences from each feature) 

perl multiple_sequence_splitter.pl refseqs.fas

# outputs:
#  mss_log                 (empty, not sure what its for...)
#  nameparsed.<GENE_NAME>  (one for each gene)
#  resultsgenecounts       (counts how many sequences exist for each gene)

# rename nameparsed gene file you are interested in as 'mamrefseq16S' and delete all other nameparsed files
# this file contains a 'high quality' reference sequence for a specific gene and will be used as a query to get more sequences

mv nameparsed.16S mamrefseq16S
rm nameparsed.*

#######################################################################################################################################

# now we need to use the reference 16S sequence as a query to get more reference homologs

# use perl script 'preprocess_fasta_database.pl' to remove from reference database:
#  using -filter_seq_length_outliers flag
#   sequences that are too long  (very long sequences cause dereplication processes to crash)
#   sequences that are too short

perl preprocess_fasta_database.pl -in mamDB.parsed -filter_seq_length_outliers -lower_length_limit 200 -upper_length_limit 32000

# output files:
#  mamDB.parsed.genome  (not sure what this file is for...)
#  mamDB.parsed.ng

# use perl script 'preprocess_fasta_database.pl' to remove from reference database:
#  using -reduce_redundency flag
#   duplicate entries  (same sequence and acccession)
#   identical sequences if coming from the same species  (retained if from different species)
# Note: this step requires software: usearch4.2.66_i86linux32

perl preprocess_fasta_database.pl -in mamDB.parsed.ng -reduce_redundency -usearch_command /home/wangxy/scripts/usearch/usearch4.2.66_i86linux32

# output files:
#  mamDB.parsed.ng.rr
#  mamDB.parsed.ng.sorted      (sorted by decreasing length)
#  mamDB.parsed.ng.sorted.nr   (no species name, just cluster numbers)
#  mamDB.parsed.ng.uclust_log  (log file containing info of each sequence)

# create a BLAST database from the parsed and processed reference database using makeblastdb

makeblastdb -in mamDB.parsed.ng.rr -dbtype nucl -parse_seqids -max_file_sz 600000000

# output files (a bunch of binary files):
#  mamDB.parsed.ng.rr.nhr
#  mamDB.parsed.ng.rr.nin
#  mamDB.parsed.ng.rr.nog
#  mamDB.parsed.ng.rr.nsd
#  mamDB.parsed.ng.rr.nsi
#  mamDB.parsed.ng.rr.nsq

# use 'high quality' reference sequence of the gene you are interested in (mamrefseq16S) as a query
# for searching in the BLAST database using blastn

blastn -task blastn -db mamDB.parsed.ng.rr -out mamDB.blastout -dust no -strand both -evalue 1e-6 -query mamrefseq16S -num_threads 1 -max_target_seqs 200000 -outfmt '6 qseqid sseqid evalue pident length sstart send qframe sframe'

# output files:
#  mamDB.blastout

# use perl script 'parse_hits.pl' to parse blast results and get full sequences only 
# my $db 				               = $ARGV[0];
# my $uclust_outfile 		       = $ARGV[1];
# my $sequence_identity_cutoff = $ARGV[2];
# my $trim_retrieved_sequences = $ARGV[3];  (0 = default, no trim, 1 = trim according to positions of longest hit, 2 = trim according to left-most and right-most positions from all hits to given database sequence)
# my $sequence_length_limit 	 = $ARGV[4];  (lower limit sequence length, may need to adjust based on taxa/gene)
# my $parse_which		 	         = $ARGV[5];  (1/2 = column 1 (query) or 2 (hit, default))
# my $blastdbcmd 			         = $ARGV[6];

perl parse_hits.pl mamDB.parsed.ng.rr mamDB.blastout 50 1 1500 2 blastdbcmd

# output files:
#  mamDB.blastout.retreived
#  mamDB.blastout.retreivedINFO

# use perl script 'species_filter.pl' to filter reference of full sequences to obtain only one sequence per species
# my $format 			         = $ARGV[1]  (2 = species_filtering_tobycoded, else species_filtering())
# open script to modify the following settings:
#  identified_species_only = 1         (???)
#  $filter_by_MRS			     = 2         (0 = take longest sequence for each species) 
#                                      (1 = take the most representative sequence for each species)
#					                             (2 = take the sequence with least ambiguous data (N's etc), these make alignment more difficult)

perl species_filter.pl mamDB.blastout.retreived 1

# output files:
#  mamDB.blastout.retreived.ID_filtered
#  species_filter_LOG

#######################################################################################################################################

# now we need to do a multiple alignment of the mined reference sequences

# rename the fasta file of unaligned reference sequences to 'refs_unaligned'

cp mamDB.blastout.retreived.ID_filtered refs_unaligned

# use clustalo to do an initial alignment

/home/wangxy/scripts/Clustal/clustalo --infile=refs_unaligned --outfile=refs_unaligned.clo --outfmt=fa --threads=4 --force --log=clustalo_logfile

# output files:
#  refs_unaligned.clo

# use perl script 'format_conversion.pl' to convert clustalo alignment from fasta to phylip format

perl format_conversion.pl refs_unaligned.clo refs_unaligned.clo.phy fasta phylip

# output files:
#  refs_unaligned.clo.phy

# use RAxML to create a guide tree (using prank to create a guide tree is very slow)

/home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -D -e 1.0 -c 4 -s refs_unaligned.clo.phy -n refs_unaligned.clo_tree -m GTRCAT -p 12345

# output files:
#  RAxML_bestTree.refs_unaligned.clo_tree
#  RAxML_info.refs_unaligned.clo_tree (logfile)

# use prank with RAxML guide tree to do a final multiple alignment of the reference sequences
# Note: because prank requires modules that cannot be installed on the server, you will have to download prank and run it on a local machine

/home/wangxy/scripts/prank/bin/prank -d=refs_unaligned -o=refs_unaligned.prank -t=RAxML_bestTree.refs_unaligned.clo_tree -DNA

# output files:
#  refs_unaligned.prank.best.fas

# finally, we need a reference phylogeny
# this is a quick method to get a nice-ish tree using a previously published tree as a backbone
# BE_mammal_supertree.nwk is a supertree of mammals from Bininda-Emonds (ask DC for source)
# outgroup that might be used in the RAxML search: -o Tachyglossus_aculeatus
# Note: may need to first remove previous files - rm BE_mammal_supertree.nwk.overlapping *.overlapping *.overlapping.phy R*ref_tree

# use perl script 'parse_taxa_overlapping_newick_and_fasta.pl' to obtain only overlapping sequences in the supertree and the prank alignment of reference sequences
# Note: make sure $include_nonoverlapping_taxa_in_fasta_file = 1 or else reference alignment will be constrained by the supertree

perl parse_taxa_overlapping_newick_and_fasta.pl BE_mammal_supertree.nwk refs_unaligned.prank.best.fas

# output files:
#  BE_mammal_supertree.nwk.overlapping
#  refs_unaligned.prank.best.fas.overlapping

# use perl script 'format_conversion.pl' to convert the overlapping sequences from the alignment of references from fasta to phylip format

perl format_conversion.pl refs_unaligned.prank.best.fas.overlapping refs_unaligned.prank.best.fas.overlapping.phy fasta phylip

# output files:
#  refs_unaligned.prank.best.fas.overlapping.phy

# use RAxML to create the reference phylogeny using the supertree as a backbone

/home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -s refs_unaligned.prank.best.fas.overlapping.phy -n ref_tree -m GTRCAT -c 25 -g BE_mammal_supertree.nwk.overlapping -p 12345

# output files:
#  RAxML_log.ref_tree
#  RAxML_info.ref_tree
#  RAxML_bestTree.ref_tree
