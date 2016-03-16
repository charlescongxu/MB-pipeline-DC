# SECTION ONE: Create a reference framework from mined NCBI data

# Note: If you already have a reference multiple sequence alignment, skip to SECTION TWO.

# change your working directory to where you want to download references
cd ~/home/wangxy/charles/data/references

# download the NCBI taxonomy database
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# unzip the NCBI taxonomy database
tar -xvzf taxdump.tar.gz

# output files:
# citations.dmp
# delnodes.dmp
# division.dmp
# gc.prt
# gencode.dmp
# merged.dmp
# names.dmp
# nodes.dmp
# readme.txt

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
# many gzipped sequence files of the target taxonomic groups

# use perl script 'parse_ncbi_tax_database.pl' to parse NCBI database files (names.dmp and nodes.dmp)
# for taxonomic information of specified taxonomic group
# 40674 = ncbi tax number for mammals
# script options: $parse_species_only= 1; ignores partially labeled IDs
perl parse_ncbi_tax_database.pl 40674

# output files:
# key_<MONTH+YEAR>_Mammalia
# parse_ncbi_tax_database_LOG

# use perl script 'create_fasta_database_from_genbank_flatfiles.pl' to parse NCBI flatfiles
# for creating a single fasta format database
# open script to modify the following settings:
#  my @genbank_divisions = "<all GenBank database divisions used (e.g. gbmam)>"
#  my $gene_specific_analysis = 0/1 (depending on if you are interested in 1 or many genes)
#  my @product_of_interests = ("<gene name (e.g. 16S)>")
#  my $print_genome_taxa = 0/1 (depending on if you want to print full genomes (e.g. D. melanogaster))
#  $upper_entry_length 		= 10000000 (?just put a big number here?)
#  $limit_taxon = 1/0 (?what does this do?)
#  $limit_taxon_name = "<name of taxon you want to limit (e.g. Insecta)>"
#  $parse_protein = 1/0 (depending on if you want to also parse proteins)
#  my $verbose = 1/0 (?what does this do?)
#  my $accession_keyfile_ID = "<accession abbreviation of taxonomic group (e.g. inv)>"
perl create_fasta_database_from_genbank_flatfiles.pl -out mamDB -outformat 1
