
#!/usr/bin/perl
#
# multiple_sequence_splitter.pl, 
# by Douglas Chesters, Chinese Academy of Sciences,
# based on multiple_sequence_splitter.rb by Benjamin Meyer (Peters et al. 2011. BMC Biology 9:55)
# re-written initially to increase speed on very large datasets
#  (does not individually query NCBI for each sequence)
# and standardizes gene names
# and will parse protien sequence
#
#
# to run:
# perl multiple_sequence_splitter.pl input_file_name
#
# input is a fasta format sequence file, fasta IDs are: >species_name_accession
# also reads genbank flatfiles which need specifying in the array: @zipped_genbank_flatfiles
# [need to compare to /home/douglas/usr_scripts/extract_CDS_from_genbank_flatfile_genome.pl]
#
#
# outputs:
# mss_log
# nameparsed.$name, for each gene name
# resultsgenecouts, counts how many sequences for each gene.
#
# 
# change log:
# 28 Jun 2014: bugfix's for name parsing
# 25 Aug 2014: option to use accession or GI for sequence entries
# 29 Oct 2014: warns the user if accession/gi seem not correctly specified
# 27 Dec 2014: option to parse translated protein sequences instead of DNA
# 29 Dec 2014: bugfix above
# 09 Mar 2016: bugfix when using just one sequence
# 
# 
# 
# 
# 
# 
# NOTES:
# 
# WARNING: optimized on insect mtgenomes.
#  currently the names used may be appropriate for eukaryote mtgenomes only ...
#  eg 'SMALLSUBUNITRIBOSOMALRNA' => '12S',
#
# cant figure out why its not properly parsing atp8 from insect mtgenomes??? 
#
# insect mtgenes will be named thus:
#  ATP6 CO1 CO2 CO3 CYTB NAD1 NAD2 NAD3 NAD4 NAD4L NAD5 NAD6
#
#
# 
# 
# 
# 
# 
# 
########################################################################################################



$infile = $ARGV[0];
$gbfile = $ARGV[1];

$excised_seq_length_limit 	= 60;
$verbose			= 0;
$accession_or_gi		= 2; # 1=accession; 2=gi number
$parse_protein			= 0; # default=0; 1=extract translated sequences from each feature 


##$number_of_attempts = 5; # to get url for genbank entry. not used anymore.


#my $genbank_division  		= "gbmam"; 	# gbpln 1..60; gbbct 1..94; gbenv 1..57; gbinv 1..33; gbmam 1..8; gbpri 1..45; gbrod 1..30; gbvrt 1..28; gbuna 1
#my $file_list = `ls $genbank_division*`;
#my @zipped_genbank_flatfiles = split /\n/ , $file_list;

#my @zipped_genbank_flatfiles = ("mam_refseq_complete_mtgenome.gb");
#my @zipped_genbank_flatfiles = ("mam_all_complete_mtgenomes.gb");
#my @zipped_genbank_flatfiles = ("wilting_complete_mtgenomes_fetched");
#@zipped_genbank_flatfiles = ("insect_refseq_complete_mtgenome.gb");
#@zipped_genbank_flatfiles = ("InsectComplMtGenom");
my @zipped_genbank_flatfiles = $gbfile;


#print "\n\nwill search for insect db entries in these files:@zipped_genbank_flatfiles\n\n";

unless($#zipped_genbank_flatfiles >= 0){die "\nerror, no ncbi files found\n"}


##########################################################################################################################################
#
#
#
#
#
##########################################################################################################################################




%gene_synonyms;
%gene_counts;
%store_all_accessions;
%full_id_for_accessions;
%record_ids_printed; 		# some entries are duplicates, just need one per species/accession/gene


# read a curated list of gene names and alternative names, for the purpose of congruence testing later:

######################
read_gene_synonyms();#
######################



###################
read_accessions();#
###################


system("rm nameparsed.*");


############################
parse_genbank_flatfiles2();#
############################



print "
countmissingaccessions:$countmissingaccessions
";

print "\n\nfin.\n\n";
print "
WARNING: currently names used may be appropriate for eukaryote mtgenomes only ...
";

exit;


##########################################################################################################################################
#
#
#
#
#
##########################################################################################################################################






sub read_accessions
{
print "sub read_accessions\n";
open(LOGFILE, ">mss_log");


my $file_as_string = "";
open(IN_FILTER , $infile) || die "\nerror 159. cant open $infile\n";
print "opened file:($infile)\n";
while (my $line = <IN_FILTER>)
	{$file_as_string .= $line};close(IN_FILTER);
my @all_lines = split />/, $file_as_string;
	# note $all_lines[0] contains nothing.

print $#all_lines , " seqs in file\n";


my $formatA = 0;my $formatB = 0;


for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];#print "entry number $each_line\nentry:$line\n";
	if($line =~ /^(.+)/)
		{
		my $speciesid = $1;my $complete_fasta_id = $speciesid;
		$speciesid =~ s/\n//g;$speciesid =~ s/\r//g;
		$line =~ s/^.+\n//;$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/[\s\t]//g;

		if($speciesid =~ /_[^_]+$/)
			{$speciesid =~ s/.+_([^_]+)$/$1/;
			}else{die "cant read accession in fasta ID:$line\n"};
		if($speciesid =~ /^\d+$/){$formatB++	
			}else{
			if($speciesid =~ /^\w+\d+$/){$formatA++}
			}

		$store_all_accessions{$speciesid}=$line;
		$full_id_for_accessions{$speciesid}=$complete_fasta_id;

		};


	}#for my $each_line(1 .. $#all_lines)


#close(LOGFILE);

my @test = keys %store_all_accessions;
print scalar @test , " accessions read\n";
unless(scalar @test >= 1){die "
error 205. did not read accessions from input file
"};

print "
formatA$formatA formatB$formatB 
";



if($formatA > $formatB && $accession_or_gi == 2)
	{print "\n\nWARNING. you specified gi numbers in the script. but looks like your fasta file uses accessions.\n\n"}
if($formatA < $formatB && $accession_or_gi == 1)
	{print "\n\nWARNING. you specified accessions in the script. but looks like your fasta file uses gi numbers.\n\n"}

#$accession_or_gi		= 2; # 1=accession; 2=gi number




}









##########################################################################################################################################
#
#
#
#
#
##########################################################################################################################################



sub fetch_id
{
my $gb_entry = shift;

#print "entry ($gb_entry)\n";

#my @split_entry = split /\s+CDS|tRNA|rRNA|misc_feature\s+/ , $gb_entry;
my @split_entry = split /\s+CDS\s+|\s+tRNA\s+|\s+rRNA\s+|\s+misc_feature\s+/ , $gb_entry;# bugfix 28/06/14

#     rRNA            1107..2705
#                     /product="16S ribosomal RNA"
#                     /note="L-rRNA"
#     tRNA            2706..2779
#                     /product="tRNA-Leu"
#                     /codon_recognized="UUR"
#     gene            2784..3761
#                     /gene="ND1"
#                     /db_xref="GeneID:18667249"
#     CDS             2784..3761
#                     /gene="ND1"



my %store_current_features  = ();
my %store_current_proteins  = ();
my $products_found = 0;

for  my $i(1 .. $#split_entry)
	{
	$feature = $split_entry[$i];#	print "\n\nFEATURE1:($feature)\n";

	$feature =~ s/\s+gene\s+.+$//s;# 28/06/14: if the next feature is a protein coding gene (see example snippit above), 
					# you will have some details of it in this feature, so rm

	my $product = "NA";my $start = "NA";my $stop	= "NA";
#	print "\nFEATURE2:($feature)\n";#<1..>1468

	if($feature =~ /^\s*[complement]*\(*[><]*(\d+)\.\.[><]*(\d+)\)*/)
		{
		$start = $1;$stop	= $2;

		#FEATURE2:(4153..4219
 		#                    /gene="tRNA-Lys"
  		#                   /product="tRNA-Lys"
   		#                  /db_xref="GeneID:19351220")
		#feature number:10 product:NA start:4153 stop:4219

     		#rRNA            5415..6728
      		#               /gene="16S rRNA"
      		#               /product="l-rRNA"
      		#               /note="16S ribosomal RNA"
      		#               /db_xref="GeneID:14412113"

		#if($feature =~ /\s+\/gene="(.+)"/)
		if($feature =~ /\s+\/gene="([\w\d\s\-]+)"/)# 28/06/14: new regex for parsing names, 
			{$product = $1}				# previously was missing lots of 16S due to space or hyphen in each label
		if($feature =~ /\s+\/note="([\w\d\s\-]+)"/)
			{$product = $1}
		if($feature =~/\s+\/product="([\w\d\s\-]+)"/)
			{$product = $1}else{}

		if($product eq "NA")
			{
		#	print LOGFILE "\ncant read gene/note/product name of this feature.\n";
		#	print  "\ncant read gene/note/product name of this feature.\n($feature)\n";
			}else{
			$product = uc($product);$product =~ s/[^A-Z0-9]//g;
			$store_current_features{$start . "_" . $stop} = $product;$products_found++;
			}

		}else{
		if($feature =~ /^\s+\(*join\(\<*\d+|^order\(\d+/)
			{}else{
			#print LOGFILE "\ncannot read feature positions.\n";
			}
		}


#	print "feature number:$i product:$product start:$start stop:$stop\n\n";

#	if($product =~ /16S/){die ""}


	};# for each feature.


my @feature_keys = keys %store_current_features;
my $return_string = "";
foreach my $key(@feature_keys)
	{
				# $store_current_features{$start . "_" . $stop} = $product
	$return_string .= 	 $store_current_features{$key} . "__" . $key . "\t";
	}

if($products_found == 0)
	{
	return("ENTRY_NOT_FOUND");
	}else{
	$return_string =~ s/\t$//;
	return($return_string);
	}

}



##########################################################################################################################################
#
#
#
#
#
##########################################################################################################################################



sub extract_seqeunce_for_this_feature
{
my $current_feature = $_[0];my $current_sequence = $_[1];my $current_accession = $_[2];

		#:cytochrome oxidase subunit 1__1_658
if($current_feature =~ /(.+)__(\d+)_(\d+)/) 
	{
	my $name = $1; my $start = $2;my $end = $3;
#	print "name:$name start:$start end:$end\n";
	my $excision_length = $end-$start;

	if($end == $start)
		{return("NO_SEQUENCE")}

	if($end < $start)
		{# presumably this could indicate reverse complement? might be needed:
		#	$entry_retrieved = reverse($entry_retrieved);
		#	$entry_retrieved =~ tr/ACGTYRacgtyr/TGCARYtgcary/;
		die "\nerror 192. accession:$current_accession feature:$current_feature\n";
		}


	if($end > length($current_sequence)){die "\nerror.\n$current_feature\nend > length(current_sequence\n"}
	my $excised_seq = substr $current_sequence , $start-1, $excision_length+1;

#	print "\n\nfeature\n$current_feature\nroiginal sequence\n$current_sequence\nexcised sequence\n$excised_seq\n";
#	unless($current_sequence eq $excised_seq){print "\nOH LOOK\n"}

	return ($excised_seq);
	}

return("NO_SEQUENCE");
}



##########################################################################################################################################
#
#
#
#
#
##########################################################################################################################################





sub read_gene_synonyms
{

print "\nsub read_gene_synonyms\n";

# what to do with 'putative' , or 'similar to' , and misspellings: CYTOCHOMEOXIDSEII
# cad is probably short for cadherin, but what is relation to e-carherin, cadherin-1, cdh1?


%gene_synonyms = (
	'12S' 			=> '12S',
	'12SRDNA' 		=> '12S',
	'SRRNA' 		=> '12S',
'12SRIBOSMALRNA' 		=> '12S',
'SMALLSUBUNITRIBOSOMALRNA' 		=> '12S',

	'16S' 			=> '16S',
	'16SRDNA' 		=> '16S',
	'16SRIBOSOMALRNA' 	=> '16S',
	'16SRRNA' 		=> '16S',
'16SRIBOSMALRNA' 		=> '16S',

	'LRRNA' 		=> '16S',
	'LARGESUBUNITRIBOSOMALRNA'=> '16S',

	'18SRRNA' 		=> '18S',
	'18SRRNAGENE' 		=> '18S',

	'28S' 			=> '28S',
	'28SRDNA' 		=> '28S',
	'28SRIBOSOMALRNA' 	=> '28S',
	'28SRIBOSOMALRNAGENE' 	=> '28S',
	'28SRRNA' 		=> '28S',
	'28SRRNAGENE' 		=> '28S',

	'58SRRNA' 		=> '58SRRNA',
	'58SRRNAGENE' 		=> '58SRRNA',
	'5SRRNA' 		=> '58SRRNA',

	'CO1' 			=> 'CO1',
	'COI' 			=> 'CO1',
	'COX1' 			=> 'CO1',
	'COXI' 			=> 'CO1',
	'CYTOCHROMEOXIDASEI' 	=> 'CO1',

	'CO2' 			=> 'CO2',
	'COII' 			=> 'CO2',
	'COX2' 			=> 'CO2',
	'COXII' 		=> 'CO2',

	'CO3' 			=> 'CO3',
	'COIII' 		=> 'CO3',
	'COX3' 			=> 'CO3',
	'COXIII' 		=> 'CO3',

	'COB' 			=> 'CYTB',
	'CYTB' 			=> 'CYTB',
	'CYTOCHROMEB' 		=> 'CYTB',

	'CYTC' 			=> 'CYTC',
	'CYTOCHROMEC' 		=> 'CYTC',

	'EF1' 			=> 'EF1',
	'EF1A' 			=> 'EF1',
	'EF1AF1' 		=> 'EF1',
	'EF1AF2' 		=> 'EF1',
	'EF1ALPHA' 		=> 'EF1',
	'EF1ALPHA100E' 		=> 'EF1',
	'EF1ALPHAF2' 		=> 'EF1',
	'ELONGATIONFACTOR1ALPHA'=> 'EF1',
	'ELONGATIONFACTOR2LIKE' => 'EF1',

	'NAD1' 			=> 'NAD1',
	'ND1' 			=> 'NAD1',
	'NADH1' 		=> 'NAD1',
	'NDI' 			=> 'NAD1',

	'NAD2' 			=> 'NAD2',
	'ND2' 			=> 'NAD2',
	'NADH2' 		=> 'NAD2',

	'NAD3' 			=> 'NAD3',
	'ND3' 			=> 'NAD3',
	'NADH3' 		=> 'NAD3',

	'NAD4' 			=> 'NAD4',
	'NAD4L' 		=> 'NAD4',
	'ND4' 			=> 'NAD4',
	'ND4L' 			=> 'NAD4',
	'NADH4' 		=> 'NAD4',

	'NAD5' 			=> 'NAD5',
	'ND5' 			=> 'NAD5',
	'NADH5' 		=> 'NAD5',

	'NAD6' 			=> 'NAD6',
	'ND6' 			=> 'NAD6',
	'NADH6' 		=> 'NAD6',

# later additions:

'CYTOCHROMEOXIDASESUBUNIT1'  	=> 'CO1',
'CYTOCHROMEOXIDASESUBUNIT2'  	=> 'CO2',
'CYTOCHROMEOXIDASESUBUNIT3'  	=> 'CO3',

'CYTOCHROMEOXIDASESUBUNITI'  	=> 'CO1',
'CYTOCHROMEOXIDASESUBUNITII'  	=> 'CO2',
'CYTOCHROMEOXIDASESUBUNITIII'  	=> 'CO3',

'CYTOCHROMECOXIDASESUBUNITI'  	=> 'CO1',
'CYTOCHROMECOXIDASESUBUNITII'  	=> 'CO2',
'CYTOCHROMECOXIDASESUBUNITIII' 	=> 'CO3',

'CYTOCHROMECOXIDASESUBUNIT3' 	=> 'CO3',

'SIMILARTOCYTOCHROMEOXIDASESUBUNIT1'  	=> 'CO1',
'CYTOCHROMECOXIDASESUBUNIT1'  	=> 'CO1',
'SIMILARTOCYTOCHROMEOXIDASESUBUNITI'  	=> 'CO1',
'CYTOCHROMECOXIDASE1'  	=> 'CO1',
'CYTOCHROMECOXIDASEI'  	=> 'CO1',
'CYTOCHROMEOXIDASECSUBUNITI'  	=> 'CO1',
'CYTOCHROMEOXIDASE1'  	=> 'CO1',
'CYTOCHROMECOXIDASEISUBUNIT'  	=> 'CO1',

'CYTOCHROMECOXIDASESUBUNIT2'  	=> 'CO2',
'CYTOCHROMECOXIDASEII'  	=> 'CO2',

'NADHDEHYDROGENASE1'  => 'NAD1',
'NADHDEHYDROGENASE2'  => 'NAD2',
'NADHDEHYDROGENASE3'  => 'NAD3',
'NADHDEHYDROGENASE4'  => 'NAD4',
'NADHDEHYDROGENASE5'  => 'NAD5',

'NADHDEHYDROGENASESUBUNIT1'=> 'NAD1',
'NADHDEHYDROGENASESUBUNIT2'=> 'NAD2',
'NADHDEHYDROGENASESUBUNIT3'=> 'NAD3',
'NADHDEHYDROGENASESUBUNIT4'=> 'NAD4',
'NADHDEHYDROGENASESUBUNIT5'=> 'NAD5',

'NADHDEHYDROGENASESUBUNIT5ND5'=> 'NAD5',
'NADHDEHYDROGENASESUBUNITI'=> 'NAD1',

'NADHDEHYDROGENASESUBUNIT6'=> 'NAD6',

# nad4 and nad4L seem to be two seperate loci, adjacent
'NADHDEHYDROGENASESUBUNIT4L'=> 'NAD4L',

'NADHSUBUNIT1'=> 'NAD1',
'NADHSUBUNIT2'=> 'NAD2',
'NADHSUBUNIT3'=> 'NAD3',
'NADHSUBUNIT4'=> 'NAD4',
'NADHSUBUNIT4L'=> 'NAD4L',
'NADHSUBUNIT5'=> 'NAD5',
'NADHSUBUNIT6'=> 'NAD6',

'NADHDEHYDROGENSESUBUNIT1'=> 'NAD1',
'NADHDEHYDROGENSESUBUNIT2'=> 'NAD2',
'NADHDEHYDROGENSESUBUNIT3'=> 'NAD3',
'NADHDEHYDROGENSESUBUNIT4'=> 'NAD4',
'NADHDEHYDROGENSESUBUNIT4L'=> 'NAD4L',
'NADHDEHYDROGENSESUBUNIT5'=> 'NAD5',
'NADHDEHYDROGENSESUBUNIT6'=> 'NAD6',

'12SRIBOSOMALRNA' => '12S',

'CYTOCHROMEOXIDASEII'  	=> 'CO2',
'CYTOCHROMEOXIDASEIII' 	=> 'CO3',

'WINGLESS'=>'WNT',
'WINGLESSPROTEIN'=>'WNT',

'CAD'=>'CAD',
'CADHERIN'=>'CAD',

'28SLARGESUBUNITRIBOSOMALRNA' => '28S',
'12SSMALLSUBUNITRIBOSOMALRNA'=> '12S',
'16SLARGESUBUNITRIBOSOMALRNA'=> '16S',
'18SRIBOSOMALRNA' => '18S',
'18SRIBOSOMALRNAGENE'=> '18S',
'18SSMALLSUBUNITRIBOSOMALRNA'=> '18S',

'INTERNALTRANSCRIBEDSPACER1ITS1' => 'ITS1',

'INTERNALTRANSCRIBEDSPACER2ITS2' => 'ITS2',

'HISTONEH3' => 'H3',
'HISTONE3'=> 'H3',

'TRANSLATIONELONGATIONFACTOR1ALPHA'=> 'EF1',
'ELONGATIONFACTOR1A'=> 'EF1',


'ATPASE6' 	=> 'ATP6',
'ATPASESUBUNIT6' 	=> 'ATP6',
'ATPSYNTHASEFOSUBUNIT6' => 'ATP6',
'ATPSYNTHASEF0SUBUNIT6' => 'ATP6',
'ATPSYNTHASESUBUNIT6' => 'ATP6',

'ATPASESUBUNIT8' 	=> 'ATP8',
'ATPSYNTHASEF0SUBUNIT8' => 'ATP8',
'ATPSYNTHASEFOSUBUNIT8' => 'ATP8',


'CADHERINLIKEPROTEIN'=>'CAD',


	);


my @gene_synonym_keys = keys %gene_synonyms;@gene_synonym_keys = sort @gene_synonym_keys;

foreach my $key(@gene_synonym_keys)
{
my $belongs_to = $gene_synonyms{$key};
#print "$key standardized to $belongs_to\n";

}



};


##############################################################################################################
#
#
#
#
#
##############################################################################################################




sub parse_genbank_flatfiles2
{
#open(LOG2, ">>DodgyDnaSeqsFound") or die "cantopen\n";
#open(MAIN_OUT, ">inv.fas") or die "cantopen\n";
#open(OUTKEY, ">accession_key.$month.$accession_keyfile_ID") || die "\nerror 404\n";
#print OUTKEY "accession\ttaxid\torder\ttaxstring\tspecies_id\tproduct\tcountry\tcoordinates\n";
#close(OUTKEY);
#print LOG2 "\n\n\nRUNNING:$date\n";

for my $zipped_infile(@zipped_genbank_flatfiles)
	{

#	##$zipped_infile = "gbinv19.seq.gz";

	my $command = "gunzip $zipped_infile";
	print "command:$command\n";
	system($command);
	$zipped_infile =~ s/\.gz//;
	$zipped_infileTEST = $zipped_infile;

	####################################
	parse_genbank_flatfileNEW($zipped_infile);#
	####################################

	my $command = "gzip $zipped_infile";
	print "command:$command\n";
	system($command);

	my @count_samples = keys %gene_counts;@count_samples = sort @count_samples;
	open(OUT_GENE_COUNTS, ">resultsgenecouts") || die "\nerror\n";
	foreach my $gene(@count_samples)
		{
		print OUT_GENE_COUNTS "$gene $gene_counts{$gene}\n";
		}
	close(OUT_GENE_COUNTS);

	print "\nfrom all files, $countnumberprinted were printed out of $counttotal\n";

	}

#close(LOG2);
#close(MAIN_OUT);

#my @speciesremoved = keys %sp_rm;
#print scalar @speciesremoved , " removed\n@speciesremoved\n";
#print "\n$countnumberprinted were printed out of $counttotal\n";

}



##############################################################################################################
#
#
#
#
#
##############################################################################################################

sub parse_genbank_flatfileNEW
{
my $inputfile = shift;
open(IN2, $inputfile) || die "error 702. cant open infile:$inputfile\n";
print "\nsub parse_genbank_flatfileNEW\nopened $inputfile\n";

my $fileasstring = "";
while(my $line = <IN2>)
	{$fileasstring .= $line}
close(IN2);
my @file_as_array = split(/\n\/\/\n/ , $fileasstring);

my $count_entries_this_flatfile=0;

my $entrycount = 0;
foreach my $entry(@file_as_array)
	{
	$entrycount++;
	#print "entry no $entrycount. length:" , length($entry) , "\n";
	$count_entries_this_flatfile++;if($count_entries_this_flatfile =~ /0000$/){print "$inputfile $count_entries_this_flatfile of $#file_as_array\n"}
	#print "\n$inputfile. from all files, $countnumberprinted were printed out of $counttotal\n";
	#if ($entry =~ /Z93710/i){print $entry}


	#########################
	parse_this_entry($entry)#	
	#########################


	}


unless($done>= 1)
	{print "\nerror 736.\n"; die};


}




##############################################################################################################
#
#
#
#
#
##############################################################################################################


sub parse_this_entry
{
my $current_entry = shift;
my $current_entry_copy = $current_entry;
$counttotal++;


#######		ACCESSION

my $accession = "UnknownAccesion";
if($current_entry =~ /ACCESSION\s+(\S+)\n/)
	{
	$accession = $1# single values for accession found
	}else{
	# multiple values for accession, first one usually matches with that found on description
	#DEFINITION  Caenorhabditis elegans DNA for 32-kDa galectin, complete cds.
	#ACCESSION   AB000802 D85885 D85886 D85887
	if($current_entry =~ /ACCESSION\s+(\S+)\s(.+)\n/)
		{$accession = $1
		}else{
		$countmissingaccessions++;
		#die "\n\nERROR, no accession\n\n$current_entry\n"
		}
	}


my $gi_number = "UnknownGI";
if($current_entry =~ /VERSION\s+\S+\s+GI\:(\d+)\n/)
	{
	# single values for accession found
	$gi_number = $1
	};


if($accession_or_gi == 2)
	{$accession = $gi_number};

	if(length($store_all_accessions{$accession})>= 1)
		{
		$done++;

		if($verbose ==1) {print "\n$accession, in list of insect entries.\n"};

		my $returned = "";
		if($parse_protein == 1)
			{
			#print "\ntest\n";
			$returned = fetch_proteins($current_entry);
			}else{

			##########################
			$returned = fetch_id($current_entry);#
			##########################

			};
		#print "test";
		if($verbose ==1) {
		print "returned:$returned\n";};

		unless($returned eq "ENTRY_NOT_FOUND")
			{
			my @each_feature = split /\t/, $returned;
			#print "each_feature:@each_feature\n";die;
			foreach my $feature(@each_feature)
				{
				#print "\nfeature:$feature\n";#die;
				my $excised_sequence = "NO_SEQUENCE";
				if($parse_protein == 1)
					{
					if($feature =~ /^.+__(.+)/)#/(.+)__(\d+)_(\d+)/ 
						{
						$excised_sequence = $1;
						}else{die "\n\nerror 804.\n\n"}
					
					}else{
					#####################################################
					$excised_sequence = extract_seqeunce_for_this_feature($feature , $store_all_accessions{$accession} , $accession);#
					#####################################################
					};


				unless($excised_sequence eq "NO_SEQUENCE")
					{
					my $name = "gene_unknown";
					if($feature =~ /^(.+)__/)#/(.+)__(\d+)_(\d+)/ 
						{
						$name = $1;#print "\tname:$name\n";
						}
					if(exists($gene_synonyms{$name}))
						{
					#	print "\tfound name ($name), changing it to standardized name ($gene_synonyms{$name})\n";
						$name = $gene_synonyms{$name};
						}else{
						#print LOGFILE "name $name not found in curated list\n";
						if($verbose ==1) {print "name $name not found in curated list\n"};
					}
					if(length($excised_sequence)>= $excised_seq_length_limit)
						{
						my $complete_fas_id = $full_id_for_accessions{$accession};
						if(exists($record_ids_printed{$complete_fas_id . $name}))
							{
							print "duplicated ID ($complete_fas_id for gene $name) come across $record_ids_printed{$complete_fas_id . $name} time before\n";
							$record_ids_printed{$complete_fas_id . $name}++;
							}else{
							open(OUT , ">>nameparsed.$name") || die "\nerror 74\n";	
							print OUT ">$complete_fas_id\n$excised_sequence\n";
							close(OUT);
							$gene_counts{$name}++;
							$record_ids_printed{$complete_fas_id . $name}++;
							}
						};
					}
				}
			}




		}else{
		#print "$accession, not in list.\n"
		};




}


##############################################################################################################
#
#
#
#
#
##############################################################################################################


sub fetch_proteins
{

my $gb_entry = shift;

#print "entry ($gb_entry)\n";
my @split_entry = split /\s+CDS\s+|\s+tRNA\s+|\s+rRNA\s+|\s+misc_feature\s+/ , $gb_entry;# bugfix 28/06/14

#my %store_current_features  = ();
my %store_current_proteins  = ();
my $products_found = 0;

for  my $i(1 .. $#split_entry)
	{
	$feature = $split_entry[$i];#	print "\n\nFEATURE1:($feature)\n";
	$feature =~ s/\s+gene\s+.+$//s;
	my $product = "NA";my $start = "NA";my $stop	= "NA";
	#print "\nFEATURE2:($feature)\n";#<1..>1468

	if($feature =~ /^\s*[complement]*\(*[><]*(\d+)\.\.[><]*(\d+)\)*/)
		{
		$start = $1;$stop	= $2;

		if($feature =~ /\s+\/gene="([\w\d\s\-]+)"/)
			{$product = $1}			
		if($feature =~ /\s+\/note="([\w\d\s\-]+)"/)
			{$product = $1}
		if($feature =~/\s+\/product="([\w\d\s\-]+)"/)
			{$product = $1}else{}

		if($product eq "NA")
			{
		#	print LOGFILE "\ncant read gene/note/product name of this feature.\n";
		#	print  "\ncant read gene/note/product name of this feature.\n($feature)\n";
			}else{
			$product = uc($product);$product =~ s/[^A-Z0-9]//g;

			my $protein_sequence = "";
			if($feature =~/\s+\/translation="(.+)"/s)# reads (AA) string over multi-lines
				{
				$protein_sequence = $1;$protein_sequence =~ s/[\s\n\r]//g;
				$protein_sequence =~ s/^([A-Z]+)\".+/$1/;
				unless($protein_sequence =~ /^[A-Z]+$/)
					{die "\nerror 944: strange characters in protein sequence:($protein_sequence)\n"}
				#print "\tprotein_sequence:($protein_sequence)\n";
				$store_current_proteins{$protein_sequence} = $product;$products_found++;
				};

			}

		};

	};# for each feature.


my $return_string = "";
my @prot_keys = keys %store_current_proteins;
foreach my $key(@prot_keys)
	{
	$return_string .= $store_current_proteins{$key} . "__" . $key . "\t";
	#print "key:($key)\nstore_current_features{key}:(" , $store_current_proteins{$key} , ")\n\n";
	};

if($products_found == 0)
	{
	return("ENTRY_NOT_FOUND");
	}else{
	$return_string =~ s/\t$//;
	return($return_string);
	}




};#sub fetch_proteins


##############################################################################################################
#
#
#
#
#
##############################################################################################################











