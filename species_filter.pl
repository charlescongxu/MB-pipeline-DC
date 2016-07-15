#!/usr/bin/perl
# 
# 
# species_filter.pl
# Douglas Chesters. 
# Chinese Academy of Sciences.
# 
# 
# 
# change log 
# 01Aug2014: for each species, select the most representative (as opposed to the longest)
# 04Aug2014: bugfix for printing identified / non-identified
# 05Aug2014: option to filter by least ambiguous bases.
# 26Dec2014: insignificant changes
# 16Jan2015: deletes subspecies names
# 24Jan2015: discard entries where cant parse binomial.
#		removes substrings for verbose species IDs
# 31Jan2015: stop printing 'nr' species 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
##################################################################################################


my $version = 1.006;


my $input 			= $ARGV[0];
my $format 			= $ARGV[1]; #if($format == 2)species_filtering_tobycoded else species_filtering();


$identified_species_only 	= 1;
$filter_by_MRS			= 0; 	# 0 = take longest sequence for each species. 
					# 1 = take the most representative sequence for each species.
					# 2 = take the sequence with least ambiguous data (N's etc). these make alignment more difficult. 

$MRS_script 			= "~/usr_scripts/most_representative_seq.pl";


$trim_off_accession 		= 0;	# must be 0 or 1

$print_screen			= 1;	# 0 = quiet

############################################################################################################



$date = localtime time;

open(LOG, ">>species_filter_LOG");
print LOG "\n$date input:$input format:$format\n";
system("rm $input.ID_filtered");


unless($print_screen == 0)
{
print "
\n ***** species_filter.pl (v $version)  *****  
";
if($filter_by_MRS == 0){print "filter_by_MRS == 0. take longest sequence for each species.\n"};
if($filter_by_MRS == 1)
	{
	print "filter_by_MRS == 1. take the most representative sequence for each species.\n";
	print "MRS_script:$MRS_script\n";
	};
if($filter_by_MRS == 2){print "filter_by_MRS == 2. take the sequence with least ambiguous data (N's etc). these make alignment more difficult.\n"};
print "trim_off_accession:$trim_off_accession\n";
print "identified_species_only:$identified_species_only\n\n";

};




if($format == 2)
	{
	species_filtering_tobycoded();
	}else{
	species_filtering();
	}

close(LOG);


unless($print_screen == 0){print "\nFIN.\n"};


exit;

##############################################################################################################
#
#
#
#
#
##############################################################################################################




sub species_filtering
	{

unless($print_screen == 0){print "\nsub species_filtering\n"};

open (IN_FILTER, $input) || die "cant open $input\n";

my $file_as_string = "";

while (my $line = <IN_FILTER>)
	{$file_as_string .= $line}
close(IN_FILTER);

my @all_lines = split />/, $file_as_string;

unless($print_screen == 0){print scalar @all_lines , " seqs in your input file.\n"};

for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];

#	if($line =~ /^([^\_]+[_\.][^\_]+)(.*)/)
	if($line =~ /(.+)/)
		{
		my $description = $1;$description =~ s/[\.\-]//g;#print "\ndescription:$description\n";
		my $speciesid = "NA";my $rest = "NA";
		my $discard_entry = 0;

		if($description =~ /^([A-Z][a-z]+_[a-z]+)_([a-z]+)(_\d+)$/)#Apis_mellifera_intermissa_693583962
			{
			$speciesid = $1; my $subspecies = $2; my $gi=$3;$rest = $gi;



			#print "speciesid:$speciesid subspecies:$subspecies gi:$gi\n";
			}elsif($description =~ /^([A-Z][a-z]+_[a-z]+)(_.+)$/)#Helicoverpa_punctigera_586947476
				{
				$speciesid = $1;$rest = $2;

				if($rest =~ s/.+_.+_.+_(.+)$/$1/)
					{
					print "\nwarning. long ID ($description) doesnt look like accession, " , 
					"will cut some of the string out: ($rest) \n"
					};
				#>Colletes_fasciatus_s_l_NML_2007_EF028541
				#>Colletes_clypearis_group_sp_ZQN_2013_KC469666

				}else{
				print "\nWARNING. cant parse binomial from ID:($description) ... discarding this entry.\n";
				$discard_entry = 1;$discarded_entries++;
				};



		
		#print "$speciesid rest:$rest\n";
		$line =~ s/^.+\n//;
		$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/\012\015?|\015\012?//g;

		unless($discard_entry == 1)
		{

		if($filter_by_MRS == 1)
			{

			$sequences{$speciesid} .= ">$speciesid$rest\n$line\n";

			}else{



			if($filter_by_MRS == 0)	# take longest
				{
				if(exists($sequences{$speciesid}))
					{
					my $len_existing = length($sequences{$speciesid});
					my $len_new = length($line);	#	print "\texisting:$len_existing new:$len_new\n";
					if($len_new >= $len_existing)
						{$sequences{$speciesid} = $line;$therest{$speciesid} = $rest;#print "\treplacing\n"
						}
					
					}else{
					$sequences{$speciesid} = $line;
					$therest{$speciesid} = $rest;
					}
				}


			if($filter_by_MRS == 2)	# take that with least ambig data
				{
				if(exists($sequences{$speciesid}))
					{
					print LOG "speciesid:$speciesid observed already\n";

					my $existing_seq = $sequences{$speciesid};my $new_seq = $line;
					#print "existing_seq:$existing_seq\nnew_seq:$new_seq\n";
					my $count_ambig_existing = 0;my $count_ambig_new = 0;
					while($existing_seq =~ s/[nyr]//i){$count_ambig_existing++};
					while($new_seq =~ s/[nyr]//i){$count_ambig_new++};
					print LOG "\tcount_ambig_existing:$count_ambig_existing count_ambig_new:$count_ambig_new\n";

					if($count_ambig_new < $count_ambig_existing)
						{$sequences{$speciesid} = $line;$therest{$speciesid} = $rest;#print "\tless ambiguous data. replacing\n"
						}

					}else{
					print LOG "speciesid:$speciesid is novel\n";
					$sequences{$speciesid} = $line;
					$therest{$speciesid} = $rest;
					
					};				
				}


			}

			};

		}else{
		print "\n\nwarning. unexpected ID.\n"
		}
	}

my @unique_species = keys %sequences;
@unique_species = sort @unique_species;
open(OUT, ">$input.ID_filtered") || die "\n\nerror\n";
close(OUT);
my $printed_to_output =0;


foreach my $sp(@unique_species)
	{


	if($filter_by_MRS == 1)
		{
		

		my $theMRS = find_MRS($sp);

			open(OUT, ">>$input.ID_filtered") || die "\n\nerror\n";
			if($trim_off_accession == 1){$theMRS =~ s/(>[A-Z][a-z]+_[a-z]+).+/$1/}
			print OUT "$theMRS";
			close OUT;
$printed_to_output++;

		}else{

		my $print_current = "";
		if($identified_species_only == 1)
			{
			if($sp =~ /^[A-Z][a-z]+[_][a-z][a-z]+$/ )#>Ceratina_n_subgen_sp_SMR_2010_GU321539
				{
				# >Euglossa_cf_variabilis_SR417_AY920314

				if( $sp =~ /^[A-Z][a-z]+[_]sp$/ || $sp =~ /^[A-Z][a-z]+[_]aff$/ || 
					$sp =~ /^[A-Z][a-z]+[_]cf$/ || $sp =~ /^[A-Z][a-z]+[_]nr$/)
					{
					$print_current = 0;
					}else{
					$print_current = 1;
		
					}
				
				#if($sp =~ /Braunsapis_nr/){print "sp:$sp\n";die "\n\n"}

				}else{
				my $print_current = 0;
				}  
			}else{
			$print_current =1;
			}

		if($print_current == 1)
			{
		#	print "sp:$sp\n";
			
			open(OUT, ">>$input.ID_filtered") || die "\n\nerror\n";
			print OUT ">$sp" , "";
			if(length($therest{$sp})>=2 && $trim_off_accession == 0)
				{print OUT $therest{$sp}}
			print OUT "\n" , $sequences{$sp}  , "\n";#	print ">$sp\n" , $sequences{$sp} , $therest{$sp} , "\n";
			close OUT;
			$printed_to_output++;

			};

		}

	}#foreach my $sp(@unique_species)



unless($print_screen == 0)
	{
print "
discarded_entries:$discarded_entries
$printed_to_output printed_to_output\n"
};


}



############################################################
#
#
#
#
#
############################################################


sub species_filtering_tobycoded
	{

print "\nspecies_filtering, assuming taxstring underscore accession\n";


open (IN_FILTER, $input) || die "cant open $input\n";

my $file_as_string = "";

while (my $line = <IN_FILTER>)
	{$file_as_string .= $line}
close(IN_FILTER);

my @all_lines = split />/, $file_as_string;

print scalar @all_lines , " seqs in file\n";

for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];

	if($line =~ /^(.+)[_]([^\_\n\r]+)/)
		{
		my $speciesid = $1;my $rest = $2;$rest =~ s/[\n\r \s]//g;unless($rest =~ /^[A-Z]+[0-9]+$/){die "\n\nwierd accession in:$rest\n"}
	#	print "$speciesid rest:$rest\n";
		$line =~ s/^.+\n//;
		$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/\012\015?|\015\012?//g;

		if(exists($sequences{$speciesid}))
			{
			my $len_existing = length($sequences{$speciesid});
			my $len_new = length($line);
		#	print "\texisting:$len_existing new:$len_new\n";
			$number_replacments{$speciesid}++;

			if($len_new >= $len_existing){$sequences{$speciesid} = $line;$therest{$speciesid} = $rest;
				#print "\treplacing\n"
				}
			}else{
			$sequences{$speciesid} = $line;
			$therest{$speciesid} = $rest;
			}

		}
	}

my @unique_species = keys %sequences;
@unique_species = sort @unique_species;
open(OUT, ">$input.ID_filtered") || die "\n\nerror\n";

my $countprinted=0;

foreach my $sp(@unique_species)
	{

	my $print_current = 1;
	if($identified_species_only == 1)
		{
		unless($sp =~ /^[A-Z][a-z]+[_][a-z]+$/){$print_current = 0}# only works with format 1
		}else{
		$print_current =1;
		}


	if($print_current == 1)
		{
		print OUT ">$sp" , "";
	$countprinted++;
		if(length($therest{$sp})>=2)
			{
			print OUT $therest{$sp}; 
			}

		print OUT "\n" , $sequences{$sp}  , "\n";
	#	print ">$sp\n" , $sequences{$sp} , $therest{$sp} , "\n";
		}

if($number_replacments{$sp}>=100){print LOG "species ID:$sp occured:$number_replacments{$sp} times in file\n"}

unless($sp =~ /^[^7]+7[a-z]+$/){$unidentified++}
	}


close(OUT);

print scalar @all_lines , " seqs in file, filtered to:$countprinted\n";
print LOG scalar @all_lines , " seqs in file, filtered to:$countprinted. unidentified:$unidentified\n";



}



############################################################
#
#
#
#
#
############################################################


sub find_MRS
{
my $species = shift;

my $data = $sequences{$species};


my @all_lines = split />/, $data;splice @all_lines , 0,1;
unless($print_screen == 0){print "species:$species, seqs:" , scalar @all_lines , "\n"};

my $MRS = "NA";

if(scalar @all_lines == 1)
	{
	
	$MRS = $data;
#	print "only 1 sequence:$MRS\n";

	}elsif(scalar @all_lines == 2)
	{
	# 2 sequences for current species, pick at random.

	$MRS = ">$all_lines[0]";
#	print "only 2 sequence, selecting 1 randomly:$MRS\n";

#	die;
	}else{

	# 3 or more sequences, find MRS

	open(OUT77, ">all_seqs_current_sp") || die "\nerror 286\n";
	print OUT77 "$data";
	close OUT77;

	my $command = "perl $MRS_script -i all_seqs_current_sp -s 120 -b blastn-2.2.28+-64bit -a 1200 -p 80";
	# -i followed by name of fasta file of DNA sequences
	# -s followed by number of sequences to sample from the input file,
	#    from which the MRS will be found
	# -b followed by the name of the command for running blast on your system
	#    e.g. blastn or /home/plantID/ncbi-blast-2.2.27+/bin/blastn or blastn-2.2.28+-64bit
	# -a followed by number, which is sequence length below which hits are ignored
	# -p follows by number, percent identity cutoff

	system("rm all_seqs_current_sp.MRS");
	system($command);

	my $test_MRS = `cat all_seqs_current_sp.MRS`;

	print "test_MRS:($test_MRS)\n";
	my @test_splitfile = split />/, $test_MRS;splice @test_splitfile , 0,1;

	if(scalar @test_splitfile == 1)
		{
		$MRS = 	$test_MRS;		
		}else{
		print "\nerror:MRS output contains more than one sequence. or none. just taking first by random.\n";
		$MRS = ">$all_lines[0]";

		}
	#die;
	}

return($MRS);

}






############################################################
#
#
#
#
#
############################################################









