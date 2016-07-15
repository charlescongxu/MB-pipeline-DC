#!/usr/bin/perl

###################################################################################################################################
#
#
#
# 	print_taxa_overlapping_newick_and_fasta.pl, 
#		Perl script that parses the taxa overlapping between a phylogeny and a file of DNA sequences 
#
#
#    	Copyright (C) 2014  Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
#
###################################################################################################################################
#
#
#
#	Takes a file with a phylogeny and a file with DNA sequences, 
#	reads the taxa on each, then prints a new phylogeny and new sequence file of only the overlapping taxa
#	Currently only works for terminals/fasta IDs in strict binomial format
# 	Script is entirely self-contained, does not load any modules for anything,
#	so you can modify any aspect for whatever your purpose.
#
#	To run:
#	perl print_taxa_overlapping_newick_and_fasta.pl [newick_file_name] [fasta_file_name]
#
#
#
#
#	change log:
# 	03JUL2014: first release version 
#	24JUL2014: option to include taxa in the fasta file that are not in the tree.
#		these are included since a raxml constraint analysis can add them
#	05AUG2014: rm boot values from input reference tree. this was stopping correct functioning.
#
#
	$script_version = 1.03;
#
#
#
#
#
###################################################################################################################################



$include_nonoverlapping_taxa_in_fasta_file = 1;
$verbose = 0;


$treefile 	= $ARGV[0];
$fastafile 	= $ARGV[1];



# check input files:

#########
check();#
#########



print "\n\n\n\tprint_taxa_overlapping_newick_and_fasta.pl (version:$script_version)
\tinput files, tree:$treefile sequences:$fastafile\n\n";





# some global variables:

%nodes;
$root_node;
%parent;
$count_terminals_from_current_node = 0;




# read the tree into memory

######################
read_file($treefile);#
######################


my @terminal_array = keys %terminals;@terminal_array = sort @terminal_array;
print scalar @terminal_array , " terminals parsed\n";


# read the sequence file and record entries that overlap with those found in the tree earlier

################
fasta_reader();#
################



my $notfound = 0;
my $found = 0;

foreach my $terminal(@terminal_array)
	{
	unless($terminal =~ /^[A-Z][a-z]+_[a-z]+$/)
		{#print "non binomial:$terminal\n"
		}

	if($terminals{$terminal} == 2)
		{
		$found++; if($verbose ==1){	print "tree terminal:$terminal is found in sequence DB\n"};
		}else{
		$notfound++; if($verbose ==1){	print "tree terminal:$terminal is NOT found in sequence DB\n"};	
		}
	}


print "\nOverlapping tree and DB:$found. 
In parsed tree but not found in DB:$notfound\n\n";




# There are probably many ways to prune a tree, of which this is just one.....
# The tree (even when unrooted) has a natural dimension, root to tips.
# Reading the tree earlier is easy from tip to root due to size of the strings,
# but once the tree is in memory, it seems more natural to travel root to tip (only one starting point).
# So traverse in this way, but ignore the members that need pruning 
# (by not printing them to the new newick string).
# However you dont want to wait until you actually reach these nodes,
# since what if the sister node also requires pruning, 
# that means you have just processed a node (the parent) that is now also a non-required node,
# and you have to go back up and thats complicated.
# So, first check each node to see if its a 'dead-end', 
# i.e. all decendends are to be pruned, 
# that way, you can then traverse the tree and know all nodes to be pruned, not just all terminals.


################################
find_deadend_nodes($root_node);#
################################


# a new newick string will be developed in which non-overlapping members are ommited, start of this string:
$new_newick_string = "($root_node);";

##################################
build_tree_structure($root_node);#
##################################


# some minor processing of the new newick tree, then print to file:
$new_newick_string =~ s/^\((.+)\)(\;)/$1$2/;
open(OUT, ">$treefile.overlapping") || die "\nerror, cant open output file ($treefile.overlapping). \nquitting.\n";
print OUT "$new_newick_string\n";
close OUT;


# then print fasta entries which are overlapping, this is a lot simpler

my @overlapping_tax = keys %store_seqs_for_overlapping_species;@overlapping_tax  =sort @overlapping_tax;

open(OUT, ">$fastafile.overlapping") || die "\nerror 100.\n";
foreach my $tax(@overlapping_tax)
	{print OUT ">$tax\n$store_seqs_for_overlapping_species{$tax}\n"}

close OUT;



print "\nfin.\n";
exit;


# end of script, subroutines are below


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_file
{
my $file = shift;

open(IN , $file) || die "\nerror 24.\n";
my $tree;my $found_tree=0;

while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /(.+\(.+)\;/)
		{
		$tree  = $1;
		$found_tree =1;
		}
	}

close(IN);

if($found_tree == 1)
	{

	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	$tree =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/g; 

	# remove regular branchlengths: 0.02048
	$tree =~ s/\:\-*\d+\.\d+//g; 

	# remove 0 length branchlengths
	$tree =~ s/\:\d+//g;

	# remove these things
	$tree =~ s/[\[\]]//g; 

	# tree from 'delayed rise of modern day mammals':
	# don':5.0,Perameles_bougainville:12.4,(Perameles_gunnii:7.7,Perameles_nasuta:7.7)'2077_Perameles*':4.7)'2074_Peramelidae*':6.1)'2070':17.7,Macrotis
	# remove these node labels
	$tree =~ s/\)\'[^\']+\'/)/g; 

	# remove boot support values (05 aug 2014). not sure why this wasnt here before
	while($tree =~ /\)\d+[\)\,]/){$tree =~ s/(\))\d+([\)\,])/$1$2/g}; 


print "\nNOTE:
  tree file has been read. comments, branchlengths and support values have been removed.\n\n";

#print "tree:$tree\n";die;


	######################
	parse_newick ($tree);#
	######################


	}else{die "\nerror 229.\n"}

}




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub parse_newick
{
my $newick_string = shift;

#if($job == 1){$tree_copy = $newick_string}


#print "newick_string:$newick_string\n";
#print OUT "$newick_string\n";

$interal_node =0;


# tree is parsed inside-out, this means we are dealing with short simpler strings
# parse one node at a time, by looking for instances of ( ....... )


while ($newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/)
	{

	# found a pair of parentheses (with no parentheses inside).
	# this string is taken, and replaced with a single node (INTERNAL_NODE_\d)
	# node is numbered (\d in variable $interal_node), 
	# the number being incremented at the end of the loop

	my $node = $1;my $nodeID = "INTERNAL_NODE_$interal_node";

	# find what is at the current node (decendents), by splitting at the commas
	
	my @child_nodes = split /\,/ , $node;
	$child_counts{$nodeID} = $#child_nodes;

	for $i(0 .. $#child_nodes)
		{
		my $bl = "NA"; 	# branch length stuff, currently not relevent
		if($child_nodes[$i] =~ /\:(.+)/)
			{$bl = $1
			}else{#die "error no banchlength"
			};
		$child_nodes[$i] =~ s/\:(.+)//;
		#$node_branchlengths{$child_nodes[$i]} = $bl;

		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }

		$nodes{$nodeID}{$i} = $child_nodes[$i];
		$parent{$child_nodes[$i]} = $nodeID;

	#	print "node:$child_nodes[$i]\tbl:$bl\n";

		# and record whether the current node is a terminal one

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/){$terminals{$child_nodes[$i]} =1}

		}
#	print "node:$interal_node\n\tchild nodes:@child_nodes\n";


	$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.

	$interal_node++;
	}


#print "newick string has been read\n";

}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub build_tree_structure
{
my $current_node = shift;
my $child_count_current_node = $child_counts{$current_node};
#print "\n\nNEW NODE:$current_node number of children:$child_count_current_node\n";


########################################################################################################################


# start string that will replace current node:
my $current_node_new_string = "";

for my $j(0 .. $child_count_current_node)
	{

	# retrieve decendends of the current node, go through each and check if it will be incorporated,

	my $child_node = $nodes{$current_node}{$j};
	if($terminals{$child_node} == 1)		# this hash is filled as: $terminals{$binomial} = 2;
		{					# current entry of fasta file ($speciesid) is from a species in the tree ($binomial)\n";
							# otherwise, it would be 1

		#print " ... is NOT recorded as overlapping\n";
		}else{
		#print " ... might be overlapping ... ";

		if($deadendnode{$child_node}==1)
			{
			#print "but is a deadend node\n";
			}else{
			#print "appending to node replacemnt string\n";

			# current child node is overlapping, so put into current string:
			$current_node_new_string .= "$child_node,";	
			}

		};
	}


########################################################################################################################


# remove an unwanted comma at end of line
$current_node_new_string =~ s/\,$//;

if($current_node_new_string =~ /\,/)
	{
	# current node has more than one member, will need to be encased in parentheses
	$current_node_new_string = "(" . $current_node_new_string . ")"
	};

#print "\tcurrent node:$current_node becomes:$current_node_new_string\n";

if($current_node_new_string eq "")
	{
	print "NOTHING\n";die "\nerror 409\n";
	}else{

	# update newick string at the current node, replace old node with the decendents

	$new_newick_string =~ s/([\(\)\,])$current_node([\(\)\,])/$1$current_node_new_string$2/;

	};


########################################################################################################################


for my $j(0 .. $child_count_current_node)
	{

	# current node is now processed. 
	# go though the list of decendents again, and if appropriate, traverse to next node

	my $child_node = $nodes{$current_node}{$j};#print "\tnode:$current_node child:$child_node\n";
	if($child_node =~ /^INTERNAL_NODE_\d+/)
		{
		if($deadendnode{$child_node}==1)
			{
			#print "\tchild_node:$child_node has no terminals to be printed\n";
			}else{
			#print "\tchild_node:$child_node, recursing\n";

			###################################
			build_tree_structure($child_node);# current child node is not terminal, and not a dead-end, so keep going tip-ward
			###################################

			}
		}else{
		#print "\tchild_node:$child_node, end of the line\n";
		}
	}




}




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub fasta_reader
{


my $file_as_string = "";
open(IN_FILTER , $fastafile) || die "\nerror. cant open\n";

print "sub fasta reader, reading:$fastafile\n";

while (my $line = <IN_FILTER>)	
	{$file_as_string .= $line};
close(IN_FILTER);


# this is a simple fasta reader which reads all the file into memory,
# therefor not suitable for very very large sequence files.

my @all_lines = split />/, $file_as_string;
my $number_of_binomials_read=0;

for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];
	if($line =~ /^(.+)/)
		{
		my $speciesid = $1;	#print "$speciesid\n";

		$line =~ s/^.+\n//;
		$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;
		$line =~ s/[\s\t]//g;

		if($speciesid =~ /^([A-Z][a-z]+_[a-z]+)/)
			{
			my $binomial = $1;$number_of_binomials_read++;
			if(exists($terminals{$binomial}))
				{
				#print "current entry of fasta file ($speciesid) is from a species in the tree ($binomial)\n";

				# current entry is overlapping, record the ID and store the DNA sequence:
				$terminals{$binomial} = 2;
				$store_seqs_for_overlapping_species{$binomial} = $line;

				}else{

				if($include_nonoverlapping_taxa_in_fasta_file == 1)
					{
					$store_seqs_for_overlapping_species{$binomial} = $line;#######################
					}
			#	print "binomial ($binomial) wasnt recorded from tree\n";
				}

			}
		}#if($line =~ /^(.+)/)
	}#for my $each_line(1 .. $#all_lines)


print "

NOTE:
  currently this script is just for use with identifiers in strict binomial format:Genus_species
  if you need a more general implementation, send me you input files and i will modify the script accordingly\n\n";


print scalar @all_lines , " seqs in file\n";

unless($number_of_binomials_read >= 2)
	{
	print "\nerror parsing binomial fasta IDs. expecting:[A-Z][a-z]+_[a-z]+\nquitting.\n\n"
	}



};




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub get_terminals_from_this_node
{

my $current_node_for_terminals 	= shift;
my $child_count2 		= $child_counts{$current_node_for_terminals};
my $current_node_new_string 	= "";

for my $j(0 .. $child_count2)
	{
	my $child_node2 = $nodes{$current_node_for_terminals}{$j};#print "\tnode:$current_node child:$child_node\n";
	if($child_node2 =~ /^INTERNAL_NODE_\d+/)
		{

		############################################
		get_terminals_from_this_node($child_node2);#
		############################################

		}else{

		# traverse through to a terminal node and check if it is overlapping

		if($terminals{$child_node2} == 2)
			{
			$count_terminals_from_current_node++;
		#	print "terminal:$child_node2 overlapping:1 ";
			}else{
		#	print "terminal:$child_node2 overlapping:0 ";
			};
		};
	}

return();

}



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub find_deadend_nodes
{

my $current_node = shift;
my $child_count_current_node = $child_counts{$current_node};

#print "\nNODE1:$current_node number of children:$child_count_current_node\n";


# for current node check all terminals that decendend from it,
# and count how many are overlapping with the sequence file

$count_terminals_from_current_node =0;

############################################
get_terminals_from_this_node($current_node);#
############################################

#print "\tcount_terminals_from_current_node:$count_terminals_from_current_node\n";

if($count_terminals_from_current_node == 0)
	{

	# if non of the terminal are overlapping, record current node as a dead-end.
	# then no point to proceed further, break from current recurse

	$deadendnode{$current_node}=1;

	}else{

	# otherwise current node has decendent terminals that are overlapping

	for my $j(0 .. $child_count_current_node)
		{
		my $child_node = $nodes{$current_node}{$j};#print "\tnode:$current_node child:$child_node\n";
		if($child_node =~ /^INTERNAL_NODE_\d+/)
			{
			#################################
			find_deadend_nodes($child_node);# not reached the tips yet, so recurse.
			#################################
			}
		}
	}



}


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub check
{


unless($treefile =~ /[\w\d]/ && $fastafile =~ /[\w\d]/)
	{
	print "\nerror in command. \nto run type:\nperl print_taxa_overlapping_newick_and_fasta.pl " ,
	"[newick_file_name] [fasta_file_name]\n\n... quitting ...\n";die;
	};


# check tree file looks like a tree file:

open(IN, $treefile) || die "\nerror. cannot open the tree file you specified ($treefile)\n\n" ,
	"check this file is named correctly and in the working directory\n";
my $count_trees =0;
while(my $line=<IN>)
	{
	if($line =~ /\(/){$count_trees++}
	}
close IN;
unless($count_trees == 1)
	{die "\n\nerror, treefile ($treefile) doesnt look right. expecting newick format and a single tree. quitting.\n\n"}


# same for seq file:

open(IN, $fastafile) || die "\nerror. cannot open the fasta file you specified ($fastafile)\n\n" ,
	"check this file is named correctly and in the working directory\n";
my $count_fas =0;
while(my $line=<IN>)
	{
	if($line =~ /^>.+/){$count_fas++}
	}
close IN;
unless($count_fas >= 2)
	{die "\n\nerror, sequence file ($fastafile) doesnt look right. expecting fasta format with multiple entries. quitting.\n\n"}


close IN;



};



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




