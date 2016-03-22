#######################################################################################################################################
### SECTION THREE: phylogenetic placement #############################################################################################
#######################################################################################################################################
### OPTION ONE: RAxML EPA #############################################################################################################
#######################################################################################################################################
# input files:

queryfile      = <name_of_file_containing_query_sequences>
aligned_refs   = refs_unaligned.prank.best.fas.overlapping
reference_tree = RAxML_bestTree.ref_tree

# make a new directory (placement, RAxML_EPA) for section 3, placement option 1 and copy input files, 'format_conversion.pl' into it

mkdir ../placement
mkdir ../placement/RAxML_EPA
cp refs_unaligned.prank.best.fas.overlapping ../placement/RAxML_EPA
cp BWL_CROP98.cluster.align ../placement/RAxML_EPA

cd ../references
cp format_conversion.pl ../placement/RAxML_EPA
cp RAxML_bestTree.ref_tree ../placement/RAxML_EPA

cd ../placement/RAxML_EPA

# we need an alignment of both the references and the queries
# we use the 'overlapping' file because it contains accessions rm, which matches the constraint tree

# concatenate aligned queries with 'overlapping' reference alignment and then convert from fasta to phylip format
cat refs_unaligned.prank.best.fas.overlapping BWL_CROP98.cluster.align > aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq
perl format_conversion.pl aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy fasta phylip

# output files:
#  aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy

# use RAxML with previously created reference phylogeny (RAxML_bestTree.ref_tree) to do phylogenetic placement of queries/references alignment
# "-f v": classify a bunch of environmental sequences into a reference tree
# -o Tachyglossus_aculeatus: use short-beaked echidna as outgroup
# Note: may need to remove files from previous runs of RAxML

/home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -f v -t RAxML_bestTree.ref_tree -m GTRCAT -n BWL_CROP98.cluster.align.raxmlEPAout -s aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy -o Tachyglossus_aculeatus

# Note: may need to remove some sequences consisting entirely of undetermined values (-) or else RAxML will not run, these are not the taxa you are looking for...

sed --in-place '/all_seqs_89to121bp_noChimera_uniques_1199/d' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy
sed --in-place '/all_seqs_89to121bp_noChimera_uniques_44259/d' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy
sed --in-place '/all_seqs_89to121bp_noChimera_uniques_23857/d' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy
sed --in-place '/all_seqs_89to121bp_noChimera_uniques_24008/d' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy
sed --in-place '/all_seqs_89to121bp_noChimera_uniques_63886/d' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy

# Also note: make sure the header of your queries/references alignment is correct, the 1st number should be the number of total sequences
#            otherwise you will get this stupid error
#            "Taxon Name too long at taxon XXX, adapt constant nmlngth in axml.h, current setting 256"

#######################################################################################################################################
### OPTION THREE: bagpipe phylo #######################################################################################################
#######################################################################################################################################

# necessary script:

bagpipe_phylo.pl
(https://sourceforge.net/projects/bagpipe/)

# input files:

aligned queries and references that overlap with supertree backbone in phylip format (from RAxML EPA placement option)
overlapping supertree backbone (from section 1 when making the reference database)

# make a new directory (bagpipe_phylo) for placement option 3 and copy input files into it

mkdir ../bagpipe_phylo

# copy nodes.dmp and names.dmp from section 1 when making the reference database

cp nodes.dmp ../placement/bagpipe_phylo
cp names.dmp ../placement/bagpipe_phylo

# use RAxML to create a constraint tree from the input files

/home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -s aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy -n BWL_CROP98.constrn -m GTRCAT -c 25 -p 12345 -g BE_mammal_supertree.nwk.overlapping

# output files:
#  aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy.reduced
#  RAxML_bestTree.BWL_CROP98.constrn
#  RAxML_info.BWL_CROP98.constrn
#  RAxML_log.BWL_CROP98.constrn
#  RAxML_result.BWL_CROP98.constrn

# copy aligned queries (output from mothur) to bagpipe directory

cp ../../query/BWL_CROP98.cluster.align .

# use perl script 'bagpipe_phylo.pl' to do phylogenetic placement using the constrained RAxML tree and the aligned queries

perl bagpipe_phylo.pl -treefile RAxML_bestTree.BWL_CROP98.constrn -seqfile BWL_CROP98.cluster.align -support 0 -node 40674
