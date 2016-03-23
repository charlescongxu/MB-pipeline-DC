#######################################################################################################################################
### SECTION THREE: phylogenetic placement #############################################################################################
#######################################################################################################################################
### OPTION ONE: RAxML EPA #############################################################################################################
#######################################################################################################################################
# input files:

BWL_CROP98.cluster.align                  = alignment of queries to reference, output from mothur
refs_unaligned.prank.best.fas.overlapping = aligned references that overlap with the supertree backbone
RAxML_bestTree.ref_tree                   = reference tree

# make a new directory (placement, RAxML_EPA) for section 3, placement option 1 and copy input files and 'format_conversion.pl' into it

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

sed '/_32350\|_35659/{N;d;}' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy > fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy

/home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -f v -t RAxML_bestTree.ref_tree -m GTRCAT -n BWL_CROP98.cluster.align.raxmlEPAout -s fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy -o Tachyglossus_aculeatus

# Also note: make sure the header of your queries/references alignment is correct, the 1st number should be the number of total sequences
# otherwise you will get this stupid error: "Taxon Name too long at taxon XXX, adapt constant nmlngth in axml.h, current setting 256"
# the number of sequences is equal to 1 less than the number of lines, which you can check using 'wc <file>'

# output files:
#  

#######################################################################################################################################
### OPTION TWO: pplacer ###############################################################################################################
#######################################################################################################################################

# input files:

BWL_CROP98.cluster.align = alignment of queries to reference, output from mothur
refs_unaligned.prank.best.fas.overlapping = aligned references that overlap with the supertree backbone
RAxML_bestTree.ref_tree  = reference tree from section two
RAxML_info.ref_tree      = info file of reference tree from section two

# make a new directory (pplacer) for placement option 2 and copy input files into it

mkdir ../pplacer
cd ../pplacer
cp ../../query/BWL_CROP98.cluster.align .
cp ../../query/refs_unaligned.prank.best.fas.overlapping .
cp ../../references/RAxML_bestTree.ref_tree .
cp ../../references/RAxML_info.ref_tree .

# concatenate aligned queries with 'overlapping' reference alignment and then convert from fasta to phylip format

cat refs_unaligned.prank.best.fas.overlapping BWL_CROP98.cluster.align > aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa

# use pplacer with reference tree and concatenated aligned queries with 'overlapping' reference alignment to do phylogenetic placement

pplacer-Linux-v1.1.alpha17/pplacer --keep-at-most 1 -t RAxML_bestTree.ref_tree -s RAxML_info.ref_tree aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa

# Note: may need to remove some sequences consisting entirely of undetermined values (-) or else RAxML will not run, these are not the taxa you are looking for...

sed '/_32350\|_35659/{N;d;}' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa > fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa

pplacer-Linux-v1.1.alpha17/pplacer --keep-at-most 1 -t RAxML_bestTree.ref_tree -s RAxML_info.ref_tree fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa

# output files:
#  fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.jplace

# guppy does various analyses of pplacer output
# here, we use guppy to fatten tree edges where queries are assigned
pplacer-Linux-v1.1.alpha17/guppy fat fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.jplace

# output files:
#  fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.xml

# here, we use guppy to make a tree with each of the reads represented as a pendant edge
pplacer-Linux-v1.1.alpha17/guppy tog fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.jplace

# output files:
#  fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.tog.tre

#######################################################################################################################################
### OPTION THREE: bagpipe phylo #######################################################################################################
#######################################################################################################################################

# necessary script:

bagpipe_phylo.pl
(https://sourceforge.net/projects/bagpipe/)

# input files:
refs_unaligned.prank.best.fas.overlapping = aligned references that overlap with the supertree backbone
BE_mammal_supertree.nwk.overlapping       = overlapping supertree backbone (from section 1 when making the reference database)
BWL_CROP98.cluster.align                  = alignment of queries to reference, output from mothur

# make a new directory (bagpipe_phylo) for placement option 3 and copy input files and 'format_conversion.pl' into it

mkdir ../bagpipe_phylo
cd ../bagpipe_phylo

cp ../../query/refs_unaligned.prank.best.fas.overlapping .
cp ../../query/BWL_CROP98.cluster.align .
cp ../../references/BE_mammal_supertree.nwk.overlapping .
cp ../../references/format_conversion.pl .

# we need an alignment of both the references and the queries
# we use the 'overlapping' file because it contains accessions rm, which matches the constraint tree

# concatenate aligned queries with 'overlapping' reference alignment and then convert from fasta to phylip format

cat refs_unaligned.prank.best.fas.overlapping BWL_CROP98.cluster.align > aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq
perl format_conversion.pl aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy fasta phylip

# output files:
#  aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy

# use RAxML to create a constraint tree from the input files

/home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -s aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy -n BWL_CROP98.constrn -m GTRCAT -c 25 -p 12345 -g BE_mammal_supertree.nwk.overlapping

# output files:
#  aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy.reduced
#  RAxML_bestTree.BWL_CROP98.constrn
#  RAxML_info.BWL_CROP98.constrn
#  RAxML_log.BWL_CROP98.constrn
#  RAxML_result.BWL_CROP98.constrn

# Note: may need to remove some sequences consisting entirely of undetermined values (-) or else RAxML will not run, these are not the taxa you are looking for...

sed '/_32350\|_35659/{N;d;}' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy > fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy

/home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -s fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy -n BWL_CROP98.constrn -m GTRCAT -c 25 -p 12345 -g BE_mammal_supertree.nwk.overlapping

# Also note: make sure the header of your queries/references alignment is correct, the 1st number should be the number of total sequences
# otherwise you will get this stupid error: "Taxon Name too long at taxon XXX, adapt constant nmlngth in axml.h, current setting 256"
# the number of sequences is equal to 1 less than the number of lines, which you can check using 'wc <file>'

# output files:
#  

# copy nodes.dmp and names.dmp from section 1 when making the reference database to bagpipe directory

cp ../../references/nodes.dmp .
cp ../../references/names.dmp .

# use perl script 'bagpipe_phylo.pl' to do phylogenetic placement using the constrained RAxML tree and the aligned queries

perl bagpipe_phylo.pl -treefile RAxML_bestTree.BWL_CROP98.constrn -seqfile BWL_CROP98.cluster.align -support 0 -node 40674

# Note: may need to remove some sequences consisting entirely of undetermined values (-) or else RAxML will not run, these are not the taxa you are looking for...

sed '/_32350\|_35659/{N;d;}' BWL_CROP98.cluster.align > fixed.BWL_CROP98.cluster.align

perl bagpipe_phylo.pl -treefile RAxML_bestTree.BWL_CROP98.constrn -seqfile fixed.BWL_CROP98.cluster.align -support 0 -node 40674

# output files:
#  RAxML_bestTree.BWL_CROP98.constrn.query_clades
