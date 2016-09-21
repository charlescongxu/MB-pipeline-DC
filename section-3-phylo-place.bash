#######################################################################################################################################
### SECTION THREE: phylogenetic placement #############################################################################################
#######################################################################################################################################
### OPTION ONE: RAxML mothur/EPA ######################################################################################################
#######################################################################################################################################
# input files:

# pwd # be in query folder
mothuralignment=NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.align  # alignment of queries to reference, output from mothur
aligned_refs=mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping # aligned references that overlap with the supertree backbone
epa_tree=RAxML_bestTree.ref_tree_epa    # reference tree

# make a new directory (placement, RAxML_EPA) for section 3, placement option 1 and copy input files and 'format_conversion.pl' into it
# cd into the query folder
mkdir ../placement
mkdir ../placement/RAxML_EPA
cp ${aligned_refs} ../placement/RAxML_EPA
cp ${mothuralignment} ../placement/RAxML_EPA

cd .. # should be in main working folder

cp ${epa_tree} placement/RAxML_EPA

cd placement/RAxML_EPA

# we need an alignment of both the references and the queries
# we use the 'overlapping' file because it contains accessions, which matches the constraint tree

# concatenate aligned queries with 'overlapping' reference alignment and then convert from fasta to phylip format

cat ${mothuralignment} ${aligned_refs} > ${mothuralignment}_${aligned_refs}.rpq
perl ../../format_conversion.pl ${mothuralignment}_${aligned_refs}.rpq ${mothuralignment}_${aligned_refs}.rpq.phy fasta phylip

# output files:
#  ${mothuralignment}_${aligned_refs}.rpq.phy

# use RAxML with previously created reference phylogeny (RAxML_bestTree.ref_tree) to do phylogenetic placement of queries/references alignment
# "-f v": classify a bunch of environmental sequences into a reference tree
# -o Tachyglossus_aculeatus: use short-beaked echidna as outgroup
# Note: may need to remove files from previous runs of RAxML

raxmlHPC-PTHREADS-SSE3 -T 2 -f v -t ${epa_tree} -m GTRCAT -n ${mothuralignment}.raxmlEPAout -s ${mothuralignment}_${aligned_refs}.rpq.phy -o Tachyglossus_aculeatus

# Note: may need to remove some sequences consisting entirely of undetermined values (-) or else RAxML will not run, these are not the taxa you are looking for...

# sed '/_32350\|_35659/{d;}' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy > fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy

# /home/wangxy/scripts/RAxML-7.2.8-ALPHA/raxmlHPC -f v -t RAxML_bestTree.ref_tree_epa -m GTRCAT -n BWL_CROP98.cluster.align.raxmlEPAout -s fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.phy -o Tachyglossus_aculeatus

# Also note: make sure the header of your queries/references alignment is correct, the 1st number should be the number of total sequences
# otherwise you will get this stupid error: "Taxon Name too long at taxon XXX, adapt constant nmlngth in axml.h, current setting 256"
# the number of sequences is equal to 1 less than the number of lines, which you can check using 'wc <file>'

# output files:
#  RAxML_originalLabelledTree.${mothuralignment}.raxmlEPAout
#  RAxML_labelledTree.${mothuralignment}.raxmlEPAout
#  RAxML_portableTree.${mothuralignment}.raxmlEPAout.jplace
#  RAxML_info.${mothuralignment}.raxmlEPAout
#  RAxML_entropy.${mothuralignment}.raxmlEPAout
#  RAxML_classificationLikelihoodWeights. ${mothuralignment}.raxmlEPAout
#  RAxML_classification.${mothuralignment}.raxmlEPAout

# use guppy from pplacer to make a tree with each of the reads represented as a pendant edge
guppy fat RAxML_portableTree.${mothuralignment}.raxmlEPAout.jplace
guppy tog RAxML_portableTree.${mothuralignment}.raxmlEPAout.jplace

# output files:
#  fat: RAxML_portableTree.${mothuralignment}.raxmlEPAout.xml  # for forester.jar (Archeopteryx) to visualise the tree
#  tog: RAxML_portableTree.${mothuralignment}.raxmlEPAout.tog.tre  # for figtree

#######################################################################################################################################
### OPTION TWO: pynast/pplacer ########################################################################################################
#######################################################################################################################################

# NOTE: pplacer only works with older versions of RAxML so you will need to remake 'RAxML_bestTree.ref_tree' and 'RAxML_info.ref_tree' with version 7.2.7
# Unless you use taxit to create a reference package first

# input files:
pynastrefalignment=mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping_NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.fasta.pynast.rpq.fa  # pynast alignment of queries to reference, then concatenated to aligned_refs
pplacer_tree=RAxML_bestTree.ref_tree_pplacer    # reference tree
pplacer_info=RAxML_info.ref_tree_pplacer # info file
# make a new directory (pplacer) for placement option 2 and copy input files into it

mkdir ../pplacer
cd ../pplacer
cp ../../query/${pynastrefalignment} ./
cp ../../${pplacer_tree} .
cp ../../${pplacer_info} .

# use taxit to create reference package
taxit create -l pplacer -P my.refpkg --aln-fasta mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping --tree-stats ${pplacer_info} --tree-file ${pplacer_tree}

# use pplacer on taxit reference package
pplacer --keep-at-most 7 -o pynast_pplacer.jplace -p -c my.refpkg ${pynastrefalignment}

# use pplacer with reference tree and concatenated aligned queries with 'overlapping' reference alignment to do phylogenetic placement

pplacer --keep-at-most 7 -t ${pplacer_tree} -s ${pplacer_info} ${pynastrefalignment}

# Note: may need to remove some sequences consisting entirely of undetermined values (-) or else RAxML will not run, these are not the taxa you are looking for...
# sed '/_32350\|_35659/{N;d;}' aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa > fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa
# /home/wangxy/scripts/pplacer/pplacer-Linux-v1.1.alpha17/pplacer --keep-at-most 7 -t RAxML_bestTree.ref_tree_pplacer -s RAxML_info.ref_tree_pplacer fixed.aligned_refs.BWL_CROP98.cluster.align.overlapping.rpq.fa

# output files:
#  ${pynastrefalignment}.jplace # but without the .fa part

# guppy does various analyses of pplacer output
# here, we use guppy to fatten tree edges where queries are assigned
~/scripts/pplacer/guppy fat  mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping_NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.fasta.pynast.rpq.jplace
# output files:
#  ${pynastrefalignment}.xml
#  xml file can be visualised using forester.jar; be sure to click on options:  colorize branches and use branch widths

# here, we use guppy to make a tree with each of the reads represented as a pendant edge
~/scripts/pplacer/guppy tog mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping_NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.fasta.pynast.rpq.jplace
# output files:
#  ${pynastrefalignment}.tog.tre

#######################################################################################################################################
### Genesis:  this gives us the table and other output that we use to call consensus taxonomies #######################################
#######################################################################################################################################

# genesis is a toolkit for working with .jplace files. Here, we use genesis to output placement results in table format.

# load gcc and GSL module necessary for genesis
module load GSL/1.15
module load gcc/4.9.1

# use genesis to on jplace files
placement_classification_table XXX.jplace output_file

# use genesis to create nexus file from .jplace file using 'visualize_placements'
# Note: genesis is currently only installed on GRACE so you will need to upload the .jplace file there using scp
~/scripts/genesis-0.8.0/bin/visualize_placements mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping_NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.fasta.pynast.rpq.jplace BWL_CROP98_pplacer.nexus

# use genesis to create a table of the top 7 taxonomic placements for each query seq (note that my pplacer run asked for 7 to be kept)
#~/scripts/genesis-0.8.0/bin/placement_classification_table XXX.jplace output_file
#pynast/pplacer
~/scripts/genesis-0.8.0/bin/placement_classification_table mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping_NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.fasta.pynast.rpq.jplace BWL_CROP98_pplacer.txt
#mothur/epa
~/scripts/genesis-0.8.0/bin/placement_classification_table RAxML_portableTree.NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.align.raxmlEPAout.jplace BWL_CROP98_raxmlepa.txt

# The output table then needs to be edited and viewed in Megan.
