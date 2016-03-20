#######################################################################################################################################
### SECTION TWO: align query sequences to reference database ##########################################################################
#######################################################################################################################################

# Note: there are three options DC listed (PYNAST, MOTHUR and SINA), only MOTHUR is implemented here

# input files:

queryfile    = <name_of_file_containing_query_sequences>
aligned_refs = refs_unaligned.prank.best.fas.overlapping

# make a new directory (query) for section 2 and copy input files into it
# make sure query file is also in this directory

mkdir ../query
cp refs_unaligned.prank.best.fas.overlapping ../query

cd ../query

# use mothur to align query sequences to the reference database
# key options you can vary are:
#  search = kmer/suffix/blast
#  align  = needleman/gotoh

/home/wangxy/scripts/mothur/mothur
align.seqs(candidate=BWL_CROP98.cluster.fasta, template=refs_unaligned.prank.best.fas.overlapping, search=kmer, align=needleman, threshold=0.85)

# output files:
#  BWL_CROP98.cluster.align
#  BWL_CROP98.cluster.align.report
#  BWL_CROP98.cluster.flip.accnos
#  refs_unaligned.prank.best.fas.8mer
#  a mothur logfile

# the output (BWL_CROP98.cluster.align) needs some processing now
# we need to replace all '.' with '-' for input into RAxML to do phylogenetic placement

sed -i 's/\./-/g' BWL_CROP98.cluster.align
