#######################################################################################################################################
### SECTION TWO: align query sequences to reference database ##########################################################################
#######################################################################################################################################

# Note: there are three options DC listed (PYNAST, MOTHUR and SINA), only PYNAST and MOTHUR are implemented here

# input files:

queryfile    = <name_of_file_containing_query_sequences>
aligned_refs = refs_unaligned.prank.best.fas.overlapping

# make a new directory (query) for section 2 and copy input files into it
# make sure query file is also in this directory

mkdir ../query
cp refs_unaligned.prank.best.fas.overlapping ../query

cd ../query

########################################################################################################################
### PYNAST #############################################################################################################
########################################################################################################################

# Note: PyNAST must be installed and used on the local linux machine since it require modules unavilable on the barcode server

scp wangxy@10.0.16.80:<PATH OF FILE> <PATH TO DOWNLOAD TO>

# delete previous pynast output file (if present)

rm *.pynast1 *.pynast2

# use PyNAST to align query sequences to the reference database
# Note: default parameter for alignment length is way too high, set at 80 for our short reads
# clustal for pw step instead of default uclust

/home/ecec/qiime_software/pynast-1.2-release/bin/pynast --pairwise_alignment_method clustal -i BWL_CROP98.cluster.fasta -t refs_unaligned.prank.best.fas.overlapping --min_pct_id=67 --min_len=80 --fasta_out_fp=refs_unaligned.prank.best.fas.overlapping.BWL_CROP98.cluster.fasta.pynast1

# pynast modifies fasta IDs so remove these extra bits using 'format_conversion.pl'

perl format_conversion.pl refs_unaligned.prank.best.fas.overlapping.BWL_CROP98.cluster.fasta.pynast1 refs_unaligned.prank.best.fas.overlapping.BWL_CROP98.cluster.fasta.pynast2 fasta fasta

# combine reference and query alignments

cat refs_unaligned.prank.best.fas.overlapping refs_unaligned.prank.best.fas.overlapping.BWL_CROP98.cluster.fasta.pynast2 > refs_unaligned.prank.best.fas.overlapping.BWL_CROP98.cluster.fasta.pynast.rpq

# upload query to reference alignment back to barcode server

scp <PATH OF FILE> wangxy@10.0.16.80:<PATH TO UPLOAD TO>

########################################################################################################################
### MOTHUR #############################################################################################################
########################################################################################################################

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
