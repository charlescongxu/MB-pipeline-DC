#######################################################################################################################################
### SECTION TWO: align query sequences to reference database ##########################################################################
#######################################################################################################################################
### Only need to run this script once per query dataset ###############################################################################
#######################################################################################################################################

# input files:

queryfile=NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.fasta # query OTUs
aligned_refs=mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping # these are the aligned reference sequences, formatted

# make a new directory (query) for section 2 and copy input files into it
# make sure query file is also in this directory

mkdir query
cp ${queryfile} query
cp ${aligned_refs} query
cd query

#######################################################################################################################################
### PYNAST for pynast/pplacer #########################################################################################################
#######################################################################################################################################

# Note: PyNAST must be installed and used on the local linux machine since it require modules unavilable on the barcode server
# scp wangxy@10.0.16.80:<PATH OF FILE> <PATH TO DOWNLOAD TO>

# delete previous pynast output file (if present)
rm *.pynast1 *.pynast2

# may need to load python 2.7.9 first
module load python/2.7.9

# use PyNAST to align query sequences to the reference database
# Note: default parameter for alignment length is way too high, set at 80 for our short reads
# clustal for pw step instead of default uclust
pynast --pairwise_alignment_method clustal -i ${queryfile} -t ${aligned_refs} --min_pct_id=67 --min_len=80 --fasta_out_fp="${aligned_refs}_${queryfile}.pynast1"

# ${queryfile}_pynast_fail.fasta holds sequences that pynast can't align

# pynast modifies fasta IDs so remove these extra bits using 'format_conversion.pl'
# or use seqtk: seqtk seq -AC "${aligned_refs}_${queryfile}.pynast1" > "${aligned_refs}_${queryfile}.pynast2"
perl ../format_conversion.pl "${aligned_refs}_${queryfile}.pynast1" "${aligned_refs}_${queryfile}.pynast2" fasta fasta

# combine reference and query alignments
cat "${aligned_refs}" "${aligned_refs}_${queryfile}.pynast2" > "${aligned_refs}_${queryfile}.pynast.rpq.fa"

########################################################################################################################################
### MOTHUR for mothur/epa ##############################################################################################################
########################################################################################################################################

# use mothur to align query sequences to the reference database
# key options you can vary are:
#  search = kmer/suffix/blast
#  align  = needleman/gotoh

mothur # starts mothur environment
align.seqs(candidate=${queryfile}, template=${aligned_refs}, search=kmer, align=needleman, threshold=0.85)
quit()

# output files:
#  mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.8mer
#  NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.flip.accnos
#  NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.align.report
#  NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.align
#  a mothur logfile

# the output (NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.align) needs some processing now
# we need to replace all '.' with '-' for input into RAxML to do phylogenetic placement

# Linux syntax
sed -i 's/\./-/g' NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.align

# macOS syntax
sed -i '' 's/\./-/g' NYM_BFCusearchMF_CROP98.cluster_Blast_Filtered.align

########################################################################################################################################
### PaPaRa #############################################################################################################################
########################################################################################################################################

#use format_conversion.pl to convert reference alignment fasta to phylip
perl format_conversion.pl mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy fasta phylip

#output:
#mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy

#The phylip file (option -s) must contain the reference alignment (RA), consistent with the reference tree (option -t).
#The FASTA file (option -q) contains the unaligned query sequences (QS)
#The output alignment will be written to papara_alignment.default.
#You can change the file suffix (i.e., “default”) by supplying a run-name with parameter -n
papara -t RAxML_bestTree.ref_tree_epa -s mamDB16S_full_unambig.ng.rr.ID_filtered.fixed.mafft.overlapping.phy -q ../../16s_seqs_nonchimeric_rep_set__20plus_allmaps_nopighumbac_mammals.fasta -n 16s_torrey

#output:
#papara_log.16s_torrey (papara log file)
#papara_alignment.16s_torrey (papara alignment file in phylip format)
#papara_quality.16s_torrey (empty??)
