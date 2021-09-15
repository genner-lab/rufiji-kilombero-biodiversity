#!/usr/bin/env sh

# https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
# install CLI as above, add to PATH

# authenticate (log in)
# first log into basespace 
bs auth --force
# list prjects 
bs list projects
bs list runs

# cd
cd ../temp-local-only

# download a project - just fastq files with bases called online
#bs download project -i 292701419 -o AS_6samples_2pools

# download a run - including all the raw data
bs download run -i 213208043 -o AS_RC_6samples_2pools_08092021

# run bcl2fastq to generate fastq from scratch
cd AS_RC_6samples_2pools_08092021
bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --adapter-stringency 0.9 --barcode-mismatches 0 --mask-short-adapter-reads 35 --minimum-trimmed-read-length 35 --sample-sheet "SampleSheetUsed.csv" --no-lane-splitting

# md5sum the regenerated fastq 
md5sum IRW_LIB1_ICE_09122019_RUN/Data/Intensities/BaseCalls/LIB1_ICE_S1_R1_001.fastq.gz
md5sum IRW_LIB1_ICE_09122019_RUN/Data/Intensities/BaseCalls/LIB1_ICE_S1_R2_001.fastq.gz

# md5sums
md5sum Data/Intensities/BaseCalls/A06_S1_R1_001.fastq.gz
md5sum Data/Intensities/BaseCalls/A06_S1_R2_001.fastq.gz
md5sum Data/Intensities/BaseCalls/A12_S2_R1_001.fastq.gz
md5sum Data/Intensities/BaseCalls/A12_S2_R2_001.fastq.gz

# 75aca1160bc5989fc2766725edd3b78b A06_S1_R1_001.fastq.gz
# cbd3fa76d0a61c53db9bde41978d21b1 A06_S1_R2_001.fastq.gz
# 0b8eebf52b1185a5595402e1bfe69851 A12_S2_R1_001.fastq.gz
# 419abc4af510cbccaa5a11605cd652cd A12_S2_R2_001.fastq.gz

# merge runs
cat Data/Intensities/BaseCalls/A06_S1_R1_001.fastq.gz Data/Intensities/BaseCalls/A12_S2_R1_001.fastq.gz > Data/Intensities/BaseCalls/SeaDNA_Tele02_Lib03v2_R1.fastq.gz
cat Data/Intensities/BaseCalls/A06_S1_R2_001.fastq.gz Data/Intensities/BaseCalls/A12_S2_R2_001.fastq.gz > Data/Intensities/BaseCalls/SeaDNA_Tele02_Lib03v2_R2.fastq.gz

# md5sums
md5sum Data/Intensities/BaseCalls/SeaDNA_Tele02_Lib03v2_R1.fastq.gz
md5sum Data/Intensities/BaseCalls/SeaDNA_Tele02_Lib03v2_R2.fastq.gz

# 1d2e263e21236cf42486e1037ba936f8 SeaDNA_Tele02_Lib03v2_R1.fastq.gz
# 1f27da3c73633ba5355639d4e683b941 SeaDNA_Tele02_Lib03v2_R2.fastq.gz
