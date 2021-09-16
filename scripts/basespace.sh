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
# to install bcl2fastq download the RPM from Illumina or use '~/Dropbox/Permanent/reinstalls/software/bcl2fastq2-v2-20-0-linux-x86-64.zip'
# just unpack and copy the whole dir into '~/Software'
cd AS_RC_6samples_2pools_08092021
bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --adapter-stringency 0.9 --barcode-mismatches 0 --mask-short-adapter-reads 35 --minimum-trimmed-read-length 35 --sample-sheet "SampleSheetUsed.csv" --no-lane-splitting

## cat the libraries together
# Kilombero_Tele02_Lib01_R1.fastq.gz
cat Data/Intensities/BaseCalls/AS_A01_S4_R1_001.fastq.gz Data/Intensities/BaseCalls/AS_A04_S3_R1_001.fastq.gz Data/Intensities/BaseCalls/AS_A05_S1_R1_001.fastq.gz > Data/Intensities/BaseCalls/Kilombero_Tele02_Lib01_R1.fastq.gz
# Kilombero_Tele02_Lib01_R2.fastq.gz
cat Data/Intensities/BaseCalls/AS_A01_S4_R2_001.fastq.gz Data/Intensities/BaseCalls/AS_A04_S3_R2_001.fastq.gz Data/Intensities/BaseCalls/AS_A05_S1_R2_001.fastq.gz > Data/Intensities/BaseCalls/Kilombero_Tele02_Lib01_R2.fastq.gz
# Kilombero_Tele02_Lib02_R1.fastq.gz
cat Data/Intensities/BaseCalls/AS_A07_S2_R1_001.fastq.gz Data/Intensities/BaseCalls/AS_A08_S6_R1_001.fastq.gz Data/Intensities/BaseCalls/AS_A11_S5_R1_001.fastq.gz > Data/Intensities/BaseCalls/Kilombero_Tele02_Lib02_R1.fastq.gz
# Kilombero_Tele02_Lib02_R2.fastq.gz
cat Data/Intensities/BaseCalls/AS_A07_S2_R2_001.fastq.gz Data/Intensities/BaseCalls/AS_A08_S6_R2_001.fastq.gz Data/Intensities/BaseCalls/AS_A11_S5_R2_001.fastq.gz > Data/Intensities/BaseCalls/Kilombero_Tele02_Lib02_R2.fastq.gz

# get md5sums for all 
md5sum Data/Intensities/BaseCalls/*.fastq.gz > Data/Intensities/BaseCalls/fastq.md5

# check md5
md5sum -c Data/Intensities/BaseCalls/fastq.md5

# count reads
seqkit -j 8 stats Data/Intensities/BaseCalls/Kilombero_Tele02_Lib01_R1.fastq.gz
seqkit -j 8 stats Data/Intensities/BaseCalls/Kilombero_Tele02_Lib02_R1.fastq.gz

# copy the reads over to the rdsf
