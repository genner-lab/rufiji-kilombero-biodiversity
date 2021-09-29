# rufiji-kilombero-biodiversity
12S eDNA biodiversity survey of the Kilombero and Rufiji rivers

```
# clone repos
git clone https://github.com/genner-lab/rufiji-kilombero-biodiversity.git
cd rufiji-kilombero-biodiversity
git clone https://github.com/genner-lab/meta-fish-pipe.git
git clone https://github.com/genner-lab/meta-fish-lib.git
git clone https://github.com/genner-lab/refseq-reflib.git

# restore R libs
Rscript -e "renv::restore()"
# get data from ncbi
#scripts/get-data.sh

# assemble fish reference library 
cd meta-fish-lib
mkdir -p reports temp/fasta-temp
cp ../temp/ncbi-key.R assets/ncbi-key.R
Rscript -e "renv::restore()"
cp ../assets/species-table.csv assets/species-table.csv
scripts/sequences-download.R -t 4 -e false
scripts/references-assemble.R -t 4 -m 12s.taberlet
scripts/qc.R -p raxmlHPC -t 1

#(and add local seqs)
#locals <- read_csv(file="assets/local-12s.csv")
#reflib.orig %>% bind_rows(locals) %>% write_csv(file="meta-fish-pipe/assets/meta-fish-lib-v243.csv")

# get refseq
cd refseq-reflib
mkdir temp references
Rscript -e "renv::restore()"
scripts/download.sh
scripts/extract.R -p tele02
scripts/annotate.R -s 42 -p tele02
cd ..
cp refseq-reflib/references/refseq208-annotated-tele02.csv meta-fish-pipe/assets/refseq208-annotated-tele02.csv

# copy across sample sheet and contam file to the pipeline lib
cp assets/sequencing-master.csv meta-fish-pipe/assets/sequencing-master.csv
cp assets/contaminants-exclude.csv meta-fish-pipe/assets/contaminants-exclude.csv

# set up pipeline
cd meta-fish-pipe
Rscript -e "renv::restore()"
scripts/session-info.sh  -r assets/refseq208-annotated-tele02.csv -c assets/meta-fish-lib-v243.csv

# set up libs
scripts/prepare-libraries.sh -p tele02 -l lib1
scripts/prepare-libraries.sh -p tele02 -l lib2

cd ..
# make symlinks LIB1
ln -s -r temp-local-only/AS_RC_6samples_2pools_08092021/Data/Intensities/BaseCalls/Kilombero_Tele02_Lib01_R1.fastq.gz meta-fish-pipe/temp/processing/tele02-lib1/fastq/R1.fastq.gz
ln -s -r temp-local-only/AS_RC_6samples_2pools_08092021/Data/Intensities/BaseCalls/Kilombero_Tele02_Lib01_R2.fastq.gz meta-fish-pipe/temp/processing/tele02-lib1/fastq/R2.fastq.gz
# make symlinks LIB2
ln -s -r temp-local-only/AS_RC_6samples_2pools_08092021/Data/Intensities/BaseCalls/Kilombero_Tele02_Lib02_R1.fastq.gz meta-fish-pipe/temp/processing/tele02-lib2/fastq/R1.fastq.gz
ln -s -r temp-local-only/AS_RC_6samples_2pools_08092021/Data/Intensities/BaseCalls/Kilombero_Tele02_Lib02_R2.fastq.gz meta-fish-pipe/temp/processing/tele02-lib2/fastq/R2.fastq.gz

# generate barcodes
cd meta-fish-pipe
scripts/generate-barcodes.R -p tele02 -l lib1 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib2 -f 18 -r 20 -m assets/sequencing-master.csv

# demultiplex
scripts/demultiplex.sh -p tele02 -l lib1 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib2 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

######
got up to here - need reflib priors
###### 

# denoise with dada2
scripts/dada2.R -p tele02 -l lib1
scripts/dada2.R -p tele02 -l lib2

# generate stats
scripts/generate-stats.sh -p tele02 -l lib1 -t 8
scripts/generate-stats.sh -p tele02 -l lib2 -t 8

# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p tele02

# assemble results
scripts/assemble-results.R -c assets/contaminants-exclude.csv
```