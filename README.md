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
echo 'ncbi.key <- "my-ncbi-key"' > assets/ncbi-key.R
Rscript -e "renv::restore()"
cp ../assets/species-table.csv assets/species-table.csv
scripts/sequences-download.R -t 4 -e true
scripts/references-assemble.R -t 4 -m 12s.taberlet



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
rm temp/duckdb
cd ..
cp refseq-reflib/references/refseq206-annotated-tele02.csv meta-fish-pipe/assets/refseq206-annotated-tele02.csv

# copy across sample sheet and contam file to the pipeline lib
cp assets/sequencing-master.csv meta-fish-pipe/assets/sequencing-master.csv
cp assets/contaminants-exclude.csv meta-fish-pipe/assets/contaminants-exclude.csv

# set up pipeline
cd meta-fish-pipe
Rscript -e "renv::restore()"
scripts/session-info.sh  -r assets/refseq206-annotated-tele02.csv -c assets/meta-fish-lib-v243.csv

# set up libs
scripts/prepare-libraries.sh -p tele02 -l lib1
scripts/prepare-libraries.sh -p tele02 -l lib2
scripts/prepare-libraries.sh -p tele02 -l lib3
scripts/prepare-libraries.sh -p tele02 -l lib4

cd ..
# make symlinks LIB1
ln -s -r temp/data/SeaDNA_Teleo02_01_S1_L001_R1_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib1/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Teleo02_01_S1_L001_R2_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib1/fastq/R2.fastq.gz
# make symlinks LIB2
ln -s -r temp/data/SeaDNA_Teleo02_02_S2_L001_R1_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib2/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Teleo02_02_S2_L001_R2_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib2/fastq/R2.fastq.gz
# make symlinks LIB3
ln -s -r temp/data/SeaDNA_Tele02_Lib03v2_R1.fastq.gz.1 meta-fish-pipe/temp/processing/tele02-lib3/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Tele02_Lib03v2_R2.fastq.gz.1 meta-fish-pipe/temp/processing/tele02-lib3/fastq/R2.fastq.gz
# make symlinks LIB4
ln -s -r temp/data/SeaDNA_Teleo02_Lib-04_S2_L001_R1_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib4/fastq/R1.fastq.gz
ln -s -r temp/data/SeaDNA_Teleo02_Lib-04_S2_L001_R2_001.fastq.gz meta-fish-pipe/temp/processing/tele02-lib4/fastq/R2.fastq.gz

# generate barcodes
cd meta-fish-pipe
scripts/generate-barcodes.R -p tele02 -l lib1 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib2 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib3 -f 18 -r 20 -m assets/sequencing-master.csv
scripts/generate-barcodes.R -p tele02 -l lib4 -f 18 -r 20 -m assets/sequencing-master.csv

# demultiplex
scripts/demultiplex.sh -p tele02 -l lib1 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib2 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib3 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
scripts/demultiplex.sh -p tele02 -l lib4 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib1
scripts/dada2.R -p tele02 -l lib2
scripts/dada2.R -p tele02 -l lib3
scripts/dada2.R -p tele02 -l lib4

# generate stats
scripts/generate-stats.sh -p tele02 -l lib1 -t 8
scripts/generate-stats.sh -p tele02 -l lib2 -t 8
scripts/generate-stats.sh -p tele02 -l lib3 -t 8
scripts/generate-stats.sh -p tele02 -l lib4 -t 8

# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p tele02

# assemble results
scripts/assemble-results.R -c assets/contaminants-exclude.csv
```