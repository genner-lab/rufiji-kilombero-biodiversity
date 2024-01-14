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
cp ../assets/exclusions.csv assets/exclusions.csv
scripts/sequences-download.R -q 2000 -t 3 -e false
scripts/references-assemble.R -t 4 -m 12s.taberlet
# backup (or copy back)
cp assets/reference-library-master.csv.gz assets/reference-library-master-nolocal.csv.gz
#cp assets/reference-library-master-nolocal.csv.gz assets/reference-library-master.csv.gz

# add local seqs and qc in R (need to make into script)
library("here")
library("tidyverse")
source(here::here("scripts/references-load-local.R"))
source(here::here("scripts/references-clean.R"))
locals <- read_csv(file=here("../assets/local-12s.csv"))
reflib.local <- reflib.cleaned %>% mutate(dbid=as.character(dbid)) %>% bind_rows(locals) %>% arrange(phylum,class,order,family,genus,sciNameValid,dbid)
reflib.local %>% write_csv(file=gzfile(here("assets/reference-library-master.csv.gz")), na="")
reflib.local %>% write_csv(file=here("../meta-fish-pipe/assets/meta-fish-lib-v259.csv"), na="")

# qc
scripts/qc.R -p ~/Software/standard-RAxML/raxmlHPC-AVX -t 1
scripts/qc.R -p raxmlHPC -t 1

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
scripts/session-info.sh  -r assets/refseq208-annotated-tele02.csv -c assets/meta-fish-lib-v259.csv

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

# jan 2022 - reanalyse with new reference library
# in reflib copy files
cd meta-fish-lib
cp ../assets/exclusions.csv assets/exclusions.csv
cp assets/reference-library-master-nolocal.csv.gz assets/reference-library-master.csv.gz
# in R add local seqs
library("here")
library("tidyverse")
source(here::here("scripts/references-load-local.R"))
source(here::here("scripts/references-clean.R"))
locals <- read_csv(file=here("../assets/local-12s.csv"))
reflib.local <- reflib.orig %>% mutate(dbid=as.character(dbid)) %>% bind_rows(locals) %>% arrange(phylum,class,order,family,genus,sciNameValid,dbid)
reflib.local %>% write_csv(file=gzfile(here("assets/reference-library-master.csv.gz")), na="")
reflib.local %>% write_csv(file=here("../meta-fish-pipe/assets/meta-fish-lib-v245.csv"), na="")
# qc
scripts/qc.R -p raxmlHPC -t 2

# in pipe
cd meta-fish-pipe
cp assets/meta-fish-lib-v245.csv temp/taxonomic-assignment/custom-reference-library.csv
scripts/taxonomic-assignment.sh -t 8 -p tele02
scripts/assemble-results.R -c assets/contaminants-exclude.csv
```

########### jan 2024 renanalyse with new ref lib (do in meta-fish-pipe) ##############

# redo reflib
```sh
scripts/references-assemble.R -t 4 -m 12s.taberlet
```

```r
cd meta-fish-lib-v259
cp ../assets/exclusions.csv assets/exclusions.csv
# in R add local seqs
library("here")
library("tidyverse")
source(here::here("scripts/references-load-local.R"))
source(here::here("scripts/references-clean.R"))
locals <- read_csv(file=here("../assets/local-12s.csv"))
reflib.local <- reflib.cleaned %>% mutate(dbid=as.character(dbid)) %>% bind_rows(locals) %>% arrange(phylum,class,order,family,genus,sciNameValid,dbid)
reflib.local %>% write_csv(file=gzfile(here("assets/reference-library-master.csv.gz")), na="")
#reflib.local %>% write_csv(file=here("../meta-fish-pipe/assets/meta-fish-lib-v259.csv"), na="")

```sh
scripts/qc.R -t 4 -v false
```

# check the exclusions
```r
library("here")
library("tidyverse")
library("ape")
library("phangorn")
tr <- ape::read.tree(file="temp/qc_v259_Jan-2024/12s.taberlet.noprimers.fas.ali.raxml.rba.raxml.bestTree")
tr <- ape::ladderize(phangorn::midpoint(tr))
exc <- read_csv(file=here("../assets/exclusions.csv"))
cols <- rep("gray20",length(tr$tip.label))
cols[str_split_fixed(tr$tip.label,"\\|",3)[,1] %in% pull(exc,dbid)] <- "hotpink"
# plot tree
pdf(file=here("temp/qc_v259_Jan-2024/12s.taberlet.references.exluded.pdf"), width=15, height=length(tr$tip.label)/10)
plot.phylo(tr, tip.col=cols, cex=0.5, font=1, label.offset=0.01, no.margin=TRUE)
dev.off()
```


# rerun assignment
```bash
# in meta-fish-pipe
cd meta-fish-pipe
cp path/to/meta-fish-lib-v259.csv temp/taxonomic-assignment/custom-reference-library.csv
scripts/prep-reflibs.R -p tele02
vsearch --threads 8 --sintax temp/taxonomic-assignment/asvs-clean-cat-relabel-derep.fasta --db temp/taxonomic-assignment/refseq-annotated.fasta --sintax_cutoff 0.7 --tabbedout temp/taxonomic-assignment/sintax-output.tsv
scripts/extract-fishes.R
makeblastdb -in temp/taxonomic-assignment/custom-reference-library.fasta -dbtype nucl -blastdb_version 5
blastn -task blastn -num_threads 8 -evalue 1000 -word_size 7 -max_target_seqs 500 -db temp/taxonomic-assignment/custom-reference-library.fasta -outfmt "6 qseqid sseqid evalue length pident nident score bitscore" -out temp/taxonomic-assignment/fish-blast.out -query temp/taxonomic-assignment/asvs-fish.fasta
echo -e 'asv\tblastDbid\tblastEvalue\tblastLength\tblastPident\tblastNident\tblastScore\tblastBitscore' | sed -e "s/-e //g" > temp/taxonomic-assignment/headers
cat temp/taxonomic-assignment/headers temp/taxonomic-assignment/fish-blast.out > temp/taxonomic-assignment/fish-blast-result.tsv
rm temp/taxonomic-assignment/fish-blast.out
rm temp/taxonomic-assignment/headers
scripts/process-blast.R
scripts/prep-epa.R
mafft --thread 8 --retree 2 --maxiterate 2 temp/taxonomic-assignment/epa/epa.input.fasta > temp/taxonomic-assignment/epa/epa.aligned.fasta
scripts/split.R
raxml-ng --parse --msa temp/taxonomic-assignment/epa/epa.references.fasta --model TN93+G --redo
raxml-ng --search --msa temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba --tree pars{1} --seed 42 --threads auto{8} --workers auto --lh-epsilon 0.01 --redo
raxml-ng --evaluate --msa temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba --tree temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba.raxml.bestTree --prefix temp/taxonomic-assignment/epa/opt --redo
epa-ng --ref-msa temp/taxonomic-assignment/epa/epa.references.fasta --tree temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba.raxml.bestTree --query temp/taxonomic-assignment/epa/epa.queries.fasta --outdir temp/taxonomic-assignment/epa --model temp/taxonomic-assignment/epa/opt.raxml.bestModel --redo --preserve-rooting off
gappa examine assign --per-query-results --jplace-path temp/taxonomic-assignment/epa/epa_result.jplace --taxon-file temp/taxonomic-assignment/epa/references.taxonomy.tsv --out-dir temp/taxonomic-assignment/epa --allow-file-overwriting
mv temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba.raxml.bestTree temp/taxonomic-assignment/epa/references.tree.nwk
rm temp/taxonomic-assignment/epa/*.log
rm temp/taxonomic-assignment/epa/*raxml*
scripts/process-epa.R
# assemble
scripts/assemble-results.R -c assets/contaminants-exclude.csv
```

```r
library("here")
library("tidyverse")
library("ape")
library("phangorn")

# read csv
tax.ass.new <- read_csv(here("results/taxonomic-assignments.csv"),guess_max=9999)
tax.ass.old <- read_csv(here("results/taxonomic-assignments-old.csv"),guess_max=9999)

# get old list 
tax.ass.old.red <- tax.ass.old |> select(asvHash,assignedName) |> rename(assignedNameOld=assignedName)

# join to new list
tax.ass.new |> 
    left_join(tax.ass.old.red,by="asvHash") |> 
    relocate(assignedNameOld,.after="assignedName") |> 
    select(asvHash,nreads,isFish,assigned,assignedName,assignedNameOld) |> 
    mutate(same=if_else(assignedName==assignedNameOld,".","DIFFERENT")) |>
    write_csv(here("results/taxonomic-assignments-differences.csv"))

# tree
tr <- ape::read.tree(file="results/12s.taberlet.noprimers.fas.ali.raxml.rba.raxml.bestTree")
tr <- ape::ladderize(phangorn::midpoint(tr))

# fix names
ass.spp <- tax.ass.new |> filter(isFish==TRUE & assigned==TRUE) |> distinct(assignedName) |> pull(assignedName) |> str_replace_all(" ","_") |> str_replace_all("'","_")

# check names
ass.spp %in% str_split_fixed(tr$tip.label,"\\|",3)[,2]

# make cols
cols <- rep("gray20",length(tr$tip.label))
cols[str_split_fixed(tr$tip.label,"\\|",3)[,2] %in% ass.spp] <- "hotpink"

# plot tree
pdf(file=here("results/12s.taberlet.noprimers.fas.ali.raxml.rba.raxml.bestTree.pdf"), width=15, height=length(tr$tip.label)/10)
plot.phylo(tr, tip.col=cols, cex=0.5, font=1, label.offset=0.01, no.margin=TRUE)
dev.off()

```

