#!/usr/bin/env sh

# make db
makeblastdb -in temp-local-only/asvs-lulu.fasta -parse_seqids -dbtype nucl

# blast 
blastn -db temp-local-only/asvs-lulu.fasta -outfmt '6 qseqid sseqid pident' -out temp-local-only/lulu-blast-matchlist.tsv -qcov_hsp_perc 80 -perc_identity 84 -query temp-local-only/asvs-lulu.fasta
