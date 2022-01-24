#!/usr/bin/env Rscript

# lib
library("here")
library("tidyverse")
library("lulu")
library("ape")
library("phangorn")
library("RColorBrewer")
library("ggtree")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")
# https://github.com/tobiasgf/lulu
options(width=180)


# NEW RAXML-NG FUN
raxml_ng <- function(file) {
        string.mafft <- paste0("mafft --thread -1 --maxiterate 2 --retree 2 ",file," > ",file,".ali")
        system(command=string.mafft,ignore.stdout=FALSE)
        string.parse <- paste0("raxml-ng --parse --msa ",file,".ali --model TN93+G --seed 42 --redo --threads auto")
        system(command=string.parse,ignore.stdout=FALSE)
        string.search <- paste0("raxml-ng --search --msa ",file,".ali.raxml.rba --tree pars{1} --lh-epsilon 0.1 --seed 42 --redo --threads auto")
        system(command=string.search,ignore.stdout=FALSE)
        rax.tr <- ape::read.tree(file=paste0(file,".ali.raxml.rba.raxml.bestTree"))
    return(rax.tr)
}


# run mptp
run_mptp <- function(file) {
        string.mptp <- paste0("mptp --ml --single --tree_file ",file," --output_file ",file,".mptp.out")
        #  --minbr 0.0001
        #system("export PATH=~/Software/mptp/bin:$PATH")
        system(command=string.mptp,ignore.stdout=FALSE)
        #x <- x
    #return(string.mptp)
}


# FUNCTION TO READ MPTP OUTPUT FILES
read_mptp <- function(file,skiplines) {
    mptp.raw <- readr::read_delim(file,skip=skiplines,delim=",",col_names="individual",show_col_types=FALSE)
    mptp.tab <- mptp.raw %>% 
        dplyr::mutate(species=ifelse(grepl(":",individual),individual,NA)) %>%
        tidyr::fill(species,.direction="down") %>%
        dplyr::filter(!grepl(":",individual)) %>%
        dplyr::mutate(species=stringr::str_replace_all(species,":| ","")) %>%
        dplyr::relocate(species,.before=individual)
    return(mptp.tab)
}
