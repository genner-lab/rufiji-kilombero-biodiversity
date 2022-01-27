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
run_mptp <- function(file,threshold,minbr) {
    string.mptp <- paste0("mptp --ml --",threshold," --minbr ",minbr," --tree_file ",file," --output_file ",file,".mptp.out")
    system(command=string.mptp,ignore.stdout=FALSE)
    #return(string.mptp)
}


# FUNCTION TO READ MPTP OUTPUT FILES
read_mptp <- function(file) {
    mptp.scan <- scan(file=file,what="character",sep="\n",quiet=TRUE)
    skiplines <- grep("Species 1:",mptp.scan)
    skiplines <- skiplines - 1
    writeLines(mptp.scan[1:skiplines])
    mptp.raw <- readr::read_delim(file,skip=skiplines,delim=",",col_names="individual",show_col_types=FALSE)
    mptp.tab <- mptp.raw %>% 
        dplyr::mutate(species=ifelse(grepl(":",individual),individual,NA)) %>%
        tidyr::fill(species,.direction="down") %>%
        dplyr::filter(!grepl(":",individual)) %>%
        dplyr::mutate(species=stringr::str_replace_all(species,":| ","")) %>%
        dplyr::relocate(species,.before=individual)
    return(mptp.tab)
}


# run swarm
run_swarm <- function(file) {
    string.swarm <- paste0("swarm --threads 4 --differences 1 --fastidious --output-file ",file,".out ",file)
    system(command=string.swarm,ignore.stdout=FALSE)
    #return(string.swarm)
    swarm.tab <- readr::read_table(file=paste0(file,".out"),col_names=FALSE)
    swarm.tab.flipped <- swarm.tab %>% 
        dplyr::mutate(swarmCluster=paste0("Swarm",str_pad(seq_along(1:nrow(swarm.tab)),pad="0",width=4))) %>%
        dplyr::relocate(swarmCluster,.before=everything()) %>%
        tidyr::pivot_longer(cols=!swarmCluster,names_to="column",values_to="asvHash",values_drop_na=TRUE) %>% 
        dplyr::select(!column) %>%
        dplyr::mutate(asvHash=str_replace_all(asvHash,"_[0-9]+",""))
    return(swarm.tab.flipped)
}
