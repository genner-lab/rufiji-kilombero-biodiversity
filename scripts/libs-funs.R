#!/usr/bin/env Rscript

# lib
library("here")
library("tidyverse")
library("lulu")
library("ape")
library("phangorn")
library("RColorBrewer")
library("ggtree")
library("parallel")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")
# https://github.com/tobiasgf/lulu
options(width=180)


# NJ TREE FUN
nj_og <- function(file) {
        string.mafft <- paste0("mafft --thread -1 --maxiterate 2 --retree 2 ",file," > ",file,".ali")
        system(command=string.mafft,ignore.stdout=FALSE)
        ali <- ape::read.FASTA(paste0(file,".ali"))
        rax.tr <- ape::njs(ape::dist.dna(ali,model="TN93",gamma=0.5,pairwise.deletion=TRUE))
    return(rax.tr)
}


# NEW RAXML-NG FUN
raxml_ng <- function(file,epsilon) {
        string.mafft <- paste0("mafft --thread -1 --maxiterate 2 --retree 2 ",file," > ",file,".ali")
        system(command=string.mafft,ignore.stdout=FALSE)
        string.parse <- paste0("raxml-ng --parse --msa ",file,".ali --model TN93+G --seed 42 --redo --threads auto")
        system(command=string.parse,ignore.stdout=FALSE)
        string.search <- paste0("raxml-ng --search --msa ",file,".ali.raxml.rba --tree pars{1} --lh-epsilon ",epsilon," --seed 42 --redo --threads auto")
        system(command=string.search,ignore.stdout=FALSE)
        rax.tr <- ape::read.tree(file=paste0(file,".ali.raxml.rba.raxml.bestTree"))
    return(rax.tr)
}

# RAXML RANDOM 100
raxml_random <- function(file,epsilon) {
        string.mafft <- paste0("mafft --thread -1 --maxiterate 2 --retree 2 ",file," > ",file,".ali")
        system(command=string.mafft,ignore.stdout=FALSE)
        string.parse <- paste0("raxml-ng --parse --msa ",file,".ali --model TN93+G --seed 42 --redo --threads auto")
        system(command=string.parse,ignore.stdout=FALSE)
        string.search <- paste0("raxml-ng --search --msa ",file,".ali.raxml.rba --tree rand{100} --lh-epsilon ",epsilon," --redo --threads auto")
        system(command=string.search,ignore.stdout=FALSE)
        rax.tr <- ape::read.tree(file=paste0(file,".ali.raxml.rba.raxml.bestTree"))
    return(rax.tr)
}


# PARSIMONY FUN
parsimony_ng <- function(file) {
        string.mafft <- paste0("mafft --thread -1 --maxiterate 2 --retree 2 ",file," > ",file,".ali")
        system(command=string.mafft,ignore.stdout=FALSE)
        string.parse <- paste0("raxml-ng --parse --msa ",file,".ali --model TN93+G --seed 42 --redo --threads auto")
        system(command=string.parse,ignore.stdout=FALSE)
        string.start <- paste0("raxml-ng --start --msa ",file,".ali.raxml.rba --tree pars{1} --seed 42 --redo --threads auto")
        system(command=string.start ,ignore.stdout=FALSE)
        string.search <- paste0("raxml-ng --evaluate --msa ",file,".ali.raxml.rba --tree ",file,".ali.raxml.rba.raxml.startTree --lh-epsilon 0.1 --seed 42 --redo --threads auto")
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
    mptp.raw <- readr::read_delim(file,skip=skiplines,delim=",",col_names="asvHash",show_col_types=FALSE)
    mptp.tab <- mptp.raw %>% 
        dplyr::mutate(mptpDelim=ifelse(grepl(":",asvHash),asvHash,NA)) %>%
        tidyr::fill(mptpDelim,.direction="down") %>%
        dplyr::filter(!grepl(":",asvHash)) %>%
        dplyr::mutate(mptpDelim=stringr::str_replace_all(mptpDelim,":","")) %>%
        dplyr::mutate(mptpDelim=stringr::str_replace_all(mptpDelim,"Species ","")) %>%
        dplyr::mutate(mptpDelim=paste0("mptp",str_pad(mptpDelim,pad="0",width=4))) %>%
        dplyr::relocate(mptpDelim,.before=asvHash)
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


# run mptp
iterate_mptp <- function(df,filter,threshold,epsilon) {
    #file.remove(list.files(here("temp-local-only"),pattern="mptp",full.names=TRUE))
    tax.sub.filtered <- df %>% 
        filter(nreads >= filter)
    tax.sub.filtered.fas <- tab2fas(df=tax.sub.filtered ,seqcol="nucleotides",namecol="asvHash")
    base.name <- paste0("temp-local-only/mptp.",str_pad(as.character(filter),pad="0",width=4),".fasta")
    write.FASTA(tax.sub.filtered.fas,here(base.name))
    #mptp.tr <- parsimony_ng(file=here(base.name))
    #mptp.tr <- raxml_ng(file=here(base.name),epsilon=epsilon)
    mptp.tr <- nj_og(file=here(base.name))
    write.tree(ladderize(midpoint(mptp.tr)),here(paste0(base.name,".nwk"))) 
    mptp.tr <- read.tree(here(paste0(base.name,".nwk")))
    min.edge.length <- min(mptp.tr$edge.length)
    run_mptp(file=here(paste0(base.name,".nwk")),threshold=threshold,minbr=min.edge.length)
    mptp.tab <- read_mptp(file=here(paste0(base.name,".nwk.mptp.out.txt")))
    mptp.species <- mptp.tab %>% 
        summarise(nClust=length(unique(mptpDelim))) %>% 
        mutate(filterThreshold=filter) %>% 
        relocate(filterThreshold,.before=nClust)
    #file.remove(list.files(here("temp-local-only"),pattern="mptp",full.names=TRUE))
    return(mptp.species)
}


# FUNCTION TO GENERATE FILTER THRESHOLD VALUES
generate_thresholds <- function(start,max) {
    i <- start 
    res <- NULL
    while(i < max) { 
        i <- i*1.5
        res <- c(res,i)
    }
    res <- ceiling(c(start,res))
    return(res)
}
