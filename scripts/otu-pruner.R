#!/usr/bin/env Rscript

# lib
library("here")
library("tidyverse")
library("lulu")
# https://github.com/tobiasgf/lulu
options(width=180)

# load the assigned taxonomy file
tax.ass <- read_csv(here("meta-fish-pipe/results/taxonomic-assignments.csv"))


# load master asv table
asvs.all <- read_csv(file=here("meta-fish-pipe/temp/taxonomic-assignment/asvs-all.csv"))

# make id
asvs.all %<>% mutate(asvCode=paste(primer,lib,asv,sep="|")) %>% 
    select(-asv) %>% 
    rename(asvHash=md5,primerSet=primer,library=lib)

# get dirs for otu tables
libs.dirs <- list.dirs(path=here("meta-fish-pipe/temp/processing"),full.names=TRUE,recursive=FALSE)

# load up all the otu tables
tables.list <- lapply(here(libs.dirs,"results/asv-table.tsv"),read_tsv)

# pivot long and merge otu tables
asv.tables.df <- purrr::map_dfr(tables.list,pivot_longer,cols=-asv,names_to="sampleHash",values_to="nreads") %>% 
    filter(nreads > 0) %>% 
    rename(asvCode=asv)

# load up table of events
events.list <- lapply(here(libs.dirs,"results/events-hashes.csv"),read_csv)

# merge otu tables
events.df <- bind_rows(events.list) %>% rename(sampleHash=hashLabel)

# join the codes table and events
asv.tables.events <- asv.tables.df %>% left_join(events.df,by="sampleHash")

# join with the hash table 
asvs.by.sample <- asv.tables.events %>% left_join(asvs.all,by=c("asvCode","primerSet","library")) 

# filter the taxonomic assignment table
tax.sub <- tax.ass %>% 
    filter(isFish==TRUE & isContaminated==FALSE) %>% 
    mutate(assignedName=if_else(assigned==FALSE,sintaxSpeciesID,assignedName)) %>%
    select(asvHash,nreads,isFish,assigned,assignedName,nucleotides) %>% 
    mutate(label=paste(assignedName,assigned,asvHash),label=str_replace_all(label," ","_"))

# filter the asv list
asvs.by.sample %>% 
    filter(asvHash %in% pull(tax.sub,asvHash)) %>% 
    select(sampleHash,nreads,asvHash) %>% 
    left_join(select(tax.sub,asvHash,label)) %>%
    select(-asvHash) %>%
    pivot_wider(names_from=sampleHash,values_from=nreads,values_fill=0)
