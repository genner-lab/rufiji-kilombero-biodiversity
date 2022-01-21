#!/usr/bin/env Rscript

# lib
library("here")
library("tidyverse")
library("lulu")
library("ape")
library("phangorn")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")
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
    select(asvHash,nreads,isFish,assigned,assignedName,blastSpeciesID,blastPident,nucleotides)# %>% 
    #mutate(label=paste(assignedName,assigned,asvHash),label=str_replace_all(label," ","_"))

# write out asv table
asvs.by.sample %>% 
    filter(asvHash %in% pull(tax.sub,asvHash)) %>% 
    filter(!grepl("Blank",eventID)) %>%
    select(sampleHash,nreads,asvHash) %>%
    pivot_wider(names_from=sampleHash,values_from=nreads,values_fill=0) %>% 
    write_csv(here("temp-local-only/asvs-table-lulu.csv"))


# write out fasta (need to use hash because max length for blast is 50 chars)
tab2fas(df=tax.sub,seqcol="nucleotides",namecol="asvHash") %>% write.FASTA(here("temp-local-only/asvs-lulu.fasta"))

# run blast to generate matchlist
system("./scripts/blast-match-list.sh")

# load otu table using base R!
asv.table <- read.csv(here("temp-local-only/asvs-table-lulu.csv"),sep=",",header=TRUE,as.is=TRUE,row.names=1)
#head(asv.table)

# load blast table using base R!
match.list <- read.table(here("temp-local-only/lulu-blast-matchlist.tsv"),sep="\t",header=FALSE,as.is=TRUE,stringsAsFactors=FALSE)
#head(match.list)

# run lulu
#?lulu
#,minimum_ratio_type="min",minimum_ratio=1,minimum_match=84,minimum_relative_cooccurence=0.95
lulu.result <- lulu(otutable=asv.table,matchlist=match.list)

# format lulu result - add identities, and assignments
lulu.result.by.asv <- lulu.result$otu_map %>% 
    rownames_to_column(var="asvHash") %>% 
    tibble() %>%
    left_join(rename(tibble(match.list),asvHash=V1,parent_id=V2,parentIdentity=V3)) %>%
    left_join(rename(select(tax.sub,asvHash,assignedName,assigned),parent_id=asvHash,parentName=assignedName,parentAssigned=assigned))

# join lulu result with assignment table
tax.sub.lulu  <- tax.sub %>% 
    left_join(lulu.result.by.asv) %>%
    select(asvHash,nreads,assigned,assignedName,blastSpeciesID,blastPident,parent_id,parentName,parentAssigned,curated,parentIdentity,rank,spread) %>%
    arrange(parentName,desc(nreads))


# write
tax.sub.lulu %>% write_csv(here("temp-local-only/taxonomic-assignments-lulu.csv"))
#read back in 
tax.sub.lulu <- read_csv(here("temp-local-only/taxonomic-assignments-lulu.csv"))

# get number parent otus
tax.sub.lulu %>% filter(curated=="parent") %>% distinct(asvHash)

# subset the merged but already assigned spp - where lulu is WRONG
tax.sub.lulu %>% 
    filter(curated=="merged" & parentName!=assignedName & blastPident>parentIdentity) %>% #assigned==TRUE & 
    arrange(desc(assigned),parent_id,desc(nreads)) %>%
    #write_csv(here("temp-local-only/taxonomic-assignments-lulu-merged.csv"))
    print(n=Inf)

# subset the merged but same assigned spp - where lulu is CORRECT
tax.sub.lulu %>% 
    filter(curated=="merged" & parentName==assignedName) %>% #assigned==TRUE & 
    #arrange(desc(assigned),parent_id,desc(nreads)) %>%
    #write_csv(here("temp-local-only/taxonomic-assignments-lulu-merged.csv"))
    print(n=Inf)


# get assigned names
lulu.assigned.spp <- tax.sub.lulu %>% 
    filter(curated=="parent" & assigned==TRUE) 

# get unassigned names
lulu.unassigned.spp <- tax.sub.lulu %>% 
    filter(curated=="parent" & assigned==FALSE)

# print unique names of assigned unassigned
lulu.assigned.spp %>% 
    distinct(parentName) %>% 
    pull(parentName)
lulu.unassigned.spp %>%
    distinct(parentName) %>%
    pull(parentName)

# 
setdiff(lulu.unassigned.spp %>%
    distinct(parentName) %>%
    pull(parentName),
    lulu.assigned.spp %>% 
    distinct(parentName) %>% 
    pull(parentName))

# pull out unassigned for blast
unassigned.fas <- tax.sub %>% filter(asvHash %in% pull(lulu.unassigned.spp,asvHash)) %>% mutate(label=paste(asvHash,assignedName,nreads,sep="_"))
tab2fas(df=unassigned.fas,seqcol="nucleotides",namecol="label") %>% write.FASTA(here("temp-local-only/unassigned-blast.fasta"))


# correct the incorrect merges and filter
tax.sub.lulu.filtered <- tax.sub.lulu %>% 
    mutate(luluMerge=curated) %>%
    mutate(luluMerge=if_else(curated=="merged" & parentName!=assignedName & blastPident>parentIdentity,"corrected",luluMerge)) %>% 
    filter(luluMerge=="parent" | luluMerge=="corrected") %>%
    filter(nreads>1) %>%
    filter(blastPident<100)# %>% #assigned!=FALSE | 


parents <- tax.sub %>% filter(asvHash %in% pull(tax.sub.lulu.filtered,asvHash)) 
parents.fas <- tab2fas(df=parents ,seqcol="nucleotides",namecol="asvHash")

# load reflib
refs.fas <- read.FASTA(here("meta-fish-pipe/temp/taxonomic-assignment/custom-reference-library.fasta"))
#refs.csv <- read_csv(here("meta-fish-pipe/temp/taxonomic-assignment/custom-reference-library.csv"))

# join 
write.FASTA(c(refs.fas,parents.fas),here("temp-local-only/mptp-parents.fasta"))


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


# run raxml
mptp.tr <- raxml_ng(file=here("temp-local-only/mptp-parents.fasta"))
#mptp.tr <- ape::read.tree(file=paste0(here("temp-local-only/mptp-parents.fasta"),".ali.raxml.rba.raxml.bestTree"))

plot(ladderize(midpoint(mptp.tr)),cex=0.001)

# write out rooted tree
write.tree(ladderize(midpoint(mptp.tr)),here("temp-local-only/mptp-parents.nwk")) 



# run mptp
mptp <- function(file) {
        string.mptp <- paste0("mptp --ml --single --tree_file ",file," --output_file ",file,".mptp.out")
        #  --minbr 0.0001
        #system("export PATH=~/Software/mptp/bin:$PATH")
        system(command=string.mptp,ignore.stdout=FALSE)
        #x <- x
    #return(string.mptp)
}

mptp(file=here("temp-local-only/mptp-parents.nwk"))


## mptp --ml --single --minbr 0.0001 --tree_file ml.haps.tr.nwk --output_file ml.haps.tr.nwk.out









##################################################
lulu.result.by.asv <- lulu.result$otu_map %>% 
    rownames_to_column(var="asvHash") %>% 
    tibble() %>%
    mutate(assignedName=pull(tax.sub,assignedName)[match(asvHash,pull(tax.sub,asvHash))]) %>%
    mutate(assigned=pull(tax.sub,assigned)[match(asvHash,pull(tax.sub,asvHash))]) %>%
    mutate(parentName=pull(tax.sub,assignedName)[match(parent_id,pull(tax.sub,asvHash))]) %>%
    mutate(parentAssigned=pull(tax.sub,assigned)[match(parent_id,pull(tax.sub,asvHash))]) %>%
    select(asvHash,assignedName,assigned,total,spread,parent_id,parentName,parentAssigned,curated,rank)

#
lulu.result.by.asv %>% 
    left_join(rename(tibble(match.list),asvHash=V1,parent_id=V2)) %>% 
    filter(assignedName!=parentName) %>% 
    arrange(parentName,desc(total)) %>% 
    print(n=Inf)


rename(tibble(match.list),asvHash=V1,parent_id=V2)



tax.sub %>% 



# load up matchlist
lulu.matchlist <- read_tsv(here("temp-local-only/lulu-blast-matchlist.tsv"),col_names=c("asv1","asv2","identity"))

# add names back 
lulu.matchlist.names <- lulu.matchlist %>% 
    mutate(label1=pull(tax.sub.names,label)[match(asv1,pull(tax.sub,asvHash))],
        label2=pull(tax.sub.names,label)[match(asv2,pull(tax.sub,asvHash))]) %>%
    select(-asv1,-asv2) %>%
    relocate(identity,.after=last_col())


lulu.result <- lulu(otutable=data.frame(asv.table.fish,row.names=1),matchlist=lulu.matchlist.names)

match.list %>% filter(V1=="d2fa7148466d7bf1678000c60330e05c" & V2=="ae4318a5d24d1264b3524c5b1642abf0")
match.list %>% filter(V2=="d2fa7148466d7bf1678000c60330e05c" & V1=="ae4318a5d24d1264b3524c5b1642abf0")


c5e48c62f125ef22221d4f236b7b65af
ae4318a5d24d1264b3524c5b1642abf0


 %>% filter(V2=="ae4318a5d24d1264b3524c5b1642abf0")
 %>% filter(V2=="cf14519ed8e5fa7af73e942989bce49e")
str(lulu.result)

lulu.result$curated_table
lulu.result$curated_count
lulu.result$otu_map
lulu.result$curated_otus
lulu.result$discarded_count
lulu.result$minimum_match
lulu.result$minimum_relative_cooccurence

lulu.result$discarded_otus


# make an annotated raw table combining assigments
lulu.result$curated_table %>% rownames_to_column(var="mother") %>% as_tibble() %>% 
    mutate(size=pull(tax.ass.df,size)[match(mother,pull(tax.ass.df,md5))],
    bestId=pull(tax.ass.df,bestId)[match(mother,pull(tax.ass.df,md5))],
    matchLength=pull(tax.ass.df,matchLength)[match(mother,pull(tax.ass.df,md5))],
    identity=pull(tax.ass.df,identity)[match(mother,pull(tax.ass.df,md5))],
    blastId=pull(tax.ass.df,blastId)[match(mother,pull(tax.ass.df,md5))]) %>% 
    select(mother,size,bestId,matchLength,identity,blastId,everything()) %>%
    arrange(bestId,desc(size)) %>%
    write_csv(path="results/otu-table-raw-annotated-lulu.csv")

#       1. ‘curated_table’ - a curated OTU table with daughters merged
#          with their matching parents.
#
#       2. ‘curated_count’ - number of curated (parent) OTUs.
#
#       3. ‘curated_otus’ - ids of the OTUs that were accepted as valid
#          OTUs.
#
#       4. ‘discarded_count’ - number of discarded (merged with parent)
#          OTUs.
#
#       5. ‘discarded_otus’ - ids of the OTUs that were identified as
#          errors (daughters) and merged with respective parents.
#
#       6. ‘runtime’ - time used by the script.
#
#       7. ‘minimum_match’ - the id threshold (minimum match % between
#          parent and daughter) for evaluating co-occurence (set by
#          user).
#
#       8. ‘minimum_relative_cooccurence’ - minimum ratio of
#          daughter-occurences explained by co-occurence with parent
#          (set by user).
#
#       9. ‘otu_map’ - information of which daughters were mapped to
#          which parents.
#
#      10. ‘original_table’ - original OTU table.#


read_csv(here("assets/local-12s.csv")) %>% distinct(sciNameValid)
