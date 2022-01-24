#!/usr/bin/env Rscript

# libs and funs
source(here("scripts/libs-funs.R"))

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
refs.csv <- read_csv(here("meta-fish-pipe/temp/taxonomic-assignment/custom-reference-library.csv"))

# join 
write.FASTA(c(refs.fas,parents.fas),here("temp-local-only/mptp-parents.fasta"))


# run raxml
mptp.tr <- raxml_ng(file=here("temp-local-only/mptp-parents.fasta"))
#mptp.tr <- ape::read.tree(file=paste0(here("temp-local-only/mptp-parents.fasta"),".ali.raxml.rba.raxml.bestTree"))
plot(ladderize(midpoint(mptp.tr)),cex=0.001)

# write out rooted tree
write.tree(ladderize(midpoint(mptp.tr)),here("temp-local-only/mptp-parents.nwk")) 


# run mptp
mptp(file=here("temp-local-only/mptp-parents.nwk"))


## mptp --ml --single --minbr 0.0001 --tree_file ml.haps.tr.nwk --output_file ml.haps.tr.nwk.out

# run func
mptp.tab <- read_mptp(file=here("temp-local-only/mptp-parents.nwk.mptp.out.txt"),skiplines=6)

# 
tax.sub.lulu.filtered.red <- tax.sub.lulu.filtered %>% 
    mutate(source="ASV") %>% 
    distinct(asvHash,nreads,assigned,parentName,source) %>% 
    mutate(label=paste(asvHash,assigned,parentName,source,nreads,sep="_"),label=str_replace_all(label," ","_")) %>% 
    rename(individual=asvHash,assignedSpecies=parentName)

refs.csv.red <- refs.csv %>% 
    distinct(source,dbid,sciNameValid) %>% 
    mutate(label=paste(dbid,sciNameValid,source,sep="_"),label=str_replace_all(label," ","_")) %>% 
    rename(individual=dbid,assignedSpecies=sciNameValid)

# annotate the mptp table
mptp.tab.annotated <- mptp.tab %>% 
    left_join(bind_rows(tax.sub.lulu.filtered.red,refs.csv.red)) %>%
    relocate(species,.after=label)


# load up tree
mptp.tr <- read.tree(here("temp-local-only/mptp-parents.nwk"))

#mptp.tr$tip.label <- pull(mptp.tab.annotated,label)[match(mptp.tr$tip.label,pull(mptp.tab.annotated,individual))]

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
set.seed(9)
cols <- sample(getPalette(n=length(unique(pull(mptp.tab.annotated,species)))))


setdiff(mptp.tr$tip.label,pull(mptp.tab.annotated,individual))
setdiff(pull(mptp.tab.annotated,individual),mptp.tr$tip.label)


# plot the tree
p <- ggtree(mptp.tr, ladderize=TRUE, color="grey20", size=0.8)
p <- p %<+% select(mptp.tab.annotated,-label) + 
    geom_tiplab(aes(fill=species),geom="label",size=3) #+
    #scale_fill_manual(values=cols)
plot(p)

ggsave(plot=p, filename="../temp/delim_all2.pdf", width=14, height=30, bg="transparent", limitsize=FALSE)

