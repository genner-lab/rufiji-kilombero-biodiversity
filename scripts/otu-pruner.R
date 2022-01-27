#!/usr/bin/env Rscript

# libs and funs
source(here::here("scripts/libs-funs.R"))


############ LOAD ASSIGNMENTS ############
############ LOAD ASSIGNMENTS ############

# load the assigned taxonomy file
tax.ass <- read_csv(here("meta-fish-pipe/results/taxonomic-assignments.csv"))

# filter the taxonomic assignment table
tax.sub <- tax.ass %>% 
    filter(isFish==TRUE & isContaminated==FALSE) %>% 
    mutate(assignedName=if_else(assigned==FALSE,sintaxSpeciesID,assignedName)) %>%
    select(asvHash,nreads,isFish,assigned,assignedName,blastSpeciesID,blastPident,nucleotides)


############ RUN LULU ASVS ############
############ RUN LULU ASVS ############

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

# write out lulu result
tax.sub.lulu %>% write_csv(here("temp-local-only/taxonomic-assignments-lulu.csv"))


############ LULU RESULTS ############
############ LULU RESULTS ############

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

# compare
setdiff(lulu.unassigned.spp %>%
    distinct(parentName) %>%
    pull(parentName),
    lulu.assigned.spp %>% 
    distinct(parentName) %>% 
    pull(parentName))

# pull out unassigned to blast manually
unassigned.fas <- tax.sub %>% filter(asvHash %in% pull(lulu.unassigned.spp,asvHash)) %>% mutate(label=paste(asvHash,assignedName,nreads,sep="_"))
tab2fas(df=unassigned.fas,seqcol="nucleotides",namecol="label") %>% write.FASTA(here("temp-local-only/unassigned-blast.fasta"))


############ PREP MPTP ############
############ PREP MPTP ############

#read back in 
tax.sub.lulu <- read_csv(here("temp-local-only/taxonomic-assignments-lulu.csv"))

# correct the incorrect merges and filter
tax.sub.lulu.filtered <- tax.sub.lulu %>% 
    mutate(luluMerge=curated) %>%
    mutate(luluMerge=if_else(curated=="merged" & parentName!=assignedName & blastPident>parentIdentity,"corrected",luluMerge))

# get distrib of reads for filter ABOVE
tax.sub.lulu.filtered %>% filter(assigned==FALSE) %>% pull(nreads) %>% quantile(probs=0.9,names=FALSE)

tax.sub.lulu.filtered <- tax.sub.lulu.filtered %>% 
    filter(luluMerge=="parent" | luluMerge=="corrected") %>%
    filter(nreads>20 & spread >= 1)

# get nucleotides of filtered asvs
parents <- tax.sub %>% filter(asvHash %in% pull(tax.sub.lulu.filtered,asvHash)) 
parents.fas <- tab2fas(df=parents ,seqcol="nucleotides",namecol="asvHash")

# load references
refs.fas <- read.FASTA(here("meta-fish-pipe/temp/taxonomic-assignment/custom-reference-library.fasta"))
refs.csv <- read_csv(here("meta-fish-pipe/temp/taxonomic-assignment/custom-reference-library.csv"))

# remove the references that had no assignments made to it
refs.csv.ass <- refs.csv %>% filter(sciNameValid %in% unique(pull(filter(tax.sub.lulu.filtered,assigned==TRUE),assignedName)))
refs.fas.red <- refs.fas[names(refs.fas) %in% pull(refs.csv.ass,dbid)]

# join and write out for mPTP
write.FASTA(parents.fas,here("temp-local-only/mptp-parents.fasta"))#c(refs.fas.red,parents.fas)

# run raxml
mptp.tr <- raxml_ng(file=here("temp-local-only/mptp-parents.fasta"))
# raxml-ng --evaluate --msa mptp-parents.fasta.ali --tree mptp-parents.fasta.ali.raxml.rba.raxml.startTree --model TN93+G --seed 42 --redo --threads auto
# mptp.tr <- ape::read.tree(file=paste0(here("temp-local-only/rax-test/mptp-parents.fasta"),".ali.raxml.bestTree"))
#plot(ladderize(midpoint(mptp.tr)),cex=0.001)

# write out rooted tree for mptp
write.tree(ladderize(midpoint(mptp.tr)),here("temp-local-only/mptp-parents.nwk")) 

# read tree back in
mptp.tr <- read.tree(here("temp-local-only/mptp-parents.nwk"))

# get distribution of edges
min.edge.length <- min(mptp.tr$edge.length)

# run mptp
run_mptp(file=here("temp-local-only/mptp-parents.nwk"),threshold="single",minbr=min.edge.length)


############ PLOT MPTP ############
############ PLOT MPTP ############

# read in mPTP data
mptp.tab <- read_mptp(file=here("temp-local-only/mptp-parents.nwk.mptp.out.txt"))

# make asv table with tip labels
tax.sub.lulu.filtered.red <- tax.sub.lulu.filtered %>% 
    mutate(source="ASV") %>% 
    distinct(asvHash,nreads,assigned,assignedName,source) %>% 
    mutate(tiplabel=paste(str_trunc(asvHash,width=10,side="right",ellipsis=""),assignedName,assigned,nreads,sep="_"),tiplabel=str_replace_all(tiplabel," ","_")) %>% 
    rename(individual=asvHash)

# make references table with tip labels
refs.csv.red <- refs.csv %>% 
    distinct(source,dbid,sciNameValid) %>% 
    mutate(tiplabel=paste(dbid,sciNameValid,sep="_"),tiplabel=str_replace_all(tiplabel," ","_")) %>% 
    rename(individual=dbid,assignedName=sciNameValid)

# annotate the mptp table
mptp.tab.annotated <- mptp.tab %>% 
    left_join(bind_rows(tax.sub.lulu.filtered.red,refs.csv.red)) %>%
    relocate(species,.after=tiplabel)


############ RUN SWARM ############
############ RUN SWARM ############

# write out fasta with abundances
parents.swarm <- parents %>% mutate(label=paste(asvHash,nreads,sep="_"))

#write out fasta
swarm.fas <- tab2fas(df=parents.swarm,seqcol="nucleotides",namecol="label")
write.FASTA(swarm.fas,here("temp-local-only/swarm.fasta"))

# run swarm
swarm.result <- run_swarm(here("temp-local-only/swarm.fasta"))

# join with mptp table
mptp.tab.annotated <- mptp.tab.annotated %>% left_join(rename(swarm.result,individual=asvHash))

############ PLOT TREE ############
############ PLOT TREE ############

# load up tree
mptp.tr <- read.tree(here("temp-local-only/mptp-parents.nwk"))

# get come random colours
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
set.seed(42)
cols <- sample(getPalette(n=length(unique(pull(mptp.tab.annotated,swarmCluster)))))

# annotate and plot the mptp tree
p <- ggtree(mptp.tr, ladderize=TRUE, color="grey20") + xlim(0,3)
p <- p %<+% mptp.tab.annotated
p <- p + geom_tiplab(aes(fill=swarmCluster,label=tiplabel),geom="label",label.size=0,label.padding=unit(0.1,"lines")) +
    geom_tippoint(aes(shape=source,color=source),size=3) +
    theme(legend.position="none") +
    scale_fill_manual(values=cols)

# save plot as pdf
ggsave(plot=p, filename=here("temp-local-only/tree.pdf"), width=nrow(mptp.tab.annotated)/30, height=nrow(mptp.tab.annotated)/2.5, units="cm", bg="transparent", limitsize=FALSE)

#nrow(mptp.tab.annotated)/2.5
