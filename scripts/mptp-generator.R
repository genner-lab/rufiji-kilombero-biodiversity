#!/usr/bin/env Rscript

############ LOAD LIBS AND FUNS ############
############ LOAD LIBS AND FUNS ############

source(here::here("scripts/libs-funs.R"))


############ LOAD ASSIGNMENT TABLE ############
############ LOAD ASSIGNMENT TABLE ############

# load the assigned taxonomy file
tax.ass <- read_csv(here("meta-fish-pipe/results/taxonomic-assignments.csv"),show_col_types=FALSE)

# filter the taxonomic assignment table
tax.sub <- tax.ass %>% 
    filter(isFish==TRUE & isContaminated==FALSE) %>% 
    mutate(assignedName=if_else(assigned==FALSE,sintaxSpeciesID,assignedName)) %>%
    mutate(tipLabel=paste(str_trunc(asvHash,width=10,side="right",ellipsis=""),assignedName,assigned,nreads,sep="_"),tiplabel=str_replace_all(tipLabel," ","_")) %>%
    mutate(swarmLabel=paste(asvHash,nreads,sep="_")) %>%
    select(asvHash,nreads,assigned,assignedName,blastSpeciesID,blastPident,nucleotides,tipLabel,swarmLabel)


############ ADD SWARM IF REQUIRED ############
############ ADD SWARM IF REQUIRED ############

#write out fasta
swarm.fas <- tab2fas(df=tax.sub,seqcol="nucleotides",namecol="swarmLabel")
write.FASTA(swarm.fas,here("temp-local-only/swarm.fasta"))

# run swarm
swarm.result <- run_swarm(here("temp-local-only/swarm.fasta"))

# join with mptp table
tax.sub <- tax.sub %>% left_join(swarm.result)

# filter by swarm cluster with greatest number seqs
tax.sub <- tax.sub %>% 
    group_by(swarmCluster) %>% 
    slice_max(order_by=nreads,n=1,with_ties=FALSE) %>% 
    ungroup()


############ READ ABUNDANCE FILTER ############
############ READ ABUNDANCE FILTER ############


# generate distrib of values
seqs <- generate_thresholds(start=1,max=100000)

# iterate mptp over values
mptp.results.list <- mcmapply(function(x) iterate_mptp(df=tax.sub,filter=x,threshold="single",epsilon=10), x=seqs, USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=2)
#file.remove(list.files(here("temp-local-only"),pattern="mptp",full.names=TRUE))

# join results
mptp.results <- bind_rows(mptp.results.list)
mptp.results %>% print(n=Inf)

# writeout
mptp.results %>% 
    mutate(method="single-njs-swarm") %>%
    relocate(method,.before=filterThreshold) %>% 
    write_csv(here("temp-local-only/delim-results-temp.csv"))

# plot
mptp.table <- read_csv(here("temp-local-only/delim-results.csv"),show_col_types=FALSE)
mptp.table %>% 
    ggplot(aes(x=filterThreshold,y=nClust,color=method)) + 
    geom_point() + 
    geom_line() + 
    scale_x_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=6) +
    geom_hline(yintercept=68,lty=2) +
    scale_y_continuous(labels=scales::comma,breaks=seq(0,2000,200)) +
    ggthemes::theme_clean()
    #
    #scale_y_continuous(trans="log10",labels=scales::comma,breaks=seq(0,2000,100))
    #scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=10))
    #scale_y_continuous(breaks=seq(0,1900,by=100))

# get n spp (= 68)
#tax.ass %>% 
#    filter(isFish==TRUE & isContaminated==FALSE & assigned==TRUE) %>%
#    distinct(assignedName)


# run raxml random
tax.sub.filtered.fas <- tab2fas(df=tax.sub ,seqcol="nucleotides",namecol="asvHash")
base.name <- paste0("temp-local-only/mptp.",str_pad(as.character(1),pad="0",width=4),".fasta")
write.FASTA(tax.sub.filtered.fas,here(base.name))
# run
rand.tr <- raxml_random(file=here(base.name),epsilon=10)




#####################################################################################################################
# filter on nreads
# NEED TO ADD SPREAD
tax.sub.filtered <- tax.sub %>% 
    filter(nreads > 10000) #  & spread >= 1

# make fasta
tax.sub.filtered.fas <- tab2fas(df=tax.sub.filtered ,seqcol="nucleotides",namecol="asvHash")

# write out for mPTP
write.FASTA(tax.sub.filtered.fas,here("temp-local-only/mptp-filtered.fasta"))

# run raxml
mptp.tr <- parsimony_ng(file=here("temp-local-only/mptp-filtered.fasta"))

# write out rooted tree for mptp
write.tree(ladderize(midpoint(mptp.tr)),here("temp-local-only/mptp-filtered.nwk")) 

# read tree back in
mptp.tr <- read.tree(here("temp-local-only/mptp-filtered.nwk"))

# get distribution of edges
min.edge.length <- min(mptp.tr$edge.length)

# run mptp
run_mptp(file=here("temp-local-only/mptp-filtered.nwk"),threshold="single",minbr=min.edge.length)

# read in mPTP data
mptp.tab <- read_mptp(file=here("temp-local-only/mptp-filtered.nwk.mptp.out.txt"))

# get num spp
mptp.species <- mptp.tab %>% summarise(nClust=length(unique(mptpDelim))) %>% pull()
return(mptp.species)
