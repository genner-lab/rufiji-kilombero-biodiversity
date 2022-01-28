#!/usr/bin/env Rscript

############ LOAD LIBS AND FUNS ############
############ LOAD LIBS AND FUNS ############

source(here::here("scripts/libs-funs.R"))


############ LOAD ASSIGNMENT TABLE ############
############ LOAD ASSIGNMENT TABLE ############

# load the assigned taxonomy file
tax.ass <- read_csv(here("meta-fish-pipe/results/taxonomic-assignments.csv"))

# filter the taxonomic assignment table
tax.sub <- tax.ass %>% 
    filter(isFish==TRUE & isContaminated==FALSE) %>% 
    mutate(assignedName=if_else(assigned==FALSE,sintaxSpeciesID,assignedName)) %>%
    mutate(tiplabel=paste(str_trunc(asvHash,width=10,side="right",ellipsis=""),assignedName,assigned,nreads,sep="_"),tiplabel=str_replace_all(tiplabel," ","_")) %>%
    select(asvHash,nreads,assigned,assignedName,blastSpeciesID,blastPident,nucleotides,tiplabel)


############ READ ABUNDANCE FILTER ############
############ READ ABUNDANCE FILTER ############


# run mptp
iterate_mptp <- function(df,threshold) {
    #file.remove(list.files(here("temp-local-only"),pattern="mptp",full.names=TRUE))
    tax.sub.filtered <- df %>% 
        filter(nreads >= threshold)
    tax.sub.filtered.fas <- tab2fas(df=tax.sub.filtered ,seqcol="nucleotides",namecol="asvHash")
    base.name <- paste0("temp-local-only/mptp.",str_pad(as.character(threshold),pad="0",width=4),".fasta")
    write.FASTA(tax.sub.filtered.fas,here(base.name))
    mptp.tr <- parsimony_ng(file=here(base.name))
    write.tree(ladderize(midpoint(mptp.tr)),here(paste0(base.name,".nwk"))) 
    mptp.tr <- read.tree(here(paste0(base.name,".nwk")))
    min.edge.length <- min(mptp.tr$edge.length)
    run_mptp(file=here(paste0(base.name,".nwk")),threshold="single",minbr=min.edge.length)
    mptp.tab <- read_mptp(file=here(paste0(base.name,".nwk.mptp.out.txt")))
    mptp.species <- mptp.tab %>% 
        summarise(nClust=length(unique(mptpDelim))) %>% 
        mutate(filterThreshold=threshold) %>% 
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

# generate distrib of values
seqs <- generate_thresholds(start=1,max=100000)


# iterate mptp over values
mptp.results.list <- mcmapply(function(x) iterate_mptp(df=tax.sub,threshold=x), x=seqs, USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=2)
#file.remove(list.files(here("temp-local-only"),pattern="mptp",full.names=TRUE))

# join results
mptp.results <- bind_rows(mptp.results.list)
mptp.results %>% print(n=Inf)

# plot
mptp.results %>% ggplot(aes(x=filterThreshold,y=nClust)) + geom_point() + geom_line() + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
mptp.results %>% write_csv(here("temp-local-only/delim-results.csv"))


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
