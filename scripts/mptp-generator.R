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

# make a basename
base.name <- here("temp-local-only/mptp.fasta")


############ FILTER BY SWARM IF REQUIRED ############
############ FILTER BY SWARM IF REQUIRED ############

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


############ FILTER BY SEQUENCE LENGTH IF REQUIRED ############
############ FILTER BY SEQUENCE LENGTH IF REQUIRED ############

# set length min/max at 6%
tax.sub <- tax.sub %>% mutate(
        length=str_length(nucleotides),
        avgLength=median(length),
        minLength=avgLength-ceiling(avgLength*0.06),
        maxLength=avgLength+ceiling(avgLength*0.06)) %>%
    mutate(include=if_else(length >= minLength & length <= maxLength,TRUE,FALSE))

# filter
tax.sub <- tax.sub %>% filter(include==TRUE)


############ FILTER BY N EVENTS  ############
############ FILTER BY N EVENTS  ############

# get samples data
asvs.by.sample <- get_samples_table() 


# filter and count samples per asv
asvs.by.sample.n <- asvs.by.sample %>% 
    filter(asvHash %in% pull(tax.sub,asvHash)) %>% 
    filter(!grepl("Blank",eventID)) %>%
    group_by(asvHash) %>%
    summarise(nEvents=length(unique(sampleHash)))

# join with table
tax.sub <- tax.sub %>% left_join(asvs.by.sample.n,by="asvHash")

# filter by reads > 1 or events > 1
tax.sub <- tax.sub %>% filter(nEvents>1 | nreads>1) %>% arrange(desc(nreads))


############ RUN RAXML FOR 100 RANDOM STARTING TREES ############
############ RUN RAXML FOR 100 RANDOM STARTING TREES ############

        # write out fasta for alignment
        tax.sub.filtered.fas <- tab2fas(df=tax.sub,seqcol="nucleotides",namecol="asvHash")
        write.FASTA(tax.sub.filtered.fas,base.name)
        
        # run raxml over 100 random starting trees
        rand.tr <- raxml_random(file=base.name,epsilon=0.1)


############ ITERATE MPTP ############
############ ITERATE MPTP ############

# get max cutoff at 0.01% of total reads
tot.reads <-  ceiling(sum(pull(tax.sub,nreads))*0.0001)

# generate distrib of values
seqs <- generate_thresholds(start=1,max=tot.reads)

# load up 100 ml trees
trs <- read.tree(paste0(base.name,".ali.raxml.rba.raxml.mlTrees"))

# run mptp in parallel over a list of trees
mptp.tab.all.trees <- mcmapply(function(x,y) mptp_parallel(df=tax.sub,base.name=base.name,tr=x,num=y,threshold="single",filt=seqs,parsespp=TRUE), x=trs, y=seq_along(trs), USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=8)
file.remove(list.files(here("temp-local-only"),pattern=".nwk",full.names=TRUE))

# join into dataframe
mptp.combined <- bind_rows(mptp.tab.all.trees)

# save dataframe
mptp.combined %>% write_csv(here("temp-local-only/delim-results.csv"))
mptp.combined %>% write_csv(here("temp/results/delim-results.csv"))

# print
writeLines("mPTP analyses complete")

# read back in if needed
mptp.combined <- read_csv(here("temp-local-only/delim-results.csv"),show_col_types=FALSE)

# get density
adj <- 1
d <- density(pull(mptp.combined,nClust),adjust=adj)
max.d <- ceiling(d$x[which.max(d$y)])
print(max.d)

# get threshold estimate
#mptp.combined %>% filter(nClust >= max.d-1 & nClust <= max.d+1) %>% pull(filterThreshold) %>% density(adjust=0.1) %>% plot()
#mptp.combined %>% filter(nClust == max.d) %>% pull(filterThreshold) %>% density(adjust=0.1) %>% plot()
#est.threshold <- mptp.combined %>% filter(nClust >= max.d-1 & nClust <= max.d+1) %>% pull(filterThreshold) %>% median()
est.threshold <- mptp.combined %>% filter(nClust == max.d) %>% pull(filterThreshold) %>% median()
print(est.threshold)

# plot gam 
mptp.combined %>% 
    ggplot(aes(x=filterThreshold,y=nClust)) + 
    geom_boxplot(aes(group=filterThreshold),outlier.shape=NA) +
    geom_jitter(height=0,alpha=0.1) + 
    geom_smooth(method=mgcv::gam,formula=y~s(x,bs="cs",k=4),alpha=0,color="seagreen") +
    scale_x_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=6) +
    geom_hline(yintercept=68,lty=2,color="cornflowerblue") +
    geom_hline(yintercept=max.d,lty=2,color="firebrick1") +
    #geom_hline(yintercept=median(pull(mptp.combined,nClust)),lty=2,color="forestgreen") +
    scale_y_continuous(labels=scales::comma_format(accuracy=1),n.breaks=12) + 
    #scale_y_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=12) + 
    ggthemes::theme_clean()
    #

# plot density
mptp.combined %>% 
    ggplot(aes(x=nClust)) + 
    geom_density(adjust=adj) + 
    geom_vline(xintercept=max.d,color="firebrick1",linetype=2) +
    scale_x_continuous(labels=scales::comma_format(accuracy=1),n.breaks=12) +
    ggthemes::theme_clean()


############ GET FREQUENCY OF DELIMITATIONS ############
############ GET FREQUENCY OF DELIMITATIONS ############

# load up 100 ml trees
trs <- read.tree(paste0(base.name,".ali.raxml.rba.raxml.mlTrees"))

# run mptp in parallel over a list of trees
mptp.tab.filter.trees <- mcmapply(function(x,y) mptp_parallel(df=tax.sub,base.name=base.name,tr=x,num=y,threshold="single",filt=est.threshold,parsespp=FALSE), x=trs, y=seq_along(trs), USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=8)

# join into dataframe
mptp.combined <- bind_rows(mptp.tab.filter.trees)

# get hashes for sets
mptp.sets <- mptp.combined %>% 
    group_by(replicate,mptpDelim) %>% 
    arrange(asvHash,.by_group=TRUE) %>% 
    mutate(set=paste(asvHash,collapse=",")) %>% 
    ungroup() %>%
    mutate(setHash=md5(set))

# get frequencies of sets
mptp.sets.freq <- mptp.sets %>% 
    distinct(replicate,setHash) %>% 
    count(setHash) %>% 
    rename(setFreq=n) %>%
    mutate(setRank=row_number(desc(setFreq)))

# join sets and asvs
mptp.sets.joined <- mptp.sets %>% 
    left_join(mptp.sets.freq) %>% 
    distinct(asvHash,setHash,setFreq,setRank) %>% 
    arrange(setRank) %>% 
    group_by(asvHash) %>%
    mutate(asvFreq=length(unique(setHash))) %>%
    ungroup()

# reduce the assignment df
tax.sub.red <- tax.sub %>% select(asvHash,nreads,assigned,assignedName,tipLabel)

# join with mptp sets and annotate
mptp.sets.ann <- mptp.sets.joined %>% left_join(tax.sub.red)

# write out
mptp.sets.ann %>% write_csv(here("temp-local-only/delimitation-probabilities.csv"))
#mptp.sets.ann <- read_csv(here("temp-local-only/delimitation-probabilities.csv",show_col_types=FALSE))


############ MAKE A CONSENSUS TREE SET ############
############ MAKE A CONSENSUS TREE SET ############

# get all trees
tree.list <- mapply(function(x) read.tree(x), x=list.files(here("temp-local-only"),pattern=".nwk$",full.names=TRUE),USE.NAMES=FALSE,SIMPLIFY=FALSE)
file.remove(list.files(here("temp-local-only"),pattern=".nwk",full.names=TRUE))

# write out trees
class(tree.list) <- "multiPhylo"
write.nexus(tree.list,file=here("temp-local-only/all.trees.nex"),translate=FALSE)

# run treeannotator to get consensus tree
all.trees.nex <- tree_annotator(file=here("temp-local-only/all.trees.nex"),heights="keep")

# select asvs with highest set freq
mptp.sets.ann.sub <- mptp.sets.ann %>% 
    group_by(asvHash) %>% 
    slice_max(setFreq,with_ties=FALSE) %>%
    ungroup()
#mptp.sets.ann %>% arrange(asvHash,desc(setFreq)) %>% print(n=Inf)
#mptp.sets.ann.sub %>% distinct(setHash) %>% print(n=Inf)

# get come random colours
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
set.seed(42)
cols <- sample(getPalette(n=length(unique(pull(mptp.sets.ann.sub,setHash)))))
print(length(unique(pull(mptp.sets.ann.sub,setHash))))

# annotate and plot the mptp tree
p <- ggtree(all.trees.nex, ladderize=TRUE, color="grey20") + xlim(0,4)
p <- p %<+% mptp.sets.ann
p <- p + geom_tiplab(aes(label=setFreq), geom="text") + 
    geom_tiplab(aes(fill=setHash,label=tipLabel),geom="label",label.size=0,label.padding=unit(0.2,"lines"),offset=0.1) + 
    theme(legend.position="none") +
    scale_fill_manual(values=cols)

# save plot as pdf
ggsave(plot=p, filename=here("temp-local-only/tree.pdf"), width=length(all.trees.nex$tip.label)/4, height=length(all.trees.nex$tip.label)/2, units="cm", bg="transparent", limitsize=FALSE)






############ VSEARCH 3% CLUSTERING ############
############ VSEARCH 3% CLUSTERING ############

tax.sub.vsearch <- tax.sub %>% mutate(vsearchLabel=paste(asvHash,nreads,sep=";size="))
tax.sub.vsearch.fas <- tab2fas(df=tax.sub.vsearch,seqcol="nucleotides",namecol="vsearchLabel")
write.FASTA(tax.sub.vsearch.fas,here("temp-local-only/vsearch.fasta"))
# run in terminal
#vsearch --cluster_size vsearch.fasta --centroids vsearch.fasta.out --id 0.97 --sizein --sizeout
    #Clusters: 843 Size min 1, max 24, avg 1.5
    #Singletons: 685, 55.6% of seqs, 81.3% of clusters


############ MESSING WITH GAMS ############
############ MESSING WITH GAMS ############

library("mgcv")
library("gratia")
#https://fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/
m0 <- mgcv::gam(nClust~s(log10(filterThreshold),bs="cs",k=4),data=mptp.combined)#,bs="cs",sp=0.5
summary(m0)
plot(m0)
#gam.check(m0)
gratia::derivatives(m0,partial_match = TRUE) %>% print(n=Inf)
gratia::appraise(m0)






############ ZZZZ OLD ZZZZ ############
############ ZZZZ OLD ZZZZ ############


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


# 
root_trees <- function(base.name,tr,num,threshold,filt) {
    base.name.rep <- here(paste0(base.name,".tr",str_pad(as.character(num),pad="0",width=3),".nwk"))
    write.tree(ladderize(midpoint(tr)),base.name.rep) 
    mptp.tr <- read.tree(base.name.rep)
    min.edge.length <- min(mptp.tr$edge.length)
    run_mptp(file=base.name.rep,threshold=threshold,minbr=min.edge.length)
    mptp.tab <- read_mptp(file=paste0(base.name.rep,".mptp.out.txt"))
    mptp.species <- mptp.tab %>% 
        summarise(nClust=length(unique(mptpDelim))) %>% 
        mutate(filterThreshold=filt,replicate=num) %>% 
        relocate(filterThreshold,.before=nClust)
    return(mptp.species)
}

tr1.mptp <- root_trees(base.name=base.name,tr=trs[[1]],num=1,threshold="single",filt=1)

trs.mptp.list <- mcmapply(function(x,y) root_trees(base.name=base.name,tr=x,num=y,threshold="single",filt=1), x=trs, y=seq_along(trs), USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=4)

mptp.results <- bind_rows(trs.mptp.list)
mptp.results %>% print(n=Inf)

#plot(mptp.results$nClust)

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
    scale_y_continuous(trans="log10") + #labels=scales::comma,breaks=seq(0,2000,200)
    ggthemes::theme_clean()
    #
    #scale_y_continuous(trans="log10",labels=scales::comma,breaks=seq(0,2000,100))
    #scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=10))
    #scale_y_continuous(breaks=seq(0,1900,by=100))

# get n spp (= 68)
#tax.ass %>% 
#    filter(isFish==TRUE & isContaminated==FALSE & assigned==TRUE) %>%
#    distinct(assignedName)