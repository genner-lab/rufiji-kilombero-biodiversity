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


############ RUN RAXML FOR 100 RANDOM STARTING TREES ############
############ RUN RAXML FOR 100 RANDOM STARTING TREES ############

# make a basename
base.name <- here("temp-local-only/mptp.fasta")

# write out fasta for alignment
tax.sub.filtered.fas <- tab2fas(df=tax.sub,seqcol="nucleotides",namecol="asvHash")
write.FASTA(tax.sub.filtered.fas,base.name)

# run raxml over 100 random starting trees
rand.tr <- raxml_random(file=base.name,epsilon=10)


############ ITERATE MPTP ############
############ ITERATE MPTP ############

# make a basename
base.name <- here("temp-local-only/mptp.fasta")

# get max cutoff at 0.01% of total reads
tot.reads <-  ceiling(sum(pull(tax.sub,nreads))*0.0001)

# generate distrib of values
seqs <- generate_thresholds(start=1,max=tot.reads)

# load up 100 ml trees
trs <- read.tree(paste0(base.name,".ali.raxml.rba.raxml.mlTrees"))

# run mptp in parallel over a list of trees
mptp.tab.all.trees <- mcmapply(function(x,y) mptp_parallel(df=tax.sub,base.name=base.name,tr=x,num=y,threshold="single",filt=seqs), x=trs, y=seq_along(trs), USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=8)
file.remove(list.files(here("temp-local-only"),pattern=".nwk",full.names=TRUE))

# join into dataframe
mptp.combined <- bind_rows(mptp.tab.all.trees)

# save dataframe
mptp.combined %>% write_csv(here("temp-local-only/delim-results.csv"))
#mptp.combined %>% write_csv(here("temp/results/delim-results.csv"))

# read back in if needed
mptp.combined <- read_csv(here("temp-local-only/delim-results.csv"),show_col_types=FALSE)

# get density
adj <- 1
d <- density(pull(mptp.combined,nClust),adjust=adj)
max.d <- ceiling(d$x[which.max(d$y)])
print(max.d)

# plot density
mptp.combined %>% 
    ggplot(aes(x=nClust)) + 
    geom_density(adjust=adj) + 
    geom_vline(xintercept=max.d,color="firebrick1",linetype=2) +
    scale_x_continuous(labels=scales::comma_format(accuracy=1),n.breaks=12) +
    ggthemes::theme_clean()

# plot gam 
mptp.combined %>% 
    ggplot(aes(x=filterThreshold,y=nClust)) + 
    geom_boxplot(aes(group=filterThreshold),outlier.shape=NA) +
    geom_jitter(height=0,alpha=0.1) + 
    geom_smooth(method=mgcv::gam,formula=y~s(x,bs="cs",k=7),alpha=0,color="seagreen") +
    scale_x_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=6) +
    geom_hline(yintercept=68,lty=2,color="cornflowerblue") +
    geom_hline(yintercept=max.d,lty=2,color="firebrick1") +
    #geom_hline(yintercept=median(pull(mptp.combined,nClust)),lty=2,color="forestgreen") +
    scale_y_continuous(labels=scales::comma_format(accuracy=1),n.breaks=12) + 
    #scale_y_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=12) + 
    ggthemes::theme_clean()
    #


############ PLOT A REPRESENTATIVE TREE ############
############ PLOT A REPRESENTATIVE TREE ############

#mptp.combined %>% filter(nClust >= max.d-1 & nClust <= max.d+1) %>% pull(filterThreshold) %>% median()
est.threshold <- mptp.combined %>% filter(nClust == max.d) %>% pull(filterThreshold) %>% median()
print(est.threshold)

# filter df at optimised threshold
tax.sub.opt <- tax.sub %>% filter(nreads>est.threshold)

# subset and write fasta
tax.sub.opt.fas <- tab2fas(df=tax.sub.opt,seqcol="nucleotides",namecol="asvHash")
write.FASTA(tax.sub.opt.fas,here("temp-local-only/opt.fasta"))

# run raxml
mptp.opt.tr <- raxml_ng(file=here("temp-local-only/opt.fasta"),epsilon=0.1)
mptp.opt.tr <- phangorn::midpoint(unroot(mptp.opt.tr))

# or load prev best ml tr
mptp.opt.tr <- read.tree(file=here("temp-local-only/mptp.fasta.ali.raxml.rba.raxml.bestTree"))
mptp.opt.tr <- phangorn::midpoint(unroot(mptp.opt.tr))
mptp.opt.tr <- drop_tips(df=tax.sub,filt=est.threshold,tree=mptp.opt.tr)


#min.edge.length <- min(mptp.opt.tr$edge.length)
min.edge.length <- 1e-07

# write out for mptp
write.tree(mptp.opt.tr,here("temp-local-only/opt.fasta.ali.raxml.rba.raxml.bestTree.nwk"))

# run mptp and read in
run_mptp(file=here("temp-local-only/opt.fasta.ali.raxml.rba.raxml.bestTree.nwk"),threshold="single",minbr=min.edge.length)
mptp.tab <- read_mptp(file=here("temp-local-only/opt.fasta.ali.raxml.rba.raxml.bestTree.nwk.mptp.out.txt"))
file.remove(list.files(here("temp-local-only"),pattern="opt.",full.names=TRUE))

# join with df
mptp.tab.annotated <- tax.sub.opt %>% left_join(mptp.tab)
mptp.tab.annotated %>% distinct(mptpDelim) %>% pull() %>% length()

# get come random colours
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
set.seed(42)
cols <- sample(getPalette(n=length(unique(pull(mptp.tab.annotated,mptpDelim)))))

# annotate and plot the mptp tree
p <- ggtree(mptp.opt.tr, ladderize=TRUE, color="grey20") + xlim(0,3)
p <- p %<+% mptp.tab.annotated
p <- p + geom_tiplab(aes(fill=mptpDelim,label=tipLabel),geom="label",label.size=0,label.padding=unit(0.1,"lines")) +
    theme(legend.position="none") +
    scale_fill_manual(values=cols)

# save plot as pdf
ggsave(plot=p, filename=here("temp-local-only/tree.pdf"), width=nrow(mptp.tab.annotated)/5, height=nrow(mptp.tab.annotated)/2.5, units="cm", bg="transparent", limitsize=FALSE)


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