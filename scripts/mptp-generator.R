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


# check
#tax.sub %>% filter(include==FALSE) %>% print(n=Inf)

# load up tree
#mptp.tr <- midpoint(read.tree(here("temp-local-only/mptp.0001.fasta.ali.raxml.rba.raxml.bestTree")))

# get come random colours
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#set.seed(42)
#cols <- sample(getPalette(n=length(unique(pull(tax.sub,include)))))

# annotate and plot the mptp tree
#p <- ggtree(mptp.tr, ladderize=TRUE, color="grey20") + xlim(0,4)
#p <- p %<+% tax.sub
#p <- p + geom_tiplab(aes(fill=include,label=tipLabel),geom="label",label.size=0,label.padding=unit(0.1,"lines")) +
#        theme(legend.position="none") +
#    scale_fill_manual(values=cols)

# save plot as pdf
#ggsave(plot=p, filename=here("temp-local-only/length-tree.pdf"), width=nrow(tax.sub)/30, height=nrow(tax.sub)/2.5, units="cm", bg="transparent", limitsize=FALSE)



############ ITERATE MPTP OVER 100 RANDOM TREES FOR 30 FILTERING POINTS ############
############ ITERATE MPTP OVER 100 RANDOM TREES FOR 30 FILTERING POINTS ############

# generate distrib of values
seqs <- generate_thresholds(start=1,max=100000)

# write out fasta for alignment
tax.sub.filtered.fas <- tab2fas(df=tax.sub ,seqcol="nucleotides",namecol="asvHash")
base.name <- paste0("temp-local-only/mptp.fasta")
write.FASTA(tax.sub.filtered.fas,here(base.name))

# run raxml over 100 random starting trees
rand.tr <- raxml_random(file=here(base.name),epsilon=10)

# load up 100 ml trees
trs <- read.tree(paste0(here(base.name),".ali.raxml.rba.raxml.mlTrees"))

# run mptp in parallel over a list of trees
mptp.tab.all.trees <- mcmapply(function(x,y) mptp_parallel(df=tax.sub,base.name=base.name,tr=x,num=y,threshold="single",filt=seqs), x=trs, y=seq_along(trs), USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=8)
#file.remove(list.files(here("temp-local-only"),pattern=".nwk",full.names=TRUE))

# join into dataframe
mptp.combined <- bind_rows(mptp.tab.all.trees)

# save dataframe
mptp.combined %>% write_csv(here("temp-local-only/delim-results.csv"))
mptp.combined %>% write_csv(here("temp/results/delim-results.csv"))

# read back in if needed
mptp.combined <- read_csv(here("temp-local-only/delim-results.csv"),show_col_types=FALSE)

# plot
mptp.combined %>% 
    ggplot(aes(x=filterThreshold,y=nClust)) + 
    geom_boxplot(aes(group=filterThreshold),outlier.shape=NA) +
    geom_jitter(height=0,alpha=0.1) + 
    geom_smooth(method=mgcv::gam,formula=y~s(x,bs="cs",k=4),alpha=0) +
    scale_x_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=6) +
    geom_hline(yintercept=68,lty=2) +
    scale_y_continuous(labels=scales::comma_format(accuracy=1),n.breaks=12) + 
    #scale_y_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=8) + 
    ggthemes::theme_clean()
    #


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