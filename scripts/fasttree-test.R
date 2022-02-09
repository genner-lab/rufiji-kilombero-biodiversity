

./FastTreeDbl -gtr -gamma -nt mptp.fasta.ali > fasttree.nwk
./FastTreeDbl -nt -gamma mptp.fasta.ali > fasttree.nwk
./FastTreeDbl -nt mptp.fasta.ali > fasttree.nwk

cp /home/rc16041/Software/FastTree/fasttree.nwk /home/rc16041/Projects/genner-lab/rufiji-kilombero-biodiversity/temp-local-only/fasttree.nwk

# load up tree
tr <- read.tree(here("temp-local-only/fasttree.nwk"))
tr <- root_on_longest_edge(tr)

# DROP TIPS FROM TREES
tr.drop <- drop_tips(df=tax.sub,filt=600,tree=tr) 

min(tr.drop$edge.length)

write.tree(tr.drop,here("temp-local-only/fasttree.drop.nwk"))

run_mptp(file=here("temp-local-only/fasttree.drop.nwk"),threshold="single",minbr=5e-20)

dlim <- parse_mptp(file=here("temp-local-only/fasttree.drop.nwk.mptp.out.txt"),filt=600,num=1,species=FALSE)

tax.sub.dlim <- tax.sub %>% left_join(dlim,by="asvHash") %>% filter(!is.na(mptpDelim))

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
set.seed(42)
cols <- sample(getPalette(n=length(unique(pull(tax.sub.dlim,mptpDelim)))))


p <- ggtree(tr.drop, ladderize=TRUE, color="grey20") + xlim(0,2)
p <- p %<+% tax.sub.dlim
p <- p + geom_tiplab(aes(fill=mptpDelim,label=tipLabel),geom="label",label.size=0,label.padding=unit(0.2,"lines"),offset=0.01) + 
    theme(legend.position="none") +
    scale_fill_manual(values=cols)


ggsave(plot=p, filename=here("temp-local-only/tree.pdf"), width=length(tr.drop$tip.label)/4, height=length(tr.drop$tip.label)/2, units="cm", bg="transparent", limitsize=FALSE)





# FUN TO GENERATE BOOTSTRAPPED ALIGNMENTS
get_bs <- function(file,reps) {
    string.mafft <- paste0("mafft --thread -1 --maxiterate 1000 --localpair ",file," > ",file,".ali")
    system(command=string.mafft,ignore.stdout=FALSE)
    string.bs <- paste0("raxml-ng --bsmsa --msa ",file,".ali --model TN93+G --seed 42 --redo --threads auto --bs-trees ",reps)
    system(command=string.bs,ignore.stdout=FALSE)
}


export PATH=~/Software/FastTree:$PATH
source(here::here("scripts/libs-funs.R"))

run_fasttree <- function(file) {
    string.ft <- paste0("FastTreeDbl -nt ",file," > ",file,".nwk")
    system(command=string.ft,ignore.stdout=FALSE)
    ft.tr <- ape::read.tree(file=paste0(file,".nwk"))
    return(ft.tr)
}

tr <- run_fasttree(file=here("temp-local-only/testing/mptp.fasta.ali.raxml.bootstrapMSA.100.phy"))

#system(command="export PATH=~/Software/FastTree:$PATH")


############ CREATE 100 BOOTSTRAP TREES ############
############ CREATE 100 BOOTSTRAP TREES ############

# write out fasta for alignment
tax.sub.filtered.fas <- tab2fas(df=tax.sub,seqcol="nucleotides",namecol="asvHash")
write.FASTA(tax.sub.filtered.fas,base.name)

# set reps (100 bs trees)
nreps <- 100

# generate alignment and generate n bootstraps
get_bs(file=base.name,reps=nreps)


############ RUN FASTTREE ON EACH BS ############
############ RUN FASTTREE ON EACH BS ############

# get paths for all bs alignments
bs.paths <- paste0(base.name,".ali.raxml.bootstrapMSA.",1:nreps,".phy")
bs.paths <- paste0(base.name,".ali.boot.",str_pad(as.character(1:100),pad="0",width=4),".fas")


# run fasttree in parallel
mcmapply(function(x) run_fasttree(file=x), x=bs.paths, USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=8)

trs <- lapply(paste0(bs.paths,".nwk"),read.tree)
class(trs) <- "multiPhylo"




# run mptp in parallel over a list of trees
mptp.tab.all.trees <- mcmapply(function(x,y) mptp_parallel(df=tax.sub,base.name=base.name,tr=x,num=y,threshold="single",filt=seqs,parsespp=TRUE), x=trs, y=seq_along(trs), USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=8)
#mptp_parallel(df=tax.sub,base.name=base.name,tr=trs[[2]],num=1,threshold="single",filt=seqs,parsespp=TRUE)

file.remove(list.files(here("temp-local-only"),pattern=".nwk",full.names=TRUE))
file.remove(list.files(here("temp-local-only"),pattern="boot",full.names=TRUE))


## try alt bootstrap
tax.sub.filtered.fas.mat <- as.matrix(read.FASTA(here("temp-local-only/mptp.fasta.ali")))


# DROP TIPS FROM TREES
boot_time <- function(ali,seed) {
    set.seed(seed)
    ali.boot <- ali[,sample(ncol(ali),replace=TRUE)]
    return(ali.boot)
}

set.seed(42)
r.seeds <- sample(9:999999,replace=FALSE,size=100)

boot.mats <- mcmapply(function(x) boot_time(ali=tax.sub.filtered.fas.mat,seed=x), x=r.seeds, USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=8)

boot.files <- paste0(base.name,".ali.boot.",str_pad(as.character(1:100),pad="0",width=4),".fas")

mcmapply(function(x,y) write.FASTA(x,file=y), x=boot.mats, y=boot.files, USE.NAMES=FALSE, SIMPLIFY=FALSE, mc.cores=1)


### run single
ft.tree <- run_fasttree(file=here("temp-local-only/mptp.fasta.ali"))


mptp.tab.all.trees <- mptp_parallel(df=tax.sub,base.name=base.name,tr=ft.tree,num=1,threshold="single",filt=seqs,parsespp=TRUE)


mptp.tab.all.trees %>% 
    ggplot(aes(x=filterThreshold,y=nClust)) + 
    geom_point() +
    geom_smooth(method=mgcv::gam,formula=y~s(x,bs="cs",k=4),alpha=0,color="seagreen") +
    scale_x_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=5) +
    geom_hline(yintercept=68,lty=2,color="cornflowerblue") +
    geom_hline(yintercept=max.d,lty=2,color="firebrick1") +
    #geom_hline(yintercept=median(pull(mptp.combined,nClust)),lty=2,color="forestgreen") +
    #scale_y_continuous(labels=scales::comma_format(accuracy=1),n.breaks=12) + 
    #scale_y_continuous(trans="log10",labels=scales::comma_format(accuracy=1),n.breaks=12) + 
    ggthemes::theme_clean()
    #


# get density
adj <- 1
d <- density(pull(mptp.tab.all.trees,nClust),adjust=adj)
max.d <- ceiling(d$x[which.max(d$y)])
print(max.d)

# plot density
mptp.tab.all.trees %>% 
    ggplot(aes(x=nClust)) + 
    geom_density(adjust=1) + 
    geom_vline(xintercept=max.d,color="firebrick1",linetype=2) +
    scale_x_continuous(labels=scales::comma_format(accuracy=1),n.breaks=12) +
    ggthemes::theme_clean()

est.threshold <- mptp.tab.all.trees %>% filter(nClust == max.d) %>% pull(filterThreshold) %>% median()
print(est.threshold)


# BEAST
# beast run1 (run from temp-local-only)
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
beast -beagle -overwrite -seed 802221 mptp.fasta.ali.run1.xml
