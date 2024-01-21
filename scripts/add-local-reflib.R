#!/usr/bin/env Rscript

# libs
library("here")
library("tidyverse")

# refs
source(here::here("scripts/references-load-local.R"))
source(here::here("scripts/references-clean.R"))
locals <- read_csv(file=here("../assets/local-12s.csv"))

# merge
reflib.local <- reflib.cleaned %>% mutate(dbid=as.character(dbid)) %>% bind_rows(locals) %>% arrange(phylum,class,order,family,genus,sciNameValid,dbid)

# write
reflib.local %>% write_csv(file=gzfile(here("assets/reference-library-master.csv.gz")), na="")
reflib.local %>% write_csv(file=here("../meta-fish-pipe/assets/meta-fish-lib-v259.csv"), na="")
