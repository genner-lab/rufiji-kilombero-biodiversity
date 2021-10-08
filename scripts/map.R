#!/usr/bin/env Rscript

# lib
library("here")
library("tidyverse")
library("ggmap")
library("ggrepel")
library("svglite")

# read in data and filter
events <- read_csv(here("assets/events-master.csv")) %>% 
    mutate(site=str_split_fixed(localitySite,"-",2)[,1]) %>%
    distinct(year,site,decimalLatitude,decimalLongitude) %>%
    arrange(year,site) %>%
    group_by(year,site) %>%
        slice(n=1) %>%
    ungroup() %>%
    mutate(year=as.factor(year))

events %>% print(n=Inf)

# register with google api at
# https://cloud.google.com/maps-platform/
# save key in temp dir
source(here("temp/google-key.R"))
ggmap::register_google(key=google.key)

# get map
map.tz <- ggmap::get_googlemap("Dakawa",maptype="hybrid",zoom=7)
#ggmap(map.tz)

#close up maps
map.tz <- ggmap::get_googlemap("Kilombero Bridge",maptype="hybrid",zoom=16)
map.tz <- ggmap::get_googlemap("Utete",maptype="hybrid",zoom=10)
map.tz <- ggmap::get_googlemap("Mngeta",maptype="hybrid",zoom=10)


# plot
p <- map.tz %>% ggmap() + 
    geom_point(data=events,aes(x=decimalLongitude,y=decimalLatitude,color=year),alpha=1,size=3,shape=20) + #
    geom_label_repel(data=events,aes(x=decimalLongitude,y=decimalLatitude,label=site,fill=year),color="white",label.padding=0.1,min.segment.length=0,size=2,max.overlaps=25,seed=42) +
    xlab("Longitude") + 
    ylab("Latitude") +
    theme_bw()
plot(p)

# save
ggsave(filename=here("temp/maps/map-sites.pdf"),plot=p)


# for fixing seq sheet [delete]
library("tidyverse")
library("here")
events <- read_csv(here("assets/events-master.csv"))
samples <- read_csv(here("temp/sequencing/libs.csv"))

samples %>% left_join(select(events,fieldID,eventID,verbatimLocality)) %>% write_csv(here("temp/sequencing/libs-events.csv"))

# quick plot of seq depth
fbs <- read_csv(here("results/fishes-by-sample.csv"))
pl <- fbs %>% 
    mutate(site=str_split_fixed(eventID,"-",4)[,3]) %>% 
    group_by(library,site,sampleHash) %>% 
    summarise(reads=sum(nReads)) %>% 
    ggplot(aes(x=site,y=reads)) + 
        geom_boxplot() + 
        geom_jitter() +
        scale_y_log10()

plot(pl)

#save 
ggsave(plot=pl,filename=here("../temp/results/Results_08-10-2021/reads.pdf"),width=14,height=6,units="in")
