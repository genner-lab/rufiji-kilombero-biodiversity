#!/usr/bin/env Rscript

# lib
library("here")
library("tidyverse")
library("ggmap")
library("ggrepel")
library("svglite")

# read in data and filter
events <- read_csv(here("assets/events-master.csv")) %>% 
    distinct(year,fieldID,decimalLatitude,decimalLongitude) %>%
    arrange(year,fieldID) %>%
    group_by(year,fieldID) %>%
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
map.tz <- ggmap::get_googlemap("Dakawa",maptype="hybrid",zoom=8)
#ggmap(map.tz)

#close up maps
map.tz <- ggmap::get_googlemap("Kilombero Bridge",maptype="hybrid",zoom=16)
map.tz <- ggmap::get_googlemap("Utete",maptype="hybrid",zoom=10)
map.tz <- ggmap::get_googlemap("Mngeta",maptype="hybrid",zoom=10)


# plot
p <- map.tz %>% ggmap() + 
    geom_point(data=events,aes(x=decimalLongitude,y=decimalLatitude,color=year),alpha=1,size=3,shape=20) + #
    geom_label_repel(data=events,aes(x=decimalLongitude,y=decimalLatitude,label=fieldID,fill=year),color="white",label.padding=0.1,min.segment.length=0,size=2,max.overlaps=25,seed=42) +
    xlab("Longitude") + 
    ylab("Latitude") +
    theme_bw()
plot(p)

# save
ggsave(filename=here("temp/maps/map-mngeta-zoom10.svg"),plot=p)


# for fixing seq sheet [delete]
library("tidyverse")
library("here")
events <- read_csv(here("assets/events-master.csv"))
samples <- read_csv(here("temp/sequencing/libs.csv"))

samples %>% left_join(select(events,fieldID,eventID,verbatimLocality)) %>% write_csv(here("temp/sequencing/libs-events.csv"))
