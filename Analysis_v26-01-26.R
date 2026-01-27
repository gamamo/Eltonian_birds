# 1) Load packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(terra)
library(patchwork)
library(tidyterra)
library(rnaturalearth)
library(ggrepel)
library(scales)
library(rgbif)
library(vegan)
library(iNEXT)
library(flextable)

# 1.1) Load the shape file of biogeographic realms ---------------------------------------------------------------
# Load and do some fixing
# This shapefile comes from Holt et al. 2013 Science - I cannot share it as it does not belongs to me

realms <- vect("newRealms.shp") 
realms <- realms[,"Realm"]
realms[realms$Realm == "Panamanian", "Realm"] <- "Neotropical"
realms[realms$Realm == "Oceanina", "Realm"] <- "Oceanian"

# 1.2) Load and wrangle datasets  ---------------------------------------------------------
## 1.2.1) AVONET ---------------------------------------------------------------------------
# Avonet can be found at Tobias et al. 2022 Ecology Letters
avonet <- read_csv("AVONET1_BirdLife.csv") 

# select species only
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1) |>
  dplyr::select(Species,Trophic.Niche)

# Merge with biogeographic realms
avonetSpShp <- avonet |> 
  dplyr::mutate(Species = Species1) |>
  dplyr::select(Species,Centroid.Latitude,Centroid.Longitude)

# convert to a shapefile
avonetSpShp <- vect(avonetSpShp, geom=c("Centroid.Longitude","Centroid.Latitude"),crs="EPSG:4326")
avonetSpShp_coords <- crds(avonetSpShp)

# extract realms
avonetSpRealm <- extract(realms,avonetSpShp)
avonetSpRealm <- cbind(avonetSpRealm, avonetSpShp_coords)

# assign species names
avonetsSpRealm <- cbind(avonetSp ,avonetSpRealm ) |> 
  dplyr::select(Species, Realm,x,y, Trophic.Niche) |> 
  dplyr::mutate(Realm = case_when(Realm == 'Panamanian' ~ 'Neotropical',
                                  TRUE ~ Realm ))

#snap the NA features to the nearest realm
avonetsSpRealmNA <- avonetsSpRealm |> 
  filter(is.na(Realm))
avonetsSpRealmNA_vect <- vect(avonetsSpRealmNA,geom=c("x","y"),crs="EPSG:4326")

# use the nearest algorithm
ac <- nearest(avonetsSpRealmNA_vect,realms)
ac <- as.data.frame(ac)
r <- as.data.frame(realms)
r$realmID <- as.numeric(rownames(r))

l <- list()
for(i in r$realmID){
  temp <- r |> filter(realmID==i) |> select(Realm) |> pull()
  l[[i]] <- ac |> filter(to_id==i) |> mutate(Realm = temp)
}
ac2 <- map_df(l, ~.x) |> select(from_x,from_y, Realm) |> 
  rename(x=from_x, y=from_y)

avonetsSpRealm2 <- left_join(avonetsSpRealm,ac2, by=c("x","y")) |> 
  mutate(Realm = coalesce(Realm.x, Realm.y)) |> 
  select(-Realm.x, -Realm.y)

avonetsSpRealm3 <- avonetsSpRealm2 |> select(Species, Trophic.Niche, Realm)

## 1.2.2) BIRDBASE -------------------------------------------------
# This data file comes from Şekercioğlu et al. 2025 Scientific Reports
birdbase <- read_csv("BIRDBASE_afterGMM.csv")
birdbaseSp <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  unite("Species",Genus:Species, sep =" ") 

birdbaseSpNO <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  unite("Species",Genus:Species, sep =" ") |> 
  filter(`Primary Diet`=="Fruit" | `Primary Diet`=="Seed" ) |> 
  dplyr::select(Species)

# 1.3) Load species interactions -----------------------------------------------------
nets <- read_csv("datalong_v08_06_2023_CLEANED.csv",locale = locale(encoding = "UTF-8")) |> 
  select(Scientific, plant_id, database, ref, interaction, lat, lon)

# import the new networks 2023-2025
new <- read_csv(here("new networks","New_Networks_210126_utf8.csv"),locale = locale(encoding = "UTF-8"))
new <- new |>  mutate(lon = str_trim(lon, side = "right")) |> 
  mutate(lon = as.numeric(lon))

# check if all species are AVES (Takes a while) - this code came from chatGPT
if(F){
res <- new |>
  distinct(Scientific) |>
  mutate(
    class = map_chr(
      Scientific,
      ~ {
        x <- name_backbone(name = .x)$class
        if (is.null(x)) NA_character_ else x
      }
    )
  )

non_aves <- res |>
  filter(is.na(class) | class != "Aves")
}

# rbind with the new networks 2023-2025
nets <- rbind(nets,new)
nets <- nets |> 
  rename(Species = Scientific)

# do small fixing
nets[nets$Species=="ï»¿Mionectes oleagineus" ,"Species"] <- "Mionectes oleagineus"
nets <- nets |>  
  filter(Species!="Birds") |> 
  filter(Species!="Parrots") |> 
  filter(Species!="Contopus spp") |> 
  filter(Species!="Pogoniulus sp")

# get the realm for each species
nets2 <- left_join(nets, avonetsSpRealm3)

#check if any species has not been assigned to a realm
spnas <- nets2 |> filter(is.na(Realm)) |> select(Species) |> distinct() 

#revise the taxonomy for the species without realm
taxo <- read_csv(here("AVONET", 'AVONET_taxonomy.csv'))

l <- list()
for (i in spnas$Species ){
  temp <- spnas |> filter(Species == i) |> pull()
 
  if(temp %in% taxo$Species2_eBird == T){
    a1 <- unique(taxo[taxo$Species2_eBird %in% temp,"Species1_BirdLife"])
    a1 <- a1 |> distinct() |> pull()
  } else {
    a1 <- NA}
  if(temp %in% taxo$eBird.species.group ==T){
    b1 <- unique(taxo[taxo$eBird.species.group %in% temp,"Species1_BirdLife"])
    b1 <- b1 |> distinct() |> pull()
  } else {
    b1 <- NA}
  if(temp %in% taxo$Species3_BirdTree  == T){
    c1 <- unique(taxo[taxo$Species3_BirdTree %in% temp,"Species1_BirdLife"])
    c1 <- c1 |> distinct() |>pull()
  }else {
    c1 <- NA}
    l[[i]] <- data.frame(Species = temp, ebird1 = a1, ebird2 = b1, birdtree = c1 )
}
taxo2 <- map_df(l, ~.x) |> distinct(Species,.keep_all = TRUE)

# fix taxonomy by comparing the other taxonomic sources - use birdtree-birdlife relationship first
temp <- nets2
for (j in taxo2$Species){
  print(j)
  n <- taxo2 |> filter(Species==j) |> select(birdtree) |> pull()
  if(is.na(n)){next}
  temp[temp$Species==j,"Species"] <- taxo2[taxo2$Species==j,"birdtree"]
}
nets3 <- temp

# assign realms again
avonetsSpRealm4 <-  avonetsSpRealm3 |> rename(Realm2 = Realm, tf=Trophic.Niche)
nets4 <- left_join(nets3, avonetsSpRealm4 )

nets4 <- nets4 |> 
group_by(Species) %>%
  fill(Trophic.Niche, Realm, .direction = "downup") %>%
  ungroup()

nets5 <- nets4 |> mutate(Realm = coalesce(Realm, Realm2))  |> 
  mutate(Trophic.Niche = coalesce(tf, Trophic.Niche)) |> 
  select(-tf, -Realm2)

spnas2 <- nets5 |> filter(is.na(Realm)) |> select(Species) |> distinct() 

# As there are still taxonomy to fix, compare now with ebird-birdlife relationship 
temp <- nets5
for (j in taxo2$Species){
  print(j)
  n <- taxo2 |> filter(Species==j) |> select(ebird1) |> pull()
  if(is.na(n)){next}
  temp[temp$Species==j,"Species"] <- taxo2[taxo2$Species==j,"ebird1"]
}
nets6 <- temp

# assign realms again
nets6 <- left_join(nets6, avonetsSpRealm4 )

nets6 <- nets6|> 
  group_by(Species) %>%
  fill(Trophic.Niche, Realm, .direction = "downup") %>%
  ungroup()

nets6 <- nets6 |> mutate(Realm = coalesce(Realm, Realm2))  |> 
  mutate(Trophic.Niche = coalesce(tf, Trophic.Niche)) |> 
  select(-tf, -Realm2)

spnas3 <- nets6 |> filter(is.na(Realm)) |> select(Species) |> distinct() 

# network file to use in the analysis with the necessary filters ------------------------------------
netsAna <- nets6 |> 
  filter(interaction > 0) 

#### 1.3) Make a map --------------------------------------------------------------
wrld <- ne_countries(continent = c("africa","south america","asia","oceania","europe",
                                   "north america"),
                     type = "countries")

nets_map <- netsAna |> 
  select(database,lat,lon ) |> 
  distinct()

nets_map <- vect(nets_map, geom=c("lon", "lat"), crs="EPSG:4326")
writeVector(nets_map, "networks_210126.shp",overwrite=T)

m1 <- ggplot() +
  geom_spatvector(data=wrld, fill=NA)+
  geom_spatvector(data=realms, aes(fill=Realm))+
  geom_spatvector(data=nets_map,size=2,color="orange",shape=21,fill="black")+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#plot.background = element_rect(fill="white",color="black"),
        legend.position = "none")
m1
ggsave(here("products","map.jpeg"), units = "cm", height = 10, width = 20)

# Cartogram
nets_cartog <- netsAna |> 
  select(database,lat,lon, interaction ) |> 
  group_by(database, lat, lon) |> 
  summarise(Interactions = sum(interaction))
nets_cartog <- vect(nets_cartog , geom=c("lon", "lat"), crs="EPSG:4326")

m2 <- ggplot() +
  geom_spatvector(data=wrld, fill=NA)+
  geom_spatvector(data=realms, aes(fill=Realm))+
  geom_spatvector(data=nets_cartog,aes(size=Interactions),
                  color="black",shape=21,fill="orange",alpha=0.8)+
  scale_fill_viridis_d()+
  scale_size_continuous(range = c(1, 10))+
  theme_minimal()+
  theme(#plot.background = element_rect(fill="white",color="black"),
        legend.position = "bottom",
        legend.box = "vertical")
m2
ggsave(here("products","map_sizes.jpeg"), units = "cm", height = 10, width = 20)

(m1/m2) +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","map_int.jpeg"), units = "cm", height = 20, width = 20)

# 2) Calculate some stats ----------------------------------------------------
## 2.0) How many networks -----------------------------------------------------
netsAna |> 
  select(database) |> 
  n_distinct()

## 2.0.1) How many bird species
netsAna |> 
  select(Species) |> 
  n_distinct()

153300/11009
153300/11589

## 2.1) How many networks per realm ---------------------------------------------
netsAna |> 
  select(Realm, database) |> 
  drop_na() |> 
  distinct() |> 
  group_by(,Realm) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup() -> Nnets

Nnets

NnetsP <- Nnets |>   
  mutate(perc=(perc = `n` / sum(`n`))) |> 
  mutate(labels = scales::percent(perc))

i0 <- NnetsP |> 
  ggplot(aes(x="", y=perc,fill = Realm)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label_repel(aes(label = labels),color="white",
             position = position_stack(vjust = 0.4),
             show.legend = FALSE,size=5)+
  scale_fill_viridis_d()+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))+
  labs(title="The proportion of seed-dispersal networks \nper biogeographical realm")
i0

#### 2.2) Number of interactions per realm -------------------------------------------
netsAna |> 
  group_by(Realm) |> 
  count() |> 
  drop_na() |> 
  distinct() |> 
  ungroup() |> 
  mutate(perc=(perc = `n` / sum(`n`))) |> 
  mutate(labels = scales::percent(perc)) |> 
  rename(int = n) -> Nints

Nints

sum(Nints$int)

# plot the number of interactions per realm
i1 <- netsAna|> 
  group_by(Realm) |> 
  drop_na() |> 
  count() |> 
  ungroup() |> 
  mutate(perc=(`n` / sum(`n`))) |> 
  mutate(labels = scales::percent(perc)) |> 

  ggplot(aes(x="", y=perc,fill = Realm)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label_repel(aes(label = labels),color="white",
             position = position_stack(vjust = 0.2),
             show.legend = FALSE,size=5)+
  scale_fill_viridis_d()+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))+
  labs(title="The number of interactions \nper biogeographical realm")
i1

# plot the normalized number of interactions per realm by number of networks
inter <- cbind(Nints,Nnets$n)
inter <- inter |> rename(nets = 'Nnets$n')

i2 <- inter|> 
  select(Realm, int,nets) |> 
  mutate(norm = int/nets) |> 
  mutate(perc=(`norm` / sum(`norm`))) |> 
  mutate(labels = scales::percent(perc,accuracy = 0.01)) |> 
  
  ggplot(aes(x="", y=perc,fill = Realm)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label_repel(aes(label = labels),color ="white",
             position = position_stack(vjust = 0.3),
             show.legend = FALSE,size=5)+
  scale_fill_viridis_d()+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))+
  labs(title="The normalized number of interactions \nper biogeographical realm")
i2

### Make a composite graph --------------------------------------------------
i0+i1+i2+
  plot_annotation(tag_levels = "a",tag_suffix = ")")+
  plot_layout(guides = 'collect')&
  theme(legend.position = "bottom",legend.text = element_text(size=12))
ggsave(here("products","fig1.jpeg"),units="cm", height = 15, width = 30)

## 2.3) Species more and less common --------------------------------------------
common <- netsAna |> 
  group_by(Species, Realm) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup()

#per realm
re <- common |> dplyr::select(Realm) |> distinct() |> pull()

l <- list()
for (i in re){
  l[[i]] <- common |> filter(Realm==i)
}
intSP <- map_df(l, ~.x)

# make a graphic with the top ten species per realm
resu<- list()
for (i in re){
  resu[[i]] <- common |> filter(Realm==i) |> 
    slice(1:10)
}
resu <- map_df(resu, ~.x)
intSP <- resu

resu |> 
  ggplot(aes(x=reorder(Species,n),y=n)) +
  geom_col() +
  facet_wrap(.~Realm,drop=T,scales="free",ncol=2)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=18),
        axis.text.y = element_text(face = "italic"),
        legend.position = "none",
        axis.ticks.y = element_blank())+
  labs(y="Number of interactions of each species",x="")+
  coord_flip()
ggsave(here("products","Fig2.jpeg"),units="cm",height = 30,width = 38)

# standardize the calculations per population size
callaghan <- read_csv(here("Abundances","all_species_summary_table.csv")) |> 
  select('Scientific name','Abundance estimate') |> 
  rename(Species = 'Scientific name', Abundances='Abundance estimate')

abu <- left_join(intSP, callaghan)
abu |> ggplot(aes(x=n,y=Abundances))+
  geom_point()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  geom_smooth(method = "lm")+
  facet_wrap(.~Realm,scales="free",ncol=2)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=14),
        legend.position = "none")+
  labs(x="Number of species interactions in seed-dispersal networks",
       y="Population abundances")
ggsave(here("products","SI_1.jpeg"),units = "cm", height = 30, width = 35)

# 3)Taxonomic gaps ------------------------------------------------------------
## 3.1) Avonet -----------------------------------------------------------
# load nets data
netsAna |> select(Species)  |> n_distinct()

netsSp <- netsAna |> 
  select(Species) |> 
  distinct() |> 
  mutate(Nets = 1)

##3.2) Get proportions ----------------------------------------------------------
### 3.2.1) Overall -----------------------------------------------------------------------
# no filter
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1) |>
  dplyr::select(Species,Trophic.Niche)

avo_nets <- left_join(avonetSp,netsSp )

# make a figure
g1 <- avo_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Nets) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc)) |> 
  mutate(Nets = factor(Nets, levels = c(0, 1),
                       labels = c("Absent in networks", "Present in networks"))) |> 
  ggplot(aes(x="", y=perc,fill = Nets)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label(aes(label = labels), color = c("white", "white"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,size=5)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))+
  scale_fill_manual(
    values = c(
      "Absent in networks" = "lightsteelblue4",  
      "Present in networks" = "darkblue"))+
  labs(title="Representation in seed-dispersal networks \ncompared to AVONET")
g1

#### By realm --------------------------------------------------------------
avonetsSpRealmNF <- avonetsSpRealm3

avo_nets_realmNF <- left_join(avonetsSpRealmNF ,netsSp ) 

avo_nets_realmNF |> 
  filter(!is.na(Realm)) |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Realm,Nets) |> 
  count() |> 
  dplyr::group_by(Realm) |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  filter(Nets==1)

g1rNF <- avo_nets_realmNF |> 
  filter(!is.na(Realm)) |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Realm,Nets) |> 
  count() |> 
  dplyr::group_by(Realm) |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2)) |> 
  filter(Nets==1) |>  
  # mutate(labels = scales::percent(proportion)) |> 
  ggplot(aes(x=reorder(Realm,-n), 
             y=proportion,
             fill = Nets,
             label=proportion)) +
  geom_col(stat="identity") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=12),
        legend.position = "none",
        axis.ticks = element_blank())+
  coord_flip()+
  geom_text(nudge_y= +3.5,
            color="black",
            size = 3,
            fontface="bold")+
  scale_y_continuous(limits=c(0,100))+
  labs(#title="Representation in seed-dispersal networks",
    x="Biogeograpihc realms",y="Proportion of species in seed-dispersal networks \nin relation to regional species pools")
g1rNF

### 3.2.2) Filter 1 ---------------------------------------------------------------------
#apply filter1
avonetSp2 <- avonet |> 
  dplyr::mutate(Species = Species1) |>
  dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche) 

avo_netsF1 <- left_join(avonetSp2,netsSp )

# make a figure
g1F1 <- avo_netsF1 |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                  TRUE ~ Nets )) |> 
  dplyr::group_by(Nets) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc)) |> 
  mutate(Nets = factor(Nets, levels = c(0, 1),
                       labels = c("Absent in networks", "Present in networks"))) |> 
ggplot(aes(x="", y=perc,fill = Nets)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label(aes(label = labels), color = c("white", "white"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,size=5)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))+
  scale_fill_manual(
    values = c(
      "Absent in networks" = "lightsteelblue4",  
      "Present in networks" = "darkblue"))+
  labs(title="Representation in seed-dispersal networks \ncompared to AVONET")
g1F1

##### (How many species out ?) -----------------------------------------------------
#how many species in
sum(netsSp$Species %in% avonetSp2$Species)

#how many species out and what is their trophic niche
outAvonet <- anti_join(netsSp,avonetSp2)
outAvonet |> select(Species) |> n_distinct()

# how many out considering the number of species in networks
73499/1523

outAvonet <- left_join(outAvonet,avonetSp)
outAvonet |> 
  group_by(Trophic.Niche) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc))

#### By realm---------------------------------------------------------------------
avonetsSpRealm4 <- avonetsSpRealm3  |> 
  dplyr::filter(Trophic.Niche == "Omnivore"| Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") 

avo_nets_realmF1 <- left_join(avonetsSpRealm4 ,netsSp ) 

avo_nets_realmF1 |> 
  filter(!is.na(Realm)) |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Realm,Nets) |> 
  count() |> 
  dplyr::group_by(Realm) |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  filter(Nets==1)

g1rF1 <- avo_nets_realmF1 |> 
  filter(!is.na(Realm)) |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Realm,Nets) |> 
  count() |> 
  dplyr::group_by(Realm) |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2)) |> 
  filter(Nets==1) |>  
 # mutate(labels = scales::percent(proportion)) |> 
  ggplot(aes(x=reorder(Realm,-n), 
             y=proportion,
             fill = Nets,
             label=proportion)) +
  geom_col(stat="identity") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=12),
        legend.position = "none",
        axis.ticks = element_blank())+
  coord_flip()+
  geom_text(nudge_y= +3.5,
            color="black",
            size = 3,
            fontface="bold")+
  scale_y_continuous(limits=c(0,100))+
  labs(#title="Representation in seed-dispersal networks",
       x="Biogeograpihc realms",y="Proportion of species in seed-dispersal networks \nin relation to regional species pools")
g1rF1


### 3.2.3) Filter 2 -------------------------------------------------------------------
#apply filter 2
avonetSpNO <- avonet |> 
  dplyr::mutate(Species = Species1) |>
  dplyr::filter(Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche)

avo_netsNO <- left_join(avonetSpNO,netsSp )

sum(netsSp$Species %in% avonetSpNO$Species)

g1NO <- avo_netsNO |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Nets) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc)) |> 
  mutate(Nets = factor(Nets, levels = c(0, 1),
                       labels = c("Absent in networks", "Present in networks"))) |> 
  ggplot(aes(x="", y=perc,fill = Nets)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label(aes(label = labels), color = c("white", "white"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,size=5)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))+
  scale_fill_manual(
    values = c(
      "Absent in networks" = "lightsteelblue4",  
      "Present in networks" = "darkblue"))+
  labs(title="Representation in seed-dispersal networks \ncompared to AVONET")
g1NO

# by realm
avonetsSpRealm4NO <- avonetsSpRealm3  |> 
  dplyr::filter( Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") 

avo_nets_realmNO <- left_join(avonetsSpRealm4NO ,netsSp ) 

avo_nets_realmNO |> 
  filter(!is.na(Realm)) |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Realm,Nets) |> 
  count() |> 
  dplyr::group_by(Realm) |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  filter(Nets==1)

g1rNO <- avo_nets_realmNO |> 
  filter(!is.na(Realm)) |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Realm,Nets) |> 
  count() |> 
  dplyr::group_by(Realm) |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2)) |> 
  filter(Nets==1) |>  
  # mutate(labels = scales::percent(proportion)) |> 
  ggplot(aes(x=reorder(Realm,-n), 
             y=proportion,
             fill = Nets,
             label=proportion)) +
  geom_col(stat="identity") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=12),
        legend.position = "none",
        axis.ticks = element_blank())+
  coord_flip()+
  geom_text(nudge_y= +3.5,
            color="black",
            size = 3,
            fontface="bold")+
  scale_y_continuous(limits=c(0,100))+
  labs(#title="Representation in seed-dispersal networks",
    x="Biogeograpihc realms",y="Proportion of species in seed-dispersal networks \nin relation to regional species pools")
g1rNO
ggsave(here("products","fig5_noOmnivores.jpeg"),units = "cm", height = 10, width = 15)

##3.3) Birdbase ---------------------------------------------------------
netsSp |> n_distinct()
birdbaseSp |>  n_distinct()

birdbaseSp <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  unite("Species",Genus:Species, sep =" ") 

birdbaseSp2 <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  filter(`Primary Diet`=="Fruit" | `Primary Diet`=="Seed" | `Primary Diet`=="Omnivore" ) |> 
  unite("Species",Genus:Species, sep =" ") 

birdbaseSpNO <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  unite("Species",Genus:Species, sep =" ") |> 
  filter(`Primary Diet`=="Fruit" | `Primary Diet`=="Seed" ) |> 
  dplyr::select(Species)

##3.4) Get proportions ----------------------------------------------------------
### 3.4.1) Overall ---------------------------------------------------------------------
bb_nets <- left_join(birdbaseSp,netsSp )

g2 <- bb_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Nets) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc)) |> 
  mutate(Nets = factor(Nets, levels = c(0, 1),
                       labels = c("Absent in networks", "Present in networks"))) |> 
  ggplot(aes(x="", y=perc,fill = Nets)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label(aes(label = labels), color = c("white", "white"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,size=5)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))+
  scale_fill_manual(
    values = c(
      "Absent in networks" = "lightsteelblue4",  
      "Present in networks" = "darkblue"))+
  labs(title="Representation in seed-dispersal networks \ncompared to BirdBase")
g2

### 3.4.2) Filter 1 ---------------------------------------------------------------------
#apply filter1
bb_netsF1 <- left_join(birdbaseSp2,netsSp )

# how many species in
sum(netsSp$Species %in% birdbaseSp2$Species)

#how many species out and what is their trophic niche
outbirdbase <- anti_join(netsSp,birdbaseSp2)
outbirdbase |> select(Species) |> n_distinct()

#how many species out and what is their trophic niche
outbirdbase <- left_join(outbirdbase,birdbaseSp )
outbirdbase |> 
  group_by(`Primary Diet`) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc))


g2F1 <- bb_netsF1 |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Nets) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc)) |> 
  mutate(Nets = factor(Nets, levels = c(0, 1),
                       labels = c("Absent in networks", "Present in networks"))) |> 
  ggplot(aes(x="", y=perc,fill = Nets)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label(aes(label = labels), color = c("white", "white"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,size=5)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))+
  scale_fill_manual(
    values = c(
      "Absent in networks" = "lightsteelblue4",  
      "Present in networks" = "darkblue"))+
  labs(title="Representation in seed-dispersal networks \ncompared to BirdBase")
g2F1

### 3.4.3) Filter 2 ---------------------------------------------------------------------
# apply filter2
bb_netsNO <- left_join(birdbaseSpNO,netsSp )

g2NO <- bb_netsNO |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  dplyr::group_by(Nets) |> 
  count() |> 
  ungroup()  |>  
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc)) |> 
  mutate(Nets = factor(Nets, levels = c(0, 1),
                       labels = c("Absent in networks", "Present in networks"))) |>
  ggplot(aes(x="", y=perc,fill = Nets)) +
  geom_col(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label(aes(label = labels), color = c("white", "white"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,size=5)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))+
  scale_fill_manual(
    values = c(
      "Absent in networks" = "lightsteelblue4",  
      "Present in networks" = "darkblue"))+
  labs(title="Representation in seed-dispersal networks \ncompared to BirdBase")
g2NO

## Make a composite graph -------------------------------------------------

((g1/g2)|g1rNF)  +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products", "fig4_new.jpeg"),units="cm", width=25, height = 15)


((g1F1/g2F1)|g1rF1)  +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products", "fig4_Omnivore.jpeg"),units="cm", width=25, height = 15)

((g1NO/g2NO)|g1rNO)  +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products", "fig4_NoOmnivore.jpeg"),units="cm", width=25, height = 15)


# 4) Accumulation curves -----------------------------------------------------
com <- netsAna |> 
  select(Species, Realm, database) |> 
  distinct() |> 
  drop_na() |> 
  mutate(presence = 1) |> 
  pivot_wider(names_from = Species, values_from = presence, values_fill = 0) |> 
  select(-database)

l <- list()
for (i in unique(com$Realm)) {
  print(i)
  temp <- com |> filter(Realm==i) |> 
    select(-Realm)
  sa <- specaccum(temp, "random") 
  sa_df <- tibble(
    sites = sa$sites,
    richness = sa$richness,
    sd = sa$sd,
    Realm = i
  )
  l[[i]] <- sa_df
}
acu <- map_df(l, ~.x)

#get the end point of each realm
lab_df <- acu %>%
  group_by(Realm) %>%
  slice_max(sites, n = 1)
lab_df[7,'richness'] <- 220
lab_df[7,'sites'] <- 170


#  plot accumulation curves
ggplot(acu, aes(sites, richness, color = Realm)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = richness - sd, ymax = richness + sd, fill = Realm),
    alpha = 0.15,
    color = NA,
    show.legend = F) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=14),
        legend.position = "none")+
  scale_color_viridis_d()+
  scale_x_continuous(expand = c(0,1))+
  labs(x = "Number of seed dispersal networks", y = "Accumulated species richness")+
  geom_text(
    data = lab_df, aes(label = Realm),hjust = 0,vjust= 1.5, size = 4,
    show.legend = FALSE)+
  coord_cartesian(clip = "off")
ggsave(here("products","fig_acum.jpeg"), units="cm", height = 15, width = 25)


# calculate estimated richness using iNext
bird <- netsAna |> 
  select(Species, Realm) |> 
  #filter(Realm=="Neotropical") |> 
  #distinct() |> 
  drop_na()

lst <- bird |> 
  count(Species,Realm) |> 
  group_by( Realm) |> 
  summarise(values = list(n), .groups = "drop") |> 
  tibble::deframe()

estimation <- iNEXT(lst, q=0, datatype="abundance") 

est_rich <- estimation$AsyEst
est_rich |> 
  filter(Diversity == "Species richness") |> 
  mutate(delta = Estimator - Observed) |> 
  flextable() %>%
  autofit() |> 
  save_as_docx(path=here("products","Estimated_richness.docx"))
  



