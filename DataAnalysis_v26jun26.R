# Assessing the Eltonian shortfall for bird species in seed dispersal networks
# Moulatlet et al.--------------------------------------------------

# 1) Load packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(terra)
library(patchwork)
library(tidyterra)
library(rnaturalearth)
library(ggrepel)
library(iNEXT)
library(flextable)


## 1.2) Load the dataset --------------------------------------------------------

netsAna <- read_csv("DataForAnalysis_v11jun26.csv")

## 1.3) Make a map of the distribution of networks and interactions ---------------------

realms <- vect("newRealms.shp") # Holt et al.
realms <- realms[,"Realm"]
realms[realms$Realm == "Panamanian", "Realm"] <- "Neotropical"
realms[realms$Realm == "Oceanina", "Realm"] <- "Oceanian"

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

# make a Cartogram with the size of the dots representing the number of interactions

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

# This will be figure 1:
(m1/m2) +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","map_int.jpeg"), units = "cm", height = 20, width = 20)

# 2) Calculate some statistics ----------------------------------------------------
## 2.0) How many networks -----------------------------------------------------
netsAna |> 
  select(database) |> 
  n_distinct()

## 2.0.1) How many bird species in total
netsAna |> 
  select(Species) |> 
  n_distinct()

153300/11009 # divided by the total species in avonet
153300/11589 # divided by the total species in birdbase

## 2.1) How many networks per biogeographic realm ---------------------------------------------
netsAna |> 
  select(Realm, database) |> 
  drop_na() |> 
  distinct() |> 
  group_by(,Realm) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup() -> Nnets

print(Nnets)

## 2.1.1) What is the proportion of networks per realm

NnetsP <- Nnets |>   
  mutate(perc=(perc = `n` / sum(`n`))) |> 
  mutate(labels = scales::percent(perc))

i0 <- NnetsP |> 
  ggplot(aes(x="", y=perc,fill = Realm)) +
  geom_bar(stat="identity",color="white") +
  coord_polar("y", start=0)+
  geom_label_repel(aes(label = labels),
                   position = position_stack(vjust = 0.4),
                   show.legend = FALSE,size=5)+
  scale_fill_viridis_d(alpha=0.7)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        text = element_text(size=14))#+
  #labs(title="The proportion of seed-dispersal networks \nper biogeographical realm")
i0

#### 2.2) What is the number of interactions per realm -------------------------------------------
netsAna |> 
  group_by(Realm) |> 
  count() |> 
  drop_na() |> 
  distinct() |> 
  ungroup() |> 
  mutate(perc=(perc = `n` / sum(`n`))) |> 
  mutate(labels = scales::percent(perc)) |> 
  rename(int = n) -> Nints

print(Nints)

# total number of interactions
sum(Nints$int)

# plot the proportional number of interactions per realm
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
  geom_label_repel(aes(label = labels),
                   position = position_stack(vjust = 0.2),
                   show.legend = FALSE,size=5)+
  scale_fill_viridis_d(alpha=0.7)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        text = element_text(size=14))#+
  #labs(title="The proportion of interactions \nper biogeographical realm")
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
  geom_label_repel(aes(label = labels),
                   position = position_stack(vjust = 0.2),
                   show.legend = FALSE,size=5)+
  scale_fill_viridis_d(alpha = 0.7)+
  guides(fill = guide_legend(title = ""))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5),
        text = element_text(size=14))#+
  #labs(title="The normalized proportion of interactions \nper biogeographical realm")
i2

### Make a composite graph --------------------------------------------------
# This will be Figure 5 of the manuscritp
i0+i1+i2+
  plot_annotation(tag_levels = "a",tag_suffix = ")")+
  plot_layout(guides = 'collect')&
  theme(legend.position = "bottom",legend.text = element_text(size=14))

#ggsave(here("products","fig1_vertical.jpeg"),units="cm", height = 30, width = 15)
ggsave(here("products","fig1.jpeg"),units="cm", height = 13, width = 25)

## 2.3) Sort species by their frequency in the networks ------------------------

common <- netsAna |> 
  group_by(Species, Realm) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup()

## 2.3.1) Sort species by their frequency in the networks per realm -----------------------
re <- common |> 
  dplyr::select(Realm) |> 
  distinct() |> 
  pull()

l <- list()
for (i in re){
  l[[i]] <- common |> filter(Realm==i)
}
intSP <- map_df(l, ~.x)

# make a graphic with the top ten species per realm
# This will be the figure 3
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

# Standardize the number of and proportion of species per population size
# The data below was obtained from Callaghan et al. 2021 paper in PNAS

callaghan <- read_csv(here("Abundances","all_species_summary_table.csv")) |> 
  select('Scientific name','Abundance estimate') |> 
  rename(Species = 'Scientific name', Abundances='Abundance estimate')

abu <- left_join(common, callaghan) |> drop_na()
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
ggsave(here("products","SI_1.jpeg"),units = "cm", height = 30, width = 20)

## 2.3.2) Compare population sizes of rare vs common species ------------------------
# This is graphic goes to SM

abu |> ggplot(aes(x=n,y=Abundances))+
  geom_point()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  geom_smooth(method = "lm")+
  #facet_wrap(.~Realm,scales="free",ncol=2)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=14),
        legend.position = "none")+
 # scale_y_log10()+
  scale_x_log10()+
  labs(x="Number of species interactions in seed-dispersal networks",
       y="Population abundances")
ggsave(here("products","SI_4.jpeg"),units = "cm", height = 15, width = 20)


## 2.4) Explore taxonomic gaps at family level -------------------------------------------
### 2.4.1) How many families in avonet? --------------------------------------------------
# Avonet cata was obtained from Tobias et al. 2022 Ecology Letters

avonet <- read_csv("AVONET1_BirdLife.csv") 

# select species only
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1, Family=Family1) |>
  dplyr::select(Species,Trophic.Niche, Family)

avonetSp |> select(Family) |> n_distinct()

### 2.4.2) How many families in birdbase? ------------------------------------------------

birdbase <- read_csv("BIRDBASE_afterGMM.csv")
birdbaseSp <- birdbase |> 
  dplyr::select(Species, Genus, `Primary Diet`,`Family HBW/BirdLife v9.1 (2024)`) |> 
  unite("Species",Genus:Species, sep =" ") 
birdbaseSp |> select(`Family HBW/BirdLife v9.1 (2024)`) |> n_distinct()

### 2.4.3) How many families in networks? -----------------------------------------------
netsAna_fa <- left_join(netsAna,avonetSp)
netsAna_fa |> select(Family) |> n_distinct()

### 2.4.4) what are the most common families? -------------------------------------------
netsAna_fa |> 
  group_by(Family) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup() |> 
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc, accuracy = 0.01)) 

### 2.4.5) what are the most common families by realm? -----------------------------------

netsAna_fa |> 
  group_by(Realm, Family) |> 
  count() |> 
  #arrange(desc(n)) |> 
  ungroup() |> 
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc, accuracy = 0.01)) 

# 3) Taxonomic gaps ------------------------------------------------------------

netsSp <- netsAna |> 
  select(Species) |> 
  distinct() |> 
  mutate(Nets = 1)

## 3.1) Compare data with Avonet -----------------------------------------------------------
### Filter 1 -----------------------------------------------------------------------------
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1, Family = Family1) |>
  dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche, Family)

avo_nets <- left_join(avonetSp,netsSp )

avo_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1) -> avonetspIN

length(unique(avonetspIN$Species))

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

##### (How many species out ?) -----------------------------------------------------
#how many species in
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1, Family = Family1) |>
  #dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche, Family)

avo_nets <- left_join(avonetSp,netsSp )

avo_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1) -> avonetspIN

length(unique(avonetspIN$Species))

#how many species out and what is their trophic niche
avo_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1) |> 
  dplyr::filter(!Trophic.Niche == "Omnivore") |> 
  dplyr::filter(!Trophic.Niche=="Frugivore") |> 
  dplyr::filter(!Trophic.Niche=="Granivore") -> pOUT

dplyr::filter(!Trophic.Niche == "Omnivore") |> 
  dplyr::filter(Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") 

length(unique(pOUT$Species))

y<- left_join(pOUT,avonetsSpRealm3)
length(y$Species)

y |> dplyr::group_by( Trophic.Niche ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2))

y |> dplyr::group_by( Realm ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2))

# R#1: It would be useful to have access (maybe as a figure) to the proportion of all birds  -----
#on each realm that are included in these networks, 
#as well as the proportions after excluding the non-frugivores and the omnivores 
#(so that they can be compared visually).

# Load realms
avonetsSpRealm3 <- read_csv("Avonet_realmFixed.csv")
avonetsSpRealmNF <- avonetsSpRealm3

# avonet with no filter
avonetSp_r1 <- avonet |> 
dplyr::mutate(Species = Species1, Family = Family1) |>
 # dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche, Family)

#how many species in networks 
avo_nets_r1 <- left_join(avonetSp_r1,netsSp ) |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1) 

length(avo_nets_r1$Species) #1518

# assign realms
x <- left_join(avo_nets_r1 ,avonetsSpRealm3) 
length(x$Species)

# plot by realm
avonetALL1 <- x |> 
  dplyr::group_by(Realm) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2)) |> 
  
  ggplot(aes(x=reorder(Realm,-n), 
             y=proportion,
             #fill = Nets,
             label=proportion)) +
  geom_col(stat="identity") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=12),
        legend.position = "none",
        axis.ticks = element_blank())+
  coord_flip()+
  geom_text(nudge_y= +5,
            color="black",
            size = 3,
            fontface="bold")+
  scale_y_continuous(limits=c(0,100))+
  labs(subtitle= "Number of species = 1523",
    x="Biogeograpihc realms",y="Proportion of network species in AVONET")
avonetALL1

#now filter our omnivors, frugivors and seed eaters

avo_nets_r1_filter <- avo_nets_r1 |> 
  dplyr::filter(!Trophic.Niche == "Omnivore") |> 
  dplyr::filter(Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") 

length(avo_nets_r1_filter$Species)

# assign realms
x <- left_join(avo_nets_r1_filter ,avonetsSpRealm3) #|>
  #distinct()
length(x$Species)

avonetFI1 <-x|> 
  filter(!is.na(Realm)) |> 
  dplyr::group_by(Realm) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2)) |> 
  
  ggplot(aes(x=reorder(Realm,-n), 
             y=proportion,
             #fill = Nets,
             label=proportion)) +
  geom_col(stat="identity") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=12),
        legend.position = "none",
        axis.ticks = element_blank())+
  coord_flip()+
  geom_text(nudge_y= +5,
            color="black",
            size = 3,
            fontface="bold")+
  scale_y_continuous(limits=c(0,100))+
  labs(subtitle= "Number of species = 516",
    x="Biogeograpihc realms",y="Proportion of network species in AVONET \nexcluding non-frugivores and the omnivores")
avonetFI1

# make a figure of species in and out after filtering

avonetALL1 + avonetFI1+
  plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(here("products", "anovetINandOUT.jpeg"), dpi=300, units = "cm", width = 25, height = 12)

#### 3.1.2) Compare data with Avonet By realm --------------------------------------------------------------

avonetsSpRealm3 <- read_csv("Avonet_realmFixed.csv")
avonetsSpRealmNF <- avonetsSpRealm3

avo_nets_realmNF <- left_join(avonetsSpRealmNF ,netsSp ) |> 
  dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") 


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

###Filter 2 ---------------------------------------------------------------------

avonetSp2 <- avonet |> 
  dplyr::mutate(Species = Species1) |>
  dplyr::filter(Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche) 

avo_netsF1 <- left_join(avonetSp2,netsSp )

length(unique(avonetSp2$Species))

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


## 3.3) Compare data with Birdbase ---------------------------------------------------------
netsSp |> n_distinct()
birdbaseSp |>  n_distinct()

# Apply filters 1 and 2

birdbaseSp <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  filter(`Primary Diet`=="Fruit" | `Primary Diet`=="Seed" | `Primary Diet`=="Omnivore" ) |> 
  unite("Species",Genus:Species, sep =" ") 

birdbaseSpNO <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  unite("Species",Genus:Species, sep =" ") |> 
  filter(`Primary Diet`=="Fruit" | `Primary Diet`=="Seed" ) |> 
  dplyr::select(Species)

##3.4) Get proportions ----------------------------------------------------------
### Filter 1---------------------------------------------------------------------
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

#### (How many species out?) ------------------------------------------
#how many species out and what is their trophic niche
outBirdbase <- anti_join(netsSp,birdbaseSp)
outBirdbase|> select(Species) |> n_distinct()


### Filter 2 ---------------------------------------------------------------------

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
ggsave(here("products", "fig4_new_r1.jpeg"),units="cm", width=25, height = 15)


((g1F1/g2F1)|g1rF1)  +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products", "fig4_Omnivore_r1.jpeg"),units="cm", width=25, height = 15)

((g1NO/g2NO)|g1rNO)  +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products", "fig4_NoOmnivore_r1.jpeg"),units="cm", width=25, height = 15)


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


## 4.1) calculate estimated richness using iNext -------------------------------
# Code build with the help of Claude Sonnet 4.6 

  bird <- netsAna |> 
  select(Species, Realm) |> 
  #filter(Realm=="Neotropical") |> 
  #distinct() |> 
  drop_na()


## 4.2) Get the total number of interactions per realm (sampling units) ---------
n_papers <- bird %>%
  group_by(Realm) %>%
  summarise(n_papers = n(), .groups = "drop")

##  4.3) Get the number of papers recording each species, per realm ----------
species_freq <- bird %>%
  group_by(Realm, Species) %>%
  summarise(freq = n(), .groups = "drop")

## iNEXT incidence_freq format: first element = total sampling units,
## remaining elements = detection frequency of each species

incidence_list <- lapply(unique(bird$Realm), function(r) {
  total <- n_papers$n_papers[n_papers$Realm == r]
  freqs <- species_freq$freq[species_freq$Realm == r]
  c(total, freqs)
})
names(incidence_list) <- unique(bird$Realm)

## 4.3) Run iNEXT -------------------------------------------------
# it takes some hours to run
out <- iNEXT(incidence_list,
             q = 0,                  # q = 0: species richness
             datatype = "incidence_freq",
             knots = 1000,           # dense grid so target values are captured
             se = TRUE,
             conf = 0.95,
             nboot = 50)

## Inspect sampling effort and coverage per realm
out$DataInfo %>% select(Assemblage, T, S.obs, SC)


## 4.4) Coverage-based rarefaction at a common completeness level ------
## Standardise all realms to the lowest observed coverage (SC = 0.98)
## so comparisons are not biased by unequal numbers of source interactions

target_SC <- 0.98

rarefied_cov <- out$iNextEst$coverage_based %>%
  filter(Method == "Rarefaction") %>%
  group_by(Assemblage) %>%
  slice(which.min(abs(SC - target_SC))) %>%
  select(Assemblage, SC, qD, qD.LCL, qD.UCL) %>%
  ungroup() %>%
  arrange(desc(qD))

print(rarefied_cov)

## Check whether any realm needed extrapolation to reach target_SC
## (i.e. observed coverage was below the target)
rarefied_cov %>%
  left_join(out$DataInfo %>% select(Assemblage, SC_obs = SC),
            by = "Assemblage") %>%
  mutate(extrapolated = SC_obs < target_SC)

## 4.5) Plot rarefaction / extrapolation curves ----------------

ggiNEXT(out, type = 1, datatype = "incidence_freq") +
  labs(title = "Sample-based rarefaction/extrapolation curves",
       x = "Number of papers",
       y = "Species richness") +
  theme_bw() +
  theme(legend.title = element_blank())

## Coverage-based curves (richness vs sample completeness)
ggiNEXT(out, type = 3, datatype = "incidence_freq") +
  labs(title = "Coverage-based rarefaction/extrapolation curves",
       x = "Sample coverage",
       y = "Species richness") +
  theme_bw() +
  theme(legend.title = element_blank())


## 4.6) Export results table -----------------------------------

write.csv(rarefied_cov, here("products","rarefied_richness_by_realm.csv"), row.names = FALSE)
