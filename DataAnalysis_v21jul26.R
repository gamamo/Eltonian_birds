# Assessing the Eltonian shortfall for bird species in seed dispersal networks
# Moulatlet et al.--------------------------------------------------

# 1) Load packages -----------------------------------------------------------
library(tidyverse)     # CRAN v2.0.0
library(here)          # CRAN v1.0.2
library(terra)         # CRAN v1.9-1
library(patchwork)     # CRAN v1.3.2
library(tidyterra)     # CRAN v1.1.0
library(rnaturalearth) # CRAN v1.2.0
library(ggrepel)       # CRAN v0.9.7
library(iNEXT)         # CRAN v3.0.2
library(flextable)     # CRAN v0.9.11
library(scales)        # CRAN v1.4.0
library(vegan)         # CRAN v2.7-3
library(sf)            # CRAN v1.1-0

## 1.2) Load the dataset --------------------------------------------------------

netsAna <- read_csv("DataForAnalysis_v21jul26_encodingFix_utf_NoDupli.csv")

## How many networks -----------------------------------------------------
netsAna |> 
  select(database) |> 
  n_distinct()

## How many unique bird species -----------------------------------------------
netsAna |> 
  select(Species) |> 
  n_distinct()

## 1.3) Make a map of the distribution of networks and interactions ---------------------

realms <- vect("newRealms.shp") # Holt et al.
realms <- realms[,"Realm"]
realms[realms$Realm == "Panamanian", "Realm"] <- "Neotropical"
realms[realms$Realm == "Oceanina", "Realm"] <- "Oceanian"

wrld <- ne_countries(continent = c("africa","south america","asia","oceania","europe",
                                   "north america"),
                     type = "countries")

m1 <- ggplot() +
  geom_spatvector(data=wrld, fill=NA)+
  geom_spatvector(data=realms, aes(fill=Realm))+
  geom_spatvector(data=nets_map,size=2,color="orange",shape=21,fill="black")+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(#plot.background = element_rect(fill="white",color="black"),
    legend.position = "none")
m1

# make a Cartogram with the size of the dots representing the number of interactions

nets_cartog <- netsAna |> 
  dplyr::select(database,lat,lon, interaction ) |> 
  dplyr::group_by(database, lat, lon) |> 
  dplyr::summarise(Interactions = sum(interaction))
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

# This will be figure 1:
(m1/m2) +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","productsToPaper","fig1.jpeg"), units = "cm", height = 20, width = 20,dpi=600)

# 2) Taxonomic gaps ------------------------------------------------------------

# Select intractions = 1
netsSp <- netsAna |> 
  select(Species) |> 
  distinct() |> 
  mutate(Nets = 1)

## 2.1) Compare data with Avonet -----------------------------------------------------------

# Load avonet
avonet <- read_csv("AVONET1_BirdLife.csv") 

# Filter 1 
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1, Family = Family1) |>
  dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche, Family)

# join avonet and networks
avo_nets <- left_join(avonetSp,netsSp )

# how many species selected?
avo_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1)

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

### 2.1.1) What happen when Psittacidae are removed from the species pool? ----
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1, Family = Family1) |>
  dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche, Family) |> 
  dplyr::filter(Family != "Psittacidae")

avo_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1)

# make a figure
g1.P <- avo_nets |> 
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
  labs(title="Representation in seed-dispersal networks \ncompared to AVONET \n (Without Psittacidae)")
g1.P

###2.1.2) How many species out ? -----------------------------------------------------
avo_nets |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1) |> 
  dplyr::filter(!Trophic.Niche == "Omnivore") |> 
  dplyr::filter(!Trophic.Niche=="Frugivore") |> 
  dplyr::filter(!Trophic.Niche=="Granivore") -> pOUT

#load Avonet realms
avonetsSpRealm3 <- read_csv("Avonet_realmFixed.csv")

y<- left_join(pOUT,avonetsSpRealm3)

#how many species out and what is their trophic niche
y |> dplyr::group_by( Trophic.Niche ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2))

#how many species out and what is their realms
y |> dplyr::group_by( Realm ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2))

#how many species out and what is their families
y |> dplyr::group_by( Family ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2)) |> 
  arrange(desc(proportion))

### 2.2) Compare data with Avonet By realm --------------------------------------------------------------

avonetsSpRealmNF <- read_csv("Avonet_realmFixed.csv")

avo_nets_realmNF <- left_join(avonetsSpRealmNF ,netsSp ) |> 
  dplyr::filter(Trophic.Niche == "Omnivore" | Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") 

# get proportions of interactions per realm
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

# make a figure of the proportion of interactions per realm
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

###2.2.1) Filter 2 ---------------------------------------------------------------------

avonetSp2 <- avonet |> 
  dplyr::mutate(Species = Species1) |>
  dplyr::filter(Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") |> 
  dplyr::select(Species,Trophic.Niche) 

avo_netsF2 <- left_join(avonetSp2,netsSp )


# make a figure of traxonomic proportions
g1NO <- avo_netsF2 |> 
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


# make a figure of the proportion of interactions per realm
avo_nets_realmNO <- left_join(avonetsSpRealmNF ,netsSp ) |> 
  dplyr::filter(Trophic.Niche=="Frugivore" | Trophic.Niche=="Granivore") 
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

g1rN0 <- avo_nets_realmNO |> 
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
g1rN0


## 2.3) Compare data with Birdbase ---------------------------------------------------------
birdbase <- read_csv("BIRDBASE_afterGMM.csv")
birdbaseSp <- birdbase |> 
  dplyr::select(Species, Genus, `Primary Diet`,`Family HBW/BirdLife v9.1 (2024)`) |> 
  unite("Species",Genus:Species, sep =" ") 
birdbaseSp |> select(`Family HBW/BirdLife v9.1 (2024)`) |> n_distinct()

netsSp |> n_distinct()
birdbaseSp |>  n_distinct()

# Apply filters 1 and 2
#filter 1
birdbaseSp <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  filter(`Primary Diet`=="Fruit" | `Primary Diet`=="Seed" | `Primary Diet`=="Omnivore" ) |> 
  unite("Species",Genus:Species, sep =" ") 

#filter 2
birdbaseSpNO <- birdbase |> dplyr::select(Species, Genus, `Primary Diet`) |> 
  unite("Species",Genus:Species, sep =" ") |> 
  filter(`Primary Diet`=="Fruit" | `Primary Diet`=="Seed" ) |> 
  dplyr::select(Species)

avonetsSpRealm3 <- read_csv("Avonet_realmFixed.csv")
avonetsSpRealmNF <- avonetsSpRealm3

# Join BirdBase and networks
bb_nets <- left_join(birdbaseSp,netsSp )

bb_nets |>dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                         TRUE ~ Nets )) |> 
  filter(Nets==1)

# make a figure of taxonomic proportions
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

### 2.3.1) How many species out? ------------------------------------------

birdbaseOut <- birdbase |> 
  dplyr::select(Species, Genus, `Primary Diet`,`Family HBW/BirdLife v9.1 (2024)`) |> 
  unite("Species",Genus:Species, sep =" ") 

bb_netsOut <- left_join(birdbaseOut,netsSp )

bb_netsOut |> 
  dplyr::mutate(Nets = case_when(is.na(Nets) ~ 0,
                                 TRUE ~ Nets )) |> 
  filter(Nets==1) |> 
  dplyr::filter(!`Primary Diet` == "Omnivore") |> 
  dplyr::filter(!`Primary Diet`=="Fruit") |> 
  dplyr::filter(!`Primary Diet`=="Seed") -> pOUT

length(unique(pOUT$Species))

y<- left_join(pOUT,avonetsSpRealm3)

#how many species out and what is their trophic niche
y |> dplyr::group_by(`Primary Diet` ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2))

#how many species out and what is their realms
y |> dplyr::group_by( Realm ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2))

#how many species out and what is their families
y |> dplyr::group_by( `Family HBW/BirdLife v9.1 (2024)` ) |> 
  count() |> 
  ungroup() |> 
  dplyr::mutate(total_n = sum(n),
                proportion = (n*100) / total_n) |> 
  mutate(proportion = round(proportion,2)) |> 
  arrange(desc(proportion))


### 2.3.2) Filter 2 ---------------------------------------------------------------------

bb_netsNO <- left_join(birdbaseSpNO,netsSp )

# make a figure of taxonomic proportions
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

#### Make a composite graph -------------------------------------------------

#Filter 1
((g1/g2)|g1rNF)  +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products","productsToPaper", "fig2.jpeg"),units="cm", width=25, height = 15,dpi=600)

#Filter 2
((g1NO/g2NO)|g1rN0)  +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
ggsave(here("products", "productsToPaper","figSI_1.jpeg"),units="cm", width=25, height = 15,dpi=600)


# 3) Geographic gaps ----------------------------------------------------
## 3.1) Sort species by their frequency in the networks ------------------------

common <- netsAna |> 
  group_by(Species, Realm) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup()

### 3.1.1) Sort species by their frequency in the networks per realm -----------------------
re <- common |> 
  dplyr::select(Realm) |> 
  distinct() |> 
  pull()

l <- list()
for (i in re){
  l[[i]] <- common |> filter(Realm==i)
}
intSP <- map_df(l, ~.x)

# make a graphic with the top ten species with more interactions per realm
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
ggsave(here("products","productsToPaper","Fig3.jpeg"),units="cm",height = 30,width = 38,dpi=600)

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
ggsave(here("products","productsToPaper","figSI_2.jpeg"),units = "cm", height = 30, width = 20,dpi=600)

### 3.1.2) Compare population sizes of rare vs common species ------------------------
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
ggsave(here("products","productsToPaper","figSI_3.jpeg"),units = "cm", height = 15, width = 20,dpi=600)

## 3.2) How many networks per biogeographic realm ---------------------------------------------
netsAna |> 
  select(Realm, database) |> 
  drop_na() |> 
  distinct() |> 
  group_by(,Realm) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup() -> Nnets

print(Nnets)

#### 3.2.1) What is the proportion of networks per realm ----------------------------------------

NnetsP <- Nnets |>   
  mutate(perc=(perc = `n` / sum(`n`))) |> 
  mutate(labels = scales::percent(perc))
print(NnetsP)

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
        text = element_text(size=14))
i0

#### 3.2.2) What is the number of interactions per realm -------------------------------------------
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
        text = element_text(size=14))
i1

#### 3.2.3) plot the normalized number of interactions per realm by number of networks ------
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
        text = element_text(size=14))
i2

##### Make a composite graph --------------------------------------------------
# This will be Figure 5 of the manuscript
i0+i1+i2+
  plot_annotation(tag_levels = "a",tag_suffix = ")")+
  plot_layout(guides = 'collect')&
  theme(legend.position = "bottom",legend.text = element_text(size=14))
ggsave(here("products","productsToPaper","fig5.jpeg"),units="cm", height = 13, width = 25,dpi=600)


## 3.3) Explore taxonomic gaps at family level -------------------------------------------
### 3.3.1) How many families in avonet? --------------------------------------------------

# select species only
avonetSp <- avonet |> 
  dplyr::mutate(Species = Species1, Family=Family1) |>
  dplyr::select(Species,Trophic.Niche, Family)

avonetSp |> select(Family) |> n_distinct()

### 3.3.2) How many families in birdbase? ------------------------------------------------

birdbaseSp <- birdbase |> 
  dplyr::select(Species, Genus, `Primary Diet`,`Family HBW/BirdLife v9.1 (2024)`) |> 
  unite("Species",Genus:Species, sep =" ") 
birdbaseSp |> select(`Family HBW/BirdLife v9.1 (2024)`) |> n_distinct()

### 3.3.3) How many families in networks? -----------------------------------------------
netsAna_fa <- left_join(netsAna,avonetSp)
netsAna_fa |> select(Family) |> n_distinct()

### 3.3.4) what are the most common families? -------------------------------------------
netsAna_fa |> 
  group_by(Family) |> 
  count() |> 
  arrange(desc(n)) |> 
  ungroup() |> 
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc, accuracy = 0.01)) 

### 3.3.5) what are the most common families by realm? -----------------------------------

netsAna_fa |> 
  group_by(Realm, Family) |> 
  count() |> 
  #arrange(desc(n)) |> 
  ungroup() |> 
  mutate(perc = `n` / sum(`n`)) |> 
  mutate(labels = scales::percent(perc, accuracy = 0.01)) 

# 4) Accumulation curves -----------------------------------------------------
# make a community table
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


### 4.1) Get the estimated richness --------------------------------------------------

bird <- netsAna |> 
  select(Species, Realm) |> 
  drop_na()

lst <- bird |> 
  count(Species,Realm) |> 
  group_by( Realm) |> 
  summarise(values = list(n), .groups = "drop") |> 
  tibble::deframe()

estimation <- iNEXT(lst, q=0, datatype="abundance") 

#### Export results table ---------------------------------------------
est_rich <- estimation$AsyEst
est_rich |> 
  filter(Diversity == "Species richness") |> 
  mutate(delta = Estimator - Observed) |> 
  flextable() %>%
  autofit() |> 
  save_as_docx(path=here("products","productsToPaper","TableS1.docx"))

####  Make a Figure ---------------------------------------------
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
ggsave(here("products","productsToPaper","fig4.jpeg"), units="cm", height = 15, width = 25,dpi=300)

## 4.2) calculate rarefaction curves  -------------------------------

  bird <- netsAna |> 
  select(Species, Realm) |> 
  drop_na()

# Get the total number of interactions per realm (sampling units) 
n_papers <- bird %>%
  group_by(Realm) %>%
  summarise(n_papers = n(), .groups = "drop")

# Get the number of papers recording each species, per realm 
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

# Run iNEXT 
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


## Coverage-based rarefaction at a common completeness level 
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

## Plot rarefaction / extrapolation curves

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


####Export results table -----------------------------------

write.csv(rarefied_cov, here("products","rarefied_richness_by_realm.csv"), row.names = FALSE)

# 5) Centroid vs species area Realm Assessement ------------------
# Load Birdlife polygons

dsn = "C:/Users/manda/Dropbox/postdocINECOL/birds distribution/BOTW/BOTW.gdb"

birds <- st_read(dsn = dsn, layer = "All_Species")


# get the list of species of the networks
netsAna |> 
  select(Species) |> 
  distinct() ->splist

# select the polygon of the network species
filterSP <- birds %>% filter(binomial %in% splist$Species)

# use some filters
filterSP.o12 <- filterSP %>% filter(origin<=2) # only native distribution

filterSP.surf <- filterSP.o12 %>% filter (st_is(., "MULTISURFACE")) #select multisurface types

filterSP.poly <- filterSP.o12 %>% filter(st_is(., "MULTIPOLYGON")) #select multipolygon types

# and save shapefiles
st_write(filterSP.surf,here("shpBirdLife", "surf.shp"),append=FALSE )
st_write(filterSP.poly,here("shpBirdLife", "poly.shp"),append=FALSE )

### Compare areas and realms ###################################################

avonetsSpRealm3 <- read_csv("Avonet_realmFixed.csv")

## Surface polygons --------------------------------

# 1. Load data
poly_b <- vect(here("shpBirdLife", "surf.shp"))  # target polygons (the ones you're assigning attributes to)
poly_b <- poly_b[poly_b$seasonl %in% c(1, 2), ]
#poly_b <- aggregate(poly_b, by="binomil")
as.data.frame(poly_b)

poly_a <- vect("newRealms.shp")  # source polygons (has the trait/attribute of interest)
poly_a[poly_a$Realm == "Panamanian", "Realm"] <- "Neotropical"
poly_a[poly_a$Realm == "Oceanina", "Realm"] <- "Oceanian"

# Replace with the actual column name holding the attribute you care about in poly_a
attr_col <- "Realm"

# Give poly_b a stable, explicit ID so we can track features through the intersection
poly_b$id_b <- 1:nrow(poly_b)
poly_b$area <- expanse(poly_b, unit="km")
as.data.frame(poly_b)

# 2. Intersect
inter <- intersect(poly_b, poly_a)
inter$InterArea <- expanse(inter, unit="km")  # area of each intersected piece

# 3. Keep the row with the largest overlapping area per poly_b polygon
inter_df <- as.data.frame(inter)
inter_df <- inter_df |> select(id_b,binomil, area,Realm,InterArea)
inter_df$prop <- (inter_df$InterArea/inter_df$area)*100


# 4. Get centroid of poly_b and find which poly_a polygon it falls in 
centroids_b <- centroids(poly_b)
extracted <- extract(poly_a, centroids_b)   # row order matches poly_b/centroids_b order

centroid_attr <- data.frame(
  id_b = poly_b$id_b,
  attr_centroid = extracted[[attr_col]]
)

# Compare and save
comparison <- merge(inter_df, centroid_attr, by = "id_b", all = TRUE)
write.csv(comparison,"centroid_comparision_surf.csv")

### Multi polygons ------------------------------

# 1. Load data
poly_b <- vect(here("shpBirdLife", "poly.shp"))  # target polygons (the ones you're assigning attributes to)
poly_b <- poly_b[poly_b$seasonl %in% c(1, 2), ] # selecting for breeding polygons
as.data.frame(poly_b)

poly_a <- vect("newRealms.shp")  # source polygons (has the trait/attribute of interest)
poly_a[poly_a$Realm == "Panamanian", "Realm"] <- "Neotropical"
poly_a[poly_a$Realm == "Oceanina", "Realm"] <- "Oceanian"

# Replace with the actual column name holding the attribute you care about in poly_a
attr_col <- "Realm"

# Give poly_b a stable, explicit ID so we can track features through the intersection
poly_b$id_b <- 1:nrow(poly_b)
poly_b$area <- expanse(poly_b, unit="km")
as.data.frame(poly_b)

# 2. Intersect 
inter <- intersect(poly_b, poly_a)
inter$InterArea <- expanse(inter, unit="km")  # area of each intersected piece

# 3. Keep the row with the largest overlapping area per poly_b polygon 
inter_df <- as.data.frame(inter)
inter_df <- inter_df |> select(id_b,binomil, area,Realm,InterArea)
inter_df$prop <- (inter_df$InterArea/inter_df$area)*100

# 4. Get centroid of poly_b and find which poly_a polygon it falls in 
centroids_b <- centroids(poly_b)
extracted <- extract(poly_a, centroids_b)   # row order matches poly_b/centroids_b order

centroid_attr <- data.frame(
  id_b = poly_b$id_b,
  attr_centroid = extracted[[attr_col]]
)

# 5. Compare 
comparison <- merge(inter_df, centroid_attr, by = "id_b", all = TRUE)
write.csv(comparison,"centroid_comparision_poly.csv")

### Compare for the paper ---------------------------------------------------

# join surface and and multi poly

surf <- read_csv("surf_polygons_cetroid_match.csv")
poly <- read_csv("poly_polygons_cetroid_match.csv")

all <- rbind(surf, poly)
all$prop <- round(all$prop,3)

length(unique(all$Species))

# Total distribution area per species = sum of DISTINCT polygon areas
total_area <- all %>%
  distinct(Species, area) %>%
  group_by(Species) %>%
  summarise(total_area = sum(area), .groups = "drop")

# Proportion of the species' total range intersecting each realm
result <- all %>%
  group_by(Species, Realm) %>%
  summarise(InterArea_sum = sum(InterArea), .groups = "drop") %>%
  left_join(total_area, by = "Species") %>%
  mutate(proportion = InterArea_sum / total_area)

result$proportion <- round(result$proportion ,3)*100

# load AVONET realms
avonetsSpRealm3 <- read_csv("Avonet_realmFixed.csv")
colnames(avonetsSpRealm3)[3] <- "Realm_avonet"
centroid <- all |> select(Species, attr_centroid)

result <- left_join(result,avonetsSpRealm3)
result$match <- NA
result$match[result$proportion >= 50 & result$Realm == result$Realm_avonet] <- "Match"
result$match[result$proportion <= 50 & result$Realm == result$Realm_avonet] <- "Match"
result$match[result$proportion <= 50 & result$Realm != result$Realm_avonet] <- "No Match"
result$match[result$proportion >= 50 & result$Realm != result$Realm_avonet] <- "No Match"

# and plot
result |> filter(proportion >=50) |>
  group_by(match) |> 
  count() |> 
  drop_na() |> 
  ungroup() |> 
  mutate(perc=(`n` / sum(`n`))) |> 
  mutate(labels = scales::percent(perc)) |> 
  ggplot(aes(x=as.factor(match), y=perc*100)) +
  geom_col(stat="identity",color="white") +
  theme_bw()+
  theme(text=element_text(size=14))+
  labs(x="",y="Class proportion")+
  scale_y_continuous(limits = c(0,100))
ggsave(here("products","Comparision_overlapRangesVSCentroid.jpeg"),dpi = 300,
       units="cm",height = 15,width = 20)
