### This is the repository of the article:

*Assessing the Eltonian shortfall for bird species in seed-dispersal networks*, in R2 review in the journal Biodiversity Informatics

### Citation: 
- Moulatlet G.M., Falcão J.,  Arregoitia L.V., Mendes S.B., Dáttilo W., and Villalobos F., Assessing the Eltonian shortfall for bird species in seed-dispersal networks. Submitted.
﻿
### Contact:
- Gabriel M. Moulatlet (mandaprogabriel@gmail.com)
﻿
### Abstract: 
- Seed dispersal interactions in mutualistic networks are important for maintaining and recovering many ecosystems; however, current empirical knowledge is uneven across regions and taxa Recent decades have seen an exponential interest in studying these interactions, especially t birds and plants seed-dispersal interactions, using ecological networks as key tools for addressing long-standing ecological questions. The use of seed dispersal networks may still represent limited information on this interaction with potential taxonomic and geographic gaps. Here, we quantify the extent of the current knowledge gap on the Avian frugivore assemblages of seed dispersal networks, describing the Eltonian shortfall for this important interaction. We compiled 538 bird-plant seed dispersal networks across all biogeographic realms and evaluated their taxonomic gap (i.e. the gap between the number of bird species recorded in seed-dispersal networks and those documented as seed dispersers in global avian ecology datasets) and geographic gap (i.e., the sampling coverage of seed-dispersal interactions within and among biogeographic regions). We found that 50% of bird dispersers in seed-dispersal networks have been classified as non-frugivorous in global avian datasets. 23%-24% of all potential seed dispersing birds are present in the networks globally, but this taxonomic gap is variable according to the species’ main food diet classification. Taxonomic gaps are larger in Australian, Oceanian and Madagascan realms, where only ~10% of the species from each biogeographic realm pool of bird seed dispersers have been recorded. We found that geographic sampling coverage and the number of interactions is largely skewed towards the Neotropics, Oriental and Palearctic realms (~60% of all networks), and that Madagascan, Oceania and Australian realms have 6% of the recorded networks. By highlighting key knowledge gaps and explicitly quantifying the size of the seed-dispersal Eltonian shortfall, we also provide future directions for more targeted data collection and also for species identification.
﻿
### Responsible for writing code: 
- Gabriel M. Moulatlet (mandaprogabriel@gmail.com)
﻿
### Folders and files:
- There is one data folder and one code file.
﻿
	* Data folder:
		* 	`DataForAnalysis_v15jul26_encodingFix_utf_NoDupli.csv`, with the following columns:
    		* `Species`: Bird species
			* `plant_id`: Plant species
   			* `database`: individual network identification
    		* `ref`: Original reference from where the network was obtained
    		* `interaction`: 1 if the bird and plant species interact with each other and 0 if not.
    		* `lat` and `long`: latitude and longitude coordinates
    		* `Trophic.niche`: Bird species diet, information obtained from AVONET
    		* `Realm`: Biogegraphic realm affiliation of each species.
    	
  		* 	AVONET and BIRDbase data can be found in the orginal publications - we are not allowed to share their data. The same applies to Holt et al. Biogeographic realms shapefile, which can downloaded from the original publication:

			- Holt, B. G., J.-P. Lessard, M. K. Borregaard, S. A. Fritz, M. B. Araújo, D. Dimitrov, P.-H. Fabre, C. H. Graham, G. R. Graves, K. A. Jønsson, and Others. 2013. An update of Wallace’s zoogeographic regions of the world. Science 339:74–78. American Association for the Advancement of Science.
			- Tobias, J. A., C. Sheard, A. L. Pigot, A. J. M. Devenish, J. Yang, F. Sayol, M. H. C. Neate-Clegg, N. Alioravainen, T. L. Weeks, R. A. Barber, and Others. 2022. AVONET: morphological, ecological and geographical data for all birds. Ecology letters 25:581–597. Wiley Online Library.
			- Şekercioğlu, Ç. H., K. D. Kittelberger, F. M. M. Mota, A. N. Buxton, N. Orton, A. DeNiro, E. R. Buechley, J. J. Horns, J. D. Blount, J. Socci, and M. H. C. Neate-Clegg. 2025. BIRDBASE: A Global Dataset of Avian Biogeography, Conservation, Ecology and Life History Traits. Sci. Data 12:1558

	* Code file:
	    * `DataAnalysis_v26jun26.R`: contains the codes to generate the figures and the results of the manuscript

* Software version:
   - R version 4.5.1 (2025-06-13 ucrt)
   - Platform: x86_64-w64-mingw32/x64
   - Running under: Windows 11 x64 (build 26200)
   - Rstudio 2025.9.1.401 (Cucumberleaf Sunflower)

### Acknowledgments:
- We are thankful to Tom Bradfer-Lawrence, Clementine Durand-Bessart, and Manuel Nogales for sharing datasets on seed-dispersal networks compiled by them. GM and JF were supported by SECIHTI postdoctoral grants. LVA, FV and WD thank the Instituto de Ecologia A.C. for continuous institutional support.

