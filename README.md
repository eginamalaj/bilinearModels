# Phylogenetic bilinear models to predict toxicity of invertebrate species

This file contains the R code published with the article:

Malaj, E., Guénard, G., Schäfer, R.B., and von der Ohe, P.C. (2016). Evolutionary patterns and physicochemical properties explain macroinvertebrate sensitivity to heavy metals. Ecological Application 26(4), 1249-1259.

Here bilinear models integrating phylogenetic information of species and physicochemical properties of metals allowed to predict species sensitivity to chemicals. Combining the molecular information (DNA sequences) of 31 invertebrate species with the physicochemical properties of six bivalent metals, we built bilinear models that explained 70–80% of the variability in species sensitivity to metals. Phylogeny was the most important component of the bilinear models, as it explained the major part of the explained variance (>40%). Predicted values from bilinear modeling were in agreement with experimental values (>50%); therefore, this approach is a good starting point to build statistical models, which can potentially predict heavy metal toxicity for untested invertebrate species based on empirical values for similar species. 

To run this analysis the following files are needed from the `/data` file:
1. `data_metal.RData` contains toxicity concentrations for metals. The dataset is compiled from USEPA ECOTOX dataset and the modelled toxicities are presented in another publication: Malaj, E., Grote, M., Schäfer, R.B., Brack, W., and von der Ohe, P.C. (2012). Physiological sensitivity of freshwater macroinvertebrates to heavy metals. Environ. Toxicol. Chem. 31(8), 1754-1764.

2. `data_tree.rda` contains phylogenetic tree for invertebrate species; 

3. `data_exposure.RData` adds different exposure durations of the toxicity tests;

4. `data_physicochem.RData` adds different water parameters (pH, hardness, temp) of the toxicity tests;

5. `data_Properties.RData` adds physicochemical properties of each metal;

