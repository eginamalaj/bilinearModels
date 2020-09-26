# Phylogenetic bilinear models to predict toxicity of invertebrate species

This file contains the R code published with the article:

**Malaj, E.**, Guénard, G., Schäfer, R.B., and von der Ohe, P.C. (2016). Evolutionary patterns and physicochemical properties explain macroinvertebrate sensitivity to heavy metals. *Ecological Application* 26(4), 1249-1259, https://doi.org/10.1890/15-0346


Guillaum Guénard and collegues created the R Package **MPSEM** which is used to generate Phylogenetic Eigenvector Maps (PEM). The main paper explaining the published methodology can be found here: Guénard, G., Legendre, P., and Peres‐Neto, P. (2013). Phylogenetic eigenvector maps (PEM): a framework to model and predict species traits. *Methods in Ecology and Evolution* 4(12), 1120-1131,  https://doi.org/10.1111/2041-210X.12111

Summarized, in this work bilinear models integrating phylogenetic information of species and metal's physicochemical properties allowed to predict species sensitivity to metals. Combining the molecular information (DNA sequences) of 31 invertebrate species with the physicochemical properties of six bivalent metals, we built bilinear models that explained up to 70% of the variability in species sensitivity to metals. Phylogeny was the most important component of the bilinear models, as it explained the major part of the explained variance (>40%). Predicted values from bilinear modeling were in agreement with experimental values (>50%); therefore, this approach is a good starting point to build statistical models, which can potentially predict heavy metal toxicity for untested invertebrate species based on empirical values from similar species. 

To run this analysis the following files are needed from the `data` file:
1. `data_metal.RData` contains toxicity concentrations for metals. The dataset is compiled from USEPA ECOTOX dataset and the modelled toxicities are presented in another publication: **Malaj, E.**, Grote, M., Schäfer, R.B., Brack, W., and von der Ohe, P.C. (2012). Physiological sensitivity of freshwater macroinvertebrates to heavy metals. *Environ. Toxicol. Chem.* 31(8), 1754-1764.

2. `data_tree.rda` contains phylogenetic tree for invertebrate species; 

3. `data_exposure.RData` adds different exposure durations of the toxicity tests;

4. `data_physicochem.RData` adds different water parameters (pH, hardness, temp) of the toxicity tests;

5. `data_Properties.RData` adds physicochemical properties of each metal.


## Phylogenetic tree  

Multiple sequence alignment was performed for each gene using Muscle version 3.8.31. Individual genes were concatenated in a single super- alignment. Phylogenetic trees were based on the analysis of nucleotides with a maximum likelihood (ML) method. The branch support of the tree was assessed by 1000 bootstrap replicates. The majority consensus tree obtained from the bootstraps (`data_tree.rda`)looks like this:

![Phylogenetic tree](https://user-images.githubusercontent.com/54320408/94322270-35e63500-ff4f-11ea-91dc-59153bb8e227.jpeg)


Statistical support values are provided above each node (as a percentage). Asterisks denote nodes that received 100% bootstrap support. Nodes with support values less than 50% were collapsed.

## Results of the modelling

To evaluate the ability of the model to make accurate predictions, a leave-one-out cross-validation procedure was followed. This included:

1. Removing one species *i* 
2. Re-building the bilinear model for the *n−1* remaining species following:

⟨Y<sub>exp</sub>⟩=(Z⊗U)⟨B⟩+⟨E⟩

> where ‹…› denotes the lexicographic concatenation (unfolding) of the columns of a matrix into a single column vector, ⊗ is the Kronecker product, **Z** is a matrix of size (*j* × *l*), where *j* represents the metal and *l* represents the physicochemical properties, **U** is the influence matrix of size (*i* × *k*), *i* being the species and *k* representing the eigenvectors, **B** is a matrix of bilinear regression coefficients, and **E** is a matrix of error terms. 


3. Predicting the LC50 value for species *i*: 

Y<sub>pred</sub> = U<sub>target</sub>BZ<sup>T</sup>

> where **B** is a bilinear coefficient matrix, **U<sub>target</sub>** is thew calculated new (target) score matrix based on the position of the untested species in the phylogenetic tree, and **Z<sup>T</sup>** is the transposed matrix of **Z**.

4. Quantifying the difference between the predicted and experimental values for species *i*

![obs_vs_pred](https://user-images.githubusercontent.com/54320408/94323902-bbb8af00-ff54-11ea-95ca-6807a8822bc4.png)


In the figure above, the relationships between the experimental and predicted median lethal concentration (log10 LC50 in μg/L) for leave-one out cross-validation (A with four metals model and B with six metal model) is given. The dashed line is the 1:1 relationship that represents the perfect prediction by the model, and the solid line is the regression line. The prediction coefficient (q<sup>2</sup>) demonstrates the accuracy of the predicted values.






