# Characterisation and modelling of the native and invaded niche of castor bean (Ricinus communis L.): the Iberian Peninsula as a case study. Script and dataset.
This repository contains the script and archives used in the study “Characterisation and modelling of the native and invaded niche of castor bean (*Ricinus communis* L.): the Iberian Peninsula as a case study” by Echarren-Lucendo et al. (2026, draft). Please, read "README" for more information

# Ricinus.R
The script is divided into four sections: 

  ## PREPARING THE DATA FOR ANALYSIS
  This part prepares climatic and spatial data for niche analysis of Ricinus communis by integrating species occurrence records with bioclimatic variables. It downloads WorldClim bioclimatic layers, extracts climate values at occurrence locations from both the invaded range (Iberian Peninsula) and the native range (northeastern Africa), and filters Iberian occurrences based on minimum temperature thresholds to exclude frost-prone sites. The script then generates exploratory maps of occurrences and climatic variables at regional and global scales. Finally, it defines the native and invaded geographic extents, extracts background (environmental) climate data for each region by cropping and masking the bioclimatic rasters, and exports both occurrence-based and background climate datasets for subsequent ecological niche modelling.

  
  ## CHARACTERIZATION ANALYSIS
  Performs the climatic niche characterization and comparison of native versus invaded populations of *Ricinus communis* using the ecospat framework. It integrates occurrence data and background environmental data from both ranges, harmonizes bioclimatic variables, and projects them into a common environmental space using PCA. The script quantifies niche overlap (e.g. Schoener’s D and I indices), visualizes niche dynamics (overlap, expansion, stability, and unfilling), and assesses centroid shifts between native and invaded niches. It further implements niche equivalency and similarity tests based on randomization to evaluate hypotheses of niche conservatism or divergence. Finally, it maps potential climatic niches within each geographic range, providing a comprehensive assessment of niche structure and dynamics associated with the invasion process.

  
  ## NICHE MODELIZATION
  This section evaluates multicollinearity among bioclimatic variables, retains a subset of non-collinear predictors for MaxEnt modelling, and post-processes MaxEnt outputs for Ricinus communis occurrences. It extracts WorldClim bioclimatic values at native and invaded occurrence points, evaluates multicollinearity through correlation matrices and Variance Inflation Factor (VIF) analyses, and identifies a reduced set of ecologically meaningful, low-correlated variables for each region. The script then generates frost-free spatial masks based on minimum temperature thresholds and applies these masks to MaxEnt suitability rasters for the Iberian Peninsula and northeastern Africa. Finally, it crops and exports cleaned suitability layers and vector products for visualization and further analysis in GIS software such as QGIS.
  
  ## BOXPLOTS
This section compares the climatic conditions occupied by native and invaded populations of *Ricinus communis* using boxplots and basic statistical analyses.

If there is any issue, contact with: echarrenlucendonicolas@gmail.com










