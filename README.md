## GEB_wsl.ebc

*wsl.ebc* corrects environmental bias from an observational dataset based on environmental stratification/clustering. A map of n clusters is generated based on input raster layers (in the article, climate predictors used in PPMs are employed as inputs). Following random equal-stratified sampling design (EBCe; Hirzel & Guisan, 2002), for each cluster, the number of observations per species (relative to all others) is artificially rescaled to the total number of observations found in the cluster presenting the highest observation density. Proportional-stratified sampling design (EBCp; Hirzel & Guisan, 2002) adds a second step: for each cluster, the number of observations per species may additionally be multiplied by the cluster's area (i.e. proportion of pixels in percentage relative to the study area), or by its logarithm (default; consensus between EBCe and EBCp). Resulting output indicates a new number of observations per cluster and species, that the function automatically sub-samples with replacement over the original observational dataset. This function may be used for presences and absences distinctively. 

Detailed documentation on how to use the corrective methods may be found here: <a href="https://8ginette8.github.io/GEB_wsl.ebc/">8ginette8.github.io/GEB_wsl.ebc</a>

## Citation

Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., & Thuiller, W. (2021). Novel methods to correct for observer and sampling bias in presence‐only species distribution models. Global Ecology and Biogeography, 30(11), 2312-2325. doi: <a href="https://doi.org/10.1111/geb.13383">10.1111/geb.13383</a>

## References


