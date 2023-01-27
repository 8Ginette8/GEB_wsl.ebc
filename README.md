## GEB_wsl.ebc

*wsl.ebc* corrects environmental bias from an observational dataset based on environmental stratification/clustering. A map of n clusters is generated based on input raster layers (in the article, climate predictors used in PPMs are employed as inputs). Following random equal-stratified sampling design (EBCe; Hirzel & Guisan, 2002), for each cluster, the number of observations per species (relative to all others) is artificially rescaled to the total number of observations found in the cluster presenting the highest observation density.

Proportional-stratified sampling design (EBCp; Hirzel & Guisan, 2002) adds a second step: for each cluster, the number of observations per species may additionally be multiplied by the cluster's area (i.e. proportion of pixels in percentage relative to the study area), or by its logarithm (default; consensus between EBCe and EBCp).

Resulting output indicates a new number of observations per cluster and species, that the function automatically sub-samples with replacement over the original observational dataset. This function may be used for presences and absences distinctively. 

Detailed documentation on how to use the corrective methods may be found here: <a href="https://8ginette8.github.io/GEB_wsl.ebc/">8ginette8.github.io/GEB_wsl.ebc</a>

## Citation

Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., & Thuiller, W. (2021). Novel methods to correct for observer and sampling bias in presence‐only species distribution models. Global Ecology and Biogeography, 30(11), 2312-2325. doi: <a href="https://doi.org/10.1111/geb.13383">10.1111/geb.13383</a>

## References

Hirzel, A. & Guisan, A. (2002) Which is the optimal sampling strategy for habitat suitability modelling. Ecological Modelling, 157, 331-341. doi: <a href="https://doi.org/10.1016/S0304-3800(02)00203-X">10.1016/S0304-3800(02)00203-X</a>

Reynolds, A.P., Richards, G., De La Iglesia, B. & Rayward-Smith, V.J. (2006) Clustering rules: A comparison of partitioning and hierarchical clustering algorithms. Journal of Mathematical Modelling and Algorithms, 5, 475-504. doi: <a href="https://doi.org/10.1007/s10852-005-9022-1">10.1007/s10852-005-9022-1</a>

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M. & Hornik, K. (2019) Cluster: cluster analysis basics and extensions. R package version, 1, 56. url: <a href="https://cran.r-project.org/web/packages/cluster/index.html">cluster - CRAN</a>

Schubert, E. & Rousseeuw, P.J. (2019) Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms. Amato G., Gennaro C., Oria V., Radovanovic M. (eds) Similarity Search and Applications. SISAP 2019. Lecture Notes in Computer Science, vol 11807. Springer, Cham., pp. 171-187. doi: <a href="https://doi.org/10.1016/j.is.2021.101804">10.1016/j.is.2021.101804</a>
