# Appendix S2: EBC function, description, examples and considerations

## 1. Description of wsl.ebc

*wsl.ebc* corrects environmental bias from an observational dataset based on environmental stratification/clustering. A map of n clusters is generated based on input raster layers (in the article, climate predictors used in PPMs are employed as inputs). Following random equal-stratified sampling design (EBCe; Hirzel & Guisan, 2002), for each cluster, the number of observations per species (relative to all others) is artificially rescaled to the total number of observations found in the cluster presenting the highest observation density.

Proportional-stratified sampling design (EBCp; Hirzel & Guisan, 2002) adds a second step: for each cluster, the number of observations per species may additionally be multiplied by the cluster's area (i.e. proportion of pixels in percentage relative to the study area), or by its logarithm (default; consensus between EBCe and EBCp).

Resulting output indicates a new number of observations per cluster and species, that the function automatically sub-samples with replacement over the original observational dataset. This function may be used for presences and absences distinctively. 

## 2. R libraries and sub-function

``` r
# Load necessary libraries 
require(raster)
require(cluster)

# Load necessary sub-function 
wsl.obs.filter = function(o.xy,a.xy=NULL,grid)
{
  # Check 'o.xy' input

  if (ncol(o.xy)!=2 || !all(colnames(o.xy) %in% c("x","y"))){
    stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
  }

    # Check 'ras' input

    if (!(class(grid) %in% c("RasterBrick","RasterStack","RasterLayer"))) {
      stop("env.layer should be of class RasterBrick, RasterStack or RasterLayer!")
    }

    # Apply simple filtering
    if (!is.null(a.xy)) {
      
      # Check 'a.xy' input

      if (ncol(a.xy)!=2 || !all(colnames(a.xy) %in% c("x","y"))){
        stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
      }

      # Position presences and absences
      posP = raster::cellFromXY(grid,o.xy)
      posA = raster::cellFromXY(grid,a.xy)

      # Remove absences where we find presences
      posA = posA[!(posA %in% posP)]

      # Extract new presences/absences and regroup
      new.pxy = raster::coordinates(grid)[unique(posP),]
      new.axy = raster::coordinates(grid)[unique(posA),]
      new.oxy = list(new.pxy,new.axy)
      names(new.oxy) = c("Presences","Absences")

    } else {
      # Apply simple filtering
      posCELL = raster::cellFromXY(grid,o.xy)
      new.oxy = raster::coordinates(grid)[unique(posCELL),]
    }

    if (class(new.oxy)[1] %in% "numeric"){
      new.oxy = matrix(new.oxy,ncol=2)
      colnames(new.oxy) = c("x","y")
    }

    return(new.oxy)
}
```

## 3. Description of wsl.ebc parameters
### Usage

``` r
wsl.ebc = function(obs=NULL,ras=NULL,pportional=TRUE,plog=TRUE,nclust=50,
           sp.specific=TRUE,sp.cor=0.5,keep.bias=TRUE,filter=FALSE,resample=FALSE,path=NULL,...)
```

### Arguments

*obs*         Object of class matrix or data frame with three columns named "sp.id" (character),
              "x" (numeric) and "y" (numeric). More than one observation per species must be referenced. If
              only one species is in the data.frame, its sole observations will be resampled (equally or
              proportionally) across the environmental clusters.

*ras*         Object of class RasterBrick, RasterStack, list of RasterLayer of desired resolution and extent.
              Used to generate the map of clusters needed to summarize the environmental space of
              the study area.

*pportional*  Logical. Should environmental stratification of observations be proportional to the
              clusters' areas? If TRUE, EBCp applies.

@param plog Logical. Should EBCp apply with a logarithm? If TRUE, a stratification consensus
#' between EBCp and EBCe applies.
#' @param nclust Number of chosen clusters. Default is 50.
#' @param sp.specific Logical. Should EBC apply only for species whose environmental bias follows
#' the overall one (i.e. the number of original species observations per cluster is correlated
#' with that of the full dataset)?
#' @param sp.cor If sp.specific = TRUE, spearman's correlation tests are by default set to 0.5;
#' i.e. species with r < 0.5 are excluded from the function outputs. If sp.cor = 0, all species
#' are kept and corrected.
#' @param keep.bias Default is TRUE. Strongly recommended to use when sp.cor = TRUE. Per species,
#' should the number of observations of the cluster in which the species was originally sampled
#' most often, be preserved? Said differently, for each species, should the cluster with the most
#' original species observations be as representative as the cluster with the most corrected
#' observations? If TRUE, after EBC applies, the number of observations per species in their
#' densest original cluster is set to that of the densest corrected cluster.
#' @param filter Logical. Should the observations be filtered according to 'ras' resolution?
#' @param resample Logical. If TRUE, the number of corrected observations per species will not
#' exceed twice the number of original observations (resample without replacement)
#' @param path Path folder where the new species observation files should be saved.
#' @param ... Additional arguments passed on to the ‘clara’ function (package 'cluster') 
#' @return The function returns in ‘path’, one text file of corrected observations (presences
#' or absences) per species. If the number of new EBC observations per species is too large,
#' sampling those randomly without replacement before model calibrations is advised.
