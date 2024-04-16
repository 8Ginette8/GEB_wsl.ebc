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

**obs**: Object of class matrix or data frame with three columns named "sp.id" (character), "x" (numeric) and "y" (numeric). More than one observation per species must be referenced. If only one species is in the data.frame, its sole observations will be resampled (equally or proportionally) across the environmental clusters.

**ras**: Object of class RasterBrick, RasterStack, list of RasterLayer of desired resolution and extent. Used to generate the map of clusters needed to summarize the environmental space of the study area.

**pportional**:  Logical. Should environmental stratification of observations be proportional to the clusters' areas? If TRUE, EBCp applies.

**plog**--> Logical. Should EBCp apply with a logarithm? If TRUE, a stratification consensus between EBCp and EBCe applies.

**nclust**--> Number of chosen clusters. Default is 50.

**sp.specific**--> Logical. Should EBC apply only for species whose environmental bias follows the overall one (i.e. the number of original species observations per cluster is correlated with that of the full dataset)?

**sp.cor** If sp.specific = TRUE, spearman's correlation tests are by default set to 0.5; i.e. species with r < 0.5 are excluded from the function outputs. If sp.cor = 0, all species are kept and corrected.

**keep.bias**--> Default is TRUE. Strongly recommended to use when sp.cor = TRUE. Per species, should the number of observations of the cluster in which the species was originally sampled most often, be preserved? Said differently, for each species, should the cluster with the most original species observations be as representative as the cluster with the most corrected observations? If TRUE, after EBC applies, the number of observations per species in their densest original cluster is set to that of the densest corrected cluster.

**filter**--> Logical. Should the observations be filtered according to 'ras' resolution?

**resample**--> Logical. If TRUE, the number of corrected observations per species will not exceed twice the number of original observations (resample without replacement)

**path**--> Path folder where the new species observation files should be saved.

**...**--> Additional arguments passed on to the ‘clara’ function (package 'cluster') 

### Return value

The function returns in ‘path’, one text file of corrected observations (presences or absences) per species. If the number of new EBC observations per species is too large, sampling those randomly without replacement before model calibrations is advised.

## 4. *wsl.ebc* function

``` r
wsl.ebc=function(obs=NULL,ras=NULL,pportional=TRUE,plog=TRUE,nclust=50,
	sp.specific=TRUE,sp.cor=0.5,keep.bias=TRUE,filter=FALSE,resample=FALSE,path=NULL,...)
{
	# Solve sp.cor = 0
	if (sp.cor %in% 0){
		sp.cor = -1
	}

	# 'nclust' must be at least equal to '2'
	if (nclust<=1 || is.null(nclust)){
		stop("'nclust' must be equal to 2 or more...!")
	}

	# 'obs' must be a data.frame of observations with column "sp.id","x" and "y"
	if (!(class(obs) %in% c("data.frame","matrix"))){
		stop("'obs' must be an object of class 'data.frame' or 'matrix'...!")
	}

	if (class(obs) %in% "matrix") {
		obs = as.data.frame(obs,stringsAsFactors=FALSE)
		obs$x = as.numeric(obs$x)
		obs$y = as.numeric(obs$y)
	}
	
	if (any(table(obs$sp.id) %in% 1)){
		stop("Each species must have more than one observation...!")
	}

	# 'ras' must be of 'brick' or 'stack' class or list of 'rasters'
	if (!(class(ras) %in% c("RasterBrick","RasterStack","list"))){
		stop("'ras' must be an object of class 'RasterBrick', 'RasterStack' or 'list'...!")
	} else if (length(ras) %in% 1 || nlayers(ras) %in% 1){
		stop("'ras' must include more than one 'RasterLayer'...!")
	}
	if (class(ras) %in% "list"){
		ras = stack(ras)
	}

	# First info message
	cat("Step 1 --> CLARA algorithm processing...","\n")
	
	# Convert raster in the right CLARA format
	id.toReplace = complete.cases(ras[])
	toClara = ras[][id.toReplace,]

	# Run CLARA and assign results to right pixels
	blocks = clara(toClara,nclust,...)
	toNew.ras = ras[[1]]
	toNew.ras[][id.toReplace] = blocks$clustering
	toNew.ras[][!id.toReplace] = NA

	# Start inventorying species observations per cluster 
	sp.ref = as.character(unique(obs$sp.id))

	# Second info message
	cat("Step 2 --> Species inventory...","\n")

	cluster.infos = list()
	cluster.tab = list()
	for (i in 1:length(sp.ref))
	{
		# Extract each species one by one & filter according to grid resolution if TRUE
		sp.target = obs[obs$sp.id %in% sp.ref[i],c("x","y")]
		if (filter){
			sp.target = wsl.obs.filter(sp.target[,c("x","y")],grid=toNew.ras)
		}

		# Extract cluster values with warnings if XY outside raster extent or assign to NA values
		clustV = extract(toNew.ras,sp.target)
		if(any(is.na(clustV))) {
			warning(paste("XY coordinates of",sp.ref[i],
				"outside 'ras' boundaries, or have extracted NAs environmental values...",sep=" "))
		}
		summaryV = table(clustV)

		# Fill a vector of length 'nclust'
		sp.vec = rep(0,nclust)
		sp.vec[1:nclust %in% as.numeric(names(summaryV))]=summaryV

		# Store infos
		cluster.infos[[i]] = cbind(sp.target,clustV)
		cluster.tab[[i]] = c(sp.ref[i],sp.vec)
	}

	# Combine and apply a sum to all columns
	cluster.all = as.data.frame(do.call("rbind",cluster.tab))
	cluster.all[,2:(nclust+1)] = sapply(cluster.all[,2:(nclust+1)],
		function(x) as.numeric(as.character(x)))
	cluster.sum = apply(cluster.all[,-1],2,sum)

	# Which cluster is max?
	plateau = max(cluster.sum)

	# Calculate area values based on n pixels per cluster (formated)
	pixels.clust = table(toNew.ras[])
	area.clust = pixels.clust/sum(pixels.clust)*100

	# Get reduction coef for each cluster specific to "plateau"
	reducoef = cluster.sum/plateau

	# Apply it to the species table and round up to the integer
	down.l = lapply(1:length(reducoef),function(x) cluster.all[,1+x]/reducoef[x])
	cluster.down = do.call("cbind",down.l)
	cluster.down = ceiling(cluster.down)

	# Weight observations by cluster areas
	if (pportional){
		cat ("...Environmental proportional stratification = TRUE...","\n")

		if (plog){
			cat ("...log...","\n")
			area.clust = log(area.clust+1)
		}

		down.l2 = lapply(1:length(reducoef),function(x) cluster.down[,x]*area.clust[x])
		cluster.down2 = do.call("cbind",down.l2)
		cluster.down = ceiling(cluster.down2)

	} else {
		cat ("...Environmental proportional stratification = FALSE...","\n")
	}

	# Replace NAs by 0 in case absolutely no informations fall in one cluster
	cluster.down[is.na(cluster.down)] = 0

	# Store
	cl.out = data.frame(x0=cluster.all[,1],cluster.down)
	cl.out[,1] = as.character(cl.out[,1])

	# Do we apply EBC only for target species?
	if (!(sp.specific)){
		cat("Step 3 --> EBC set to non species-specific...","\n")
	} else {
		cat("Step 3 --> EBC set to species-specific...","\n")

		# Tchao species whose environmental bias is not correlated with the general one
		spcl.temp=
		lapply(1:nrow(cluster.all),function(x) {
			sp.sumo = as.numeric(cluster.all[x,-1])
			cor.ebc = cor(sp.sumo,cluster.sum)
			if (cor.ebc>=sp.cor | is.na(cor.ebc)){
				return(cl.out[x,])
			}
		})
		# Format cluster.infos & bind results
		ebc.null = sapply(spcl.temp,is.null)
		cluster.infos = cluster.infos[!ebc.null]
		cl.out = do.call("rbind",spcl.temp)
	}

	# Should the initial environmental bias be preserved?
	if (!keep.bias){
		cat("Step 4 --> Initial environmental bias per species is not preserved...","\n")
	} else {
		cat("Step 4 --> Initial environmental bias per species is preserved...","\n")

		# Keeping initial biased environment
		spcl.temp2 =
		lapply(1:nrow(cl.out),function(x) {
			spcl.out = cl.out[x,]
			o.clust = cluster.all[cluster.all[,1] %in% spcl.out[,1],]
			maxo = which.max(o.clust[,-1])[1]
			spcl.out[,maxo+1] = max(spcl.out[,-1],na.rm=TRUE)
			return(spcl.out)
		})
		# Format cluster.infos & bind results
		cl.out = do.call("rbind",spcl.temp2)
	}

	# Fourth infos message
	cat("Step 5 --> Environmental Bias Correction (EBC) starting...","\n")
	
	# Write new observation files per species in choosen 'path'
	for (i in 1:nrow(cl.out))
	{
		# Extract number of obs. per cluster to sample
		sp.infos = cl.out[i,]
		sp.save = as.character(sp.infos[1])

		# Extract original observations from species
		cl.i = cluster.infos[[i]]

		# Fourth infos message
		cat("Processing ",sp.save,"...","\n")

		# Extract obs. according to each cluster
		Cselect=list()
		for (j in 1:nclust) {

			# Extract observations in clusters
			tg = cl.i[cl.i[,"clustV"] %in% j,c("x","y"),drop=FALSE]

			# Sample with replacements
			n.samp = as.numeric(sp.infos[-1][j])
			do.samp = tg[sample(1:nrow(tg),n.samp,replace=TRUE),c("x","y")]

			# Store with cluster ID depending on the output
			if (nrow(do.samp) %in% 0){
				next
			} else if (class(do.samp)[1] %in% "numeric") {
        		do.samp = matrix(do.samp,ncol=2)
        		colnames(do.samp) = c("x", "y")
    		}

			Cselect[[j]] = data.frame(cluster.id=paste0("clust.",j),do.samp,stringsAsFactors=FALSE)
		}

		# resample = TRUE ?
		new.obs = do.call("rbind",Cselect)
		if (resample & nrow(new.obs) > nrow(cl.i)*2) {
			new.obs = new.obs[sample(1:nrow(new.obs),nrow(cl.i)*2,replace=FALSE),]
		}

		# Save
		out.sp = paste0(path,"/",Sys.Date(),"_obs_corrected_",sp.save,".txt")
		write.table(new.obs,out.sp)
	}
}
```

## 5. Examples

``` r
# Load data
load("exrst.RData")
load("xy_ppm.RData")

# Run EBC function with the log consensus
# sp.specific = TRUE --> Select corrected XY outputs in a species-specific manner
# keep.bias = TRUE   --> Preserve initial observer bias of each species
wsl.ebc(obs = xy.ppm,
        ras = rst[[1:5]],
        pportional = TRUE,
        plog = TRUE,
        nclust = 50,
        sp.specific = TRUE,
        sp.cor = 0.5,
        filter = TRUE,
        keep.bias = TRUE,
        resample= = FALSE,
	path = getwd())
```
``` r
## Step 1 --> CLARA algorithm processing... 
## Step 2 --> Species inventory... 
## ...Environmental proportional stratification = TRUE... 
## ...log... 
## Step 3 --> EBC set to species-specific... 
## Step 4 --> Initial environmental bias per species is preserved... 
## Step 5 --> Environmental Bias Correction (EBC) starting... 
## Processing  sp1 ... 
## Processing  sp2 ... 
## Processing  sp3 ... 
## Processing  sp4 ... 
## Processing  sp5 ... 
## Processing  sp6 ... 
## Processing  sp7 ... 
## Processing  sp8 ... 
## Processing  sp9 ... 
## Processing  sp10 ...
```
``` r
# Open corrected observations
files = list.files(getwd())
target.files = files[grep("_obs_corrected_",files)]
correct.obs = lapply(target.files, function(x) obs=read.table(paste0(getwd(),"/",x)))
correct.obs = do.call("rbind",correct.obs)

# Extract with corrected and uncorrected observations environmental values of 'rst'
correct.env = extract(rst,correct.obs[,-1])
uncorrect.env = extract(rst,xy.ppm[,-1])

# Sampling of the environment (e.g. second layer) and of environmental values extracted with observations
f.correct.env = correct.env[,2]
f.uncorrect.env = uncorrect.env[,2]
first.env = rst[[2]][!is.na(rst[[2]])]
samp.env = first.env[sample(1:length(first.env),2e3)]
samp.vc=f.correct.env[sample(1:length(f.correct.env),2e3)]
samp.vnc=f.uncorrect.env[sample(1:length(f.uncorrect.env),2e3)]
 
# Plot distribution histograms for each 
par(mfrow = c(2,2))
hist(samp.env,breaks = seq(min(samp.env),max(samp.env),length.out=11)) # Cluster distribution across the region
hist(samp.vc,breaks = seq(min(samp.vc),max(samp.vc),length.out=11))  # Corrected obs. density (OD) per cluster
hist(samp.vnc,breaks = seq(min(samp.vnc),max(samp.vnc),length.out=11)) # Original/uncorrected OD per cluster
```
![Untitled](https://github.com/8Ginette8/GEB_wsl.ebc/assets/43674773/287e6d38-3189-4c12-89e4-47f69c9e638e)

## 6. Implementation of wsl.ebc function in our analysis
### Chosen parameters

*pportional*         TRUE (EBCp) and FALSE (EBCe)

*plog*                  FALSE

*nclust*                100

*sp.specific*         FALSE

*sp.cor*                NULL

*keep.bias*          FALSE

*filter*                   FALSE 

### Comments on the Discussion

**In the Discussion**: EBC should (1) only be applied if environmental bias is apparent, (2) ideally implemented in a species-specific manner, and (3) preserve the influence of the original environmental cluster within which the species was most sampled.

**Comment**: While EBC limitation (1) must be adressed by the user, EBC limitations (2) and (3) are adressed in wsl.ebc via the logical and correlation parameters sp.specific and sp.cor, and the logical parameter keep.bias respectively. (2) and (3) were not included in our analysis in order to simplify model and bias correction comparisons. Furthermore, keep.bias would only be efficient in preserving per species the environmental influence of its most sampled original cluster, if the selected number of clusters remains reasonable (by default, nclust = 50); i.e. summarizing not too precisely the environmental space of the study area, yet, large enough to account for the its environmental complexity. 

## 7. To cite

Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., & Thuiller, W. (2021). Novel methods to correct for observer and sampling bias in presence-only species distribution models. Global Ecology and Biogeography, 30, 2312–2325. https://doi.org/10.1111/geb.13383

## 8. Literature

 Hirzel, A. & Guisan, A. (2002) Which is the optimal sampling strategy for habitat suitability modelling. Ecological Modelling, 157, 331-341.

Reynolds, A.P., Richards, G., De La Iglesia, B. & Rayward-Smith, V.J. (2006) Clustering rules: A comparison of partitioning and hierarchical clustering algorithms. Journal of Mathematical Modelling and Algorithms, 5, 475-504.

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M. & Hornik, K. (2019) Cluster: cluster analysis basics and extensions. R package version, 1, 56.

Schubert, E. & Rousseeuw, P.J. (2019) Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms. Amato G., Gennaro C., Oria V., Radovanovic M. (eds) Similarity Search and Applications. SISAP 2019. Lecture Notes in Computer Science, vol 11807. Springer, Cham., pp. 171-187. 
