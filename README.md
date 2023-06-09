# InterpolationIDF.jl

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

This repository provides functions for the interpolation of precipitation extremes on a large domain toward *Intensity-Duration-Frequency* (IDF) curve construction at unmonitored locations developed in Jalbert *et al.* (2022).

## Installation

The following julia command will install the package:

```julia
julia> import Pkg
julia> Pkg.add(url="https://github.com/JuliaExtremes/InterpolationIDF.jl")
```

## Data

This version uses the hourly reanalyzed precipitations from the Regional Deterministic Reforecast System (*RDRS v2.1*, Gasset *et al.*, 2021) as the spatial covariate, available [here](https://github.com/julemai/CaSPAr/wiki/Available-products) and IDF data from ECCC's Engineering Climate Datasets can be found [here](https://collaboration.cmc.ec.gc.ca/cmc/climate/Engineer_Climate/IDF/).

## Tutorial

This quick tutorial presents the main features of InterpolationIDF.jl. More details can be found [here](https://github.com/jojal5/Publications/blob/master/JalbertGenestPerreault2022/IDF_interpolation%20-%20RDRS%20v2.1.ipynb).

### 1. Data loading and preparation

The steps of data preparation are as follows:
- Loading the gridded spatial covariate.
- Loading the station list where the IDF curves are available in Canada.
- Organizing the data and selecting the station if more than one lie in the same grid cell.

#### 1.1. Loading the gridded spatial covariate

The precipitation from *RDRS v2.1* are loaded for the 78 474 cells in a rectangular grid of $319 \times 246$ centered in Eastern Canada.

```julia
gridded_cov = CSV.read(GRIDDED_COV_PATH, copycols=true, DataFrame)  

m₁, m₂ = 319, 246
m = m₁*m₂

lat = reshape(gridded_cov[:,1], m₁, m₂)
lon = reshape(gridded_cov[:,2], m₁, m₂)
pr = reshape(gridded_cov[:,3], m₁, m₂);
```

#### 1.2. Preparing the station data list

Then, the station list where the IDF curves are available in Canada has to be loaded and processed. Only the stations located in the set PROVINCES are retained and if more than one station lie in the same grid cell, we keep only the one with the longest period of record.

Provided the paths for the right files, the duration and the provinces list, the function `get_station_list` will prepare the data and return a *DataFrame* with the information and data of each station.

```julia
station_list = get_station_list(GRIDDED_COV_PATH, STATION_INFO_PATH, STATION_DATA_PATH, DURATION, PROVINCES)
```

### 2. Latent iGMRF definition

The grid structure for the latent iGMRFs is defined using the `iGMRF` function from the package *[GMRF.jl](https://github.com/jojal5/GMRF.jl)*. The set of grid cells where the observations are is defined in $S$ and the remaining grid cells are defined in $\bar{S}$. Conditional independant subsets are also defined to improve MCMC performance.

```julia
G = GMRF.iGMRF(m₁,m₂,1,1).G
```

Then, we create a datastructure using the average daily precipitation from *RDRS v2.1* and the station elevation as the spatial covariates with the function `create_datastructure` :

```julia
datastructure = create_datastructure(G, station_list, m₁, m₂, log.(gridded_cov[:,:pr]), Float64.(station_list.Elevation))
```
It returns a `datastructure` type :

```
DataStructure
Y :		Vector{Vector{Float64}}[318]
X₁ᵢ :		Matrix{Float64}[(318, 2)]
X₂ᵢ :		Matrix{Float64}[(318, 2)]
G :
	GridStructure
	gridSize :	(319, 246)
	nbs :		Vector{Vector{Int64}}[78474]

S :
	GridPointStructure
	V :		Vector{Int64}[318]
	CondIndSubset :		Vector{Vector{Int64}}[2]
	CondIndIndex :		Vector{Vector{Int64}}[2]

S̄ :
	GridPointStructure
	V :		Vector{Int64}[78156]
	CondIndSubset :		Vector{Vector{Int64}}[2]
	CondIndIndex :		Vector{Vector{Int64}}[2]
```

### 3. Markov Chain Monte Carlo for parameter estimation

Markov Chain Monte Carlo procedure to generate a random sample of the posterior distribution can be executed using the function `mcmc`. Note that only the parameter values at the grid cells containing observations are retained to save memory. The values at the remaining grid cells can be interpolated offline.

```julia
C = mcmc(datastructure, niter=NITER, warmup=WARMUP, thin=THIN) 
```

It will return an object of type "Chains" from the package *MambaLite.jl*. To trace the chain and the density plot, the function `trace_plot` can be used.

### 4. Interpolation of the parameters at unmonitored locations

To interpolate the values at the remaining grid cells, the function `interpolation` can be used. Note that considering the size of the domain, it can be recommended to extrapolate the parameters to a single point of interest rather than to the whole domain.

### 5. Cross validation

Useful validation criterions to implement a cross validation can be obtained with the functions `get_ω̂` and `get_ω̄`.

## References

* Jalbert, J., Genest, C., and Perreault, L. (2022). Interpolation of precipitation extremes on a large domain toward IDF curve construction at unmonitored locations. *Journal of Agricultural, Biological and Environmental Statistics*.
* Gasset, N., Fortin, V., Dimitrijevic, M., Carrera, M., Bilodeau, B., Muncaster, R., Gaborit, E., Roy, G., Pentcheva, N., Bulat, M., Wang, X., Pavlovic, R., Lespinas, F., Khedhaouiria, D., and Mai, J. (2021). A 10 km North American precipitation and land-surface reanalysis based on the GEM atmospheric model. *Hydrology and Earth System Sciences*, 25(9):4917–4945.
