# InterpolationIDF.jl

This repository provides functions for the interpolation of precipitation extremes on a large domain toward *Intensity-Duration-Frequency* (IDF) curve construction at unmonitored locations developed in Jalbert *et al.* (2022).

## Installation

The following julia command will install the package:

```julia
julia> import Pkg
julia> Pkg.add(url="https://github.com/JuliaExtremes/InterpolationIDF.jl")
```

:warning: **The unregisterred package GMRF.jl is also required.** To install it, run in the Julia package manager the following command: 

```julia
julia> Pkg.add(url="https://github.com/jojal5/GMRF.jl")
```    

## Data

This version uses the hourly reanalyzed precipitations from the Regional Deterministic Reforecast System (Gasset *et al.*, 2021) as the spatial covariate.


## Tutorial

### 1. Data loading and preparation

`idf_load`

### 2. Latent iGMRF definition

`GMRF.iGMRF`, `create_datastructure`

### 3. Markov Chain Monte Carlo for parameter estimation

`mcmc`, `trace_plot`

### 4. Computing a sample of the GEV parameters posterior distribution

`getGEV`

### 5. Cross validation

`get_ω̂`, `get_ω̄`

## References

* Jalbert, J., Genest, C., and Perreault, L. (2022). Interpolation of precipitation extremes on a large domain toward IDF curve construction at unmonitored locations. *Journal of Agricultural, Biological and Environmental Statistics*.
* Gasset, N., Fortin, V., Dimitrijevic, M., Carrera, M., Bilodeau, B., Muncaster, R., Gaborit, E., Roy, G., Pentcheva, N., Bulat, M., Wang, X., Pavlovic, R., Lespinas, F., Khedhaouiria, D., and Mai, J. (2021). A 10 km North American precipitation and land-surface reanalysis based on the GEM atmospheric model. *Hydrology and Earth System Sciences*, 25(9):4917–4945.
