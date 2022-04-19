# vizfrac3D


This repository contains the source code to conduct 3D fracture statistics and generate 3D stochastic DFN models. 


## Content
- main.py
Script that load the interpreted data obtained from outcrop immersive interpretation, compute and plot fracture statistiscs, and generate the deterministic and stochastic DFN models
- frac3D/file_handling.py
Module to load strike/dip data from files and load Mosis XP project files
- frac3D/fracture_clustering.py
Python module to cluster fracture plane data with spherical k-means and the elbow method using Fisher statistics
- frac3d/fracture_generator.py
Python module for stochastic 3D fracture generation using Fisher statistics and spacing betweeen fractures
- frac3D/fracture_stats.py
Python module with statitics functions to obtain directional statistics (Fisher) and other fracture related statistics
- frac3D/geometry.py
Python module with functions and classes to support 3D geometry operation (point distances, ray/triangle mesh intersection, object definitions, etc)
- frac3D/plot.py
Python module with functions to plot statistical graphs and the DFN models (fracture planes and cell fracture intensity model)

## Installation

This Python script requires a Python 3 environment and the following installed libraries as seen in the file requirements.txt.

The code has been tested using packages of:  
- Python (version 3.7)
- numpy
- matplotlib
- mplstereonet
- open3D
- colorsys
- spherical_kde
- sklearn
- scipy
- numba
- math
- time

## Usage
Running the code `main.py` will perform the statistics and DFN models generation. The results of data file and figures can be found in the folder "output\\". 


## Credits	
This work is credited to the [Vizlab | X-Reality and GeoInformatics Lab](http://vizlab.unisinos.br/) and the following developers:	[Ademir Marques Junior](https://www.researchgate.net/profile/Ademir_Junior) and [Graciela dos Reis Racolte](https://www.researchgate.net/profile/Graciela-Racolte).


## License

All source code is made available under a BSD 3-clause license. You can freely use and modify the code, without warranty, so long as you provide attribution to the authors. See LICENSE.md for the full license text.

## How to cite

If you find our work useful in your research please consider citing one of our papers:

```bash
To be published
```

