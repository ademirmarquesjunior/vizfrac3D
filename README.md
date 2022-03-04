# Immersive Virtual Reality for Geological Interpretation and 3D Fracture Network Modeling

by Graciela Racolte<sup>1,2</sup>, Ademir Marques Jr<sup>1,2</sup>, Leonardo Scalco<sup>1,2</sup>, Tiago Siqueira de Miranda<sup>3</sup>, Daniel Capella Zanotta<sup>1,2</sup>, Caroline Cazarin <sup>4</sup>, Luiz Gonzaga Jr<sup>1,2</sup>, Mauricio Roberto Veronez<sup>1,2</sup>

<sup>1</sup>Vizlab | X-Reality and Geoinformatics Lab. Unisinos University. São Leopoldo, RS, Brazil
  
<sup>2</sup>Graduate Program in Applied Computing at UNISINOS, São Leopoldo, RS, Brazil

<sup>3</sup>Federal University of Pernambuco, Recife, Brazil
  
<sup>4</sup>CENPES - PETROBRAS, Rio de Janeiro, Brazil
  


Corresponding author: Graciela Racolte, email: gracielarr@unisinos.br


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

## Dependencies
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

## Running the files
Running the code `main.py` will perform the statistics and DFN models generation. The results of data file and figures can be found in the folder "output\\". 


## Credits	
This work is credited to the [Vizlab | X-Reality and GeoInformatics Lab](http://vizlab.unisinos.br/) and the following developers:	[Ademir Marques Junior](https://www.researchgate.net/profile/Ademir_Junior) and [Graciela dos Reis Racolte](https://www.researchgate.net/profile/Graciela-Racolte).


## License


All source code is made available under a BSD 3-clause license. You can freely use and modify the code, without warranty, so long as you provide attribution to the authors. See LICENSE.md for the full license text.

The manuscript text is not open source. The authors reserve the rights to the article content, which is currently submitted for publication in the Computers & Geosciences.

