# <img src="https://github.com/ademirmarquesjunior/VizFrac3D/blob/main/output/logo.png" width="400" alt="Segmented image">

Fracture modeling plays an important role to understand the fluid flow in carbonate reservoirs. However, these reservoirs are often hard to access which makes it difficult to carry out direct surveys. Alternatively, analogue outcrops are a convenient alternative to studying carbonate reservoirs. Furthermore, the fracture characterization to generate Discrete Fracture Networks (DFNs) can take advantage of the analogue outcrop studies through Virtual Outcrop Models (VOMs), acquired by Unmanned Aerial Vehicles (UAV) and digital photogrammetry. 
The stochastic DFN generation is an important step in reservoir modeling as it brings more representative data to the process and has long been studied. However, optimizations concerning automatizing some of the steps necessary to its generation like data clustering are still open to advancements. In this sense, this work aims the fracture data clustering and the definition of the number of clusters when gathering data for the stochastic process, developing an Elbow method for spherical data and a balanced K-means, both based on Fisher statistics.
Regarding the clustering balance, our method achieved a lower standard deviation between sets while maintaining the Fisher values greater to obtain fracture sets with lower dispersion during the stochastic generation.

<img src="https://github.com/ademirmarquesjunior/VizFrac3D/blob/main/output/vizfrac.png" width="600" alt="Segmented image">



## Package overview
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
This work is credited to the [Vizlab | X-Reality and GeoInformatics Lab](http://vizlab.unisinos.br/) and the following developers:	[Ademir Marques Junior](https://www.researchgate.net/profile/Ademir_Junior) and [Graciela Racolte](https://www.researchgate.net/profile/Graciela-Racolte).


## License

All source code is made available under a BSD 3-clause license. You can freely use and modify the code, without warranty, so long as you provide attribution to the authors. See LICENSE.md for the full license text.

## How to cite

If you find our work useful in your research please consider citing one of the following papers:

```bash
@ARTICLE{9794644,
  author={Racolte, Graciela and Marques, Ademir and Scalco, Leonardo and Tonietto, Leandro and Zanotta, Daniel Capella and Cazarin, Caroline Lessio and Gonzaga, Luiz and Veronez, Maurício Roberto},
  journal={IEEE Access}, 
  title={Spherical K-Means and Elbow Method Optimizations With Fisher Statistics for 3D Stochastic DFN From Virtual Outcrop Models}, 
  year={2022},
  volume={10},
  number={},
  pages={63723-63735},
  doi={10.1109/ACCESS.2022.3182332}}
```

