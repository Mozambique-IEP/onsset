onsset : Open Source Spatial Electrification Tool
=================================================

[![PyPI version](https://badge.fury.io/py/onsset.svg)](https://badge.fury.io/py/onsset)
[![Build Status](https://travis-ci.com/OnSSET/onsset.svg?branch=master)](https://travis-ci.com/OnSSET/onsset)
[![Coverage Status](https://coveralls.io/repos/github/OnSSET/onsset/badge.svg?branch=master)](https://coveralls.io/github/OnSSET/onsset?branch=master)
[![Documentation Status](https://readthedocs.org/projects/onsset/badge/?version=latest)](https://onsset.readthedocs.io/en/latest/?badge=latest)

# Scope

This repository contains the source code of the Open Source Spatial Electrification Tool
([OnSSET](http://www.onsset.org/)) adapted for the Mozambique IEP.

# Installation

### Requirements

OnSSET requires Python > 3.8 with the following packages installed:
- numpy
- pandas
- rasterio
- jupyter
- seaborn
- alphashape
- geojson
- click-log
- trimesh
- numba

### Install from GitHub

Download or clone the repository which includes all of the codes and the environment file, then open Anaconda prompt 
and run the following commands:

First, browse to the location where the downloaded and extracted files are found on your computer:
```
cd ...
```

Next, run the following command to install all the required packages:

```
conda env create --name moz_onsset_env --file gep_onsset_env.yml
```

Once the installation is complete, you can run the following two commands to activate the environment and start a 
Jupyter Notebook:

```
conda activate moz_onsset_env

jupyter notebook
```

## Contact
For more information regarding the tool, its functionality and implementation
please visit https://www.onsset.org.
