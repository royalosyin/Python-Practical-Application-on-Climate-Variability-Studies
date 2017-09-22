# Python & Practical Application on Climate Variability Studies
Main objective of this tutorial is the transference of know-how in practical applications and management of statistical tools commonly used to explore meteorological time series, focusing on applications to study issues related with the climate variability and climate change. This tutorial starts with some basic statistic for time series analysis as estimation of means, anomalies, standard deviation, correlations, arriving the estimation of particular climate indexes (Niño 3), detrending single time series and decomposition of time series, filtering, interpolation of climate variables on regular or irregular grids, leading modes of climate variability (EOF or HHT), signal processing in the climate system (spectral and wavelet analysis). In addition, this tutorial also deals with different data formats such as CSV, NetCDF, Binary, and matlab'mat, etc. It is assumed that you have basic knowledge and understanding of statistics and Python.

## Generic libraries for scientific analysis
The default Python library for dealing with large arrays of numeric data (e.g. four dimensional latitude/longitude/altitude/time data arrays) is numpy, while the netCDF data format is commonly used in Atmospheric science and Oceanography as it is convinient to store various variables of many dimensions. The default for reading and writing netCDF files is netCDF4 (the capability to read and write text files, including .csv, is built into numpy).
### What is NetCDF?
NetCDF is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data.
NetCDF was developed and is maintained at Unidata. Unidata provides data and software tools for use in geoscience education and research. The NetCDF homepage may be found at http://www.unidata.ucar.edu/software/netcdf/. The NetCDF source-code is hosted at GitHub, and may be found directly at http://github.com/Unidata/netcdf-c.
### How to deal with NetCDF and other data with Python?
we mainly use netCDF4-python, NumPy and SciPy to process NetCDF and other data formats.
* netCDF4-python
> netcdf4-python is a Python interface to the netCDF C library. netCDF version 4 has many features not found in earlier versions of the library and is implemented on top of HDF5. This module can read and write files in both the new netCDF 4 and the old netCDF 3 format, and can create files that are readable by HDF5 clients. The API modelled after Scientific.IO.NetCDF, and should be familiar to users of that module (see more http://unidata.github.io/netcdf4-python/).
* NumPy
> NumPy is the fundamental package for scientific computing with Python (see more http://www.numpy.org/). 
It contains among other things:
a powerful N-dimensional array object
sophisticated (broadcasting) functions
tools for integrating C/C++ and Fortran code
useful linear algebra, Fourier transform, and random number capabilities
Besides its obvious scientific uses, NumPy can also be used as an efficient multi-dimensional container of generic data. Arbitrary data-types can be defined. This allows NumPy to seamlessly and speedily integrate with a wide variety of databases.
* SciPy
> SciPy is a collection of mathematical algorithms and convenience functions built on the Numpy extension of Python (https://www.scipy.org/index.html). It adds significant power to the interactive Python session by providing the user with high-level commands and classes for manipulating and visualizing data. With SciPy an interactive Python session becomes a data-processing and system-prototyping environment rivaling systems such as MATLAB, IDL, Octave, R-Lab, and SciLab.
>The additional benefit of basing SciPy on Python is that this also makes a powerful programming language available for use in developing sophisticated programs and specialized applications. Scientific applications using SciPy benefit from the development of additional modules in numerous niches of the software landscape by developers across the world. Everything from parallel programming to web and data-base subroutines and classes have been made available to the Python programmer. All of this power is available in addition to the mathematical libraries in SciPy.
## Data Sources
We will use the data publicly available as possible.
The data are mainly downloaded from https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.derived.surfaceflux.html
## Visualization
### Matplotlib
Matplotlib is a Python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms. Matplotlib can be used in Python scripts, the Python and IPython shell, the jupyter notebook, web application servers, and four graphical user interface toolkits.
Matplotlib tries to make easy things easy and hard things possible. You can generate plots, histograms, power spectra, bar charts, errorcharts, scatterplots, etc., with just a few lines of code. For a sampling, see the screenshots, thumbnail gallery, and examples directory
For simple plotting the pyplot module provides a MATLAB-like interface, particularly when combined with IPython. For the power user, you have full control of line styles, font properties, axes properties, etc, via an object oriented interface or via a set of functions familiar to MATLAB users.
See more from https://matplotlib.org/
### Basemap
Basemap is a great tool for creating maps using python in a simple way. It’s a matplotlib extension, so it has got all its features to create data visualizations, and adds the geographical projections and some datasets to be able to plot coast lines, countries, and so on directly from the library.
Basemap has got some documentation, but some things are a bit more difficult to find. I started this documentation to extend a little the original documentation and examples, but it grew a little, and now covers many of the basemap possibilities.
See more from https://basemaptutorial.readthedocs.io/en/latest/
## More
We will mainly apply these mostly generic libraries to carry out data analysis step by step. The procedures or steps are common to the atmospheric and ocean sciences. Although other advanced libraries will simpilfy the procedures, the underlying ideas should be the same in essense.
We will also introduce more other libraries such as xarray and iris in the following parts.
