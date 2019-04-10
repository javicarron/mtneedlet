.. MTNeedlet documentation master file, created by
   sphinx-quickstart on Wed Mar 13 16:06:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MTNeedlet's documentation!
=====================================

MTNeedlet is a small Python package created to filter spherical (Healpix) maps and analyse the maxima population.
It has been developed with the CMB in mind, but it can be applied to other spherical maps. It pivots around three 
basic steps:

#. The calculation of several types of needlets and their possible use to filter maps.
#. The detection of maxima (or minima) on spherical maps, their visualization and basic analysis.
#. The multiple testing approach in order to detect anomalies in the maxima population of the maps
   with respect to the expected behaviour for a random Gaussian map.

The original code was created to reproduce the results on [1]_ and later used to produce the results on [2]_.
This software heavily relies on Healpy (the Healpix implementation for Python) to eficiently deal with 
spherical maps.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
    installation
   needlets
   maxima
   mt
   tools
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. [1] Cheng, D., Cammarota, V., Fantaye, Y., Marinucci, D., & Schwartzman, A. (2016). Multiple testing of local maxima for detection of peaks on the (celestial) sphere. *Bernoulli, in press*, arXiv preprint arXiv:1602.08296.
.. [2] Carron Duque, J., Buzzelli, A., Fantaye, Y., Marinucci, D., Schwartzman, A., & Vittorio, N. (2019). Point Source Detection and False Discovery Rate Control on CMB Maps. arXiv preprint arXiv:1902.06636.
    
