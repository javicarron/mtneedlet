
Maxima
======

In *MTNeedlet*, maxima (or minima) are stored as a pd.DataFrame that contains at least two columns: "Pixel number" and
"Intensity". It can also contains other columns with the information of the location on the sky, or the *pvalue* asociated
to that intensity for an expected maxima distribution.

It is important to note that the information of the pixel number depends on the particular resolution (``nside``) of the 
map where the maxima where calculated, and therefore it is recommended to use another way to locate maxima if different resolutions
are going to be used (for example, the coordenates in the sky.

There are three types of functions regarding maxima:

#. A function to find the maxima (or minima) of a spherical map.
#. Functions to locate and analyse the maxima.
#. A function to visualize the population of the maxima.

.. currentmodule:: maxima

Maxima detection
----------------

.. autosummary::
   :toctree: functions

   hotspot


Analysis
--------

.. autosummary::
   :toctree: functions

   locate_maxima
   max_inmask
   max_threshold


Maxima visualization
--------------------

.. autosummary::
   :toctree: functions

   plot_maxima
