
Needlets
========

In *MTNeedlet*, needlets are defined by their filtering function :math:`b ({\ell\over {B^j}})`. This is represented by 
a np.ndarray of dimension 1, with the values of the function for each :math:`\ell`, starting by :math:`\ell=0`.


There are three types of functions regarding needlets:

#. Functions to create needlets from their :math:`B` and :math:`j` parameters.
#. Functions to filter a spherical map.
#. Functions to visualize and compare different needlets.


.. currentmodule:: needlets

Needlet creation
----------------

.. autosummary::
   :toctree: functions

   standardneedlet
   mexicanneedlet


Map filtering
-------------

.. autosummary::
   :toctree: functions

   filtermap
   filtermap_fromalm


Needlet visualization
---------------------

.. autosummary::
   :toctree: functions

   plot_bl
   plot_blprofile
