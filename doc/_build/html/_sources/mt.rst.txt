Multiple Testing
================

Multiple Testing is implemented to analyse the characteristics of the maxima population as a whole.
First, the intensity of the maxima is compared to the expected distribution :math:`f`, and
the probability of a maxima happening at random in a Gaussian map is computed (as opposed to the 
maximum being produced by a point source): this is the *pvalue*. These can be studied individually
or as a population. For the latter, we implement the Benjamini-Hochberg procedure of multiple testing.

There are three types of functions regarding *pvalues* and multiple testing:

#. Functions to calculate the theoretical distribution of maxima :math:`f`.
#. Functions to calculate the *pvalue*.
#. A functions to apply the multiple testing approach.

.. currentmodule:: mt

Theoretical distribution :math:`f`
----------------------------------

.. autosummary::
   :toctree: functions

   f_fromks
   f_fromcl
   f_fromSW


*pvalues*
---------

.. autosummary::
   :toctree: functions

   pvalues
   max_getpvalue


Multiple testing
----------------

.. autosummary::
   :toctree: functions

   benjamini_hochberg

