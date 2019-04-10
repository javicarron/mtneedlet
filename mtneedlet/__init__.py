'''Small package to treat with needlet filtering of spherical maps,
analysis of the maxima of the maps and multiple testing tools over
these maxima.
'''
#test

from .needlets import (plot_bl,
                      plot_blprofile,
                      standardneedlet,
                      mexicanneedlet,
                      filtermap_fromalm,
                      filtermap)

from .maxima import (hotspot,
                    locate_maxima,
                    max_inmask,
                    max_threshold,
                    plot_maxima)

from .mt import (f_fromks,
                f_fromcl,
                f_fromSW,
                pvalues,
                max_getpvalue,
                benjamini_hochberg)

from .tools import maximum_info

from .version import __version__

