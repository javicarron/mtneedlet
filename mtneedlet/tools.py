
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
# from scipy.special import gamma as gafun
# from scipy.stats import norm
# import scipy.integrate as integrate
# import os
import warnings
# import subprocess
import pandas as pd
# import numpy.ma as ma


# try:
#     import camb
#     usecamb=True
# except ModuleNotFoundError:
#     print('Module "camb" not found, please input the Cl you want to use')
#     usecamb=False


def _set_text(fig,text,hh,**kwargs):
    '''Auxiliar function to print text in an image'''
    fig.text(0.05, 1-hh, text,
        verticalalignment='bottom', horizontalalignment='left',
        transform=fig.transFigure, fontsize=15, **kwargs)
    
    
def maximum_info(maxima,n,mapp,nmax=None):
    '''Show a plot with a maximum and some information about it.
    
    Utility function to easily view the basic propierties of a specific maximum.
    
    Parameters
    ----------
    maxima : pd.DataFrame
        Pandas DataFrame with the information of the maxima. Has to contain at least the columns 'Pixel number' and 'Intensity'.
        If it contains 'pvalue', it will use it. Otherwise, the p-value will be set to ``np.nan``.
    n : int
        Row of the ``maxima`` with the information of the maximum to be shown.
    mapp : np.ndarray
        Healpix map to be shown in the image.
    nmax : int or None, optional
        Total number of maxima detected in the map. If ``None`` is introduced, it uses the length of the ``maxima`` table. Since this
        can affect the results (if ``maxima`` does not contain all the maxima), the lines affected are colored grey and the function 
        raises a warning.
    
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the information about the maximum.
    '''
    pixn=maxima['Pixel number'].iloc[n]
    intensity=maxima['Intensity'].iloc[n]
    try:
        pval=maxima['pvalue'].iloc[n]
    except KeyError:
        pval=np.nan
    nrows=1
    nside=hp.get_nside(mapp)
    (lon,lat)=hp.pix2ang(ipix=pixn,nside=nside,lonlat=True)

    fig,axs=plt.subplots(nrows,3, figsize=(15,6*nrows),squeeze=False)
    for ax in axs.flatten():
        ax.axis('off')
    _set_text(fig,f"Pixel number: {pixn}",0.1/nrows)
    _set_text(fig,f"Intensity: {intensity:.2f}",0.2/nrows)
    _set_text(fig,f"Longitude [deg]: {lon:.3f}",0.3/nrows)
    _set_text(fig,f"Latitude [deg]: {lat:.3f}",0.4/nrows)

    _set_text(fig,f"Index: {n+1}",0.52/nrows)
    _set_text(fig,f"pvalue: {pval:.2}",0.62/nrows)
    if nmax !=None:
        _set_text(fig,f"Expected number: {pval*nmax:.2}",0.72/nrows)
    else:
        _set_text(fig,f"Expected number: {pval*maxima.shape[0]:.2}",0.72/nrows,alpha=0.5)
        warnings.warn('The number of maxima nmax has not been introduced. Using the size of maxima instead. If this table does not cointain all the maxima in the map, the result of "Expected number" is not reliable.')
    
    _set_text(fig,f"Intensity in the map: {mapp[pixn]:.2f}",0.85/nrows)
    _set_text(fig,f"Std of the map: {np.std(mapp):.2f}",0.95/nrows)

    hp.mollview(mapp,sub=[nrows,3,3],margins=(0,0.12/nrows,0.01,0))
    hp.projplot(lon,lat,'o',color='orange',lonlat=True,fillstyle='none',markersize=8)
    hp.gnomview(mapp,rot=(lon,lat),sub=[nrows,3,2],margins=(0.02,0.02/nrows,0.02,0.02/nrows))
    
    return(fig,axs)