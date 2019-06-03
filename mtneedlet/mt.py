
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gafun
from scipy.stats import norm
import scipy.integrate as integrate
# import os
# import warnings
# import subprocess
import pandas as pd
# import numpy.ma as ma

 
def f_fromks(k1,k2):
    '''Get the theoretical maxima distribution f, from the values of k_1, k_2.
    
    For more information, see [1]_ or [2]_ (in the latter, k_1 corresponds to (\kappa_j)^2 and k_2 to (\eta_j)^2).
    
    Parameters
    ----------
    k1 : float
        The parameter k_1.
    k2 : float
        The parameter k_2.

    Returns
    -------
    function
        Theoretical distribution of the maxima, f.
        
    
    .. [1] Cheng, D., Cammarota, V., Fantaye, Y., Marinucci, D., & Schwartzman, A. (2016). Multiple testing of local maxima for detection of peaks on the (celestial) sphere. *Bernoulli, in press*, arXiv preprint arXiv:1602.08296.
    .. [2] Carron Duque, J., Buzzelli, A., Fantaye, Y., Marinucci, D., Schwartzman, A., & Vittorio, N. (2019). Point Source Detection and False Discovery Rate Control on CMB Maps. arXiv preprint arXiv:1902.06636.
    '''
    def f(x):
        fx=((2.*np.sqrt(3.+k1))/(2.+k1*np.sqrt(3.+k1))* 
        ((k1+k2*(x**2.-1.))*norm.pdf(x)*norm.cdf((x*np.sqrt(k2))/(np.sqrt(2.+k1-k2))) +
        np.sqrt(k2*(2.+k1-k2))/(2.*np.pi)*x*np.exp((-(2.+k1)*x**2.)/(2.*(2.+k1-k2))) +
        np.sqrt(2./(np.pi*(3.+k1-k2)))*np.exp((-(3.+k1)*x**2.)/(2.*(3.+k1-k2)))*
        norm.cdf((np.sqrt(k2)*x)/np.sqrt((2.+k1-k2)*(3.+k1-k2)))
        ))
        return(fx)
    return(f)


def f_fromcl(cls):
    '''Get the theoretical maxima distribution f, from the angular power spectra C_l of a map.
    
    For more information, see [1]_ or [2]_ (in the latter, k_1 corresponds to (\kappa_j)^2 and k_2 to (\eta_j)^2).
    
    Parameters
    ----------
    cls : np.ndarray
        Angular power spectrum of the map.

    Returns
    -------
    function
        Theoretical distribution of the maxima, f.
        
    
    .. [1] Cheng, D., Cammarota, V., Fantaye, Y., Marinucci, D., & Schwartzman, A. (2016). Multiple testing of local maxima for detection of peaks on the (celestial) sphere. *Bernoulli, in press*, arXiv preprint arXiv:1602.08296.
    .. [2] Carron Duque, J., Buzzelli, A., Fantaye, Y., Marinucci, D., Schwartzman, A., & Vittorio, N. (2019). Point Source Detection and False Discovery Rate Control on CMB Maps. arXiv preprint arXiv:1902.06636.
    '''
    ls=np.arange(len(cls))
    c1=np.sum((2.*ls+1.)*ls*(ls+1.)/(8.*np.pi)*cls)
    c2=np.sum((2.*ls+1.)*ls*(ls+1.)*(ls-1.)*(ls+2.)/(32.*np.pi)*cls)
    k1=c1/c2
    k2=c1**2./c2
    f=f_fromks(k1,k2)
    return(f)
    
    
def _getcs(gamma,n,p=1.):
    '''Returns c_{p,2n}(gamma)
    '''
    return(2.**(gamma/2. - 2. - n - 2.*p) * gafun(1.- gamma/2. + n +2.*p))


def f_fromSW(j,B,gamma=2.5,p=1):
    '''Get the theoretical maxima distribution for a Sachs-Wolfe-like spectra filtered with a Mexican needlet.
    
    For more information, see [1]_.
    
    Parameters
    ----------
    j : float
        Frequency ``j`` of the Mexican needlet
    B : float
        Parameter ``B`` of the Mexican needlet.
    gamma : float, optional
        Parameter ``gamma`` for the Sach-Wolfe spectra. In the CMB, 2 < gamma < 3.
    p : int, optional
        Order of the Mexican needlet
    
    Returns
    -------
    function
        Theoretical distribution of the maxima, f.
        
    
    .. [1] Cheng, D., Cammarota, V., Fantaye, Y., Marinucci, D., & Schwartzman, A. (2016). Multiple testing of local maxima for detection of peaks on the (celestial) sphere. *Bernoulli, in press*, arXiv preprint arXiv:1602.08296.
    '''
    cp0=_getcs(gamma,0.,p)
    cp2=_getcs(gamma,1.,p)
    cp4=_getcs(gamma,2.,p)
    k1=4.*(cp2/cp4)*(B**(-2.*j))
    k2=2.*(cp2)**2./(cp0*cp4)
    f=f_fromks(k1,k2)
    return(f)


def _pvalue(x,f,**kwargs):
    '''Not vectorised version, accept only one x'''
    return(integrate.quad(f,x,np.inf,**kwargs))


def pvalues(xvec,f,returnerror=False,**kwargs):
    '''Calculate the p-values for a certain maxima distribution f, of diferent values of the intensity.
    
    Parameters
    ----------
    xvec : np.ndarray
        Array with the intensities where the p-values will be calculated.
    f : function
        Theoretical distribution of the maxima, f.
    returnerror : bool, optional
        If ``True``, the output will contain both the p-value and the error of the integration.
    **kwargs : dict
        Additional arguments will be passed to np.integrate.quad for the calculation of the p-value.
        
    Returns
    -------
    np.ndarray
        Array containing the p-values calculated. If ``returnerror`` is ``True``, it will also containg the error of the integration.
    '''
    pvals=np.vectorize(_pvalue)
    if returnerror == True:
        return(pvals(xvec,f,**kwargs))
    else:
        return(pvals(xvec,f,**kwargs)[0])
    
    
def _max_getpvalue_exact(vec_max,f,**kwargs):
    '''Calculate the exact p-values for the given intensities.'''
    return(pvalues(vec_max,f,**kwargs))


def _max_getpvalue_approx(vec_max,f,step=0.05,**kwargs):
    '''Calculate the approximated p-values for the given intensities. It interpolates from exact values calculated at intervals of ``step``.'''
    xvec=np.arange(vec_max.min()-step,vec_max.max()+step,step)
    f_inxvec=pvalues(xvec,f,**kwargs)
    approx_pvalues=np.interp(vec_max,xvec,f_inxvec)
    return(approx_pvalues)


def max_getpvalue(maxima,f,n_exact=1000,step=0.05,correct=True,**kwargs):
    '''Get the p-values for the maxima for a given expected distribution.
    
    Parameters
    ----------
    maxima : pd.DataFrame
        Pandas DataFrame with the information of the maxima. Has to contain at least the column 'Intensity'.
    f : function
        Theoretical distribution of the maxima, f.
    n_exact : int, optional
        Number of maxima where the p-value is computed directly. After the ``n_exact`` most intense values, the p-values are computed 
        in an array with step ``step`` and then interpolated. This introduces a small error in the computation but greatly speeds up the
        procedure.
    step : float, optional
        Step of the array where the p-values are exactly computed after the first ``n_exact`` maxima.
    correct : bool, optional
        If ``True``, check that the pvalues are decreasing with intensity and, if it is not the case, change the value accordingly.
        
    Returns
    -------
    pd.DataFrame
        The same Pandas DataFrame as the input ``maxima``, but with an additional column ``pvalues``.
    '''
    maxima=maxima.sort_values(by='Intensity',ascending=False)

    max_exact=maxima.iloc[:n_exact]
    max_approx=maxima.iloc[n_exact:]
    
    pval_exact=_max_getpvalue_exact(max_exact['Intensity'],f,**kwargs)
    pval_approx=_max_getpvalue_approx(max_approx['Intensity'],f,step=step,**kwargs)
    
    pvals=np.concatenate((pval_exact,pval_approx))
    maxima['pvalue']=pvals
    
    if correct:
        for icorrect in np.argwhere(maxima['pvalue'].diff()<0.)[:,0]:       #correct the pvalues 
            maxima['pvalue'].iloc[icorrect]=maxima['pvalue'].iloc[icorrect-1]
    
    return(maxima)


def benjamini_hochberg(maxima,alpha,plot=False):
    '''Select a subset of the maxima to be candidates to Point Source.
    
    It applies the procedure in multiple testing called Benjamini-Hochberg procedure.
    
    Parameters
    ----------
    maxima : pd.DataFrame
        Pandas DataFrame with the information of the maxima. Has to contain at least the column 'pvalue'.
    alpha : float
        Parameter alpha of the Benjamini-Hochberg procedure. Should be between 0 and 1.
    plot : bool, optional
        If ``True``, the function prints the expected number of sources for the intensities of the maxima versus the actual 
        number of maxima at that intensity, along with the threshold (expected*alpha). Two plots are given: one with all the 
        detections plus 5; and one with all the maxima.
        
    Returns
    -------
    pd.DataFrame
        The same Pandas DataFrame as the input ``maxima``, but only with the candidates to be point sources.
    '''
    if 'pvalue' not in maxima:
        raise ValueError('maxima does not have a "pvalue" column, try running max_getpvalue(maxima,f) to calculate it.')
    pval=maxima['pvalue']
    expect=pval*pval.size
    limit=np.arange(1,pval.size+1)*alpha
    
    if np.sum(expect<limit) == 0:
        detect=0
    else:
        detect=np.argwhere(expect<limit)[-1][0]+1
    print(f'{detect} points have been reported as detections (alpha={alpha})')
    
    if plot == True:
        fig,[ax1,ax2]=plt.subplots(1,2,figsize=(10,5))
        ax1.plot(np.arange(detect+5)+1,expect[:detect+5],label='Measured')
        ax1.plot(np.arange(detect+5)+1,limit[:detect+5],label=f'Expected * {alpha}')
        lim1y=ax1.get_ylim()
        lim1x=ax1.get_xlim()
        ax1.plot(np.arange(detect+6),np.arange(detect+6),'--',color='grey',linewidth=1,label='Expected')
        ax1.set_ylim(lim1y)
        ax1.set_xlim(lim1x)
        ax1.legend()
        ax2.plot(np.arange(pval.size)+1,expect,label='Measured')
        ax2.plot(np.arange(pval.size)+1,limit,label=f'Expected * {alpha}')
        lim2y=ax2.get_ylim()
        lim2x=ax2.get_xlim()
        ax2.plot(np.arange(pval.size+1),np.arange(pval.size+1),'--',color='grey',linewidth=1,label='Expected')
        ax2.set_ylim(lim2y)
        ax2.set_xlim(lim2x)
        ax2.legend()
        plt.tight_layout()
        
    return(maxima[:detect])

