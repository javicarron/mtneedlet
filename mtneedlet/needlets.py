
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


def plot_bl(bl,newfigure=True,label=None):
    '''Plot the needlet filter b(l) in multipole space.
    
    Parameters
    ----------
    bl : np.ndarray
        A numpy array containing the values of a needlet filter b(l), starting with l=0. It
        can also be an array containing several needlet filters.
    newfigure : bool, optional
        If ``True``, the image will be created in a new Figure.
    label : str or list of strings, optional
        Label to show in the legend of the figure.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot.
    '''
    if newfigure:
        fig=plt.figure()
    else:
        fig=plt.gcf()
    ax=plt.gca()
    if bl.ndim ==1:
        bl=np.array(bl,ndmin=2)
        label=[label]
    if label is None:
        label=[None]*bl.shape[0]
    for ii in np.arange(bl.shape[0]):
        ax.plot(np.arange(len(bl[ii])),bl[ii], label=label[ii])
    ax.set_xlabel(r'$\ell$', fontsize=14)
    ax.set_ylabel(r'$b(\ell)$', fontsize=14)
    if label[0] is not None:
        plt.legend()
    return(fig)


def plot_blprofile(bl,newfigure=False,unit='min',label=None,**kwargs):
    '''Plot the profile of the needlet in pixel space.
    
    Plot the value of the needlet as a function of the angular distance between a given point and the center of the needlet
    
    Parameters
    ----------
    bl : np.ndarray
        A numpy array containing the values of a needlet filter b(l), starting with l=0.
    newfigure : bool, optional
        If ``True``, the image will be created in a new Figure.
    unit : str, optional
        The unit is one of the following: ``'min'`` for arcminutes (default), ``'sec'`` for arcseconds,
        ``'deg'`` for degrees or ``'rad'`` for radians
    label : str, optional
        Label to show in the legend of the figure.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot.
    '''
    if newfigure:
        fig=plt.figure()
    else:
        fig=plt.gcf()
        
    if bl.ndim ==1:
        bl=np.array(bl,ndmin=2)
        label=[label]
    if label is None:
        label=[None]*bl.shape[0]
    ls=np.arange(bl.shape[1])
    for ii in np.arange(bl.shape[0]):
        leg=np.polynomial.legendre.Legendre(bl[ii]*(ls*2.+1.)/(4.*np.pi))
        x_prev=np.arange(-180.*60.,180.*60.)
        y_prev=leg(np.cos(np.deg2rad(x_prev/60.)))
        limit=np.min((x_prev[np.argwhere(np.abs(y_prev/leg(1))>0.001)][-1,0]*1.5,180.*60.))
        x=np.linspace(-limit,limit,300)/60.

        if unit=='min':
            plt.plot(x*60.,leg(np.cos(np.deg2rad(x))),label=label[ii],**kwargs)
        elif unit=='deg':
            plt.plot(x,leg(np.cos(np.deg2rad(x))),label=label[ii],**kwargs)
        elif unit=='sec':
            plt.plot(x*3600.,leg(np.cos(np.deg2rad(x))),label=label[ii],**kwargs)
        elif unit=='rad':
            plt.plot(np.deg2rad(x),leg(np.cos(np.deg2rad(x))),label=label[ii],**kwargs)
            
    if unit=='min':
        plt.xlabel(r"$\alpha (x,\xi )$ [min]", fontsize=14)
    elif unit=='deg':
        plt.xlabel(r"$\alpha (x,\xi )$ [deg]", fontsize=14)
    elif unit=='sec':
        plt.xlabel(r"$\alpha (x,\xi )$ [sec]", fontsize=14)
    elif unit=='rad':
        plt.xlabel(r"$\alpha (x,\xi )$ [rad]", fontsize=14)
    else:
        raise ValueError('"unit" not recognised, valid units are "min", "deg", "sec" and "rad".')

    plt.ylabel(r'$\Psi_j\,(\langle x,\xi \rangle)$', fontsize=14)
    if label[0] is not None:
        plt.legend()
    plt.tight_layout()
    return(fig)


def __f_need(t):
    '''Auxiliar function f to define the standard needlet'''
    if t <= -1.:
        return(0.)
    elif t >= 1.:
        return(0.)
    else:
        return(np.exp(1./(t**2.-1.)))
    
def __psi(u):
    '''Auxiliar function psi to define the standard needlet'''
    return(integrate.quad(__f_need,-1,u)[0]/integrate.quad(__f_need,-1,1)[0])

def __phi(q,B):
    '''Auxiliar function phi to define the standard needlet'''
    B=float(B)
    if q < 0.:
        raise ValueError('The multipole should be a non-negative value')
    elif q <= 1./B:
        return(1.)
    elif q >= 1.:
        return(0)
    else:
        return(__psi(1.-(2.*B/(B-1.)*(q-1./B))))
    
def __b2_need(xi,B):
    '''Auxiliar function b^2 to define the standard needlet'''
    b2=__phi(xi/B,B)-__phi(xi,B)
    return(np.max([0.,b2]))  ## np.max in order to avoid negative roots due to precision errors

def standardneedlet(B,j,lmax):
    '''Return the needlet filter b(l) for a standard needlet with parameters ``B`` and ``j``.
    
    Parameters
    ----------
    B : float
        The parameter B of the needlet, should be larger that 1.
    j : int or np.ndarray
        The frequency j of the needlet. Can be an array with the values of ``j`` to be calculated.
    lmax : int
        The maximum value of the multipole l for which the filter will be calculated (included).

    Returns
    -------
    np.ndarray
        A numpy array containing the values of a needlet filter b(l), starting with l=0. If ``j`` is
        an array, it returns an array containing the filters for each frequency.
        
    Note
    ----
    Standard needlets are always normalised by construction: the sum (in frequencies ``j``) of the
    squares of the filters will be 1 for all multipole ell.
    '''
    ls=np.arange(lmax+1)
    j=np.array(j,ndmin=1)
    needs=[]
    bl2=np.vectorize(__b2_need)

    for jj in j:
        xi=(ls/B**jj)
        bl=np.sqrt(bl2(xi,B))
        needs.append(bl)
        
    return(np.squeeze(needs))


def mexicanneedlet(B,j,lmax,p=1,normalised=True):
    '''Return the needlet filter b(l) for a Mexican needlet with parameters ``B`` and ``j``.
    
    Parameters
    ----------
    B : float
        The parameter B of the needlet, should be larger that 1.
    j : int or np.ndarray
        The frequency j of the needlet. Can be an array with the values of ``j`` to be calculated.
    lmax : int
        The maximum value of the multipole l for which the filter will be calculated (included).
    p : int
        Order of the Mexican needlet.
    normalised : bool, optional
        If ``True``, the sum (in frequencies ``j``) of the squares of the filters will be 1 for all multipole ell.

    Returns
    -------
    np.ndarray
        A numpy array containing the values of a needlet filter b(l), starting with l=0. If ``j`` is
        an array, it returns an array containing the filters for each frequency.
    '''
    ls=np.arange(lmax+1)
    j=np.array(j,ndmin=1)
    needs=[]
    if normalised != True:
        for jj in j:
            u=(ls/B**jj)
            bl=u**(2.*p)*np.exp(-u**2.)
            needs.append(bl)
    else:
        K=np.zeros(lmax+1)
        jmax=np.log(5.*lmax)/np.log(B)
        for jj in np.arange(1,jmax):
            u=(ls/B**jj)
            bl=u**2.*np.exp(-u**2.)
            K=K+bl**2.
            if np.isin(jj,j):
                needs.append(bl)
        needs=needs/np.mean(K[int(lmax/3):int(2*lmax/3)])
    return(np.squeeze(needs))
    

def filtermap_fromalm(alm,bl,nside,returnalm=False):
    '''Filter the map (given by its alm) with an input needlet.
    
    Parameters
    ----------
    alm : np.ndarray
        The alm of a Healpy map of the sky, as given by hp.map2alm()
    bl : np.ndarray
        A numpy array containing the values of a needlet filter b(l), starting with l=0.
    nside : int
        The ``nside`` of the input map. The output map will also have this ``nside``.
    returnalm : bool, optional
        If ``True``, the function will also return the ``alm`` for the *filtered* map.

    Returns
    -------
    np.ndarray or [np.ndarray,np.ndarray]
        If ``returnalm=False``, filtered map in the Healpix format. If ``returnalm=True``, a list with the filtered map and its corresponding ``alm``.
    '''
    filtered_alm=hp.almxfl(alm,bl)/(np.sqrt(12)*nside)
    betas=hp.alm2map(filtered_alm,nside);
    if returnalm == True:
        return([betas,filtered_alm])
    else:
        return(betas)
    
    
def filtermap(mapp,bl,returnalm=False):
    '''Filter the map with an input needlet.
    
    Parameters
    ----------
    mapp : np.ndarray
        A Healpy map of the sky.
    bl : np.ndarray
        A numpy array containing the values of a needlet filter b(l), starting with l=0.
    returnalm : bool, optional
        If ``True``, the function will also return the ``alm`` for the *filtered* map.

    Returns
    -------
    np.ndarray or [np.ndarray,np.ndarray]
        If ``returnalm=False``, filtered map in the Healpix format. If ``returnalm=True``, a list with the filtered map and its corresponding ``alm``.
    '''
    alm=hp.map2alm(mapp)
    nside=hp.get_nside(mapp)
    return(filtermap_fromalm(alm,bl,nside,returnalm=returnalm))

    