# nicaea.py module

# Martin Kilbinger 2016
# See http://cosmostat.org/nicaea


# Compability with python2.x for x>6
from __future__ import print_function


import re
import sys
import os
import glob

import numpy as np
from math import *
import pylab as plt
from bisect import bisect_left, bisect_right
from itertools import cycle

import warnings

import mkstuff as stuff
from mkstuff.mkplot import make_log_ticks

try:
    import pyfits as fits
except ImportError:
    pass
    try:
        from astropy.io import fits
    except ImportError:
        raise

from astropy.io import ascii
from astropy import units as u



class cl_theo:
    """Theoretical predictions of tomographic convergence power spectrum."""

    def __init__(self, ell, cl, nzbin=1, name='', dimensionless=False, dim=2):
        """Initialize power spectrum class.

        Parameters
        ----------
        ell: array of double
            angular Fourier modes
        cl: nell x nzbin*(nzbin)/2 -dimensional array of double
            tomographic power spectra
        name: string
            attach string to class, default empty string
        dimensionless: boolean
            if False (default), input cl's are raw P_kappa(ell), if True: input cl's
            are 'dimensionless', ell (ell+1) / (2pi) C(ell)
        dim: integer
            space dimension, 2 (default) or 3.

        Returns
        -------
        dat: class cl_theo
            Data
        """

        self.ell  = np.array(ell)
        self.cl   = np.array(cl)
        self.name = name
        self.dim  = dim

        if dim != 2 and dim != 3:
            raise ValueError('Invalid dimension dim={}'.format(self.dim))


        self.nzcorr = self.cl.shape[0]
        self.n      = self.ell.shape[0]
        if nzbin != -1:
            if self.nzcorr != int(nzbin * (nzbin+1) / 2):
                raise ValueError('Invalid number of correlations {} and number of z-bins {}'.
                                 format(self.nzcorr, nzbin))
            self.nzbin = nzbin
        else:
            # TODO calculate nzbin(nzcorr)
            raise ValueError('implement nzbin(nzcorr)!')
        if self.cl.shape[1] != self.n:
            raise ValueError('Incompatible length of ell ({}) and cl ({})'.format(self.cl.shape[1], len(cl)))

        if dimensionless is True:
            self.cl_transform(ell_factor_action='divide')


    def get_cl_ij(self, i_bin, j_bin):
        """Return power spectrum correlation corresponding to (i_bin, j_bin)

        Parameters
        ----------
        i_bin, j_bin: integers
            redshift bins

        Returns
        -------
        cl_ij: numpy array
            power spectrum for bins (i, j)
        """

        pass



    def get_cl_n(self, n, dimensionless=False):
        """Return nth power spectrum correlation

        Parameters
        ----------
        n: integers
            redshift correlation index
        dimensionless: boolean, default False
            if True, return dimensionless power spectrum

        Returns
        -------
        cl_n: numpy array
            nth power spectrum
        """

        if n >= self.nzcorr:
            raise ValueError('Requested C_l index {} has to be smaller than nzcorr={}'.format(n, self.nzcorr))

        if dimensionless is True:
            return self.cl[n] * self.get_ell_factor()
        else:
            return self.cl[n]



    def set_cl_n(self, cl, n, dimensionless=False):

        q = cl
        if dimensionless is True:
            q = q * self.ell_factor()
        self.cl[n] = q


    def get_min_max_cl(self, dimensionless=False):
        """Return minumum and maximum values of the Cl's.

        Parameters
        ----------
        dimensionless: boolean, default False
            if True, return dimensionless power spectrum

        Returns
        -------
        cl_min, cl_max: doubles
            min and max values
        """

        if dimensionless is True:
            q = self.cl * self.get_ell_factor()
        else:
            q = self.cl

        return q.min(), q.max()


    def get_ell_factor(self):
        """Return factor to transform between dimensional and dimensionless power spectrum.
    
        Parameters
        ----------
        None

        Returns
        -------
        ell_factor: numpy array of doubles
            dim=2: ell * (ell+1) / (2*pi), dim=3: k^3 / (2 pi^2)
        """
        
        if self.dim == 2:
            return self.ell * (self.ell+1) / (2 * np.pi)
        elif self.dim == 3:
            return self.ell**3 / (2 * np.pi**2)


    def cl_transform(self, ell_factor_action='multiply'):
        """Transforms power spectra from dimensional to dimensionless.
           Private method, should only be called once
           in __init__!

        Parameters
        ----------
        ell_factor_action: string
            'multiply' (default) or 'divide', to multiply with or divide by factor.
        dim: integer
            space dimension, 2 (default) or 3.

        Returns
        -------
        None
        """

        ell = self.ell
        f   = self.get_ell_factor()
        for n in range(self.nzcorr):
            cl = self.get_cl_n(n)
            if ell_factor_action == 'multiply':
                cl = cl * f
            elif ell_factor_action == 'divide':
                cl = cl / f
            else:
                raise ValueError('Invalid string value ell_factor=\'{}\''.format(ell_factor))
            self.set_cl_n(cl, n)
            
        return


    def print(self, fname_out=None, dimensionless=True):

        f = open(fname_out, 'w')
        for i in range(len(self.ell)):
            print(self.ell[i], file=f, end=' ')
            for n in range(self.nzcorr):
                cl_n_i = self.get_cl_n(n)[i]
                if dimensionless is True:
                    cl_n_i = cl_n_i * self.get_ell_factor()[i]
                print(cl_n_i, file=f, end=' ')
            print('', file=f)
        f.close()



    def plot(self, fname_out, dimensionless=True, log_y=True, savefig=True, range_x=None, range_y=None, ylabel=None, label=None, ls_start=0):

        #lines  = ['s', 'o', 'D', 'p', '^', '*']
        lines  = ['-', '--', ':', '-.']
        #lines = ['-']
        colors = ['red', 'blue', 'green', 'magenta', 'black', 'orange']

        if savefig is True:
            plt.clf()

        linecycler  = cycle(lines)
        colorcycler = cycle(colors)
        for c in range(ls_start):
           next(linecycler)
           next(colorcycler) 

        plt.xlabel('$\\ell$')
        ylabel_base = 'C_{{ij}}(\\ell)'
        if dimensionless is True:
            my_ylabel = '$\\ell (\\ell+1)/(2\\pi) {}$'.format(ylabel_base)
            f = self.get_ell_factor()
        else:
            my_ylabel = '${}$'.format(ylabel_base)
            f = 1

        if ylabel is not None:
            my_ylabel = ylabel

        plt.ylabel(my_ylabel)

        log_ell = np.log10(self.ell)
        n = 0
        for i_bin in range(0, self.nzbin):
            for j_bin in range(i_bin, self.nzbin):

                next_line  = '{}'.format(next(linecycler))
                next_color = next(colorcycler)
                if label is None:
                    my_label = '({}, {})'.format(i_bin, j_bin)  
                else:
                    my_label = label
                cl         = self.get_cl_n(n, dimensionless=dimensionless) 
                if log_y is True:
                    y = np.log10(cl)
                else:
                    y = cl
                plt.plot(log_ell, y, next_line, label=label, color=next_color)

                n += 1

        if savefig is True:
            plt.legend(loc='best')

        font = {'size'   : 12,
                'family' : 'sans-serif'}
        plt.rc('text', usetex=False)
        plt.rc('font', **font)


        if range_x is not None:
            log_xmin = np.log10(range_x[0])
            log_xmax = np.log10(range_x[1])
        else:
            log_xmin = log_ell[0]
            log_xmax = log_ell[-1]
        make_log_ticks('x', log_xmin, log_xmax)

        if log_y is True:
            log_min_max =  np.log10(self.get_min_max_cl(dimensionless=dimensionless))
            make_log_ticks('y', log_min_max[0], log_min_max[1])
        else:
            plt.ticklabel_format(style='sci', axis='y', scilimits=(-5,5))

        if range_y is not None:
            plt.ylim(range_y[0], range_y[1])

        if self.name is not None and savefig is True:
            plt.title(self.name)

        if savefig is True:
            plt.savefig('{}.pdf'.format(fname_out))




def create_link(dst, path_cosmo):
    """Create link to file 'path_cosmo/dst'.

    Parameters
    ----------

    dst: string
        Destination file name
    path_cosmo: string
        Directory name

    Returns
    -------
    None
    """

    if os.path.isfile(dst):
        pass
    else:
        src = '{}/{}'.format(path_cosmo, dst)
        if os.path.isfile(src):
            os.symlink(src, dst)
        else:
            raise IOError('File {} not found at dir {}'.format(dst, path_cosmo))



def create_links_to_cosmo(path_to_cosmo):
    """Create links to all cosmo files.

    Parameters
    ----------
    path_to_cosmo: string
        Directory name

    Returns
    -------
    None
    """

    files = ['cosmo.par', 'cosmo_lens.par', 'nofz.par']
    for dst in files:
        create_link(dst, path_to_cosmo)

    files = glob.glob('{}/nofz_*'.format(path_to_cosmo))
    for dst in files:
        create_link(os.path.basename(dst), path_to_cosmo)


