"""

:Name: populationWL.py

:Date: 2015, 2019

:Author: Martin Kilbinger, <martin.kilbinger@cea.fr>

:Description: Some useful classes and functions for population weak lensing (lensing by foreground objects such as galaxies, clusters, or voids)

:Package: mkstuff

"""


import glob
import os
import re

import numpy as np

from astropy.io import ascii, fits
from astropy.table import Table, Column

import matplotlib
matplotlib.use("Agg")
import pylab as plt

from mkstuff import athena
from mkstuff import mkstuff



class cluster_data:
    """Class for cluster data

    Parameters
    ----------

    name: string
        Cluster name
    ra: float
        Cluster right ascension
    dec: float
        Cluster declination
    z: float
        Cluster redshift
    M: float
        Cluster mass
    unit: string
        Unit of ra and dec, one of 'rad', 'deg', 'arcmin', 'arcsec'
    scale: float, optional, default=1
        scale parameter, e.g. radius
    """

    def __init__(self, name, ra, dec, z, M, unit, scale=1):

        self.n   = len(name)
        if self.n != len(ra) or self.n != len(dec) or self.n != len(z):
            mkstuff.error('Unequal lengths for cluster input data')

        self.name = [str(n) for n in name]
        self.ra   = athena.unit_to_rad(ra, unit)
        self.dec  = athena.unit_to_rad(dec, unit)
        self.z    = z
        self.M    = M
        self.unit = unit
        self.scale = scale


    def print_min_max(self, x_min, x_max, string):
        """Prints x_min and x_max together with string.
        """

        print('{} = {: .3f} ... {: .3f} {}'.format(string, x_min, x_max, 'rad'), end='')
        if self.unit != 'rad':
            print(' ({:g} ... {:g} {})'.format(athena.rad_to_unit(x_min, self.unit),
                                        athena.rad_to_unit(x_max, self.unit), self.unit), end='')
        print('')


    def print_extend(self):
        """Prints extends in ra, dec, and z to the screen
        """
        ra_min, ra_max = min(self.ra), max(self.ra)
        dec_min, dec_max = min(self.dec), max(self.dec)
        z_min, z_max = min(self.z), max(self.z)

        self.print_min_max(ra_min, ra_max, ' RA')
        self.print_min_max(dec_min, dec_max, 'DEC')
        print('  z = {:.2f} ... {:.2f}'.format(z_min, z_max))


    def get_index(self, name, verbose=False):
        """Returns index of cluster with name <name> (-1 if not found).
        """

        try:
            idx = self.name.index(name)
        except ValueError:
            if verbose == True:
                print('Cluster \'{}\' not found in list'.format(name))
            idx = -1
    
        return idx


class nofz:
    """Redshift information, a la nicaea
    """

    def __init__(self, fname=None, zdata=None, zbin=None, err_if_empty=False):

        if fname is None or fname == '':

            if zdata is None and zbin is None:

                if err_if_empty == True:
                    mkstuff.error('No redshift info given, not creating empty nofz')

                # No nofz input given, create empty nofz
                self.z    = None
                self.n    = None
                self.nz   = None
                self.type = 'none'
            elif zdata is None or zbin is None:
                mkstuff.error('Only one of zbin and zdata defined')
            else:
                # Input histogram
                self.n, self.z = np.histogram(zdata, zbin)

                # Right-most input bin nofz.z[-1] is correctly interpreted as last bin edge.
                # Returned n has one entry less than z.
                self.n = np.append(self.n, [0])
                self.nz   = self.z.shape[0] - 1
                self.type = 'hist'

        else:

            if zdata is not None or zbin is not None:
                mkstuff.error('Both n(z) input file name and zdata/zbin defined')
            else:
                # Read nofz from file
                dat, n, head = mkstuff.read_table(fname, count_comment=True)
                if re.search('hist', head[0]) is not None:
                    self.z    = dat[:,0]
                    self.n    = dat[:,1]
                    self.nz   = dat.shape[0] - 1
                    self.type = 'hist'
                else:
                    mkstuff.error('nofz: only \'hist\' type implemented')

        self.is_valid()


    def is_valid(self):
        """Checks whether nofz has valid format and content.
         """

        if self.n[-1] != 0:
            mkstuff.error('Wrong histogram format, frequency for rightmost bin corner has ' +
                              'to be 0, found instead {}'.format(self.n[-1]))
        if sum(self.n) == 0:
            mkstuff.error('sum of n(z) is zero')
        if max(self.z) == 0:
            mkstuff.error('All redshifts are zero')

        return True



def get_cluster_data(options):
    """Read cluster catalogue and return cluster info. Compact form, using options structure.
    """

    clusters = get_cluster_data_narg(options.cluster_cat, options.cluster_name, options.cluster_ra, options.cluster_dec,
            options.cluster_z, cluster_M=options.cluster_M, unit=options.unit, cluster_scale=options.cluster_scale, verbose=options.verbose)

    return clusters



def get_cluster_data_narg(cluster_cat, cluster_name, cluster_ra, cluster_dec, cluster_z, cluster_M=None,
                          cluster_scale=None, unit=None, verbose=False):
    """Read cluster catalogue and return cluster info. See get_cluster_data (compact form).
    """

    extension = os.path.splitext(cluster_cat)[1]
    if re.search('fit', extension):
        format = 'fits'
        hdu  = fits.open(cluster_cat)
        dat  = hdu[1].data
        if unit is None:
            unit = athena.get_unit_fits(hdu[1].header, verbose=verbose)
        else:
            unit = unit
    else:
        format = 'ascii'
        dat  = ascii.read(cluster_cat)
        if unit is None:
            unit = athena.get_unit(cluster_ra, verbose=verbose)
        else:
            unit = unit

    name = mkstuff.get_data_col(dat, cluster_name, format=format)
    ra   = mkstuff.get_data_col(dat, cluster_ra, format=format)
    dec  = mkstuff.get_data_col(dat, cluster_dec, format=format)
    z    = mkstuff.get_data_col(dat, cluster_z, format=format, action=None)
    if z is None:
        z = [-1] * len(dec)

    if cluster_M is not None:
        M = mkstuff.get_data_col(dat, cluster_M, format=format, action=None)
    else:
        M = None
    if M is None:
        M = [-1] * len(dec)

    if cluster_scale is not None:
        scale = mkstuff.get_data_col(dat, cluster_scale, format=format, action=None)
    else:
        scale = [1] * len(dec)

    clusters = cluster_data(name, ra, dec, z, M, unit, scale=scale)

    if verbose:
        #print('Extend of catalogue {}'.format(cluster_cat))
        ##clusters.print_extend()
        #print('Coordinates in {}'.format(unit))
        pass

    return clusters



def plot_histograms(dir, X, hist_type='logM', out_dir='plots', out_name_ext='', verbose=False):
    """Plot cluster mass and redshift histograms"""

    if max(X) < 0:
        # X is array of [-1] if no corresponding column in cluster data found (see get_cluster_data)
        if verbose is True:
            print('Columns for hist_type \'{}\' not valid/not found, skipping plot_histogram'.format(hist_type))
        return

    if hist_type == 'logM':
        # Mass
        x      = np.log10(X)
        Mman, Mexp = mkstuff.frexp10(X.mean())
        text   = 'mean $M = {} \\times 10^{{{}}}$'.format('%.1f' % Mman, Mexp)
        xlabel = '$\\log M$'
        create_histogram(dir, out_dir, x, xlabel, text, hist_type, name_ext=out_name_ext)

    elif hist_type == 'z':
        # Redshift
        x = X
        text = 'mean $z$ = {}'.format('%.2f' % X.mean())
        xlabel = '$z$'
        create_histogram(dir, out_dir, x, xlabel, text, hist_type, name_ext=out_name_ext)

    elif hist_type == 'SNR':
        # Signal-to-noise
        x = X
        text      = 'mean SNR = {}'.format('%.2f' % X.mean())
        xlabel    = 'SNR{}'.format(out_name_ext)
        create_histogram(dir, out_dir, x, xlabel, text, hist_type, name_ext=out_name_ext)

    else:
        mkstuff.error('Wrong hist_type \'{}\''.format(hist_type))



def create_histogram(dir, out_dir, x, xlabel, text, name_base, name_ext=''):

    plt.clf()
    ax = plt.subplot(1, 1, 1)
    plt.hist(x)

    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.text(0.95, 0.95, text, ha='right', transform=ax.transAxes)
    plt.title('{} clusters'.format(len(x)))

    dst_dir = '{}/{}'.format(dir, out_dir)
    mkstuff.mkdir_p(dst_dir)
    fig_name = '{}/hist_{}{}.pdf'.format(dst_dir, name_base, name_ext)
    print('Saving figure \'{}\''.format(fig_name))
    plt.savefig(fig_name)



def signal_to_noise(fname, col=1, method='var', fname_cov=''):
    """Return signal to noise ratio for WL measurements in file fname.

    Parameters
    ----------
    fname: string
        Input file name
    col: integer
        Column number for signal, default=1
    method: string
        Noise method, 'var' (default), 'cov'
    fname_cov: string
        Covariance file, optional

    Returns
    -------
    SNR: float
        signa to noise ratio
    """

    if method == 'var':
        dat  = ascii.read(fname)
        y    = dat[dat.keys()[col]]  # gt or DeltaSigma
        dy   = dat['sqrt_Dcor']
        snr  = athena.signal_to_noise_w2_var(y, dy)
    else:
        mkstuff.error('Method \'{}\' not (yet) implemented'.format(method))
 
    return snr


def get_n_files(fname='wgl_1.txt', in_dir='.'):
    """Return the number of cross-correlation ascii files in
       subdirectories.
    """

    files = glob.glob('{}/*/{}'.format(in_dir, fname))

    return len(files)



def stack(stack, n_cl=-1, fname='wgl_1.txt', subtype='wgl', verbose=False):
    """Stack shear profiles of clusters. Call 'mean_allcolumns.pl'.
       (Formerly in wl_cluster.py.)

    Parameters
    ----------
    stack    : string
        Stack type: 'physical', 'angular'
    n_cl     : integer
        Number of clusters for checking, default=-1. If < 0 no check is performed.
    fname    : string
        File name of wgl output file to be stacked
    subtype  : string
        'wgl' (default), 'wgl_res', 'DeltaSigma'
    verbose  : bool
        Verbose output

    Returns
    -------
    None
     """

    if stack is None:
        if verbose is True:
            print('No stacking')
        return

    # TODO: Only use clusters with non-zero profile, modify -t option.

    in_name = fname

    if n_cl < 0:
        t_flag = ''
    else:
        t_flag = '-t {}'.format(n_cl)

    if subtype == 'wgl' or subtype == 'DeltaSigma':
        w_flags = '-w 3'            # Weight column (w)
    elif subtype == 'wgl_res':
        w_flags = '-w 2 -f 1/x^2'   # rms column, use as w=1/rms^2
    else:
        mkstuff.error('Unknown subtype \'{}\''.format(subtype))

    cmd = 'mean_allcolumns.pl -k {} {} -r {}'.format(t_flag, w_flags, in_name)

    n, out_mgs, err_msg = mkstuff.run_cmd(cmd, verbose=verbose)

    if n == 0:
        mkstuff.error('No  file \'{}\' found for stacking, exiting populationWL.stack with error'.format(in_name))

    mkstuff.mkdir_p('stack')
    os.rename('{}.mean'.format(in_name), 'stack/{}.mean'.format(in_name))
    os.rename('{}.rms'.format(in_name), 'stack/{}.rms'.format(in_name))
    os.rename('{}.meanrms'.format(in_name), 'stack/{}.meanrms'.format(in_name))

    return n



def get_bg_cat_zmin(dir, lensing_cat_base, lensing_ext, z):
    """Return background (lensing) catalogue with zmin > z.

    Parameters
    ----------
    dir: string
        Directory of catalogues
    lensing_cat_base: string
        Catalogue base name
    lensing_ext: string
        Catalogue extension
    z: float
        Minimum redshift

    Returns
    -------
    cat_lens: string
        Catalogue name. '' if no match found.
    """

    # Find all lensing catalogues
    base    = '{}/{}_zmin_'.format(dir, lensing_cat_base)
    pattern = '{}*.{}'.format(base, lensing_ext)
    files = glob.glob(pattern)
    if len(files) == 0:
        pattern = '{}.{}'.format(dir, lensing_cat_base, lensing_ext)
        print(pattern)
        if os.path.isfile(pattern):
            return pattern
        else:
            mkstuff.error('No matching lensing catalogue slices \'{}\' found'.format(pattern))

    # Extract zmin from file name
    min_zmin = 99
    cat_lens = ''
    for f in files:
        m = re.findall('{}(.*).{}'.format(base, lensing_ext), f)
        if len(m) == 0:
            continue
        this_zmin = float(m[0])
        if this_zmin >= z:
            if min_zmin > this_zmin:
                min_zmin = this_zmin
                cat_lens = f

    return cat_lens
