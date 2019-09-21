"""

:Name: athena.py

:Author: Martin Kilbinger, <martin.kilbinger@cea.fr>

:Package: mkstuff

"""

# Compability with python2.x for x>6
from __future__ import print_function



import re

import sys
import os

import numpy as np
from math import *
import pylab as plt
from bisect import bisect_left, bisect_right

import warnings

from mkstuff import mkstuff
from mkstuff import mkplot

from astropy.io import fits
from astropy.io import ascii
from astropy import units as u


#################
### Constants ###
#################

# The following stuff is from pallas.py

# TODO: transform to class
units       = {'rad'    : 1,
               'deg'    : np.pi / 180.0,
               'arcmin' : np.pi / 180.0 / 60.0,
               'arcsec' : np.pi / 180.0 / 60.0 / 60.0,
               'none'   : 1}
unit_default = 'arcmin'


class col_descr:
    """Defines a pair (name, num) to describe a column name and number.
    """

    def __init__(self, name, num):
        self.name = name
        self.num  = num


class in_cols():
    """ Names and numbers of input columns. Names are used for fits (and ascii_col_names),
        numbers for ascii files.
    """

    # Angular scale
    scale   = col_descr('theta', 0)

    # Fourier scale
    scale_F = col_descr('ell', 0)

    # Signal
    signal1 = {'xi'    : col_descr('xi_p', 1),
               'gl'    : col_descr('g_t', 1),
               'w'     : col_descr('w', 1),
               'xi_res': col_descr('xi_p_resample', 1),
               'Pb'    : col_descr('P_E', 1),
              }

    # Second signal (only for some)
    signal2 = {'xi'    : col_descr('xi_m', 2),
               'xi_res': col_descr('xi_m_resample', 2),
               'Pb'    : col_descr('P_B', 2),
              }

    # Parity-violating signal (not for 'w')
    parity  = {'xi'  : col_descr('xi_x', 3),
               'gl'  : col_descr('g_x', 2),
               'Pb'  : col_descr('P_EB', 3),
              }

    error1  = {'Pb'  : col_descr('dP_E', 4),
              }

    error2  = {'Pb'  : col_descr('dP_B', 5),
              }

    error_parity = {'Pb': col_descr('dP_EB', 6)
              }

    # Scale boundaries
    scale_lo = {'Pb': col_descr('ell_lower', 7),
               }
              
    scale_hi = {'Pb': col_descr('ell_upper', 8),
               }
              
    # Resample rms
    rms1 = {'xi_res': col_descr('rms_p_resample', 3)}
    rms2 = {'xi_res': col_descr('rms_m_resample', 4)}

    xi_names       = ['theta', 'xi_p', 'xi_m']
    xi_indices     = [0, 1, 2]

    xi_names_opt   = ['xi_x']
    xi_indices_opt = [3]

    xi_res_names   = ['theta', 'xi_p_resample', 'xi_m_resample', 'rms_p_resample', 'rms_m_resample'] 
    xi_res_indices = [0, 1, 2, 3, 4]

    gl_names       = ['theta', 'g_t']
    gl_indices     = [0, 1]

    gl_names_opt   = ['g_x']
    gl_indices_opt = [2]

    w_names        = ['theta', 'w']
    w_indices      = [0, 1]

    Pb_names        = ['ell', 'P_E', 'P_B', 'P_EB']
    Pb_indices      = [0, 1, 2, 3]

    Pb_names_err   = ['dP_E', 'dP_B', 'dP_EB']
    Pb_indices_err = [4, 5, 6]



class xi_data:
    """Correlation function data
    
    Parameters
    ----------
    theta: array of double
        angular scales [rad]
    xip, xim, xix: array of double
        Two-point shear correlation function components
    file_format: string
        file format
    force_reg: bool, optional, default=False
        if True, forces angular binning to be on regular grid
    verbose: bool, optional, default=False
        verbose output if True
    """

    def __init__(self, theta, xip, xim, xix, file_format, force_reg=False, verbose=False):
        self.theta       = theta
        self.xip         = xip
        self.xim         = xim
        self.xix         = xix
        self.length      = len(theta)
        self.file_format = file_format

        self.set_binning(force_reg, verbose)
        if verbose == True:
            print('Binning type of 2PCF is \'{0}\''.format(self.binning))


    def set_binning(self, force_reg=False, verbose=False):
        """Sets irregular (recommended), linear, or logarithmic binning
        """

        if (self.length < 2): error('Data vector has length {0}, has to be larger than 2'.format(self.length))

        if force_reg == False:
            self.binning, self.delta = 'ireg', 0.0
            return

        eps     = 1.0e-2

        dtheta1 = self.theta[1] - self.theta[0]
        dtheta2 = self.theta[self.length-1] - self.theta[self.length-2]
        drel = abs(dtheta2 - dtheta1)/self.theta[1]
        if drel < eps:
            self.binning, self.delta = 'lin', dtheta2
            return

        dlogtheta1 = np.log10(self.theta[1]) - np.log10(self.theta[0])
        dlogtheta2 = np.log10(self.theta[self.length-1]) - np.log10(self.theta[self.length-2])
        drel = abs(dlogtheta2 - dlogtheta1)
        if  drel < eps:
            self.binning, self.delta = 'log', dlogtheta2
            return

        if verbose == True:
            print('Bins seem not to be regular, falling back to \'ireg\'')
        self.binning, self.delta = 'ireg', 0.0

    def get_binning(self):
        """ Returns the binning type and delta (0 if 'ireg')
        """

        if self.binning is None or self.delta is None:
            error('Binning not set')
        return self.binning, self.delta

    def dtheta_theta(self, i):
        """Returns dtheta_i * theta_i
        """

        if self.binning == 'lin':
            # dtheta * theta
            d = self.delta * self.theta[i]
        elif self.binning == 'log':
            # d log theta * theta^2 = d theta / theta * theta^2 = dtheta * theta
            d = self.delta * self.theta[i] * self.theta[i]
        elif self.binning == 'ireg':
            if i != self.length-1:
                d = (self.theta[i+1] - self.theta[i]) * self.theta[i]
            else:
                d = (self.theta[i] - self.theta[i-1]) * self.theta[i]
        else:
            mkstuff.error('Invalid binning type \'{}\''.format(self.binning))

        return d

    def gather_component(self, typ, append=False):
        """Return information on correlation function component.

        Parameters
        ----------

        typ: string
            'p', 'm', or 'x' for xi_p, xi_m, xi_x
        append: bool, optional, default=False
            if True, only append correlation function value, do not
            include scale

        Returns
        -------
        dat: 1- or 2-d array
            array of data columns
        """

        dat = []
        if not append:
            dat.append(self.theta)

        if typ == 'p':
            corr = self.xip
        elif typ == 'm':
            corr = self.xim
        elif typ =='x':
            corr = self.xix
        else:
            raise ValueError('Invalid correlation type \'{}\''.format(typ))

        dat.append(corr)

        return dat


    def write_component_ascii(self, output_path, typ):
        """Write one of the correlation function components to an ascii file.

        Parameters
        ----------
        output_path: string
            output file path
        typ: string
            'p', 'm', or 'x' for xi_p, xi_m, xi_x

        Returns
        -------
        None
        """

        dat = self.gather_component(typ, append=False)
        f = open(output_path, 'w')
        write_xi_data(f, dat)
        f.close()


class xi_data_tomo:
    """Tomographic correlation function data

    Parameters
    ----------
    xi: list of class xi_data
        Two-point shear correlation function for different redshift combinations
    verbose: bool, optional, default=False
        verbose output if True
    """

    def __init__(self, xi, verbose=False):

        # Number of redshift correlations
        self.nzcorr = len(xi)

        # Infer number of redshift bins, solve quadratic equation
        dnzbin = 0.5*(np.sqrt(8*self.nzcorr + 1) - 1)
        nzbin = int(dnzbin + 0.5)
        if dnzbin != nzbin:
            raise ValueError('Number of correlations nc={} has no integer solution nz for nc=nz(nz+1)/2'.format(self.nzcorr))

        self.nzbin = nzbin

        xi_all = []

        for xi_c in xi:
            xi_all.append(xi_c)

        self.xi_all = xi_all


    # TODO: check consistency between xis, e.g. theta.

    @classmethod
    def from_lists(self, theta, xip, xim, xix, file_format, force_reg=False, verbose=False):
        """Create tomographic correlation function from list of components.
           Probably not needed.

        Parameters
        ----------
        theta: array of double
            angular scales [rad]
        xip, xim, xix: list of arrays of double
            Two-point shear correlation function component for different redshift combinations
        file_format: string
            file format
        force_reg: bool, optional, default=False
            if True, forces angular binning to be on regular grid
        verbose: bool, optional, default=False
            verbose output if True
        """

        self.nzcorr = len(xip)
        self.xi_all = []

        for nz in range(self.nzcorr):
            xi = xi_data(theta, xip[nz], xim[nz], xix[nz], file_format, force_reg=force_reg, verbose=verbose)
            self.xi_all.append(xi)


    def get_xi(self, zc):
        """Return correlation function for a given redshift combination

        Parameters
        ----------
        zc: int
            redshift correlation index
        """

        return self.xi_all[zc]


    def get_theta(self):
        """Return angular scale array
        """

        return self.xi_all[0].theta


    def get_file_format(self):
        """Return file format string
        """

        return self.xi_all[0].file_format


    def write_component_ascii(self, output_path, typ):
        """Write one of the correlation function components to an ascii file.

        Parameters
        ----------
        output_path: string
            output file path
        typ: string
            'p', 'm', or 'x' for xi_p, xi_m, xi_x
        append: bool, optional, default=False
            if True, only append correlation function value, do not
            write scale or newline

        Returns
        -------
        None
        """

        xi = self.get_xi(0)
        dat = xi.gather_component(typ, append=False)

        for zc in range(1, self.nzcorr):
            xi = self.get_xi(zc)
            xi_comp = xi.gather_component(typ, append=True)
            dat.append(xi_comp[0])
        f = open(output_path, 'w')
        write_xi_data(f, dat)
        f.close()


def get_iz_jz(nzbin, zcorr):
    """Return redshift bins (iz, jz) corresponding to redshift correlation
       #zcorr
    """

    icorr = 0
    for iz in range(nzbin):
        for jz in range(iz, nzbin):
            if icorr == zcorr:
                return iz, jz
            icorr = icorr + 1

    raise ValueError('No redshift bins found for zcorr={}'.format(zcorr))



class pkappa_data:
    """Data for derived second-order functions: P_kappa as fct. of ell, the band-power spectrum P_i,
       and smooth real-space functions such as <M_ap^2>(theta).
       From pallas.py.
    """

    def __init__(self, ell_min, ell_max, band=False, error=False, Nell=None, unit_out=None):
        """Initialise and return a pkappa_data object. 

        Parameters
        ----------
        ell_min: float
            minimum ell mode or angular scale[rad]
        ell_max: float
            maximum ell mode or angular scale[rad]
        band: boolean
            If True, create band-power or angular bin limits on output (see *write*). Default is False.
        error: boolean
            If True, include errors, e.g. from resampling. Default is False.
        Nell: integer
            Number of logarithmic modes or bins. If None (default), Nell=(ell_max-ell_min) linear modes are created.
        unit_out: string
            Output unit for angular bins. If None (default), the standard unit 'arcmin' is used..

        Returns
        -------
        pkappa: class pkappa_data
            Data
        """

        if Nell is None:
            self.ell  = np.arange(ell_min, ell_max)
        else:
            self.ell  = np.logspace(np.log10(ell_min), np.log10(ell_max), Nell)

        self.pE   = np.zeros(np.shape(self.ell))
        self.pB   = np.zeros(np.shape(self.ell))
        self.pEB  = np.zeros(np.shape(self.ell))
        self.Nell = len(self.ell)
        self.band = band
        if unit_out is None:
            self.unit_out = 'none'
        else:
            self.unit_out = unit_out
        if error == True:
            self.dpE  = np.zeros(np.shape(self.ell))
            self.dpB  = np.zeros(np.shape(self.ell))
            self.dpEB = np.zeros(np.shape(self.ell))
        else:
            self.dpE  = None
            self.dpB  = None
            self.dpEB = None



    def write(self, fname, file_format, header=None, verbose=False):
        """Writes the power spectrum/smooth function to the ascii or fits
           file 'fname.[txt|fits]'
        """

        if file_format == 'ascii':
            ffname = fname + '.txt'
            self.write_ascii(ffname, header)
        elif file_format == 'fits':
            ffname = fname + '.fits'
            self.write_fits(ffname, header)
        else:
            error('Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(file_format))

        if verbose == True:
            print('Writing output file \'{0}\''.format(ffname))


    def write_ascii(self, fname, header=None):
        """Writes the power spectrum/smooth function to the ascii file 'fname'
        """

        f = open(fname, 'w')
        if header is not None:
            if self.unit_out is not 'none':
                my_header = header.split()
                my_header[1] = '{}[{}]'.format(my_header[1], self.unit_out)
                header = ' '.join(my_header)
            f.write(header + '\n')

        for j in range(self.Nell):
            scale = rad_to_unit(self.ell[j], self.unit_out)
            f.write('{0:10.3f} {1: 12.5e} {2: 12.5e} {3: 12.5e}'.format(scale, self.pE[j], self.pB[j], self.pEB[j]))
            if self.dpE is not None:
                f.write('{1: 12.5e} {2: 12.5e} {3: 12.5e}'.format(scale, self.dpE[j], self.dpB[j], self.dpEB[j]))
            if self.band == True:
                ell_l, ell_u = self.ell_l_u(j)
                ell_l = rad_to_unit(ell_l, self.unit_out)
                ell_u = rad_to_unit(ell_u, self.unit_out)
                f.write('{0:10.3f} {1:10.3f}'.format(ell_l, ell_u))
            f.write('\n')

        f.close()


    def write_fits(self, fname, header=None):
        """Writes the power spectrum/smooth function to the ascii file 'fname'
        """

        n_col = 4
        if self.dpE is not None:
            n_col += 3
        if self.band is True:
            n_col += 2
        if header is None:
            my_header = ['col_{}'.format(x) for x in range(n_col)]
        else:
            my_header = header.replace('#', '').split()

        scales = rad_to_unit(self.ell, self.unit_out)
        cols = []
        i = 0
        cols.append(fits.Column(name=my_header[i], format='E', array=scales))
        i += 1
        cols.append(fits.Column(name=my_header[i], format='E', array=self.pE))
        i += 1
        cols.append(fits.Column(name=my_header[i], format='E', array=self.pB))
        i += 1
        cols.append(fits.Column(name=my_header[i], format='E', array=self.pEB))

        if self.dpE is not None:
            i += 1
            cols.append(fits.Column(name=my_header[i], format='E', array=self.dpE))
            i += 1
            cols.append(fits.Column(name=my_header[i], format='E', array=self.dpB))
            i += 1
            cols.append(fits.Column(name=my_header[i], format='E', array=self.dpEB))

        if self.band is True:
            ell_l = []
            ell_u = []
            for j in range(self.Nell):
                my_ell_l, my_ell_u = self.ell_l_u(j)
                my_ell_l = rad_to_unit(my_ell_l, self.unit_out)
                my_ell_u = rad_to_unit(my_ell_u, self.unit_out)
                ell_l.append(my_ell_l)
                ell_u.append(my_ell_u)

            i += 1
            cols.append(fits.Column(name=my_header[i], format='E', array=ell_l))
            i += 1
            cols.append(fits.Column(name=my_header[i], format='E', array=ell_u))

        fcols = fits.ColDefs(cols)

        if sys.version_info >=  (3,3):
            hdu = fits.BinTableHDU.from_columns(fcols)
        else:
            hdu = fits.new_table(fcols)
        if self.unit_out is not 'none':
            hdu.header.append(card=('UNITS', self.unit_out, 'Coordinate units'))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            hdu.writeto(fname, clobber=True)


    def ell_l_u(self, j):
        """Returns the lower and upper band ell
        """

        dlogell = log(self.ell[1]) - log(self.ell[0])
        return self.ell[j] / sqrt(exp(dlogell)), self.ell[j] * sqrt(exp(dlogell))



def unit_to_rad(theta, unit):
    return theta * units[unit]



def rad_to_unit(theta, unit):
    return theta / units[unit]



def get_unit(field, verbose=False, exit=False):
    """Returns unit for angular scales, from header line with format e.g. 'theta[UNIT]'.
       If exit=True, exists with error. If exit=False, returns unit_default.
    """

    for u in units.keys():
        if u in field:
            if verbose == True:
                print('Unit of input angular scales = {0}'.format(u))
            return u

    if exit == True:
        mkstuff.error('No unit found, exiting')

    if verbose == True:
        warnings.warn('No unit for angular scales found, assuming \'{0}\''.format(unit_default))

    return unit_default



def get_unit_fits(header, verbose=False, exit=False):

    if 'UNITS' in header.keys():
        unit = header['UNITS']
    else:
        if exit == True:
            mkstuff.error('No unit found, exiting')
        if verbose == True:
            mkstuff.warning('No \'UNITS\' keyword in fits header found, assuming default unit for angular scales \'{0}\''.format(unit_default))
        return unit_default

    # Go through valid units and compare to assigned one (from header)
    for u in units.keys():
        if u == unit:
            if verbose == True:
                print('Unit of input angular scales = {0}'.format(u))
            return u

    if verbose == True:
        mkstuff.warning('Unit from fits header (\'UNITS\' = \'{0}\') not a valid unit, reverting to default unit \'{1}\''.format(unit, unit_default))
        return unit_default



def get_athena_path(version='C'):
    """Return path to athena.

    Parameters
    ----------
    version: string
        One of ['C' (default)]

    Returns
    -------
    path: string
        /path/to/athena
    """

    if version == 'C':
        return '{0}/athena_git_Euclid/bin'.format(os.environ['HOME'])
    else:
        mkstuff.error('Wrong athena verions {}'.format(version))



# Writes 'athena' config file to path
def write_config_tree(path='config_tree', galcat1=None, galcat2=None, sformat='standard', ncol=None, col_names=None,
	        	      wcorr=1, scoord_input='arcmin', scoord_output=None, thmin=0.05, thmax=120, nth=20,
                      bintype='LOG', radec=0, oath=0.03, serror=None, nresample='0 0'):

    if galcat2 is None:
        galcat2 = '-'

    if scoord_output is None:
        scoord_output = scoord_input

    if serror is None:
        serror    = 'none'
        nresample = '0 0'

    if type(nresample) is int:
        nresample = '{} {}'.format(nresample, nresample)

    conf = open(path, 'w')

    print('GALCAT1        = ' + galcat1, file=conf)
    print('GALCAT2        = ' + galcat2, file=conf)
    print('WCORR          = ' + str(wcorr), file=conf)
    if wcorr == 2:
        print('SWCORR_SUBTYPE  = nn_2d', file=conf)
    print('SFORMAT        = ' + sformat, file=conf)
    if sformat == 'fits':
    	print('NCOL           = {}'.format(ncol), file=conf)
    	print('COL_NAMES      = {}'.format(col_names), file=conf)
    print('SCOORD_INPUT   = ' + str(scoord_input), file=conf)
    print('SCOORD_OUTPUT  = ' + str(scoord_output), file=conf)
    print('THMIN          = ' + str(thmin), file=conf)
    print('THMAX          = ' + str(thmax), file=conf)
    print('NTH            = ' + str(nth), file=conf)
    print('BINTYPE        = ' + bintype, file=conf)
    print('RADEC          = ' + str(radec), file=conf)
    print('OATH           = ' + str(oath), file=conf)
    print('SERROR         = ' + str(serror), file=conf)
    print('NRESAMPLE      = ' + nresample, file=conf)

    conf.close()



def run_athena_list(config_names, suffixes=None, n_cpu=1, version='C', verbose=False, run=True):
    """Run athena.

    Parameters
    ----------
    config_names: list of strings
        list of config names
    suffixes: list of strings
	list of suffix names, default empty, uses config_names instead
    n_cpu: integer
        number of CPUs, default=1
    version: string
        One of ['C' (default)]
    verbose: boolean
        verbose output, default False
    run: boolean
	if False, do not run athena. Default = True
    
    Returns
    -------
    None
    """

    athena_path = get_athena_path(version=version)
    if verbose:
        q_ath = ''
    else:
        q_ath = ' -q'


    cmds = []
    for i, cn in enumerate(config_names):
        if os.path.isfile(cn):
            if suffixes is None or len(suffixes) == 0:
                suff = cn
            else:
                suff = suffixes[i]
            cmds.append('{}/athena -c {} --out_suf {}{}'.format(athena_path, cn, suff, q_ath))

    if verbose is True:
        print('Running {} athena jobs on {} CPUs'.format(len(cmds), n_cpu))

    res = mkstuff.run_cmd_list(cmds, ncpu=n_cpu, verbose=verbose, run=run)
    if res != 0:
        mkstuff.error('At least one athena run returned error code {}'.format(res))



def check_fields(header, strings, indices, exit, verbose=False):
    """Prints error (warning) message if header (array of strings) does not contain necessary (optional)
       strings.
    """

    if header is None:
        if verbose is True:
            mkstuff.warning('No header found, using default columns')
        if exit is True:
            return False

    for i in range(len(strings)):
        if len(header) <= indices[i] or not strings[i] in header[indices[i]]:
            msg = 'Field \'{0}\' not found in header (expected in field #{1})'.format(strings[i], indices[i])
            if exit is True:
                mkstuff.error(msg)
            else:
                if verbose is True:
                    mkstuff.warning('{}, continuing'.format(msg))
                return False

    return True



def get_header_columns_from_fits_table(hdu):
    """Checks hdu for table format and returns header column values.

    Parameters
    ----------
    hdu: hdu
        FITS hdu

    Returns
    -------
    header_col_vals: array of strings
        header column names, or None if hdu not FITS table
    """

    hdu_type = type(hdu)
    if hdu_type == fits.hdu.table.BinTableHDU or hdu_type == fits.hdu.table.TableHDU:
        header_all_keys = hdu.header.keys()
        header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]
        header_col_vals = [hdu.header.get(s) for s in header_col_keys]
        return header_col_vals
    else:
        return None


def get_file_format(file_format, fname):
    """Return file_format if not None, or 'fits' if extension is '.fits' or '.FITS', and 'ascii' otherwise.
    """

    if file_format is None:
        # Determine format using file extension
        extension = os.path.splitext(fname)[1]
        if extension == '.fits' or extension == '.FITS':
            file_format = 'fits'
        else:
            # Default
            file_format = 'ascii'

    elif file_format != 'fits' and file_format != 'ascii':
       print(file_format)
       error('get_file_format: Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(file_format))

    return file_format


def read_xi(fname, file_format, force_reg, verbose=False):
    """Reads a shear-shear correlation function file (ascii or fits).
    """

    my_file_format = get_file_format(file_format, fname)

    if verbose:
        print("File format is {0}".format(my_file_format))

    if my_file_format == 'fits':
        theta, xip, xim, xix = read_xi_fits(fname, force_reg, verbose)
    elif my_file_format == 'ascii':
        theta, xip, xim, xix = read_xi_ascii(fname, force_reg, verbose)
    else:
       mkstuff.error('read_xi: Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(my_file_format))

    xi = xi_data(theta, xip, xim, xix, my_file_format, force_reg=force_reg, verbose=verbose)

    return xi



def read_xi_ascii(fname, force_reg, verbose=False):
    """Read a shear-shear correlation function ascii file.
    """

    if verbose is True: print('Reading input xi file \'{0}\''.format(fname))

    xi, ncomment, header  = mkstuff.read_table(fname, True)

    if xi.shape[1] <= 2:
        mkstuff.error('Input file \'{0}\' has only {1} columns. Minimum columns are (theta, xi+, xi-)'.format(fname, xi.shape[1]))

    if ncomment > 0:
        fields = header[0].replace('#', '').split()
    else:
        fields = None

    # Required fields
    check_fields(fields, in_cols.xi_names, in_cols.xi_indices, True, verbose)

    theta = xi[:, in_cols.scale.num]
    xip   = xi[:, in_cols.signal1['xi'].num]
    xim   = xi[:, in_cols.signal2['xi'].num]

    # Optional fields
    has_xix = check_fields(fields, in_cols.xi_names_opt, in_cols.xi_indices_opt, False, verbose)
    if has_xix == True:
        xix  = xi[:, in_cols.parity['xi'].num]
    else:
        xix  = None

    # Units of angular scales
    if fields is not None:
        unit = get_unit(fields[in_cols.scale.num], verbose)
    else:
        unit = unit_default
    theta  = unit_to_rad(theta, unit)

    return theta, xip, xim, xix


def read_xi_fits(fname, force_reg, verbose=False):
    """Read a shear-shear correlation function file in fits format.
    """

    hdulist = fits.open(fname)
    for ihdu, hdu in enumerate(hdulist):

        header_col_vals = get_header_columns_from_fits_table(hdu)
        if header_col_vals is not None:

            # Required fields
            check_fields(header_col_vals, in_cols.xi_names, in_cols.xi_indices, True, verbose)
            theta = hdu.data[in_cols.scale.name]
            xip   = hdu.data[in_cols.signal1['xi'].name]
            xim   = hdu.data[in_cols.signal2['xi'].name]

            # Optional fields
            has_xix = check_fields(header_col_vals, in_cols.xi_names_opt, in_cols.xi_indices_opt, False, verbose)
            if has_xix == True:
                xix = hdu.data[in_cols.parity['xi'].name]
            else:
                xix = None

            # Units of angular scales
            unit  = get_unit_fits(hdu.header, verbose)
            theta = unit_to_rad(theta, unit)

            break


    hdulist.close()

    return theta, xip, xim, xix


def read_xi_rms_resample_fits(fname, force_reg, verbose=False, stop=True):
    """Read a shear-shear correlation function resample (containing mean and rms) file in fits format.

    Parameters
    ----------
    fname: string
        input file name
    force_reg: boolean
        if True, forces regular binning
    verbose: boolean
        Verbose output, default=False
    stop: boolean
        If True (default), exit with error. If False, return None

    Returns
    -------
    theta: array of double, length N_theta
        angular scales
    xi: array of double, length N_resample x N_theta
        resampled correlation values
    resample_type: string
        resample type, 'bootstrap', 'jackknife'
    """

    hdulist = fits.open(fname)
    nhdu    = len(hdulist)
    theta   = None
    for hdu in range(nhdu):
        hdu_type = type(hdulist[hdu])

        # TODO: Use get_header_columns_from_fits_table
        if hdu_type == fits.hdu.table.BinTableHDU or hdu_type == fits.hdu.table.TableHDU:

            if 'OBJECT' in hdulist[hdu].header and hdulist[hdu].header.get('OBJECT') == 'xi.resample':
                if verbose == True:
                    print('Table found in fits file {0} (hdu #{1})'.format(input, hdu))
                header_all_keys = hdulist[hdu].header.keys()
                header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]
                header_col_vals = [hdulist[hdu].header.get(s) for s in header_col_keys]

                # Required fields
                check_fields(header_col_vals, in_cols.xi_res_names, in_cols.xi_res_indices, True, verbose)
                theta     = hdulist[hdu].data[in_cols.scale.name]
                xi_p_res  = hdulist[hdu].data[in_cols.signal1['xi_res'].name]
                xi_m_res  = hdulist[hdu].data[in_cols.signal2['xi_res'].name]
                rms_p_res = hdulist[hdu].data[in_cols.rms1['xi_res'].name]
                rms_m_res = hdulist[hdu].data[in_cols.rms2['xi_res'].name]

                # Units of angular scales
                unit  = get_unit_fits(hdulist[hdu].header, verbose)
                theta = unit_to_rad(theta, unit)

                break


    hdulist.close()

    if theta is None:
        msg = 'No hdu with resample information (OBJECT \'xi.resample\') found in file {}'.format(fname)
        if stop is True:
            mkstuff.error(msg)
        else:
            if verbose is True:
                print('{}, continuing'.format(msg))
            return None

    return theta, xi_p_res, xi_m_res, rms_p_res, rms_m_res



def read_xi_resample_fits(fname, force_reg, check_OBJECT=None, verbose=False):
    """Read a shear-shear correlation function file with individual resample results, in fits format (xi_p or xi_m).

    Parameters
    ----------
    fname: string
        input file name
    force_reg: boolean
        if True, forces regular binning
    check_OBJECT: string
        If not None, checks whether the given string equals the value of the OBJECT key in
        fits header, and returns with error if not. Default: None
    verbose: boolean
        Verbose output, default=False

    Returns
    -------
    theta: array of double, length N_theta
        angular scales
    xi: array of double, length N_resample x N_theta
        resampled correlation values
    resample_type: string
        resample type, 'bootstrap', 'jackknife'
    """

    hdulist = fits.open(fname)
    nhdu    = len(hdulist)
    for hdu in range(nhdu):
        hdu_type = type(hdulist[hdu])
    
        # TODO: Use get_header_columns_from_fits_table
        if hdu_type == fits.hdu.table.BinTableHDU or hdu_type == fits.hdu.table.TableHDU:
            if verbose == True:
                print('Table found in fits file {0} (hdu #{1})'.format(input, hdu))
            header_all_keys = hdulist[hdu].header.keys()
            header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]
            header_col_vals = [hdulist[hdu].header.get(s) for s in header_col_keys]

            if check_OBJECT is not None:
                if 'OBJECT' in hdulist[hdu].header:
                    val = hdulist[hdu].header.get('OBJECT')
                    if val != check_OBJECT:
                        mkstuff.error('Found in fits header \'OBJECT = {}\', not equal to required \'\''.format(val, check_OBJECT))
                else:
                        mkstuff.error('Key \'OBJECT\' not found in header')

            if 'TYPE' in hdulist[hdu].header:
                resample_type = hdulist[hdu].header.get('TYPE')
                if resample_type != 'bootstrap' and resample_type != 'jackknife':
                    mkstuff.error('Unknown resample type \'{}\' found in fits file (keyword \'TYPE\')'.format(resample_type))

            check_fields(header_col_vals, [in_cols.xi_names[0]], [in_cols.xi_indices[0]], True, verbose)

            # Angular scales
            theta  = hdulist[hdu].data[in_cols.scale.name]
            unit  = get_unit_fits(hdulist[hdu].header, verbose)
            theta = unit_to_rad(theta, unit)

            # Correlation data
            xi = []
            for i in range(1, len(header_col_vals)):
                key = header_col_vals[i]
                xi_this_res = hdulist[hdu].data[key]
                xi.append(xi_this_res)

            break

    return theta, xi, resample_type



def read_pkappa_band_fits(fname, verbose=False, stop=False):
    """Read from file the band-power spectrum P_i.

    Parameters
    ----------
    fname: string
        file name
    verbose: boolean
        verbose mode if True, default False
    stop: boolean
        If True (default), exit with error. If False, return None

    Returns
    -------
    pb: class pkappa_data
        data
    """

    try:
        hdulist = fits.open(fname)
    except IOError:
        msg = 'File \'{}\' does not exist'.format(fname)
        mkstuff.error(msg, stop=stop, verbose=verbose)
        if stop is False:
            return None

    for ihdu, hdu in enumerate(hdulist):

        header_col_vals = get_header_columns_from_fits_table(hdu)
        if header_col_vals is not None:
            check_fields(header_col_vals, in_cols.Pb_names, in_cols.Pb_indices, True, verbose)

            if verbose is True:
                print('Table with band-power spectrum data found at hdu #{}'.format(ihdu))

            ell  = hdu.data[in_cols.scale_F.name]
            P_E  = hdu.data[in_cols.signal1['Pb'].name]
            P_B  = hdu.data[in_cols.signal2['Pb'].name]
            P_EB = hdu.data[in_cols.parity['Pb'].name]
            ell_lower = hdu.data[in_cols.scale_lo['Pb'].name]
            ell_upper = hdu.data[in_cols.scale_hi['Pb'].name]

            has_errors = check_fields(header_col_vals, in_cols.Pb_names_err, in_cols.Pb_indices_err, False, verbose)
            if has_errors is True:
                dP_E  = hdu.data[in_cols.error1['Pb'].name]
                dP_B  = hdu.data[in_cols.error2['Pb'].name]
                dP_EB = hdu.data[in_cols.error_parity['Pb'].name]
            else:
                dP_E  = None
                dP_B  = None
                dP_EB = None

            break

    hdulist.close()

    if ell is None:
        pb = None
    else:
        pb = pkappa_data(ell[0], ell[1], band=True, error=True, Nell=len(ell))
        pb.ell  = ell
        pb.pE   = P_E
        pb.pB   = P_B
        pb.pEB  = P_EB
        pb.dpE  = dP_E
        pb.dpB  = dP_B
        pb.dpEB = dP_EB

    return pb


def write_xi_data(f, dat):

    scale = rad_to_unit(dat[0], unit_default)
    offset = 1
    nzcorr = len(dat) - offset
    for j in range(len(scale)):
        f.write('{:15.8f}'.format(scale[j]))
        for zc in range(offset, nzcorr + offset):
            f.write(' {:15.8e}'.format(dat[zc][j]))
        f.write('\n')


################
### Plotting ###
################

def prepare_plot(theta, title, xlabel, ylabel, ylog=False):
    """Common routine for plot preparation
    """

    f_err = np.log10(1.02)
    f_lim = 1.1

    plt.clf()

    log_x = np.log10(theta)
    mkplot.make_log_ticks('x', min(log_x), max(log_x))

    #plt.xlim(theta[0] / f_lim, theta[-1] * f_lim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    return f_err, f_lim, log_x



def get_theta_xlabel(scale_mode, wgl, unit_force=None):
    """Return angular scale column and xlabel string, from file wgl 

    Parameters
    ----------
    scale_mode: string
        'physical' or 'angular'
    wgl: Class table
        wgl file content, e.g. obtained from ascii.read
    unit_force: string, default=None
        If not None, use unit_force for label

    Returns
    -------
    theta: array of double
        content of scale column
    xlabel: string
        x-axis label
    """

    if scale_mode == "angular":
        x        = ['theta']
        for key in units:
            x.append('theta[{}]'.format(key))
    elif scale_mode == "physical":
        x        = ['R[Mpc]', 'R']

    theta = None
    for myx in x:
        if myx in wgl.keys():
            theta  = wgl[myx]
            xlabel = myx
            break

    if unit_force is not None:
        pattern = re.compile('\[.*\]')
        xlabel = pattern.sub('[{}]'.format(unit_force), xlabel)

    if theta is None:
        mkstuff.error('Column for angular scale not found')

    return theta, xlabel



def get_col_label(subtype, unit=''):
    """Return column names and plot labels.

    Parameters
    ----------
    subtype: string
        'wgl', 'wgl_res', 'DeltaSigma'
    unit: string
        '', 'M_sun/pc2'

    Returns
    -------
    None
    col: array of strings
        column names
    label: array of strings
        labels
    ylabel: string
        y-axis label
    """

    if unit == '' and subtype == 'DeltaSigma':
        unit = 'M_sun/pc2'

    if unit != '':
        unit_str   = '[{}]'.format(unit)
        unit_label =  '[{}]'.format(u.Unit(unit).to_string('latex_inline'))

    if subtype == 'wgl':
        col    = ['g_t', 'g_x', 'sqrt_Dcor', 'sqrt_Dcor']
        label  = ['$\\gamma_{\\rm t}$', '$\\gamma_\\times$']
        ylabel = 'Tangential and cross-shear'
    elif subtype == 'wgl_res':
        col    = ['g_t_resample', 'g_t_rms_resample']
        label  = ['$\\gamma_{\\rm t}$']
        ylabel = 'Tangential shear'
    elif subtype == 'DeltaSigma':
        col    = ['Delta_Sigma_t{}'.format(unit_str), 'Delta_Sigma_x{}'.format(unit_str), 'sqrt_Dcor', 'sqrt_Dcor']
        label  = ['$\\Delta \\Sigma_{{\\rm t}}$ {}'.format(unit_label), '$\\Delta \\Sigma_\\times$ {}'.format(unit_label)]
        ylabel = 'Projected overdensity'
    else:
        mkstuff.error('Unknown subtype \'{}\''.format(subtype))

    return col, label, ylabel


def signal_to_noise_w2_var(y, dy):
    """Return signal to noise for w2 type correlation using the variance.
    """

    snr = np.sqrt(sum(y**2 / dy**2))

    return snr


def get_index_range(theta, scales=''):
    """Return min and max indices for *theta* that are within the range specified by *scales*.

    Parameters
    ----------
    theta: array of float
        x-axis data (scales)
    scales_SNR: string with separator '_' or ' '.
        min and max scales. If None (default) return full range.

    Returns
    -------
    imin, imax: integers
        Min and Max index.
    """

    if scales is None:
        imin, imax = 0, len(theta)-1
    else:
        smin, smax = [float(i) for i in mkstuff.my_string_split(scales, num=2, stop=True)]

        # Index of value in theta to the right of smin. Can be 0.
        imin = bisect_left(theta, smin)
        # Index of value in theta to the left of xmin. Can be len(theta).
        imax = bisect_right(theta, smax)-1

    return imin, imax



def plot_corr_w2(in_name='wgl', out_name=None, scale_mode='angular', subtype='wgl', SNR=False, scales_SNR=None,
                 unit='', with0=False, unit_force=None, ylog=False):
    """Plot correlation of type which=2 (shear-position).

    Parameters
    ----------
    in_name: string
        input file name
    out_name: string
        input file name, with '.pdf' added. If None (default), out_name = in_name
    scale_mode: string
        'angular' (default) or 'physical'
    subtype: string
        'wgl', 'wgl_res', 'DeltaSigma'
    SNR: boolean
        If True, calculates and adds as text to plot signa-to-noise ratio.
    scales_SNR: string, two numbers with separator '_' or ' '
        Minimum and maximum scale for SNR calculation. If None (default) use full range.
    unit: string
        '', 'M_sun/pc2'
    with0: bool
        If True, adds the y=0 line
    unit_force: string, default=None
        If not None, use unit_force for label
    ylog: boolean
        If True, logarithmic y-axis. Default is False.

    Returns
    -------
    None
    """

    col, label, ylabel = get_col_label(subtype, unit=unit)

    plot_corr_w2_internal(in_name, col, label, ylabel, scale_mode=scale_mode, SNR=SNR, scales_SNR=scales_SNR, with0=with0, unit_force=unit_force, ylog=ylog)

    if out_name is None:
        out_name = in_name
    plt.savefig('{}.pdf'.format(out_name))




def plot_corr_w2_internal(in_path, col, label, ylabel, scale_mode='angular', SNR=False, scales_SNR=None,
                          with0=False, unit_force=None, ylog=False):
    """Plot correlation of type which=2 (shear-position), internal version. Called by plot_corr_w2.

    Parameters
    ----------
    in_path: string
        input file path
    col: array of string
        column labels
    label: array of strings
        labels
    ylabel: string
        y-axis label
    scale_mode: string
        'angular' (default) or 'physical'
    SNR: boolean
        If True, calculates and adds as text to plot signa-to-noise ratio.
    scales_SNR: string, two numbers with separator '_' or ' '
        Minimum and maximum scale for SNR calculation. If None (default) use full range.
    with0: bool
        If True, adds the y=0 line
    unit_force: string, default=None
        If not None, use unit_force for label
    ylog: boolean
        If True, logarithmic y-axis. Default is False.

    Returns
    -------
    None
    """

    w2 = ascii.read('{}'.format(in_path))

    theta, xlabel = get_theta_xlabel(scale_mode, w2, unit_force=unit_force)

    ww = []
    for c in col:
        ww.append(w2[c])

    plot_corr_w2_col(theta, ww, xlabel, ylabel, label, in_path, SNR=SNR, scales_SNR=scales_SNR, with0=with0, ylog=ylog)



def plot_corr_w2_col(theta, w2, xlabel, ylabel, label, title='', SNR=False, scales_SNR=None, with0=False, ylog=False):
    """Plot columns theta versus w2, containing data and error column(s). Called by plot_corr_w2_internal.

    Parameters
    ----------
    theta: array of length N of float
        x-axis data
    w2: array NxNcol of float
        y-axis data, sets of data and error columns
    xlabel: string
        x-axis label
    ylabel: string
        y-axis label
    label: array of length 2 of string
        labels for the two data sets (w2)
    title: string
        title string (default: empty)
    SNR: boolean
        If True, calculates and adds as text to plot signa-to-noise ratio.
    scales_SNR: string, two numbers with separator '_' or ' '
        Minimum and maximum scale for SNR calculation. If None (default) use full range.
    with0: bool
        If True, adds the y=0 line. Default: False
    ylog: bool
        If True, plot log of y-axis, default: False

    Returns
    -------
    None
    """

    color = ['b', 'g']
    msize = [5, 4]
    sign  = [-1, 1]
    f_err, f_lim, log_x = prepare_plot(theta, title, xlabel, ylabel)

    if len(w2) == 4:
        yerr = [w2[2], w2[3]]
        Ndata = 2
    elif len(w2) == 2:
        yerr = [w2[1]]
        Ndata = 1
    else:
        mkstuff.error('Wrong number {} of columns in w2'.format(len(w2)))

    if ylog is True:

        for i in range(0, Ndata):

            # dlogy = dy/y
            yerr[i] = yerr[i] / w2[i]

            id_pos = np.where(w2[i]>0)
            id_neg = np.where(w2[i]<=0)
            y_pos  = np.log10(w2[i][id_pos])
            y_neg  = np.log10(-w2[i][id_neg])
            if i == 0:
                mkplot.make_log_ticks('y', min(y_pos)-2, max(y_pos)+1)

            plt.plot(log_x[id_pos] + sign[i]*f_err, y_pos, '{}.'.format(color[i]), marker='o', markersize=msize[i], label=label[i])
            plt.plot(log_x[id_neg] + sign[i]*f_err, y_neg, '{}o'.format(color[i]), marker='o', mfc='none', markersize=msize[i], label='')
            plt.errorbar(log_x[id_pos] + sign[i]*f_err, y_pos, yerr=yerr[i][id_pos], fmt='none',
                ecolor='{}'.format(color[i]), elinewidth=1, label='')
            plt.errorbar(log_x[id_neg] + sign[i]*f_err, y_neg, yerr=yerr[i][id_neg], fmt='none',
                ecolor='{}'.format(color[i]), elinewidth=0.5, label='')

    else:

        for i in [0, Ndata-1]:
            plt.plot(log_x + sign[i]*f_err, w2[i], '{}'.format(color[i]), linewidth=1, label=label[i])
            plt.errorbar(log_x + sign[i]*f_err, w2[i], yerr=yerr[i], fmt='none', ecolor='{}'.format(color[i]), label='')

        if with0 is True:
            plt.plot([min(log_x), max(log_x)], [0, 0], 'k', linewidth=1, label=None)

    plt.legend(loc='upper right', numpoints=1, frameon=False)

    if SNR is True:

        imin, imax = get_index_range(theta, scales_SNR)

        ax = plt.subplot(1, 1, 1)
        if len(w2) == 4:
            snr   = signal_to_noise_w2_var(w2[0][imin:imax+1], w2[2][imin:imax+1])
            snr_x = signal_to_noise_w2_var(w2[1][imin:imax+1], w2[3][imin:imax+1])
            plt.text(0.95, 0.7, 'SNR_x = {}'.format('%.2f' % snr_x), ha='right', transform=ax.transAxes)
        elif len(w2) == 2:
            snr   = signal_to_noise_w2_var(w2[0][imin:imax+1], w2[1][imin:imax+1])
        plt.text(0.95, 0.75, 'SNR_t = {}'.format('%.2f' % snr), ha='right', transform=ax.transAxes)
        ymin, ymax = plt.ylim()
        dx2 = (log_x[1] - log_x[0])/2
        plt.plot([log_x[imin]-dx2, log_x[imin]-dx2], [ymin, ymax], linestyle='dashed', color='black')
        plt.plot([log_x[imax]+dx2, log_x[imax]+dx2], [ymin, ymax], linestyle='dashed', color='black')

    # From plot_wgl.py:plot_stacked_wgl
    mkplot.update_lim('x', max_new = np.log10(max(theta) * f_lim))

    
 
