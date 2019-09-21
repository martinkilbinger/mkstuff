#!/usr/bin/env python

"""

:Name: ascii2fits.py

:Description: Transforms an ascii (table) file into FITS.

:Author: Martin Kilbinger

:Date: 2014

:Package: mkstuff

"""

# Compability with python2.x for x>6
from __future__ import print_function


import sys
from optparse import OptionParser

from astropy.io import ascii
from astropy.table import Table, Column
import random

import mkstuff

try:
    import pyfits as fits
except ImportError:
    pass
    try:
        from astropy.io import fits
    except ImportError:
        error("Could not import pyfits or astropy.io/fits library")


import mkstuff



def read_ascii_file(name, no_header, verbose=False):
    """Read ascii file and return data, header.
    """

    if no_header is True:
        data = ascii.read(name, format='no_header', delimiter='\s')
    else:
        data = ascii.read(name)
        if verbose:
            print(data.keys())

    header = data.keys()

    if verbose:
        print('Header')
        print(header)

    return data, header



def get_header(in_header, user_header):
    """Return either data header or user-defined header given on command line
    """

    if user_header is None:
        header = in_header
    else:
        header = user_header.split()

    # Are header and data dimensions equal?
    ncol = len(in_header)
    if ncol != len(header):
        mkstuff.error('Data has {0} columns, header on command line {1}.'.format(ncol, len(header)))

    return header



def every(data, every):
    """Return a random sample of 1/every objects
    """

    n      = len(data)
    sample = random.sample(range(0, n), n/every)
    data_ev = data[sample]

    return data_ev



def write_fits_file(data, in_header, out_header, name, verbose=False):
    """Write fits file.
    """

    ncol = len(out_header)

    if verbose:
        print('Guessing column type info, from first-row entries:')

    cols  = []
    for i in range(ncol):

        entry = data[in_header[i]][0]
        if isinstance(entry, float):
            form = 'E'
        elif isinstance(entry, int):
            form = 'I'
        else:
            form = '100A'

        if verbose:
            print('  col {}: {} -> \'{}\''.format(i, entry, form, end=''))

        new_col = fits.Column(name = out_header[i], format=form, array=data[in_header[i]])
        cols.append(new_col)

    coldefs = fits.ColDefs(cols)
    hdu     = fits.BinTableHDU.from_columns(coldefs)
    if verbose:
        print('Writing FITS file \'{}\''.format(name))
    hdu.writeto(name, clobber=True)



def main(argv=None):
    """Main program
    """

    # Command line options
    usage  = "%prog [OPTIONS]"

    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--input', dest='input', type='string', help='input xi FITS file name')
    parser.add_option('-o', '--output', dest='output', type='string', default=None,
                    help='output base name (default = <INPUT>.fits)')
    parser.add_option('-e', '--every', dest='every', type='int', default=1,
                    help='Write only a fraction of 1/EVERY objects (default = 1')
    parser.add_option('-n', '--no_header', dest='no_header', action='store_true', default=False,
                    help='No header in input file')
    parser.add_option('-H', '--header', dest='header', type='string',
                    help='Header string \'HEADER\' (white-spaced list; default: input file header)')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='Verbose')


    options, args = parser.parse_args()

    see_help = 'See option \'-h\' for help.'

    if options.input is None:
        print('Input ascii file not given (use option \'-i\'). ' + see_help, file=sys.stderr)
        return

    if options.output is None:
        options.output = '{0}.fits'.format(options.input)


    if options.verbose == True:
        print('Reading input ascii file \'{}\''.format(options.input))
    in_ascii, in_header = read_ascii_file(options.input, options.no_header)

    if options.every > 1:
        if options.verbose == True:
            print('Selection randomly 1 in {} objects'.format(options.every))
        in_ascii = every(in_ascii, options.every)

    out_header = get_header(in_header, options.header)

    if options.verbose == True:
        print('Writing ascii file {}'.format(options.output))

    write_fits_file(in_ascii, in_header, out_header, options.output, verbose=options.verbose)


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))


