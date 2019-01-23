#!/usr/bin/env python

  
"""

:Name: fits2ascii.py

:Description: Prints content of FITS file to ascii.

:Author: Martin Kilbinger

:Date: 2014

:Version: 1.0 (2014)

:Package: mkstuff

"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys

from astropy.io import fits

from optparse import OptionParser



def read_and_print_fits_file(input, output=None, verbose=False, n_row=-1, start=0, every=-1, out_hdu=-1, header=2, line_count=False):
    """Reads fits file and prints content to file/stdout
    """

    keys = ['OBJECT', 'TYPE', 'NRESAMPLE', 'UNITS', 'EXTNAME']
    

    hdulist = fits.open(input, mode='denywrite', memmap=True, do_not_scale_image_data=True)
    nhdu    = len(hdulist)

    for hdu in range(nhdu):

        if out_hdu != -1 and out_hdu != hdu:
            continue

        hdu_type = type(hdulist[hdu])
        if hdu_type == fits.hdu.table.BinTableHDU or hdu_type == fits.hdu.table.TableHDU:

            if output is None:
                fout = sys.stdout
            else:
                if out_hdu == -1:
        	        fout = open('{}_{}.txt'.format(output, hdu), 'w')
                else:
        	        fout = open('{}.txt'.format(output), 'w')
 
            if verbose:
                print('Table found in fits file {0} (hdu #{1})'.format(input, hdu), file=sys.stderr)
            header_all_keys = hdulist[hdu].header.keys()
            header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]

            # Print some header keys 
            if header == 2:
                for key in keys:
                    if key in hdulist[hdu].header:
                        print('# {0} = {1} ({2})'.format(key, hdulist[hdu].header.get(key), hdulist[hdu].header.comments[key]), file=fout)

            # Print header
            if header != 0:
                print('#', end=' ', file=fout)
                for key in header_col_keys:
                    header_col_val = hdulist[hdu].header.get(key)
                    print('{0:>18}'.format(header_col_val), end='', file=fout)
                print('', file=fout)

            # Print data
            nrows = hdulist[hdu].header.get('NAXIS2')
            if line_count == True:
                print('# number of lines = {} = {:.2e}'.format(nrows, nrows))
                return

            if n_row > 0:
                n_max = min(n_row, nrows)
            else:
                n_max = nrows

            for n in range(start, start+n_max):

                if every!=-1 and n % every != 0:
                    continue

                print('   ', end='', file=fout)
                for k, key in enumerate(header_col_keys):
                    header_col_val = hdulist[hdu].header.get(key)
                    val   = hdulist[hdu].data[header_col_val][n]

                    # Print according to format string in header
                    tform = 'TFORM{}'.format(k+1)
                    code  = hdulist[hdu].header[tform]
                    if code == 'E' or code == 'D':
                        print('{0: 17.10g} '.format(val), end='', file=fout)
                    elif code == 'I' or code == 'K':
                        print('{0:18d} '.format(val), end='', file=fout)
                    else:
                        print('{0:7s} '.format(str(val)), end='', file=fout)
                print('', file=fout)

            if output is not None:
                fout.close()

    hdulist.close()



def main(argv=None):
    """Main program
    """

    # Command line options
    usage  = "%prog [OPTIONS]"

    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--input', dest='input', type='string', help='input xi FITS file name')
    parser.add_option('-l', '--line_count', dest='line_count', action='store_true', default=False,
                       help='Only output number of lines')
    parser.add_option('-o', '--output', dest='output', type='string', default=None,
                      help='output base name (default = stdout), to which is added \'_<hdu>.txt\', or \'.txt\' if --hdu option is given')
    parser.add_option('-n', '--n_row', dest='n_row', type='int', default=-1, help='print only first n_row rows')
    parser.add_option('-e', '--every', dest='every', type='int', default=-1, help='print only every EVERY rows')
    parser.add_option('-s', '--start', dest='start', type='int', default=-1, help='start at row START')
    parser.add_option('', '--hdu', dest='hdu', type='int', default=-1, help='print only hdu HDU')
    parser.add_option('-H', '--header', dest='header', type='int', default=2, help='Prints full (2; default), minimal (1) or no header (0)')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true',
                      help='Verbose')


    options, args = parser.parse_args()

    see_help = 'See option \'-h\' for help.'

    if options.input is None:
        print('Input file not given (use option \'-i\'). ' + see_help, file=sys.stderr)
        return

    if options.start < 0:
        options.start = 0

    read_and_print_fits_file(options.input, output=options.output, verbose=options.verbose, n_row=options.n_row,
            start=options.start, every=options.every, out_hdu=options.hdu, header=options.header, line_count=options.line_count)



if __name__ == "__main__":
    sys.exit(main(sys.argv))


