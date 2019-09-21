"""

:Name: mkstuff.py

:Author: Martin Kilbinger, <martin.kilbinger@cea.fr>

:Package: mkstuff

"""


# Compability with python2.x for x>6
from __future__ import print_function, unicode_literals

import re
import numpy as np
import math
import os
import errno

import subprocess
import shlex
import sys

from scipy import stats

from optparse import IndentedHelpFormatter


r"""
=======
Classes
=======
"""

class param:
    """General class to store (default) variables
    """

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        print(self.__dict__)

    def var_list(self, **kwds):
        return vars(self)


class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


###########
### I/O ###
###########


def read_table(fname, count_comment=False, type='array'):
    """Reads the file 'fname', skips comment lines
    (starting with '#') and returns a numpy matrix.
    If count_comment=True, also returns number of comment
    lines and the header.
    Return type is ndarray (for type = 'array', 'ndarray')
    or matrix (for type ='matrix').
    """

    data = []
    header = []
    ncomment = 0
    f = open(fname, "r")
    all_lines = f.readlines()
    f.close()
    for line in all_lines:
        m = re.search("#", line)
        if m is None:
            data.append([])
            all_el = line.split()
            for el in all_el:
                data[-1].append(float(el))
        else:
            ncomment += 1
            header.append(line)

    mdata = np.matrix(data)
    if type == 'array' or type == 'ndarray':
        mdata = mdata.getA()

    if count_comment == True:
        return mdata, ncomment, header
    else:
        return mdata



def get_column(name, header, verbose=True):
    """Returns the column number for field 'name', from ascii header
    """

    header_fields = header[0].replace('#', '').split()
    for i in range(len(header_fields)):
        if name == header_fields[i]:
            return i
    if verbose is True:
       warning('Field {0} not found in header'.format(name)) 

    return -1



def get_data_col(data, col_name, action='exit', format='fits'):
    """Return the column with name *col_name* from file *data*.

    Parameters
    ----------
    data: table file struct
        Contains table information and data
    col_name: string
        Column name
    action: string
       *action* one of 'exit', 'warn', 'none'/None.
    format: string
        one of 'fits', 'ascii'

    Returns
    -------
    column: array
        Data array corresponding to column *col_name*
    """

    found = True

    if format == 'fits':
        if not col_name in data.names:
            found = False
    elif format == 'ascii':
        if not col_name in data.keys():
            found = False
    else:
        error('Unsupported file format {}'.format(format))

    if found == False:
        msg = 'Column \'{}\' not found in {} file'.format(col_name, format)
        if action == 'exit':
            error(msg)
        elif action == 'warn':
            warning(msg)
            return None
        else:
            return None

    return data[col_name]



def mkdir_p(path, verbose=False):
    """Create diretorcy by calling os.makedirs. Emulate shell function 'mkdir -p':
       If path already exists, returns without raising an error.

    Parameters
    ----------
    path: string
        Directory name
    verbose: boolean
        Verbose mode, default False

    Returns
    -------
    None
    
    """

    if verbose is True:
        print('Creating directory \'{}\''.format('{}'.format(path)))

    try:
        os.makedirs(str(path))
    #except OSError, exc: # Python <2.5
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise



def chdir(path, verbose=False):
    """Change directory by calling os.chdir. Emulate shell function 'cd'.

    Parameters
    ----------
    path: string
        Directory name
    verbose: boolean
        Verbose mode, default False

    Returns
    -------
    None
    """

    if verbose is True:
        print('Changing to directory \'{}\''.format(path))

    os.chdir(path)



def ln_s(orig, new, verbose=False, force=False):
    """Create symbolic link.

    Parameters:
    -----------
    orig: string
        Name of original file
    new: string
        Name of new, to be created, link
    verbose: bool
        Verbose output
    force: bool
        If True, link creation is forced even if file exists

    Returns:
    --------
    None
    """

    if os.path.isfile(orig) or os.path.isdir(orig):
        if os.path.isfile(new) or os.path.isdir(new):
            if force == False:
                if verbose:
                    print('File \'{}\' exists, skipping link creation...'.format(new))
            else:
                if verbose:
                    print('File \'{}\' exists, deleting file and create new link...'.format(new))
                os.remove(new)
                os.symlink(orig, new)
        else:
            if verbose:
                print('Creating link \'{}\' <- \'{}\''.format(orig, new))
            os.symlink(orig, new)
    else:
        if verbose:
            print('Original file \'{}\' does not exist, skipping...'.format(orig))



def write_matrix(a, name):
    """Writes the matrix a to the file name.
    """

    Nx, Ny = a.shape
    f = open(name, 'w')
    for x in range(Nx):
        for y in range(Ny):
            #print(a[x,y] + ' ', file=f, end='')
            print('{} '.format(a[x,y]), file=f, end='')
        print('', file=f)
    f.close()



def number_of_lines(fname, count_comment_lines=False):
    """Return number of lines in file.

    Parameters
    ----------
    fname: string
        Input file name
    count_comment_lines: boolean
        If False (default), do not include comment lines in count

    Returns
    -------
    num_lines: integer
        Number of lines in file (0 for empty file)
    """

    if count_comment_lines is True:
        num_lines = sum(1 for line in open(fname))
    else:
        num_lines = 0
        pattern   = '\s*#' # Comment with 0 or more preceding whitespace characters
        for line in open(fname):
            if re.match(pattern, line):   # Look for match at beginning of string
                continue
            num_lines += 1

    return num_lines



###################
### Job control ###
###################

def log_command(argv, name=None, close_no_return=True):
    """Write command with arguments to a file or stdout.
       Choose name = 'sys.stdout' or 'sys.stderr' for output on sceen.

    Parameters
    ----------
    argv: array of strings
        Command line arguments
    name: string
        Output file name (default: 'log_<command>')
    close_no_return: bool
        If True (default), close log file. If False, keep log file open
        and return file handler

    Returns
    -------
    log: filehandler
        log file handler (if close_no_return is False)
    """

    if name is None:
        name = 'log_' + os.path.basename(argv[0])

    if name == 'sys.stdout':
        f = sys.stdout
    elif name == 'sys.stderr':
        f = sys.stderr
    else:
        f = open(name, 'w')

    for a in argv:

        # Quote argument if special characters
        if ']' in a or ']' in a:
            a = '\"{}\"'.format(a)

        print(a, end='', file=f)
        print(' ', end='', file=f)

    print('', file=f)

    if close_no_return == False:
        return f

    if name != 'sys.stdout' and name != 'sys.stderr':
        f.close()



def error(str, val=1, stop=True, verbose=True):
    """Print message str and exits program with code val.
       See [ABC:]covest.py for 2.7 version.

    Parameters
    ----------
    str: string
        message
    val: integer
        exit value, default=1
    stop: boolean
        stops program if True (default), continues if False
    verbose: boolean
        verbose output if True (default)

    Returns
    -------

    None
    """

    if verbose is True:
        print_color('red', str, file=sys.stderr, end='')

    if stop is False:
        if verbose is True:
            print_color('red', ', continuing', file=sys.stderr)
    else:
        if verbose is True:
            print_color('red', '', file=sys.stderr)
        sys.exit(val)



def warning(str):
    """Prints message to stderr
    """

    error('Warning: ' + str, val=None, stop=False)



def run_cmd_list(cmd_list, ncpu=1, run=True, verbose=True, stop=False, file_list=None):
    """Run list of commands cmd_list on ncpu CPUs in parallel, calling
       run_cmd.

    Parameters
    ----------
    cmd_list: array of strings
        list of commands
    ncpu: int
        Number of CPUs, default = 1
    run: bool
        If True (default), run commands. run=False is for testing and debugging purpose
    verbose: bool
        If True (default), verbose output
    stop: bool
        If False (default), do not stop after command exits with error.
    file_list: array of strings
        If file_list[i] exists, cmd_list[i] is not run. Default value is None

    Returns
    -------
    sum_ex: int
        Sum of exit codes of all commands
    """

    njob = len(cmd_list)

    if njob == 0:
        return 0

    if verbose is True and njob > 1:
        if ncpu == 1:
            str_ncpu = 'CPU'
        else:
            str_ncpu = 'CPUs'
        print('Running {} jobs on {} {} in parallel'.format(njob, ncpu, str_ncpu))

    naverage = ncpu
    nchunk   = int(njob / ncpu)
    if njob % ncpu != 0:
        nchunk += 1
    files    = []
    idx      = []

    i_start = 0
    i_stop  = 0
    for i in range(nchunk):
        idx.append([])

        i_stop += naverage
        if i_stop > njob:
            i_stop = njob

        idx[i].extend(np.arange(i_start, i_stop))

        i_start = i_stop

    if verbose is True:
        print('idx    = {}, max/min njob / cpu = {}/{}'.format(len(idx), len(idx[0]), len(idx[-1])))

    sum_ex = 0
    for i in idx:
        if file_list is None:
            files = None
        else:
            files = [file_list[j] for j in i]
        res, out_list, err_list = run_cmd([cmd_list[j] for j in i], run=run, verbose=verbose, stop=stop,
                          file_list=files, parallel=True)
        sum_ex = sum_ex + res

    return sum_ex



def run_cmd(cmd_list, run=True, verbose=True, stop=False, parallel=True, file_list=None, devnull=False):
    """Run shell command or a list of commands using subprocess.Popen().

    Parameters
    ----------

    cmd_list: string, or array of strings
        list of commands
    run: bool
        If True (default), run commands. run=False is for testing and debugging purpose
    verbose: bool
        If True (default), verbose output
    stop: bool
        If False (default), do not stop after command exits with error.
    parallel: bool
        If True (default), run commands in parallel, i.e. call subsequent comands via
        subprocess.Popen() without waiting for the previous job to finish.
    file_list: array of strings
        If file_list[i] exists, cmd_list[i] is not run. Default value is None
    devnull: boolean
        If True, all output is suppressed. Default is False.

    Returns
    -------
    sum_ex: int
        Sum of exit codes of all commands
    """

    if type(cmd_list) is not list:
        cmd_list = [cmd_list]

    if verbose is True and len(cmd_list) > 1:
        print('Running {} commands, parallel = {}'.format(len(cmd_list), parallel))


    ex_list   = []
    pipe_list = []
    out_list  = []
    err_list  = []
    for i, cmd in enumerate(cmd_list):

        ex = 0
        out = ''
        err = ''

        if run is True:

            # Check for existing file
            if file_list is not None and os.path.exists(file_list[i]):
                if verbose is True:
                    print_color('blue', 'Skipping command \'{}\', file \'{}\' exists'.format(cmd, file_list[i]))
            else:
                if verbose is True:
                        print_color('green', 'Running command \'{0}\''.format(cmd))

                # Run command
                try:
                    cmds = shlex.split(cmd)
                    if devnull is True:
                        pipe = subprocess.Popen(cmds, stdout=subprocess.DEVNULL)
                    else:
                        # MKDEBUG NEW 17/04/2018
                        pipe = subprocess.Popen(cmds, stdout=subprocess.PIPE)
                        # See https://www.endpoint.com/blog/2015/01/28/getting-realtime-output-using-python
                        while True:
                            output = pipe.stdout.readline().decode('UTF-8')
                            if output == '' and pipe.poll() is not None:
                                break
                            if output:
                                print(output.strip())
                            ex = pipe.poll()

                    if parallel is False:
                        # Wait for process to terminate
                        pipe.wait()

                    pipe_list.append(pipe)

                    # If process has not terminated, ex will be None
                    #ex = pipe.returncode
                except OSError as e:
                    print_color('red', 'Error: {0}'.format(e.strerror))
                    ex = e.errno

                    check_error_stop([ex], verbose=verbose, stop=stop)

        else:
            if verbose is True:
                print_color('yellow', 'Not running command \'{0}\''.format(cmd))

        ex_list.append(ex)
        out_list.append(out)
        err_list.append(err)


    if parallel is True:
        for i, pipe in enumerate(pipe_list):
            pipe.wait()

            # Update exit code list
            ex_list[i] = pipe.returncode


    s = check_error_stop(ex_list, verbose=verbose, stop=stop)

    #return s
    return s, out_list, err_list



def check_error_stop(ex_list, verbose=True, stop=False):
    """Check error list and stop if one or more are != 0 and stop=True

    Parameters
    ----------
    ex_list: list of integers
        List of exit codes
    verbose: boolean
        Verbose output, default=True
    stop: boolean
        If False (default), does not stop program

    Returns
    -------
    s: integer
        sum of absolute values of exit codes
    """

    if ex_list is None:
        s = 0
    else:
        s = sum([abs(i) for i in ex_list])


    # Evaluate exit codes
    if s > 0:
        n_ex = sum([1 for i in ex_list if i != 0])
        if verbose is True:
            if len(ex_list) == 1:
                print_color('red', 'The last command returned sum|exit codes|={}'.format(s), end='')
            else:
                print_color('red', '{} of the last {} commands returned sum|exit codes|={}'.format(n_ex, len(ex_list), s), end='')
        if stop is True:
            print_color('red', ', stopping')
        else:
            print_color('red', ', continuing')

        if stop is True:
            sys.exit(s)


    return s



def do_job(options, job, verbose=True):
    """Returns true if job #job is within range given in options.
       From ~/bin/Great3/run_sex_moments.py.
    """

    res =  (options.job_first <= job) and (options.job_last >= job)

    if verbose is True:
        if res is True:
            print_color('green', 'Running job  #{}: {}'.format(job, options.jobs[job]))
        else:
            print_color('blue', 'Skipping job #{}: {}'.format(job, options.jobs[job]))

    return res


def job_description(jobs):
    """Return description of jobs for usage string.
    """

    job_last = len(jobs) - 1

    if job_last <= 0:
        return '', -1

    job_list = 'Jobs:\n'
    for j, descr in enumerate(jobs):
        job_list += ' {:3d} {}\n'.format(j, descr)

    return job_list, job_last



def get_verbose_flag(verbose):
    """Return verbose flag ' -v' if verbose=True, and empty string otherwise,
       for programs executed as shell commands.

    Parameters
    ----------
    verbose: bool
        Verbose mode

    Returns
    -------
    verbose_flag: string
        Verbose flag for shell command
    """

    if verbose is True:
        verbose_flag = ' -v'
    else:
        verbose_flag = ''

    return verbose_flag


############
### Maths ##
############


# Returns density value levels of L corresponding
# to confidence levels cl
def get_density_levels(L, cl=None):
    if cl == None: cl   = [0.6827, 0.9545, 0.9973]

    Ls   = np.sort(L, axis=None)[::-1]
    Lsum = Ls.sum()
    N    = Ls.shape[0]

    cum  = 0
    levels = np.zeros(shape=len(cl)) - 1
    for i in range(N):
        cum = cum + Ls[i]
        for l in range(len(cl)):
            if cum >= cl[l] * Lsum and levels[l] < 0:
                levels[l] = Ls[i]
    return levels


# x can be list or np.array
def mean(x):
    return float(sum(x))/len(x)


def corr_coeff(a):
    """Return correlation matrix (correlation coefficient)
       of matrix a.

    Parameters
    ----------
    a: matrix of float
        input matrix

    Returns
    -------
    ra: matrix of float
        correlation matrix

    Raises
    ------
    ValueError
        if input is not matrix
    ZeroDivisionError
        if an input diagonal element is zero
    """


    try:
        Nx, Ny = a.shape
    except ValueError as err:
        print('Input is not a 2D matrix')
        raise

    if any(np.diag(a) == 0):
        raise ZeroDivisionError('At least one input diagonal element is zero')

    ra = np.zeros(shape=a.shape) 
    for i in range(Nx):
        for j in range(Ny):
            ra[i,j] = a[i,j] / np.sqrt(a[i,i] * a[j,j])
    return ra


def frexp10(x):
    """Return the mantissa and exponent of x, as pair (m, e).
       m is a float and e is an int, such that x = m * 10.0**e.
       See math.frexp()

        >>> mkstuff.frexp10(1240)
        (1.24, 3)

    """

    if x == 0: return (0, 0)
    try:
        l = math.log10(abs(x))
    except:
        print('Error with math.log10(|' + (str(x)) + '|)')
        return None, None
    if l < 1: l = l - 1 + 1e-10
    exp = int(l)
    return x / 10**exp, exp 


def exp2string_nice(x):
    """Return nicely formatted string in scientific notation.
    """

    m, e = frexp10(x)

    if m == 1:
        s = '10^{}'.format(e)
    else:
        s = '{:g}\\times10^{}'.format(m, e)

    return s


def draw_from_hist(bins, values, n):
    """Draws n random numbers from a histogram (bins, values).
       The np.arrays bins and values must have the same length.
     """

    v = values / sum(values)
    cdf = np.add.accumulate(v)
    rd  = np.digitize(np.random.random_sample(n), cdf)

    return bins[rd]



def get_bin_idx(values, x, which='lower', mode=None):
    """Returns index i such that bins[i] <= x < bin[i+1].
       Also returns out-of-range flag that is set to one if
       x is not in the range of values
       if which='right' or 'higher', returns i+1.
       If mode == 'stop_on_error', stops if x is out of range.
       Uses bisection.
    """

    i_low  = 0
    i_high = len(values) - 1

    if x < values[i_low]:

        out_of_range = 1
        i = 0

    elif x >= values[i_high]:

        out_of_range = 1
        i = i_high

    else:

        out_of_range = 0
        while (i_high - i_low > 1):
            i = (i_high + i_low) >> 1
            if values[i] > x:
                i_high = i
            else:
                i_low = i

        if which == 'lower' or which == 'left':
            i = i_low
        elif which == 'higher' or which == 'right':
            i = i_high
        else:
            error('Wrong mode which={0}'.format(which))

    if out_of_range == 1 and mode == 'stop_on_error':
        error('Value {0} out of range {1} .. {2}'.format(x, values[0], values[-1]))

    return i, out_of_range



def is_number(str, type=None):
    """Returns None if str is not a number. Returns type(str) otherwise.
       type can be 'int(eger)' or 'float'
    """

    try:
        f = float(str)
        is_num = True
    except ValueError:
        is_num = False

    if is_num == False:
        return None

    if type is None:
        return True

    elif type == 'integer' or type == 'int':
        try:
            f_int = int(str)
            return f_int
        except ValueError:
            return None

    elif type == 'float':
        return f


def column(matrix, i):
    """Returns i-th column from two-dimensional list matrix
    """
    return [row[i] for row in matrix]



###################
### Coordinates ###
###################


coord = {'rad'    : 1,
         'deg'    : np.pi / 180.0,
         'arcmin' : np.pi / 180.0 / 60.0,
         'arcsec' : np.pi / 180.0 / 60.0 / 60.0,
         'none'   : 1}


def coord_to_rad(theta, c):
    """Transform angle theta[c] to theta[rad]
    """

    return theta * coord[c]



def rad_to_coord(theta, c):
    """Transform angle theta[rad] to theta[c]
    """

    return theta / coord[c]



#############
### Other ###
#############


# From http://code.activestate.com/recipes/577279-generate-list-of-numbers-from-hyphenated-and-comma/
def hyphen_range(s, first=None, last=None):
    """ creates a list with each integer from a complex range string like "1-9,12, 15-20,23".
        The wildcard '*' returns [first, last+1].

    >>> list(hyphen_range('1-9,12, 15-20,23'))
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16, 17, 18, 19, 20, 23]

    >>> list(hyphen_range('1-9,12, 15-20,2-3-4'))
    Traceback (most recent call last):
        ...
    ValueError: format error in 2-3-4
    """

    res = []

    if s == '*':
        if first is None or last is None:
            raise ValueError('Wildcard range requires first and last of index list')
        res = list(range(first, last+1))

    else:
        for x in s.split(','):
            elem = x.split('-')
            if len(elem) == 1: # a number
                res.append(int(elem[0]))
            elif len(elem) == 2: # a range inclusive
                start, end = map(int, elem)
                for i in range(start, end+1):
                    res.append(i)
            else: # more than one hyphen
                raise ValueError('format error in %s' % x)

    return res


def my_string_split(string, num=-1, sep_force=None, verbose=False, stop=False):
    """Split a *string* into a list of strings. Choose as separator
        the first in the list [space, underscore] that occurs in the string.
        (Thus, if both occur, use space.)

    Parameters
    ----------
    string: string
        Input string
    num: int
        Required length of output list of strings, -1 if no requirement.
    sep_force: string, optional, default=None
        if not None, use this separator instead of ' ' or '_'
    verbose: bool, optional, default=False
        Verbose output
    stop: bool, optional, default=False
        Stop programs with error if True, return None and continues otherwise

    Raises
    ------
    MyError
        If number of elements in string and num are different, for stop=True

    Returns
    -------
    list_str: string, array()
        List of string on success, and None if failed.
    """

    if string is None:
        return None

    if sep_force:
        has_sep = string.find(sep_force)
        if has_sep != -1:
            sep = sep_force
        else:
            # string has neither, consists of one element
            if num == -1 or num == 1:
                # one-element string is ok
                sep = None
            else:
                error('Separator \'{} \' not found in string \'{}\', cannot split'.format(string))

    else:
        has_space      = string.find(' ')
        has_underscore = string.find('_')

        if has_space != -1:
            # string has white-space
            sep = ' '
        else:
            if has_underscore != -1:
            # string has no white-space but underscore
                sep = '_'
            else:
                if num == -1 or num == 1:
                    sep = None
                else:
                    error('Neither \' \' nor \'_\' found in string \'{}\', cannot split'.format(string))
 
    #res = string.split(sep=sep) # python v>=3?
    res = string.split(sep)

    if num != -1 and num != len(res) and stop==True:
        raise MyError('String \'{}\' has length {}, required is {}'.format(string, len(res), num))

    return res


def print_color(color, txt, file=sys.stdout, end='\n'):
    """Print text with color. If not supported, print standard text.

    Parameters
    ----------
    color: string
        color name
    txt: string
        message
    file: file handler
        output file handler, default=sys.stdout
    end: string
        end string, default='\n'

    Returns
    -------
    None
    """

    try:
        import colorama
        colors = {'red'    : colorama.Fore.RED,
                  'green'  : colorama.Fore.GREEN,
                  'blue'   : colorama.Fore.BLUE,
                  'yellow' : colorama.Fore.YELLOW,
                  'black'  : colorama.Fore.BLACK,
                 }

        if colors[color] is None:
            col = colorama.Fore.BLACK
        else:
            col = colors[color]

        print(col + txt + colors['black'] + '', file=file, end=end)

    except ImportError:
        print(txt, file=file, end=end)

