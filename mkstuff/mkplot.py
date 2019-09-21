"""

:Name: mkplot.py

:Author: Martin Kilbinger, <martin.kilbinger@cea.fr>

:Package: mkstuff

"""

from mkstuff import *
import pylab as plt


################
### Plotting ###
################



### Colors ###

#darkred = RGB(100, 30, 30)


# From http://www.scipy.org/Cookbook/Matplotlib/CustomLogLabels
def log_10_product(x, pos):
    """The two args are the value and tick position.
    Label ticks with the product of the exponentiation"""
    res = '$%1i$' % (x)
    return res


def lin_product(x, pos):
    thres = 1
    if abs(x) < thres:
        m, e = frexp10(x)
        if (m, e) == (0, 0): res = '$0\,$'
        else:
            if m == 1: res = '$10^{%g}$' % (e)
            else: res = '$%g \\cdot 10^{%g}$' % (m, e)
    else:
        #res = str(x)
        res = '{0:g}'.format(x)
    return res


def no_product(x, pos):
    thres = 1
    if abs(x) < thres:
        m, e = frexp10(x)
        if (m, e) == (0, 0): res = '$0\,$'
        else: res = '$10^{%g}$' % e
    else:
        res = str(x)
    return res



def set_axis(axis, linlog, ax=None):
    """Axis setting to linear or logarithmic.
       Seems not to work any more (python 3 problem?).
       See also make_log_ticks.
    """

    if ax is None:
        ax = plt.subplot(111, axes='equal')

    formatter = plt.FuncFormatter(lin_product)

    if axis == 'x':
        print(ax)
        if linlog == 'log': ax.set_xscale(linlog, nonposx='clip')
        ax.xaxis.set_major_formatter(formatter)
    if axis == 'y':
        if linlog == 'log': ax.set_yscale(linlog, nonposy='clip')
        ax.yaxis.set_major_formatter(formatter)
        
    return ax



def set_square_plot():
    """Set axis ratio such that a square plot is created.
    To avoid a too large bounding box in the pdf figure of the resulting plot, use
    the options 'bbox_inches = 'tight'' in plt.savefig.

    Parameters
    ----------
    None

    Retunrs
    -------
    None
    """

    plt.axes().set_aspect((plt.xlim()[1] - plt.xlim()[0]) / (plt.ylim()[1] - plt.ylim()[0]))



def set_equal_axes_plot():
    """Set axis ratio such that equal dx and dy are drawn equal in the plot.

    Parameters
    ----------
    None

    Retunrs
    -------
    None
    """

    plt.gca().set_aspect('equal', adjustable='box')



def update_lim(xy, min_new=None, max_new=None):
    """Updates x- or y-axis limit.
    """

    if xy == 'x':
        min_cur, max_cur = plt.xlim()
    elif xy == 'y':
        min_cur, max_cur = plt.ylim()
    else:
        mkstuff.error('Wrong axis type {}'.format(xy))

    if min_new is None:
        min_new = min_cur
    if max_new is None:
        max_new = max_cur

    if xy == 'x':
        plt.xlim(min(min_cur, min_new), max(max_cur, max_new))
        min_upd, max_upd = plt.xlim()
    elif xy == 'y':
        plt.ylim(min(min_cur, min_new), max(max_cur, max_new))
        min_upd, max_upd = plt.ylim()

    return min_upd, max_upd

    

def make_log_ticks(xy, log_min, log_max, first_call=True):
    """Creates logarithmic ticks.
    """ 

    # Nearest integer for axis range
    log_ax_min   = np.floor(log_min)
    log_ax_max   = np.ceil(log_max)

    # Set (and update) limits
    if xy == 'x':
        if first_call == True:
            min_cur = log_ax_min
            max_cur = log_ax_max
        else:
            min_cur, max_cur = plt.xlim()
        plt.xlim(min(min_cur, log_ax_min), max(max_cur, log_ax_max))
        min_new, max_new = plt.xlim()
    elif xy == 'y':
        if first_call == True:
            min_cur = log_ax_min
            max_cur = log_ax_max
        else:
            min_cur, max_cur = plt.ylim()
        plt.ylim(min(min_cur, log_ax_min), max(max_cur, log_ax_max))
        min_new, max_new = plt.ylim()


    # Set ticks and their labels
    ticks        = np.arange(min_new, max_new + 1, 1)
    ticklabels   = []
    for t in ticks:
        ticklabels.append('$10^{{{:.0f}}}$'.format(t))

    minor_ticks = []
    for i in range(int(log_ax_min), int(log_ax_max)):
        for j in range(2, 10):
            minor_ticks.append(i + np.log10(j))

    if xy == 'x':
        plt.xticks(ticks, ticklabels)
        plt.gca().set_xticks(minor_ticks, minor=True)
    elif xy == 'y':
        plt.yticks(ticks, ticklabels)
        plt.gca().set_yticks(minor_ticks, minor=True)
    else:
        mkstuff.error('Wrong axis name {}'.format(xy))



def axes_limits_sym(lim=1):
    """Sets symmetrix axes limits.
    """

    plt.xlim((-lim, lim))
    plt.ylim((-lim, lim))


def axes_no_offset():
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)


def plot_one(name, col, style='', drawstyle='default', label=None, linewidth=1):
    dat, nc, header  = read_table(name, count_comment=True)
    x    = dat[:,0]
    y    = dat[:,col]
    my_plot = plt.plot(x, y, style, drawstyle=drawstyle, label=label, linewidth=linewidth)
    return header, my_plot


def plot_one_err(name, col, col_err, color=None):
    dat, nc, header  = read_table(name, count_comment=True)
    X    = dat[:,0]
    Y    = dat[:,col]
    DY   = dat[:,col_err]
    x    = np.array(X)[:,0]
    y    = np.array(Y)[:,0]
    dy   = np.array(DY)[:,0]
    plt.errorbar(x, y, yerr=dy, fmt=None, ecolor=color)


def imshow(img, scale_val=1, ax=None, *args, **kwargs):
    """Tweaked version of imshow, to also print z-value of image
    """

    if ax is None:
         ax = plt.gca()
    im = ax.imshow(img, *args, **kwargs)
    numrows, numcols = img.shape

    # Need lambda function to access to img data
    # TODO: account for extent dimensions
    ax.format_coord = lambda x, y: ('x=%1.2f    y=%1.2f    '%(x, y), 'x=%1.2f    y=%1.2f    z=%1.4g    '%(x, y, img[int(x+0.5),int(y+0.5)]))[x>=-0.5 and x<=numcols-0.5 and y>=-0.5 and y<=numrows-0.5]

    ax.figure.canvas.draw()
    return im


def shear_plot(x, y, gamma1, gamma2, fwidth=1, scale=1.0, every=1, gamma_min=0.01, title=None, out_name=None, unit=None):
    """Create stick plot of shear data.
       See https://github.com/drphilmarshall/Pangloss/blob/c9e4efb51d51cf7f27b10fd9563bf2da86adde23/pangloss/shearmap.py
       for a more professional code. See also yorick/lensing.i:shearmap2.

    Parameters
    ----------
    x: array of doubles
        x coordinates
    y: array of doubles
        y coordinates
    gamma1: array of doubles
        gamma1 values
    gamam2: array of doubles
        gamma2 values
    fwidth: double, optional
        width factor for  sticks, default=1
    scale: double, optional
        scale factor for sticks, default=1
    every: integer
        plot a fraction of 1/every objects, default=1
    gamma_min: double
        minimum shear module plotted, default=0.01
    title: string, optional
        title string, default None
    out_name: string, optional
        output name, default None
    unit: string, optional
        coordinate units, for axis labels, default None

    Returns
    -------
    None
    """

    # Prepare plot if output name given
    if out_name is not None:
        plt.clf()
        plt.xlabel('ra [{}]'.format(unit))
        plt.ylabel('dec [{}]'.format(unit))
        if title is not None:
            plt.title(title)

    # TODO: pixellise data, plot mean per pixel instead of every key

    # Plot only some objects
    if every > 1:
        ev_list = np.arange(0, len(x), every)
        x       = x[ev_list]
        y       = y[ev_list]
        gamma1  = gamma1[ev_list]
        gamma2  = gamma2[ev_list]

    mod_gamma = np.sqrt(gamma1*gamma1 + gamma2*gamma2)

    # Do not plot very small sticks
    ind_min   = np.where(mod_gamma<gamma_min)
    mod_gamma = np.delete(mod_gamma, ind_min)
    x         = np.delete(x, ind_min)
    y         = np.delete(y, ind_min)
    gamma1    = np.delete(gamma1, ind_min)
    gamma2    = np.delete(gamma2, ind_min)

    phi_gamma = np.arctan2(gamma2, gamma1) / 2.0

    dx        = scale * mod_gamma * np.cos(phi_gamma)
    dy        = scale * mod_gamma * np.sin(phi_gamma)

    # 0.005 is recommended in quiver help
    width     = (plt.xlim()[1] - plt.xlim()[0]) * 0.005

    #plt.figure(figsize=(70, 70))
    plt.quiver(x, y, dx, dy, headwidth=0, width=fwidth*width, headlength=0, headaxislength=0, minshaft=0, minlength=0, pivot='middle')
    set_equal_axes_plot()

    if out_name is not None:
        import plotly
        from plotly.tools import FigureFactory as FF
        print('FF.create_quiver')
        fig = FF.create_quiver(x, y, dx, dy, arrow_scale=0.00001, auto_open=False)
        print('plotly create plot')
        plotly.potly.plot(fig, filename='out_name')
        print('done')

    if out_name is not None:
        plt.savefig(out_name, dpi=900)



def shear_wheel():
    """Create 'wheel' plot of shear, see Kilbinger et al. (2015) Fig. 2.
    """

    N = 13
    x = 0.5 * np.cos(np.arange(N)/(N-1) * 2 *np.pi)
    y = 0.5 * np.sin(np.arange(N)/(N-1) * 2 *np.pi)
    gamma1 = x
    gamma2 = y

    shear_plot(x, y, gamma1, gamma2, fwidth=0.5)

    plt.savefig('wheel.pdf')



