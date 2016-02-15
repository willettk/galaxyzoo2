# ppgplot_numpy_spb.py

# imports ppgplot functions and sets up useful variables

from ppgplot import *
import numpy as N
import os
from os.path import join as pjoin
import pyfits
from glob import glob
from gaussian import gaussian2d
import sys
#import nr as nr
import rostat as rostat
from numpy import median
from scipy.stats import scoreatpercentile
from checkarray import checkarray

pgplot_extras_path = '/Users/spb/Software/pgplot_extras'

# Point style to use for plots
pointStyles = {"dot": 17, "point": 1, "square": 0, "plus": 2, "asterisk": 3,
               "circle": 22, "cross": 5, "triangle": 7, "star": 12, "tiny": -1,
               "filled-square": 16, "filled-circle": 17, "filled-star": 18,
               "filled-triangle": 13,
               "arrow-up": 30, "arrow-down": 31,
               "arrow-left": 28, "arrow-right": 29}
def pgxpt(x, y, s):
    x = N.asarray(x)
    y = N.asarray(y)
    ch1 = None
    if s.endswith('filled-square'):
        ch1 = pgqch()
	pgsch(1.25*ch1)
    if s.startswith('small'):
	ch = pgqch()
	pgsch(0.75*ch)
	pgpt(x, y, pointStyles[s.replace('small', '')])
	pgsch(ch)
    elif s.startswith('tiny'):
	ch = pgqch()
	pgsch(0.5*ch)
	pgpt(x, y, pointStyles[s.replace('tiny', '')])
	pgsch(ch)
    else:
        pgpt(x, y, pointStyles[s])
    if not ch1 is None:
        pgsch(ch1)


# PGPLOT linestyles
lineStyles = {"solid": 1, "dashed": 2, "dash-dot": 3, "dotted": 4}
def pgxsls(s):  pgsls(lineStyles[s])

def pgxsci(s):  pgsci(colourIndices[s])

# OLD COLOUR SETUP ---------------------------------------------------
# # PGPLOT colour indices
# colourIndices = {"white": 0, "black": 1, "red": 2, "green": 3, "blue":4,
# 		 "brightcyan": 5, "magenta": 6, "yellow": 7, "orange": 8,
# 		 "lime": 9, "seafoam": 10, "cyan": 11, "purple": 12,
# 		 "strawberry": 13,
#                  "darkgray":14,"lightgray":15,
# 		 "darkred":16, "darkgreen":17, "darkblue":18,
# 		 "lightcyan":19,
# 		 "mocha": 20}
# #colourIndices = {"black": 1, "red": 1, "green": 1, "blue":1, "magenta": 1, "cyan": 1, "purple": 1, "darkgray":14,"lightgray":15}
# def pg_more_colours():
#     pgscr(16, 0.75, 0.0, 0.0)
#     pgscr(17, 0.0, 0.6, 0.0)
#     pgscr(18, 0.0, 0.0, 0.75)
#     pgscr(19, 0.5, 1.0, 1.0)
#     pgscr(20, 0.5, 0.25, 0.0)
# --------------------------------------------------------------------
# NEW COLOUR SETUP ---------------------------------------------------
colourIndices = {"white": 0, "black":1}
def pg_more_colours():
    from colours import colours
    ci = 2
    for c in colours:
        if c not in ['white', 'black']:
            r, g, b = colours[c]
            pgscr(ci, r, g, b)
            colourIndices[c] = ci
            ci += 1
# --------------------------------------------------------------------

# Default character height
ch = 1.3
# Default line width
lw = 3
# Use roman font
cf = 2

def pgerrbx(x, y, xerr, cap=1.0):
    x = checkarray(x)
    y = checkarray(y)
    xerr = checkarray(xerr)
    pgerrb(5, x, y, xerr, cap)

def pgerrby(x, y, yerr, cap=1.0):
    x = checkarray(x)
    y = checkarray(y)
    yerr = checkarray(yerr)
    pgerrb(6, x, y, yerr, cap)

def pgsetup(nx=1, ny=1):
    anx = abs(nx)
    any = abs(ny)
    if anx > 1 or any > 1:
	pgpap(anx*5.0, float(any)/anx)
	pgsubp(nx, ny)
    else:
	pgpap(0, 1.0)
    pg_more_colours()
    pgsci(colourIndices["black"])
    pgslw(lw); pgsch(ch); pgscf(cf)
    #pgeras()

def pgaqt():
    pgopen('/aqt')

def trend_bins(x, y, xlow=None, xhigh=None, xbinwidth=None, nmin=100,
               lowpc=2.5, highpc=97.5, usemedian=True):
    if xlow is None:  xlow = scoreatpercentile(x, 1)
    if xhigh is None:  xhigh = scoreatpercentile(x, 99)
    if xbinwidth is None:  xbinwidth = (xhigh - xlow) * 100*nmin / len(x)
    x_bin = N.arange(xlow, xhigh+xbinwidth/2.0, xbinwidth)
    n_bin = len(x_bin)
    xx_bin = N.zeros(n_bin, N.float) - 99999
    y_bin = N.zeros((3, n_bin), N.float) - 99999
    ok = N.ones(n_bin, N.bool)
    for i, xb in enumerate(x_bin):
        inbin = (x >= xb - 0.5*xbinwidth) & (x < xb + 0.5*xbinwidth)
        x_inbin = x[inbin]
        y_inbin = y[inbin]
        if len(y_inbin) > nmin:
            xx_bin[i] = median(x_inbin)
            if usemedian:
                y_bin[0, i] = median(y_inbin)
            else:
                y_bin[0, i] = N.mean(y_inbin)
            y_bin[1, i] = scoreatpercentile(y_inbin, lowpc)
            y_bin[2, i] = scoreatpercentile(y_inbin, highpc)
        else:
            ok[i] = False
    #x_bin = x_bin[ok]
    xx_bin = xx_bin[ok]
    y_bin = y_bin[:,ok]
    return xx_bin, y_bin

def bin_array(d, nbin, low, high, weights=None):
    n = len(d)
    step = (high-low)/float(nbin)
    x_bin = N.arange(nbin) * step + low + step/2.0
    if weights is None:
        d_bin = N.zeros(nbin)
        w = N.ones(n)
    else:
        d_bin = N.zeros(nbin, N.float)
        w = weights
    for i in range(n):
	bin_index = int((d[i] - low)/step)
	if 0 <= bin_index < nbin:
	    d_bin[bin_index] += w[i]
    return x_bin, d_bin

def bin_array_2d(x, y, nxbin, xlow, xhigh, nybin, ylow, yhigh,
		 completeness=None):
    nx = len(x)
    ny = len(y)
    if nx != ny:
	print 'Error: len(x) != len(y)'
	return
    xstep = float(xhigh-xlow)/nxbin
    ystep = float(yhigh-ylow)/nybin
    x_bin = N.arange(nxbin) * xstep + xlow + xstep/2.0
    y_bin = N.arange(nybin) * ystep + ylow + ystep/2.0
    d_bin = N.zeros((nybin, nxbin), N.float)
    for k in range(nx):
	jbin_index = int((x[k] - xlow)/xstep)
	ibin_index = int((y[k] - ylow)/ystep)
	if completeness is None:
	    c = 1
	else:
	    c = completeness[k]
	if 0 <= jbin_index < nxbin and 0 <= ibin_index < nybin:
	    d_bin[ibin_index, jbin_index] += 1.0/c
    return x_bin, y_bin, d_bin

def ave_array_2d(x, y, z, nxbin, xlow, xhigh, nybin, ylow, yhigh,
		 completeness=None):
    nx = len(x)
    ny = len(y)
    if nx != ny:
	print 'Error: len(x) != len(y)'
	return
    xstep = float(xhigh-xlow)/nxbin
    ystep = float(yhigh-ylow)/nybin
    x_bin = N.arange(nxbin) * xstep + xlow + xstep/2.0
    y_bin = N.arange(nybin) * ystep + ylow + ystep/2.0
    d_bin = N.zeros((nybin, nxbin), N.float)
    z_bin = N.zeros((nybin, nxbin), N.float)
    for k in range(nx):
	jbin_index = int((x[k] - xlow)/xstep)
	ibin_index = int((y[k] - ylow)/ystep)
	if completeness is None:
	    c = 1
	else:
	    c = completeness[k]
	if 0 <= jbin_index < nxbin and 0 <= ibin_index < nybin:
	    d_bin[ibin_index, jbin_index] += 1.0/c
            z_bin[ibin_index, jbin_index] += z[k]
    z_bin /= d_bin
    N.putmask(z_bin, d_bin < 1, 0.0)
    return x_bin, y_bin, z_bin


def min_array_2d(x, y, z, nxbin, xlow, xhigh, nybin, ylow, yhigh,
                 listoutput=False):
    nx = len(x)
    ny = len(y)
    if nx != ny:
	print 'Error: len(x) != len(y)'
	return
    xstep = float(xhigh-xlow)/nxbin
    ystep = float(yhigh-ylow)/nybin
    x_bin = N.arange(nxbin) * xstep + xlow + xstep/2.0
    y_bin = N.arange(nybin) * ystep + ylow + ystep/2.0
    d_bin = N.zeros((nybin, nxbin), N.float) + 1e100
    if listoutput:
        xl = N.zeros(nxbin*nybin, N.float)
        yl = N.zeros(nxbin*nybin, N.float)
        dl = N.zeros(nxbin*nybin, N.float)
    for k in range(nx):
	jbin_index = int((x[k] - xlow)/xstep)
	ibin_index = int((y[k] - ylow)/ystep)
	if (0 <= jbin_index < nxbin and 0 <= ibin_index < nybin
            and d_bin[ibin_index, jbin_index] > z[k]):
	    d_bin[ibin_index, jbin_index] = z[k]
            if listoutput:
                xl[jbin_index*nxbin+ibin_index] = (jbin_index+0.5)*xstep+xlow
                yl[jbin_index*nxbin+ibin_index] = (ibin_index+0.5)*ystep+ylow
                dl[jbin_index*nxbin+ibin_index] = z[k]
    if listoutput:
        return xl, yl, dl
    else:
        return x_bin, y_bin, d_bin


def bin_array_2d_smooth(x, y, nxbin, xlow, xhigh, nybin, ylow, yhigh,
			smooth_fwhm=10, completeness=None):
    nx = len(x)
    ny = len(y)
    if nx != ny:
	print 'Error: len(x) != len(y)'
	return
    xstep = float(xhigh-xlow)/nxbin
    ystep = float(yhigh-ylow)/nybin
    fwhm = smooth_fwhm
    sigma = fwhm/2.3548
    rmax = 4 * sigma
    size = int(round(rmax)) * 2 + 1
    gc = (size-1)/2
    g = N.asarray(gaussian2d((size, size), 1.0, gc, gc, sigma, sigma))
    x_bin = N.arange(nxbin) * xstep + xlow + xstep/2.0
    y_bin = N.arange(nybin) * ystep + ylow + ystep/2.0
    d_bin = N.zeros((nybin+2*gc, nxbin+2*gc), N.float)
    percent_prev = 0
    print '%:',
    for k in range(nx):
	percent = int((100.0*k)/nx)
	if percent/10 != percent_prev/10:
	    print percent,
	    percent_prev = percent
	    sys.stdout.flush()
	jbin_index = int((x[k] - xlow)/xstep)
	ibin_index = int((y[k] - ylow)/ystep)
	if completeness is None:
	    c = 1
	else:
	    c = completeness[k]
	if 0 <= jbin_index < nxbin and 0 <= ibin_index < nybin:
	    d_bin[ibin_index:ibin_index+2*gc+1,
		  jbin_index:jbin_index+2*gc+1] += g/c
    print
    d_bin = N.pi*(smooth_fwhm/2.0)**2 * d_bin[gc:-gc, gc:-gc]
    return x_bin, y_bin, d_bin


def bin_array_2d_super_smooth(x, y, nxbin, xlow, xhigh, nybin, ylow, yhigh,
			      n_per_smooth_fwhm=100,
			      completeness=None):
    nx = len(x)
    ny = len(y)
    if nx != ny:
	print 'Error: len(x) != len(y)'
	return
    xstep = float(xhigh-xlow)/nxbin
    ystep = float(yhigh-ylow)/nybin
    fwhm_list = (2.0 * N.sqrt(n_per_smooth_fwhm/(N.pi*N.arange(1,101))))
    sigma_list = fwhm_list/2.3548
    rmax = 4 * sigma_list[-1]
    size = int(round(rmax)) * 2 + 1
    gc = (size-1)/2
    g_list = []
    for sigma in sigma_list:	
	g = N.asarray(gaussian2d((size, size), 1.0, gc, gc, sigma, sigma))
	g_list.append(g)
    x_bin, y_bin, d_bin = bin_array_2d(x, y, nxbin, xlow, xhigh,
				       nybin, ylow, yhigh,
				       completeness)
    N.putmask(d_bin, d_bin < 1, 1)
    d_max = d_bin.max()
    gindex_array = (100*d_bin / d_max).astype(N.int)
    N.putmask(gindex_array, gindex_array < 0, 0)
    N.putmask(gindex_array, gindex_array > 99, 99)
    gindex_array_sorted = N.sort(N.ravel(gindex_array))
    gindex_min = gindex_array_sorted[0]
    gindex_max = gindex_array_sorted[-1]
    print 'Smoothing FWHM varies between',
    print '%3.2f and %3.2f bins'%(fwhm_list[gindex_min], fwhm_list[gindex_max])
    d_bin_smooth = N.zeros((nybin+2*gc, nxbin+2*gc), N.float)
    percent_prev = 0
    print '%:',
    for k in range(nx):
	percent = int((100.0*k)/nx)
	if percent/10 != percent_prev/10:
	    print percent,
	    percent_prev = percent
	    sys.stdout.flush()
	jbin_index = int((x[k] - xlow)/xstep)
	ibin_index = int((y[k] - ylow)/ystep)
	if 0 <= jbin_index < nxbin and 0 <= ibin_index < nybin:
	    gindex = gindex_array[ibin_index, jbin_index]
	    d_bin_smooth[ibin_index:ibin_index+2*gc+1,
			 jbin_index:jbin_index+2*gc+1] += g_list[gindex]
    print
    d_bin_smooth = N.pi*(fwhm_list[gindex_min]/2.0)**2 * d_bin_smooth[gc:-gc, gc:-gc]
    return x_bin, y_bin, d_bin_smooth, fwhm_list[gindex_min], fwhm_list[gindex_max]


def get_colour_table(name):
    filename = pjoin(pgplot_extras_path, 'tables', name+'.fits')
    rgb = pyfits.getdata(filename)
    return rgb

def get_colour_ramp(name):
    filename = pjoin(pgplot_extras_path, 'ramps', name+'.fits')
    i = pyfits.getdata(filename)
    return i

def get_colour_table_names():
    fileglob = pjoin(pgplot_extras_path, 'tables', '*.fits')
    filelist = glob(fileglob)
    names = []
    for f in filelist:
	n = f.split('/')[-1]
	n = n.replace('.fits', '')
	names.append(n)
    return names

def get_colour_ramp_names():
    fileglob = pjoin(pgplot_extras_path, 'ramps', '*.fits')
    filelist = glob(fileglob)
    names = []
    for f in filelist:
	n = f.split('/')[-1]
	n = n.replace('.fits', '')
	names.append(n)
    return names

def setup_colour_table(table='standard', ramp='ramp',
		       contra=1, bright=0.5, flip_colours=False,
                       inset_ends=0.0):
    l = get_colour_ramp(ramp)[::2]
    if flip_colours:
	flip = -1
    else:
	flip = 1
    r, g, b = get_colour_table(table)[:,::flip*2]
    nc = len(r)
    pgscir(32, 32+nc)
    pgctab(l, r, g, b, nc, contra, bright)

def plot_colour_schemes():
    pgopen('colour_schemes.ps/cps')
    ramps = get_colour_ramp_names()
    tables = get_colour_table_names()
    contra = 1
    bright = 0.5
    nstep = 128
    for ramp in ramps:
	l = get_colour_ramp(ramp)
	pgsetup()
	pgsci(1)
	pgpap(0, 2)
	pgslw(lw-1)
	pgsch(0.8*ch)
	pgsvp(0.2, 0.9, 0.1, 0.9)
	pgswin(-2, nstep+1, 0, 1.5*len(tables))
	pglab('', '', ramp)
	for i, t in enumerate(tables):
	    pgsci(1)
	    pgptxt(-nstep/100.0, 1.5*i+0.25, 0.0, 1.0, t)
	    x = N.arange(nstep)
	    y = N.ones(nstep) * 1.5*i + 0.5
	    setup_colour_table(t, ramp, flip_colours=True)
	    for j in range(nstep):
		pgsci(32+j)
		pgpt(x[j:j+1], y[j:j+1], pointStyles['filled-square'])
    pgclos()

def draw_circle(xc, yc, dx, dy, d):
    r = d/2.0
    x = (N.arange(0,101) / 100.0 - 0.5)*2.0
    y = N.sqrt(1-x**2)
    x = dx*x*r
    y = dy*y*r
    pgline(xc + x, yc + y)
    pgline(xc + x[::-1], yc - y[::-1])

def optimal_bin_width(x):
    # A first-order attempt at selecting an optimum histogram bin-width
    # from Scott (1979) as given in Wand (1997).
    n = len(x)
    sigma = x.std()
    lq, hq = rostat.quartiles(x)
    iqr = hq-lq
    sigma_hat = min(sigma, iqr/1.349)
    #print 'n =', n
    #print 'sigma =', sigma
    #print 'iqr =', iqr
    #print 'sigma_hat =', sigma_hat
    h = 3.49 * sigma_hat * n**(-1.0/3.0)
    #print 'bin_width =', h
    return h

def pgbinc(x, d):
    pgbin(x, d, True)

def ps2gif(f):
    """Convert ps to gif."""
    f = f.replace('/cps', '')
    f = f.replace('/ps', '')
    os.system('convert -density 576 -geometry 25% -rotate 90 ' +
	      '%s %s &'%(f, f.replace('.ps', '.gif')))
