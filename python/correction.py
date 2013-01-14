from ppgplot_spb import *
import pyfits
import numpy as N
import cosmology
from math import log10, e, exp, pi
import sys
from scipy.interpolate.fitpack import *
from scipy.interpolate import interp1d
from scipy.optimize import *
from scipy.integrate import *
from scipy.stats import scoreatpercentile
from scipy.special.orthogonal import chebyt
from scipy.ndimage import *
from astrotools import calc_ang_dist
import cPickle as pickle
import time
#import fitmc
import gc
import sys
from combine_fits_tables import combine_fits_tables
#import get_zoo_data

data_path = '/Users/spb/Work/projects/GalaxyZoo/dr6z7fiddleup/'
data_file = 'dr6z7_zoo_selected_counts.fits'
plots_path = data_path+'plots/'

#plots_path = '../plots/bad_seeing/'
#data_path = '/Users/spb/Work/projects/GalaxyZoo/bad_seeing/'
#data_file = 'dr6z7_zoo_selected_counts_bad_seeing.fits'

#plots_path = '../plots/colour_likelihoods/'
#data_path = '/Users/spb/Work/projects/GalaxyZoo/colour_likelihoods/'
#data_file = 'dr6z7_zoo_selected_counts_colour_likelihoods.fits'

min_ratio = -2.0
max_ratio = 2.0
range_ratio = max_ratio - min_ratio
unknown_ratio = -999
min_ratio_flag = -99
max_ratio_flag = 99
eps = 0.000001
fg_ratio = max_ratio + 0.1*range_ratio
bg_ratio = min_ratio - 0.2*range_ratio
fg_count = 3.5
bg_count = -0.75

oversamp = 10

eg_bins = [(5, 15), (6, 11), (8, 8), (14, 4)]

default_scheme = ('WEIGHTS', 'CSUW')

N.random.seed()

# free whatever memory we can
gc.collect()

def bintype(s):
    if s in ('FUKUGITA', 'MOSES'):
        return s
    else:
        return 'ZOO'


def set_sizes(p=None, s=default_scheme[1], auto=False):
    dataset = bintype(s)
    if dataset == 'ZOO':  dataset = 'UW'
    global min_zbin, max_zbin, min_magbin, max_magbin, min_sizebin, max_sizebin
    if auto:
        # BEWARE WITH AUTO - procedure assumes min is zero 
        min_zbin, max_zbin,\
        min_magbin, max_magbin,\
        min_sizebin, max_sizebin = N.ravel([N.sort(i)[:,[0, -1]] for i
                        in N.where(p['COUNTS_ALL_%s'%dataset].data > 5)])
    else:
        min_zbin = 0
        max_zbin = 24
        if dataset == 'FUKUGITA':
            min_sizebin = 0      # MUST BE ZERO
            max_sizebin = 11
            min_magbin = 0       # MUST BE ZERO
            max_magbin = 13
        else:
            min_sizebin = 0      # MUST BE ZERO
            max_sizebin = 23
            min_magbin = 0       # MUST BE ZERO
            max_magbin = 27
            

def plot_ratio(zbin, min_count=10, scheme=default_scheme,
               scale=True, fiddle=None, twiddle=None):
    f = data_path+data_file
    p = pyfits.open(f)
    d_el = p['%s_EL_%s'%scheme].data
    d_sp = p['%s_SP_%s'%scheme].data
    set_sizes(p, scheme[1])
    zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = p['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = p['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    if fiddle is not None:
        fiddlez = 10**(0.4 * fiddle * (zbins.field('CENTRE') - 0.1))
        d_el_fiddle = N.array([d_el[i] * N.sqrt(f) for i, f in enumerate(fiddlez)], d_el.dtype)
        d_sp_fiddle = N.array([d_sp[i] / N.sqrt(f) for i, f in enumerate(fiddlez)], d_sp.dtype)
        d_tot = d_el + d_sp
        d_fiddle = d_el_fiddle + d_sp_fiddle
        f_fiddle = d_tot / d_fiddle
        d_el_fiddle *= f_fiddle
        d_sp_fiddle *= f_fiddle
        d_el = d_el_fiddle
        d_sp = d_sp_fiddle
    if twiddle is not None and twiddle != 0:
        twiddlez = twiddle * (zbins.field('CENTRE') - 0.1)
        d_el = do_twiddle(d_el, twiddlez)
        d_sp = do_twiddle(d_sp, twiddlez)
    el = d_el[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    sp = d_sp[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    zsel = zbins[zbin].field('CENTRE')
    dmod = cosmology.dmod_flat(zsel)
    ang_scale = cosmology.ang_scale_flat(zsel)
    mlim = 17.77 - dmod
    mask = (sp + el) < min_count
    if scheme[0] == 'COUNTS':
        N.putmask(el, el < 1, -1.0)
        N.putmask(sp, sp < 1, -1.0)
    ratio = el.astype(N.float)/sp
    ratio = N.log10(ratio)
    count = el + sp
    N.putmask(count, count < 1, eps)
    count = N.log10(count)
    if scheme[0] == 'COUNTS':
        N.putmask(ratio, el < 1, min_ratio)
        N.putmask(ratio, sp < 1, max_ratio)
    N.putmask(ratio, mask, unknown_ratio)
    ratio = ratio.transpose()
    count = count.transpose()
    zrange = '%6.4f < z < %6.4f'%(zbins[zbin].field('MIN'),
                                  zbins[zbin].field('MAX'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    SBapp = 23.0
    SBlim_size = N.arange(200) / 10.0
    SBlim_mag = (SBapp - dmod - 2.5*N.log10(6.283185*(SBlim_size/ang_scale)**2))
    # plot ratio
    pgend()
    if scale:
        fname = 'ratio_grid_zbin_%02i.ps/cps'%zbin
    else:
        fname = 'ratio_grid_zbin_%02i_noscale.ps/cps'%zbin
    pgopen(plots_path+'%s_%s/'%scheme+fname)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', '')
    pgimag_s(hires(ratio), fg_ratio, bg_ratio)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
#     for imag, mag in enumerate(magvals):
# 	for isize, size in enumerate(sizevals):
# 	    if ratio[isize, imag] > min_ratio+eps and ratio[isize, imag] < max_ratio-eps:
# 		pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
#                       '%.1f'%(ratio[isize, imag]))
# 		#pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
# 		#       '%.0f %.0f'%(el[imag,isize], sp[imag,isize]))
    pgxsls('dotted')
    pgslw(3*lw)
    pgline(N.array([magbins[min_magbin].field('MIN'),
                    magbins[max_magbin].field('MAX')]),
           N.array([ang_scale, ang_scale]))
    pgline(N.array([mlim, mlim]),
           N.array([sizebins[min_sizebin].field('MIN'),
                    sizebins[max_sizebin].field('MAX')]))
    pgline(SBlim_mag, SBlim_size)
    pgxsls('solid') 
    pgslw(lw)
    pgsch(1.2*ch)
    pgstbg(1)
    pgptxt(-19.0, sizebins[max_sizebin].field('MAX')-1, 0.0, 0.5, zrange)
    pgstbg(-1)
    pgxsci('black')
    pgslw(lw)
    pgsch(ch)
    if scale:
        pgwedg_s(bg_ratio, fg_ratio, 'RI',
                 'log\d10\u(n\del\u/n\dsp\u)', 1.0, 4.0)
        if scheme[0] == 'COUNTS':
            pgmtxt('R', 2.1, 1-0.1/1.3, 0.5, '\(0903)')
            pgxsci('white')
            pgmtxt('R', 2.1, 0.2/1.3, 0.5, '\(0903)')
    pgxsci('black')
    pgsch(ch)
    # plot counts
    pgend()
    if scale:
        fname = 'count_grid_zbin_%02i.ps/cps'%zbin
    else:
        fname = 'count_grid_zbin_%02i_noscale.ps/cps'%zbin
    pgopen(plots_path+'%s_%s/'%scheme+fname)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', '')
    pgimag_s(hires(count), fg_count, bg_count)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
#     for imag, mag in enumerate(magvals):
# 	for isize, size in enumerate(sizevals):
# 	    if count[isize, imag] > 0:
# 		pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
#                       '%.1f'%(count[isize, imag]))
    pgxsls('dotted')
    pgslw(3*lw)
    pgline(N.array([magbins[min_magbin].field('MIN'),
                    magbins[max_magbin].field('MAX')]),
           N.array([ang_scale, ang_scale]))
    pgline(N.array([mlim, mlim]),
           N.array([sizebins[min_sizebin].field('MIN'),
                    sizebins[max_sizebin].field('MAX')]))
    pgline(SBlim_mag, SBlim_size)
    pgxsls('solid') 
    pgslw(lw)
    pgsch(1.2*ch)
    pgstbg(1)
    pgptxt(-19.0, sizebins[max_sizebin].field('MAX')-1, 0.0, 0.5, zrange)
    pgstbg(-1)
    pgxsci('black')
    pgslw(lw)
    pgsch(ch)
    if scale:
        pgwedg_s(bg_count, fg_count, 'RI',
                 'log\d10\u(n\del\u + n\dsp\u)', 1.0, 4.0)
    pgend()
    p.close()


def plot_ratio_quality_compare(zbin, qtype='seeing',
                               min_count=10, scheme=default_scheme,
                               scale=True, fiddle=None, twiddle=None):
    ratioq = {}
    for q in ('bad', 'good'):
        data_path = '/Users/spb/Work/projects/GalaxyZoo/%s_%s/'%(q, qtype)
        data_file = 'dr6z7_zoo_selected_counts_%s_%s.fits'%(q, qtype)
        f = data_path+data_file
        p = pyfits.open(f)
        d_el = p['%s_EL_%s'%scheme].data
        d_sp = p['%s_SP_%s'%scheme].data
        set_sizes(p, scheme[1])
        zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
        magbins = p['MR_%s_BINS'%bintype(scheme[1])].data
        sizebins = p['R50_KPC_%s_BINS'%bintype(scheme[1])].data
        if fiddle is not None:
            fiddlez = 10**(0.4 * fiddle * (zbins.field('CENTRE') - 0.1))
            d_el_fiddle = N.array([d_el[i] * N.sqrt(f) for i, f in enumerate(fiddlez)], d_el.dtype)
            d_sp_fiddle = N.array([d_sp[i] / N.sqrt(f) for i, f in enumerate(fiddlez)], d_sp.dtype)
            d_tot = d_el + d_sp
            d_fiddle = d_el_fiddle + d_sp_fiddle
            f_fiddle = d_tot / d_fiddle
            d_el_fiddle *= f_fiddle
            d_sp_fiddle *= f_fiddle
            d_el = d_el_fiddle
            d_sp = d_sp_fiddle
        if twiddle is not None and twiddle != 0:
            twiddlez = twiddle * (zbins.field('CENTRE') - 0.1)
            d_el = do_twiddle(d_el, twiddlez)
            d_sp = do_twiddle(d_sp, twiddlez)
        el = d_el[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
        sp = d_sp[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
        zsel = zbins[zbin].field('CENTRE')
        dmod = cosmology.dmod_flat(zsel)
        ang_scale = cosmology.ang_scale_flat(zsel)
        mlim = 17.77 - dmod
        mask = (sp + el) < min_count
        N.putmask(el, el < 1, -1.0)
        N.putmask(sp, sp < 1, -1.0)
        ratio = el.astype(N.float)/sp
        ratio = N.log10(ratio)
        N.putmask(ratio, el < 1, min_ratio)
        N.putmask(ratio, sp < 1, max_ratio)
        N.putmask(ratio, mask, unknown_ratio)
        ratioq[q] = ratio.transpose()
    ratioq = ratioq['bad']-ratioq['good']

    zrange = '%6.4f < z < %6.4f'%(zbins[zbin].field('MIN'),
                                  zbins[zbin].field('MAX'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    SBapp = 23.0
    SBlim_size = N.arange(200) / 10.0
    SBlim_mag = (SBapp - dmod - 2.5*N.log10(6.283185*(SBlim_size/ang_scale)**2))
    # plot ratio
    pgend()
    if scale:
        fname = 'ratio_grid_zbin_%02i_%s.ps/cps'%(zbin, qtype)
    else:
        fname = 'ratio_grid_zbin_%02i_%s_noscale.ps/cps'%(zbin, qtype)
    pgopen(plots_path+'%s_%s/'%scheme+fname)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', '')
    pgimag_s(hires(ratioq), 1.0, -1.0)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
#     for imag, mag in enumerate(magvals):
# 	for isize, size in enumerate(sizevals):
# 	    if ratioq[isize, imag] > min_ratio+eps and ratio[isize, imag] < max_ratio-eps:
# 		pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
#                       '%.1f'%(ratioq[isize, imag]))
# 		#pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
# 		#       '%.0f %.0f'%(el[imag,isize], sp[imag,isize]))
    pgxsls('dotted')
    pgslw(3*lw)
    pgline(N.array([magbins[min_magbin].field('MIN'),
                    magbins[max_magbin].field('MAX')]),
           N.array([ang_scale, ang_scale]))
    pgline(N.array([mlim, mlim]),
           N.array([sizebins[min_sizebin].field('MIN'),
                    sizebins[max_sizebin].field('MAX')]))
    pgline(SBlim_mag, SBlim_size)
    pgxsls('solid') 
    pgslw(lw)
    pgsch(1.2*ch)
    pgstbg(1)
    pgptxt(-19.0, sizebins[max_sizebin].field('MAX')-1, 0.0, 0.5, zrange)
    pgstbg(-1)
    pgxsci('black')
    pgslw(lw)
    pgsch(ch)
    if scale:
        pgwedg_s(-1.0, 1.0, 'RI',
                 'log\d10\u(n\del\u/n\dsp\u)', 1.0, 4.0)
    pgxsci('black')
    pgsch(ch)
    pgend()
    p.close()


def plot_ratio_zbins(min_count=10, scheme=default_scheme, fiddle=None, twiddle=None):
    zbins = pyfits.getdata(data_path+data_file,
                           'REDSHIFT_%s_BINS'%bintype(scheme[1]))
    for zbin in zbins.field('BIN'):
	plot_ratio(zbin, min_count, scheme, scale=True, fiddle=fiddle, twiddle=twiddle)
    for zbin in [3, 8]:
	plot_ratio(zbin, min_count, scheme, scale=False, fiddle=fiddle, twiddle=twiddle)


def determine_baseline(min_count=50, scheme=default_scheme, fiddle=None, twiddle=None):
    f = data_path+data_file
    p = pyfits.open(f)
    d_el = p['%s_EL_%s'%scheme].data
    d_sp = p['%s_SP_%s'%scheme].data
    set_sizes(p, scheme[1])
    n_magbin = max_magbin - min_magbin + 1
    n_sizebin = max_sizebin - min_sizebin + 1
    ratio_baseline = N.zeros((n_magbin, n_sizebin), N.float) + unknown_ratio
    counts_baseline = N.zeros((n_magbin, n_sizebin, 2), N.int) + unknown_ratio
    zbins = pyfits.getdata(f, 'REDSHIFT_%s_BINS'%bintype(scheme[1]))
    magbins = p['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = p['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    mag_step = (magbins[0].field('MAX') - magbins[0].field('MIN'))
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    if fiddle is not None:
        fiddlez = 10**(0.4 * fiddle * (zbins.field('CENTRE') - 0.1))
        print 'fiddle:', fiddle
        d_el_fiddle = N.array([d_el[i] * N.sqrt(f) for i, f in enumerate(fiddlez)], d_el.dtype)
        d_sp_fiddle = N.array([d_sp[i] / N.sqrt(f) for i, f in enumerate(fiddlez)], d_sp.dtype)
        d_tot = d_el + d_sp
        d_fiddle = d_el_fiddle + d_sp_fiddle
        f_fiddle = d_tot / d_fiddle
        N.putmask(f_fiddle, d_fiddle < 1.0, 0.0)
        d_el_fiddle *= f_fiddle
        d_sp_fiddle *= f_fiddle        
        rbefore = d_el.sum(1).sum(1)/d_sp.sum(1).sum(1)
        d_el = d_el_fiddle
        d_sp = d_sp_fiddle
        rafter = d_el.sum(1).sum(1)/d_sp.sum(1).sum(1)
        print 'Redshift bin centres:', zbins.field('CENTRE')
        print 'Overall el/sp ratio before fiddle:', rbefore
        print 'Overall el/sp ratio after fiddle:', rafter
    if twiddle is not None and twiddle != 0:
        rbefore = d_el.sum(1).sum(1)/d_sp.sum(1).sum(1)
        twiddlez = twiddle * (zbins.field('CENTRE') - 0.1)
        d_el = do_twiddle(d_el, twiddlez)
        d_sp = do_twiddle(d_sp, twiddlez)
        rafter = d_el.sum(1).sum(1)/d_sp.sum(1).sum(1)
        print 'Redshift bin centres:', zbins.field('CENTRE')
        print 'mag bin centres:', magbins.field('CENTRE')
        print 'size bin centres:', sizebins.field('CENTRE')
        print 'Overall el/sp ratio before twiddle:', rbefore
        print 'Overall el/sp ratio after twiddle:', rafter
    elz = d_el[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    spz = d_sp[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    for iz in zbins.field('BIN'):
        zsel = zbins.field('CENTRE')[iz]
	el = elz[iz]
	sp = spz[iz]
        # ignore cells with low counts
        mask = (sp + el) < min_count
        # ignore cells too close to selection limits
        dmod = cosmology.dmod_flat(zsel)
        ang_scale = cosmology.ang_scale_flat(zsel)
        mlim = 17.77 - dmod
        SBapp = 23.0
        for i in range(mask.shape[0]):
            for j in range(mask.shape[1]):
                mask[i,j] = mask[i,j] | (magvals[i] > (mlim - 1.0))
                #mask[i,j] = mask[i,j] | (magvals[i] > -19.0)
                mask[i,j] = mask[i,j] | (sizevals[j] < (2.0*ang_scale))
                #SBlim_mag = (SBapp - dmod - 2.5*N.log10(6.283185*(sizevals[j]/ang_scale)**2))
                SB = (magvals[i] + dmod + 
                      2.5*N.log10(6.283185*(sizevals[j]/ang_scale)**2))
                mask[i,j] = mask[i,j] | (SB > (SBapp - 1.0))
        if scheme[0] == 'COUNTS':
            # ell/sp ratio - mask zero counts for division
            N.putmask(el, el < 1, -1.0)
            N.putmask(sp, sp < 1, -1.0)
        ratio = el.astype(N.float)/sp
        ratio = N.log10(ratio)
        if scheme[0] == 'COUNTS':
            N.putmask(ratio, el < 1, min_ratio_flag)
            N.putmask(ratio, sp < 1, max_ratio_flag)
        N.putmask(ratio, mask, unknown_ratio)
	select = N.logical_not(mask)
        if scheme[0] == 'COUNTS':
            # consider cells with no ell or sp as empty
            empty = ((ratio_baseline <= min_ratio_flag+eps) |
                     (ratio_baseline >= max_ratio_flag-eps))
        else:
            empty = ratio_baseline <= unknown_ratio+eps
        select &= empty
	ratio_baseline[select] = ratio[select]
	counts_baseline[select] = N.transpose([el[select], sp[select]])
    pyfits.writeto(data_path+'%s_%s/'%scheme+'ratio_baseline.fits',
		   ratio_baseline, clobber=True)    
    pyfits.writeto(data_path+'%s_%s/'%scheme+'counts_baseline.fits',
		   counts_baseline, clobber=True)    
    p.close()


def plot_ratio_baseline(type='data', scheme=default_scheme):
    f = data_path+data_file
    p = pyfits.open(f)
    set_sizes(p, scheme[1])
    zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = p['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = p['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    p.close()
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    if bintype(scheme[1]) == 'MOSES':
        ratio_zoo = pyfits.getdata(data_path+'COUNTS_CSDW/'+'ratio_baseline.fits')
    if type=='fit' or type=='fitsmooth':
	ratio = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.fits')
    elif type=='fitres':
	ratio = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_fitres.fits')
    elif type=='fitnormres':
	ratio = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_fitnormres.fits')
    elif type=='sigma':
	ratio = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_sigma.fits')
    elif type=='filled':
	ratio = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_filled.fits')
    else:
	ratio = data
    #ratio = ratio[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    #ratio_zoo = ratio_zoo[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_sp = ratio >= max_ratio_flag-eps
    mask_ell = N.logical_and(ratio <= min_ratio_flag+eps,
                             ratio > unknown_ratio+eps)
    mask_ok = N.logical_and(data > min_ratio_flag+eps,
                            data < max_ratio_flag-eps)
    if bintype(scheme[1]) == 'MOSES':
        mask_zoo = ((ratio_zoo > min_ratio_flag+eps) &
                    (ratio_zoo < max_ratio_flag-eps))
    #N.putmask(ratio, mask_sp, max_ratio)
    #N.putmask(ratio, mask_ell, min_ratio)
    if type not in ('fit'):
        N.putmask(ratio, N.logical_not(mask_ok), min_ratio_flag)
    ratio = ratio.transpose()
    # selection limits at z=0.1
    zsel = 0.085
    dmod = cosmology.dmod_flat(zsel)
    ang_scale = cosmology.ang_scale_flat(zsel)
    mlim = 17.77 - dmod
    SBapp = 23.0
    SBlim_size = N.arange(200) / 10.0
    SBlim_mag = (SBapp - dmod - 2.5*N.log10(6.283185*(SBlim_size/ang_scale)**2))
    SBappfid = 20.25
    SBfid_mag = (SBappfid - dmod - 2.5*N.log10(6.283185*(SBlim_size/ang_scale)**2))
    pgend()
    fg = fg_ratio
    bg = bg_ratio
    if type=='fit':
	pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid_fit.ps/cps')
	lab = 'n\del\u/n\dsp\u baseline fit'
    elif type=='fitres':
	pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid_fitres.ps/cps')
	lab = 'n\del\u/n\dsp\u baseline fit residuals'
	fg = 1.0
	bg = -1.0
    elif type=='fitnormres':
	pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid_fitnormres.ps/cps')
	lab = 'n\del\u/n\dsp\u baseline fit residuals/sigma'
	fg = 5.0
	bg = -5.0
    elif type=='sigma':
	pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid_sigma.ps/cps')
	lab = 'n\del\u/n\dsp\u baseline sigma'
	fg = 0.5
	bg = -0.05
    elif type=='filled':
	pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid_filled.ps/cps')
	lab = 'n\del\u/n\dsp\u baseline filled'
    else:
	pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid.ps/cps')
	lab = 'n\del\u/n\dsp\u baseline'
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', '')
    pgimag_s(hires(ratio), fg, bg)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
#     for imag, mag in enumerate(magvals):
# 	for isize, size in enumerate(sizevals):
# 	    if (ratio[isize, imag] > min_ratio+eps 
#                  and ratio[isize, imag] < max_ratio-eps):
# 		pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
# 		       '%.1f'%(ratio[isize, imag]))
    pgslw(3*lw)
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
    if bintype(scheme[1]) == 'MOSES':
        pgxsci('cyan')
        mask_zoo = -N.transpose(mask_zoo).astype(N.float)
        pgcont_s(mask_zoo.astype(N.float), 1, N.array([-0.5]))
    pgxsci('white')
    pgxsls('dotted')
    pgslw(3*lw)
    pgline(N.array([magbins[min_magbin].field('MIN'),
                    magbins[max_magbin].field('MAX')]),
           N.array([ang_scale, ang_scale]))
    pgline(N.array([mlim, mlim]),
           N.array([sizebins[min_sizebin].field('MIN'),
                    sizebins[max_sizebin].field('MAX')]))
    pgline(SBlim_mag, SBlim_size)
    pgxsls('dashed')
    pgslw(3*lw/2)    
    pgline(SBfid_mag, SBlim_size)
    pgxsls('solid')
    pgsch(ch)
    if type == 'fit':
        egx = N.array([magvals[i[0]] for i in eg_bins])
        egy = N.array([sizevals[i[1]] for i in eg_bins])
        pgxpt(egx, egy, 'dot')
    pgxsci('black')
    pgslw(lw)
    pgwedg_s(bg, fg, 'RI', lab, 1.0, 4.0)
    pgend()


def plot_ratio_baseline_fit(nover=5, scheme=default_scheme):
    f = data_path+data_file
    p = pyfits.open(f)
    set_sizes(p, scheme[1])
    magbins = p['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = p['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    mag_step = (magbins[0].field('MAX') - magbins[0].field('MIN'))
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    p.close()
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = N.logical_and(data > min_ratio_flag+eps,
                            data < max_ratio_flag-eps)
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    magvals_over = N.arange(magvals[0] - (1-1.0/nover)*mag_step/2.0,
			    magvals[-1] + mag_step - (1-1.0/nover)*mag_step/2.0,
			    mag_step/nover)
    sizevals_over = N.arange(sizevals[0] - (1-1.0/nover)*size_step/2.0,
			     sizevals[-1] + size_step - (1-1.0/nover)*size_step/2.0,
			     size_step/nover)
    p, perr = pickle.load(file(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.pickle'))
    f = ratio_function(p, magvals_over, sizevals_over)
    f = f.transpose()
    pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid_fitsmooth.ps/cps')
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', 'n\del\u/n\dsp\u baseline fit smooth')
    fg = fg_ratio
    bg = bg_ratio
    pgimag_s(f, fg, bg)
    pgxsci('white')
    pgsch(0.3*ch/nover)
    pgslw(lw)
    if nover == 1:
        for imag, mag in enumerate(magvals_over):
            for isize, size in enumerate(sizevals_over):
                pgptxt(mag, size-0.15*size_step/nover, 0.0, 0.5,
                       '%.1f'%(f[isize, imag]))
    pgslw(2*lw)
    pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
    pgend()
    ferr = ratio_function_err(p, perr, magvals_over, sizevals_over)
    ferr = ferr.transpose()
    pgopen(plots_path+'%s_%s/'%scheme+'ratio_baseline_grid_sigma_fitsmooth.ps/cps')
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', 'n\del\u/n\dsp\u baseline fit sigma smooth')
    fg = fg_ratio
    bg = bg_ratio
    pgimag_s(ferr, fg, bg)
    pgxsci('white')
    pgsch(0.3*ch/nover)
    pgslw(lw)
    if nover == 1:
        for imag, mag in enumerate(magvals_over):
            for isize, size in enumerate(sizevals_over):
                pgptxt(mag, size-0.15*size_step/nover, 0.0, 0.5,
                       '%.1f'%(ferr[isize, imag]))
    pgslw(2*lw)
    pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
    pgend()


def fit_ratio_baseline(scheme=default_scheme):
    f = data_path+data_file
    pf = pyfits.open(f)
    set_sizes(pf, scheme[1])
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    sigma = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_sigma.fits')
    ratio = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #pinit = N.array((8.0, 2.0, -23.0, 0.4, 0.4, -4.5, 3.0))
    #pinit = N.array([6.0, 1.5, -22.8, 0.2, 0.2, -3.0, 1.25])
    pinit = N.array([-0.3, 1.7, -22.8, 0.3, 0.3, -2.6, 1.5, 1.0, 1.0])
    chi2r = 0.0
    niter = 0
#     sys = sigma[sigma > 0.0].mean()/10.0
#     oldchi2r = -100
#     while abs(chi2r - 1.0) >  0.001 and abs(chi2r - oldchi2r) > 0.01:
#         print 'niter =', niter
#         print 'sys =', sys
#         niter += 1
#         weights = 1.0/(sigma**2 + sys**2)
#         #N.putmask(weights, weights > 0.001, 0.1)
#         #fout = file('../xyz', 'w')
#         #for imag, mag in enumerate(magvals):
#         #    for isize, size in enumerate(sizevals):
#         #        if weights[imag, isize] > 0.00001:
#         #            fout.write('%f %f %f\n'%(mag, size, ratio[imag, isize]))
#         #fout.close()
#         p, perr = bootstrap_fmin(ratio_minfunc, pinit,
#                                  magvals, sizevals, ratio, weights, 10)
#         oldchi2r = chi2r
#         chi2r = ratio_minfunc(p, magvals, sizevals, ratio, weights)
#         sys = sys * chi2r**(2/3.)
#         print 'param =', p
#         print 'param errors =', perr
#         print 'chi2r =', chi2r
    sys = 0.0
    weights = 1.0/(sigma**2 + sys**2)
    # unweighted
    weights = weights.ravel()
    weights = 0.0 * sigma + weights.compress(weights > 0.000001).mean()
    N.putmask(weights, sigma > 98, 0.0)
    N.putmask(weights, sigma < -98, 0.0)
    p, perr = bootstrap_fmin(ratio_minfunc, pinit,
                             magvals, sizevals, ratio, weights)
    #oldchi2r = chi2r
    chi2r = ratio_minfunc(p, magvals, sizevals, ratio, weights)
    #sys = sys * chi2r**(2/3.)
    print 'param =', p
    print 'param errors =', perr
    print 'chi2r =', chi2r
#     print 'Anneal:'
#     lower = p-N.absolute(p)*0.1
#     upper = p+N.absolute(p)*0.1
#     p, retval = anneal(ratio_minfunc, p,
#                        (magvals, sizevals, ratio, weights),
#                        schedule = 'boltzmann',
#                        lower=lower, upper=upper,
#                        maxiter=100)
#     print 'param =', p
#     print 'reval =', retval
#     s = ratio_minfunc(p, magvals, sizevals, ratio, weights)
#     print 'chi2 = %4.3f'%s
    ratio_baseline_fit = ratio_function(p, magvals, sizevals)
    pickle.dump((p, perr), file(data_path+'%s_%s/'%scheme+
                                'ratio_baseline_fit.pickle', 'w'))
    res = ratio-ratio_baseline_fit
    normres = res/sigma
    N.putmask(normres, weights < 0.000001, unknown_ratio)
    pyfits.writeto(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.fits',
		   ratio_baseline_fit, clobber=True)    
    pyfits.writeto(data_path+'%s_%s/'%scheme+'ratio_baseline_fitres.fits',
		   res, clobber=True) 
    pyfits.writeto(data_path+'%s_%s/'%scheme+'ratio_baseline_fitnormres.fits',
		   normres, clobber=True) 
    pf.close()


def determine_ratio_baseline_sigma(scheme=default_scheme):
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'counts_baseline.fits')
    el = data[:,:,0].astype(N.float)
    sp = data[:,:,1].astype(N.float)
    mask_el = el < 1
    mask_sp = sp < 1
    mask_all = N.logical_and(mask_el, mask_sp)
    mask = N.logical_or(mask_el, mask_sp)
    ratio = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    sum = (el + sp).astype(N.float)
    product = (el * sp).astype(N.float)
    N.putmask(product, mask, 1.0)
    N.putmask(sum, mask, 1.0)
    sigma = (log10(e))**2 * sum / product
    sigma = N.sqrt(sigma)
    N.putmask(sigma, mask, unknown_ratio)
    pyfits.writeto(data_path+'%s_%s/'%scheme+'ratio_baseline_sigma.fits',
		   sigma, clobber=True)    


def ratio_intfunc(y, x, p):
    f = ratio_function(p, N.array([x]), N.array([y]))[0]
    return f[0]

def ratio_function2(p, x, y):
    #p = (3.5, -21, -1.5, 0.05, -0.2, -18, 3.0, -2.0)
    x0, x1, x2, x3, x4, x5, y2, y3 = p
    y0 = x0*N.exp((x-x1)/x2)
    y1 = x3*(x - x5)**2 + x4*(x - x5)
    z = N.zeros((len(x), len(y)), N.float)
    for i in range(len(x)):
	z[i] = y2/(1.0 + N.exp((y-y0[i])/y1[i])) + y3
    return z

def ratio_function3(p, x, y):
    #p = (0.1, 2, 2, 0.001, -0.01, 0.5, 3.0, -2.0)
    x0, x1, x2, x3, x4, x5, y2, y3 = p
    y0 = x0*x**2 + x1*x + x2
    y1 = x3*x**2 + x4*x + x5
    z = N.zeros((len(x), len(y)), N.float)
    for i in range(len(x)):
	z[i] = y2/(1.0 + N.exp((y-y0[i])/y1[i])) + y3
    return z

def ratio_function4(p, x, y):
    #p = (2.0, 14.8, 0.5, 4.0, 3.0, -2.0)
    x0, x1, x2, x3, y2, y3 = p
    y0 = x0/20**x1 * (-x)**x1
    y1 = x2/20**x3 * (-x)**x3
    z = N.zeros((len(x), len(y)), N.float)
    for i in range(len(x)):
	z[i] = y2/(1.0 + N.exp((y-y0[i])/y1[i])) + y3
#     print x
#     print y
#     print p
#     print y0, y1
#     print z
    return z


def ratio_function5(p, x, y):
    #p = (-1.0E+04, -2.2E+01, 6.1E-01, -7.3E+00, 3.3E-01, 1.0E+04)
    a, b, c, d, e, f = p
    z = N.zeros((len(x), len(y)), N.float)
    for i in range(len(x)):
	z[i] = a / ((1.0 + N.exp(b - c*x[i])) * (1.0 + N.exp(d - e*y))) + f
    return z


def ratio_function6(p, x, y):
    #p = (7.0, 1.5, -23.5, 0.5, 0.3, -3.5, 1.5)
    a, b, c, d, e, f, g = p
    z = N.zeros((len(x), len(y)), N.float)
    x0 = a*(b**-y) + c
    x1 = d + e*(x0-c)
    #x1 = d + e**(x0-c)
    for i in range(len(x)):
	z[i] = f / (1.0 + N.exp((x0 - x[i])/x1)) + g
    return z

def ratio_function(p, x, y):
    a, b, c, d, e, f, g, h, i = p
    z = N.zeros((len(x), len(y)), N.float)
    x0 = b**-(a + h*y**i) + c
    x1 = d + e*(x0-c)
    for i in range(len(x)):
	z[i] = f / (1.0 + N.exp((x0 - x[i])/x1)) + g
    return z

def ratio_function_err(p, perr, x, y, nmc=100):
    ffit = ratio_function(p, x, y)
    mc = N.zeros(ffit.shape+(nmc,), N.float)
    for i in range(nmc):
        pmc = p + N.random.normal(0.0, 1.0, len(perr)) * perr
        mc[:,:,i] = ratio_function(pmc, x, y)
    sigma = mc.std(axis=-1)
#     sigma = N.array([[(scoreatpercentile(mc[i,j,:], 84.0) -
#                       scoreatpercentile(mc[i,j,:], 16.0))/2.0
#                      for j in range(ffit.shape[1])]
#                     for i in range(ffit.shape[0])])
    return sigma

def ratio_minfunc(p, x, y, z, w):
    f = ratio_function(p, x, y)
    r2 = (z - f)**2
    # prefer positive residuals (function is more of a lower limit)
    #N.putmask(r2, f > z, 2*r2)
    df = (w > 0).ravel().astype(N.int).sum() - len(p)
    s = (w * r2).sum() / df
    return s

def ratio_binfunc(p, x, y, xstep, ystep):
    f = N.zeros((len(x), len(y)))
    for j in range(len(y)):
	def gfun(x):  return y[j]-0.5*ystep
	def hfun(x):  return y[j]+0.5*ystep
	for i in range(len(x)):
	    print i, j
	    f[i,j] = dblquad(ratio_intfunc,
			     x[i]-0.5*xstep, x[i]+0.5*xstep,
			     gfun, hfun, (p,), 0.001, 0.01)[0]
    return f


def determine_baseline_correction(min_count=10, scheme=default_scheme,
                                  nmc=100, fromfit=True, fill=True,
                                  forcezero=True, fiddle=None, twiddle=None):
    pf = pyfits.open(data_path+data_file)
    if fromfit:
        p, perr = pickle.load(file(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.pickle'))
    set_sizes(pf, scheme[1])
    d_el = pf['%s_EL_%s'%scheme].data
    d_sp = pf['%s_SP_%s'%scheme].data
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    if fiddle is not None:
        fiddlez = 10**(0.4 * fiddle * (zbins.field('CENTRE') - 0.1))
        d_el_fiddle = N.array([d_el[i] * N.sqrt(f) for i, f in enumerate(fiddlez)], d_el.dtype)
        d_sp_fiddle = N.array([d_sp[i] / N.sqrt(f) for i, f in enumerate(fiddlez)], d_sp.dtype)
        d_tot = d_el + d_sp
        d_fiddle = d_el_fiddle + d_sp_fiddle
        f_fiddle = d_tot / d_fiddle
        d_el_fiddle *= f_fiddle
        d_sp_fiddle *= f_fiddle
        d_el = d_el_fiddle
        d_sp = d_sp_fiddle
    if twiddle is not None and twiddle != 0:
        twiddlez = twiddle * (zbins.field('CENTRE') - 0.1)
        d_el = do_twiddle(d_el, twiddlez)
        d_sp = do_twiddle(d_sp, twiddlez)
    elz = d_el[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    spz = d_sp[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    if fromfit:
        ffit = ratio_function(p, magvals, sizevals)
        ffit_sigma = ratio_function_err(p, perr, magvals, sizevals, nmc)
    else:
        ffit = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
        ffit_sigma = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_sigma.fits')
    #ffit_sigma = ffit * 0.0
    correction = N.zeros(elz.shape[:3], N.float)
    correction_sigma = N.zeros(elz.shape[:3], N.float)
    for iz, z in enumerate(zbins):
	el = elz[iz]
	sp = spz[iz]
        mask_baseline = ffit < min_ratio_flag + eps
        mask_min = ((el + sp) < min_count) | mask_baseline
        if scheme[0] == 'COUNTS':
            mask_el = el < 1
            mask_sp = sp < 1
            mask = mask_min | mask_sp | mask_el
            N.putmask(sp, mask_sp, 1.0)
            N.putmask(el, mask_el, 1.0)
        else:
            mask = mask_min
	ratio = el.astype(N.float)/sp
        sum = (el + sp).astype(N.float)
        product = (el * sp).astype(N.float)
        N.putmask(product, mask, 1.0)
        N.putmask(sum, mask, 1.0)
        sigma = (log10(e))**2 * sum / product
        sigma = N.sqrt(sigma)
        N.putmask(sigma, mask, unknown_ratio)
	ratio = N.log10(ratio)
        if scheme[0] == 'COUNTS':
            N.putmask(ratio, mask_sp, max_ratio_flag)
            N.putmask(ratio, mask_el, min_ratio_flag)
            N.putmask(ratio, N.logical_and(mask_sp, mask_el), unknown_ratio)
        else:
            N.putmask(ratio, mask, unknown_ratio)
	correction[iz] = ratio - ffit
        #correction_sigma[iz] = N.sqrt(sigma**2 + ffit_sigma**2)
        correction_sigma[iz] = sigma  # assuming baseline "truth"
        if scheme[0] == 'COUNTS':
            N.putmask(correction[iz], mask_el, 0.0)
            N.putmask(correction[iz], mask_sp, 0.0)
	N.putmask(correction[iz], mask_min, -99999)
        N.putmask(correction_sigma[iz], mask, -99999)
        #if iz == 5:
        #    print ratio
        #    print ffit
        #    print correction[iz]
        #if N.any(mask != True):
        #    cx = correction[iz].ravel()
        #    mx = N.logical_not(mask.ravel())
        #    print '%5.3f: %5.3f'%(z.field('CENTRE'), cx[mx].mean())
        #    #print cx[mx]
        if forcezero:
            N.putmask(correction[iz],
                      (correction[iz] < 0.0) & (correction[iz] > -9999), 0.0)
        if fill:
            correction_filled = correction[iz] * 0.0
            nx = len(magvals)
            ny = len(sizevals)
            xlist = N.arange(nx)
            ylist = N.arange(ny)
            xlist = xlist.reshape((len(xlist), 1))
            masked = ((correction[iz] > (max_ratio_flag-eps)) |
                      (correction[iz] < (min_ratio_flag+eps)))
            #print 'nx, ny =', nx, ny
            for xi in range(nx):
                for yi in range(ny):
                    if masked[xi,yi]:
                        #print 'xi, yi =', xi, yi
                        dx = xlist - xi
                        dy = ylist - yi
                        d = N.sqrt(dx**2 + dy**2)
                        N.putmask(d, masked, nx*ny)
                        d = d.ravel()
                        dmin = d.min()
                        di = N.argsort(d)
                        dsorted = d[di]
                        di = di[dsorted < 1.2*dmin]
                        #print dmin, len(di)
                        xdi = di / ny
                        ydi = di - xdi*ny
                        #print xdi, ydi
                        #print correction[iz, xdi, ydi]
                        correction_filled[xi, yi] = correction[iz, xdi, ydi].mean()
                        #print correction_filled[xi, yi]
                    else:
                        correction_filled[xi, yi] = correction[iz, xi, yi]
            correction[iz] = correction_filled
    pyfits.writeto(data_path+'%s_%s/'%scheme+'baseline_correction.fits',
		   correction, clobber=True)
    pyfits.writeto(data_path+'%s_%s/'%scheme+'baseline_correction_sigma.fits',
		   correction_sigma, clobber=True)


def plot_baseline_correction(zbin, scheme=default_scheme,
                             scale=True, min_count=10,
                             fromfit=True, fiddle=None, twiddle=None):
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    zsel = zbins[zbin].field('CENTRE')
    zrange = '%6.4f < z < %6.4f'%(zbins[zbin].field('MIN'),
                                  zbins[zbin].field('MAX'))
    d = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction.fits')
    dsel = d[zbin]
    #dsel = dsel[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    dsel = dsel.transpose()
    if fromfit:
        data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.fits')
    else:
        data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    f = data_path+data_file
    p = pyfits.open(f)
    d_el = p['%s_EL_%s'%scheme].data
    d_sp = p['%s_SP_%s'%scheme].data
    set_sizes(p, scheme[1])
    if fiddle is not None:
        zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
        fiddlez = 10**(0.4 * fiddle * (zbins.field('CENTRE') - 0.1))
        d_el_fiddle = N.array([d_el[i] * N.sqrt(f) for i, f in enumerate(fiddlez)], d_el.dtype)
        d_sp_fiddle = N.array([d_sp[i] / N.sqrt(f) for i, f in enumerate(fiddlez)], d_sp.dtype)
        d_tot = d_el + d_sp
        d_fiddle = d_el_fiddle + d_sp_fiddle
        f_fiddle = d_tot / d_fiddle
        d_el_fiddle *= f_fiddle
        d_sp_fiddle *= f_fiddle
        d_el = d_el_fiddle
        d_sp = d_sp_fiddle
    if twiddle is not None and twiddle != 0:
        twiddlez = twiddle * (zbins.field('CENTRE') - 0.1)
        d_el = do_twiddle(d_el, twiddlez)
        d_sp = do_twiddle(d_sp, twiddlez)
    elz = d_el[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    spz = d_sp[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = ((data > min_ratio_flag+eps) &
               (data < max_ratio_flag-eps) &
               (elz+spz >= min_count))
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    mask_view = N.transpose(elz+spz >= 1.0)
    binstruct = generate_binary_structure(2, 2)
    mask_view = binary_opening(mask_view, binstruct)
    mask_view = N.logical_not(mask_view)
    N.putmask(dsel, mask_view, -999)
    pgend()
    if scale:
        fname = 'baseline_correction_zbin_%02i.ps/cps'%zbin
    else:
        fname = 'baseline_correction_zbin_%02i_noscale.ps/cps'%zbin
    pgopen(plots_path+'%s_%s/'%scheme+fname)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', ' ')
    pgimag_s(hires(dsel), 2.0, -1.0)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
#     for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
# 	for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
# 	    if dsel[isize, imag] > -98 and dsel[isize, imag] < 98:
#                 pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
#                        '%.1f'%(dsel[isize, imag]))
    pgslw(3*lw)
    pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
    dmod = cosmology.dmod_flat(zsel)
    mlim = 17.77 - dmod
    pgxsls('dotted')  
    pgline(N.array([mlim, mlim]),
           N.array([sizebins[min_sizebin].field('MIN'),
                    sizebins[max_sizebin].field('MAX')]))
    pgxsls('solid')  
    pgslw(lw)
    pgsch(1.2*ch)
    #pgstbg(1)
    pgptxt(-19.0, sizebins[max_sizebin].field('MAX')-1, 0.0, 0.5, zrange)
    #pgstbg(-1)
    pgxsci('black')
    pgslw(lw)
    pgsch(ch)
    if scale:
        pgwedg_s(-1.0, 2.0, 'RI',
                  'log\d10\u(n\del\u/n\dsp\u) adjustment', 1.0, 4.0)
    pgend()

def plot_baseline_correction_sigma(zbin, scheme=default_scheme,
                                   min_count=10, scale=True,
                                   fromfit=True, fiddle=None, twiddle=None):
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    zsel = zbins[zbin].field('CENTRE')
    zrange = '%6.4f < z < %6.4f'%(zbins[zbin].field('MIN'),
                                  zbins[zbin].field('MAX'))
    d = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction_sigma.fits')
    dsel = d[zbin]
    #dsel = dsel[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    dsel = dsel.transpose()
    if fromfit:
        data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.fits')
    else:
        data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    f = data_path+data_file
    p = pyfits.open(f)
    d_el = p['%s_EL_%s'%scheme].data
    d_sp = p['%s_SP_%s'%scheme].data
    set_sizes(p, scheme[1])
    if fiddle is not None:
        zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
        fiddlez = 10**(0.4 * fiddle * (zbins.field('CENTRE') - 0.1))
        d_el_fiddle = N.array([d_el[i] * N.sqrt(f) for i, f in enumerate(fiddlez)], d_el.dtype)
        d_sp_fiddle = N.array([d_sp[i] / N.sqrt(f) for i, f in enumerate(fiddlez)], d_sp.dtype)
        d_tot = d_el + d_sp
        d_fiddle = d_el_fiddle + d_sp_fiddle
        f_fiddle = d_tot / d_fiddle
        d_el_fiddle *= f_fiddle
        d_sp_fiddle *= f_fiddle
        d_el = d_el_fiddle
        d_sp = d_sp_fiddle
    if twiddle is not None and twiddle != 0:
        twiddlez = twiddle * (zbins.field('CENTRE') - 0.1)
        d_el = do_twiddle(d_el, twiddlez)
        d_sp = do_twiddle(d_sp, twiddlez)
    elz = d_el[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    spz = d_sp[zbin, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = ((data > min_ratio_flag+eps) &
               (data < max_ratio_flag-eps) &
               (elz+spz >= min_count))
    pgend()
    if scale:
        fname = 'baseline_correction_sigma_zbin_%02i.ps/cps'%zbin
    else:
        fname = 'baseline_correction_sigma_zbin_%02i_noscale.ps/cps'%zbin
    pgopen(plots_path+'%s_%s/'%scheme+fname)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    pglab('M\dr\u', 'R\d50\u (kpc)', 'n\del\u/n\dsp\u correction '+zrange)
    pgimag_s(hires(dsel), 1.5, -0.5)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
#     for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
# 	for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
# 	    if dsel[isize, imag] > -98 and dsel[isize, imag] < 98:
#                 pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
#                        '%.1f'%(dsel[isize, imag]))
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    pgslw(2*lw)
    pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
    pgxsls('solid') 
    pgxsci('black')
    pgslw(lw)
    pgsch(ch)
    if scale:
        pgwedg_s(-0.5, 1.5, 'RI',
                  'log\d10\u(n\del\u/n\dsp\u) adjustment uncertainty', 1.0, 4.0)
    pgend()


def plot_baseline_correction_zbins(scheme=default_scheme, min_count=10,
                                   fromfit=True, fiddle=None, twiddle=None):
    zbins = pyfits.getdata(data_path+data_file,
                           'REDSHIFT_%s_BINS'%bintype(scheme[1]))
    for zbin in zbins.field('BIN'):
	plot_baseline_correction(zbin, scheme=scheme, min_count=min_count,
                                 fromfit=fromfit, fiddle=fiddle, twiddle=twiddle)
        plot_baseline_correction_sigma(zbin, scheme=scheme, min_count=min_count,
                                       fromfit=fromfit, fiddle=fiddle, twiddle=twiddle)
    for zbin in [3, 8]:
	plot_baseline_correction(zbin, scheme=scheme, min_count=min_count,
                                 scale=False, fromfit=fromfit, fiddle=fiddle, twiddle=twiddle)
        plot_baseline_correction_sigma(zbin, scheme=scheme, min_count=min_count,
                                       scale=False, fromfit=fromfit, fiddle=fiddle, twiddle=twiddle)


def fit_baseline_correction_zfunc(scheme=default_scheme):
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    d = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction.fits')
    dsel = d
    #dsel = d[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    derr = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction_sigma.fits')
    derrsel = derr
    #derrsel = derr[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = N.logical_and(data > min_ratio_flag+eps,
                            data < max_ratio_flag-eps)
    pgend()
    pgopen(plots_path+'%s_%s/'%scheme+'baseline_correction_zfunc.ps/cps')
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    if (scheme[0] =='WEIGHTS'):
        corrmax = 0.999999
    else:
        corrmax = 3.0
    pgenv(zbins[0].field('MIN'), zbins[-1].field('MAX'), -0.3, corrmax)
    pglab('z', 'log\d10\u(n\del\u/n\dsp\u) adjustment', '')
    correction_zfit = N.zeros(data.shape + (5,), N.float) - 99999
    N.random.seed(12345)
    i = 36
    for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
	i = 36+8*imag
	for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
            #r = N.random.uniform(0.0, 1.0)
            #if r > 1.0/5:
            #    continue
	    if mask_ok[imag, isize]:
                if (imag, isize) not in eg_bins:
                    continue
                print imag, isize, mag, size
		zselect = N.where((dsel[:, imag, isize] > -98) & 
                                  (derr[:, imag, isize] > -98))[0]
		if (len(zselect) > 1):
		    #pgsci(i)
		    zmin = zselect.min()
		    zmax = zselect.max() + 1
                    zb = zbins[zselect].field('CENTRE')
                    mean_z = zb.mean()
                    zbfit = zb #-mean_z
#  		    fit = fitmc.fitmc(zbfit, zbfit*0.0,
#                                       dsel[zselect, imag, isize],
#                                       derrsel[zselect, imag, isize])
#                     a, siga, b, sigb, corab, afit, bfit = fit
#  		    correction_zfit[imag, isize] = [a, b, siga, sigb, corab]
#                     pgxsls('solid') 
                    #pgline(zbins[zselect].field('CENTRE'),
                    #       N.asarray(a+b*zbfit))
                    #pgxsls('dashed')
                    #for j in range(20):
                    #    ra = N.random.normal(0.0, 1.0)
                    #    r = N.random.normal(0.0, 1.0)
                    #    rb = corab*ra + (1-corab)*r
                    #    pgline(zb, N.asarray((a+siga*ra) + (b+sigb*rb)*zbfit))
                    #    #pgline(zb, N.asarray((afit[j])+ (bfit[j])*zbfit))
                    pgxsls('solid') 
                    pgerrby(zb,
                            dsel[zselect, imag, isize],
                            derrsel[zselect, imag, isize])
                    pgslw(2*lw)
                    pgline(zb,
                           dsel[zselect, imag, isize])
                    pgslw(lw)
                    pgsch(2.0*ch)
                    for j in range(len(zb)):
                        corr = dsel[zselect, imag, isize][j:j+1]
                        col = int(32 + 128 * (corr + 1.0) / 3.0)
                        pgsci(col)
                        pgxpt(zb[j:j+1], corr, 'dot')
                    pgsch(ch)
                    pgxsci('black')
    pgxsci('black')
    pgend()
    pyfits.writeto(data_path+'%s_%s/'%scheme+'baseline_correction_zfit.fits',
		   correction_zfit, clobber=True)

def plot_baseline_correction_zfit(scheme=default_scheme):
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = N.logical_and(data > -98, data < 98)
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    zfit = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction_zfit.fits')
    #zfit = zfit[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    a = zfit[:,:,0]
    b = zfit[:,:,1]
    a = a.transpose()
    b = b.transpose()
    for f, l, bg, fg in ((a, 'int', -2.0, 1.0), (b, 'slope', 0.0, 50.0)):
	fout = file('../z_%s'%l, 'w')
	lab = 'n\del\u/n\dsp\u baseline correction zfit '+l
	pgend()
	pgopen(plots_path+'%s_%s/'%scheme+'baseline_correction_zfit_%s.ps/cps'%l)
	pgsetup()
	pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
	setup_colour_table('smooth2', 'ramp', flip_colours=False)
	pgsch(1.2*ch)
        pgenv(magbins[min_magbin].field('MIN'),
              magbins[max_magbin].field('MAX'),
              sizebins[min_sizebin].field('MIN'),
              sizebins[max_sizebin].field('MAX'))
	pglab('M\dr\u', 'R\d50\u (kpc)', lab)
	pgimag_s(f, fg, bg)
	pgxsci('white')
	pgsch(0.3*ch)
        pgslw(lw)
	for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
	    for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
		if f[isize, imag] > -98:
		    pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
			   '%.1f'%(f[isize, imag]))
		    if b[isize, imag] > 3:
			fout.write('%f %f %f\n'%(mag, size, f[isize, imag]))
	fout.close()
	pgslw(2*lw)
	pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
    pgend()


def fit_baseline_correction_zfit(scheme=default_scheme):
    zfit = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction_zfit.fits')
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    a, b, siga, sigb, corab = [zfit[:,:,i] for i in range(zfit.shape[2])]
    mask_ok = (b > 3).astype(N.float)
    aw = mask_ok * 1.0/siga**2
    bw = mask_ok * 1.0/sigb**2
    ap0 = (-0.5, 0.0, 0.0)
    ap, aperr = bootstrap_fmin(simplified_linear_minfunc, ap0,
                               magvals, sizevals, a, aw)
    #ap = fmin_powell(simplified_linear_minfunc, ap0,
    #	              (magvals, sizevals, a, mask_ok))
    #afit = simplified_linear_function(ap, magvals_sel, sizevals_sel)
    print ap
    print aperr
    bp0 = (15.0, 0.0, 0.0)
    bp, bperr = bootstrap_fmin(simplified_linear_minfunc, bp0,
                        magvals, sizevals, b, bw)
    #bp = fmin_powell(simplified_linear_minfunc, bp0,
    #	              (magvals, sizevals, b, mask_ok))
    #bfit = simplified_linear_function(bp, magvals_sel, sizevals_sel)
    print bp
    print bperr
    p = (ap, bp)
    pickle.dump(p,
		file(data_path+'%s_%s/'%scheme+'baseline_correction_fit.pickle', 'w'))
    #C = d + e*x + f*y +  (g + h*x + i*y) * z
    print
    print 'correction = %5.3f + %5.3f * mag + %5.3f * size + '%tuple(ap)
    print '(%5.3f + %5.3f * mag + %5.3f * size) * z'%tuple(bp)
    print
    print 'zero correction = (%5.3f + %5.3f * mag + %5.3f * size) * z'%tuple(bp)


def fit_baseline_correction(scheme=default_scheme):
    bc = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction.fits')
    bc_sigma = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction_sigma.fits')
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    zvals = zbins.field('CENTRE')
    mask_ok = (bc > -98).astype(N.float)
    w = mask_ok * 1.0/bc_sigma**2
    # unweighted
    w = w.ravel()
    w = mask_ok * w.compress(w > 0.000001).mean()
    #p0 = (-0.5, -0.1, 0.1, 85.0, 3.7, -0.1)
    #p0 = (4.0, 0.2, 0.16, 17000.0, -0.2, -0.6, -0.002, -0.006)
    # k = a + b*z + (c + d*z)*x + (e + f*z)*y + (g + h*z)*y**2
    #p0 = (9.0, 0.0, 0.5, -0.0, 0.6, 5.0, -0.1, 0.0)
    # k = a**(-x-b) + c*exp(-(y-d)**2/e**2)
    #p0 = (1.5, -20.0, 0.5, 3.0, 2.0)
    #k = a0 / (1.0 + N.exp((x0 - b0)/c0)) * exp(-(y-d0)**2/e0**2)
    #p0 = (2.0, 0.0, -16.0, -50.0, -1.0, 0.0, 0.5, 50.0, 3.0, 0.0)
    #p0 = (1.6, -6.0, -17.0, -48.04, -1.15, 0.48, 1.7, 39.0, -0.2, 63.0)
    #p0 = (9.0, -37.9, -14.1, -34.0, -1.35, -5.8, 2.1, 33.7, 0.4, 41.6)
    #p0 = (1.0, 1.0, -0.3, 2.0, 2.0)
    #pb, pberr = pickle.load(file(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.pickle'))
    #a, b, c, d, e, f, g, h, i = pb
    #p0 = (21.0, d, e, 2.0)
    nx = 3
    ny = 3
    nz = 3
    n = (nx, ny, nz)
    p0 = N.array([0.0]*4**3)
    p111 =  [  3.450e-01, 7.901e-01, -1.147e-02, -8.152e-03,
               7.850e-04,-3.405e-02,  7.094e-04,  3.394e-03]
    #p0[:8] = p111
    p222 = [9.39263257e-01,  6.90985164e+00, -3.46143214e+01,  5.44391854e-02,
            3.58124371e-01, -1.15336256e+00, -1.33863367e-02,  5.54497224e-02,
            2.70694291e-02,  3.45291135e-02, -2.00184139e-01,  8.79786798e-01,
            -1.35846536e-03, -1.17594488e-02,  4.36795026e-02,  5.84882172e-04,
            -2.35966297e-03, -1.59468740e-03, -1.45900418e-03,  6.66513593e-03,
            -3.53799258e-02,  8.51293770e-05,  6.28972439e-04, -1.65027554e-03,
            -2.54165643e-05,  1.00024846e-04, -3.33031329e-04]
    p0[:27] = p222
    #p, perr = bootstrap3_fmin(bc_minfunc, p0,
    #                          magvals, sizevals, zvals, bc, w)
    make_perms(n)
    p = fmin_powell(bc_minfunc, p0, (magvals, sizevals, zvals, bc, n, w),
                    ftol = 0.01)
    #p = p0
    perr = [0*i for i in p]
    pickle.dump((n, p, perr),
      file(data_path+'%s_%s/'%scheme+'baseline_correction_fullfit.pickle', 'w'))
    #C = d + e*x + f*y +  (g + h*x + i*y) * z
    print
    #print ('correction = %5.3f + %5.3f * mag + %5.3f * size + '+
    #       '(%5.3f + %5.3f * mag + %5.3f * size) * z')%tuple(p)
    print


def fit_baseline_correction_spline(scheme=default_scheme):
    bc = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction.fits')
    bc_sigma = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction_sigma.fits')
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    zvals = zbins.field('CENTRE')
    mask_ok = (bc > -98).astype(N.float)
    weights = mask_ok * 1.0/bc_sigma**2
    # unweighted
    weights = weights.ravel()
    weights = mask_ok * weights.compress(weights > 0.000001).mean()
    tcklist = []
    xz = []
    yz = []
    bz = []
    wz = []
    zz = []
    nx = len(magvals)
    ny = len(sizevals)
    xlist = N.arange(nx)
    ylist = N.arange(ny)
    xlist = xlist.reshape((len(xlist), 1))
    for i, z in enumerate(zvals):
        gc.collect()
        print
        print z
        x = []
        y = []
        b = []
        w = []
        for xi in range(nx):
            for yi in range(ny):
                if weights[i,xi,yi] > 0.000001:
                    x.append(magvals[xi])
                    y.append(sizevals[yi])
                    b.append(bc[i,xi,yi])
                    w.append(weights[i,xi,yi])
                else:
                    dx = xlist - xi
                    dy = ylist - yi
                    d = N.sqrt(dx**2 + dy**2)
                    N.putmask(d, bc[i] < -98, nx*ny)
                    di = N.argmin(d)
                    xdi = di / len(sizevals)
                    ydi = di - xdi*len(sizevals)
                    dmin = d[xdi, ydi]
                    if dmin <=  100:
                        x.append(magvals[xi])
                        y.append(sizevals[yi])
                        b.append(bc[i,xdi,ydi])
                        w.append(weights[i,xdi,ydi]/dmin)
        #xz.extend(x)
        #yz.extend(y)
        #bz.extend(b)
        #wz.extend(w)
        #zz.extend([z]*len(w))
        if len(x) < 10: continue
        tx = N.array([-23.875, -23.875, -23.875, -23.875,
                      -20.5, 
                      -17.125, -17.125, -17.125, -17.125])
        ty = N.array([0.25,   0.25,   0.25,   0.25,
                      6.0,
                      11.75,  11.75,  11.75,  11.75])

        tck, fp, ier, msg = bisplrep(x, y, b, w,
                                     xb = magvals[0], xe = magvals[-1],
                                     yb = sizevals[0], ye = sizevals[-1],
                                     kx = 3, ky = 3,
                                     task = -1,
                                     #s = len(x),
                                     nxest = len(x)/3,
                                     nyest = len(y)/3,
                                     #tx = None, ty = None,
                                     tx = tx, ty = ty,
                                     full_output = 1, quiet = 0)
        print tck
        print fp
        print ier
        print msg
        tcklist.append(tck)
#     print 'ALL'
#     out = splprep([xz, yz, bz], wz, k = 3,
#                   task = 0, s = 100, full_output = 1)
#     tcku, fp, ier, msg = out
#     tck, u = tcku
#     print u
#     print tck
#     print fp
#     print ier
#     print msg
    pickle.dump(tcklist,
        file(data_path+'%s_%s/'%scheme+'baseline_correction_fit_spline.pickle', 'w'))
    #pickle.dump((tck, u),
    #    file(data_path+'%s_%s/'%scheme+'baseline_correction_fit_spline_all.pickle', 'w'))
    print



def bootstrap_fmin(f, p0, x, y, z, mask, nboot=50):
    plist = N.zeros((nboot, len(p0)), N.float)
    ndata = len(mask.ravel().nonzero()[0])
    for i in range(nboot):
        bootmask = N.zeros(mask.shape)
        while bootmask.sum() < ndata:
            rx = int(N.random.uniform(0, x.shape))
            ry = int(N.random.uniform(0, y.shape))
            if mask[rx,ry] > 0.000001:
                bootmask[rx,ry] += 1
        bootmask *= mask
        plist[i] = fmin_powell(f, p0, (x, y, z, bootmask))
    p = fmin_powell(f, p0, (x, y, z, mask))
    pmed = N.median(plist, axis=0)
    perr = N.array([(scoreatpercentile(plist[:,k], 84.0) -
            scoreatpercentile(plist[:,k], 16.0))/2.0 for k in range(len(p))])
    return p, perr

def bootstrap3_fmin(f, p0, x, y, z, k, mask, nboot=25):
    plist = N.zeros((nboot, len(p0)), N.float)
    ndata = len(mask.ravel().nonzero()[0])
    for i in range(nboot):
        bootmask = N.zeros(mask.shape)
        while bootmask.sum() < ndata:
            rx = int(N.random.uniform(0, x.shape))
            ry = int(N.random.uniform(0, y.shape))
            rz = int(N.random.uniform(0, z.shape))
            if mask[rz,rx,ry] > 0.000001:
                bootmask[rz,rx,ry] += 1
        bootmask *= mask
        plist[i] = fmin_powell(f, p0, (x, y, z, k, bootmask))
    p = fmin_powell(f, p0, (x, y, z, k, mask))
    pmed = N.median(plist, axis=0)
    perr = N.array([(scoreatpercentile(plist[:,k], 84.0) -
            scoreatpercentile(plist[:,k], 16.0))/2.0 for k in range(len(p))])
    return p, perr

def simplified_linear_function(p, x, y):
    a, b, c = p
    z = N.zeros((len(x), len(y)), N.float)
    for i in range(len(x)):
	z[i] = a + b*x[i] + c*y
    return z

def simplified_linear_minfunc(p, x, y, z, w=None):
    f = simplified_linear_function(p, x, y)
    r2 = (z - f)**2
    if w is None:
	s = r2.sum()
    else:
        df = (w > 0).ravel().astype(N.int).sum() - len(p)
	s = (w * r2).sum() / df
    return s

def bc_function_grid3(p, x, y, z):
    a, b, c, d, e, f, g, h = p
    d = 0.0
    k = N.zeros((len(z), len(x), len(y)), N.float)
    for i in range(len(x)):
        for j in range(len(y)):
            k[:,i,j] = a + b*z + (c + d*z)*x[i] + (e + f*z)*y[j] + (g + h*z)*y[j]**2
    return k

def bc_function_grid4(p, x, y, z):
    a, b, c, d, e = p
    k = N.zeros((len(z), len(x), len(y)), N.float)
    for i in range(len(x)):
        for j in range(len(y)):
            k[:,i,j] = a**(x[i]-b) * c*N.exp(-(y[j]-d)**2/e**2)
    return k


def bc_function_grid5(p, x, y, z):
    a1, a2, b1, b2, c1, c2, d1, d2, e1, e2 = p
    a0 = a1 + a2*z
    b0 = b1 + b2*z
    c0 = c1 + c2*z
    d0 = d1 + d2*z
    e0 = e1 + e2*z
    k = N.zeros((len(z), len(x), len(y)), N.float)
    for i in range(len(x)):
        for j in range(len(y)):
            k[:,i,j] = (a0 / (1.0 + N.exp((x[i] - b0)/c0)) * 
                        N.exp(-(y[j]-d0)**2/e0**2))
    return k

def bc_function_grid6(p, x, y, z):
    a, b, c, d, e = p
    dmod = N.array([cosmology.dmod_flat(iz) for iz in z])
    ang_scale = N.array([cosmology.ang_scale_flat(iz) for iz in z])
    mlim = 17.77 - dmod
    a0 = a
    b0 = b + mlim
    c0 = c
    d0 = d * ang_scale
    e0 = e
    k = N.zeros((len(z), len(x), len(y)), N.float)
    for i in range(len(x)):
        for j in range(len(y)):
            k[:,i,j] = (a0 / (1.0 + N.exp((x[i] - b0)/c0)) * 
                        N.exp(-(y[j]-d0)**2/e0**2))
    return k

def bc_function_grid7(p, x, y, z, pb):
    a, b, c, d, e, f, g, h, i = pb
    q, r, s, t = p
    dmod = N.array([cosmology.dmod_flat(iz) for iz in z])
    ang_scale = N.array([cosmology.ang_scale_flat(iz) for iz in z])
    mlim = 17.77 - dmod
    k = N.zeros((len(z), len(x), len(y)), N.float)
    for ii in range(len(x)):
        for j in range(len(y)):
            x0 = b**-(a + h*y[j]**i) + c# + mlim + q
            x1 = r + s*(x0-c)
            k[:,ii,j] = t / (1.0 + N.exp((x0 - x[ii])/x1))
    return k

def bc_function_grid8(p, x, y, z, n):
    # superposition of orthogonal polynomials on intervals
    # -24 < x < -17
    # 0 < y < 12
    # 0 < z < 0.25
    nx, ny, nz = n
    f = chebyt
    normx = (x + 20.5)/3.5
    normy = (y - 6.0)/6.0
    normz = (z - 0.125)/0.125
    kx = N.zeros(len(x), N.float)
    for i in range(nx):
        kx += p[i] * f(i).weight_func(normx)*f(i)(normx)
    ky = N.zeros(len(y), N.float)
    for i in range(ny):
        ky += p[nx+i] * f(i).weight_func(normy)*f(i)(normy)
    kz = N.zeros(len(z), N.float)
    for i in range(nz):
        kz += p[nx+ny+i] * f(i).weight_func(normz)*f(i)(normz)
    kz = kz.reshape((len(z),1,1))
    kx = kx.reshape((1,len(x),1))
    ky = ky.reshape((1,1,len(y)))
    k = kx*ky*kz
    return k

def bc_function_grid(p, x, y, z, n):
    k = N.zeros((len(z), len(x), len(y)), N.float)
    for zii, zi in enumerate(z):
        for xii, xi in enumerate(x):
            k[zii, xii] = poly(xi, y, zi, p)
    return k


def bc_function2(p, x, y, z):
    a, b, c, d, e, f = p
    k = a + d*z + (b + e*z)*x + (c + f*z)*y
    return k

def bc_function3(p, x, y, z):
    a, b, c, d, e, f, g, h = p
    k = a + b*z + (c + d*z)*x + (e + f*z)*y + (g + h*z)*y**2
    return k

def bc_function4(p, x, y, z):
    a, b, c, d, e = p
    k = a**(x-b) * c*exp(-(y-d)**2/e**2)
    return k

def bc_function5(p, x, y, z):
    a1, a2, b1, b2, c1, c2, d1, d2, e1, e2 = p
    a0 = a1 + a2*z
    b0 = b1 + b2*z
    c0 = c1 + c2*z
    d0 = d1 + d2*z
    e0 = e1 + e2*z
    k = a0 / (1.0 + N.exp((x - b0)/c0)) * exp(-(y-d0)**2/e0**2)
    return k

def bc_function6(p, x, y, z, pb):
    a, b, c, d, e, f, g, h, i = pb
    q, r, s, t = p
    dmod = N.array([cosmology.dmod_flat(iz) for iz in z])
    ang_scale = N.array([cosmology.ang_scale_flat(iz) for iz in z])
    mlim = 17.77 - dmod
    x0 = b**-(a + h*y**i) + c + mlim + q
    x1 = r + s*(x0-c)
    k = t / (1.0 + N.exp((x0 - x)/x1))
    return k

def bc_function7(p, x, y, z, n):
    # superposition of orthogonal polynomials on intervals
    # -24 < x < -17
    # 0 < y < 12
    # 0 < z < 0.25
    nx, ny, nz = n
    f = chebyt
    normx = (x + 20.5)/3.5
    normy = (y - 6.0)/6.0
    normz = (z - 0.125)/0.125
    kx = N.zeros(len(x), N.float)
    for i in range(nx):
        kx += p[i] * f(i).weight_func(normx)*f(i)(normx)
    ky = N.zeros(len(y), N.float)
    for i in range(ny):
        ky += p[nx+i] * f(i).weight_func(normy)*f(i)(normy)
    kz = N.zeros(len(z), N.float)
    for i in range(nz):
        kz += p[nx+ny+i] * f(i).weight_func(normz)*f(i)(normz)
    k = kx*ky*kz
    return k

def bc_function(p, x, y, z, n):
    k = poly(x, y, z, p)
    return k

def make_perms(orders):
    global perms
    perms = []
    for x in range(orders[0]+1):
        for y in range(orders[1]+1):
            for z in range(orders[2]+1):
                p = (x, y, z)
                if p not in perms:
                    perms.append(p)
    score = []
    ordered_perms = []
    for o in range(max(orders), -1, -1):
        for p in perms:
            if o in p and p not in ordered_perms:
                ordered_perms.append(p)
    print ordered_perms[::-1]

def poly(x, y, z, coeffs):
    psum = 0.0
    for i, p in enumerate(perms):
        psum += coeffs[i] * x**p[0] * y**p[1] * z**p[2]
    return psum


def bc_minfunc2(p, x, y, z, k, w=None):
    f = bc_function_grid(p, x, y, z)
    r2 = (k - f)**2
    if w is None:
	s = r2.sum()
    else:
        df = (w > 0).ravel().astype(N.int).sum() - len(p)
	s = (w * r2).sum() / df
    print p
    print s
    return s

def bc_minfunc3(p, x, y, z, k, pb, w=None):
    f = bc_function_grid(p, x, y, z, pb)
    r2 = (k - f)**2
    if w is None:
	s = r2.sum()
    else:
        df = (w > 0).ravel().astype(N.int).sum() - len(p)
	s = (w * r2).sum() / df
    print p
    print s
    return s

def bc_minfunc(p, x, y, z, k, n, w=None):
    f = bc_function_grid(p, x, y, z, n)
    r2 = (k - f)**2
    if w is None:
	s = r2.sum()
    else:
        df = (w > 0).ravel().astype(N.int).sum() - len(p)
	s = (w * r2).sum() / df
    print p
    print s
    sys.stdout.flush()
    return s



def plot_baseline_correction_zfitfit(nover=1, scheme=default_scheme):
    ap, bp = pickle.load(file(data_path+'%s_%s/'%scheme+'baseline_correction_fit.pickle'))
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    mag_step = (magbins[0].field('MAX') - magbins[0].field('MIN'))
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = N.logical_and(data > -98, data < 98)
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    zfit = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction_zfit.fits')
    #zfit = zfit[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    zfitb = zfit[:,:,1]
    mask_zfit = zfitb > 3
    mask_zfit = -N.transpose(mask_zfit).astype(N.float)
    magvals_over = N.arange(magvals[min_magbin] - (1-1.0/nover)*mag_step/2.0,
			    magvals[max_magbin] + mag_step - (1-1.0/nover)*mag_step/2.0,
			    mag_step/nover)
    sizevals_over = N.arange(sizevals[min_sizebin] - (1-1.0/nover)*size_step/2.0,
			     sizevals[max_sizebin] + size_step - (1-1.0/nover)*size_step/2.0,
			     size_step/nover)
    a = simplified_linear_function(ap, magvals_over, sizevals_over)
    a = a.transpose()
    b = simplified_linear_function(bp, magvals_over, sizevals_over)
    b = b.transpose()
    for f, l, bg, fg in ((a, 'int', -2.0, 1.0), (b, 'slope', 0.0, 50.0)):
	lab = 'n\dell\u/n\dsp\u baseline correction zfit %s fit'%l
	pgend()
	pgopen(plots_path+'%s_%s/'%scheme+'baseline_correction_zfitfit_%s.ps/cps'%l)
	pgsetup()
	pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
	setup_colour_table('smooth2', 'ramp', flip_colours=False)
	pgsch(1.2*ch)
        pgenv(magbins[min_magbin].field('MIN'),
              magbins[max_magbin].field('MAX'),
              sizebins[min_sizebin].field('MIN'),
              sizebins[max_sizebin].field('MAX'))
	pglab('M\dr\u', 'R\d50\u (kpc)', lab)
	pgimag_s(f, fg, bg)
	pgxsci('white')
	pgsch(0.3*ch)
        pgslw(lw)
	for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
	    for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
		if f[isize, imag] > -98:
		    pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
			   '%.1f'%(f[isize, imag]))
	pgslw(2*lw)
	pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
        pgxsci('black')
        pgcont_s(mask_zfit.astype(N.float), 1, N.array([-0.5]))
    pgend()


def baseline_correction_zfitfit(z, mag, size, zero=False, scheme=default_scheme):
    ap, bp = pickle.load(file(data_path+'%s_%s/'%scheme+'baseline_correction_fit.pickle'))
    if zero:
	a = 0.0
    else:
	a = simplified_linear_function(ap, mag, size)
    b = simplified_linear_function(bp, mag, size)
    correction = a + b*z
    return correction

def plot_baseline_correction_fit(zbin, zero=False, scheme=default_scheme):
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    zsel = zbins[zbin].field('CENTRE')
    corr = baseline_correction_zfitfit(zsel, magvals, sizevals, zero, scheme)
    corr = corr.transpose()
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = N.logical_and(data > -98, data < 98)
    pgend()
    if zero:
	fname = 'baseline_correction_fit_zero_zbin_%02i.ps/cps'
    else:
	fname = 'baseline_correction_fit_zbin_%02i.ps/cps'
    pgopen(plots_path+'%s_%s/'%scheme+fname%zbin)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    zrange = '%6.4f < z < %6.4f'%(zbins[zbin].field('MIN'),
                                  zbins[zbin].field('MAX'))
    if zero:
	lab = 'n\dell\u/n\dsp\u correction zero fit '+zrange
    else:
	lab = 'n\dell\u/n\dsp\u correction fit '+zrange
    pglab('M\dr\u', 'R\d50\u (kpc)', lab)
    pgimag_s(hires(corr), 3.0, -2.0)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
    for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
	for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
	    if corr[isize, imag] > -98 and corr[isize, imag] < 98:
		pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
		       '%.1f'%(corr[isize, imag]))
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    pgslw(2*lw)
    pgcont_s(mask_ok.astype(N.float), 1, N.array([-0.5]))
    pgxsls('solid') 
    pgxsci('black')
    pgsch(ch)
    pgend()

def plot_baseline_correction_fit_spline(zbin, zero=False, scheme=default_scheme):
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    zsel = zbins[zbin].field('CENTRE')
    tck = pickle.load(file(data_path+'%s_%s/'%scheme+'baseline_correction_fit_spline.pickle'))[zbin]
    corr = bisplev(magvals, sizevals, tck)
    corr = corr.transpose()
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    mask_ok = N.logical_and(data > -98, data < 98)
    bc = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction.fits')
    bc_ok = bc[zbin] > -98
    pgend()
    if zero:
	fname = 'baseline_correction_fit_spline_zero_zbin_%02i.ps/cps'
    else:
	fname = 'baseline_correction_fit_spline_zbin_%02i.ps/cps'
    pgopen(plots_path+'%s_%s/'%scheme+fname%zbin)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    zrange = '%6.4f < z < %6.4f'%(zbins[zbin].field('MIN'),
                                  zbins[zbin].field('MAX'))
    if zero:
	lab = 'n\dell\u/n\dsp\u correction zero fit '+zrange
    else:
	lab = 'n\dell\u/n\dsp\u correction fit '+zrange
    pglab('M\dr\u', 'R\d50\u (kpc)', lab)
    pgimag_s(hires(corr), 3.0, -2.0)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
    for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
	for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
	    if corr[isize, imag] > -98 and corr[isize, imag] < 98:
		pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
		       '%.1f'%(corr[isize, imag]))
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    bc_ok = -N.transpose(bc_ok).astype(N.float)
    pgslw(2*lw)
    #pgcont_s(mask_ok, 1, N.array([-0.5]))
    pgcont_s(bc_ok, 1, N.array([-0.5]))
    pgxsls('solid') 
    pgxsci('black')
    pgsch(ch)
    pgend()


def plot_baseline_correction_fit_zbins(zero=False, scheme=default_scheme):
    zbins = pyfits.getdata(data_path+data_file,
                           'REDSHIFT_%s_BINS'%bintype(scheme[1]))
    for zbin in zbins.field('BIN'):
	plot_baseline_correction_fit(zbin, zero, scheme)

def plot_baseline_correction_fullfit(zbin, zero=False, scheme=default_scheme,
                                     scale=True):
    pf = pyfits.open(data_path+data_file)
    set_sizes(pf, scheme[1])
    zbins = pf['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = pf['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = pf['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    size_step = (sizebins[0].field('MAX') - sizebins[0].field('MIN'))
    magvals = magbins[min_magbin:max_magbin+1].field('CENTRE')
    sizevals = sizebins[min_sizebin:max_sizebin+1].field('CENTRE')
    zsel = zbins[zbin].field('CENTRE')
    bc = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction.fits')
    bc = bc[zbin]
    mask_bc = bc > -98
    np, p, perr = pickle.load(file(data_path+'%s_%s/'%scheme+'baseline_correction_fullfit.pickle'))
#     print ('correction = ' +
#            '(%5.3f +- %5.3f) + '%(p[0], perr[0]) +
#            '(%5.3f +- %5.3f) * mag + '%(p[1], perr[1]) +
#            '(%5.3f +- %5.3f) * size + '%(p[2], perr[2]) +
#            '((%5.3f +- %5.3f) + '%(p[3], perr[3]) +
#            '(%5.3f +- %5.3f) * mag + '%(p[4], perr[4]) +
#            '(%5.3f +- %5.3f) * size) * z'%(p[5], perr[5]))
    #pb, pberr = pickle.load(file(data_path+'%s_%s/'%scheme+'ratio_baseline_fit.pickle'))
    corr = bc_function_grid(p, magvals, sizevals, N.array([zsel]), np)
    corr = corr[0]
    corr = corr.transpose()
    data = pyfits.getdata(data_path+'%s_%s/'%scheme+'ratio_baseline.fits')
    #data = data[min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
    mask_ok = N.logical_and(data > -98, data < 98)
    pgend()
    if zero:
	fname = 'baseline_correction_fit_zero_zbin_%02i.ps/cps'
    else:
	fname = 'baseline_correction_fullfit_zbin_%02i.ps/cps'
    if not scale:
        fname = fname.replace('.ps', '_noscale.ps')
    pgopen(plots_path+'%s_%s/'%scheme+fname%zbin)
    pgsetup()
    pgpap(0, float(max_sizebin-min_sizebin+1)/(max_magbin-min_magbin+1))
    setup_colour_table('smooth2', 'ramp', flip_colours=False)
    pgsch(1.2*ch)
    pgenv(magbins[min_magbin].field('MIN'),
          magbins[max_magbin].field('MAX'),
          sizebins[min_sizebin].field('MIN'),
          sizebins[max_sizebin].field('MAX'))
    zrange = '%6.4f < z < %6.4f'%(zbins[zbin].field('MIN'),
                                  zbins[zbin].field('MAX'))
    if zero:
	lab = 'log\d10\u(n\dell\u/n\dsp\u correction zero fit)'
    else:
	lab = 'log\d10\u(n\dell\u/n\dsp\u correction fit)'
    pglab('M\dr\u', 'R\d50\u (kpc)', '')
    pgimag_s(hires(corr), 3.0, -2.0)
    pgbox('BC', 0, 0, 'BC', 0, 0)
    pgxsci('white')
    pgsch(0.3*ch)
    pgslw(lw)
#     for imag, mag in enumerate(magvals[min_magbin:max_magbin+1]):
# 	for isize, size in enumerate(sizevals[min_sizebin:max_sizebin+1]):
# 	    if corr[isize, imag] > -98 and corr[isize, imag] < 98:
# 		pgptxt(mag, size-0.15*size_step, 0.0, 0.5,
# 		       '%.1f'%(corr[isize, imag]))
    mask_ok = -N.transpose(mask_ok).astype(N.float)
    mask_bc = -N.transpose(mask_bc).astype(N.float)
    pgslw(2*lw)
    pgcont_s(mask_ok, 1, N.array([-0.5]))
    pgcont_s(mask_bc, 1, N.array([-0.5]))
    pgxsls('solid') 
    pgslw(lw)
    pgsch(1.2*ch)
    pgstbg(1)
    pgptxt(-19.0, sizebins[max_sizebin].field('MAX')-1, 0.0, 0.5, zrange)
    pgstbg(-1)
    pgxsci('black')
    pgsch(ch)
    if scale:
        pgwedg_s(-2.0, 3.0, 'RI', lab, 1.0, 4.0)
    pgend()

def plot_baseline_correction_fullfit_zbins(zero=False, scheme=default_scheme):
    zbins = pyfits.getdata(data_path+data_file,
                           'REDSHIFT_%s_BINS'%bintype(scheme[1]))
    for zbin in zbins.field('BIN'):
	plot_baseline_correction_fullfit(zbin, zero, scheme)


def hires(x):
    x_hires = N.zeros([s*oversamp for s in x.shape])
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):    
            x_hires[i*oversamp:(i+1)*oversamp,
                        j*oversamp:(j+1)*oversamp] = x[i, j]
    return x_hires


def do_all(start_at=0, scheme=('WEIGHTS', 'CSUW'), fiddle=None, twiddle=None):
    if (scheme[1] in ['FUKUGITA']):
        print 'min_count: 10'
        min_count = 10
    else:
        min_count = 30
    if start_at < 1:
	plot_ratio_zbins(scheme=scheme, fiddle=fiddle, twiddle=twiddle)
        plot_ratio(3, scheme=scheme, scale=False, fiddle=fiddle, twiddle=twiddle)
        plot_ratio(8, scheme=scheme, scale=False, fiddle=fiddle, twiddle=twiddle)
    if start_at < 2:
	determine_baseline(min_count=min_count, scheme=scheme, fiddle=fiddle, twiddle=twiddle)
	plot_ratio_baseline('data', scheme=scheme)
	determine_ratio_baseline_sigma(scheme=scheme)
	plot_ratio_baseline('sigma', scheme=scheme)
        fit_ratio_baseline(scheme=scheme)
	plot_ratio_baseline('fit', scheme=scheme)
	plot_ratio_baseline('fitres', scheme=scheme)
	plot_ratio_baseline('fitnormres', scheme=scheme)
	plot_ratio_baseline_fit(scheme=scheme)
    if start_at < 3:
	determine_baseline_correction(min_count=min_count, scheme=scheme,
                                      fromfit=False, fill=False, fiddle=fiddle, twiddle=twiddle)
	plot_baseline_correction_zbins(scheme=scheme, min_count=min_count, fiddle=fiddle, twiddle=twiddle)
	fit_baseline_correction_zfunc(scheme=scheme)  # just for plot
    if start_at < 4:
        pass
	#fit_baseline_correction(scheme=scheme, fiddle=fiddle, twiddle=twiddle)  # full fit
	#plot_baseline_correction_fullfit_zbins(zero=False, scheme=scheme, fiddle=fiddle, twiddle=twiddle)
	#plot_baseline_correction_fullfit(3, zero=False, scheme=scheme,
        #                                 scale=False, fiddle=fiddle, twiddle=twiddle)
	#plot_baseline_correction_fullfit(8, zero=False, scheme=scheme,
        #                                 scale=False, fiddle=fiddle, twiddle=twiddle)
	#plot_baseline_correction_fit_zbins(zero=True, scheme=scheme, fiddle=fiddle, twiddle=twiddle)
    if start_at < 5:
        correct_dr6z7(fiddle=fiddle, twiddle=twiddle)
    #if start_at < 6:
    #    limited_selection()


def do_fiddle_twiddle():
    global data_path, plots_path
    data_path = '/Users/spb/Work/projects/GalaxyZoo/dr6z7fiddleup/'
    plots_path = data_path+'plots/'
    do_all_fit(fiddle=+4.0)
    data_path = '/Users/spb/Work/projects/GalaxyZoo/dr6z7fiddledn/'
    plots_path = data_path+'plots/'
    do_all_fit(fiddle=-4.0)
    data_path = '/Users/spb/Work/projects/GalaxyZoo/dr6z7twiddleup/'
    plots_path = data_path+'plots/'
    do_all_fit(fiddle=+4.0)
    data_path = '/Users/spb/Work/projects/GalaxyZoo/dr6z7twiddledn/'
    plots_path = data_path+'plots/'
    do_all_fit(fiddle=-4.0)
    data_path = '/Users/spb/Work/projects/GalaxyZoo/dr6z7defaultfit/'
    plots_path = data_path+'plots/'
    do_all_fit()


def do_all_fit(start_at=0, scheme=default_scheme, fiddle=None, twiddle=None):
    if (scheme[1] in ['FUKUGITA']):
        print 'min_count: 10'
        min_count = 10
    else:
        min_count = 50
    if start_at < 1:
	plot_ratio_zbins(scheme=scheme, fiddle=fiddle, twiddle=twiddle)
    if start_at < 2:
	determine_baseline(min_count=min_count, scheme=scheme, fiddle=fiddle, twiddle=twiddle)
	plot_ratio_baseline('data', scheme=scheme)
	determine_ratio_baseline_sigma(scheme=scheme)
	plot_ratio_baseline('sigma', scheme=scheme)
	fit_ratio_baseline(scheme=scheme)
	plot_ratio_baseline('fit', scheme=scheme)
	plot_ratio_baseline('fitres', scheme=scheme)
	plot_ratio_baseline('fitnormres', scheme=scheme)
	plot_ratio_baseline_fit(scheme=scheme)
    if start_at < 3:
	determine_baseline_correction(scheme=scheme, fromfit=True, fiddle=fiddle, twiddle=twiddle)
	plot_baseline_correction_zbins(scheme=scheme)
	fit_baseline_correction_zfunc(scheme=scheme)	
	plot_baseline_correction_zfit(scheme=scheme)
    if start_at < 4:
	fit_baseline_correction_zfit(scheme=scheme)
	plot_baseline_correction_zfitfit(1, scheme=scheme)
	plot_baseline_correction_fit_zbins(zero=False, scheme=scheme)
	plot_baseline_correction_fit_zbins(zero=True, scheme=scheme)
    if start_at < 5:
        correct_dr6z7(fiddle=fiddle, twiddle=twiddle)


def do_all_all(comparisons=False, quality=False):
    scheme_list = [('WEIGHTS', 'CSUW'),
                   ('COUNTS', 'CSUW'),
                   ('WEIGHTS', 'CSDW'),
                   ('COUNTS', 'CSDW'),                   
                   ('WEIGHTS', 'CSW1'),
                   ('COUNTS', 'CSW1'),
                   ('WEIGHTS', 'CSW2'),
                   ('COUNTS', 'CSW2')]
                   #('COUNTS', 'DW'), ('COUNTS', 'UW'),
                   #('WEIGHTS', 'DW'), ('WEIGHTS', 'UW')]
    if comparisons:
        scheme_list += [('COUNTS', 'MOSES'), ('COUNTS', 'FUKUGITA')]

    if quality:
        quality_list = ['good_seeing', 'bad_seeing',
                        'good_sky', 'bad_sky']
    else:
        quality_list = []

    global plots_path, data_path, data_file

    for scheme in scheme_list:
        gc.collect()
        print scheme
        do_all(scheme=scheme)
    
    for q in quality_list:
        if q is None:
            plots_path = '../plots/'
            data_path = '/Users/spb/Work/projects/GalaxyZoo/'
            data_file = 'dr6z7_zoo_selected_counts.fits'
        else:
            plots_path = '../plots/%s/'%q
            data_path = '/Users/spb/Work/projects/GalaxyZoo/%s/'%q
            data_file = 'dr6z7_zoo_selected_counts_%s.fits'%q

        for scheme in scheme_list[:2]:
            gc.collect()
            print scheme, q
            do_all(scheme=scheme)

def correct_dr6z7(subsample=False, fiddle=None, twiddle=None):
    if subsample:
        pinfile = data_path+'%s/dr6z7_zoo_selected_counts_%s.fits'%(subsample, subsample)
    else:
        pinfile = data_path+'dr6z7_zoo_selected_counts.fits'
    p = pyfits.open(pinfile)
    pcols = p['DATA'].columns
    d = p['DATA'].data
    cols = []
    for c in pcols:
	cols.append(pyfits.Column(name=c.name, format=c.format,
				  array=d.field(c.name)))
    if subsample and ('colour_likelihoods' in subsample):
        scheme_list = [('WEIGHTS', 'BALDRYCOLOUR'),
                       ('WEIGHTS', 'JIMCOLOUR'),
                       ('WEIGHTS', 'JIMMORPH')]
    else:
        scheme_list = [('WEIGHTS', 'CSUW')]
                   #    ('COUNTS', 'CSUW'),
                   #    ('COUNTS', 'CSDW')]
                   #('WEIGHTS', 'CSDW'),
                   #('WEIGHTS', 'CSW1'),
                   #('COUNTS', 'CSW1'),
                   #('WEIGHTS', 'CSW2'),
                   #('COUNTS', 'CSW2')]
                   #('COUNTS', 'DW'), ('COUNTS', 'UW'),
                   #('WEIGHTS', 'DW'), ('WEIGHTS', 'UW')]
    for scheme in scheme_list:
        s = scheme[1]
        set_sizes(p, s)
        if ((min_magbin != 0) or (min_sizebin != 0)):
            print 'WARNING! (min_magbin != 0) or (min_sizebin != 0)'
        d_el = p['%s_EL_%s'%scheme].data
        d_sp = p['%s_SP_%s'%scheme].data
        zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
        if fiddle is not None:
            fiddlez = 10**(0.4 * fiddle * (zbins.field('CENTRE') - 0.1))
            d_el_fiddle = N.array([d_el[i] * N.sqrt(f) for i, f in enumerate(fiddlez)], d_el.dtype)
            d_sp_fiddle = N.array([d_sp[i] / N.sqrt(f) for i, f in enumerate(fiddlez)], d_sp.dtype)
            d_tot = d_el + d_sp
            d_fiddle = d_el_fiddle + d_sp_fiddle
            f_fiddle = d_tot / d_fiddle
            d_el_fiddle *= f_fiddle
            d_sp_fiddle *= f_fiddle
            d_el = d_el_fiddle
            d_sp = d_sp_fiddle
        if twiddle is not None and twiddle != 0:
            twiddlez = twiddle * (zbins.field('CENTRE') - 0.1)
            d_el = do_twiddle(d_el, twiddlez)
            d_sp = do_twiddle(d_sp, twiddlez)
        elz = d_el[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
        spz = d_sp[:, min_magbin:max_magbin+1, min_sizebin:max_sizebin+1]
	rb = elz.astype(N.float)/spz
        bc = pyfits.getdata(data_path+'%s_%s/'%scheme+'baseline_correction.fits')
        p_el = d.field('P_EL_%s'%s)
        p_sp = d.field('P_SP_%s'%s)
        #mag = d.field('MR')
        #size = d.field('R50_KPC')
        #z = d.field('REDSHIFT')
        #p, perr = pickle.load(file(data_path+'%s_%s/'%scheme+'baseline_correction_fullfit.pickle'))
        #corr = bc_function(p, mag, size, z)
        ks = bintype(scheme[1])
        magbin = d.field('MR_%s_BIN'%ks)
        N.putmask(magbin, magbin > max_magbin, max_magbin)
        N.putmask(magbin, magbin < min_magbin, min_magbin)
        sizebin = d.field('R50_KPC_%s_BIN'%ks)
        N.putmask(sizebin, sizebin > max_sizebin, max_sizebin)
        N.putmask(sizebin, sizebin < min_sizebin, min_sizebin)
        zbin = d.field('REDSHIFT_%s_BIN'%ks)
        N.putmask(zbin, zbin > max_zbin, max_zbin)
        N.putmask(zbin, zbin < min_zbin, min_zbin)
        corr = bc[zbin, magbin, sizebin]
        corr = N.maximum(corr, 0.0)
        corr = N.minimum(corr, 5.0)
        corr = 10**corr
        ratio_bin = rb[zbin, magbin, sizebin]
        ratio_bin =  N.minimum(ratio_bin, 100000.0)
        ratio_bin =  N.maximum(ratio_bin, 1/100000.0)
        ratio_meas = p_el/p_sp
        ratio_meas =  N.minimum(ratio_meas, 100000.0)
        ratio_meas =  N.maximum(ratio_meas, 1/100000.0)
        if subsample and ('colour_likelihoods' in subsample):
            print 'Using colour likelihoods correction'
            p_sp_bin_corr = 1.0/(ratio_bin/corr + 1.0)
            p_el_bin_corr = 1.0 - p_sp_bin_corr
            p_sp_bin = 1.0/(ratio_bin + 1.0)
            p_el_bin = 1.0 - p_sp_bin
            delta_p_el = p_el_bin_corr - p_el_bin
            delta_p_sp = p_sp_bin_corr - p_sp_bin
            p_el_corr = N.where(p_el > p_sp,
                                p_el + delta_p_el, p_el)
            p_sp_corr = N.where(p_el > p_sp,
                                p_sp + delta_p_sp, p_sp)
            p_x = 0.0
        else:
            p_x = 1.0 - p_el - p_sp
            ratio_corr = ratio_meas / corr
            p_sp_corr = 1.0/(ratio_corr + 1 + p_x/p_sp)
            p_el_corr = 1.0/(1.0/ratio_corr + 1 + p_x/p_el)
            N.putmask(p_el_corr, p_el < eps, 0.0)
            N.putmask(p_el_corr, p_sp < eps, 1.0-p_x)
            N.putmask(p_sp_corr, p_el < eps, 1.0-p_x)
            N.putmask(p_sp_corr, p_sp < eps, 0.0)
        cols.append(pyfits.Column(name='CORR_%s_%s'%scheme,
                                  format='E', array=corr))
        cols.append(pyfits.Column(name='RATIO_BIN_%s_%s'%scheme,
                                  format='E', array=ratio_bin))
        cols.append(pyfits.Column(name='P_EL_CORR_%s_%s'%scheme,
                                  format='E', array=p_el_corr))
        cols.append(pyfits.Column(name='P_SP_CORR_%s_%s'%scheme,
                                  format='E', array=p_sp_corr))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    outfile = data_path+'dr6z7_zoo_corrected.fits'
    if subsample and ('colour_likelihoods' in subsample):
        outfile = outfile.replace('.fits', '_jim.fits')
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    p.close()


def limited_selection(zlimit=0.085):
    infile = data_path+'dr6z7_zoo_corrected.fits'
    p = pyfits.open(infile)
    dmod = cosmology.dmod_flat(zlimit)
    ang_scale = cosmology.ang_scale_flat(zlimit)
    mlim = 17.77 - dmod
    SBapp = 23.0
    d = p['DATA'].data
    z = d.field('REDSHIFT')
    mag = d.field('MR')
    size = d.field('R50_KPC')
    SB = (mag + dmod + 2.5*N.log10(6.283185*(size/ang_scale)**2))
    select = z < zlimit
    select &= mag < mlim
    select &= size > ang_scale
    select &= SB < SBapp
    p['DATA'].data = p['DATA'].data[select]
    outfile = data_path+'dr6z7_zoo_limited_%03i.fits'%int(zlimit*1000)
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    p.writeto(outfile)
    p.close()

def SB_size_excluded(zlimit=0.085):
    infile = '../dr6z7_zoo_corrected_vmax_spclass_ba.fits'
    p = pyfits.open(infile)
    dmod = cosmology.dmod_flat(zlimit)
    ang_scale = cosmology.ang_scale_flat(zlimit)
    mlim = 17.77 - dmod
    SBapp = 23.0
    d = p['DATA'].data
    z = d.field('REDSHIFT')
    mag = d.field('MR')
    size = d.field('R50_KPC')
    SB = (mag + dmod + 2.5*N.log10(6.283185*(size/ang_scale)**2))
    select = z < zlimit
    select &= mag < mlim
    select &= (size <= ang_scale) | (SB > SBapp)
    p['DATA'].data = p['DATA'].data[select]
    outfile = '../dr6z7_zoo_SB_size_excluded_%03i.fits'%int(zlimit*1000)
    file_exists = os.path.isfile(outfile)
    if file_exists:
       os.remove(outfile)
    p.writeto(outfile)
    p.close()


def select_odd_galaxies(lim, zbin=5, nodd=8, scheme=default_scheme):
    import webbrowser
    f = data_path+data_file
    p = pyfits.open(f)
    d = p['DATA'].data
    zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    zsel = zbins.field('CENTRE')[zbin]
    print 'z: %5.3f - %5.3f'%(zbins.field('MIN')[zbin],
                              zbins.field('MAX')[zbin])
    ks = bintype(scheme[1])
    redshift_bin = d.field('REDSHIFT_%s_BIN'%ks)
    select = redshift_bin == zbin
    d = d[select]
    dmod = cosmology.dmod_flat(zsel)
    ang_scale = cosmology.ang_scale_flat(zsel)
    mlim = 17.77 - dmod
    SBapp = 23.0
    close_to_limit = {}
    close_to_limit['mag'] = ((d.field('MR') > (mlim - 0.5)) &
                          (d.field('MR') < mlim)).nonzero()[0]
    close_to_limit['ang'] = ((d.field('R50_KPC') < (1.1*ang_scale)) &
                          (d.field('R50_KPC') > ang_scale)).nonzero()[0]
    SB = (d.field('MR') + dmod + 
          2.5*N.log10(6.283185*(d.field('R50_KPC')/ang_scale)**2))
    close_to_limit['SB'] = ((SB > (SBapp - 0.5)) & (SB < SBapp)).nonzero()[0]
    ids = d.field('objid')
    print 'Close to %s limit:'%lim
    idlist = []
    newwindow = 1
    while True:
        r = int(N.random.uniform(0, len(close_to_limit[lim])))
        i = close_to_limit[lim][r]
        id = ids[i]
        if id not in idlist:
            idlist.append(id)
            print id
            url = 'http://cas.sdss.org/dr6/en/tools/explore/obj.asp?id='+str(id)
            webbrowser.open(url, newwindow)
            newwindow = 0
        if len(idlist) >= nodd: break


def for_others(other):
    if other == 'Mehri':
        names = ('OBJID', 'RA', 'DEC', 'REDSHIFT',
                 'U', 'G', 'R', 'I', 'Z', 
                 'UERR', 'GERR', 'RERR', 'IERR', 'ZERR',
                 'U_MODEL', 'G_MODEL', 'R_MODEL', 'I_MODEL', 'Z_MODEL', 
                 'UERR_MODEL', 'GERR_MODEL', 'RERR_MODEL',
                 'IERR_MODEL', 'ZERR_MODEL',
                 'MU', 'MG', 'MR', 'MI', 'MZ',
                 'MUERR', 'MGERR', 'MRERR', 'MIERR', 'MZERR',
                 'MU_MODEL', 'MG_MODEL', 'MR_MODEL', 'MI_MODEL', 'MZ_MODEL',
                 'MUERR_MODEL', 'MGERR_MODEL', 'MRERR_MODEL',
                 'MIERR_MODEL', 'MZERR_MODEL',
                 'CUR', 'LOGMSTAR_BALDRY06', 'REDNESS_BALDRY06',
                 'R50_KPC',
                 'C4_CLUSTERID', 'C4_REDSHIFT', 'C4_MASS', 'C4_VELDISP', 
                 'C4_RVIR',  'C4_CLOSEST_ID', 'C4_CLOSEST_DISTANCE',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW')
    elif other == 'Bob':
        names = ('OBJID', 'RA', 'DEC',
                 'REDSHIFT', 'REDSHIFTERR', 'REDSHIFTCONF',
                 'U', 'G', 'R', 'I', 'Z', 
                 'UERR', 'GERR', 'RERR', 'IERR', 'ZERR',
                 'U_MODEL', 'G_MODEL', 'R_MODEL', 'I_MODEL', 'Z_MODEL', 
                 'UERR_MODEL', 'GERR_MODEL', 'RERR_MODEL',
                 'IERR_MODEL', 'ZERR_MODEL',
                 'FLAGS', 'STATUS', 'INSIDEMASK', 'PRIMTARGET', 'CALIBSTATUS_U',
                 'MU', 'MG', 'MR', 'MI', 'MZ',
                 'MUERR', 'MGERR', 'MRERR', 'MIERR', 'MZERR',
                 'MU_MODEL', 'MG_MODEL', 'MR_MODEL', 'MI_MODEL', 'MZ_MODEL',
                 'MUERR_MODEL', 'MGERR_MODEL', 'MRERR_MODEL',
                 'MIERR_MODEL', 'MZERR_MODEL',
                 'CUR', 'LOGMSTAR_BALDRY06', 'REDNESS_BALDRY06',
                 'R50_ARCSEC', 'R50_KPC',
                 'C4_CLUSTERID', 'C4_RA', 'C4_DEC', 'C4_REDSHIFT',
                 'C4_MASS', 'C4_VELDISP',
                 'C4_RVIR', 'C4_CLOSEST_CLUSTERID', 
                 'C4_CLOSEST_REDSHIFT', 'C4_CLOSEST_RA',
                 'C4_CLOSEST_DEC', 'C4_CLOSEST_MASS', 'C4_CLOSEST_VELDISP',
                 'C4_CLOSEST_RVIR', 'C4_CLOSEST_DISTANCE',
                 'C4_NORMCLOSEST_CLUSTERID', 
                 'C4_NORMCLOSEST_REDSHIFT', 'C4_NORMCLOSEST_RA',
                 'C4_NORMCLOSEST_DEC', 'C4_NORMCLOSEST_MASS',
                 'C4_NORMCLOSEST_VELDISP',
                 'C4_NORMCLOSEST_RVIR', 'C4_NORMCLOSEST_DISTANCE',
                 'COMOVING_DISTANCE', 
                 'IVAN_DENSITY', 'IVAN_MIN_DENSITY', 'IVAN_MAX_DENSITY',
                 'VAC_SIGMA_5', 'VAC_SIGMA_10', 
                 'P_EL_CSUW', 'P_SP_CSUW', 'P_DK_CSUW', 'P_MG_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW')
    elif other == 'Ramin':
        names = ('OBJID', 'RA', 'DEC', 'REDSHIFT',
                 'U', 'G', 'R', 'I', 'Z', 
                 'UERR', 'GERR', 'RERR', 'IERR', 'ZERR',
                 'MU', 'MG', 'MR', 'MI', 'MZ',
                 'MUERR', 'MGERR', 'MRERR', 'MIERR', 'MZERR',
                 'CUR', 'LOGMSTAR_BALDRY06', 'REDNESS_BALDRY06',
                 'R50_ARCSEC', 'R50_KPC',
                 'P_EL_CSUW', 'P_SP_CSUW', 'P_DK_CSUW', 'P_MG_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW',
                 'REDSHIFTERR', 'REDSHIFTCONF',
                 'FLAGS', 'STATUS', 'INSIDEMASK', 'PRIMTARGET', 'CALIBSTATUS_U',
                 'CW_SPCLASS', 'ACW_SPCLASS', 'EDGEON_SPCLASS',
                 'BA_ISOPHOTE_R', 'BA_MOMENTS_R')
    elif other == 'S0':
        names = ('OBJID',
                 'REDSHIFT',
                 'MU', 'MG', 'MR', 'MI', 'MZ',
                 'CUR', 'LOGMSTAR_BALDRY06', 'REDNESS_BALDRY06', 'R50_KPC', 
                 'C4_NORMCLOSEST_MASS', 'C4_NORMCLOSEST_VELDISP',
                 'C4_NORMCLOSEST_RVIR', 'C4_NORMCLOSEST_DISTANCE',
                 'IVAN_DENSITY', 'IVAN_MIN_DENSITY', 'IVAN_MAX_DENSITY',
                 'P_EL_CSUW', 'P_SP_CSUW', 'P_DK_CSUW', 'P_MG_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW')
    elif other == 'Kris':
        names = ('OBJID', 'RA', 'DEC', 'REDSHIFT',
                 'P_EL_CSUW', 'P_SP_CSUW', 'P_DK_CSUW', 'P_MG_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW')
    elif other == 'Moein':
        names = ('OBJID', 'RA', 'DEC', 'REDSHIFT',
                 'U', 'G', 'R', 'I', 'Z',
                 'UERR', 'GERR', 'RERR', 'IERR', 'ZERR',
                 'U_MODEL', 'G_MODEL', 'R_MODEL', 'I_MODEL', 'Z_MODEL',
                 'UERR_MODEL', 'GERR_MODEL', 'RERR_MODEL',
                 'IERR_MODEL', 'ZERR_MODEL',
                 'MU', 'MG', 'MR', 'MI', 'MZ',
                 'MUERR', 'MGERR', 'MRERR', 'MIERR', 'MZERR',
                 'MU_MODEL', 'MG_MODEL', 'MR_MODEL', 'MI_MODEL', 'MZ_MODEL',
                 'MUERR_MODEL', 'MGERR_MODEL', 'MRERR_MODEL',
                 'MIERR_MODEL', 'MZERR_MODEL',
                 'CUR', 'LOGMSTAR_BALDRY06', 'REDNESS_BALDRY06',
                 'R50_ARCSEC', 'R50_KPC', 'R90_ARCSEC',
                 'P_EL_CSUW', 'P_SP_CSUW', 'P_DK_CSUW', 'P_MG_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW',
                 'CW_SPCLASS', 'ACW_SPCLASS', 'EDGEON_SPCLASS',
                 'IVAN_DENSITY', 'IVAN_MIN_DENSITY', 'IVAN_MAX_DENSITY',
                 'BA_ISOPHOTE_R', 'BA_MOMENTS_R')
    elif other == 'Karen':
        names = ('OBJID', 'RA', 'DEC', 'REDSHIFT',
                 'U', 'G', 'R', 'I', 'Z',
                 'UERR', 'GERR', 'RERR', 'IERR', 'ZERR',
                 'MU', 'MG', 'MR', 'MI', 'MZ',
                 'MUERR', 'MGERR', 'MRERR', 'MIERR', 'MZERR',
                 'MU_MODEL', 'MG_MODEL', 'MR_MODEL', 'MI_MODEL', 'MZ_MODEL',
                 'MUERR_MODEL', 'MGERR_MODEL', 'MRERR_MODEL',
                 'MIERR_MODEL', 'MZERR_MODEL',
                 'LOGMSTAR_BALDRY06', 'REDNESS_BALDRY06',
                 'R50_ARCSEC', 'R50_KPC', 'R90_ARCSEC', 'R90_KPC',
                 'P_EL_CSUW', 'P_SP_CSUW', 'P_DK_CSUW', 'P_MG_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW',
                 'CW_SPCLASS', 'ACW_SPCLASS', 'EDGEON_SPCLASS',
                 'IVAN_DENSITY', 'IVAN_MIN_DENSITY', 'IVAN_MAX_DENSITY',
                 'C4_NORMCLOSEST_MASS', 'C4_NORMCLOSEST_VELDISP',
                 'C4_NORMCLOSEST_RVIR', 'C4_NORMCLOSEST_DISTANCE',
                 'INDR5PLATE',
                 'BA_ISOPHOTE_G', 'BA_ISOPHOTE_R', 'BA_ISOPHOTE_I',
                 'BA_MOMENTS_G', 'BA_MOMENTS_R', 'BA_MOMENTS_I',
                 'ISOA_G', 'ISOA_R', 'ISOA_I',
                 'ISOB_G', 'ISOB_R', 'ISOB_I')
    elif other == 'GANDALF_redsp':
        names = ('OBJID',)
    elif other == 'GANDALF_bluesp':
        names = ('OBJID',)
    elif other == 'GANDALF_redsp_full':
        names = None
    elif other == 'GANDALF_bluesp_full':
        names = None
    elif other == 'dr_table2':
        names = ('OBJID','NVOTE',
                 'P_EL_UW', 'P_CW_UW', 'P_ACW_UW', 'P_EDGEON_UW',
                 'P_DK_UW', 'P_MG_UW', 'P_SP_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW',
                 'P_DK_CORR_WEIGHTS_CSUW', 'P_MG_CORR_WEIGHTS_CSUW')
    else:
        print 'Who for?'
        return
    if other in ('S0', 'Moein'):
        pinfile = data_path+'dr6z7_zoo_limited_085_spclass_ba.fits'
    elif other.startswith('dr_'):
        pinfile = data_path+'dr6z7_zoo_corrected_vmax_spclass.fits'
    else:
        pinfile = data_path+'dr6z7_zoo_corrected_vmax_spclass_ba.fits'
    p = pyfits.open(pinfile)
    pcols = p['DATA'].columns
    d = p['DATA'].data
    if other == 'Mehri':
        select = d.field('C4_CLUSTERID') > 0
        d = d[select]
    elif other.startswith('GANDALF_redsp'):
        select0, select1, select2, select3, select4, select5, select6 = redsp_select(d)
        samplenum = N.zeros(len(d), N.int)
        samplenum[select1.nonzero()[0]] = 1
        samplenum[select2.nonzero()[0]] = 2
        samplenum[select3.nonzero()[0]] = 3
        samplenum[select4.nonzero()[0]] = 4
        samplenum[select5.nonzero()[0]] = 5
        samplenum[select6.nonzero()[0]] = 6
        selectall = select1 | select2 | select3 | select4 | select5 | select6
        d = d[selectall]
        samplenum = samplenum[selectall]
        print 'Total:', len(d)
    elif other.startswith('GANDALF_bluesp'):
        select0, select1, select2, select3, select4, select5, select6 = redsp_select(d, blue=True)
        samplenum = N.zeros(len(d), N.int)
        samplenum[select1.nonzero()[0]] = 1
        samplenum[select2.nonzero()[0]] = 2
        samplenum[select3.nonzero()[0]] = 3
        samplenum[select4.nonzero()[0]] = 4
        samplenum[select5.nonzero()[0]] = 5
        samplenum[select6.nonzero()[0]] = 6
        selectall = select1 | select2 | select3 | select4 | select5 | select6
        d = d[selectall]
        samplenum = samplenum[selectall]
        print 'Total:', len(d)
    sp = d.field('P_SP_CORR_COUNTS_CSUW') > 0.8
    el = d.field('P_EL_CORR_COUNTS_CSUW') > 0.8
    uc = N.logical_not(sp | el)
    if names is None:
        names = d.names
    cols = []
    for n in names:
        for c in pcols:
            if c.name == n:
                cols.append(pyfits.Column(name=c.name, format=c.format,
                                          array=d.field(c.name)))
    if not other.startswith('GANDALF_'):
        cols.append(pyfits.Column(name='SPIRAL', format='I', array=sp))
        cols.append(pyfits.Column(name='ELLIPTICAL', format='I', array=el))
        cols.append(pyfits.Column(name='UNCERTAIN', format='I', array=uc))
    else:
        cols.append(pyfits.Column(name='SAMPLENUM', format='I', array=samplenum))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    outfile = data_path+'dr6z7_zoo_for_%s.fits'%other
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    p.close()
    if other == 'GANDALF_redsp_full':
        t1 = '../dr6z7_zoo_for_GANDALF_redsp_full.fits'
        t2 = '../dr6z7_zoo_for_GANDALF_redsp_spec.fits'
        tout = '../dr6z7_zoo_for_GANDALF_redsp_full_spec.fits'
        file_exists = os.path.isfile(t2)
        if file_exists:
            t = [pyfits.open(i) for i in (t1, t2)]
            combine_fits_tables(t, tout)
    if other == 'GANDALF_bluesp_full':
        t1 = '../dr6z7_zoo_for_GANDALF_bluesp_full.fits'
        t2 = '../dr6z7_zoo_for_GANDALF_bluesp_spec.fits'
        tout = '../dr6z7_zoo_for_GANDALF_bluesp_full_spec.fits'
        file_exists = os.path.isfile(t2)
        if file_exists:
            t = [pyfits.open(i) for i in (t1, t2)]
            combine_fits_tables(t, tout)

def dr_table2(alt=False):
    if not alt:
        names = ('OBJID', 'RA', 'DEC', 'NVOTE',
                 'P_EL_UW', 'P_CW_UW', 'P_ACW_UW', 'P_EDGEON_UW',
                 'P_DK_UW', 'P_MG_UW', 'P_SP_CSUW',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW')
    else:
        names = ('OBJID',
                 'P_EL_CORR_WEIGHTS_CSUW', 'P_SP_CORR_WEIGHTS_CSUW')
    pinfile = data_path+'dr6z7_zoo_corrected_vmax_spclass.fits'
    p = pyfits.open(pinfile)
    pcols = p['DATA'].columns
    d = p['DATA'].data
    sp = d.field('P_SP_CORR_COUNTS_CSUW') > 0.8
    el = d.field('P_EL_CORR_COUNTS_CSUW') > 0.8
    uc = N.logical_not(sp | el)
    cols = []
    for n in names:
        for c in pcols:
            if c.name == n:
                if c.name == 'NVOTE':
                    format = 'I'
                else:
                    format = c.format
                cols.append(pyfits.Column(name=n, format=format,
                                          array=d.field(n)))
    cols.append(pyfits.Column(name='SPIRAL', format='I', array=sp))
    cols.append(pyfits.Column(name='ELLIPTICAL', format='I', array=el))
    cols.append(pyfits.Column(name='UNCERTAIN', format='I', array=uc))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    if not alt:
        outfile = data_path+'dr6z7_zoo_for_dr_table2.fits'
    else:
        outfile = data_path+'dr6z7_zoo_for_dr_table2alt.fits'
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    p.close()

def dr_table3(alt=False):
    z = pyfits.getdata(data_path+'dr6z7_zoo_corrected_vmax_spclass.fits')
    zobjid = z.OBJID
    cols = []
    p = pyfits.open('../data/data080101/uw_stage2.spb.fits')
    pcols = p[1].columns
    d = p[1].data
    objid = d.field('objid').astype(N.int64)
    print 'zobjid is sorted:', N.all(zobjid == N.sort(zobjid))
    print 'objid is sorted:', N.all(objid == N.sort(objid))
    #idx = zobjid.searchsorted(objid)
    #ok = (idx < len(zobjid))
    #zobjid[idx[ok]] = objid[ok]
    idx = N.searchsorted(zobjid, objid)
    # remove indices for unmatched ids
    # first mark those at end
    unmatched = idx >= len(zobjid)
    idx[unmatched] = -999
    # then mark the rest
    unmatched = zobjid[idx] != objid
    idx[unmatched] = -999
    # create array of indices with matches
    matched = N.where(idx >= 0)[0]
    objidmatched = objid[matched]
    unmatched = N.where(idx < 0)[0]
    objidunmatched = objid[unmatched]
    select = unmatched
    if not alt:
        d = d[select]
    colnames = {'objid':  'OBJID',        
                'nobs':   'NVOTE',        
                'histo1': 'P_EL_UW',      
                'histo2': 'P_CW_UW',      
                'histo3': 'P_ACW_UW',     
                'histo4': 'P_EDGEON_UW',  
                'histo5': 'P_DK_UW',      
                'histo6': 'P_MG_UW' }
    for c in pcols:
        if c.name in colnames.keys():
            if c.name == 'NVOTE':
                format = 'I'
            else:
                format = c.format
            cols.append(pyfits.Column(name=colnames[c.name], format=format,
                                      array=d.field(c.name)))
    p = pyfits.open('../data/data080101/cs_uw_stage2.spb.fits')
    pcols = p[1].columns
    d = p[1].data
    d = d[select]
    colnames = {'histo4': 'P_SP_CSUW'}
    for n, m in colnames.items():
        for c in pcols:
            if c.name == n:
                cols.append(pyfits.Column(name=m, format=c.format,
                                          array=d.field(c.name)))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    if not alt:
        outfile = data_path+'dr6z7_zoo_for_dr_table3.fits'
    else:
        outfile = data_path+'dr6z7_zoo_for_dr_table3alt.fits'
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    p.close()


def redsp_select(d, blue=False):
        MR = d.field('MR')
        z = d.field('REDSHIFT')
        red = d.field('REDNESS_BALDRY06') > 0
        sp = d.field('P_SP_CORR_WEIGHTS_CSUW') > 0.8
        if blue:
            # if blue=True, redsp actually selects blue spirals!!!
            redsp = N.logical_not(red) & sp
        else:
            redsp = red & sp
        select0 = redsp.copy()
        MR0 = MR[select0]
        z0 = z[select0]
        print 'all red spirals:'
        print '  n = %i'%len(MR0)
        print '  Mr: min = %.2f, max = %.2f'%(MR0.min(), MR0.max())
        print '  z: min = %.2f, max = %.2f'%(z0.min(), z0.max())
        select6 = redsp.copy()
        select6 &= z >= 0.02
        select6 &= z < 0.03
        select6 &= MR < -18.0
        select6 &= MR >= -18.75
        MR6 = MR[select6]
        z6 = z[select6]
        print 'sample 1:'
        print '  n = %i'%len(MR6)
        print '  Mr: min = %.2f, max = %.2f'%(MR6.min(), MR6.max())
        print '  z: min = %.2f, max = %.2f'%(z6.min(), z6.max())
        select1 = redsp.copy()
        select1 &= z >= 0.025
        select1 &= z < 0.045
        select1 &= MR < -18.75
        select1 &= MR >= -19.5
        MR1 = MR[select1]
        z1 = z[select1]
        print 'sample 2:'
        print '  n = %i'%len(MR1)
        print '  Mr: min = %.2f, max = %.2f'%(MR1.min(), MR1.max())
        print '  z: min = %.2f, max = %.2f'%(z1.min(), z1.max())
        select2 = redsp.copy()
        select2 &= z >= 0.03
        select2 &= z < 0.06
        select2 &= MR < -19.5
        select2 &= MR >= -20.25
        MR2 = MR[select2]
        z2 = z[select2]
        print 'sample 3:'
        print '  n = %i'%len(MR2)
        print '  Mr: min = %.2f, max = %.2f'%(MR2.min(), MR2.max())
        print '  z: min = %.2f, max = %.2f'%(z2.min(), z2.max())
        select3 = redsp.copy()
        select3 &= z >= 0.035
        select3 &= z < 0.08
        select3 &= MR < -20.25
        select3 &= MR >= -21.0
        MR3 = MR[select3]
        z3 = z[select3]
        print 'sample 4:'
        print '  n = %i'%len(MR3)
        print '  Mr: min = %.2f, max = %.2f'%(MR3.min(), MR3.max())
        print '  z: min = %.2f, max = %.2f'%(z3.min(), z3.max())
        select4 = redsp.copy()
        select4 &= z >= 0.04
        select4 &= z < 0.085
        select4 &= MR < -21.0
        select4 &= MR >= -21.75
        MR4 = MR[select4]
        z4 = z[select4]
        print 'sample 5:'
        print '  n = %i'%len(MR4)
        print '  Mr: min = %.2f, max = %.2f'%(MR4.min(), MR4.max())
        print '  z: min = %.2f, max = %.2f'%(z4.min(), z4.max())
        select5 = redsp.copy()
        select5 &= z >= 0.055
        select5 &= z < 0.085
        select5 &= MR < -21.75
        select5 &= MR >= -24.0
        MR5 = MR[select5]
        z5 = z[select5]
        print 'sample 6:'
        print '  n = %i'%len(MR5)
        print '  Mr: min = %.2f, max = %.2f'%(MR5.min(), MR5.max())
        print '  z: min = %.2f, max = %.2f'%(z5.min(), z5.max())
        print 'total of six redsp samples: n = %i'%(len(MR6)+len(MR1)+len(MR2)+len(MR3)+len(MR4)+len(MR5))
        pgopen('../plots/redsp_selection.ps/cps')
        pgsetup()
        pgenv(0.0, 0.1, -17, -24)
        pglab('z', 'M\dr\u', '')
        pgxsci('black')
        pgxpt(z0, MR0, 'point') 
        pgxsci('grape')
        pgxpt(z5, MR5, 'dot') 
        pgxsci('orange')
        pgxpt(z4, MR4, 'dot') 
        pgxsci('blue')
        pgxpt(z3, MR3, 'dot') 
        pgxsci('green')
        pgxpt(z2, MR2, 'dot') 
        pgxsci('red')
        pgxpt(z1, MR1, 'dot') 
        pgxsci('aqua')
        pgxpt(z6, MR6, 'dot') 
        pgclos()
        return select0, select6, select1, select2, select3, select4, select5

def for_lumfn(plotz=False, all=False):
    #pinfile = '../dr6z7_zoo_corrected_vmax_spclass.fits'
    pinfile = '../dr6z7_zoo_corrected_jim.fits'
    p_zoo = pyfits.open(pinfile)
    d_zoo = p_zoo['DATA'].data
    z = d_zoo.field('REDSHIFT')
    # select only objects with z<0.2
    select = z < 0.2
    z = z[select]
    d_zoo = d_zoo[select]
    id_zoo_sortargs = d_zoo.field('OBJID').argsort()
    d_zoo = d_zoo[id_zoo_sortargs]
    id_zoo = d_zoo.field('OBJID')
    print '%i objects in Zoo table'%len(id_zoo)
    # read Will and Jim's DR5 sample DR6 objids
    jinfile = '../../GZ_lumfn/jim_positions_dr6.fits'
    j = pyfits.open(jinfile)
    jcols = j[1].columns
    d_jim = j[1].data
    jfullinfile = '../../GZ_lumfn/FINAL.allz.gal.RAsorted.uncertainties.fits'
    jfull = pyfits.open(jfullinfile)
    jfullcols = jfull[1].columns
    d_jimfull = jfull[1].data
    id_jim = d_jim.field('DR6OBJID')
    z_jim = d_jimfull.field('REDSHIFT')
    r_jim = d_jimfull.field('R')
    print d_jimfull.names
    if plotz:
        pgopen('/aqt')
        pgsetup()
        pgenv(0.0, 0.55, 10, 21)
        pgxpt(z_jim, r_jim, 'point')
        pgclos()
    print "%i objects in Jim's table"%len(id_jim)
    print "%i objects in Jim's full table"%len(d_jimfull)
    print "Jim's tables match:", N.all((d_jimfull.field('REDSHIFT') - d_jim.field('z')) < 0.0001)
    print id_zoo
    print id_jim
    i_zoo = N.searchsorted(id_zoo, id_jim)
    # remove indices for unmatched ids
    # first mark those at end
    unmatched = i_zoo >= len(id_zoo)
    i_zoo[unmatched] = -999
    # then mark the rest
    unmatched = id_zoo[i_zoo] != id_jim
    i_zoo[unmatched] = -999
    # create array of indices with matches
    matched = N.where(i_zoo >= 0)[0]
    print "%i objects in Jim's table match the Zoo table"%len(matched)
    # then cut data arrays down to just matches
    d_jim = d_jim[matched]
    d_jimfull = d_jimfull[matched]
    cgr_jim = d_jimfull.field('absG') - d_jimfull.field('absR')
    mr_jim = d_jimfull.field('absR')
    i_zoo = i_zoo[matched]
    d_zoo = d_zoo[i_zoo]
    id_zoo = d_zoo.field('OBJID')
    p_sp_raw = d_zoo.field('P_SP_CSUW')
    p_el_raw = d_zoo.field('P_EL_CSUW')
    p_sp = d_zoo.field('P_SP_CORR_WEIGHTS_CSUW')
    p_el = d_zoo.field('P_EL_CORR_WEIGHTS_CSUW')
    p_sp_baldrycolour = d_zoo.field('P_SP_BALDRYCOLOUR')
    p_el_baldrycolour = d_zoo.field('P_EL_BALDRYCOLOUR')
    p_sp_jimcolour = d_zoo.field('P_SP_JIMCOLOUR')
    p_el_jimcolour = d_zoo.field('P_EL_JIMCOLOUR')
    p_sp_jimmorph = d_zoo.field('P_SP_JIMMORPH')
    p_el_jimmorph = d_zoo.field('P_EL_JIMMORPH')
    p_sp_jimmorphcorr = d_zoo.field('P_SP_JIMMORPHCORR')
    p_el_jimmorphcorr = d_zoo.field('P_EL_JIMMORPHCORR')
    p_sp_corr_baldrycolour = d_zoo.field('P_SP_CORR_WEIGHTS_BALDRYCOLOUR')
    p_el_corr_baldrycolour = d_zoo.field('P_EL_CORR_WEIGHTS_BALDRYCOLOUR')
    p_sp_corr_jimcolour = d_zoo.field('P_SP_CORR_WEIGHTS_JIMCOLOUR')
    p_el_corr_jimcolour = d_zoo.field('P_EL_CORR_WEIGHTS_JIMCOLOUR')
    p_sp_corr_jimmorph = d_zoo.field('P_SP_CORR_WEIGHTS_JIMMORPH')
    p_el_corr_jimmorph = d_zoo.field('P_EL_CORR_WEIGHTS_JIMMORPH')
    r50 = d_zoo.field('R50_ARCSEC')
    ci = d_zoo.field('R50_ARCSEC')/d_zoo.field('R90_ARCSEC')
    redness = d_zoo.field('REDNESS_BALDRY06')
    mstar = d_zoo.field('LOGMSTAR_BALDRY06')
    gc.collect()
    #sp = p_sp > 0.8
    #el = p_el > 0.8
    #uc = N.logical_not(sp | el)
    cols = []
    if all:
        prefix = 'JIM_'
    else:
        prefix = ''
    for c in jfullcols:
        cols.append(pyfits.Column(name=prefix+c.name.upper(), format=c.format,
                                  array=d_jimfull.field(c.name)))
    if all:
        for c in p_zoo['DATA'].columns:
            cols.append(pyfits.Column(name=c.name.upper(), format=c.format,
                                  array=d_zoo.field(c.name)))
    else:
        cols.append(pyfits.Column(name='DR6OBJID', format='K', array=id_zoo))
        cols.append(pyfits.Column(name='ZOO_REDSHIFT', format='E', array=z))
        cols.append(pyfits.Column(name='P_EL', format='E', array=p_el))
        cols.append(pyfits.Column(name='P_SP', format='E', array=p_sp))
        cols.append(pyfits.Column(name='P_EL_RAW', format='E', array=p_el_raw))
        cols.append(pyfits.Column(name='P_SP_RAW', format='E', array=p_sp_raw))
        cols.append(pyfits.Column(name='P_RED_BALDRYCOLOUR', format='E',
                                  array=p_el_baldrycolour))
        cols.append(pyfits.Column(name='P_BLUE_BALDRYCOLOUR', format='E',
                                  array=p_sp_baldrycolour))
        cols.append(pyfits.Column(name='P_RED_JIMCOLOUR', format='E',
                                  array=p_el_jimcolour))
        cols.append(pyfits.Column(name='P_BLUE_JIMCOLOUR', format='E',
                                  array=p_sp_jimcolour))
        cols.append(pyfits.Column(name='P_EL_JIMMORPH', format='E',
                                  array=p_el_jimmorph))
        cols.append(pyfits.Column(name='P_SP_JIMMORPH', format='E',
                                  array=p_sp_jimmorph))
        cols.append(pyfits.Column(name='P_EL_JIMMORPHCORR', format='E',
                                  array=p_el_jimmorphcorr))
        cols.append(pyfits.Column(name='P_SP_JIMMORPHCORR', format='E',
                                  array=p_sp_jimmorphcorr))
        cols.append(pyfits.Column(name='P_RED_CORR_BALDRYCOLOUR', format='E',
                                  array=p_el_corr_baldrycolour))
        cols.append(pyfits.Column(name='P_BLUE_CORR_BALDRYCOLOUR', format='E',
                                  array=p_sp_corr_baldrycolour))
        cols.append(pyfits.Column(name='P_RED_CORR_JIMCOLOUR', format='E',
                                  array=p_el_corr_jimcolour))
        cols.append(pyfits.Column(name='P_BLUE_CORR_JIMCOLOUR', format='E',
                                  array=p_sp_corr_jimcolour))
        cols.append(pyfits.Column(name='P_EL_CORR_JIMMORPH', format='E',
                                  array=p_el_corr_jimmorph))
        cols.append(pyfits.Column(name='P_SP_CORR_JIMMORPH', format='E',
                                  array=p_sp_corr_jimmorph))
        cols.append(pyfits.Column(name='INV_CONCENTRATION', format='E', array=ci))
        cols.append(pyfits.Column(name='REDNESS', format='E', array=redness))
        cols.append(pyfits.Column(name='MSTAR', format='E', array=mstar))
        cols.append(pyfits.Column(name='CGR', format='E', array=cgr_jim))
        cols.append(pyfits.Column(name='MR', format='E', array=mr_jim))
        cols.append(pyfits.Column(name='R50_ARCSEC', format='E', array=r50))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    if all:
        outfile = '../../GZ_lumfn/zoo_for_lumfn_all.fits'
    else:
        outfile = '../../GZ_lumfn/zoo_for_lumfn.fits'
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    j.close()

def for_ari():
    pinfile = '../dr6z7_zoo_corrected.fits'
    d_zoo = pyfits.getdata(pinfile, 'DATA')
    #select = d_zoo.field('REDSHIFT') > 0.001
    #select &= d_zoo.field('REDSHIFT') < 1.0
    #select &= d_zoo.field('INDR4PLATE').astype(N.bool)
    select = N.ones(len(d_zoo), N.bool)
    z = d_zoo.field('REDSHIFT')[select]
    id_zoo = d_zoo.field('OBJID')[select]
    ra = d_zoo.field('RA')[select].astype(N.float64)
    dec = d_zoo.field('DEC')[select].astype(N.float64)
    p_sp = d_zoo.field('P_SP_CORR_WEIGHTS_CSUW')[select]
    p_el = d_zoo.field('P_EL_CORR_WEIGHTS_CSUW')[select]
    #ci = d_zoo.field('R50_ARCSEC')[select]/d_zoo.field('R90_ARCSEC')[select]
    #redness = d_zoo.field('REDNESS_BALDRY06')[select]
    #mstar = d_zoo.field('LOGMSTAR_BALDRY06')[select]
    del d_zoo
    gc.collect()
    print '%i objects in Zoo table'%len(id_zoo)
    print 'id_zoo is sorted:', N.all(id_zoo == N.sort(id_zoo))
    # read Ari's sample positions
    jinfile = '../../GZ_inclination/maller08.fits'
    j = pyfits.open(jinfile)
    jcols = j[1].columns
    d_ari = j[1].data
    id_ari = d_ari.field('vacid')
    ra_ari = d_ari.field('ra').astype(N.float64)
    dec_ari = d_ari.field('dec').astype(N.float64)
    n_ari = len(d_ari)
    print "%i objects in Ari's table"%n_ari
    # for each of Ari's objects calculate the angular distance
    # to each object in the gz table
    ari_matched = []
    zoo_matched = []
    DPIBY2 = 0.5 * pi
    convRA = pi / 180.0
    convDEC = pi / 180.0
    theta_ari = dec_ari*convDEC + DPIBY2
    sintheta_ari = N.sin(theta_ari)
    costheta_ari = N.cos(theta_ari)
    theta = dec*convDEC + DPIBY2
    sintheta = N.sin(theta)
    costheta = N.cos(theta)
    ra_ari *= convRA
    ra *= convRA
    for i in range(n_ari):
        if i%100 == 0: print i,
        sys.stdout.flush()
        cosgamma = (sintheta_ari[i] * sintheta * N.cos(ra_ari[i]-ra) + 
                    costheta_ari[i] * costheta)
        jmin = cosgamma.argmax()
        mindist = N.arccos(min(cosgamma[jmin], 1.0))
        while mindist >= 2*pi:
            print 'Warning mindist >= 2*pi'
            mindist -= 2*pi
        while mindist < 0.0:
            print 'Warning mindist < 0.0'
            mindist *= -1
        if mindist < 0.0001:
            ari_matched.append(i)
            zoo_matched.append(jmin)
    print "%i objects in Ari's table match the Zoo table"%len(ari_matched)
    # then cut data arrays down to just matches
    d_ari = d_ari[ari_matched]
    id_zoo = id_zoo[zoo_matched]
    z = z[zoo_matched]
    p_sp = p_sp[zoo_matched]
    p_el = p_el[zoo_matched]
    #ci = ci[zoo_matched]
    #redness = redness[zoo_matched]
    #mstar = mstar[zoo_matched]
    gc.collect()
    #sp = p_sp > 0.8
    #el = p_el > 0.8
    #uc = N.logical_not(sp | el)
    cols = []
    for c in jcols:
        cols.append(pyfits.Column(name=c.name.upper(), format=c.format,
                                  array=d_ari.field(c.name)))
    cols.append(pyfits.Column(name='DR6OBJID', format='K', array=id_zoo))
    cols.append(pyfits.Column(name='ZOO_REDSHIFT', format='E', array=z))
    cols.append(pyfits.Column(name='P_EL', format='E', array=p_el))
    cols.append(pyfits.Column(name='P_SP', format='E', array=p_sp))
    #cols.append(pyfits.Column(name='INV_CONCENTRATION', format='E', array=ci))
    #cols.append(pyfits.Column(name='REDNESS', format='E', array=redness))
    #cols.append(pyfits.Column(name='MSTAR', format='E', array=mstar))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    outfile = '../../GZ_inclination/zoo_for_ari.fits'
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    j.close()

def assess_correction(f='b.dat'):
    # for the blue etypes
    pinfile = '../dr6z7_zoo_corrected.fits'
    data_corrected = pyfits.getdata(pinfile)
    b_objid = N.array([int(i) for i in file('../%s'%f)])
    nb_objid = N.array([int(i) for i in file('../%s'%f)]) 
    objids = data_corrected.OBJID
    b_select = []
    for o in b_objid:
        b_select.extend(N.where(objids==o)[0])
    n = len(b_select)
    b_p_el = data_corrected.P_EL_CSUW[b_select]
    b_p_el_c = data_corrected.P_EL_CORR_COUNTS_CSUW[b_select]
    b_ratio = b_p_el_c/b_p_el
    n_lt = {}
    n_lt_c = {}
    pclist = (80, 75, 70, 65, 60)
    for pc in pclist:
        n_lt[pc] = len(N.where(b_p_el < (pc/100.0))[0])/float(n)
        n_lt_c[pc] = len(N.where(b_p_el_c < (pc/100.0))[0])/float(n)
    print '%s'%f
    print 'p_el_c / p_el = %.3f (stdev = %.3f)'%(b_ratio.mean(), b_ratio.std())
    print '%18s  %8s  %8s'%('fraction less than',
                            'p_el', 'p_el_c')
    for pc in pclist:
        print '%18.2f  %8.3f  %8.3f'%(pc/100.0, n_lt[pc], n_lt_c[pc])
    

def gz_conf_for_lumfn():
    scheme_list = [('WEIGHTS', 'BALDRYCOLOUR'),
                   ('WEIGHTS', 'JIMCOLOUR'),
                   ('WEIGHTS', 'JIMMORPH')]
    for scheme in scheme_list:
        gz_confidence_abs_bins(thresh=None,
                               label=scheme[1].lower(),
                               scheme=scheme, jim=True)

def gz_confidence_abs_bins(thresh=None,  # p threshold, or None for greater
                           label='greater',
                           scheme=default_scheme,
                           jim=False):
    f = data_path+data_file
    p = pyfits.open(f)
    zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])].data
    magbins = p['MR_%s_BINS'%bintype(scheme[1])].data
    sizebins = p['R50_KPC_%s_BINS'%bintype(scheme[1])].data
    if jim:
        pinfile = data_path+'../dr6_zoo_corrected_jim.fits'
    else:
        pinfile = data_path+'dr6z7_zoo_corrected.fits'
    data = pyfits.getdata(pinfile)
    f_unclass = N.zeros(len(data), N.float)
    f_misclass = N.zeros(len(data), N.float)
    avcorr = N.zeros(len(data), N.float)
    stdcorr = N.zeros(len(data), N.float)
    f_unclass_grid = N.zeros((len(zbins),len(magbins),len(sizebins)), N.float)
    f_misclass_grid = N.zeros((len(zbins),len(magbins),len(sizebins)), N.float)
    avcorr_grid = N.zeros((len(zbins),len(magbins),len(sizebins)), N.float)
    stdcorr_grid = N.zeros((len(zbins),len(magbins),len(sizebins)), N.float)
    print 'iz =',
    for izi, iz in enumerate(zbins.field('BIN')):
        print iz,
        sys.stdout.flush()
        z_selection = data.field('REDSHIFT_ZOO_BIN') == iz
        if not N.any(z_selection):  continue
        for imagi, imag in enumerate(magbins.field('BIN')):
            mag_selection = z_selection & (data.field('MR_ZOO_BIN') == imag)
            if not N.any(mag_selection):  continue
            #print 'mag =', mag
            for isizei, isize in enumerate(sizebins.field('BIN')):
                #print 'size =', size
                selection = mag_selection & (data.field('R50_KPC_ZOO_BIN') == isize)
                #print iz, imag, isize, len(selection.nonzero()[0])
                if not N.any(selection):  continue
                d = data[selection]
                n = float(len(d))
                p_sp_raw = d.field('P_SP_%s'%scheme[1])
                p_el_raw = d.field('P_EL_%s'%scheme[1])
                if thresh is None:
                    s_raw = p_sp_raw >= p_el_raw
                    e_raw = N.logical_not(s_raw)
                else:
                    s_raw = p_sp_raw >= thresh
                    e_raw = p_el_raw >= thresh
                p_sp_corr = d.field('P_SP_CORR_%s_%s'%scheme)
                p_el_corr = d.field('P_EL_CORR_%s_%s'%scheme)
                if thresh is None:
                    s_corr = p_sp_corr >= p_el_corr
                    e_corr = N.logical_not(s_corr)
                    unclass = N.zeros(n, N.bool)
                else:
                    s_corr = p_sp_corr >= thresh
                    e_corr = p_el_corr >= thresh
                    unclass = N.logical_not(s_corr | e_corr)
                f_unclass_i = len(N.where(unclass)[0])/n
                f_unclass[selection] = f_unclass_i
                misclass = (s_raw & e_corr) | (e_raw & s_corr)
                f_misclass_i = len(N.where(misclass)[0])/n
                f_misclass[selection] = f_misclass_i
                avcorr_i = N.mean(p_sp_corr - p_sp_raw)
                avcorr[selection] = avcorr_i
                stdcorr_i = N.std(p_sp_corr - p_sp_raw)
                stdcorr[selection] = stdcorr_i
                f_unclass_grid[izi, imagi, isizei] = f_unclass_i
                f_misclass_grid[izi, imagi, isizei] = f_misclass_i
                avcorr_grid[izi, imagi, isizei] = avcorr_i
                stdcorr_grid[izi, imagi, isizei] = stdcorr_i
    cols = []
    cols.append(pyfits.Column(name='OBJID',
                              format='K', array=data.field('objid')))
    cols.append(pyfits.Column(name='F_UNCLASS_ABS_BIN',
                              format='E', array=f_unclass))
    cols.append(pyfits.Column(name='F_MISCLASS_ABS_BIN',
                              format='E', array=f_misclass))
    cols.append(pyfits.Column(name='AVCORR_ABS_BIN',
                              format='E', array=avcorr))
    cols.append(pyfits.Column(name='STDCORR_ABS_BIN',
                              format='E', array=stdcorr))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    if jim:
        outfile = data_path+'../dr6_zoo_confidence_abs_bin_%s.fits'%label
    else:
        outfile = data_path+'dr6z7_zoo_confidence_abs_bin_%s.fits'%label
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    hdu_primary = pyfits.PrimaryHDU()
    hdu_zbins = p['REDSHIFT_%s_BINS'%bintype(scheme[1])]
    hdu_zbins.name = 'zbins'
    hdu_magbins = p['MR_%s_BINS'%bintype(scheme[1])]
    hdu_magbins.name = 'magbins'
    hdu_sizebins = p['R50_KPC_%s_BINS'%bintype(scheme[1])]
    hdu_sizebins.name = 'sizebins'
    hdu_f_unclass = pyfits.ImageHDU(f_unclass_grid)
    hdu_f_unclass.name = 'f_unclass'
    hdu_f_misclass = pyfits.ImageHDU(f_misclass_grid)
    hdu_f_misclass.name = 'f_misclass'
    hdu_avcorr = pyfits.ImageHDU(avcorr_grid)
    hdu_avcorr.name = 'avcorr'
    hdu_stdcorr = pyfits.ImageHDU(stdcorr_grid)
    hdu_stdcorr.name = 'stdcorr'
    hdulist = pyfits.HDUList([hdu_primary, hdu_zbins, hdu_magbins,
                              hdu_sizebins, hdu_f_unclass,
                              hdu_f_misclass, hdu_avcorr, hdu_stdcorr])
    if jim:
        outfile = data_path+'../dr6_zoo_conf_abs_bin_%s_grids.fits'%label
    else:
        outfile = data_path+'dr6z7_zoo_conf_abs_bin_%s_grids.fits'%label
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    p.close()


def gz_confidence_app_bins(thresh=None, 
                           label='greater',
                           scheme=default_scheme):
    pinfile = data_path+'dr6z7_zoo_corrected.fits'
    data = pyfits.getdata(pinfile)
    mag = data.field('R')
    mag_low = scoreatpercentile(mag, 0.1)
    mag_high = scoreatpercentile(mag, 99.9)
    mag_delta = 0.125
    #mag_delta = 1.0
    magbins = N.arange(mag_low-3*mag_delta, mag_high+3*mag_delta+eps, mag_delta)
    size = data.field('R50_ARCSEC')
    size_low = scoreatpercentile(size, 0.1)
    size_high = scoreatpercentile(size, 99.9)
    size_delta = 0.25
    #size_delta = 2.0
    sizebins = N.arange(size_low-3*size_delta, size_high+3*size_delta+eps, size_delta)
    f_unclass = N.zeros(len(data), N.float)
    f_misclass = N.zeros(len(data), N.float)
    avcorr = N.zeros(len(data), N.float)
    stdcorr = N.zeros(len(data), N.float)
    print len(magbins), magbins
    print len(sizebins), sizebins
    print len(magbins) * len(sizebins)
    f_unclass_grid = N.zeros((len(magbins),len(sizebins)), N.float)
    f_misclass_grid = N.zeros((len(magbins),len(sizebins)), N.float)
    avcorr_grid = N.zeros((len(magbins),len(sizebins)), N.float)
    stdcorr_grid = N.zeros((len(magbins),len(sizebins)), N.float)
    for imagi, imag in enumerate(magbins):
        mag_selection = (mag >= imag) & (mag < (imag+mag_delta))
        if not N.any(mag_selection):  continue
        print 'imag =', imag
        for isizei, isize in enumerate(sizebins):
            selection = (mag_selection & (size >= isize) & 
                         (size < (isize+size_delta)))
            if not N.any(selection):  continue
            d = data[selection]
            n = float(len(d))
            p_sp_raw = d.field('P_SP_%s'%scheme[1])
            p_el_raw = d.field('P_EL_%s'%scheme[1])
            if thresh is None:
                s_raw = p_sp_raw >= p_el_raw
                e_raw = N.logical_not(s_raw)
            else:
                s_raw = p_sp_raw >= thresh
                e_raw = p_el_raw >= thresh
            p_sp_corr = d.field('P_SP_CORR_%s_%s'%scheme)
            p_el_corr = d.field('P_EL_CORR_%s_%s'%scheme)
            if thresh is None:
                s_corr = p_sp_corr >= p_el_corr
                e_corr = N.logical_not(s_corr)
                unclass = N.zeros(n, N.bool)
            else:
                s_corr = p_sp_corr >= thresh
                e_corr = p_el_corr >= thresh
                unclass = N.logical_not(s_corr | e_corr)
            f_unclass_i = len(N.where(unclass)[0])/n
            f_unclass[selection] = f_unclass_i
            misclass = (s_raw & e_corr) | (e_raw & s_corr)
            f_misclass_i = len(N.where(misclass)[0])/n
            f_misclass[selection] = f_misclass_i
            avcorr_i = N.mean(p_sp_corr - p_sp_raw)
            avcorr[selection] = avcorr_i
            stdcorr_i = N.std(p_sp_corr - p_sp_raw)
            stdcorr[selection] = stdcorr_i
            f_unclass_grid[imagi, isizei] = f_unclass_i
            f_misclass_grid[imagi, isizei] = f_misclass_i
            avcorr_grid[imagi, isizei] = avcorr_i
            stdcorr_grid[imagi, isizei] = stdcorr_i
    cols = []
    cols.append(pyfits.Column(name='OBJID',
                              format='K', array=data.field('objid')))
    cols.append(pyfits.Column(name='F_UNCLASS_APP_BIN',
                              format='E', array=f_unclass))
    cols.append(pyfits.Column(name='F_MISCLASS_APP_BIN',
                              format='E', array=f_misclass))
    cols.append(pyfits.Column(name='AVCORR_APP_BIN',
                              format='E', array=avcorr))
    cols.append(pyfits.Column(name='STDCORR_APP_BIN',
                              format='E', array=stdcorr))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    outfile = data_path+'dr6z7_zoo_confidence_app_bin_%s.fits'%label
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)
    hdu_primary = pyfits.PrimaryHDU()
    hdu_magbins = pyfits.ImageHDU(magbins)
    hdu_magbins.name = 'magbins'
    hdu_sizebins = pyfits.ImageHDU(sizebins)
    hdu_sizebins.name = 'sizebins'
    hdu_f_unclass = pyfits.ImageHDU(f_unclass_grid)
    hdu_f_unclass.name = 'f_unclass'
    hdu_f_misclass = pyfits.ImageHDU(f_misclass_grid)
    hdu_f_misclass.name = 'f_misclass'
    hdu_avcorr = pyfits.ImageHDU(avcorr_grid)
    hdu_avcorr.name = 'avcorr'
    hdu_stdcorr = pyfits.ImageHDU(stdcorr_grid)
    hdu_stdcorr.name = 'stdcorr'
    hdulist = pyfits.HDUList([hdu_primary, hdu_magbins,
                              hdu_sizebins, hdu_f_unclass,
                              hdu_f_misclass, hdu_avcorr, hdu_stdcorr])
    outfile = data_path+'dr6z7_zoo_conf_app_bin_%s_grids.fits'%label
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)

def apply_gz_conf_app_bins(label='greater'):
    pinfile = data_path+'dr6_photometry_all.fits'
    data = pyfits.getdata(pinfile)
    mag = data.field('R')
    mag_delta = 0.125
    size = data.field('petroR50_r')
    size_delta = 0.25
    f_unclass = N.zeros(len(data), N.float) + 1.1
    f_misclass = N.zeros(len(data), N.float) + 1.1
    avcorr = N.zeros(len(data), N.float) + 1.1
    stdcorr = N.zeros(len(data), N.float) + 1.1
    outfile = data_path+'dr6z7_zoo_conf_app_bin_%s_grids.fits'%label
    grids = pyfits.open(outfile)
    magbins = grids['magbins'].data
    sizebins = grids['sizebins'].data
    f_unclass_grid = grids['f_unclass'].data
    f_misclass_grid = grids['f_misclass'].data
    avcorr_grid = grids['avcorr'].data
    stdcorr_grid = grids['stdcorr'].data
    for imagi, imag in enumerate(magbins):
        mag_selection = (mag >= imag) & (mag < (imag+mag_delta))
        if not N.any(mag_selection):  continue
        print 'imag =', imag
        for isizei, isize in enumerate(sizebins):
            selection = (mag_selection & (size >= isize) & 
                         (size < (isize+size_delta)))
            if not N.any(selection):  continue
            #d = data[selection]
            #n = float(len(d))
            f_unclass[selection] = f_unclass_grid[imagi, isizei]
            f_misclass[selection] = f_misclass_grid[imagi, isizei]
            avcorr[selection] = avcorr_grid[imagi, isizei]
            stdcorr[selection] = stdcorr_grid[imagi, isizei]
    cols = []
    cols.append(pyfits.Column(name='OBJID',
                              format='K', array=data.field('objid')))
    cols.append(pyfits.Column(name='F_UNCLASS_APP_BIN',
                              format='E', array=f_unclass))
    cols.append(pyfits.Column(name='F_MISCLASS_APP_BIN',
                              format='E', array=f_misclass))
    cols.append(pyfits.Column(name='AVCORR_APP_BIN',
                              format='E', array=avcorr))
    cols.append(pyfits.Column(name='STDCORR_APP_BIN',
                              format='E', array=stdcorr))
    tbhdu=pyfits.new_table(cols)
    tbhdu.name = 'DATA'
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.append(tbhdu)
    outfile = data_path+'dr6_zoo_all_conf_app_bin_%s.fits'%label
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdulist.writeto(outfile)

def dr_table4():
    abs_greater = pyfits.open(data_path+
                              'dr6z7_zoo_confidence_abs_bin_greater.fits')
    abs_clean = pyfits.open(data_path+
                            'dr6z7_zoo_confidence_abs_bin_clean.fits')
    print 'OK:', (abs_clean[1].data.OBJID == abs_greater[1].data.OBJID).all()
    cols = [abs_clean[1].columns[0]]
    for c in abs_clean[1].columns[1:]:
        cc = c.copy()
        cc.name = 'CLEAN_'+cc.name
        cols.append(cc)
    for c in abs_greater[1].columns[1:]:
        cc = c.copy()
        cc.name = 'GREATER_'+cc.name
        cols.append(cc)
    hdu = pyfits.new_table(cols)
    outfile = data_path+'dr6z7_zoo_for_dr_table4.fits'
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdu.writeto(outfile)

def dr_table8():
    app_greater = pyfits.open(data_path+
                              'dr6_zoo_all_conf_app_bin_greater.fits')
    app_clean = pyfits.open(data_path+
                            'dr6_zoo_all_conf_app_bin_clean.fits')
    print 'OK:', (app_clean[1].data.OBJID == app_greater[1].data.OBJID).all()
    cols = [app_clean[1].columns[0]]
    for c in app_clean[1].columns[1:]:
        cc = c.copy()
        cc.name = 'CLEAN_'+cc.name
        cols.append(cc)
    for c in app_greater[1].columns[1:]:
        cc = c.copy()
        cc.name = 'GREATER_'+cc.name
        cols.append(cc)
    hdu = pyfits.new_table(cols)
    outfile = data_path+'dr6z7_zoo_for_dr_table8.fits'
    file_exists = os.path.isfile(outfile)
    if file_exists:
	os.remove(outfile)
    hdu.writeto(outfile)

def dr_tables(create=False):
    import fits2csv
    if create:
        dr_table2()
        #dr_table2(True)
        dr_table3()
        #dr_table3(True)
        dr_table4()
    for i in ('2', '3', '4'):
        fits2csv.fits2latex(data_path+'dr6z7_zoo_for_dr_table%s.fits'%i,
                     data_path+'../data_release/dr6z7_zoo_for_dr_table%s.tex'%i)
        fits2csv.fits2csv_round(data_path+'dr6z7_zoo_for_dr_table%s.fits'%i,
            data_path+'../data_release/dr6z7_zoo_for_dr_table%s.csv'%i, round=3)

def do_twiddle(d, twiddle):
    # shift grid by twiddle bins
    for i in range(len(d)):
        di = d[i]
        n = di.shape[0]
        interp = interp1d(N.arange(n), di, axis=0, copy=True,
                          bounds_error=False, fill_value=0)
        di[:] = interp(N.arange(n)+twiddle[i])
    return d

def do_twiddle_old(d, twiddle):
    # shift grid by twiddle bins
    if twiddle > 0:
        d[:,0:-twiddle,:] = N.array(d[:,twiddle:,:])
        d[:,-twiddle:,:] = 0
    else:
        twiddle = -twiddle
        d[:,twiddle:,:] = N.array(d[:,:-twiddle,:])
        d[:,:twiddle,:] = 0
    return d
