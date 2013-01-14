import gc
import numpy as np
import numpy.ma as ma
import cPickle as pickle
import time
import warnings
import cosmology
import galaxyzoo2 as gz2

from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.stats import scoreatpercentile
from scipy.optimize import fmin_powell
from scipy.signal import convolve2d
from astropy.io import fits as pyfits
from astroML.plotting import hist

gz_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/'
fits_path = gz_path+'fits/'
plots_path = gz_path+'plots/'
gz2table_data_file = fits_path+'gz2table.fits'
gz2both_data_file = fits_path+'gz2main_table_sample.fits'

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
prob_err = 0.01

gc.collect()
plt.ion()

def weighted_bulge_prominence(prob_arr):
    n_responses = 4
    if len(prob_arr) != n_responses:
        raise ValueError('Array must be of length %s' % n_responses)

    # Check to make sure probabilities are reasonable. Slightly higher than one since
    # values may be truncated due to rounding in data conversion

    if np.sum(prob_arr) > 1. + prob_err:
        warnings.warn('Caution: probabilities in array sum to greater than 1')
        print np.sum(prob_arr)

    weights = np.arange(n_responses)/np.float(n_responses - 1)
    wbp = np.sum(weights * prob_arr)
    return wbp

def compute_wbp(data):
    p0 = data['t05_bulge_prominence_a10_no_bulge_weighted_fraction']
    p1 = data['t05_bulge_prominence_a11_just_noticeable_weighted_fraction']
    p2 = data['t05_bulge_prominence_a12_obvious_weighted_fraction']
    p3 = data['t05_bulge_prominence_a13_dominant_weighted_fraction']
    wbp_val = np.zeros(len(data),dtype='float')
    for index,(arr) in enumerate(zip(p0,p1,p2,p3)):
        wbp_val[index] = weighted_bulge_prominence(arr)

    return wbp_val

def create_binned_data(min_classifications=20,
                       #min_prob=0.8,
                       zmin = 0.01, zmax = 0.16, zstep = 0.01,
                       magmin = -24, magmax = -16, mag_numbins = 20,
                       sizemin = 0, sizemax = 15, size_numbins = 20,
                       plot=False
                       ):

    p = pyfits.open(gz2both_data_file)

    gzdata = p[1].data

    # Find galaxies with sufficient classifications on bulge prominence question

    taskcount = gzdata['t05_bulge_prominence_total_weight'] > min_classifications
    bulge = gzdata[taskcount]

    print ''
    print '%s bulge galaxies above threshold'%len(bulge)
    print ''

    bulgestat = compute_wbp(bulge)

    lobulge = bulge[bulgestat < 0.5]
    hibulge = bulge[bulgestat >= 0.5]

    datalist = [lobulge, hibulge]
    namelist = ['lobulge','hibulge']

    if plot:
        fig = plt.figure(1,(14,8))
        fig.clf()

    for index, (mdata, mname) in enumerate(zip(datalist,namelist)):

        redshift = mdata['REDSHIFT']
        mr = mdata['PETROMAG_MR']
        r50_kpc = mdata['PETROR50_R_KPC']
        r90 = mdata['PETROR90_R']

        # Mask the data where no values exist

        redshift_mask = np.isfinite(redshift)

        # Add masks for magnitude, size, and surface brightness limits

        dmod = cosmology.dmod_flat(redshift[redshift_mask])
        ang_scale = cosmology.ang_scale_flat(redshift[redshift_mask])
        absmag_lim = 17.77 - dmod
        surfacebrightness_app_lim = 23.0

        r90_kpc = r90[redshift_mask] * ang_scale

        magnitude_mask = mr[redshift_mask] > (absmag_lim - 1.0)
        size_mask = r90_kpc < (3.0 * ang_scale)
        surfacebrightness_mask = (mr[redshift_mask] + dmod + 2.5*np.log10(6.283185*(r50_kpc[redshift_mask]/ang_scale)**2)) > (surfacebrightness_app_lim - 1.0)
        totalmask = magnitude_mask | size_mask | surfacebrightness_mask

        print '%s galaxies removed from %s sample due to magnitude cutoff'%(np.sum(magnitude_mask.astype(int)),mname)
        print '%s galaxies removed from %s sample due to angular size cutoff'%(np.sum(size_mask.astype(int)),mname)
        print '%s galaxies removed from %s sample due to surface brightness cutoff'%(np.sum(surfacebrightness_mask.astype(int)),mname)
        print '%s percent of the total'%(np.sum(totalmask.astype(float))/len(redshift_mask) * 100)
        print ' '

        z_masked = redshift[redshift_mask]
        m_masked = mr[redshift_mask]
        r_masked = r50_kpc[redshift_mask]

        zbins = np.arange(zmin, zmax, zstep)
        magbins = np.linspace(magmin, magmax, mag_numbins)
        sizebins = np.linspace(sizemin, sizemax, size_numbins)
        h, edges = np.histogramdd((z_masked, m_masked, r_masked),bins=(zbins,magbins,sizebins))

        if plot:
            ax = fig.add_subplot(1,2,index)
            imarr = np.squeeze(h[3:4,:,:])                # Pick an example redshift range to slice on
            im = ax.imshow(imarr, extent=(magmin,magmax,sizemin,sizemax), interpolation="nearest", origin='lower')
            ax.set_title(mname)
            ax.set_xlabel(r'$M_R [mag]$')
            ax.set_ylabel(r'$R_{50} [kpc]$')
            ax.set_aspect('auto')
            cb = plt.colorbar(im)

            ax.set_xlim(min(magbins),max(magbins))
            ax.set_ylim(min(sizebins),max(sizebins))

        edges_redshift = edges[0]
        edges_mag = edges[1]
        edges_size = edges[2]

        centers_redshift = edges_redshift[:-1] + (edges_redshift[1]-edges_redshift[0])/2.
        centers_mag = edges_mag[:-1] + (edges_mag[1]-edges_mag[0])/2.
        centers_size = edges_size[:-1] + (edges_size[1]-edges_size[0])/2.

        col1_edge = pyfits.Column(name = 'edges', format='E', array=edges_redshift)
        col2_edge = pyfits.Column(name = 'edges', format='E', array=edges_mag)
        col3_edge = pyfits.Column(name = 'edges', format='E', array=edges_size)

        col1_centers = pyfits.Column(name = 'centers', format='E', array=centers_redshift)
        col2_centers = pyfits.Column(name = 'centers', format='E', array=centers_mag)
        col3_centers = pyfits.Column(name = 'centers', format='E', array=centers_size)

        primary_hdu = pyfits.PrimaryHDU(h)
        hdulist = pyfits.HDUList([primary_hdu])

        tb1_hdu = pyfits.new_table(pyfits.ColDefs([col1_edge]))
        tb1_hdu.name = 'REDSHIFT_BIN_EDGES'
        tb2_hdu = pyfits.new_table(pyfits.ColDefs([col2_edge]))
        tb2_hdu.name = 'MR_BIN_EDGES'
        tb3_hdu = pyfits.new_table(pyfits.ColDefs([col3_edge]))
        tb3_hdu.name = 'R50_KPC_BIN_EDGES'
        tb4_hdu = pyfits.new_table(pyfits.ColDefs([col1_centers]))
        tb4_hdu.name = 'REDSHIFT_BIN_CENTERS'
        tb5_hdu = pyfits.new_table(pyfits.ColDefs([col2_centers]))
        tb5_hdu.name = 'MR_BIN_CENTERS'
        tb6_hdu = pyfits.new_table(pyfits.ColDefs([col3_centers]))
        tb6_hdu.name = 'R50_KPC_BIN_CENTERS'

        hdulist.append(tb1_hdu)
        hdulist.append(tb2_hdu)
        hdulist.append(tb3_hdu)
        hdulist.append(tb4_hdu)
        hdulist.append(tb5_hdu)
        hdulist.append(tb6_hdu)

        hdulist.writeto(fits_path+'%s_'%mname+'binned_counts.fits',clobber=True)    

    p.close()

    return None

def determine_ratio_baseline(min_galperbin=50,linear=False):

    # Load the GZ2 sample data

    p_lobulge = pyfits.open(fits_path+'lobulge_binned_counts.fits')
    p_hibulge = pyfits.open(fits_path+'hibulge_binned_counts.fits')

    d_lobulge = p_lobulge[0].data.astype(int)
    d_hibulge = p_hibulge[0].data.astype(int)

    # Bin sizes are set when FITS file is created; same for both sets

    zbins = p_lobulge['REDSHIFT_BIN_CENTERS'].data['centers']
    magbins = p_lobulge['MR_BIN_CENTERS'].data['centers']
    sizebins = p_lobulge['R50_KPC_BIN_CENTERS'].data['centers']

    n_magbin = len(magbins)
    n_sizebin = len(sizebins)

    # Empty 2-D arrays to store the ratio and number of galaxies for elliptical/spiral ratio

    ratio_baseline = np.zeros((n_magbin, n_sizebin), np.float) + unknown_ratio
    counts_baseline = np.zeros((n_magbin, n_sizebin, 2), np.int) + unknown_ratio

    # Trim the data so it only includes counts in the relevant magnitude and size bins

    lobulgez = d_lobulge
    hibulgez = d_hibulge

    # Loop over each slice in redshift space. Cells without entries will be filled in the first available loop.

    for z_index, zsel in enumerate(zbins[1:]):
        lobulge = lobulgez[z_index, :, :]
        hibulge = hibulgez[z_index, :, :]
     
        # Create mask for cells with low counts
        mask = (lobulge + hibulge) < min_galperbin
     
        if linear:
            ratio = hibulge.astype(float)/lobulge
        else:
            ratio = np.log10(hibulge.astype(float)/lobulge)

        # Mask galaxies outside of the selection limits

        np.putmask(ratio, mask, unknown_ratio)        # If ratio is masked or infinite, replace it with the value ``unknown ratio''
        np.putmask(ratio, np.isinf(ratio), unknown_ratio)      
        select = np.logical_not(mask)                # Invert the mask so True = good values in select

        empty = ratio_baseline <= unknown_ratio+eps        # Second mask; "True" is where master array has no values yet

        # Combine empty and select masks

        select &= empty                                # Union is where master array has no values AND ratio is determined. 

        # Populate the 2-D empty arrays with the el/sp ratio for all non-masked cells

        ratio_baseline[select] = ratio[select]
        counts_baseline[select] = np.transpose([lobulge[select], hibulge[select]])

    ratio_baseline_masked = ma.masked_less_equal(ratio_baseline, unknown_ratio)
    counts_baseline_masked = ma.masked_less_equal(counts_baseline, unknown_ratio)

    # Write the results to FITS files

    pyfits.writeto(fits_path+'local_ratio_baseline_bulge.fits',
                   ratio_baseline, clobber=True)    
    pyfits.writeto(fits_path+'local_counts_baseline_bulge.fits',
                   counts_baseline, clobber=True)    

    pickle.dump(ratio_baseline_masked, open(gz_path+'local_ratio_baseline_bulge_masked.pkl','wb')) 
    pickle.dump(counts_baseline_masked, open(gz_path+'local_counts_baseline_bulge_masked.pkl','wb')) 

    # Close PyFITS objects

    p_lobulge.close()
    p_hibulge.close()

    return None

def determine_ratio_baseline_sigma(plot=False):
    data = pyfits.getdata(fits_path+'local_counts_baseline_bulge.fits')
    el = data[:,:,0].astype(np.float)
    sp = data[:,:,1].astype(np.float)
    mask_el = el < 1
    mask_sp = sp < 1
    mask_all = np.logical_and(mask_el, mask_sp)
    mask = np.logical_or(mask_el, mask_sp)
    ratio = pyfits.getdata(fits_path+'local_ratio_baseline_bulge.fits')
    count_sum = (el + sp).astype(np.float)
    count_product = (el * sp).astype(np.float)
    np.putmask(count_product, mask, 1.0)
    np.putmask(count_sum, mask, 1.0)
    sigma = (np.log10(np.e))**2 * count_sum / count_product
    sigma = np.sqrt(sigma)
    np.putmask(sigma, mask, unknown_ratio)

    sigma_masked = ma.masked_less_equal(sigma, unknown_ratio)

    if plot:
        fig = plt.figure(2)
        fig.clf()
        ax = fig.add_subplot(111)
        im = ax.imshow(sigma_masked, interpolation='nearest', origin='lower')
        cb = plt.colorbar(im)
        ax.set_aspect('auto')
        ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
        ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)

    pyfits.writeto(fits_path+'local_ratio_baseline_bulge_sigma.fits',
           sigma, clobber=True)    

    return None

def plot_ratio_baseline(plot_mag_lims=False):

    f_ratio = fits_path+'local_ratio_baseline_bulge.fits'
    f_counts = fits_path+'lobulge_binned_counts.fits'
    p_ratio = pyfits.open(f_ratio)
    p_counts = pyfits.open(f_counts)

    zbins = p_counts['REDSHIFT_BIN_CENTERS'].data['centers']
    magbins = p_counts['MR_BIN_CENTERS'].data['centers']
    sizebins = p_counts['R50_KPC_BIN_CENTERS'].data['centers']

    zedges = p_counts['REDSHIFT_BIN_EDGES'].data['edges']
    magedges = p_counts['MR_BIN_EDGES'].data['edges']
    sizeedges = p_counts['R50_KPC_BIN_EDGES'].data['edges']

    ratio_baseline_masked = pickle.load(open(gz_path+'local_ratio_baseline_bulge_masked.pkl','rb'))
    counts_baseline_masked = pickle.load(open(gz_path+'local_counts_baseline_bulge_masked.pkl','rb'))

    sumcounts = np.sum(counts_baseline_masked, axis=2)

    titlenames = ('counts','ratio')
    labelnames = (r'$N_{lobulge} + N_{hibulge}$',r'$log_{10}(N_{lobulge}/N_{hibulge})$')

    fig = plt.figure(3,(14,8))
    fig.clf()
    vrange = [(np.min(sumcounts),np.max(sumcounts)),(min_ratio,max_ratio)]

    for index, data in enumerate((sumcounts, ratio_baseline_masked)):
        ax = fig.add_subplot(1,2,index+1,aspect=1)
        cmap = cm.jet
        cmap.set_bad('k')
        #imextent=(magbins[0],magbins[-1],sizebins[0],sizebins[-1])
        imextent=(magedges[0],magedges[-1],sizeedges[0],sizeedges[-1])
        im = ax.imshow(data, 
                       vmin = vrange[index][0],
                       vmax = vrange[index][1],
                       extent=imextent,
                       interpolation='nearest',origin='lower')
        ax.set_title('hibulge/lobulge '+titlenames[index])
        ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
        ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        cb = plt.colorbar(im,orientation='vertical')
        cb.set_label(labelnames[index],fontsize=16)

        SBapp = 23.0
        appmag_lim = 17.0
        SBlim_size = np.arange(200) / 10.0
        SBlim_mag = (SBapp - cosmology.dmod_flat(np.mean(zbins))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(zbins)))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(zbins))
        absmag_lim_loz = appmag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = appmag_lim - cosmology.dmod_flat(0.25)
        r50_lim = cosmology.ang_scale_flat(np.mean(zbins))
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        if plot_mag_lims:
            ax.plot(np.zeros(2)+absmag_lim, np.array([sizeedges[0],sizeedges[-1]]),'w--')
            ax.plot(np.zeros(2)+absmag_lim_loz, np.array([sizeedges[0],sizeedges[-1]]),'g--')
            ax.plot(np.zeros(2)+absmag_lim_hiz, np.array([sizeedges[0],sizeedges[-1]]),'g--')

        ax.plot(np.array([magedges[0],magedges[-1]]), np.zeros(2)+r50_lim,'w--')
        plt.show()
        plt.draw()

    p_ratio.close()
    p_counts.close()

    return None

def bootstrap_fmin(f, p0, x, y, z, mask, nboot=50):
    plist = np.zeros((nboot, len(p0)), np.float)                                # Zero array for number of tries, number of parameters
    ndata = len(mask.ravel().nonzero()[0])                                      # number of non-masked cells in the ratio data
    for i in range(nboot):
        bootmask = np.zeros(mask.shape)                                         # Grid of zeroes same size as mask, data
        while bootmask.sum() < ndata:                                           # Loop until grid is filled at data locations
            rx = int(np.random.uniform(0, x.shape))                             # Pick random indices from the bin arrays
            ry = int(np.random.uniform(0, y.shape))                             # 
            if mask[rx,ry] > 0.000001:                                          # If existing data mask in that cell is non-zero
                bootmask[rx,ry] += 1                                            #   then put a 1 in bootmask cell
                                                                                # Preserves total weight, but assigns random weights
                                                                                #   to each data cell. Perturbs to get multiple guesses. 
        bootmask *= mask                                                        # Multiply new weights by original mask
        plist[i] = fmin_powell(f, p0, (x, y, z, bootmask),disp=0)                      # Fit parameters based on new mask passed to function
    p = fmin_powell(f, p0, (x, y, z, mask))                                     # Fit parameters based on original (unweighted) mask
    pmed = np.median(plist, axis=0)                                             # Find the median best fit parameters for each variation in the list
    perr = np.array([(scoreatpercentile(plist[:,k], 84.0) -                     # The error is the mean of the 84th percentile parameter
            scoreatpercentile(plist[:,k], 16.0))/2.0 for k in range(len(p))])   #   and the 16th percentile parameter. Ah. Looping is only to estimate error.
    return p, perr

def ratio_function(p, x, y):
    a, b, c, d, e, f, g, h, i = p
    z = np.zeros((len(x), len(y)), np.float)
    x0 = b**-(a + h*y**i) + c
    x1 = d + e*(x0-c)
    for i in range(len(x)):
        z[i] = f / (1.0 + np.exp((x0 - x[i])/x1)) + g
    return z

def ratio_minfunc(p, x, y, z, w):
    f = ratio_function(p, x, y)
    r2 = (z - f)**2                                     # Difference between model and data
    df = (w > 0).ravel().astype(np.int).sum() - len(p)  # Degrees of freedom: N_cells - N_parameters
    s = (w * r2).sum() / df                             # Chi-squared statistic
    return s

def fit_ratio_baseline(nboot=50, plot=False):

    f_ratio = fits_path+'local_ratio_baseline_bulge.fits'
    p_ratio = pyfits.open(f_ratio)
    f_counts = fits_path+'lobulge_binned_counts.fits'
    p_counts = pyfits.open(f_counts)

    zbins = p_counts['REDSHIFT_BIN_CENTERS'].data['centers']
    magbins = p_counts['MR_BIN_CENTERS'].data['centers']
    sizebins = p_counts['R50_KPC_BIN_CENTERS'].data['centers']

    zedges = p_counts['REDSHIFT_BIN_EDGES'].data['edges']
    magedges = p_counts['MR_BIN_EDGES'].data['edges']
    sizeedges = p_counts['R50_KPC_BIN_EDGES'].data['edges']

    sigma = pyfits.getdata(fits_path+'local_ratio_baseline_bulge_sigma.fits')
    ratio = pyfits.getdata(fits_path+'local_ratio_baseline_bulge.fits')
    
    # Steven's and my initial parameters
    pinit_sb = np.array([-0.3, 1.7, -22.8, 0.3, 0.3, -2.6, 1.5, 1.0, 1.0])
    pinit_kw = np.array([-3.5, 1.9, -22.8, 0.3, 0.3, -4.4, 2.2, 1.1, 1.0])
    pinit = pinit_kw

    chi2r = 0.0
    niter = 0
    sys = 0.0

    weights = 1.0/(sigma**2 + sys**2)
    # Equal weights for ALL bins (empty + counts) - unweighted.
    weights = weights.ravel()
    weights = 0.0 * sigma + weights.compress(weights > 0.000001).mean()

    # Mask with equal weights for all the bins with counts in them
    np.putmask(weights, sigma > 98, 0.0)
    np.putmask(weights, sigma < -98, 0.0)

    # Try normalizing the weights

    weights /= weights.max()

    p, perr = bootstrap_fmin(ratio_minfunc, pinit,
                             magbins, sizebins, ratio, weights, nboot)
    chi2r = ratio_minfunc(p, magbins, sizebins, ratio, weights)

    print 'param =', p
    print 'param errors =', perr
    print 'param errors (%) =', perr/p
    print 'chi squared =', chi2r

    ratio_baseline_fit = ratio_function(p, magbins, sizebins)
    pickle.dump((p, perr), file(fits_path+'ratio_baseline_bulge_fit.pickle', 'w'))
    res = ratio - ratio_baseline_fit
    normres = res/sigma
    np.putmask(normres, weights < 0.000001, unknown_ratio)

    if plot:
        gz2.plot_ratio_function(p,magbins,sizebins)

    pyfits.writeto(fits_path+'local_ratio_baseline_bulge_fit.fits',
                   ratio_baseline_fit, clobber=True)    
    pyfits.writeto(fits_path+'local_ratio_baseline_bulge_fitres.fits',
                   res, clobber=True) 
    pyfits.writeto(fits_path+'local_ratio_baseline_bulge_fitnormres.fits',
                   normres, clobber=True) 

    p_ratio.close()
    p_counts.close()

    return None

def plot_ratio_baseline_fit(fitbins=19, kernel_size = 3, plot_contour=False):

    ratio_baseline_fit = pyfits.getdata(fits_path+'local_ratio_baseline_bulge_fit.fits')

    f_counts = fits_path+'lobulge_binned_counts.fits'
    p_counts = pyfits.open(f_counts)

    zbins = p_counts['REDSHIFT_BIN_CENTERS'].data['centers']
    magbins = p_counts['MR_BIN_CENTERS'].data['centers']
    sizebins = p_counts['R50_KPC_BIN_CENTERS'].data['centers']

    zedges = p_counts['REDSHIFT_BIN_EDGES'].data['edges']
    magedges = p_counts['MR_BIN_EDGES'].data['edges']
    sizeedges = p_counts['R50_KPC_BIN_EDGES'].data['edges']

    ratio_baseline_masked = pickle.load(open(gz_path+'local_ratio_baseline_bulge_masked.pkl','rb'))
    ratio_baseline_fit_masked = ma.array(ratio_baseline_fit,mask=ratio_baseline_masked.mask)

    kernel = np.ones((kernel_size,kernel_size),dtype=int)
    convarr = np.logical_not(ratio_baseline_masked.mask)
    bw = convolve2d(np.logical_not(ratio_baseline_masked.mask).astype(int),
                    kernel,mode='same').astype(bool).astype(int) * -1

    fig = plt.figure(5)
    fig.clf()
    ax = fig.add_subplot(111)

    # Plot the best fit function as a semitransparent layer

    pfit, pfit_err = pickle.load(file(fits_path+'ratio_baseline_bulge_fit.pickle', 'r'))
    magbinsplot = np.linspace(magbins[0],magbins[-1],fitbins)
    sizebinsplot = np.linspace(sizebins[0],sizebins[-1],fitbins)
    fitarray_fn = ratio_function(pfit, magbinsplot, sizebinsplot)

    fit_extent=(magedges[0],magedges[-1],sizeedges[0],sizeedges[-1])
    imf = ax.imshow(fitarray_fn,
                   alpha=1.0,
                   extent = fit_extent,
                   vmin = min_ratio, vmax = max_ratio, 
                   interpolation='nearest', origin='lower')
    cb = plt.colorbar(imf)
    cb.set_label(r'log$_{10}(n_{el}/n_{sp})$',fontsize=20)

    # Plot the local masked baseline relation as opaque layer

    masked_extent = fit_extent
    cmap = cm.jet
    cmap.set_bad(alpha=0.0)
    im = ax.imshow(ratio_baseline_masked, 
                   alpha=0.5,
                   extent = masked_extent,
                   vmin = min_ratio, vmax=max_ratio,
                   interpolation='nearest',origin='lower')

    # Plot the outline of the local region. 
        # This plots the contour, but not in the right region (or quite the right size)
    
    if plot_contour:
        cont = plt.contour(magbins, sizebins, bw, 
                           levels = [-1], 
                           color='w', 
                           linewidths=3, linestyles='dashed')

        zc = cont.collections[0]
        zc.set_color('w')
 
    SBapp = 23.0
    appmag_lim = 17.0
    SBlim_size = np.arange(200) / 10.0
    SBlim_mag = (SBapp - cosmology.dmod_flat(np.mean(zbins))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(zbins)))**2))
    absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(zbins))
    r50_lim = cosmology.ang_scale_flat(np.mean(zbins))
    ax.autoscale(False)
    ax.plot(SBlim_mag, SBlim_size,'w--')
    ax.plot(np.zeros(2)+absmag_lim, np.array([sizeedges[0],sizeedges[-1]]),'w--')
    ax.plot(np.array([magedges[0],magedges[-1]]), np.zeros(2)+r50_lim,'w--')

    # Set general axes properties

    ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
    ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
    ax.set_aspect('auto')
    rc(('xtick','ytick'), labelsize=12)

    plt.show()

    return None

def do_all(min_galperbin=25, nboot=1, min_classifications=15, fitbins = 1000):
    
    tstart = time.time()

    create_binned_data(min_classifications)
    determine_ratio_baseline(min_galperbin,linear=False)
    plot_ratio_baseline()
    determine_ratio_baseline_sigma()
    fit_ratio_baseline(nboot,plot=True)
    plot_ratio_baseline_fit(fitbins)

    tend = time.time()
    print 'Time elapsed for do_all: %i seconds' % (tend - tstart)


