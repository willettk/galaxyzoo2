"""

This code is intended to reduce data from Galaxy Zoo 2, and to 
reproduce the tables, figures, and data in Willett et al. (in prep). It
relies on the aggregate GZ2 data tables created by Steven Bamford, as well
as other metadata retrieved from the SDSS CasJobs server. 

"""

import gc
import cPickle as pickle
import numpy as np
import numpy.ma as ma
import time
import subprocess
import os
import warnings

from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
from random import sample
from itertools import chain
from scipy.stats import scoreatpercentile
from scipy.stats import tsem
from scipy.optimize import fmin_powell
from scipy.signal import convolve2d
from astropy.io import fits as pyfits
from astroML.plotting import hist as histML

# From Steven Bamford's Python modules
import cosmology

# Set paths for I/O
gz_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/'
fits_path_main = gz_path+'fits/'
fits_path_task = fits_path_main+'tasks/'
pkl_path = gz_path+'pickle/'
plots_path = gz_path+'plots/'
dropbox_figs_path = gz_path+'gz2dropbox/figures/'

# Locations of input FITS data files
gz2table_data_file = fits_path_main+'gz2table.fits'
gz2_photoz_data_file = fits_path_main+'gz2_photoz_table_sample.fits'
gz2_coadd2_data_file = fits_path_main+'gz2_coadd2_table_sample.fits'
gz2_maglim_data_file = fits_path_main+'gz2_original_extra_s82norm_r17_table_sample.fits'
gz2_full_data_file = fits_path_main+'gz2_original_extra_s82norm_table_sample.fits'
gz2_both_data_file = fits_path_main+'gz2main_table_sample.fits'
gz2_stripe82_data_file = fits_path_main+'gz2_stripe82_normal.fits'

# Global parameters for data reduction
min_ratio = -2.0
max_ratio = 2.0
range_ratio = max_ratio - min_ratio
unknown_ratio = -999
no_correction_applied = unknown_ratio + 10
eps = 0.000001
SBlim_app = 23.0
appmag_lim = 17.0
SBlim_size = np.arange(200) / 10.0

# Collect memory and run plotting interactively
gc.collect()
plt.ion()

def bin_data_idl(task_dict, zmin = 0.00, zmax = 0.26, zstep = 0.01,
                 magmin = -24, magmax = -16, magstep = 0.25,
                 sizemin = 0, sizemax = 15, sizestep = 0.5):

    """ Bins the galaxies in three dimensions: absolute magnitude,
        physical size, and redshift. Only bins galaxies with
        spectroscopic redshifts. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    zmin : float or integer
        Minimum redshift value for binning galaxies

    zmax : float or integer
        Maximum redshift value for binning galaxies

    zstep : float or integer
        Redshift interval for binning galaxies

    magmin : float or integer
        Minimum absolute magnitude value for binning galaxies

    magmax : float or integer
        Maximum absolute magnitude value for binning galaxies

    magstep : float or integer
        Absolute magnitude interval for binning galaxies

    sizemin : float or integer
        Minimum size [kpc] value for binning galaxies

    sizemax : float or integer
        Maximum size [kpc] value for binning galaxies

    sizestep : float or integer
        Size interval [kpc] for binning galaxies

    Returns
    -------
    None

    Notes
    -------
    Calls the IDL function `bin_gz2_from_python.pro`. Uses IDL since
        it was unclear how to bin in three dimensions while tracking
        indexes of metadata in pure Python. 

    """

    tstart = time.time()
    p = pyfits.open(gz2_maglim_data_file)
    gzdata = p[1].data
    p.close()

    gzdata_hasredshift = gzdata[np.isfinite(gzdata['REDSHIFT'])]

    redshift = gzdata_hasredshift['REDSHIFT']
    mr = gzdata_hasredshift['PETROMAG_MR']
    r50_kpc = gzdata_hasredshift['PETROR50_R_KPC']
    r90 = gzdata_hasredshift['PETROR90_R']

    # Add masks for magnitude, size, and surface brightness limits. Remember that these are masks, 
    # so values of True are NOT good galaxies (hence the logical_not later on).  

    dmod = cosmology.dmod_flat(redshift)
    ang_scale = cosmology.ang_scale_flat(redshift)
    absmag_lim = appmag_lim - dmod
    r90_kpc = r90 * ang_scale

    absmag_padding = 0.0
    sb_padding = 0.0

    taskcount_mask = gzdata_hasredshift[task_dict['task_name_count']] < task_dict['min_classifications']
    magnitude_mask = mr > (absmag_lim - absmag_padding)
    size_mask = r90_kpc < (3.0 * ang_scale)
    s82_maglim_mask = (gzdata_hasredshift['PETROMAG_R'] - gzdata_hasredshift['EXTINCTION_R']) > appmag_lim

    """
    Removed the surface brightness requirement - no selection algorithms for the GZ2 sample should depend
    on this, and I'm worried that the signs in the equation below are not correct. 

    surfacebrightness_mask = (mr + dmod + \
        2.5*np.log10(6.283185*(r50_kpc/ang_scale)**2)) > (SBlim_app - sb_padding)
    """

    if task_dict['var_def'] not in ('task01','task06'):
        if isinstance(task_dict['dependent_tasks'],tuple):
            nprev = len((task_dict['dependent_tasks']))
            arr = np.reshape(np.ravel([gzdata_hasredshift[task] for task in task_dict['dependent_tasks']]),(nprev,len(gzdata_hasredshift)))
            probarr = np.product(arr,axis=0)
        elif isinstance(task_dict['dependent_tasks'],str):
            nprev = 1
            probarr = gzdata_hasredshift[task_dict['dependent_tasks']]

        taskprob_mask = (probarr <= 0.5**nprev)
    else:
        taskprob_mask = np.zeros(len(gzdata_hasredshift),bool)

    totalmask = magnitude_mask | size_mask | taskcount_mask | taskprob_mask | s82_maglim_mask

    totalmask = magnitude_mask | size_mask | taskcount_mask | taskprob_mask | s82_maglim_mask

    print ' '
    print '%7i galaxies removed from sample due to absolute magnitude cutoff'%(np.sum(magnitude_mask.astype(int)))
    print '%7i galaxies removed from sample due to apparent magnitude cutoff'%(np.sum(s82_maglim_mask.astype(int)))
    print '%7i galaxies removed from sample due to angular size cutoff'%(np.sum(size_mask.astype(int)))
    print '%7i galaxies removed from sample due to minimum number (%i) of classifications'%(np.sum(taskcount_mask.astype(int)),task_dict['min_classifications'])
    #print '%7i galaxies removed from sample due to surface brightness cutoff'%(np.sum(surfacebrightness_mask.astype(int)))
    print '%7i galaxies removed from sample due to task probability cutoff'%(np.sum(taskprob_mask.astype(int)))
    print '%6.2f percent of the total (%i galaxies) is kept for %s' % ((1. - np.sum(totalmask.astype(float))/len(gzdata_hasredshift)) * 100., np.sum(np.logical_not(totalmask).astype(int)),task_dict['var_def'])
    print ' '

    gzgood = gzdata_hasredshift[np.logical_not(totalmask)]

    # Write the masked GZ2 data to a FITS file to be binned by IDL script

    pyfits.writeto(fits_path_task+'%s_data_for_idl.fits' % task_dict['var_def'],
                   gzgood, clobber=True)    

    # Open both the root file and a temporary file to contain variables + IDL script

    root_script_file_name = 'bin_gz2_from_python.pro'
    new_script_file_name = 'temp_'+root_script_file_name
    gzidl_path = gz_path + 'idl/'
    f1 = open(gzidl_path+root_script_file_name,'r')
    f2 = open(gzidl_path+ new_script_file_name,'w')

    # Pass variables that need to be set in the IDL script

    f2.write('zmin = %6.3f \n' % zmin)
    f2.write('zmax = %6.3f \n' % zmax)
    f2.write('zstep = %6.3f \n' % zstep)
    f2.write('magmin = %6.3f \n' % magmin)
    f2.write('magmax = %6.3f \n' % magmax)
    f2.write('magstep = %6.3f \n' % magstep)
    f2.write('sizemin = %6.3f \n' % sizemin)
    f2.write('sizemax = %6.3f \n' % sizemax)
    f2.write('sizestep = %6.3f \n' % sizestep)
    f2.write('task_names_wf = "%s" \n' % ', '.join(task_dict['task_names_wf']))
    f2.write('var_def = "%s" \n' % task_dict['var_def'])

    # Insert the rest of the IDL script
    f2.write(f1.read())

    f1.close()
    f2.close()

    # Call the IDL script from Python

    print "Calling the IDL script"
    subp = subprocess.Popen("idl -quiet -e '.run %s'" % new_script_file_name,
            stderr = subprocess.PIPE, stdout = subprocess.PIPE, shell = True)

    # Make sure the FITS file finishes I/O by checking on the existence of a tiny dummy file

    idlfilecheck = True
    while idlfilecheck:
        idlfilecheck = not os.path.isfile(gz_path+'idlfilecreated')
        time.sleep(0.5)

    os.remove(gz_path+'idlfilecreated')

    # End of timer

    tend = time.time()
    print 'Time elapsed to bin GZ2 likelihood data for %s: %i seconds' % (task_dict['var_def'], tend - tstart)

    return None

def get_bins(task_dict):

    """ Retrive the bin edges and centers from `bin_data_idl`

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    Returns
    -------
    centers_redshift, centers_mag, centers_size: float, float, float
        Centers of the bins in redshift, absolute magnitude, and size
    
    edges_redshift, edges_mag, edges_size center : float, float, float
        Edges of the bins in redshift, absolute magnitude, and size

    """

    idl_binned = pyfits.open(fits_path_task+'%s_idlbinned.fits' % (task_dict['var_def']))
    bindata = idl_binned[1].data

    edges_redshift = np.squeeze(bindata['EDGES_REDSHIFT'])
    edges_mag = np.squeeze(bindata['EDGES_MAG'])
    edges_size = np.squeeze(bindata['EDGES_SIZE'])
    centers_redshift = np.squeeze(bindata['CENTERS_REDSHIFT'])
    centers_mag = np.squeeze(bindata['CENTERS_MAG'])
    centers_size = np.squeeze(bindata['CENTERS_SIZE'])

    return centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size 

def determine_ratio_baseline(task_dict, vartop = 0, varbot = 1):

    """ Determine the morphology ratio baseline for a single
        GZ2 task. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    vartop: int
        Index for the variable in the numerator of the baseline
        morphology ratio. 

    varbot: int
        Index for the variable in the denominator of the baseline
        morphology ratio. 

    Returns
    -------
    None

    """

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load the GZ2 sample data

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var_str[vartop],bintype))
        p_var2 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var_str[varbot],bintype))
        d_var1 = p_var1[0].data.astype(int)                         # Data is pre-binned
        d_var2 = p_var2[0].data.astype(int)

        zbins = p_var1['REDSHIFT_BIN_CENTERS'].data['centers']
        magbins = p_var1['MR_BIN_CENTERS'].data['centers']
        sizebins = p_var1['R50_KPC_BIN_CENTERS'].data['centers']
        
        p_var1.close()
        p_var2.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path_task+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path_task+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        if task_dict['reverse']:
            d_var1 = np.squeeze(d_allvar[varbot,:,:,:])
            d_var2 = np.squeeze(d_allvar[vartop,:,:,:])
        else:
            d_var1 = np.squeeze(d_allvar[vartop,:,:,:])
            d_var2 = np.squeeze(d_allvar[varbot,:,:,:])

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    # Bin sizes are set when FITS file is created

    cubeshape = d_allvar.shape
    n_magbin = cubeshape[2]
    n_sizebin = cubeshape[3]

    # Empty 2-D arrays to store the ratio and number of galaxies for chosen variables

    ratio_baseline = np.zeros((n_magbin, n_sizebin), np.float) + unknown_ratio
    counts_baseline = np.zeros((n_magbin, n_sizebin, 2), np.int) + unknown_ratio
    redshift_baseline = np.zeros((n_magbin, n_sizebin), np.float) + unknown_ratio

    # Loop over each slice in redshift space. Cells without entries will be filled in the first available loop.

    for z_index, zsel in enumerate(centers_redshift):
        var1 = d_var1[z_index, :, :]
        var2 = d_var2[z_index, :, :]
     
        # Create mask for cells with low counts
        mask = (var2 + var1) < task_dict['min_galperbin']
     
        # Compute the elliptical/spiral ratio for the entire array

        if task_dict['ratio_type'] is 'linear':
            ratio = var1.astype(float)/var2
        else:
            ratio = np.log10(var1.astype(float)/var2)

        # Mask galaxies outside of the selection limits

        np.putmask(ratio, mask, unknown_ratio)        # If ratio is masked, replace it with the value ``unknown ratio''
        np.putmask(ratio, np.isinf(ratio), unknown_ratio)      
        select = np.logical_not(mask)                # Invert the mask so True = good values in select
        empty = ratio_baseline <= unknown_ratio + eps        # Second mask; "True" is where master array has no values yet
        
        redshift_this_slice = np.zeros((n_magbin, n_sizebin), np.float) + zsel
        np.putmask(redshift_this_slice, mask, unknown_ratio)
        np.putmask(redshift_this_slice, np.isinf(redshift_this_slice), unknown_ratio)      
        select_redshift = np.logical_not(mask)
        empty_redshift = redshift_baseline <= unknown_ratio + eps        

        """
        Combine empty and select masks 
        Union is where master array has no values AND ratio is determined. 
        """

        select &= empty                               
        select_redshift &= empty_redshift

        # Populate the 2-D empty arrays with the el/sp ratio for all non-masked cells

        ratio_baseline[select] = ratio[select]
        counts_baseline[select] = np.transpose([var1[select], var2[select]])
        redshift_baseline[select_redshift] = redshift_this_slice[select_redshift]

    ratio_baseline_masked = ma.masked_less_equal(ratio_baseline, unknown_ratio)
    counts_baseline_masked = ma.masked_less_equal(counts_baseline, unknown_ratio)
    redshift_baseline_masked = ma.masked_less_equal(redshift_baseline, unknown_ratio)

    # Write the results to FITS files

    pyfits.writeto(fits_path_task+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot),
                   ratio_baseline, clobber=True)    
    pyfits.writeto(fits_path_task+'%s_r%s%s_local_counts_baseline.fits' % (var_def,vartop,varbot),
                   counts_baseline, clobber=True)    
    pyfits.writeto(fits_path_task+'%s_r%s%s_local_redshift_baseline.fits' % (var_def,vartop,varbot),
                   redshift_baseline, clobber=True)    

    pickle.dump(ratio_baseline_masked, open(pkl_path+'%s_r%s%s_local_ratio_baseline_masked.pkl' % (var_def,vartop,varbot),'wb')) 
    pickle.dump(counts_baseline_masked, open(pkl_path+'%s_r%s%s_local_counts_baseline_masked.pkl' % (var_def,vartop,varbot),'wb')) 
    pickle.dump(redshift_baseline_masked, open(pkl_path+'%s_r%s%s_local_redshift_baseline_masked.pkl' % (var_def,vartop,varbot),'wb')) 

    return None

def determine_ratio_baseline_sigma(task_dict, plot=False, vartop=0, varbot=1):


    """ Determine the uncertainty in the morphology ratio baseline 
        for a single GZ2 task. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    plot: bool
        When `True`, plot and save the uncertainty in the morphology
        ratio. 

    vartop: int
        Index for the variable in the numerator of the baseline
        morphology ratio. 

    varbot: int
        Index for the variable in the denominator of the baseline
        morphology ratio. 

    Returns
    -------
    None

    """

    var_def = task_dict['var_def']

    data = pyfits.getdata(fits_path_task+'%s_r%s%s_local_counts_baseline.fits' % (var_def,vartop,varbot))
    var1 = data[:,:,0].astype(np.float)
    var2 = data[:,:,1].astype(np.float)
    mask_var1 = var1 < 1
    mask_var2 = var2 < 1
    mask_all = np.logical_and(mask_var1, mask_var2)
    mask = np.logical_or(mask_var1, mask_var2)
    ratio = pyfits.getdata(fits_path_task+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot))
    count_sum = (var1 + var2).astype(np.float)
    count_product = (var1 * var2).astype(np.float)
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
        im = ax.imshow(sigma_masked.T, interpolation='nearest', origin='lower')
        cb = plt.colorbar(im)
        ax.set_aspect('auto')
        ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
        ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_title('%s baseline ratio sigma' % var_def,fontsize=22)

    pyfits.writeto(fits_path_task+'%s_r%s%s_local_ratio_baseline_sigma.fits' % (var_def,vartop,varbot),
           sigma, clobber=True)    

    return None

def plot_ratio_baseline(task_dict,
                        unset_vrange = False,
                        ratio_vrange = (min_ratio,max_ratio),
                        plot_mag_lims=False,
                        vartop=0, varbot=1
                        ):

    """ Plot the morphology ratio baseline for a single
        GZ2 task. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    unset_vrange : bool
        When `True`, min and max of the colormap are set to the 
        range of the data.

    ratio_vrange: tuple
        Sets the min and max of the colormap. Only applies if
        `unset_vrange` = `False`.

    plot_mag_lims : bool
        When `True`, overplot the default absolute magnitude limits
        of the sample.

    vartop: int
        Index for the variable in the numerator of the baseline
        morphology ratio. 

    varbot: int
        Index for the variable in the denominator of the baseline
        morphology ratio. 

    Returns
    -------
    ax:
        matplotlib Axes object

    """

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    f_ratio = fits_path_task+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot)
    p_ratio = pyfits.open(f_ratio)

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    ratio_baseline_masked = pickle.load(open(pkl_path+'%s_r%s%s_local_ratio_baseline_masked.pkl' % (var_def,vartop,varbot),'rb'))
    counts_baseline_masked = pickle.load(open(pkl_path+'%s_r%s%s_local_counts_baseline_masked.pkl' % (var_def,vartop,varbot),'rb'))
    sumcounts = np.sum(counts_baseline_masked, axis=2)

    if task_dict['ratio_type'] is 'linear':
        label_prefix = ''
    else:
        label_prefix = 'log_{10}'

    titlenames = ('counts','ratio')
    labelnames = (r'$N_{%s} + N_{%s}$' % (var_str[vartop],var_str[varbot]),r'$%s(N_{%s}/N_{%s})$' % (label_prefix,var_str[vartop],var_str[varbot]))

    if unset_vrange:
        ratio_vrange = (np.min(ratio_baseline_masked),np.max(ratio_baseline_masked))

    vrange = [(0,np.max(sumcounts)),ratio_vrange]

    fig = plt.figure(3,(14,8))
    fig.clf()

    for index, data in enumerate((sumcounts, ratio_baseline_masked)):
        ax = fig.add_subplot(1,2,index+1,aspect=1)
        cmap = cm.jet
        cmap.set_bad('k')
        imextent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1])
        im = ax.imshow(data.T, 
                       vmin = vrange[index][0],
                       vmax = vrange[index][1],
                       extent=imextent,
                       interpolation='nearest',
                       origin='lower'
                       )
        ax.set_title('%s/%s ' % (var_str[vartop],var_str[varbot])+titlenames[index])
        ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
        ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        cb = plt.colorbar(im,orientation='vertical')
        cb.set_label(labelnames[index],fontsize=16)

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
        absmag_lim_loz = appmag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = appmag_lim - cosmology.dmod_flat(0.25)
        size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        if plot_mag_lims:
            ax.axvline(absmag_lim,color='w',linestyle='dashed')
            ax.axvline(absmag_lim_loz,color='g',linestyle='dashed')
            ax.axvline(absmag_lim_hiz,color='g',linestyle='dashed')

        ax.axhline(size_1arcsec, color='w', linestyle='dashed')
        #plt.show()
        #plt.draw()

    fig.savefig(plots_path+'%s_r%s%s_ratio_baseline.png' % (var_def,vartop,varbot), dpi=200)
    p_ratio.close()

    return ax

def plot_ratio_baseline_redshift(task_dict,vartop=0,varbot=1):

    """ Plot the morphology baseline at all redshift slices
        in the binned data. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    vartop: int
        Index for the variable in the numerator of the baseline
        morphology ratio. 

    varbot: int
        Index for the variable in the denominator of the baseline
        morphology ratio. 

    Returns
    -------
    None

    Notes
    -------


    """

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load the GZ2 sample data

    """
    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var1_str,bintype))
        p_var2 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var2_str,bintype))
        d_var1 = p_var1[0].data.astype(int)                         # Data is pre-binned
        d_var2 = p_var2[0].data.astype(int)

        zbins = p_var1['REDSHIFT_BIN_CENTERS'].data['centers']
        magbins = p_var1['MR_BIN_CENTERS'].data['centers']
        sizebins = p_var1['R50_KPC_BIN_CENTERS'].data['centers']
        
        p_var1.close()
        p_var2.close()
    """
    assert bintype is 'rawlikelihood', \
        "task_dict['bintype'] must be 'rawlikelihood' to run this module"

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path_task+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path_task+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        centers_redshift, centers_mag, centers_size, edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    hmask_hypercube = np.ma.masked_equal(d_allvar,0.)

    fig = plt.figure(1)
    fig.clf()

    if len(centers_redshift) < 24:
        maxzind = -1
    else:
        maxzind = 24
    for idx,z in enumerate(centers_redshift[:maxzind]):
        ax = fig.add_subplot(6,4,idx+1)
        cmap = cm.jet
        cmap.set_bad('k')
        if task_dict['ratio_type'] is 'linear':
            label_prefix = ''
            #imarr = np.squeeze(hmask_hypercube[0,idx,:,:]/np.sum(hmask_hypercube[:,idx,:,:],axis=0))
            imarr = np.squeeze(hmask_hypercube[vartop,idx,:,:]/hmask_hypercube[varbot,idx,:,:])
        else:
            label_prefix = 'log_{10}'
            #imarr = np.squeeze(np.log10(hmask_hypercube[0,idx,:,:]/np.sum(hmask_hypercube[:,idx,:,:],axis=0)))
            imarr = np.squeeze(np.log10(hmask_hypercube[vartop,idx,:,:]/hmask_hypercube[varbot,idx,:,:]))

        imarr_mask = np.ma.masked_equal(imarr,0)
        im = ax.imshow(imarr_mask.T, 
                       extent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1]),
                       interpolation="nearest", origin='lower',
                       vmin = -2., vmax = 2., 
                       cmap=cmap)
        ax.set_aspect('auto')
        #cb = plt.colorbar(im)
        rc(('xtick','ytick'), labelsize=10)

        ax.set_xlim(min(edges_mag),max(edges_mag))
        ax.set_ylim(min(edges_size),max(edges_size))

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(z)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(z))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(edges_redshift[idx])
        size_1arcsec = cosmology.ang_scale_flat(edges_redshift[idx])
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        ax.axvline(absmag_lim,color='w',linestyle='dashed')
        ax.axhline(size_1arcsec, color='w', linestyle='dashed')

    cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r'$%s(N_{%s}/N_{%s})$' % 
                 (label_prefix,var_str[vartop],var_str[varbot]), 
                 fontsize=20)
    fig.text(0.5,0.05,r'$M_R [mag]$',ha='center',fontsize=20)
    fig.text(0.05,0.5,r'$R_{50} [kpc]$',va='center', rotation='vertical',fontsize=20)
    fig.text(0.5,0.94,
             '%s ratio per redshift bin for GZ2' % task_dict['var_def'], 
             fontsize=22, ha='center')

    fig.savefig(plots_path+'%s_ratio_redshift.png' % task_dict['var_def'], dpi=200)

    return None

def bootstrap_fmin(f, p0, x, y, z, mask, funcname, nboot):

    """ Run optimization function on data and permute weights
        on fit to estimate the uncertainty. 

    Parameters
    ----------
    f : str
        Name of the minimization function.

    p0 : numpy.array
        Starting parameters to which to fit data.

    x : numpy.array
        x-axis of the data grid

    y : numpy.array
        y-axis of the data grid

    z : numpy.array
        2-D array of the data values to be fit

    mask: numpy.masked_array
        Boolean array giving values of data that are not fit

    funcname: str
        Name of functional form to which data is fit. Options:

        * `sb` - original, 9-parameter equation from GZ1

        * `tilt - function varying linearly in both x and y 

    nboot : int
        Number of tries to permute weights. Higher numbers give
        better accuracy on the estimated uncertainty. 


    Returns
    -------
    p : tuple
        Best fit parameters of the function based on original,
        unweighted bask

    perr : tuple
        Error in the fit parameters based on the mean of the
        16th and 84th percentiles. 

    Notes
    -------


    """

    plist = np.zeros((nboot, len(p0)), np.float)                                # Zero array for number of tries, number of parameters
    ndata = len(mask.ravel().nonzero()[0])                                      # number of non-masked cells in the ratio data
    for i in range(nboot):
        bootmask = np.zeros(mask.shape)                                         # Grid of zeros same size as mask, data
        while bootmask.sum() < ndata:                                           # Loop until grid is filled at data locations
            rx = int(np.random.uniform(0, x.shape))                             # Pick random indices from the bin arrays
            ry = int(np.random.uniform(0, y.shape))                             # 
            if mask[rx,ry] > 0.000001:                                          # If existing data mask in that cell is non-zero
                bootmask[rx,ry] += 1                                            #   then put a 1 in bootmask cell
                                                                                # Preserves total weight, but assigns random weights
                                                                                #   to each data cell. Perturbs to get multiple guesses. 
        bootmask *= mask                                                        # Multiply new weights by original mask
        plist[i] = fmin_powell(f, p0, (x, y, z, bootmask,funcname),disp=0)      # Fit parameters based on new mask passed to function
    p = fmin_powell(f, p0, (x, y, z, mask,funcname))                            # Fit parameters based on original (unweighted) mask
    pmed = np.median(plist, axis=0)                                             # Find the median best fit parameters for each variation in the list
    perr = np.array([(scoreatpercentile(plist[:,k], 84.0) -                     # The error is the mean of the 84th percentile parameter
            scoreatpercentile(plist[:,k], 16.0))/2.0 for k in range(len(p))])   #   and the 16th percentile parameter. Ah. Looping is only to estimate error.
    return p, perr

def ratio_function(p,x,y,funcname):

    """ Return values for the functional form of the 
        morphology ratio.

    Parameters
    ----------
    p0 : numpy.array
        Starting parameters to which to fit data.

    x : numpy.array
        x-axis of the data grid

    y : numpy.array
        y-axis of the data grid

    funcname: str
        Name of functional form to which data is fit. Options:

        * `sb` - original, 9-parameter equation from GZ1

        * `tilt - function varying linearly in both x and y 

    Returns
    -------
    z : numpy.array
        2-D array with values of funcname(x,y)

    Notes
    -------


    """

    assert funcname in ('sb','tilt'), \
        'funcname must be either "tilt" (for linear stretch in mag,size) or "sb" (original fn. from Bamford+09)'

    if funcname is 'sb':

        a, b, c, d, e, f, g, h, i = p
        z = np.zeros((len(x), len(y)), np.float)
        x0 = b**(-(a + h*y**i)) + c
        x1 = d + e*(x0-c)
        for i in range(len(x)):
            z[i] = f / (1.0 + np.exp((x0 - x[i])/x1)) + g

    if funcname is 'tilt':

        z = np.zeros((len(x),len(y)),np.float)
        
        for i in np.arange(len(x)):
            z[i] = p[0]*(x[i] - p[1]) + p[2]*(y - p[3]) + p[4]

    return z

def ratio_minfunc(p, x, y, z, w, funcname):

    """ Compute chi-squared statistic for the fit to weighted
        data using `funcname`. 

    Parameters
    ----------
    p0 : numpy.array
        Starting parameters to which to fit data.

    x : numpy.array
        x-axis of the data grid

    y : numpy.array
        y-axis of the data grid

    z : numpy.array
        2-D array of the data values to be fit

    w: numpy.array
        float array giving weights of each data point in z

    funcname: str
        Name of functional form to which data is fit. Options:

        * `sb` - original, 9-parameter equation from GZ1

        * `tilt - function varying linearly in both x and y 

    Returns
    -------
    s : float
        chi-squared statistic from fit to data

    Notes
    -------


    """

    f = ratio_function(p, x, y, funcname)
    r2 = (z - f)**2                                     # Difference between model and data
    df = (w > 0).ravel().astype(np.int).sum() - len(p)  # Degrees of freedom: N_cells - N_parameters
    s = (w * r2).sum() / df                             # Chi-squared statistic
    return s

def fit_ratio_baseline(task_dict, nboot=50, plot=False, unset_vrange=False, vartop=0, varbot=1):

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']
    funcname = task_dict['funcname']

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    sigma = pyfits.getdata(fits_path_task+'%s_r%s%s_local_ratio_baseline_sigma.fits' % (var_def,vartop,varbot))
    ratio = pyfits.getdata(fits_path_task+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot))
    
    # Steven's and my initial parameters

    if funcname is 'sb':
        pinit_sb = np.array([-0.3, 1.7, -22.8, 0.3, 0.3, -2.6, 1.5, 1.0, 1.0])
        pinit_kw = np.array([-3.5, 1.9, -22.8, 0.3, 0.3, -4.4, 2.2, 1.1, 1.0])
        pinit = pinit_kw

    if funcname is 'tilt':
        pinit = np.array([0.5,-20.0,0.5,7.0,0.1])

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
                             centers_mag, centers_size, ratio, weights, funcname, nboot)
    chi2r = ratio_minfunc(p, centers_mag, centers_size, ratio, weights,funcname)

    print ' '
    print 'Number of bootstrap iterations: %i' % nboot
    print 'param = %s' % p
    print 'param errors = %s' % perr
    print 'param errors (percent) = %s' % (perr/p)
    print 'chi2r = %f' % chi2r
    print ' '

    ratio_baseline_fit = ratio_function(p, centers_mag, centers_size, funcname)
    pickle.dump((p, perr), file(pkl_path+'%s_r%s%s_ratio_baseline_fit.pkl' % (var_def,vartop,varbot), 'w'))
    res = ratio - ratio_baseline_fit
    normres = res/sigma
    np.putmask(normres, weights < 0.000001, unknown_ratio)

    if plot:
        plot_ratio_function(p,centers_mag,centers_size, unset_vrange,funcname)

    pyfits.writeto(fits_path_task+'%s_r%s%s_local_ratio_baseline_fit.fits' % (var_def,vartop,varbot),
                   ratio_baseline_fit, clobber=True)    
    pyfits.writeto(fits_path_task+'%s_r%s%s_local_ratio_baseline_fitres.fits' % (var_def,vartop,varbot),
                   res, clobber=True) 
    pyfits.writeto(fits_path_task+'%s_r%s%s_local_ratio_baseline_fitnormres.fits' % (var_def,vartop,varbot),
                   normres, clobber=True) 

    return None

def plot_ratio_baseline_fit(task_dict, 
                            fitbins=1000, 
                            unset_vrange = False, 
                            vrange = (min_ratio,max_ratio), 
                            match_databins = False, 
                            plot_local_transparent = True,
                            kernel_size = 3, 
                            plot_contour=False,
                            vartop=0,varbot=1
                            ):

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    ratio_baseline_fit = pyfits.getdata(fits_path_task+'%s_r%s%s_local_ratio_baseline_fit.fits' % (var_def,vartop,varbot))

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    ratio_baseline_masked = pickle.load(open(pkl_path+'%s_r%s%s_local_ratio_baseline_masked.pkl' % (var_def,vartop,varbot),'rb'))
    ratio_baseline_fit_masked = ma.array(ratio_baseline_fit,mask=ratio_baseline_masked.mask)

    kernel = np.ones((kernel_size,kernel_size),dtype=int)
    convarr = np.logical_not(ratio_baseline_masked.mask)
    bw = convolve2d(np.logical_not(ratio_baseline_masked.mask).astype(int),
                    kernel,mode='same').astype(bool).astype(int) * -1

    # Plot the best fit function as an opaque layer

    pfit, pfit_err = pickle.load(file(pkl_path+'%s_r%s%s_ratio_baseline_fit.pkl' % (var_def,vartop,varbot), 'r'))
    if match_databins:
        magbinsplot,sizebinsplot = centers_mag, centers_size
    else:
        magbinsplot,sizebinsplot = np.linspace(edges_mag[0],edges_mag[-1],fitbins),np.linspace(edges_size[0],edges_size[-1],fitbins)

    fitarray_fn = ratio_function(pfit, magbinsplot, sizebinsplot,task_dict['funcname'])
    fit_extent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1])

    if task_dict['ratio_type'] is 'linear':
        label_prefix = ''
    else:
        label_prefix = 'log_{10}'

    fig = plt.figure(5)
    fig.clf()
    ax = fig.add_subplot(111)
    if unset_vrange:
        vrange = (np.min(fitarray_fn),np.max(fitarray_fn))

    imf = ax.imshow(fitarray_fn.T,
                   alpha=1.0,
                   extent = fit_extent,
                   vmin = vrange[0], vmax = vrange[1], 
                   interpolation='nearest', origin='lower')
    cb = plt.colorbar(imf)
    cb.set_label(r'$%s(n_{%s}/n_{%s})$' % (label_prefix,var_str[vartop],var_str[varbot]),fontsize=20)

    # Plot the local masked baseline relation as semi-transparent layer

    if plot_local_transparent:
        masked_extent = fit_extent
        cmap_gray = cm.gray
        cmap_gray.set_bad(alpha=0.0)
        im = ax.imshow(ratio_baseline_masked.T, 
                       alpha=0.5,
                       cmap = cmap_gray,
                       extent = masked_extent,
                       vmin = vrange[0], vmax = vrange[1], 
                       interpolation='nearest',
                       origin='lower')

    # Plot the outline of the local region. 
        # This plots the contour, but not in the right region (or quite the right size)
    
    if plot_contour:
        cont = plt.contour(centers_mag, centers_size, bw, 
                           levels = [-1], 
                           color='w', 
                           linewidths=3, linestyles='dashed')

        zc = cont.collections[0]
        zc.set_color('w')
 
    SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
    absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
    size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
    ax.autoscale(False)
    ax.plot(SBlim_mag, SBlim_size,'w--')
    ax.axhline(size_1arcsec, color='w', linestyle='dashed')
    ax.axvline(absmag_lim,  color='w', linestyle='dashed')

    # Set general axes properties

    ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
    ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
    ax.set_title('Baseline ratio data + fit',fontsize=22)
    ax.set_aspect('auto')
    rc(('xtick','ytick'), labelsize=12)

    #plt.show()
    #plt.draw()

    fig.savefig(plots_path+'%s_r%s%s_ratio_baseline_fit.png' % (var_def,vartop,varbot), dpi=200)

    return None

def plot_ratio_function(p=np.array([-3.5, 1.9, -22.8, 0.3, 0.3, -4.4, 2.2, 1.1, 1.0]),
                        x = np.linspace(-24,-16, 50),
                        y = np.linspace(0,15,50),
                        unset_vrange = False,
                        funcname='sb',
                        fignum=4,
                        vrange = (min_ratio,max_ratio)):

    # The effect of the nine parameters:

    # a - horizontal shift (along mag axis)
    # b - curvature tilt
    # c - vertical shift (along size axis)
    # d - gradient of the curve (smaller values = very sharp transition)
    # e - vertical flare at more negative magnitudes
    # f - multiplicative scaling of the cell values
    # g - additive scaling of the cell values
    # h - curvature strength (multiplicative)
    # i - curvature strength (exponential)

    fit = ratio_function(p,x,y,funcname)

    if unset_vrange:
        vrange = (np.min(fit),np.max(fit))

    fig = plt.figure(fignum)
    fig.clf()
    ax = fig.add_subplot(111)
    im = ax.imshow(fit.T, 
        extent = (min(x),max(x),min(y),max(y)), 
        vmin = vrange[0], vmax = vrange[1],
        interpolation='nearest', 
        origin='lower')
    cb = plt.colorbar(im)
    ax.set_aspect('auto')
    ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
    ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
    ax.set_title('Ratio function',fontsize=22)

    #plt.show()
    #plt.draw()

    fig.savefig(plots_path+'ratio_function.png', dpi=200)

    return None

def get_task_dict(task):

    # Define dictionaries that define which GZ2 parameters to measure

    assert task in ('smooth','edgeon','bar','spiral','bulge',\
        'rounded','odd','odd_feature','bulge_shape','arms_winding','arms_number'), \
        "Task must be a string in 'smooth','edgeon','bar','spiral','bulge','rounded','odd','odd_feature','bulge_shape','arms_winding','arms_number')"

    # Default dictionary

    task_dict = {'task_name_count' : 't01_smooth_or_features_total_weight',
                 'task_names_wf' : ('t01_smooth_or_features_a01_smooth_weighted_fraction',
                                    't01_smooth_or_features_a02_features_or_disk_weighted_fraction',
                                    't01_smooth_or_features_a03_star_or_artifact_weighted_fraction'),
                 'var_str' : ('el','sp','ar'),
                 'var_def' : 'task01',
                 'task_str' : task,
                 'wp_lim' : 0.5,
                 'wp_type' : 'binary',
                 'ratio_type' : 'log',
                 'bintype': 'rawlikelihood', 
                 'min_prob': 0.8,
                 'min_classifications': 30, 
                 'min_galperbin': 25, 
                 'direct': True, 
                 'reverse' : False,
                 'funcname' : 'tilt'
                 }

    # Change entries for various GZ2 tasks

    if task is 'smooth':
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['direct'][0,1] = False
        task_dict['funcname'] = 'sb'
        pass
 
    if task is 'edgeon':
        task_dict['task_name_count'] = 't02_edgeon_total_weight'
        task_dict['task_names_wf'] = ('t02_edgeon_a04_yes_weighted_fraction',
                                      't02_edgeon_a05_no_weighted_fraction')
        task_dict['var_str'] = ('edgeon','notedgeon')
        task_dict['var_def'] = 'task02'
        task_dict['min_classifications'] = 15
        task_dict['min_galperbin'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
 
    if task is 'bar':
        task_dict['task_name_count'] = 't03_bar_total_weight'
        task_dict['task_names_wf'] = ('t03_bar_a06_bar_weighted_fraction',
                                      't03_bar_a07_no_bar_weighted_fraction')
        task_dict['var_str'] = ('bar','nobar')
        task_dict['var_def'] = 'task03'
        task_dict['min_classifications'] = 10
        task_dict['min_galperbin'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
 
    if task is 'spiral':
        task_dict['task_name_count'] = 't04_spiral_total_weight'
        task_dict['task_names_wf'] = ('t04_spiral_a08_spiral_weighted_fraction',
                                      't04_spiral_a09_no_spiral_weighted_fraction')
        task_dict['var_str'] = ('spiral','nospiral')
        task_dict['var_def'] = 'task04'
        task_dict['min_classifications'] = 10
        task_dict['min_galperbin'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
 
    if task is 'odd':
        task_dict['task_name_count'] = 't06_odd_total_weight'
        task_dict['task_names_wf'] = ('t06_odd_a14_yes_weighted_fraction',
                                      't06_odd_a15_no_weighted_fraction')
        task_dict['var_str'] = ('odd','notodd')
        task_dict['var_def'] = 'task06'
        task_dict['ratio_type'] = 'log'
        task_dict['min_classifications'] = 30
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
                     
    if task is 'bulge':
        task_dict['task_name_count'] = 't05_bulge_prominence_total_weight'
        task_dict['task_names_wf'] = ('t05_bulge_prominence_a10_no_bulge_weighted_fraction',
                                      't05_bulge_prominence_a11_just_noticeable_weighted_fraction',
                                      't05_bulge_prominence_a12_obvious_weighted_fraction',
                                      't05_bulge_prominence_a13_dominant_weighted_fraction')
        task_dict['var_str'] = ('nobulge','justnoticeable','obvious','dominant')
        task_dict['wp_lim'] = 0.5
        task_dict['min_classifications'] = 10
        task_dict['wp_type'] = 'weighted'
        task_dict['ratio_type'] = 'log'
        task_dict['var_def'] = 'task05'
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['direct'][0,3]   = True
        task_dict['direct'][3,0]   = True
        task_dict['direct'][2,3]   = True
        task_dict['direct'][3,2]   = True

    if task is 'rounded':
        task_dict['task_name_count'] = 't07_rounded_total_weight'
        task_dict['task_names_wf'] = ('t07_rounded_a16_completely_round_weighted_fraction',
                                      't07_rounded_a17_in_between_weighted_fraction',
                                      't07_rounded_a18_cigar_shaped_weighted_fraction')
        task_dict['var_str'] = ('round','inbetween','cigar')
        task_dict['wp_lim'] = 0.3
        task_dict['wp_type'] = 'weighted'
        task_dict['ratio_type'] = 'log'
        task_dict['var_def'] = 'task07'
        task_dict['min_classifications'] = 15
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a01_smooth_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['funcname']='sb'
        task_dict['direct'][0,1] = False
        task_dict['direct'][0,2] = False
        task_dict['direct'][1,2] = False

    if task is 'arms_winding':
        task_dict['task_name_count'] = 't10_arms_winding_total_weight'
        task_dict['task_names_wf'] = ('t10_arms_winding_a28_tight_weighted_fraction',
                                      't10_arms_winding_a29_medium_weighted_fraction',
                                      't10_arms_winding_a30_loose_weighted_fraction')
        task_dict['var_str'] = ('tight','medium','loose')
        task_dict['var_def'] = 'task10'
        task_dict['wp_lim'] = 0.4
        task_dict['wp_type'] = 'weighted'
        task_dict['min_classifications'] = 5
        task_dict['ratio_type'] = 'log'
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction','t04_spiral_a08_spiral_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)

    if task is 'arms_number':
        task_dict['task_name_count'] = 't11_arms_number_total_weight'
        task_dict['task_names_wf'] = ('t11_arms_number_a31_1_weighted_fraction',
                                      't11_arms_number_a32_2_weighted_fraction',
                                      't11_arms_number_a33_3_weighted_fraction',
                                      't11_arms_number_a34_4_weighted_fraction',
                                      't11_arms_number_a36_more_than_4_weighted_fraction',
                                      't11_arms_number_a37_cant_tell_weighted_fraction')
        task_dict['var_str'] = ('1','2','3','4','5+','cant_tell')
        task_dict['var_def'] = 'task11'
        task_dict['wp_lim'] = 0.25
        task_dict['wp_type'] = 'weighted'
        task_dict['ratio_type'] = 'log'
        task_dict['min_classifications'] = 5
        task_dict['min_galperbin'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction','t04_spiral_a08_spiral_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['direct'][0,1] = False
        task_dict['direct'][0,5] = False
        task_dict['direct'][1,0] = False
        task_dict['direct'][1,2] = False
        task_dict['direct'][1,3] = False
        task_dict['direct'][1,4] = False
        task_dict['direct'][1,5] = False
        task_dict['direct'][2,1] = False
        task_dict['direct'][2,5] = False
        task_dict['direct'][3,1] = False
        task_dict['direct'][3,5] = False
        task_dict['direct'][4,1] = False
        task_dict['direct'][4,5] = False
        task_dict['direct'][5,0] = False
        task_dict['direct'][5,1] = False
        task_dict['direct'][5,2] = False
        task_dict['direct'][5,3] = False
        task_dict['direct'][5,4] = False

    if task is 'bulge_shape':
        task_dict['task_name_count'] = 't09_bulge_shape_total_weight'
        task_dict['task_names_wf'] = ('t09_bulge_shape_a25_rounded_weighted_fraction',
                                      't09_bulge_shape_a26_boxy_weighted_fraction',
                                      't09_bulge_shape_a27_no_bulge_weighted_fraction')
        task_dict['var_str'] = ('rounded','boxy','nobulge')
        task_dict['var_def'] = 'task09'
        task_dict['wp_type'] = 'trinary'
        task_dict['min_prob'] = 0.5
        task_dict['min_classifications'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a04_yes_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['direct'][0,2] = False
        task_dict['direct'][2,0] = False

    if task is 'odd_feature':
        task_dict['task_name_count'] = 't08_odd_feature_total_weight'
        task_dict['task_names_wf'] = ('t08_odd_feature_a19_ring_weighted_fraction',
                                      't08_odd_feature_a20_lens_or_arc_weighted_fraction',
                                      't08_odd_feature_a21_disturbed_weighted_fraction',
                                      't08_odd_feature_a22_irregular_weighted_fraction',
                                      't08_odd_feature_a23_other_weighted_fraction',
                                      't08_odd_feature_a24_merger_weighted_fraction',
                                      't08_odd_feature_a38_dust_lane_weighted_fraction')
        task_dict['var_str'] = ('ring','lens','disturbed','irregular','other','merger','dustlane')
        task_dict['var_def'] = 'task08'
        task_dict['min_prob'] = 0.5
        task_dict['min_classifications'] = 5
        task_dict['min_galperbin'] = 10
        task_dict['dependent_tasks'] = ('t06_odd_a14_yes_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)

    return task_dict

def run_task(task_dict,
             nboot = 1, 
             fitbins = 1000,
             plot = True, 
             unset_vrange = True
             ):
    
    tstart = time.time()

    warnings.filterwarnings('ignore', category=UserWarning,    append=True)
    warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

    assert type(task_dict) is dict, \
        "First argument in run_task must be a dictionary -- use galaxyzoo2.get_task_dict()"

    # Execute all the tasks to bin data, find baseline morphology, fit the baseline, 
    # and adjust the raw vote fractions

    bin_data_idl(task_dict)

    ntask = len(task_dict['task_names_wf'])
    for idx1, vartop in enumerate(np.arange(ntask)):
        setdiff = np.concatenate((np.arange(vartop),np.arange(ntask - (vartop+1)) + (vartop+1) ))
        for idx2, varbot in enumerate(setdiff):

            determine_ratio_baseline(task_dict,vartop,varbot)
            determine_ratio_baseline_sigma(task_dict,plot,vartop,varbot)
            fit_ratio_baseline(task_dict,nboot,plot,unset_vrange,vartop,varbot)
            determine_baseline_correction(task_dict,vartop,varbot)

            if plot:
                plot_ratio_baseline(task_dict,unset_vrange,vartop=vartop,varbot=varbot)
                plot_ratio_baseline_fit(task_dict,fitbins,unset_vrange,vartop=vartop,varbot=varbot)
                plot_ratio_baseline_redshift(task_dict,vartop,varbot)
                plot_baseline_correction(task_dict,vartop=vartop,varbot=varbot)


    adjust_probabilities(task_dict)
    if plot:
        plot_galaxy_counts(task_dict)
        plot_type_fractions(task_dict)

    warnings.resetwarnings()

    tend = time.time()
    print 'Time elapsed for run_task on %s: %i seconds' % (task_dict['task_str'], tend - tstart)

    return None

def run_all_tasks():

    """
    Old method (skipping the Pythonic binning):
        Took 88 minutes to run at nboot = 50, although I think ipython may have stopped when computer went to sleep. 
        Takes 164 seconds to run at nboot=1. 

    New method (including the IDL-called binning):
        Takes 164 seconds to run at nboot=1. 
    """

    tstart = time.time()

    tasklist = ['smooth','edgeon','bar','spiral','odd',
                'bulge','rounded','arms_winding','arms_number',
                'bulge_shape','odd_feature']

    for task in tasklist:
        td = get_task_dict(task)
        run_task(td)

    tend = time.time()
    print 'Time elapsed for run_all_tasks: %i seconds' % (tend - tstart)

    return None

def determine_baseline_correction(task_dict, vartop = 0, varbot = 1):

    # Load best fit data

    var_def  = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var_str[vartop],bintype))
        p_var2 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var_str[varbot],bintype))
        d_var1 = p_var1[0].data.astype(int)                         # Data is pre-binned
        d_var2 = p_var2[0].data.astype(int)

        zbins = p_var1['REDSHIFT_BIN_CENTERS'].data['centers']
        magbins = p_var1['MR_BIN_CENTERS'].data['centers']
        sizebins = p_var1['R50_KPC_BIN_CENTERS'].data['centers']
        
        p_var1.close()
        p_var2.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path_task+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path_task+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        d_var1 = np.squeeze(d_allvar[vartop,:,:,:])
        d_var2 = np.squeeze(d_allvar[varbot,:,:,:])

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    ratio_baseline_masked = pickle.load(file(pkl_path+'%s_r%s%s_local_ratio_baseline_masked.pkl' % (var_def,vartop,varbot),'rb')) 

    var1_cube = d_var1
    var2_cube = d_var2

    total_cube = var1_cube + var2_cube

    # Take the ratio (log or linear) of the two cubes

    if task_dict['ratio_type'] is 'linear':
        ratio = var1_cube.astype(float)/var2_cube
    else:
        ratio = np.log10(var1_cube.astype(float)/var2_cube)

    # Mask for insufficient counts

    count_mask = total_cube < task_dict['min_galperbin']
    finite_mask = np.isinf(ratio)

    mask = count_mask | finite_mask

    # Compute correction

    if task_dict['direct'][vartop,varbot]:
        correction = ratio - ratio_baseline_masked
    else:
        pfit_params, pfit_params_err = pickle.load(file(pkl_path+'%s_r%s%s_ratio_baseline_fit.pkl' % (var_def,vartop,varbot), 'r'))
        pfit = ratio_function(pfit_params, centers_mag, centers_size,task_dict['funcname'])
        correction = ratio - pfit

    np.putmask(correction, mask, unknown_ratio)
    correction_masked = ma.masked_less_equal(correction, unknown_ratio)
    ratio_masked = ma.array(ratio, mask=correction_masked.mask)

    return correction, correction_masked, ratio_masked

def plot_baseline_correction(task_dict,vrange=(min_ratio,max_ratio),vartop=0,varbot=1):

    c,cmasked, ratio_masked = determine_baseline_correction(task_dict, vartop, varbot)

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var_str[vartop],bintype))
        p_var2 = pyfits.open(fits_path_task+'%s_%s_binned_%s.fits' % (var_def,var_str[varbot],bintype))
        d_var1 = p_var1[0].data.astype(int)                         # Data is pre-binned
        d_var2 = p_var2[0].data.astype(int)

        zbins = p_var1['REDSHIFT_BIN_CENTERS'].data['centers']
        edges_mag = p_var1['MR_BIN_EDGES'].data['edges']
        edges_size = p_var1['R50_KPC_BIN_EDGES'].data['edges']

        p_var1.close()
        p_var2.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path_task+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path_task+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        d_var1 = np.squeeze(d_allvar[vartop,:,:,:])
        d_var2 = np.squeeze(d_allvar[varbot,:,:,:])

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    zstep = centers_redshift[1] - centers_redshift[0]

    fig = plt.figure(7)
    fig.clf()
    cmap = cm.jet
    cmap.set_bad('k')

    for idx,z in enumerate(centers_redshift[:24]):
        
        ax = fig.add_subplot(6,4,idx+1)
        im = ax.imshow(cmasked[idx,:,:].T,
                       extent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1]),
                       vmin=vrange[0],vmax=vrange[1],
                       interpolation='nearest',origin='lower')
        ax.set_aspect('auto')

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(z)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(z))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(z)
        size_1arcsec = cosmology.ang_scale_flat(z)
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        ax.axhline(size_1arcsec, color='w', linestyle='dashed')
        ax.axvline(absmag_lim,  color='w', linestyle='dashed')
        rc(('xtick','ytick'), labelsize=6)

    cax = fig.add_axes([0.93, 0.1, 0.03, 0.8])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('baseline correction')
    fig.text(0.5,0.95,'Baseline correction per redshift bin for %s' % var_def, fontsize=22, ha='center')

    fig.savefig(plots_path+'%s_r%s%s_baseline_correction.png' % (var_def,vartop,varbot), dpi=200)

    return None

def plot_baseline_correction_slice(task_dict, bslice=5,
                                   vmin = -2.0, vmax = 2.0,
                                   vartop = 0, varbot=1, 
                                   smoothfunc = False,
                                   savefig=False):

    c,cmasked, ratio_masked = determine_baseline_correction(task_dict)

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    zstep = edges_redshift[1] - edges_redshift[0]

    fig = plt.figure(9)
    fig.clf()
    cmap = cm.jet
    cmap.set_bad('k')

    # Data

    zslice = ratio_masked[bslice,:,:]

    ax1 = fig.add_subplot(131)
    cmap.set_bad('k')
    im1 = ax1.imshow(zslice.T,
                     extent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1]),
                     vmin = vmin, vmax = vmax,
                     interpolation='nearest',origin='lower')
    ax1.set_title('Data ratio, %5.3f < z < %5.3f' % (centers_redshift[bslice]-zstep/2.,centers_redshift[bslice]+zstep/2.))
    ax1.set_aspect('auto')
    cb = fig.colorbar(im1)
    cb.set_label('ratio')

    SBlim_mag = (SBlim_app - cosmology.dmod_flat(0.07)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(0.07))**2))
    absmag_lim = appmag_lim - cosmology.dmod_flat(0.07)
    size_1arcsec = cosmology.ang_scale_flat(0.07)
    ax1.autoscale(False)
    ax1.plot(SBlim_mag, SBlim_size,'w--')
    ax1.axhline(size_1arcsec, color='w', linestyle='dashed')
    ax1.axvline(absmag_lim,  color='w', linestyle='dashed')
    ax1.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
    ax1.set_xlabel(r'$M_R [mag]$',fontsize=22)
    plt.show()
    plt.draw()

    # Function

    pfit, pfit_err = pickle.load(file(pkl_path+'%s_r%s%s_ratio_baseline_fit.pkl' % (var_def,vartop,varbot), 'r'))
    if smoothfunc:
        fit = ratio_function(pfit, np.linspace(edges_mag[0],edges_mag[-1],1000),np.linspace(edges_size[0],edges_size[-1],1000),task_dict['funcname'])
    else:
        fit = ratio_function(pfit, centers_mag, centers_size,task_dict['funcname'])
    fit_extent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1])

    ax2 = fig.add_subplot(132)
    cmap.set_bad(alpha=1.0)
    im2 = ax2.imshow(fit.T, 
                     alpha=1.0,
                     cmap = cmap,
                     extent = fit_extent,
                     vmin = vmin, vmax = vmax,
                     interpolation='nearest', 
                     origin='lower')
    cb = plt.colorbar(im2)
    ax2.set_title('Baseline fit')
    ax2.set_xlabel(r'$M_R [mag]$',fontsize=22)

    zslice_opaque = ma.copy(zslice)
    zslice_opaque.mask = np.logical_not(zslice_opaque.mask)
    zslice_opaque[np.logical_not(zslice_opaque.mask)] = 255.
    cmap_gray = cm.gray
    cmap_gray.set_bad(alpha=0.0)
    im2a = ax2.imshow(zslice_opaque.T, 
                      alpha=0.5,
                      extent = fit_extent,
                      vmin = vmin, vmax = vmax,
                      cmap = cmap_gray,
                      interpolation='nearest',
                      origin='lower')

    ax2.set_aspect('auto')
    ax2.autoscale(False)

    # Correction

    ax3 = fig.add_subplot(133)
    cmap.set_bad('k',alpha=1.0)
    im3 = ax3.imshow(cmasked[bslice,:,:].T,
                     extent=fit_extent,
                     vmin = vmin, vmax = vmax,
                     interpolation='nearest',origin='lower')
    ax3.set_title('Correction')
    cb = fig.colorbar(im3)
    cb.set_label('baseline correction')

    SBlim_mag = (SBlim_app - cosmology.dmod_flat(0.07)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(0.07))**2))
    absmag_lim = appmag_lim - cosmology.dmod_flat(0.07)
    size_1arcsec = cosmology.ang_scale_flat(0.07)
    ax3.plot(SBlim_mag, SBlim_size,'w--')
    ax3.axhline(size_1arcsec, color='w', linestyle='dashed')
    ax3.axvline(absmag_lim,  color='w', linestyle='dashed')
    ax3.set_xlabel(r'$M_R [mag]$',fontsize=22)
    ax3.set_xlim(fit_extent[:2])
    ax3.set_ylim(fit_extent[2:])
    ax3.autoscale(False)
    ax3.set_aspect('auto')

    if savefig:
        fig.savefig(plots_path+'%s_r%s%s_baseline_correction_slice%03i.png' % (task_dict['var_def'],vartop,varbot,bslice), dpi=200)

    return None

def ratio_adj(p_el,p_sp,correction,ratio_type):
    
    assert ratio_type in ('linear','log'), \
        "Ratio type must be defined and either 'linear' or 'log'"

    if ratio_type is 'log':
        ratio_adj = (p_el / p_sp) / (10.**(correction))

    if ratio_type is 'linear':
        ratio_adj = (p_el / p_sp) / correction

    return ratio_adj

def p_x(p_el,p_sp):

    p_x = 1. - p_el - p_sp

    return p_x

def p_el_adj(p_el,p_sp,correction,ratio_type):

    denom_el = 1./ratio_adj(p_el,p_sp,correction,ratio_type) + (p_x(p_el,p_sp) / p_el) + 1.
    p_el_adj = 1./denom_el

    return p_el_adj

def p_sp_adj(p_el,p_sp,correction,ratio_type):

    denom_sp = ratio_adj(p_el,p_sp,correction,ratio_type) + (p_x(p_el,p_sp) / p_sp) + 1.
    p_sp_adj = 1./denom_sp

    return p_sp_adj

def p_adj_new(f, corr):

    if corr.ndim > 0:
        assert len(f) - 1 == len(corr), \
            'There must be exactly one more element in the vote fraction \
            (%i) than in the correction (%i) array' % (len(f),len(corr))
    else:
        assert len(f) == 2, \
            'There must be only two vote fractions supplied for a binary correction'

    K = 10.**corr
    termlist = K * (f[1:].astype(float)/ f[0].astype(float))

    denom = np.sum(termlist) + 1.

    p_adj = 1. / denom

    return p_adj

def adjust_probabilities(task_dict, stripe82=False, photoz=False):

    timestart2 = time.time()

    # Load in the raw probabilities from the GZ2 catalog

    file_str = ''
    if stripe82:
        p = pyfits.open(fits_path_task+'%s_data_for_idl_stripe82.fits' % task_dict['var_def'])
        file_str = 'stripe82_'
    if photoz:
        p = pyfits.open(gz2_photoz_table_sample.fits)
        file_str = 'photoz_'
    else:
        p = pyfits.open(fits_path_task+'%s_data_for_idl.fits' % task_dict['var_def'])
        file_str = ''

    gzdata = p[1].data
    p.close()

    # Load in the bin sizes

    if task_dict['bintype'] is 'rawlikelihood':
        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    # Take each galaxy, find which bin it corresponds to, adjust probability

    if photoz:
        redshift = gzdata['photoz']
        mr = (gzdata['PETROMAG_R'] - gzdata['EXTINCTION_R']) - cosmology.dmod_flat(redshift)
        r50_kpc = gzdata['PETROR50_R'] * cosmology.ang_scale_flat(redshift)
    else:
        redshift = gzdata['REDSHIFT']
        mr = gzdata['PETROMAG_MR']
        r50_kpc = gzdata['PETROR50_R_KPC']

    ntask = len(task_dict['task_names_wf'])
    ngals = len(redshift)

    zstep = edges_redshift[1] - edges_redshift[0]
    magstep = edges_mag[1] - edges_mag[0]
    sizestep = edges_size[1] - edges_size[0]

    corr_test, corrmasked_test, ratiomasked_test = determine_baseline_correction(task_dict, 0, 1)
    cshape = corr_test.shape

    allcorr        = np.zeros((ntask,ntask-1,cshape[0],cshape[1],cshape[2]),dtype=float)
    allcorrmasked  = np.zeros_like(allcorr)
    allratiomasked = np.zeros_like(allcorr)

    # Loop over combinations of tasks

    for idx1, var1 in enumerate(np.arange(ntask)):
        setdiff = np.concatenate((np.arange(var1),np.arange(ntask - (var1+1)) + (var1+1) ))
        for idx2, var2 in enumerate(setdiff):
            corr, corrmasked, ratiomasked = determine_baseline_correction(task_dict, var1, var2)
            allcorr[idx1,idx2,:,:,:] = corr
            allcorrmasked[idx1,idx2,:,:,:] = corrmasked
            allratiomasked[idx1,idx2,:,:,:] = ratiomasked

    p_raw = np.zeros((ntask,ngals),dtype=float)
    for idx,task in enumerate(task_dict['task_names_wf']):
        p_raw[idx,:] = gzdata[task]

    p_task_counts = gzdata[task_dict['task_name_count']]

    # Start adjusting

    p_adj = np.zeros((ntask,ngals),dtype=float)
    corrarr = np.zeros((ntask-1,ntask-1,ngals),dtype=float) + unknown_ratio
    corrbinarr = np.zeros((3,ntask-1,ntask-1,ngals),dtype=int) + unknown_ratio
    corrgoodarr = np.zeros(ngals,dtype=bool)

    outsidebinscount = 0
    corrcount = 0
    zerocorrcount = 0
    unknowncorrcount = 0
    vfzerocount = 0
    vfnonzerocount = 0

    # Loop over all galaxies in the GZ catalog with redshifts

    for idx, (z_temp,m_temp,s_temp) in enumerate(zip(redshift,mr,r50_kpc)):

        f = p_raw[:,idx]
        p_temp = np.zeros_like(f)
        zbin = np.abs(z_temp >= edges_redshift).argmin() - 1
        mbin = np.abs(m_temp >= edges_mag).argmin() - 1
        sbin = np.abs(s_temp >= edges_size).argmin() - 1

        if -1 not in (zbin,mbin,sbin):
            applycorr = np.squeeze(allcorrmasked[:,:,zbin,mbin,sbin])

            return applycorr

            # Check to see if any corrections exist for this bin
            if np.sum((np.isnan(applycorr)) & (applycorr == unknown_ratio)) == len(applycorr.ravel()):
                unknowncorrcount += 1
                p_temp = f
            else:

                # Bins with masked corrections or NaN set to none

                applycorr[applycorr == unknown_ratio] = 0.
                applycorr[np.isnan(applycorr)] = 0.

                # See if there are any bins within 1.2 times the distance to the nearest bin center?

                # Check if any corrections are left to be applied

                if np.sum(applycorr) == 0.:
                    zerocorrcount +=1
                    p_temp = f
                    '''
                    corrarr[:,idx] = unknown_ratio
                    corrbinarr[:,:,idx] = unknown_ratio
                    '''
                else:

                    # Apply corrections to vote fractions

                    corrcount += 1
                    corrgoodarr[idx] = True
                    for fidx, votefrac in enumerate(f):
                        votefracs = np.array(list(chain.from_iterable(([votefrac], f[:fidx], f[fidx+1:]))))
                        ac = applycorr[fidx,:] if applycorr.ndim > 1 else np.array(applycorr[fidx])
                        if votefracs[0] == 0.:
                            vfzerocount += 1
                            p_temp[fidx] = votefracs[0]
                        else:
                            vfnonzerocount += 1
                            p_temp[fidx] = p_adj_new(votefracs,ac)

                        #corrarr[fidx,:,idx] = applycorr[fidx,:]
                        #corrbinarr[:,fidx,:,idx] = np.resize(np.array((zbin,mbin,sbin)),(3,3))

                p_adj[:,idx] = p_temp
          
        else: # if galaxy does not appear in the correction area. Shouldn't happen.
            outsidebinscount += 1
            p_adj[:,idx] = f
            '''
            corrarr[:,idx] = unknown_ratio
            corrbinarr[:,:,idx] = unknown_ratio
            '''

    p_both = np.array((p_raw,p_adj))

    pickle.dump(p_both, open(pkl_path+'%s_%sp_adj.pkl' % (task_dict['var_def'], file_str),'wb')) 
    pickle.dump(corrarr, open(pkl_path+'%s_%scorrarr.pkl' % (task_dict['var_def'], file_str),'wb')) 
    pickle.dump(corrbinarr, open(pkl_path+'%s_%scorrbinarr.pkl' % (task_dict['var_def'], file_str),'wb')) 
    pickle.dump(corrgoodarr, open(pkl_path+'%s_%scorrgoodarr.pkl' % (task_dict['var_def'], file_str),'wb')) 

    timeend2 = time.time()
    print ' '
    print '%7i galaxies had correction bins outside the binned volume' % outsidebinscount
    print '%7i galaxies had NaN corrections for all pairs of variables' % unknowncorrcount
    print '%7i galaxies had no non-zero corrections in any bin' % zerocorrcount
    print '%7i galaxies had finite corrections available' % corrcount
    print 'Of those: '
    print '   %4.1f percent of tasks were not adjusted because the raw vote fraction was 100 percent' % \
        ((float(vfzerocount)/(corrcount * ntask)) * 100.)
    print '   %4.1f percent of tasks had a non-zero correction applied' % \
        ((float(vfnonzerocount)/(corrcount * ntask)) * 100.)
    print ' '
    print '%7i total galaxies' % p_raw.shape[1]
    print ' '
    print 'Time elapsed to adjust probabilities: %i seconds' % (timeend2 - timestart2)
    print ' '

    return p_raw,p_adj

def plot_galaxy_counts(task_dict,min_classifications=20,magstep=0.25,sizestep=0.5, vmax=1000):

    p = pyfits.open(gz2_full_data_file)
    gzall = p[1].data
    p.close()

    gzdata = gzall[(gzall[task_dict['task_name_count']] > min_classifications) & np.isfinite(gzall['REDSHIFT'])]

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)
    zstep = centers_redshift[1] - centers_redshift[0]
    magmin,magmax = edges_mag.min(),edges_mag.max()
    sizemin,sizemax = edges_size.min(),edges_size.max()

    cmap = cm.jet
    cmap.set_bad('k')

    fig = plt.figure(20)
    fig.clf()

    for idx,z in enumerate(centers_redshift[:24]):
        
        zind = (gzdata['REDSHIFT'] >= z) & (gzdata['REDSHIFT'] < (z+zstep))
        mag = gzdata['PETROMAG_MR'][zind]
        size = gzdata['PETROR50_R_KPC'][zind]

        h,xedges,yedges = np.histogram2d(mag,size,bins=(np.arange(magmin,magmax,magstep),np.arange(sizemin,sizemax,sizestep)))
        hmask = np.ma.masked_equal(h,0)

        ax = fig.add_subplot(6,4,idx+1)
        if h.ndim > 1:
            im = ax.imshow(hmask.T,extent=(magmin,magmax,sizemin,sizemax),
                           vmin=0,vmax=vmax,
                           interpolation='nearest',origin='lower')
        #ax.set_title('%s < z < %s' % (z-zstep/2.,z+zstep/2.))
        ax.set_aspect('auto')

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(z)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(z))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(z)
        size_1arcsec = cosmology.ang_scale_flat(z)
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        ax.axhline(size_1arcsec, color='w', linestyle='dashed')
        ax.axvline(absmag_lim,  color='w', linestyle='dashed')

    cax = fig.add_axes([0.93, 0.1, 0.03, 0.8])
    fig.colorbar(im, cax=cax)
    fig.text(0.5,0.95,'Total galaxy counts per redshift bin for %s' % task_dict['var_def'], fontsize=22, ha='center')
    
    #plt.show()
    #plt.draw()

    fig.savefig(plots_path+'%s_galaxy_counts.png' % task_dict['var_def'], dpi=200)

    return fig

def plot_type_fractions(task_dict, zlo = 0.01, zhi=0.085, zwidth=0.02, stripe82=False):

    # Load data 

    if stripe82:
        s82_str = 'stripe82_'
        p = pyfits.open(gz2_stripe82_data_file)
    else:
        s82_str = ''
        p = pyfits.open(fits_path_task+'%s_data_for_idl.fits' % task_dict['var_def'])

    gzdata = p[1].data
    p.close()

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    redshift = gzdata['REDSHIFT']
    mr = gzdata['PETROMAG_MR']
    task_counts = gzdata[task_dict['task_name_count']]

    zplotbins = np.arange(min(edges_redshift[0],zlo),edges_redshift[-1],zwidth)

    ntask = len(task_dict['var_str'])

    p_both = pickle.load(open(pkl_path+'%s_%sp_adj.pkl' % (task_dict['var_def'],s82_str),'rb')) 
    corrgoodarr = pickle.load(open(pkl_path+'%s_%scorrgoodarr.pkl' % (task_dict['var_def'],s82_str),'rb')) 
    p_raw = np.squeeze(p_both[0,:,:])
    p_adj = np.squeeze(p_both[1,:,:])

    assert len(gzdata) == p_adj.shape[1], \
        'Size of metadata (%i) does not match the adjusted probabilities (%i)' % (len(gzdata),p_adj.shape[1])

    # Empty arrays for the binned type fractions

    p_raw_typefrac            = np.zeros((ntask,len(zplotbins)))
    p_adj_typefrac            = np.zeros_like(p_raw_typefrac)
    p_raw_typefrac_maglim     = np.zeros_like(p_raw_typefrac)
    p_adj_typefrac_maglim     = np.zeros_like(p_raw_typefrac)
    p_raw_typefrac_err        = np.zeros_like(p_raw_typefrac)
    p_adj_typefrac_err        = np.zeros_like(p_raw_typefrac)
    p_raw_typefrac_err_maglim = np.zeros_like(p_raw_typefrac)
    p_adj_typefrac_err_maglim = np.zeros_like(p_raw_typefrac)

    """
    Take the mean of the raw likelihoods for each galaxy per redshift bin
    """

    maglimval = appmag_lim - cosmology.dmod_flat(zhi)

    # Create the two samples
    
    n_all = np.sum((task_counts > task_dict['min_classifications']) & corrgoodarr)
    n_maglim = np.sum((mr < maglimval) & (task_counts > task_dict['min_classifications']) & corrgoodarr)

    # Loop over redshift bin and find the type fractions at each slice

    for idx,zbin in enumerate(zplotbins):

        gals_in_bin = (redshift >= zbin) & \
                      (redshift < (zbin+zwidth)) & \
                      corrgoodarr & \
                      (task_counts > task_dict['min_classifications'])

        p_raw_typefrac[:,idx] = np.mean(p_raw[:,gals_in_bin],axis=1)
        p_adj_typefrac[:,idx] = np.mean(p_adj[:,gals_in_bin],axis=1)
        p_raw_typefrac_err[:,idx] = np.std(p_raw[:,gals_in_bin],axis=1)
        p_adj_typefrac_err[:,idx] = np.std(p_adj[:,gals_in_bin],axis=1)

        gals_in_bin_maglim = (redshift >= zbin) & \
                      (redshift < (zbin+zwidth)) & \
                      (mr < maglimval) & \
                      corrgoodarr & \
                      (task_counts > task_dict['min_classifications'])

        p_raw_typefrac_maglim[:,idx] = np.mean(p_raw[:,gals_in_bin_maglim],axis=1)
        p_adj_typefrac_maglim[:,idx] = np.mean(p_adj[:,gals_in_bin_maglim],axis=1)
        p_raw_typefrac_err_maglim[:,idx] = np.std(p_raw[:,gals_in_bin_maglim],axis=1)
        p_adj_typefrac_err_maglim[:,idx] = np.std(p_adj[:,gals_in_bin_maglim],axis=1)

    pickle.dump(p_raw_typefrac,        open(pkl_path+'%s_%sraw_typefrac.pkl' % (task_dict['var_def'],s82_str),'wb')) 
    pickle.dump(p_adj_typefrac,        open(pkl_path+'%s_%sadj_typefrac.pkl' % (task_dict['var_def'],s82_str),'wb')) 
    pickle.dump(p_raw_typefrac_maglim, open(pkl_path+'%s_%sraw_typefrac_maglim.pkl' % (task_dict['var_def'],s82_str),'wb')) 
    pickle.dump(p_adj_typefrac_maglim, open(pkl_path+'%s_%sadj_typefrac_maglim.pkl' % (task_dict['var_def'],s82_str),'wb')) 

    # Plot results

    fig = plt.figure(12, (12,7))
    fig.clf()
    legend_raw_str = []
    legend_adj_str = []

    ax1 = fig.add_subplot(121)

    font = FontProperties()
    font.set_size(12)

    colorarr = ['r','b','m','g','c','y','k'][::-1]
    for idx,task_str in enumerate(task_dict['var_str']):
        linecolor = colorarr.pop()
        ax1.plot(zplotbins, p_raw_typefrac[idx,:], color=linecolor, linestyle='-' ,linewidth=2)
        legend_raw_str.append('raw %s' % task_str)

    colorarr = ['r','b','m','g','c','y','k'][::-1]
    for idx,task_str in enumerate(task_dict['var_str']):
        linecolor = colorarr.pop()
        ax1.plot(zplotbins, p_adj_typefrac[idx,:], color=linecolor, linestyle='--' ,linewidth=4)
        legend_adj_str.append('adj %s' % task_str)

        #ax1.errorbar(zplotbins, p_raw_typefrac[idx,:], yerr = p_raw_typefrac_err[idx,:], 
        #   color=linecolor, linestyle='-' ,linewidth=2,elinewidth=1)
        #ax1.errorbar(zplotbins, p_adj_typefrac[idx,:], yerr = p_adj_typefrac_err[idx,:], 
        #   color=linecolor, linestyle='-' ,linewidth=4,elinewidth=1)


    plt.legend(legend_raw_str, 'upper left', shadow=True, fancybox=True, prop=font)

    ax1.set_xlim(-0.01,0.25)
    ax1.set_ylim(-0.01,1.2)
    ax1.set_xlabel('redshift')
    ax1.set_ylabel('fraction')
    ax1.set_title('All GZ2 galaxies with redshifts, %s' % task_dict['var_def'])
    ax1.text(0.15,1.0,'%i galaxies' % n_all)

    ax2 = fig.add_subplot(122)

    colorarr = ['r','b','m','g','c','y','k'][::-1]
    for idx,task_str in enumerate(task_dict['var_str']):
        linecolor = colorarr.pop()
        ax2.plot(zplotbins, p_raw_typefrac_maglim[idx,:], color=linecolor, linestyle='-' ,linewidth=2)

    colorarr = ['r','b','m','g','c','y','k'][::-1]
    for idx,task_str in enumerate(task_dict['var_str']):
        linecolor = colorarr.pop()
        ax2.plot(zplotbins, p_adj_typefrac_maglim[idx,:], color=linecolor, linestyle='--' ,linewidth=4)

        #ax2.errorbar(zplotbins, p_raw_typefrac_maglim[idx,:], yerr = p_raw_typefrac_err_maglim[idx,:], 
        #   color=linecolor, linestyle='-' ,linewidth=2,elinewidth=1)
        #ax2.errorbar(zplotbins, p_adj_typefrac_maglim[idx,:], yerr = p_adj_typefrac_err_maglim[idx,:], 
        #   color=linecolor, linestyle='-' ,linewidth=4,elinewidth=1)

    #ax2.axvline(zlo, color='k', linestyle='--')
    ax2.axvline(zhi, color='k', linestyle='--')

    #plt.legend(legend_raw_str + legend_adj_str, 'upper left', shadow=True, fancybox=True, prop=font)
    plt.legend(legend_raw_str, 'upper left', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.2)
    ax2.set_xlabel('redshift')
    ax2.set_ylabel('fraction')
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)
    ax2.text(0.15,1.0,'%i galaxies' % n_maglim)

    fig.savefig(plots_path+'%s_%stype_fractions_new.png' % (task_dict['var_def'],s82_str), dpi=200)

    return None

def plot_all_baselines(paperplot=False):

    fig = plt.figure(13)
    fig.clf()

    titlenames = ('Smooth or features','Edge-on','Bar', 'Spiral structure','Bulge prominence','Odd',
                  'Rounded','Odd feature','Bulge shape','Arms winding','Arms number')

    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')
    odd = get_task_dict('odd')
    odd_feature = get_task_dict('odd_feature')
    arms_winding = get_task_dict('arms_winding')
    arms_number = get_task_dict('arms_number')
    bulge = get_task_dict('bulge')
    bulge_shape = get_task_dict('bulge_shape')
    rounded = get_task_dict('rounded')

    tasklist = (smooth,edgeon,bar,spiral,bulge,odd,rounded,odd_feature,bulge_shape,arms_winding,arms_number)
    nx_plot, ny_plot = 3,4
    left_plot = (0,4,8)
    bottom_plot = (8,9,10,11)

    if paperplot:
        tasklist,titlenames = tasklist[:4],titlenames[:4]
        nx_plot, ny_plot = 2,2
        left_plot = (0,2)
        bottom_plot = (2,3)

    for idx,task_dict in enumerate(tasklist):

        var_def = task_dict['var_def']
        bintype = task_dict['bintype']

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        ratio_baseline_masked = pickle.load(open(pkl_path+'%s_r01_local_ratio_baseline_masked.pkl' % var_def,'rb'))

        if task_dict['ratio_type'] is 'linear':
            label_prefix = ''
        else:
            label_prefix = 'log_{10}'

        ax = fig.add_subplot(nx_plot,ny_plot,idx+1,aspect=1)

        cmap = cm.jet
        cmap.set_bad('k')
        imextent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1])
        im = ax.imshow(ratio_baseline_masked.T, 
                       vmin = -1.5,
                       vmax = 1.5,
                       #vmin = np.min(ratio_baseline_masked),
                       #vmax = np.max(ratio_baseline_masked),
                       extent=imextent,
                       interpolation='nearest',
                       origin='lower'
                       )
        ax.set_title(titlenames[idx])
        if idx in bottom_plot:
            ax.set_xlabel(r'$M_R [mag]$',fontsize=16)
        if idx in left_plot:
            ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        if paperplot or idx==10:
            cb = plt.colorbar(im,orientation='vertical')
            cb.set_label(r'$%s(N_{%s}/N_{%s})$' % (label_prefix,task_dict['var_str'][0],task_dict['var_str'][1]),fontsize=16)

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
        absmag_lim_loz = appmag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = appmag_lim - cosmology.dmod_flat(0.25)
        size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        ax.axhline(size_1arcsec, color='w', linestyle='dashed')


    #fig.tight_layout()
    fig.savefig(dropbox_figs_path+'gz2_baselines.eps', dpi=200)

    return None

def plot_all_type_fractions(zlo = 0.01, zhi=0.085, stripe82=False, paperplot=False):

    fig = plt.figure(14)
    fig.clf()

    titlenames = ('Smooth or features','Edge-on','Bar', 'Spiral structure','Bulge prominence','Odd',
                  'Rounded','Odd feature','Bulge shape','Arms winding','Arms number')

    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')
    odd = get_task_dict('odd')
    odd_feature = get_task_dict('odd_feature')
    arms_winding = get_task_dict('arms_winding')
    arms_number = get_task_dict('arms_number')
    bulge = get_task_dict('bulge')
    bulge_shape = get_task_dict('bulge_shape')
    rounded = get_task_dict('rounded')

    tasklist = (smooth,edgeon,bar,spiral,bulge,odd,rounded,odd_feature,bulge_shape,arms_winding,arms_number)
    nx_plot, ny_plot = 3,4
    left_plot = (0,4,8)
    bottom_plot = (8,9,10,11)

    if paperplot:
        tasklist,titlenames = tasklist[:4],titlenames[:4]
        nx_plot, ny_plot = 2,2
        left_plot = (0,2)
        bottom_plot = (2,3)

    s82_str = 'stripe82_' if stripe82 else ''

    for idx,task_dict in enumerate(tasklist):

        var_def = task_dict['var_def']
        ntask = len(task_dict['var_str'])

        p_raw_typefrac = pickle.load(open(pkl_path+'%s_%sraw_typefrac.pkl' % (task_dict['var_def'],s82_str),'rb')) 
        p_adj_typefrac = pickle.load(open(pkl_path+'%s_%sadj_typefrac.pkl' % (task_dict['var_def'],s82_str),'rb')) 
        p_raw_typefrac_maglim = pickle.load(open(pkl_path+'%s_%sraw_typefrac_maglim.pkl' % (task_dict['var_def'],s82_str),'rb')) 
        p_adj_typefrac_maglim = pickle.load(open(pkl_path+'%s_%sadj_typefrac_maglim.pkl' % (task_dict['var_def'],s82_str),'rb')) 

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        zwidth = 0.02
        zplotbins = np.arange(max(edges_redshift[0],zlo),edges_redshift[-1],zwidth)
        maglimval = appmag_lim - cosmology.dmod_flat(zhi)

        ax1 = fig.add_subplot(nx_plot,ny_plot,idx+1,aspect=1)
        font = FontProperties()
        legend_raw_str = []
        legend_adj_str = []

        colorarr = ['r','b','m','g','c','y','k'][::-1]
        for td_idx,task_str in enumerate(task_dict['var_str']):
            linecolor = colorarr.pop()
            ax1.plot(zplotbins, p_raw_typefrac_maglim[td_idx,:], color=linecolor, linestyle='-' ,linewidth=2)
            legend_raw_str.append('raw %s' % task_str)

        colorarr = ['r','b','m','g','c','y','k'][::-1]
        for td_idx,task_str in enumerate(task_dict['var_str']):
            linecolor = colorarr.pop()
            ax1.plot(zplotbins, p_adj_typefrac_maglim[td_idx,:], color=linecolor, linestyle='--' ,linewidth=4)
            legend_adj_str.append('adj %s' % task_str)

        ax1.axvline(zlo, color='k', linestyle='--')
        ax1.axvline(zhi, color='k', linestyle='--')

#        plt.legend(legend_raw_str, 'upper right', shadow=True, fancybox=True, prop=font)

        ax1.set_xlim(-0.01,0.19)
        ax1.set_ylim(-0.01,1.01)
        if idx in bottom_plot:
            ax1.set_xlabel('redshift')
        if idx in left_plot:
            ax1.set_ylabel('fraction')
        ax1.set_aspect('auto')
        ax1.set_title(titlenames[idx])
        rc('xtick', labelsize=10)
        rc('ytick', labelsize=10)

    fig.tight_layout()
    fig.savefig(dropbox_figs_path+'gz2_%stype_fractions.eps' % s82_str, dpi=200)

    return None

def posterplot_debiasing(zlimval=0.085,stripe82=False):

    fig = plt.figure(15, (22,13))
    fig.clf()

    plt.subplots_adjust(left=0.04,right=0.96,top=0.95,bottom=0.08)

    taskstrings=('smooth','edgeon','bar','spiral')
    titlenames = ('Smooth or features/disk','Edge-on diskgalaxy','Bar', 'Spiral structure')
    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')

    tasklist = (smooth,edgeon,bar,spiral)

    for idx,task_dict in enumerate(tasklist):

        var_def = task_dict['var_def']
        var1_str,var2_str = task_dict['var_str']
        bintype = task_dict['bintype']

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        ratio_baseline_masked = pickle.load(open(pkl_path+'%s_local_ratio_baseline_masked.pkl' % var_def,'rb'))

        if task_dict['ratio_type'] is 'linear':
            label_prefix = ''
        else:
            label_prefix = 'log_{10}'

        ax = fig.add_subplot(2,4,idx+1,aspect=1)

        cmap = cm.jet
        cmap.set_bad('k')
        imextent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1])
        im = ax.imshow(ratio_baseline_masked.T, 
                       vmin = np.min(ratio_baseline_masked),
                       vmax = np.max(ratio_baseline_masked),
                       extent=imextent,
                       interpolation='nearest',
                       origin='lower'
                       )
        ax.set_title(titlenames[idx])
        ax.set_xlabel(r'$M_R [mag]$',fontsize=14)
        if idx is 0:
            ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        cb = plt.colorbar(im,orientation='vertical')
        cb.set_label(r'$%s(N_{%s}/N_{%s})$' % (label_prefix,var1_str,var2_str),fontsize=18)

        plt.xticks(np.arange(edges_mag[0],edges_mag[-1],2))

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
        absmag_lim_loz = appmag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = appmag_lim - cosmology.dmod_flat(0.25)
        size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        ax.axhline(size_1arcsec, color='w', linestyle='dashed')

        var_def = task_dict['var_def']
        var1_str = task_dict['var_str'][0]
        var2_str = task_dict['var_str'][1]

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        zwidth = 0.02
        zplotbins = np.arange(edges_redshift[0],edges_redshift[-1],zwidth)
        maglimval = appmag_lim - cosmology.dmod_flat(zlimval)

        s82_str = ''
        if stripe82:
            s82_str = 'stripe82_'

        p_el_raw_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'rb')) 
        p_sp_raw_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'rb')) 
        p_el_adj_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'rb')) 
        p_sp_adj_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'rb')) 

        ax2 = fig.add_subplot(2,4,idx+1 + 4)
        font = FontProperties()

        ax2.plot(zplotbins, p_el_raw_typefrac_maglim, color='r', linestyle='-' ,linewidth=2)
        ax2.plot(zplotbins, p_sp_raw_typefrac_maglim, color='b', linestyle='--',linewidth=2)
        ax2.plot(zplotbins, p_el_adj_typefrac_maglim, color='r', linestyle='-' ,linewidth=4)
        ax2.plot(zplotbins, p_sp_adj_typefrac_maglim, color='b', linestyle='--',linewidth=4)

        ax2.axvline(0.01, color='k', linestyle='--')
        ax2.axvline(zlimval, color='k', linestyle='--')

        plt.legend((var1_str,var2_str), 'center right', shadow=True, fancybox=True, prop=font)

        ax2.set_xlim(-0.01,0.25)
        ax2.set_ylim(-0.01,1.01)
        ax2.set_xlabel('redshift')
        if idx is 0:
            ax2.set_ylabel('GZ2 vote fraction',fontsize=22)

    fig.savefig(dropbox_figs_path+'gz2_posterplot_debiasing.png', dpi=200)

    """
    # Fill the sixth and empty plot with the legend

    ax2 = fig.add_subplot(2,3,6)
    font = FontProperties()

    ax2.plot(zplotbins, p_el_raw_typefrac_maglim-100, color='r', linestyle='-' ,linewidth=2)
    ax2.plot(zplotbins, p_sp_raw_typefrac_maglim-100, color='b', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_el_adj_typefrac_maglim-100, color='r', linestyle='-' ,linewidth=4)
    ax2.plot(zplotbins, p_sp_adj_typefrac_maglim-100, color='b', linestyle='--',linewidth=4)

    plt.legend(('raw '+r'$r_{high}$','raw '+r'$r_{low}$','debiased '+r'$r_{high}$','debiased '+r'$r_{low}$'), 'center', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.01)
    """

    return None

def poster_table():

    taskstrings=('smooth','edgeon','bar','spiral','odd')
    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')
    odd = get_task_dict('odd')

    tasklist = (smooth,edgeon,bar,spiral,odd)

    p = pyfits.open(gz2_maglim_data_file)
    gzdata = p[1].data
    p.close()
    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    zhi = 0.085
 
    for task_dict in tasklist:

        var_def = task_dict['var_def']
        var1_str,var2_str = task_dict['var_str']
        bintype = task_dict['bintype']
        
        p_el_adj_values = pickle.load(open(pkl_path+'%s_%s_adj.pkl' % (var_def,var1_str),'rb')) 
        p_sp_adj_values = pickle.load(open(pkl_path+'%s_%s_adj.pkl' % (var_def,var2_str),'rb')) 
        corrarr         = pickle.load(open(pkl_path+'%s_corrarr.pkl' % (var_def),'rb')) 
        task_counts = gzdata_withz[task_dict['task_name_count']]
        
        corrgood1 = corrarr > unknown_ratio
        corrgood = corrgood1

        redshift = gzdata_withz['REDSHIFT']
        mr = gzdata_withz['PETROMAG_MR']
        r50_kpc = gzdata_withz['PETROR50_R_KPC']
        
        if var_def is 'task02':
            corrarr2 = pickle.load(open(pkl_path+'%s_corrarr.pkl' % ('task01'),'rb')) 
            corrgood2 = corrarr2 > unknown_ratio
            corrgood = corrgood1 | corrgood2
            dep_task_wf = pickle.load(open(pkl_path+'%s_%s_adj.pkl' % ('task01','sp'),'rb'))[corrgood]
        if var_def in ('task03','task04'):
            corrarr2 = pickle.load(open(pkl_path+'%s_corrarr.pkl' % ('task02'),'rb')) 
            corrgood2 = corrarr2 > unknown_ratio
            corrgood = corrgood1 | corrgood2
            dep_task_wf = pickle.load(open(pkl_path+'%s_%s_adj.pkl' % ('task02','notedgeon'),'rb'))[corrgood]

        p_el_adj_values = p_el_adj_values[corrgood]
        p_sp_adj_values = p_sp_adj_values[corrgood]
        task_counts_adj = task_counts[corrgood]
        
        taskprob = 0.5
        threshprobarr = (0.5,0.8)
        hilim = 25
        lolim = 10

        if var_def in ('task02','task03','task04'):

            for threshprob in threshprobarr:

                lim = lolim
                print ' '
                print '%s, thresh prob = %4.1f' % (task_dict['var_def'],threshprob)
                assert len(task_counts_adj) == len(dep_task_wf), \
                    'Task counts != dep task'
                assert len(task_counts_adj) == len(p_el_adj_values), \
                    'Task counts != p_el_adj'
                a01 = np.sum((task_counts_adj >= lim) & (dep_task_wf >= taskprob) & (p_el_adj_values >= threshprob))
                a02 = np.sum((task_counts_adj >= lim) & (dep_task_wf >= taskprob) & (p_sp_adj_values >= threshprob))
                t01 = np.sum((task_counts_adj >= lim) & (dep_task_wf >= taskprob) )
                a01_wf = p_el_adj_values[task_counts_adj >= lim]
                a02_wf = p_sp_adj_values[task_counts_adj >= lim]
                a01_wf_ind = np.isfinite(a01_wf)
                a02_wf_ind = np.isfinite(a02_wf)
                print var2_str, '%6i' % a02,'%5.2f' % (float(a02)/t01 * 100.), \
                    '%8i' % np.sum(a02_wf[a02_wf_ind]), '%5.2f' % (np.sum(a02_wf[a02_wf_ind])/len(a02_wf_ind) * 100.)
                print var1_str, '%6i' % a01,'%5.2f' % (float(a01)/t01 * 100.), \
                    '%8i' % np.sum(a01_wf[a01_wf_ind]), '%5.2f' % (np.sum(a01_wf[a01_wf_ind])/len(a01_wf_ind) * 100.)

                if var_def is 'task02':
                    el,sp=p_el_adj_values,p_sp_adj_values

        else:

            for threshprob in threshprobarr:

                lim = hilim
                print ' '
                print '%s, thresh prob = %4.1f' % (task_dict['var_def'],threshprob)
                assert len(task_counts_adj) == len(p_el_adj_values), \
                    'Task counts != p_el_adj'
                a01 = np.sum((task_counts_adj >= lim) & (p_el_adj_values >= threshprob))
                a02 = np.sum((task_counts_adj >= lim) & (p_sp_adj_values >= threshprob))
                t01 = np.sum(task_counts_adj >= lim)
                a01_wf = p_el_adj_values[task_counts_adj >= lim]
                a02_wf = p_sp_adj_values[task_counts_adj >= lim]
                a01_wf_ind = np.isfinite(a01_wf)
                a02_wf_ind = np.isfinite(a02_wf)
                print var2_str, '%6i' % a02,'%5.2f' % (float(a02)/t01 * 100.), \
                    '%8i' % np.sum(a02_wf[a02_wf_ind]), '%5.2f' % (np.sum(a02_wf[a02_wf_ind])/len(a02_wf_ind) * 100.)
                print var1_str, '%6i' % a01,'%5.2f' % (float(a01)/t01 * 100.), \
                    '%8i' % np.sum(a01_wf[a01_wf_ind]), '%5.2f' % (np.sum(a01_wf[a01_wf_ind])/len(a01_wf_ind) * 100.)

    return None

def confidence_measures(task_dict):

    """
    Compute the confidence measures for GZ2 from Lintott et al. (2011)
    """

    var_def = task_dict['var_def']
    centers_redshift, centers_mag, centers_size, edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    corrarr         = pickle.load(open(pkl_path+'%s_corrarr.pkl' % (var_def),'rb')) 
    corrbinarr      = pickle.load(open(pkl_path+'%s_corrbinarr.pkl' % (var_def),'rb')) 
    corrgood = corrarr > unknown_ratio

    good_corr     = corrarr[corrgood]
    good_corrbins = corrbinarr[:,corrgood]

    deltap = []
    sigmap = []

    for idx_z, z, in enumerate(centers_redshift):
        temp_z = good_corrbins[0] == idx_z
        for idx_m, m in enumerate(centers_mag):
            temp_m = good_corrbins[1] == idx_m
            for idx_s, s in enumerate(centers_size):
                temp_s = good_corrbins[2] == idx_s
                p_thisbin = good_corr[temp_z & temp_m & temp_s]
                deltap.append(np.mean(p_thisbin))
                sigmap.append(np.std(p_thisbin))

    return deltap, sigmap 

def plot_confidence_measures(task_dict):

    corrarr         = pickle.load(open(pkl_path+'%s_corrarr.pkl' % (task_dict['var_def']),'rb')) 
    """
    deltap, sigmap = confidence_measures(task_dict)
    """

    fig = plt.figure(16)
    fig.clf()
    ax = fig.add_subplot(111)
    data = corrarr[corrarr > unknown_ratio]
    histML(data, bins=25, ax=ax, histtype='step', color='b',weights=np.zeros_like(data) + 1./data.size, range=(0,1.5))
    #ax.set_xlabel(r'$\langle \Delta p \rangle$')
    ax.set_xlabel(r'$C$')
    ax.set_ylabel('fraction')
    ax.set_title('Corrections for %s' % task_dict['var_def'])
    
    #plt.show()

    return None 

def save_adjusted_probabilities(photoz=False,stripe82=False):

    # Load in the raw probabilities for all GZ2 galaxies with redshifts

    file_str = ''

    if stripe82:

        file_str += 's82_'

        p = pyfits.open(gz2_coadd2_data_file)
        gzdata_all = p[1].data
        p.close()

        if photoz:
            """
            DOES NOT WORK - the PhotoZ table in the Stripe82 database is not populated.
            Use photo-Z from the DR7 table that match the original objIDs.
            Skipped for the moment.
            """
            gzdata = gzdata_all[np.logical_not(np.isfinite(gzdata_all['REDSHIFT']))]
            redshift = gzdata['photoz']
            mr = (gzdata['PETROMAG_R'] - gzdata['EXTINCTION_R']) - cosmology.dmod_flat(redshift)
            r50_kpc = gzdata['PETROR50_R'] * cosmology.ang_scale_flat(redshift)

            file_str += 'photoz_'
        else:
            gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]
            redshift = gzdata['REDSHIFT']
            mr = gzdata['PETROMAG_MR']
            r50_kpc = gzdata['PETROR50_R_KPC']

    else:
        if photoz:
            p = pyfits.open(gz2_photoz_data_file)
            gzdata = p[1].data
            p.close()

            redshift = gzdata['photoz']
            mr = (gzdata['PETROMAG_R'] - gzdata['EXTINCTION_R']) - cosmology.dmod_flat(redshift)
            r50_kpc = gzdata['PETROR50_R'] * cosmology.ang_scale_flat(redshift)

            file_str += 'photoz_'
        else:
            p = pyfits.open(gz2_maglim_data_file)
            gzdata_all = p[1].data
            p.close()

            gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]
            redshift = gzdata['REDSHIFT']
            mr = gzdata['PETROMAG_MR']
            r50_kpc = gzdata['PETROR50_R_KPC']

    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')
    odd = get_task_dict('odd')
    odd_feature = get_task_dict('odd_feature')
    arms_winding = get_task_dict('arms_winding')
    arms_number = get_task_dict('arms_number')
    bulge = get_task_dict('bulge')
    bulge_shape = get_task_dict('bulge_shape')
    rounded = get_task_dict('rounded')

    tasklist = (smooth,edgeon,bar,spiral,odd,odd_feature,arms_winding,arms_number,bulge,bulge_shape,rounded)

    allprobs = np.zeros((1,len(gzdata)))
    allprobs[0,:] = np.arange(len(gzdata))

    for task_dict in tasklist:

        timestart2 = time.time()

        # Loop over all tasks

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        # Take each galaxy, find which bin it corresponds to, adjust probability

        ntask = len(task_dict['task_names_wf'])
        ngals = len(redshift)

        zstep = edges_redshift[1] - edges_redshift[0]
        magstep = edges_mag[1] - edges_mag[0]
        sizestep = edges_size[1] - edges_size[0]

        corr_test, corrmasked_test, ratiomasked_test = determine_baseline_correction(task_dict, 0, 1)
        cshape = corr_test.shape

        allcorr        = np.zeros((ntask,ntask-1,cshape[0],cshape[1],cshape[2]),dtype=float)
        allcorrmasked  = np.zeros_like(allcorr)
        allratiomasked = np.zeros_like(allcorr)

        # Loop over combinations of tasks

        for idx1, var1 in enumerate(np.arange(ntask)):
            setdiff = np.concatenate((np.arange(var1),np.arange(ntask - (var1+1)) + (var1+1) ))
            for idx2, var2 in enumerate(setdiff):
                corr, corrmasked, ratiomasked = determine_baseline_correction(task_dict, var1, var2)
                allcorr[idx1,idx2,:,:,:] = corr
                allcorrmasked[idx1,idx2,:,:,:] = corrmasked
                allratiomasked[idx1,idx2,:,:,:] = ratiomasked

        p_raw = np.zeros((ntask,ngals),dtype=float)
        for idx,task in enumerate(task_dict['task_names_wf']):
            p_raw[idx,:] = gzdata[task]

        p_task_counts = gzdata[task_dict['task_name_count']]

        # Start adjusting

        p_adj = np.zeros((ntask,ngals),dtype=float)
        corrarr = np.zeros((ntask-1,ntask-1,ngals),dtype=float) + unknown_ratio
        corrbinarr = np.zeros((3,ntask-1,ntask-1,ngals),dtype=int) + unknown_ratio
        corrgoodarr = np.zeros(ngals,dtype=bool)

        outsidebinscount = 0
        corrcount = 0
        zerocorrcount = 0
        unknowncorrcount = 0
        vfzerocount = 0
        vfnonzerocount = 0

        # Loop over all galaxies in the GZ catalog with redshifts

        for idx, (z_temp,m_temp,s_temp) in enumerate(zip(redshift,mr,r50_kpc)):

            f = p_raw[:,idx]
            p_temp = np.zeros_like(f)
            zbin = np.abs(z_temp >= edges_redshift).argmin() - 1
            mbin = np.abs(m_temp >= edges_mag).argmin() - 1
            sbin = np.abs(s_temp >= edges_size).argmin() - 1

            if -1 not in (zbin,sbin,mbin):
                applycorr = np.squeeze(allcorrmasked[:,:,zbin,mbin,sbin])

                # Check to see if any corrections exist for this bin
                if np.sum((np.isnan(applycorr)) & (applycorr == unknown_ratio)) == len(applycorr.ravel()):
                    unknowncorrcount += 1
                    p_temp = f
                else:

                    # Bins with masked corrections or NaN set to none

                    applycorr[applycorr == unknown_ratio] = 0.
                    applycorr[np.isnan(applycorr)] = 0.

                    # Check if any corrections are left to be applied

                    if np.sum(applycorr) == 0.:
                        zerocorrcount +=1
                        p_temp = f
                        '''
                        corrarr[:,idx] = unknown_ratio
                        corrbinarr[:,:,idx] = unknown_ratio
                        '''
                    else:

                        # Apply corrections to vote fractions

                        corrcount += 1
                        corrgoodarr[idx] = True
                        for fidx, votefrac in enumerate(f):
                            votefracs = np.array(list(chain.from_iterable(([votefrac], f[:fidx], f[fidx+1:]))))
                            ac = applycorr[fidx,:] if applycorr.ndim > 1 else np.array(applycorr[fidx])
                            if votefracs[0] == 0.:
                                vfzerocount += 1
                                p_temp[fidx] = votefracs[0]
                            else:
                                vfnonzerocount += 1
                                p_temp[fidx] = p_adj_new(votefracs,ac)

                            #corrarr[fidx,:,idx] = applycorr[fidx,:]
                            #corrbinarr[:,fidx,:,idx] = np.resize(np.array((zbin,mbin,sbin)),(3,3))

                    p_adj[:,idx] = p_temp
              
            else: # if galaxy does not appear in the correction area. Shouldn't happen.
                outsidebinscount += 1
                p_adj[:,idx] = f
                '''
                corrarr[:,idx] = unknown_ratio
                corrbinarr[:,:,idx] = unknown_ratio
                '''

        p_both = np.array((p_raw,p_adj))

        timeend2 = time.time()
        print ' '
        print 'For %s: %s' % (task_dict['var_def'],file_str)
        print ' '
        print '%7i galaxies had correction bins outside the binned volume' % outsidebinscount
        print '%7i galaxies had NaN corrections for all pairs of variables' % unknowncorrcount
        print '%7i galaxies had no non-zero corrections in any bin' % zerocorrcount
        print '%7i galaxies had finite corrections available' % corrcount
        print 'Of those: '
        print '   %4.1f percent of tasks were not adjusted because the raw vote fraction was 100 percent' % \
            ((float(vfzerocount)/(corrcount * ntask)) * 100.)
        print '   %4.1f percent of tasks had a non-zero correction applied' % \
            ((float(vfnonzerocount)/(corrcount * ntask)) * 100.)
        print ' '
        print '%7i total galaxies' % p_raw.shape[1]
        print ' '
        print 'Time elapsed to adjust probabilities: %i seconds' % (timeend2 - timestart2)
        print ' '

        allprobs = np.append(allprobs,p_raw,axis=0)
        allprobs = np.append(allprobs,p_adj,axis=0)

    pickle.dump(allprobs,open(pkl_path+'save_%sadjusted_probabilities.pkl' % file_str,'wb'))

    return allprobs

def gz1_comparison():

    # Load the GZ1 and GZ2 debiased data. Matched against each other in TOPCAT

    gz1_gz2_file = fits_path_main+'gz1_gz2_debiased.fits'

    p = pyfits.open(gz1_gz2_file)
    data = p[1].data
    p.close()

    gz1_raw_el = data['P_EL']
    gz1_raw_sp = data['P_CS']
    gz2_raw_el = data['t01_smooth_or_features_a01_smooth_weighted_fraction']
    gz2_raw_sp = data['t01_smooth_or_features_a02_features_or_disk_weighted_fraction']
    gz1_adj_el = data['P_EL_DEBIASED']
    gz1_adj_sp = data['P_CS_DEBIASED']
    gz2_adj_el = data['t01_smooth_or_features_a01_smooth_debiased']
    gz2_adj_sp = data['t01_smooth_or_features_a02_features_or_disk_debiased']
    gz2_raw_spiralstructure = data['t04_spiral_a08_spiral_weighted_fraction']
    gz2_adj_spiralstructure = data['t04_spiral_a08_spiral_debiased']

    gz1_raw_mg = data['P_MG']
    gz2_raw_odd = data['t06_odd_a14_yes_weighted_fraction']
    gz2_raw_mg = data['t08_odd_feature_a24_merger_weighted_fraction']
    gz2_adj_mg = data['t08_odd_feature_a24_merger_debiased']

    # Find the differences between GZ1/GZ2 raw classifications

    data1_el = gz1_raw_el - gz2_raw_el
    data1_sp = gz1_raw_sp - gz2_raw_sp
    data1_el_clean = (gz1_raw_el >= 0.8) & (gz2_raw_el >= 0.8)
    data1_sp_clean = (gz1_raw_sp >= 0.8) & (gz2_raw_sp >= 0.8)

    fig = plt.figure(17)
    fig.clf()

    ax1 = fig.add_subplot(121)
    histML(data1_el, bins=25, ax=ax1, histtype='step', color='r',weights=np.zeros_like(data1_el) + 1./data1_el.size, range=(-1.,1.), linewidth=2, linestyle='dashed')
    histML(data1_sp, bins=25, ax=ax1, histtype='step', color='b',weights=np.zeros_like(data1_sp) + 1./data1_sp.size, range=(-1.,1.), linewidth=2, linestyle='dashed')
    histML(data1_el[data1_el_clean], bins=25, ax=ax1, histtype='step', color='r',weights=np.zeros_like(data1_el[data1_el_clean]) + 1./data1_el[data1_el_clean].size, range=(-1.,1.), linewidth=1, linestyle='solid')
    histML(data1_sp[data1_sp_clean], bins=25, ax=ax1, histtype='step', color='b',weights=np.zeros_like(data1_sp[data1_sp_clean]) + 1./data1_sp[data1_sp_clean].size, range=(-1.,1.), linewidth=1, linestyle='solid')
    ax1.set_xlabel(r'$f_{GZ1} - f_{GZ2}$')
    ax1.set_ylabel('fraction of total sample')
    ax1.set_xlim(-1.0,1.0)
    ax1.set_ylim(0,0.5)
    ax1.set_title('raw votes')

    plt.legend(('el','sp','el > 0.8','sp > 0.8'), 'upper right', shadow=True, fancybox=True)
    
    ax2 = fig.add_subplot(122)
    data2_el = gz1_adj_el - gz2_adj_el
    data2_sp = gz1_adj_sp - gz2_adj_sp
    data2_el_clean = (gz1_adj_el >= 0.8) & (gz2_adj_el >= 0.8)
    data2_sp_clean = (gz1_adj_el >= 0.8) & (gz2_adj_el >= 0.8)
    histML(data2_el, bins=25, ax=ax2, histtype='step', color='r',weights=np.zeros_like(data2_el) + 1./data2_el.size, range=(-1.,1.), linewidth=2, linestyle='dashed')
    histML(data2_sp, bins=25, ax=ax2, histtype='step', color='b',weights=np.zeros_like(data2_sp) + 1./data2_sp.size, range=(-1.,1.), linewidth=2, linestyle='dashed')
    histML(data2_el[data2_el_clean], bins=25, ax=ax2, histtype='step', color='r',weights=np.zeros_like(data2_el[data2_el_clean]) + 1./data2_el[data2_el_clean].size, range=(-1.,1.), linewidth=1, linestyle='solid')
    histML(data2_sp[data2_sp_clean], bins=25, ax=ax2, histtype='step', color='b',weights=np.zeros_like(data2_sp[data2_sp_clean]) + 1./data2_sp[data2_sp_clean].size, range=(-1.,1.), linewidth=1, linestyle='solid')
    ax2.set_xlabel(r'$p_{GZ1} - p_{GZ2}$')
    ax2.set_ylabel('fraction of total sample')
    ax2.set_xlim(-1.0,1.0)
    ax2.set_ylim(0,0.5)
    ax2.set_title('debiased votes')
    
    plt.legend(('el','sp','el > 0.8','sp > 0.8'), 'upper right', shadow=True, fancybox=True)

    fig.savefig(dropbox_figs_path+'gz1_gz2.eps', dpi=200)

    # Try Steven's "trumpet-style" plot

    fig = plt.figure(18,(10,5))
    fig.clf()

    ax1 = fig.add_subplot(121)
    H, xedges, yedges = np.histogram2d(gz1_raw_sp - gz2_raw_sp, gz1_raw_sp, 
                                       range=[[-1., 1.], [0., 1.]], bins=(50, 50))
    xcens = xedges[:-1] + (xedges[1]-xedges[0])/2.
    ycens = yedges[:-1] + (yedges[1]-yedges[0])/2.
    CS = plt.contourf(xcens,ycens,np.log10(H).T,15,cmap=plt.cm.rainbow)
    #cb = plt.colorbar(CS,orientation='vertical')
    #cb.set_label('log '+r'$N_{gal}$',fontsize=16)
    ax1.set_xlabel(r'$f_{GZ1} - f_{GZ2}$',fontsize=20)
    ax1.set_ylabel(r'$f_{sp,GZ1}$',fontsize=20)
    ax1.set_xlim(-1,1)
    ax1.set_ylim(0,1)
    ax1.set_title('raw votes')

    ax3 = fig.add_subplot(122)
    H, xedges, yedges = np.histogram2d(gz1_adj_sp - gz2_adj_sp, gz1_adj_sp, 
                                       range=[[-1., 1.], [0., 1.]], bins=(50, 50))
    xcens = xedges[:-1] + (xedges[1]-xedges[0])/2.
    ycens = yedges[:-1] + (yedges[1]-yedges[0])/2.
    CS = plt.contourf(xcens,ycens,np.log10(H).T,15,cmap=plt.cm.rainbow)
    cb = plt.colorbar(CS,orientation='vertical')
    cb.set_label('log '+r'$N_{gal}$',fontsize=16)
    ax3.set_xlabel(r'$f_{GZ1} - f_{GZ2}$',fontsize=20)
    ax3.set_ylabel(r'$f_{sp,GZ1}$',fontsize=20)
    ax3.set_ylim(0,1)
    ax3.set_xlim(-1,1)
    ax3.set_title('adj votes')

    fig.savefig(dropbox_figs_path+'gz1_gz2_trumpet.eps', dpi=200)

    gz1_el_clean_flag = (data['ELLIPTICAL'] == 1)
    gz1_el_clean_prob = (gz1_adj_el >= 0.8)
    flag_clean_raw_el = (gz2_raw_el[gz1_el_clean_flag] >= 0.8)
    flag_clean_adj_el = (gz2_adj_el[gz1_el_clean_flag] >= 0.8)
    flag_great_raw_el = (gz2_raw_el[gz1_el_clean_flag] >= 0.5)
    flag_great_adj_el = (gz2_adj_el[gz1_el_clean_flag] >= 0.5)
    prob_clean_raw_el = (gz2_raw_el[gz1_el_clean_prob] >= 0.8)
    prob_clean_adj_el = (gz2_adj_el[gz1_el_clean_prob] >= 0.8)
    prob_great_raw_el = (gz2_raw_el[gz1_el_clean_prob] >= 0.5)
    prob_great_adj_el = (gz2_adj_el[gz1_el_clean_prob] >= 0.5)

    gz1_sp_clean_flag = (data['SPIRAL'] == 1)
    gz1_sp_clean_prob = (gz1_adj_sp >= 0.8)
    flag_clean_raw_fd = (gz2_raw_sp[gz1_sp_clean_flag] >= 0.8)
    flag_clean_adj_fd = (gz2_adj_sp[gz1_sp_clean_flag] >= 0.8)
    flag_great_raw_fd = (gz2_raw_sp[gz1_sp_clean_flag] >= 0.5)
    flag_great_adj_fd = (gz2_adj_sp[gz1_sp_clean_flag] >= 0.5)
    prob_clean_raw_fd = (gz2_raw_sp[gz1_sp_clean_prob] >= 0.8)
    prob_clean_adj_fd = (gz2_adj_sp[gz1_sp_clean_prob] >= 0.8)
    prob_great_raw_fd = (gz2_raw_sp[gz1_sp_clean_prob] >= 0.5)
    prob_great_adj_fd = (gz2_adj_sp[gz1_sp_clean_prob] >= 0.5)

    flag_clean_raw_sp = (gz2_raw_spiralstructure[gz1_sp_clean_flag] >= 0.8)
    flag_clean_adj_sp = (gz2_adj_spiralstructure[gz1_sp_clean_flag] >= 0.8)
    flag_great_raw_sp = (gz2_raw_spiralstructure[gz1_sp_clean_flag] >= 0.5)
    flag_great_adj_sp = (gz2_adj_spiralstructure[gz1_sp_clean_flag] >= 0.5)
    prob_clean_raw_sp = (gz2_raw_spiralstructure[gz1_sp_clean_prob] >= 0.8)
    prob_clean_adj_sp = (gz2_adj_spiralstructure[gz1_sp_clean_prob] >= 0.8)
    prob_great_raw_sp = (gz2_raw_spiralstructure[gz1_sp_clean_prob] >= 0.5)
    prob_great_adj_sp = (gz2_adj_spiralstructure[gz1_sp_clean_prob] >= 0.5)

    gz1_mg_clean_prob = (gz1_raw_mg >= 0.6)
    prob_clean_raw_mg = (gz2_raw_mg[gz1_mg_clean_prob] >= 0.50)
    prob_clean_adj_mg = (gz2_adj_mg[gz1_mg_clean_prob] >= 0.50)
    prob_great_raw_mg = (gz2_raw_mg[gz1_mg_clean_prob] >= 0.25)
    prob_great_adj_mg = (gz2_adj_mg[gz1_mg_clean_prob] >= 0.25)


    print ' '
    print '%6i elliptical galaxies flagged in GZ1 sample.' % np.sum(gz1_el_clean_flag)
    print '   %5i (%4.1f) have raw      f_el > 0.8 in the GZ2 sample.' % (np.sum(flag_clean_raw_el),np.sum(flag_clean_raw_el,dtype=float)/np.sum(gz1_el_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_el > 0.8 in the GZ2 sample.' % (np.sum(flag_clean_adj_el),np.sum(flag_clean_adj_el,dtype=float)/np.sum(gz1_el_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have raw      f_el > 0.5 in the GZ2 sample.' % (np.sum(flag_great_raw_el),np.sum(flag_great_raw_el,dtype=float)/np.sum(gz1_el_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_el > 0.5 in the GZ2 sample.' % (np.sum(flag_great_adj_el),np.sum(flag_great_adj_el,dtype=float)/np.sum(gz1_el_clean_flag,dtype=float)*100)
    print '%6i elliptical galaxies with p_el > 0.8 in GZ1 sample.' % np.sum(gz1_el_clean_prob)
    print '   %5i (%4.1f) have raw      f_el > 0.8 in the GZ2 sample.' % (np.sum(prob_clean_raw_el),np.sum(prob_clean_raw_el,dtype=float)/np.sum(gz1_el_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_el > 0.8 in the GZ2 sample.' % (np.sum(prob_clean_adj_el),np.sum(prob_clean_adj_el,dtype=float)/np.sum(gz1_el_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have raw      f_el > 0.5 in the GZ2 sample.' % (np.sum(prob_great_raw_el),np.sum(prob_great_raw_el,dtype=float)/np.sum(gz1_el_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_el > 0.5 in the GZ2 sample.' % (np.sum(prob_great_adj_el),np.sum(prob_great_adj_el,dtype=float)/np.sum(gz1_el_clean_prob,dtype=float)*100)
    print ' '
    print ' '
    print '%6i spiral galaxies flagged in GZ1 sample.' % np.sum(gz1_sp_clean_flag)
    print '   %5i (%4.1f) have raw      f_fd > 0.8 in the GZ2 sample.' % (np.sum(flag_clean_raw_fd),np.sum(flag_clean_raw_fd,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_fd > 0.8 in the GZ2 sample.' % (np.sum(flag_clean_adj_fd),np.sum(flag_clean_adj_fd,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have raw      f_fd > 0.5 in the GZ2 sample.' % (np.sum(flag_great_raw_fd),np.sum(flag_great_raw_fd,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_fd > 0.5 in the GZ2 sample.' % (np.sum(flag_great_adj_fd),np.sum(flag_great_adj_fd,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '%6i spiral galaxies with p_cs > 0.8 in GZ1 sample.' % np.sum(gz1_sp_clean_prob)
    print '   %5i (%4.1f) have raw      f_fd > 0.8 in the GZ2 sample.' % (np.sum(prob_clean_raw_fd),np.sum(prob_clean_raw_fd,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_fd > 0.8 in the GZ2 sample.' % (np.sum(prob_clean_adj_fd),np.sum(prob_clean_adj_fd,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have raw      f_fd > 0.5 in the GZ2 sample.' % (np.sum(prob_great_raw_fd),np.sum(prob_great_raw_fd,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_fd > 0.5 in the GZ2 sample.' % (np.sum(prob_great_adj_fd),np.sum(prob_great_adj_fd,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print ' '
    print '%6i spiral galaxies flagged in GZ1 sample.' % np.sum(gz1_sp_clean_flag)
    print '   %5i (%4.1f) have raw      f_sp > 0.8 in the GZ2 sample.' % (np.sum(flag_clean_raw_sp),np.sum(flag_clean_raw_sp,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_sp > 0.8 in the GZ2 sample.' % (np.sum(flag_clean_adj_sp),np.sum(flag_clean_adj_sp,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have raw      f_sp > 0.5 in the GZ2 sample.' % (np.sum(flag_great_raw_sp),np.sum(flag_great_raw_sp,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_sp > 0.5 in the GZ2 sample.' % (np.sum(flag_great_adj_sp),np.sum(flag_great_adj_sp,dtype=float)/np.sum(gz1_sp_clean_flag,dtype=float)*100)
    print '%6i spiral galaxies with p_cs > 0.8 in GZ1 sample.' % np.sum(gz1_sp_clean_prob)
    print '   %5i (%4.1f) have raw      f_sp > 0.8 in the GZ2 sample.' % (np.sum(prob_clean_raw_sp),np.sum(prob_clean_raw_sp,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_sp > 0.8 in the GZ2 sample.' % (np.sum(prob_clean_adj_sp),np.sum(prob_clean_adj_sp,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have raw      f_sp > 0.5 in the GZ2 sample.' % (np.sum(prob_great_raw_sp),np.sum(prob_great_raw_sp,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_sp > 0.5 in the GZ2 sample.' % (np.sum(prob_great_adj_sp),np.sum(prob_great_adj_sp,dtype=float)/np.sum(gz1_sp_clean_prob,dtype=float)*100)
    print ' '
    print '%6i mergers with p_mg > 0.6 in GZ1 sample.' % np.sum(gz1_mg_clean_prob)
    print '   %5i (%4.1f) have raw      f_mg > 0.50 in the GZ2 sample.' % (np.sum(prob_clean_raw_mg),np.sum(prob_clean_raw_mg,dtype=float)/np.sum(gz1_mg_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_mg > 0.50 in the GZ2 sample.' % (np.sum(prob_clean_adj_mg),np.sum(prob_clean_adj_mg,dtype=float)/np.sum(gz1_mg_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have raw      f_mg > 0.25 in the GZ2 sample.' % (np.sum(prob_great_raw_mg),np.sum(prob_great_raw_mg,dtype=float)/np.sum(gz1_mg_clean_prob,dtype=float)*100)
    print '   %5i (%4.1f) have debiased p_mg > 0.25 in the GZ2 sample.' % (np.sum(prob_great_adj_mg),np.sum(prob_great_adj_mg,dtype=float)/np.sum(gz1_mg_clean_prob,dtype=float)*100)
    print ' '

    return None

def make_tables(makefits=True,stripe82=False,photoz=False,latex=False,imagelist=False):

    """
    Make tables for final data products; includes both FITS files (for site) and formatted LaTeX table samples (for paper)
    """

    file_str = ''

    if stripe82:

        p = pyfits.open(fits_path_main+'gz2_coadd2_table_sample.fits')
        gzdata_all = p[1].data
        p.close()

        gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]

        file_str += 's82_'

    else:
        if photoz:
            p = pyfits.open(fits_path_main+'gz2_photoz_table_sample.fits')
            gzdata = p[1].data
            p.close()

            file_str += 'photoz_'

        else:
            p = pyfits.open(gz2_maglim_data_file)
            gzdata_all = p[1].data
            p.close()
            
            gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]
    
    allprobs = pickle.load(open(pkl_path+'save_%sadjusted_probabilities.pkl' % file_str,'rb')) 

    assert len(gzdata) == allprobs.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'

    raw_t01_a01 = allprobs[ 1,:]
    raw_t01_a02 = allprobs[ 2,:]
    raw_t01_a03 = allprobs[ 3,:]
    adj_t01_a01 = allprobs[ 4,:]
    adj_t01_a02 = allprobs[ 5,:]
    adj_t01_a03 = allprobs[ 6,:]
    raw_t02_a04 = allprobs[ 7,:]
    raw_t02_a05 = allprobs[ 8,:]
    adj_t02_a04 = allprobs[ 9,:]
    adj_t02_a05 = allprobs[10,:]
    raw_t03_a06 = allprobs[11,:]
    raw_t03_a07 = allprobs[12,:]
    adj_t03_a06 = allprobs[13,:]
    adj_t03_a07 = allprobs[14,:]
    raw_t04_a08 = allprobs[15,:]
    raw_t04_a09 = allprobs[16,:]
    adj_t04_a08 = allprobs[17,:]
    adj_t04_a09 = allprobs[18,:]
    raw_t06_a14 = allprobs[19,:]
    raw_t06_a15 = allprobs[20,:]
    adj_t06_a14 = allprobs[21,:]
    adj_t06_a15 = allprobs[22,:]
    raw_t08_a19 = allprobs[23,:]
    raw_t08_a20 = allprobs[24,:]
    raw_t08_a21 = allprobs[25,:]
    raw_t08_a22 = allprobs[26,:]
    raw_t08_a23 = allprobs[27,:]
    raw_t08_a24 = allprobs[28,:]
    raw_t08_a38 = allprobs[29,:]
    adj_t08_a19 = allprobs[30,:]
    adj_t08_a20 = allprobs[31,:]
    adj_t08_a21 = allprobs[32,:]
    adj_t08_a22 = allprobs[33,:]
    adj_t08_a23 = allprobs[34,:]
    adj_t08_a24 = allprobs[35,:]
    adj_t08_a38 = allprobs[36,:]
    raw_t10_a28 = allprobs[37,:]
    raw_t10_a29 = allprobs[38,:]
    raw_t10_a30 = allprobs[39,:]
    adj_t10_a28 = allprobs[40,:]
    adj_t10_a29 = allprobs[41,:]
    adj_t10_a30 = allprobs[42,:]
    raw_t11_a31 = allprobs[43,:]
    raw_t11_a32 = allprobs[44,:]
    raw_t11_a33 = allprobs[45,:]
    raw_t11_a34 = allprobs[46,:]
    raw_t11_a36 = allprobs[47,:]
    raw_t11_a37 = allprobs[48,:]
    adj_t11_a31 = allprobs[49,:]
    adj_t11_a32 = allprobs[50,:]
    adj_t11_a33 = allprobs[51,:]
    adj_t11_a34 = allprobs[52,:]
    adj_t11_a36 = allprobs[53,:]
    adj_t11_a37 = allprobs[54,:]
    raw_t05_a10 = allprobs[55,:]
    raw_t05_a11 = allprobs[56,:]
    raw_t05_a12 = allprobs[57,:]
    raw_t05_a13 = allprobs[58,:]
    adj_t05_a10 = allprobs[59,:]
    adj_t05_a11 = allprobs[60,:]
    adj_t05_a12 = allprobs[61,:]
    adj_t05_a13 = allprobs[62,:]
    raw_t09_a25 = allprobs[63,:]
    raw_t09_a26 = allprobs[64,:]
    raw_t09_a27 = allprobs[65,:]
    adj_t09_a25 = allprobs[66,:]
    adj_t09_a26 = allprobs[67,:]
    adj_t09_a27 = allprobs[68,:]
    raw_t07_a16 = allprobs[69,:]
    raw_t07_a17 = allprobs[70,:]
    raw_t07_a18 = allprobs[71,:]
    adj_t07_a16 = allprobs[72,:]
    adj_t07_a17 = allprobs[73,:]
    adj_t07_a18 = allprobs[74,:]

    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')
    odd = get_task_dict('odd')
    odd_feature = get_task_dict('odd_feature')
    arms_winding = get_task_dict('arms_winding')
    arms_number = get_task_dict('arms_number')
    bulge = get_task_dict('bulge')
    bulge_shape = get_task_dict('bulge_shape')
    rounded = get_task_dict('rounded')

    # Set the flags for each tasks based on their individual paths and parameters

    if stripe82:
        smooth['min_classifications'] = 10
        edgeon['min_classifications'] = 5
        bar['min_classifications'] = 5
        spiral['min_classifications'] = 5
        odd['min_classifications'] = 10
        odd_feature['min_classifications'] = 5
        arms_winding['min_classifications'] = 5
        arms_number['min_classifications'] = 5
        bulge['min_classifications'] = 5
        bulge_shape['min_classifications'] = 5
        rounded['min_classifications'] = 5

    goodt01 = (gzdata['t01_smooth_or_features_total_weight'] >= smooth['min_classifications'])
    goodt02 = (gzdata['t02_edgeon_total_weight'] >= edgeon['min_classifications']) & (gzdata[edgeon['dependent_tasks']] >= 0.5)
    goodt03 = (gzdata['t03_bar_total_weight'] >= bar['min_classifications']) & (gzdata[bar['dependent_tasks'][0]] >= 0.5) & (gzdata[bar['dependent_tasks'][1]] >= 0.5)
    goodt04 = (gzdata['t04_spiral_total_weight'] >= spiral['min_classifications']) & (gzdata[spiral['dependent_tasks'][0]] >= 0.5) & (gzdata[spiral['dependent_tasks'][1]] >= 0.5)
    goodt05 = (gzdata['t05_bulge_prominence_total_weight'] >= bulge['min_classifications']) & (gzdata[bulge['dependent_tasks'][0]] >= 0.5) & (gzdata[bulge['dependent_tasks'][1]] >= 0.5)
    goodt06 = (gzdata['t06_odd_total_weight'] >= odd['min_classifications'])
    goodt07 = (gzdata['t07_rounded_total_weight'] >= rounded['min_classifications']) & (gzdata[rounded['dependent_tasks']] >= 0.5)
    goodt08 = (gzdata['t08_odd_feature_total_weight'] >= odd_feature['min_classifications']) & (gzdata[odd_feature['dependent_tasks']] >= 0.5)
    goodt09 = (gzdata['t09_bulge_shape_total_weight'] >= bulge_shape['min_classifications']) & (gzdata[bulge_shape['dependent_tasks'][0]] >= 0.5) & (gzdata[bulge_shape['dependent_tasks'][1]] >= 0.5)
    goodt10 = (gzdata['t10_arms_winding_total_weight'] >= arms_winding['min_classifications']) & (gzdata[arms_winding['dependent_tasks'][0]] >= 0.5) & (gzdata[arms_winding['dependent_tasks'][1]] >= 0.5) & (gzdata[arms_winding['dependent_tasks'][2]] >= 0.5)
    goodt11 = (gzdata['t11_arms_number_total_weight'] >= arms_number['min_classifications']) & (gzdata[arms_number['dependent_tasks'][0]] >= 0.5) & (gzdata[arms_number['dependent_tasks'][1]] >= 0.5) & (gzdata[arms_number['dependent_tasks'][2]] >= 0.5)

    cleanthresh_t01 = 0.8
    cleanthresh_t02 = 0.5
    cleanthresh_t03 = 0.5
    cleanthresh_t04 = 0.8
    cleanthresh_t05 = 0.8
    cleanthresh_t06 = 0.8
    cleanthresh_t07 = 0.8
    cleanthresh_t08 = 0.8
    cleanthresh_t09 = 0.8
    cleanthresh_t10 = 0.8
    cleanthresh_t11 = 0.8

    flag_t01_a01 = ((allprobs[ 4,:] >= cleanthresh_t01) & goodt01).astype(int)
    flag_t01_a02 = ((allprobs[ 5,:] >= cleanthresh_t01) & goodt01).astype(int)
    flag_t01_a03 = ((allprobs[ 6,:] >= cleanthresh_t01) & goodt01).astype(int)
    flag_t02_a04 = ((allprobs[ 9,:] >= cleanthresh_t02) & goodt02).astype(int)
    flag_t02_a05 = ((allprobs[10,:] >= cleanthresh_t02) & goodt02).astype(int)
    flag_t03_a06 = ((allprobs[13,:] >= cleanthresh_t03) & goodt03).astype(int)
    flag_t03_a07 = ((allprobs[14,:] >= cleanthresh_t03) & goodt03).astype(int)
    flag_t04_a08 = ((allprobs[17,:] >= cleanthresh_t04) & goodt04).astype(int)
    flag_t04_a09 = ((allprobs[18,:] >= cleanthresh_t04) & goodt04).astype(int)
    flag_t06_a14 = ((allprobs[21,:] >= cleanthresh_t06) & goodt06).astype(int)
    flag_t06_a15 = ((allprobs[22,:] >= cleanthresh_t06) & goodt06).astype(int)
    flag_t08_a19 = ((allprobs[30,:] >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a20 = ((allprobs[31,:] >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a21 = ((allprobs[32,:] >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a22 = ((allprobs[33,:] >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a23 = ((allprobs[34,:] >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a24 = ((allprobs[35,:] >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a38 = ((allprobs[36,:] >= cleanthresh_t08) & goodt08).astype(int)
    flag_t10_a28 = ((allprobs[40,:] >= cleanthresh_t10) & goodt10).astype(int)
    flag_t10_a29 = ((allprobs[41,:] >= cleanthresh_t10) & goodt10).astype(int)
    flag_t10_a30 = ((allprobs[42,:] >= cleanthresh_t10) & goodt10).astype(int)
    flag_t11_a31 = ((allprobs[49,:] >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a32 = ((allprobs[50,:] >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a33 = ((allprobs[51,:] >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a34 = ((allprobs[52,:] >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a36 = ((allprobs[53,:] >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a37 = ((allprobs[54,:] >= cleanthresh_t11) & goodt11).astype(int)
    flag_t05_a10 = ((allprobs[59,:] >= cleanthresh_t05) & goodt05).astype(int)
    flag_t05_a11 = ((allprobs[60,:] >= cleanthresh_t05) & goodt05).astype(int)
    flag_t05_a12 = ((allprobs[61,:] >= cleanthresh_t05) & goodt05).astype(int)
    flag_t05_a13 = ((allprobs[62,:] >= cleanthresh_t05) & goodt05).astype(int)
    flag_t09_a25 = ((allprobs[66,:] >= cleanthresh_t09) & goodt09).astype(int)
    flag_t09_a26 = ((allprobs[67,:] >= cleanthresh_t09) & goodt09).astype(int)
    flag_t09_a27 = ((allprobs[68,:] >= cleanthresh_t09) & goodt09).astype(int)
    flag_t07_a16 = ((allprobs[72,:] >= cleanthresh_t07) & goodt07).astype(int)
    flag_t07_a17 = ((allprobs[73,:] >= cleanthresh_t07) & goodt07).astype(int)
    flag_t07_a18 = ((allprobs[74,:] >= cleanthresh_t07) & goodt07).astype(int)

    
    if makefits:
        objid = pyfits.Column(name = 'objid', format='K', array=gzdata['objid'])
        sample = pyfits.Column(name = 'sample', format='A10', array=gzdata['sample'])
        ra = pyfits.Column(name = 'ra', format='E', array=gzdata['RA'])
        dec = pyfits.Column(name = 'dec', format='E', array=gzdata['DEC'])
        n_votes = pyfits.Column(name = 'total_classifications', format='I', array=gzdata['t01_smooth_or_features_total_count'])
        n_class = pyfits.Column(name = 'total_votes', format='I', array=gzdata['total_count'])

        col_t01_a01_ct = pyfits.Column(name='t01_smooth_or_features_a01_smooth_count',                             format='I', array = gzdata['t01_smooth_or_features_a01_smooth_count'])
        col_t01_a01_wc = pyfits.Column(name='t01_smooth_or_features_a01_smooth_weight',                            format='I', array = gzdata['t01_smooth_or_features_a01_smooth_count'])
        col_t01_a01_fr = pyfits.Column(name='t01_smooth_or_features_a01_smooth_fraction',                          format='E', array = gzdata['t01_smooth_or_features_a01_smooth_weighted_fraction'])
        col_t01_a01_wf = pyfits.Column(name='t01_smooth_or_features_a01_smooth_weighted_fraction',                 format='E', array = gzdata['t01_smooth_or_features_a01_smooth_weighted_fraction'])
        col_t01_a01_db = pyfits.Column(name='t01_smooth_or_features_a01_smooth_debiased',                          format='E', array = adj_t01_a01)
        col_t01_a01_fl = pyfits.Column(name='t01_smooth_or_features_a01_smooth_flag',                              format='I', array = flag_t01_a01)
        col_t01_a02_ct = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_count',                   format='I', array = gzdata['t01_smooth_or_features_a02_features_or_disk_count'])
        col_t01_a02_wc = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_weight',                  format='I', array = gzdata['t01_smooth_or_features_a02_features_or_disk_count'])
        col_t01_a02_fr = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_fraction',                format='E', array = gzdata['t01_smooth_or_features_a02_features_or_disk_weighted_fraction'])
        col_t01_a02_wf = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_weighted_fraction',       format='E', array = gzdata['t01_smooth_or_features_a02_features_or_disk_weighted_fraction'])
        col_t01_a02_db = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_debiased',                format='E', array = adj_t01_a02)
        col_t01_a02_fl = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_flag',                    format='I', array = flag_t01_a02)
        col_t01_a03_ct = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_count',                   format='I', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_count'])
        col_t01_a03_wc = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_weight',                  format='I', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_count'])
        col_t01_a03_fr = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_fraction',                format='E', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_weighted_fraction'])
        col_t01_a03_wf = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_weighted_fraction',       format='E', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_weighted_fraction'])
        col_t01_a03_db = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_debiased',                format='E', array = adj_t01_a03)
        col_t01_a03_fl = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_flag',                    format='I', array = flag_t01_a03)

        col_t02_a04_ct = pyfits.Column(name='t02_edgeon_a04_yes_count',                   format='I', array = gzdata['t02_edgeon_a04_yes_count'])
        col_t02_a04_wc = pyfits.Column(name='t02_edgeon_a04_yes_weight',                  format='I', array = gzdata['t02_edgeon_a04_yes_count'])
        col_t02_a04_fr = pyfits.Column(name='t02_edgeon_a04_yes_fraction',                format='E', array = gzdata['t02_edgeon_a04_yes_weighted_fraction'])
        col_t02_a04_wf = pyfits.Column(name='t02_edgeon_a04_yes_weighted_fraction',       format='E', array = gzdata['t02_edgeon_a04_yes_weighted_fraction'])
        col_t02_a04_db = pyfits.Column(name='t02_edgeon_a04_yes_debiased',                format='E', array = adj_t02_a04)
        col_t02_a04_fl = pyfits.Column(name='t02_edgeon_a04_yes_flag',                    format='I', array = flag_t02_a04)
        col_t02_a05_ct = pyfits.Column(name='t02_edgeon_a05_no_count',                    format='I', array = gzdata['t02_edgeon_a05_no_count'])
        col_t02_a05_wc = pyfits.Column(name='t02_edgeon_a05_no_weight',                   format='I', array = gzdata['t02_edgeon_a05_no_count'])
        col_t02_a05_fr = pyfits.Column(name='t02_edgeon_a05_no_fraction',                 format='E', array = gzdata['t02_edgeon_a05_no_weighted_fraction'])
        col_t02_a05_wf = pyfits.Column(name='t02_edgeon_a05_no_weighted_fraction',        format='E', array = gzdata['t02_edgeon_a05_no_weighted_fraction'])
        col_t02_a05_db = pyfits.Column(name='t02_edgeon_a05_no_debiased',                 format='E', array = adj_t02_a05)
        col_t02_a05_fl = pyfits.Column(name='t02_edgeon_a05_no_flag',                     format='I', array = flag_t02_a05)

        col_t03_a06_ct = pyfits.Column(name='t03_bar_a06_bar_count',                      format='I', array = gzdata['t03_bar_a06_bar_count'])
        col_t03_a06_wc = pyfits.Column(name='t03_bar_a06_bar_weight',                     format='I', array = gzdata['t03_bar_a06_bar_count'])
        col_t03_a06_fr = pyfits.Column(name='t03_bar_a06_bar_fraction',                   format='E', array = gzdata['t03_bar_a06_bar_weighted_fraction'])
        col_t03_a06_wf = pyfits.Column(name='t03_bar_a06_bar_weighted_fraction',          format='E', array = gzdata['t03_bar_a06_bar_weighted_fraction'])
        col_t03_a06_db = pyfits.Column(name='t03_bar_a06_bar_debiased',                   format='E', array = adj_t03_a06)
        col_t03_a06_fl = pyfits.Column(name='t03_bar_a06_bar_flag',                       format='I', array = flag_t03_a06)
        col_t03_a07_ct = pyfits.Column(name='t03_bar_a07_no_bar_count',                   format='I', array = gzdata['t03_bar_a07_no_bar_count'])
        col_t03_a07_wc = pyfits.Column(name='t03_bar_a07_no_bar_weight',                  format='I', array = gzdata['t03_bar_a07_no_bar_count'])
        col_t03_a07_fr = pyfits.Column(name='t03_bar_a07_no_bar_fraction',                format='E', array = gzdata['t03_bar_a07_no_bar_weighted_fraction'])
        col_t03_a07_wf = pyfits.Column(name='t03_bar_a07_no_bar_weighted_fraction',       format='E', array = gzdata['t03_bar_a07_no_bar_weighted_fraction'])
        col_t03_a07_db = pyfits.Column(name='t03_bar_a07_no_bar_debiased',                format='E', array = adj_t03_a07)
        col_t03_a07_fl = pyfits.Column(name='t03_bar_a07_no_bar_flag',                    format='I', array = flag_t03_a07)

        col_t04_a08_ct = pyfits.Column(name='t04_spiral_a08_spiral_count',                    format='I', array = gzdata['t04_spiral_a08_spiral_count'])
        col_t04_a08_wc = pyfits.Column(name='t04_spiral_a08_spiral_weight',                   format='I', array = gzdata['t04_spiral_a08_spiral_count'])
        col_t04_a08_fr = pyfits.Column(name='t04_spiral_a08_spiral_fraction',                 format='E', array = gzdata['t04_spiral_a08_spiral_weighted_fraction'])
        col_t04_a08_wf = pyfits.Column(name='t04_spiral_a08_spiral_weighted_fraction',        format='E', array = gzdata['t04_spiral_a08_spiral_weighted_fraction'])
        col_t04_a08_db = pyfits.Column(name='t04_spiral_a08_spiral_debiased',                 format='E', array =     adj_t04_a08)
        col_t04_a08_fl = pyfits.Column(name='t04_spiral_a08_spiral_flag',                     format='I', array =    flag_t04_a08)
        col_t04_a09_ct = pyfits.Column(name='t04_spiral_a09_no_spiral_count',                 format='I', array = gzdata['t04_spiral_a09_no_spiral_count'])
        col_t04_a09_wc = pyfits.Column(name='t04_spiral_a09_no_spiral_weight',                format='I', array = gzdata['t04_spiral_a09_no_spiral_count'])
        col_t04_a09_fr = pyfits.Column(name='t04_spiral_a09_no_spiral_fraction',              format='E', array = gzdata['t04_spiral_a09_no_spiral_weighted_fraction'])
        col_t04_a09_wf = pyfits.Column(name='t04_spiral_a09_no_spiral_weighted_fraction',     format='E', array = gzdata['t04_spiral_a09_no_spiral_weighted_fraction'])
        col_t04_a09_db = pyfits.Column(name='t04_spiral_a09_no_spiral_debiased',              format='E', array =     adj_t04_a09)
        col_t04_a09_fl = pyfits.Column(name='t04_spiral_a09_no_spiral_flag',                  format='I', array =    flag_t04_a09)

        col_t05_a10_ct = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_count',                    format='I', array = gzdata['t05_bulge_prominence_a10_no_bulge_count'])
        col_t05_a10_wc = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_weight',                   format='I', array = gzdata['t05_bulge_prominence_a10_no_bulge_count'])
        col_t05_a10_fr = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_fraction',                 format='E', array = gzdata['t05_bulge_prominence_a10_no_bulge_weighted_fraction'])
        col_t05_a10_wf = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_weighted_fraction',        format='E', array = gzdata['t05_bulge_prominence_a10_no_bulge_weighted_fraction'])
        col_t05_a10_db = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_debiased',                 format='E', array =     adj_t05_a10)
        col_t05_a10_fl = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_flag',                     format='I', array =    flag_t05_a10)
        col_t05_a11_ct = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_count',             format='I', array = gzdata['t05_bulge_prominence_a11_just_noticeable_count'])
        col_t05_a11_wc = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_weight',            format='I', array = gzdata['t05_bulge_prominence_a11_just_noticeable_count'])
        col_t05_a11_fr = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_fraction',          format='E', array = gzdata['t05_bulge_prominence_a11_just_noticeable_weighted_fraction'])
        col_t05_a11_wf = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_weighted_fraction', format='E', array = gzdata['t05_bulge_prominence_a11_just_noticeable_weighted_fraction'])
        col_t05_a11_db = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_debiased',          format='E', array =     adj_t05_a11)
        col_t05_a11_fl = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_flag',              format='I', array =    flag_t05_a11)
        col_t05_a12_ct = pyfits.Column(name='t05_bulge_prominence_a12_obvious_count',                     format='I', array = gzdata['t05_bulge_prominence_a12_obvious_count'])
        col_t05_a12_wc = pyfits.Column(name='t05_bulge_prominence_a12_obvious_weight',                    format='I', array = gzdata['t05_bulge_prominence_a12_obvious_count'])
        col_t05_a12_fr = pyfits.Column(name='t05_bulge_prominence_a12_obvious_fraction',                  format='E', array = gzdata['t05_bulge_prominence_a12_obvious_weighted_fraction'])
        col_t05_a12_wf = pyfits.Column(name='t05_bulge_prominence_a12_obvious_weighted_fraction',         format='E', array = gzdata['t05_bulge_prominence_a12_obvious_weighted_fraction'])
        col_t05_a12_db = pyfits.Column(name='t05_bulge_prominence_a12_obvious_debiased',                  format='E', array =     adj_t05_a12)
        col_t05_a12_fl = pyfits.Column(name='t05_bulge_prominence_a12_obvious_flag',                      format='I', array =    flag_t05_a12)
        col_t05_a13_ct = pyfits.Column(name='t05_bulge_prominence_a13_dominant_count',                    format='I', array = gzdata['t05_bulge_prominence_a13_dominant_count'])
        col_t05_a13_wc = pyfits.Column(name='t05_bulge_prominence_a13_dominant_weight',                   format='I', array = gzdata['t05_bulge_prominence_a13_dominant_count'])
        col_t05_a13_fr = pyfits.Column(name='t05_bulge_prominence_a13_dominant_fraction',                 format='E', array = gzdata['t05_bulge_prominence_a13_dominant_weighted_fraction'])
        col_t05_a13_wf = pyfits.Column(name='t05_bulge_prominence_a13_dominant_weighted_fraction',        format='E', array = gzdata['t05_bulge_prominence_a13_dominant_weighted_fraction'])
        col_t05_a13_db = pyfits.Column(name='t05_bulge_prominence_a13_dominant_debiased',                 format='E', array =     adj_t05_a13)
        col_t05_a13_fl = pyfits.Column(name='t05_bulge_prominence_a13_dominant_flag',                     format='I', array =    flag_t05_a13)

        col_t06_a14_ct = pyfits.Column(name='t06_odd_a14_yes_count',                    format='I', array = gzdata['t06_odd_a14_yes_count'])
        col_t06_a14_wc = pyfits.Column(name='t06_odd_a14_yes_weight',                   format='I', array = gzdata['t06_odd_a14_yes_count'])
        col_t06_a14_fr = pyfits.Column(name='t06_odd_a14_yes_fraction',                 format='E', array = gzdata['t06_odd_a14_yes_weighted_fraction'])
        col_t06_a14_wf = pyfits.Column(name='t06_odd_a14_yes_weighted_fraction',        format='E', array = gzdata['t06_odd_a14_yes_weighted_fraction'])
        col_t06_a14_db = pyfits.Column(name='t06_odd_a14_yes_debiased',                 format='E', array =     adj_t06_a14)
        col_t06_a14_fl = pyfits.Column(name='t06_odd_a14_yes_flag',                     format='I', array =    flag_t06_a14)
        col_t06_a15_ct = pyfits.Column(name='t06_odd_a15_no_count',                     format='I', array = gzdata['t06_odd_a15_no_count'])
        col_t06_a15_wc = pyfits.Column(name='t06_odd_a15_no_weight',                    format='I', array = gzdata['t06_odd_a15_no_count'])
        col_t06_a15_fr = pyfits.Column(name='t06_odd_a15_no_fraction',                  format='E', array = gzdata['t06_odd_a15_no_weighted_fraction'])
        col_t06_a15_wf = pyfits.Column(name='t06_odd_a15_no_weighted_fraction',         format='E', array = gzdata['t06_odd_a15_no_weighted_fraction'])
        col_t06_a15_db = pyfits.Column(name='t06_odd_a15_no_debiased',                  format='E', array =     adj_t06_a15)
        col_t06_a15_fl = pyfits.Column(name='t06_odd_a15_no_flag',                      format='I', array =    flag_t06_a15)


        col_t07_a16_ct = pyfits.Column(name='t07_rounded_a16_completely_round_count',             format='I', array = gzdata['t07_rounded_a16_completely_round_count'])
        col_t07_a16_wc = pyfits.Column(name='t07_rounded_a16_completely_round_weight',            format='I', array = gzdata['t07_rounded_a16_completely_round_count'])
        col_t07_a16_fr = pyfits.Column(name='t07_rounded_a16_completely_round_fraction',          format='E', array = gzdata['t07_rounded_a16_completely_round_weighted_fraction'])
        col_t07_a16_wf = pyfits.Column(name='t07_rounded_a16_completely_round_weighted_fraction', format='E', array = gzdata['t07_rounded_a16_completely_round_weighted_fraction'])
        col_t07_a16_db = pyfits.Column(name='t07_rounded_a16_completely_round_debiased',          format='E', array = adj_t07_a16)
        col_t07_a16_fl = pyfits.Column(name='t07_rounded_a16_completely_round_flag',              format='I', array = flag_t07_a16)
        col_t07_a17_ct = pyfits.Column(name='t07_rounded_a17_in_between_count',                   format='I', array = gzdata['t07_rounded_a17_in_between_count'])
        col_t07_a17_wc = pyfits.Column(name='t07_rounded_a17_in_between_weight',                  format='I', array = gzdata['t07_rounded_a17_in_between_count'])
        col_t07_a17_fr = pyfits.Column(name='t07_rounded_a17_in_between_fraction',                format='E', array = gzdata['t07_rounded_a17_in_between_weighted_fraction'])
        col_t07_a17_wf = pyfits.Column(name='t07_rounded_a17_in_between_weighted_fraction',       format='E', array = gzdata['t07_rounded_a17_in_between_weighted_fraction'])
        col_t07_a17_db = pyfits.Column(name='t07_rounded_a17_in_between_debiased',                format='E', array = adj_t07_a17)
        col_t07_a17_fl = pyfits.Column(name='t07_rounded_a17_in_between_flag',                    format='I', array = flag_t07_a17)
        col_t07_a18_ct = pyfits.Column(name='t07_rounded_a18_cigar_shaped_count',                 format='I', array = gzdata['t07_rounded_a18_cigar_shaped_count'])
        col_t07_a18_wc = pyfits.Column(name='t07_rounded_a18_cigar_shaped_weight',                format='I', array = gzdata['t07_rounded_a18_cigar_shaped_count'])
        col_t07_a18_fr = pyfits.Column(name='t07_rounded_a18_cigar_shaped_fraction',              format='E', array = gzdata['t07_rounded_a18_cigar_shaped_weighted_fraction'])
        col_t07_a18_wf = pyfits.Column(name='t07_rounded_a18_cigar_shaped_weighted_fraction',     format='E', array = gzdata['t07_rounded_a18_cigar_shaped_weighted_fraction'])
        col_t07_a18_db = pyfits.Column(name='t07_rounded_a18_cigar_shaped_debiased',              format='E', array = adj_t07_a18)
        col_t07_a18_fl = pyfits.Column(name='t07_rounded_a18_cigar_shaped_flag',                  format='I', array = flag_t07_a18)

        col_t08_a19_ct = pyfits.Column(name='t08_odd_feature_a19_ring_count',                             format='I', array = gzdata['t08_odd_feature_a19_ring_count'])
        col_t08_a19_wc = pyfits.Column(name='t08_odd_feature_a19_ring_weight',                            format='I', array = gzdata['t08_odd_feature_a19_ring_count'])
        col_t08_a19_fr = pyfits.Column(name='t08_odd_feature_a19_ring_fraction',                          format='E', array = gzdata['t08_odd_feature_a19_ring_weighted_fraction'])
        col_t08_a19_wf = pyfits.Column(name='t08_odd_feature_a19_ring_weighted_fraction',                 format='E', array = gzdata['t08_odd_feature_a19_ring_weighted_fraction'])
        col_t08_a19_db = pyfits.Column(name='t08_odd_feature_a19_ring_debiased',                          format='E', array = adj_t08_a19)    
        col_t08_a19_fl = pyfits.Column(name='t08_odd_feature_a19_ring_flag',                              format='I', array = flag_t08_a19)  
        col_t08_a20_ct = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_count',                      format='I', array = gzdata['t08_odd_feature_a20_lens_or_arc_count'])
        col_t08_a20_wc = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_weight',                     format='I', array = gzdata['t08_odd_feature_a20_lens_or_arc_count'])
        col_t08_a20_fr = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_fraction',                   format='E', array = gzdata['t08_odd_feature_a20_lens_or_arc_weighted_fraction'])
        col_t08_a20_wf = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_weighted_fraction',          format='E', array = gzdata['t08_odd_feature_a20_lens_or_arc_weighted_fraction'])
        col_t08_a20_db = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_debiased',                   format='E', array = adj_t08_a20)  
        col_t08_a20_fl = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_flag',                       format='I', array = flag_t08_a20)
        col_t08_a21_ct = pyfits.Column(name='t08_odd_feature_a21_disturbed_count',                        format='I', array = gzdata['t08_odd_feature_a21_disturbed_count'])
        col_t08_a21_wc = pyfits.Column(name='t08_odd_feature_a21_disturbed_weight',                       format='I', array = gzdata['t08_odd_feature_a21_disturbed_count'])
        col_t08_a21_fr = pyfits.Column(name='t08_odd_feature_a21_disturbed_fraction',                     format='E', array = gzdata['t08_odd_feature_a21_disturbed_weighted_fraction'])
        col_t08_a21_wf = pyfits.Column(name='t08_odd_feature_a21_disturbed_weighted_fraction',            format='E', array = gzdata['t08_odd_feature_a21_disturbed_weighted_fraction'])
        col_t08_a21_db = pyfits.Column(name='t08_odd_feature_a21_disturbed_debiased',                     format='E', array = adj_t08_a21)      
        col_t08_a21_fl = pyfits.Column(name='t08_odd_feature_a21_disturbed_flag',                         format='I', array = flag_t08_a21)    
        col_t08_a22_ct = pyfits.Column(name='t08_odd_feature_a22_irregular_count',                        format='I', array = gzdata['t08_odd_feature_a22_irregular_count'])
        col_t08_a22_wc = pyfits.Column(name='t08_odd_feature_a22_irregular_weight',                       format='I', array = gzdata['t08_odd_feature_a22_irregular_count'])
        col_t08_a22_fr = pyfits.Column(name='t08_odd_feature_a22_irregular_fraction',                     format='E', array = gzdata['t08_odd_feature_a22_irregular_weighted_fraction'])
        col_t08_a22_wf = pyfits.Column(name='t08_odd_feature_a22_irregular_weighted_fraction',            format='E', array = gzdata['t08_odd_feature_a22_irregular_weighted_fraction'])
        col_t08_a22_db = pyfits.Column(name='t08_odd_feature_a22_irregular_debiased',                     format='E', array = adj_t08_a22)    
        col_t08_a22_fl = pyfits.Column(name='t08_odd_feature_a22_irregular_flag',                         format='I', array = flag_t08_a22)  
        col_t08_a23_ct = pyfits.Column(name='t08_odd_feature_a23_other_count',                            format='I', array = gzdata['t08_odd_feature_a23_other_count'])
        col_t08_a23_wc = pyfits.Column(name='t08_odd_feature_a23_other_weight',                           format='I', array = gzdata['t08_odd_feature_a23_other_count'])
        col_t08_a23_fr = pyfits.Column(name='t08_odd_feature_a23_other_fraction',                         format='E', array = gzdata['t08_odd_feature_a23_other_weighted_fraction'])
        col_t08_a23_wf = pyfits.Column(name='t08_odd_feature_a23_other_weighted_fraction',                format='E', array = gzdata['t08_odd_feature_a23_other_weighted_fraction'])
        col_t08_a23_db = pyfits.Column(name='t08_odd_feature_a23_other_debiased',                         format='E', array = adj_t08_a23)       
        col_t08_a23_fl = pyfits.Column(name='t08_odd_feature_a23_other_flag',                             format='I', array = flag_t08_a23)     
        col_t08_a24_ct = pyfits.Column(name='t08_odd_feature_a24_merger_count',                           format='I', array = gzdata['t08_odd_feature_a24_merger_count'])
        col_t08_a24_wc = pyfits.Column(name='t08_odd_feature_a24_merger_weight',                          format='I', array = gzdata['t08_odd_feature_a24_merger_count'])
        col_t08_a24_fr = pyfits.Column(name='t08_odd_feature_a24_merger_fraction',                        format='E', array = gzdata['t08_odd_feature_a24_merger_weighted_fraction'])
        col_t08_a24_wf = pyfits.Column(name='t08_odd_feature_a24_merger_weighted_fraction',               format='E', array = gzdata['t08_odd_feature_a24_merger_weighted_fraction'])
        col_t08_a24_db = pyfits.Column(name='t08_odd_feature_a24_merger_debiased',                        format='E', array = adj_t08_a24)     
        col_t08_a24_fl = pyfits.Column(name='t08_odd_feature_a24_merger_flag',                            format='I', array = flag_t08_a24)   
        col_t08_a38_ct = pyfits.Column(name='t08_odd_feature_a38_dust_lane_count',                        format='I', array = gzdata['t08_odd_feature_a38_dust_lane_count'])
        col_t08_a38_wc = pyfits.Column(name='t08_odd_feature_a38_dust_lane_weight',                       format='I', array = gzdata['t08_odd_feature_a38_dust_lane_count'])
        col_t08_a38_fr = pyfits.Column(name='t08_odd_feature_a38_dust_lane_fraction',                     format='E', array = gzdata['t08_odd_feature_a38_dust_lane_weighted_fraction'])
        col_t08_a38_wf = pyfits.Column(name='t08_odd_feature_a38_dust_lane_weighted_fraction',            format='E', array = gzdata['t08_odd_feature_a38_dust_lane_weighted_fraction'])
        col_t08_a38_db = pyfits.Column(name='t08_odd_feature_a38_dust_lane_debiased',                     format='E', array = adj_t08_a38)   
        col_t08_a38_fl = pyfits.Column(name='t08_odd_feature_a38_dust_lane_flag',                         format='I', array = flag_t08_a38) 

        col_t09_a25_ct = pyfits.Column(name='t09_bulge_shape_a25_rounded_count',                  format='I', array = gzdata['t09_bulge_shape_a25_rounded_count'])
        col_t09_a25_wc = pyfits.Column(name='t09_bulge_shape_a25_rounded_weight',                 format='I', array = gzdata['t09_bulge_shape_a25_rounded_count'])
        col_t09_a25_fr = pyfits.Column(name='t09_bulge_shape_a25_rounded_fraction',               format='E', array = gzdata['t09_bulge_shape_a25_rounded_weighted_fraction'])
        col_t09_a25_wf = pyfits.Column(name='t09_bulge_shape_a25_rounded_weighted_fraction',      format='E', array = gzdata['t09_bulge_shape_a25_rounded_weighted_fraction'])
        col_t09_a25_db = pyfits.Column(name='t09_bulge_shape_a25_rounded_debiased',               format='E', array = adj_t09_a25)
        col_t09_a25_fl = pyfits.Column(name='t09_bulge_shape_a25_rounded_flag',                   format='I', array = flag_t09_a25)
        col_t09_a26_ct = pyfits.Column(name='t09_bulge_shape_a26_boxy_count',                     format='I', array = gzdata['t09_bulge_shape_a26_boxy_count'])
        col_t09_a26_wc = pyfits.Column(name='t09_bulge_shape_a26_boxy_weight',                    format='I', array = gzdata['t09_bulge_shape_a26_boxy_count'])
        col_t09_a26_fr = pyfits.Column(name='t09_bulge_shape_a26_boxy_fraction',                  format='E', array = gzdata['t09_bulge_shape_a26_boxy_weighted_fraction'])
        col_t09_a26_wf = pyfits.Column(name='t09_bulge_shape_a26_boxy_weighted_fraction',         format='E', array = gzdata['t09_bulge_shape_a26_boxy_weighted_fraction'])
        col_t09_a26_db = pyfits.Column(name='t09_bulge_shape_a26_boxy_debiased',                  format='E', array = adj_t09_a26)
        col_t09_a26_fl = pyfits.Column(name='t09_bulge_shape_a26_boxy_flag',                      format='I', array = flag_t09_a26)
        col_t09_a27_ct = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_count',                 format='I', array = gzdata['t09_bulge_shape_a27_no_bulge_count'])
        col_t09_a27_wc = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_weight',                format='I', array = gzdata['t09_bulge_shape_a27_no_bulge_count'])
        col_t09_a27_fr = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_fraction',              format='E', array = gzdata['t09_bulge_shape_a27_no_bulge_weighted_fraction'])
        col_t09_a27_wf = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_weighted_fraction',     format='E', array = gzdata['t09_bulge_shape_a27_no_bulge_weighted_fraction'])
        col_t09_a27_db = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_debiased',              format='E', array = adj_t09_a27)
        col_t09_a27_fl = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_flag',                  format='I', array = flag_t09_a27)

        col_t10_a28_ct = pyfits.Column(name='t10_arms_winding_a28_tight_count',                   format='I', array = gzdata['t10_arms_winding_a28_tight_count'])
        col_t10_a28_wc = pyfits.Column(name='t10_arms_winding_a28_tight_weight',                  format='I', array = gzdata['t10_arms_winding_a28_tight_count'])
        col_t10_a28_fr = pyfits.Column(name='t10_arms_winding_a28_tight_fraction',                format='E', array = gzdata['t10_arms_winding_a28_tight_weighted_fraction'])
        col_t10_a28_wf = pyfits.Column(name='t10_arms_winding_a28_tight_weighted_fraction',       format='E', array = gzdata['t10_arms_winding_a28_tight_weighted_fraction'])
        col_t10_a28_db = pyfits.Column(name='t10_arms_winding_a28_tight_debiased',                format='E', array = adj_t10_a28)
        col_t10_a28_fl = pyfits.Column(name='t10_arms_winding_a28_tight_flag',                    format='I', array = flag_t10_a28)
        col_t10_a29_ct = pyfits.Column(name='t10_arms_winding_a29_medium_count',                  format='I', array = gzdata['t10_arms_winding_a29_medium_count'])
        col_t10_a29_wc = pyfits.Column(name='t10_arms_winding_a29_medium_weight',                 format='I', array = gzdata['t10_arms_winding_a29_medium_count'])
        col_t10_a29_fr = pyfits.Column(name='t10_arms_winding_a29_medium_fraction',               format='E', array = gzdata['t10_arms_winding_a29_medium_weighted_fraction'])
        col_t10_a29_wf = pyfits.Column(name='t10_arms_winding_a29_medium_weighted_fraction',      format='E', array = gzdata['t10_arms_winding_a29_medium_weighted_fraction'])
        col_t10_a29_db = pyfits.Column(name='t10_arms_winding_a29_medium_debiased',               format='E', array = adj_t10_a29)
        col_t10_a29_fl = pyfits.Column(name='t10_arms_winding_a29_medium_flag',                   format='I', array = flag_t10_a29)
        col_t10_a30_ct = pyfits.Column(name='t10_arms_winding_a30_loose_count',                   format='I', array = gzdata['t10_arms_winding_a30_loose_count'])
        col_t10_a30_wc = pyfits.Column(name='t10_arms_winding_a30_loose_weight',                  format='I', array = gzdata['t10_arms_winding_a30_loose_count'])
        col_t10_a30_fr = pyfits.Column(name='t10_arms_winding_a30_loose_fraction',                format='E', array = gzdata['t10_arms_winding_a30_loose_weighted_fraction'])
        col_t10_a30_wf = pyfits.Column(name='t10_arms_winding_a30_loose_weighted_fraction',       format='E', array = gzdata['t10_arms_winding_a30_loose_weighted_fraction'])
        col_t10_a30_db = pyfits.Column(name='t10_arms_winding_a30_loose_debiased',                format='E', array = adj_t10_a30)
        col_t10_a30_fl = pyfits.Column(name='t10_arms_winding_a30_loose_flag',                    format='I', array = flag_t10_a30)

        col_t11_a31_ct = pyfits.Column(name='t11_arms_number_a31_1_count',                        format='I', array = gzdata['t11_arms_number_a31_1_count'])
        col_t11_a31_wc = pyfits.Column(name='t11_arms_number_a31_1_weight',                       format='I', array = gzdata['t11_arms_number_a31_1_count'])
        col_t11_a31_fr = pyfits.Column(name='t11_arms_number_a31_1_fraction',                     format='E', array = gzdata['t11_arms_number_a31_1_weighted_fraction'])
        col_t11_a31_wf = pyfits.Column(name='t11_arms_number_a31_1_weighted_fraction',            format='E', array = gzdata['t11_arms_number_a31_1_weighted_fraction'])
        col_t11_a31_db = pyfits.Column(name='t11_arms_number_a31_1_debiased',                     format='E', array = adj_t11_a31)    
        col_t11_a31_fl = pyfits.Column(name='t11_arms_number_a31_1_flag',                         format='I', array = flag_t11_a31)  
        col_t11_a32_ct = pyfits.Column(name='t11_arms_number_a32_2_count',                        format='I', array = gzdata['t11_arms_number_a32_2_count'])
        col_t11_a32_wc = pyfits.Column(name='t11_arms_number_a32_2_weight',                       format='I', array = gzdata['t11_arms_number_a32_2_count'])
        col_t11_a32_fr = pyfits.Column(name='t11_arms_number_a32_2_fraction',                     format='E', array = gzdata['t11_arms_number_a32_2_weighted_fraction'])
        col_t11_a32_wf = pyfits.Column(name='t11_arms_number_a32_2_weighted_fraction',            format='E', array = gzdata['t11_arms_number_a32_2_weighted_fraction'])
        col_t11_a32_db = pyfits.Column(name='t11_arms_number_a32_2_debiased',                     format='E', array = adj_t11_a32)  
        col_t11_a32_fl = pyfits.Column(name='t11_arms_number_a32_2_flag',                         format='I', array = flag_t11_a32)
        col_t11_a33_ct = pyfits.Column(name='t11_arms_number_a33_3_count',                        format='I', array = gzdata['t11_arms_number_a33_3_count'])
        col_t11_a33_wc = pyfits.Column(name='t11_arms_number_a33_3_weight',                       format='I', array = gzdata['t11_arms_number_a33_3_count'])
        col_t11_a33_fr = pyfits.Column(name='t11_arms_number_a33_3_fraction',                     format='E', array = gzdata['t11_arms_number_a33_3_weighted_fraction'])
        col_t11_a33_wf = pyfits.Column(name='t11_arms_number_a33_3_weighted_fraction',            format='E', array = gzdata['t11_arms_number_a33_3_weighted_fraction'])
        col_t11_a33_db = pyfits.Column(name='t11_arms_number_a33_3_debiased',                     format='E', array = adj_t11_a33)      
        col_t11_a33_fl = pyfits.Column(name='t11_arms_number_a33_3_flag',                         format='I', array = flag_t11_a33)    
        col_t11_a34_ct = pyfits.Column(name='t11_arms_number_a34_4_count',                        format='I', array = gzdata['t11_arms_number_a34_4_count'])
        col_t11_a34_wc = pyfits.Column(name='t11_arms_number_a34_4_weight',                       format='I', array = gzdata['t11_arms_number_a34_4_count'])
        col_t11_a34_fr = pyfits.Column(name='t11_arms_number_a34_4_fraction',                     format='E', array = gzdata['t11_arms_number_a34_4_weighted_fraction'])
        col_t11_a34_wf = pyfits.Column(name='t11_arms_number_a34_4_weighted_fraction',            format='E', array = gzdata['t11_arms_number_a34_4_weighted_fraction'])
        col_t11_a34_db = pyfits.Column(name='t11_arms_number_a34_4_debiased',                     format='E', array = adj_t11_a34)    
        col_t11_a34_fl = pyfits.Column(name='t11_arms_number_a34_4_flag',                         format='I', array = flag_t11_a34)  
        col_t11_a36_ct = pyfits.Column(name='t11_arms_number_a36_more_than_4_count',              format='I', array = gzdata['t11_arms_number_a36_more_than_4_count'])
        col_t11_a36_wc = pyfits.Column(name='t11_arms_number_a36_more_than_4_weight',             format='I', array = gzdata['t11_arms_number_a36_more_than_4_count'])
        col_t11_a36_fr = pyfits.Column(name='t11_arms_number_a36_more_than_4_fraction',           format='E', array = gzdata['t11_arms_number_a36_more_than_4_weighted_fraction'])
        col_t11_a36_wf = pyfits.Column(name='t11_arms_number_a36_more_than_4_weighted_fraction',  format='E', array = gzdata['t11_arms_number_a36_more_than_4_weighted_fraction'])
        col_t11_a36_db = pyfits.Column(name='t11_arms_number_a36_more_than_4_debiased',           format='E', array = adj_t11_a36)       
        col_t11_a36_fl = pyfits.Column(name='t11_arms_number_a36_more_than_4_flag',               format='I', array = flag_t11_a36)     
        col_t11_a37_ct = pyfits.Column(name='t11_arms_number_a37_cant_tell_count',                format='I', array = gzdata['t11_arms_number_a37_cant_tell_count'])
        col_t11_a37_wc = pyfits.Column(name='t11_arms_number_a37_cant_tell_weight',               format='I', array = gzdata['t11_arms_number_a37_cant_tell_count'])
        col_t11_a37_fr = pyfits.Column(name='t11_arms_number_a37_cant_tell_fraction',             format='E', array = gzdata['t11_arms_number_a37_cant_tell_weighted_fraction'])
        col_t11_a37_wf = pyfits.Column(name='t11_arms_number_a37_cant_tell_weighted_fraction',    format='E', array = gzdata['t11_arms_number_a37_cant_tell_weighted_fraction'])
        col_t11_a37_db = pyfits.Column(name='t11_arms_number_a37_cant_tell_debiased',             format='E', array = adj_t11_a37)     
        col_t11_a37_fl = pyfits.Column(name='t11_arms_number_a37_cant_tell_flag',                 format='I', array = flag_t11_a37)   


        primary_hdu = pyfits.PrimaryHDU()
        hdulist = pyfits.HDUList([primary_hdu])
        
        tb1_hdu = pyfits.new_table([\
                                    objid,  
                                    sample,
                                    ra,
                                    dec,
                                    n_class,
                                    n_votes,
                                    col_t01_a01_ct, 
                                    col_t01_a01_wc, 
                                    col_t01_a01_fr, 
                                    col_t01_a01_wf, 
                                    col_t01_a01_db, 
                                    col_t01_a01_fl, 
                                    col_t01_a02_ct, 
                                    col_t01_a02_wc, 
                                    col_t01_a02_fr, 
                                    col_t01_a02_wf, 
                                    col_t01_a02_db, 
                                    col_t01_a02_fl, 
                                    col_t01_a03_ct, 
                                    col_t01_a03_wc, 
                                    col_t01_a03_fr, 
                                    col_t01_a03_wf, 
                                    col_t01_a03_db, 
                                    col_t01_a03_fl, 
                                    col_t02_a04_ct, 
                                    col_t02_a04_wc, 
                                    col_t02_a04_fr, 
                                    col_t02_a04_wf, 
                                    col_t02_a04_db, 
                                    col_t02_a04_fl, 
                                    col_t02_a05_ct, 
                                    col_t02_a05_wc, 
                                    col_t02_a05_fr, 
                                    col_t02_a05_wf, 
                                    col_t02_a05_db, 
                                    col_t02_a05_fl, 
                                    col_t03_a06_ct, 
                                    col_t03_a06_wc, 
                                    col_t03_a06_fr, 
                                    col_t03_a06_wf, 
                                    col_t03_a06_db, 
                                    col_t03_a06_fl, 
                                    col_t03_a07_ct, 
                                    col_t03_a07_wc, 
                                    col_t03_a07_fr, 
                                    col_t03_a07_wf, 
                                    col_t03_a07_db, 
                                    col_t03_a07_fl, 
                                    col_t04_a08_ct, 
                                    col_t04_a08_wc, 
                                    col_t04_a08_fr, 
                                    col_t04_a08_wf, 
                                    col_t04_a08_db, 
                                    col_t04_a08_fl, 
                                    col_t04_a09_ct, 
                                    col_t04_a09_wc, 
                                    col_t04_a09_fr, 
                                    col_t04_a09_wf, 
                                    col_t04_a09_db, 
                                    col_t04_a09_fl, 
                                    col_t05_a10_ct, 
                                    col_t05_a10_wc, 
                                    col_t05_a10_fr, 
                                    col_t05_a10_wf, 
                                    col_t05_a10_db, 
                                    col_t05_a10_fl, 
                                    col_t05_a11_ct, 
                                    col_t05_a11_wc, 
                                    col_t05_a11_fr, 
                                    col_t05_a11_wf, 
                                    col_t05_a11_db, 
                                    col_t05_a11_fl, 
                                    col_t05_a12_ct, 
                                    col_t05_a12_wc, 
                                    col_t05_a12_fr, 
                                    col_t05_a12_wf, 
                                    col_t05_a12_db, 
                                    col_t05_a12_fl, 
                                    col_t05_a13_ct, 
                                    col_t05_a13_wc, 
                                    col_t05_a13_fr, 
                                    col_t05_a13_wf, 
                                    col_t05_a13_db, 
                                    col_t05_a13_fl, 
                                    col_t06_a14_ct, 
                                    col_t06_a14_wc, 
                                    col_t06_a14_fr, 
                                    col_t06_a14_wf, 
                                    col_t06_a14_db, 
                                    col_t06_a14_fl, 
                                    col_t06_a15_ct, 
                                    col_t06_a15_wc, 
                                    col_t06_a15_fr, 
                                    col_t06_a15_wf, 
                                    col_t06_a15_db, 
                                    col_t06_a15_fl, 
                                    col_t07_a16_ct, 
                                    col_t07_a16_wc, 
                                    col_t07_a16_fr, 
                                    col_t07_a16_wf, 
                                    col_t07_a16_db, 
                                    col_t07_a16_fl, 
                                    col_t07_a17_ct, 
                                    col_t07_a17_wc, 
                                    col_t07_a17_fr, 
                                    col_t07_a17_wf, 
                                    col_t07_a17_db, 
                                    col_t07_a17_fl, 
                                    col_t07_a18_ct, 
                                    col_t07_a18_wc, 
                                    col_t07_a18_fr, 
                                    col_t07_a18_wf, 
                                    col_t07_a18_db, 
                                    col_t07_a18_fl, 
                                    col_t08_a19_ct, 
                                    col_t08_a19_wc, 
                                    col_t08_a19_fr, 
                                    col_t08_a19_wf, 
                                    col_t08_a19_db, 
                                    col_t08_a19_fl, 
                                    col_t08_a20_ct, 
                                    col_t08_a20_wc, 
                                    col_t08_a20_fr, 
                                    col_t08_a20_wf, 
                                    col_t08_a20_db, 
                                    col_t08_a20_fl, 
                                    col_t08_a21_ct, 
                                    col_t08_a21_wc, 
                                    col_t08_a21_fr, 
                                    col_t08_a21_wf, 
                                    col_t08_a21_db, 
                                    col_t08_a21_fl, 
                                    col_t08_a22_ct, 
                                    col_t08_a22_wc, 
                                    col_t08_a22_fr, 
                                    col_t08_a22_wf, 
                                    col_t08_a22_db, 
                                    col_t08_a22_fl, 
                                    col_t08_a23_ct, 
                                    col_t08_a23_wc, 
                                    col_t08_a23_fr, 
                                    col_t08_a23_wf, 
                                    col_t08_a23_db, 
                                    col_t08_a23_fl, 
                                    col_t08_a24_ct, 
                                    col_t08_a24_wc, 
                                    col_t08_a24_fr, 
                                    col_t08_a24_wf, 
                                    col_t08_a24_db, 
                                    col_t08_a24_fl, 
                                    col_t08_a38_ct, 
                                    col_t08_a38_wc, 
                                    col_t08_a38_fr, 
                                    col_t08_a38_wf, 
                                    col_t08_a38_db, 
                                    col_t08_a38_fl, 
                                    col_t09_a25_ct, 
                                    col_t09_a25_wc, 
                                    col_t09_a25_fr, 
                                    col_t09_a25_wf, 
                                    col_t09_a25_db, 
                                    col_t09_a25_fl, 
                                    col_t09_a26_ct, 
                                    col_t09_a26_wc, 
                                    col_t09_a26_fr, 
                                    col_t09_a26_wf, 
                                    col_t09_a26_db, 
                                    col_t09_a26_fl, 
                                    col_t09_a27_ct, 
                                    col_t09_a27_wc, 
                                    col_t09_a27_fr, 
                                    col_t09_a27_wf, 
                                    col_t09_a27_db, 
                                    col_t09_a27_fl, 
                                    col_t10_a28_ct, 
                                    col_t10_a28_wc, 
                                    col_t10_a28_fr, 
                                    col_t10_a28_wf, 
                                    col_t10_a28_db, 
                                    col_t10_a28_fl, 
                                    col_t10_a29_ct, 
                                    col_t10_a29_wc, 
                                    col_t10_a29_fr, 
                                    col_t10_a29_wf, 
                                    col_t10_a29_db, 
                                    col_t10_a29_fl, 
                                    col_t10_a30_ct, 
                                    col_t10_a30_wc, 
                                    col_t10_a30_fr, 
                                    col_t10_a30_wf, 
                                    col_t10_a30_db, 
                                    col_t10_a30_fl, 
                                    col_t11_a31_ct, 
                                    col_t11_a31_wc, 
                                    col_t11_a31_fr, 
                                    col_t11_a31_wf, 
                                    col_t11_a31_db, 
                                    col_t11_a31_fl, 
                                    col_t11_a32_ct, 
                                    col_t11_a32_wc, 
                                    col_t11_a32_fr, 
                                    col_t11_a32_wf, 
                                    col_t11_a32_db, 
                                    col_t11_a32_fl, 
                                    col_t11_a33_ct, 
                                    col_t11_a33_wc, 
                                    col_t11_a33_fr, 
                                    col_t11_a33_wf, 
                                    col_t11_a33_db, 
                                    col_t11_a33_fl, 
                                    col_t11_a34_ct, 
                                    col_t11_a34_wc, 
                                    col_t11_a34_fr, 
                                    col_t11_a34_wf, 
                                    col_t11_a34_db, 
                                    col_t11_a34_fl, 
                                    col_t11_a36_ct, 
                                    col_t11_a36_wc, 
                                    col_t11_a36_fr, 
                                    col_t11_a36_wf, 
                                    col_t11_a36_db, 
                                    col_t11_a36_fl, 
                                    col_t11_a37_ct, 
                                    col_t11_a37_wc, 
                                    col_t11_a37_fr, 
                                    col_t11_a37_wf, 
                                    col_t11_a37_db, 
                                    col_t11_a37_fl, 
                                    ])
        tb1_hdu.name = 'GZ2DATA'
        
        hdulist.append(tb1_hdu)
        hdulist.writeto(fits_path_main+'gz2table_%sdebiased.fits' % file_str,clobber=True)    

    if latex:
        print " "
        for idx,gzline in enumerate(gzdata[:10]):
            print "%s & %8s & %3i & %3i & %3i & %5.1f & %5.3f & %5.3f & %5.3f & %i & %3i & %5.1f & %5.3f & %5.3f & %5.3f & %i \\\\" % \
                  (gzline['objid'], \
                  gzline['sample'], 
                  #gzline['RA'],gzline['DEC'],
                  gzline['t01_smooth_or_features_total_count'], gzline['total_count'],
                  gzline['t01_smooth_or_features_a01_smooth_count'],   gzline['t01_smooth_or_features_a01_smooth_weight'],
                  gzline['t01_smooth_or_features_a01_smooth_fraction'],gzline['t01_smooth_or_features_a01_smooth_weighted_fraction'],
                  adj_t01_a01[idx], flag_t01_a01[idx],
                  gzline['t01_smooth_or_features_a02_features_or_disk_count'],   gzline['t01_smooth_or_features_a02_features_or_disk_weight'],
                  gzline['t01_smooth_or_features_a02_features_or_disk_fraction'],gzline['t01_smooth_or_features_a02_features_or_disk_weighted_fraction'],
                  adj_t01_a02[idx], flag_t01_a02[idx])
        print " "

    if imagelist:
        gzdata_flagged = gzdata[flag_t01_a01.astype(bool)]
        print " "
        print "objid      ra         dec"
        for idx,line in enumerate(gzdata_flagged[:20]):
            print "%s & %6.3f & %6.3f" % (line['objid'],line['RA'],line['DEC'])
        print " "
    
    return None

def objlist(gzdata,str):

    ra = gzdata['RA'][gzdata[str+'_flag'].astype(bool)]
    dec = gzdata['DEC'][gzdata[str+'_flag'].astype(bool)]
    
    print ' '
    print 'ra dec'
    ll = len(ra) if len(ra) < 25 else 25
    for coo in sample(zip(ra,dec),ll):
        print "%6.3f, %6.3f" % (coo[0],coo[1])
    print ' '

    return None
