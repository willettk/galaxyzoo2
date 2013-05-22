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
import os
import matplotlib
import urllib

from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
from random import sample

from itertools import chain
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
paper_figures_path = gz_path+'datapaper/figures/'

# Locations of input FITS data files
gz2table_data_file = fits_path_main+'gz2table.fits'
gz2sample_data_file = fits_path_main+'gz2sample.fits'
gz2_photoz_data_file = fits_path_main+'gz2_photoz_table_sample.fits'
gz2_coadd1_data_file = fits_path_main+'gz2_coadd1_table_sample.fits'
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
SBlim_size = np.arange(200) / 10.0

appmag_lim_main = 17.0
appmag_lim_s82_normal = 17.7
appmag_lim_s82_coadd = appmag_lim_main + 2.0

# Collect memory and run plotting interactively
gc.collect()
plt.ion()

vt_stripe82_coadd = 10
vt_mainsample = 20

s82_dict = {'normal':{'file':gz2_stripe82_data_file,'str':'s82_normal_','maglim':appmag_lim_s82_normal,'vt':vt_mainsample},
            'coadd1':{'file':gz2_coadd1_data_file,  'str':'s82_coadd1_','maglim':appmag_lim_s82_coadd,'vt':vt_stripe82_coadd},
            'coadd2':{'file':gz2_coadd2_data_file,  'str':'s82_coadd2_','maglim':appmag_lim_s82_coadd,'vt':vt_stripe82_coadd}}

def bin_data_idl(task_dict, zmin = 0.00, zmax = 0.26, zstep = 0.01,
                 magmin = -24, magmax = -16, magstep = 0.25,
                 sizemin = 0, sizemax = 15, sizestep = 0.5,
                 stripe82 = False, depth = 'normal', 
                 whichmethod = 'new'):

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

    whichmethod : str
        'old' is for the threshold method of individual vote fractions (not systematically chosen)
        'new' uses the same threshold of 20 votes per task and selects vote fraction from
            previous dependent question (preferred)

    Returns
    -------
    None

    Notes
    -------
    Calls the IDL function `bin_gz2_from_python.pro`. Uses IDL since
        it was unclear how to bin in three dimensions while tracking
        indexes of metadata in pure Python. 

    """

    import subprocess

    tstart = time.time()

    if stripe82:
        file_str = s82_dict[depth]['str']
        datafile = s82_dict[depth]['file']
        vote_threshold = s82_dict[depth]['vt']
        mag_lim = s82_dict[depth]['maglim']
    else:
        datafile = gz2_maglim_data_file
        mag_lim = appmag_lim_main
        vote_threshold = vt_mainsample
        file_str = ''

    p = pyfits.open(datafile)
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
    absmag_lim = mag_lim - dmod
    r90_kpc = r90 * ang_scale

    absmag_padding = 0.0
    sb_padding = 0.0

    magnitude_mask = mr > (absmag_lim - absmag_padding)
    size_mask = r90_kpc < (3.0 * ang_scale)
    maglim_mask = (gzdata_hasredshift['PETROMAG_R'] - gzdata_hasredshift['EXTINCTION_R']) > mag_lim

    """
    Removed the surface brightness requirement - no selection algorithms for the GZ2 sample should depend
    on this, and I'm worried that the signs in the equation below are not correct. 

    surfacebrightness_mask = (mr + dmod + \
        2.5*np.log10(6.283185*(r50_kpc/ang_scale)**2)) > (SBlim_app - sb_padding)
    """


    # Old method -- 0.50 probability for each task, multiplied down the tree

    if whichmethod is 'old':
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

        taskcount_mask = gzdata_hasredshift[task_dict['task_name_count']] < task_dict['min_classifications']

        task_mask = taskcount_mask | taskprob_mask


    # New method: number of votes for current task + vote fraction for previous dependent task
    # from gz2_task_histograms.pro; developed with BDS, 25 Mar 2013
    
    if whichmethod is 'new':

        if task_dict['var_def'] not in ('task01','task06'):
            probarr_prevtask = gzdata_hasredshift[task_dict['dependent_tasks'][-1]] if isinstance(task_dict['dependent_tasks'],tuple) else gzdata_hasredshift[task_dict['dependent_tasks']]
            votearr_thistask = gzdata_hasredshift[task_dict['task_name_count']]

            taskprob_mask = (probarr_prevtask < task_dict['vf_prev'][str(vote_threshold)]) & (votearr_thistask < vote_threshold)
        else:
            votearr_thistask = gzdata_hasredshift[task_dict['task_name_count']]
            taskprob_mask = votearr_thistask < vote_threshold

        task_mask = taskprob_mask


    totalmask = magnitude_mask | size_mask | task_mask | maglim_mask

    print ' '
    print '%7i galaxies removed from sample due to absolute magnitude cutoff'%(np.sum(magnitude_mask.astype(int)))
    print '%7i galaxies removed from sample due to apparent magnitude cutoff'%(np.sum(maglim_mask.astype(int)))
    print '%7i galaxies removed from sample due to angular size cutoff'%(np.sum(size_mask.astype(int)))
    #print '%7i galaxies removed from sample due to surface brightness cutoff'%(np.sum(surfacebrightness_mask.astype(int)))
    #print '%7i galaxies removed from sample due to minimum number (%i) of classifications'%(np.sum(taskcount_mask.astype(int)),task_dict['min_classifications'])
    #print '%7i galaxies removed from sample due to task probability cutoff'%(np.sum(taskprob_mask.astype(int)))
    print '%7i galaxies removed from sample due to task probability and vote cutoffs'%(np.sum(task_mask.astype(int)))
    print '%6.2f percent of the total (%i galaxies) is kept for %s' % ((1. - np.sum(totalmask.astype(float))/len(gzdata_hasredshift)) * 100., np.sum(np.logical_not(totalmask).astype(int)),task_dict['var_def'])
    print ' '

    gzgood = gzdata_hasredshift[np.logical_not(totalmask)]

    # Write the masked GZ2 data to a FITS file to be binned by IDL script

    pyfits.writeto(fits_path_task+'%s_%sdata_for_idl.fits' % (task_dict['var_def'],file_str),
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
    f2.write('file_str = "%s" \n' % file_str)

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

def get_bins(task_dict,stripe82=False,depth='normal'):

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

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

    idl_binned = pyfits.open(fits_path_task+'%s_%sidlbinned.fits' % (task_dict['var_def'],file_str))
    bindata = idl_binned[1].data

    edges_redshift = np.squeeze(bindata['EDGES_REDSHIFT'])
    edges_mag = np.squeeze(bindata['EDGES_MAG'])
    edges_size = np.squeeze(bindata['EDGES_SIZE'])
    centers_redshift = np.squeeze(bindata['CENTERS_REDSHIFT'])
    centers_mag = np.squeeze(bindata['CENTERS_MAG'])
    centers_size = np.squeeze(bindata['CENTERS_SIZE'])

    return centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size 

def determine_ratio_baseline(task_dict, vartop = 0, varbot = 1, stripe82 = False, depth='normal'):

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

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    # Load the GZ2 sample data

    p_allvar = pyfits.open(fits_path_task+'%s_%sidlbinned.fits' % (var_def,file_str))
    d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
    d_var1 = np.squeeze(d_allvar[vartop,:,:,:])
    d_var2 = np.squeeze(d_allvar[varbot,:,:,:])

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict,stripe82,depth)

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

    pyfits.writeto(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline.fits' % (var_def,vartop,varbot,file_str),
                   ratio_baseline, clobber=True)    
    pyfits.writeto(fits_path_task+'%s_r%s%s_%slocal_counts_baseline.fits' % (var_def,vartop,varbot,file_str),
                   counts_baseline, clobber=True)    
    pyfits.writeto(fits_path_task+'%s_r%s%s_%slocal_redshift_baseline.fits' % (var_def,vartop,varbot,file_str),
                   redshift_baseline, clobber=True)    

    pickle.dump(ratio_baseline_masked, open(pkl_path+'%s_r%s%s_%slocal_ratio_baseline_masked.pkl' % (var_def,vartop,varbot,file_str),'wb')) 
    pickle.dump(counts_baseline_masked, open(pkl_path+'%s_r%s%s_%slocal_counts_baseline_masked.pkl' % (var_def,vartop,varbot,file_str),'wb')) 
    pickle.dump(redshift_baseline_masked, open(pkl_path+'%s_r%s%s_%slocal_redshift_baseline_masked.pkl' % (var_def,vartop,varbot,file_str),'wb')) 

    return ratio_baseline

def determine_ratio_baseline_sigma(task_dict, plot=False, vartop=0, varbot=1, stripe82=False, depth='normal'):


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

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

    var_def = task_dict['var_def']

    data = pyfits.getdata(fits_path_task+'%s_r%s%s_%slocal_counts_baseline.fits' % (var_def,vartop,varbot,file_str))
    var1 = data[:,:,0].astype(np.float)
    var2 = data[:,:,1].astype(np.float)
    mask_var1 = var1 < 1
    mask_var2 = var2 < 1
    mask_all = np.logical_and(mask_var1, mask_var2)
    mask = np.logical_or(mask_var1, mask_var2)
    ratio = pyfits.getdata(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline.fits' % (var_def,vartop,varbot,file_str))
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
        ax.set_xlabel(r'$M_r [mag]$',fontsize=22)
        ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_title('%s baseline ratio sigma' % var_def,fontsize=22)

    pyfits.writeto(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline_sigma.fits' % (var_def,vartop,varbot,file_str),
           sigma, clobber=True)    

    return None

def plot_ratio_baseline(task_dict,
                        unset_vrange = False,
                        ratio_vrange = (min_ratio,max_ratio),
                        plot_mag_lims=False,
                        vartop=0, varbot=1,
                        stripe82 = False, depth='normal'):

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

    if stripe82:
        file_str = s82_dict[depth]['str']
        mag_lim = s82_dict[depth]['maglim']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    f_ratio = fits_path_task+'%s_r%s%s_%slocal_ratio_baseline.fits' % (var_def,vartop,varbot,file_str)
    p_ratio = pyfits.open(f_ratio)

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict,stripe82,depth)

    ratio_baseline_masked = pickle.load(open(pkl_path+'%s_r%s%s_%slocal_ratio_baseline_masked.pkl' % (var_def,vartop,varbot,file_str),'rb'))
    counts_baseline_masked = pickle.load(open(pkl_path+'%s_r%s%s_%slocal_counts_baseline_masked.pkl' % (var_def,vartop,varbot,file_str),'rb'))
    sumcounts = np.sum(counts_baseline_masked, axis=2)

    if task_dict['ratio_type'] is 'linear':
        label_prefix = ''
    else:
        label_prefix = 'log_{10}'

    titlenames = ('counts','ratio')
    labelnames = (r'$N_{%s} + N_{%s}$' % (var_str[vartop],var_str[varbot]),r'$%s(N_{%s}/N_{%s})$' % (label_prefix,var_str[vartop],var_str[varbot]))

    if unset_vrange:
        ratio_vrange = (np.min(ratio_baseline_masked),np.max(ratio_baseline_masked))

    vranges = [(0,np.max(sumcounts)),ratio_vrange]

    fig = plt.figure(3,(14,8))
    fig.clf()

    for index, data in enumerate((sumcounts, ratio_baseline_masked)):
        ax = fig.add_subplot(1,2,index+1,aspect=1)
        cmap = cm.jet
        cmap.set_bad('k')
        imextent=(edges_mag[0],edges_mag[-1],edges_size[0],edges_size[-1])
        im = ax.imshow(data.T, 
                       vmin = vranges[index][0],
                       vmax = vranges[index][1],
                       extent=imextent,
                       interpolation='nearest',
                       origin='lower'
                       )
        ax.set_title('%s/%s ' % (var_str[vartop],var_str[varbot])+titlenames[index])
        ax.set_xlabel(r'$M_r [mag]$',fontsize=22)
        ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        cb = plt.colorbar(im,orientation='vertical')
        cb.set_label(labelnames[index],fontsize=16)

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
        absmag_lim = mag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
        absmag_lim_loz = mag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = mag_lim - cosmology.dmod_flat(0.25)
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

    fig.savefig(plots_path+'%s_r%s%s_%sratio_baseline.png' % (var_def,vartop,varbot,file_str), dpi=200)
    p_ratio.close()

    return ax

def plot_ratio_baseline_redshift(task_dict,vartop=0,varbot=1,stripe82 = False,depth='normal'):

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

    if stripe82:
        file_str = s82_dict[depth]['str']
        mag_lim = s82_dict[depth]['maglim']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    # Load the GZ2 sample data

    p_allvar = pyfits.open(fits_path_task+'%s_%sidlbinned.fits' % (var_def,file_str))
    d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
    centers_redshift, centers_mag, centers_size, edges_redshift, edges_mag, edges_size = get_bins(task_dict,stripe82,depth)

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
        absmag_lim = mag_lim - cosmology.dmod_flat(edges_redshift[idx])
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
    fig.text(0.5,0.05,r'$M_r [mag]$',ha='center',fontsize=20)
    fig.text(0.05,0.5,r'$R_{50} [kpc]$',va='center', rotation='vertical',fontsize=20)
    fig.text(0.5,0.94,
             '%s ratio per redshift bin for GZ2' % task_dict['var_def'], 
             fontsize=22, ha='center')

    fig.savefig(plots_path+'%s_%sratio_redshift.png' % (task_dict['var_def'],file_str), dpi=200)

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

    from scipy.stats import scoreatpercentile
    from scipy.optimize import fmin_powell
 
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

def fit_ratio_baseline(task_dict, nboot=50, plot=False, unset_vrange=False, vartop=0, varbot=1, stripe82 = False, depth='normal'):

    """ Fit the morphology ratio with an analytic function
        and estimate the uncertainty in each bin. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    nboot : int
        Number of tries to permute weights. Higher numbers give
        better accuracy on the estimated uncertainty. 

    plot : bool
        When `True`, calls `plot_ratio_function` and plots
        the data and overlaid best fit. 

    unset_vrange : bool
        When `True`, min and max of the colormap are set to the 
        range of the data.

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

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    funcname = task_dict['funcname']

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict,stripe82,depth)

    sigma = pyfits.getdata(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline_sigma.fits' % (var_def,vartop,varbot,file_str))
    ratio = pyfits.getdata(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline.fits' % (var_def,vartop,varbot,file_str))
    
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
    print 'Task %2s, responses %i and %i' % (var_def[-2:], vartop,varbot)
    print ' '
    print 'Number of bootstrap iterations: %i' % nboot
    print 'param = %s' % p
    print 'param errors = %s' % perr
    print 'param errors (percent) = %s' % (perr/p)
    print 'chi2r = %f' % chi2r
    print ' '

    ratio_baseline_fit = ratio_function(p, centers_mag, centers_size, funcname)
    pickle.dump((p, perr), file(pkl_path+'%s_r%s%s_%sratio_baseline_fit.pkl' % (var_def,vartop,varbot,file_str), 'w'))
    res = ratio - ratio_baseline_fit
    normres = res/sigma
    np.putmask(normres, weights < 0.000001, unknown_ratio)

    if plot:
        plot_ratio_function(p,centers_mag,centers_size, unset_vrange,funcname)

    pyfits.writeto(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline_fit.fits' % (var_def,vartop,varbot,file_str),
                   ratio_baseline_fit, clobber=True)    
    pyfits.writeto(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline_fitres.fits' % (var_def,vartop,varbot,file_str),
                   res, clobber=True) 
    pyfits.writeto(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline_fitnormres.fits' % (var_def,vartop,varbot,file_str),
                   normres, clobber=True) 

    return ratio_baseline_fit

def plot_ratio_baseline_fit(task_dict, 
                            fitbins=1000, 
                            unset_vrange = False, 
                            ratio_vrange = (min_ratio,max_ratio), 
                            match_databins = False, 
                            plot_local_transparent = True,
                            kernel_size = 3, 
                            plot_contour=False,
                            vartop=0,varbot=1,
                            stripe82=False, depth='normal'
                            ):

    """ Plot the best fit to morphology ratio using parameters
        computed in `fit_ratio_baseline`. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    fitbins : int
        Number of bins along x,y axes on which to plot function

    unset_vrange : bool
        When `True`, min and max of the colormap are set to the 
        range of the data.

    ratio_vrange: tuple
        Sets the min and max of the colormap. Only applies if
        `unset_vrange` = `False`.

    match_databins: bool
        When `True`, use the same bin size for the function grid
        as for the original data. 

    plot_local_transparent : bool
        When `True`, overplot the local baseline relation as
        a semi-transparent grid. 

    kernel_size : int
        Size of the kernel used to convolve the data and fit
        a boundary to the outer edge.

    plot_contour : bool
        When `True`, overplot a contour around the outer edge
        of the data as a dashed line. 

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

    from scipy.signal import convolve2d

    if stripe82:
        file_str = s82_dict[depth]['str']
        mag_lim = s82_dict[depth]['maglim']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    ratio_baseline_fit = pyfits.getdata(fits_path_task+'%s_r%s%s_%slocal_ratio_baseline_fit.fits' % (var_def,vartop,varbot,file_str))

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict,stripe82,depth)

    ratio_baseline_masked = pickle.load(open(pkl_path+'%s_r%s%s_%slocal_ratio_baseline_masked.pkl' % (var_def,vartop,varbot,file_str),'rb'))
    ratio_baseline_fit_masked = ma.array(ratio_baseline_fit,mask=ratio_baseline_masked.mask)

    kernel = np.ones((kernel_size,kernel_size),dtype=int)
    convarr = np.logical_not(ratio_baseline_masked.mask)
    bw = convolve2d(np.logical_not(ratio_baseline_masked.mask).astype(int),
                    kernel,mode='same').astype(bool).astype(int) * -1

    # Plot the best fit function as an opaque layer

    pfit, pfit_err = pickle.load(file(pkl_path+'%s_r%s%s_%sratio_baseline_fit.pkl' % (var_def,vartop,varbot,file_str), 'r'))
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
        ratio_vrange = (np.min(fitarray_fn),np.max(fitarray_fn))

    imf = ax.imshow(fitarray_fn.T,
                   alpha=1.0,
                   extent = fit_extent,
                   vmin = ratio_vrange[0], vmax = ratio_vrange[1], 
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
                       vmin = ratio_vrange[0], vmax = ratio_vrange[1], 
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
    absmag_lim = mag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
    size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
    ax.autoscale(False)
    ax.plot(SBlim_mag, SBlim_size,'w--')
    ax.axhline(size_1arcsec, color='w', linestyle='dashed')
    ax.axvline(absmag_lim,  color='w', linestyle='dashed')

    # Set general axes properties

    ax.set_xlabel(r'$M_r [mag]$',fontsize=22)
    ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
    ax.set_title('Baseline ratio data + fit',fontsize=22)
    ax.set_aspect('auto')
    rc(('xtick','ytick'), labelsize=12)

    #plt.show()
    #plt.draw()

    fig.savefig(plots_path+'%s_r%s%s_%sratio_baseline_fit.png' % (var_def,vartop,varbot,file_str), dpi=200)

    return None

def plot_ratio_function(p=np.array([-3.5, 1.9, -22.8, 0.3, 0.3, -4.4, 2.2, 1.1, 1.0]),
                        x = np.linspace(-24,-16, 50),
                        y = np.linspace(0,15,50),
                        unset_vrange = False,
                        funcname='sb',
                        ratio_vrange = (min_ratio,max_ratio)):

    """ Plot `ratio_function` on a 2D grid. 

    Parameters
    ----------
    p : numpy.array
        Parameters of `ratio_function` f[x,y;p]

    x : numpy.array
        x-axis of the data grid

    y : numpy.array
        y-axis of the data grid

    unset_vrange : bool
        When `True`, min and max of the colormap are set to the 
        range of the data.

    funcname: str
        Name of functional form to which data is fit. Options:

        * `sb` - original, 9-parameter equation from GZ1

        * `tilt - function varying linearly in both x and y 

    ratio_vrange: tuple
        Sets the min and max of the colormap. Only applies if
        `unset_vrange` = `False`.

    Returns
    -------
    None

    Notes
    -------

    # Effects of nine parameters in funcname=`sb`:

    # a - horizontal shift (along mag axis)
    # b - curvature tilt
    # c - vertical shift (along size axis)
    # d - gradient of the curve (smaller values = very sharp transition)
    # e - vertical flare at more negative magnitudes
    # f - multiplicative scaling of the cell values
    # g - additive scaling of the cell values
    # h - curvature strength (multiplicative)
    # i - curvature strength (exponential)

    """

    fit = ratio_function(p,x,y,funcname)

    if unset_vrange:
        ratio_vrange = (np.min(fit),np.max(fit))

    fig = plt.figure(4)
    fig.clf()
    ax = fig.add_subplot(111)
    im = ax.imshow(fit.T, 
        extent = (min(x),max(x),min(y),max(y)), 
        vmin = ratio_vrange[0], vmax = ratio_vrange[1],
        interpolation='nearest', 
        origin='lower')
    cb = plt.colorbar(im)
    ax.set_aspect('auto')
    ax.set_xlabel(r'$M_r [mag]$',fontsize=22)
    ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
    ax.set_title('Ratio function',fontsize=22)

    #plt.show()
    #plt.draw()

    fig.savefig(plots_path+'ratio_function.png', dpi=200)

    return None

def get_task_dict(task):

    """ Return dictionary of pre-computed parameters for
        reduction of each GZ2 task.

    Parameters
    ----------
    task : str
        Name of the GZ2 task. Possible options:

        * 'smooth'        Task 01; smooth, features/disk, or star/artifact
        * 'edgeon'        Task 02; edge-on disk or not
        * 'bar'           Task 03; bar or not
        * 'spiral'        Task 04; spiral structure or not
        * 'bulge'         Task 05; bulge prominence in disks
        * 'odd'           Task 06; anything odd?
        * 'rounded'       Task 07; ellipticity of smooth galaxies
        * 'odd_feature'   Task 08; ring, lens/arc, disturbed, irregular,
                                   other, merger, or dust lane
        * 'bulge_shape'   Task 09; rounded, boxy, or no bulge in edge-on disks
        * 'arms_winding'  Task 10; tightness (pitch angle) of spiral arms
        * 'arms_number'   Task 11; number of visible spiral arms


    Returns
    -------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. 

    Notes
    -------

    Parameters in dictionary:

      'task_name_count' : 
          name of the variable in gz2table corresponding to
          total number of weighted number of votes for that task
      
      'task_names_wf' : tuple
          Contains names of variables in gz2table corresponding to
          giving weighted vote fractions for each response

      'var_str' : tuple
          Contains string abbreviations for each reponse. 
          Used in filenames and labels in plots (legends,axes)

      'var_def' : str
          Abbreviation with task name and number (eg, 'task01'). 
          Used in filenames and labels in plots (legends,axes)
      
      'task_str' : str 
          Name of the task variable.
          Used in filenames and labels in plots (legends,axes)

      'ratio_type' : str
          Specifies form of the morphology ratio. 
          Options are either 'log' (default) or 'linear'.
      
      'min_prob' : float
          gives minimum vote fraction for clean 
          identification. Currently not in use, but could be set in
          `make_tables`. 

      'min_classifications' : int 
          minimum number of classifications per
          galaxy to be considered in establishing baseline ratio
          or included in clean sample. 
      
      'min_galperbin' : int
          minimum number of galaxies in each (M,R,z) bin to be
          included as "well-sampled" and fit with analytical function.

      'direct' : numpy.array(dtype=bool)
          Where `True`, correction for a given pair of responses is
          determined by directly taking the difference between the
          baseline ratio and the data. Where `False`, it takes the 
          difference between the analytically fit function and the
          data. 

      'funcname': str
          Name of functional form to which data is fit; to be passed to
          `ratio_function`. Options:

          * `sb` - original, 9-parameter equation from GZ1

          * `tilt - function varying linearly in both x and y 

      'dependent_tasks' : tuple
          Strings giving the variable names from gz2table for the
          weighted fractions for the responses that must be selected
          in order for this task to be viewed. For Task 02, for example,
          the user must have selected "features or disk" as the response
          to Task 01. Used to determine inclusion of galaxies in both
          baseline morphology ratio and clean samples. Not a key for
          `smooth` or `odd` tasks.  

      'vf_prev' : dict
          Vote fraction threshold for previous dependent task in tree.
          Not set for Tasks 01 or 06.
          Keys to dictionary are number of votes, for previous task, 
          as determined by find99.pro

    """

    assert task in ('smooth','edgeon','bar','spiral','bulge',\
        'rounded','odd','odd_feature','bulge_shape','arms_winding','arms_number'), \
        "Task must be a string in 'smooth','edgeon','bar','spiral','bulge','rounded','odd','odd_feature','bulge_shape','arms_winding','arms_number')"


    # Default dictionary keys

    task_dict = {'task_name_count' : 't01_smooth_or_features_total_weight',
                 'task_names_wf' : ('t01_smooth_or_features_a01_smooth_weighted_fraction',
                                    't01_smooth_or_features_a02_features_or_disk_weighted_fraction',
                                    't01_smooth_or_features_a03_star_or_artifact_weighted_fraction'),
                 'var_str' : ('el','sp','ar'),
                 'var_def' : 'task01',
                 'task_str' : task,
                 'ratio_type' : 'log',
                 'min_prob': 0.8,
                 'min_classifications': 30, 
                 'min_galperbin': 20, 
                 'direct' : np.ones((3,3),bool),
                 'funcname' : 'tilt'
                 }

    # Set parameters for various GZ2 tasks

    if task is 'smooth':
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['direct'][0,1] = False
        task_dict['funcname'] = 'sb'
 
    if task is 'edgeon':
        task_dict['task_name_count'] = 't02_edgeon_total_weight'
        task_dict['task_names_wf'] = ('t02_edgeon_a04_yes_weighted_fraction',
                                      't02_edgeon_a05_no_weighted_fraction')
        task_dict['var_str'] = ('edgeon','notedgeon')
        task_dict['var_def'] = 'task02'
        task_dict['min_classifications'] = 15
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['vf_prev'] = {'10':0.227,'20':0.430}
 
    if task is 'bar':
        task_dict['task_name_count'] = 't03_bar_total_weight'
        task_dict['task_names_wf'] = ('t03_bar_a06_bar_weighted_fraction',
                                      't03_bar_a07_no_bar_weighted_fraction')
        task_dict['var_str'] = ('bar','nobar')
        task_dict['var_def'] = 'task03'
        task_dict['min_classifications'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['vf_prev'] = {'10':0.519,'20':0.715}
 
    if task is 'spiral':
        task_dict['task_name_count'] = 't04_spiral_total_weight'
        task_dict['task_names_wf'] = ('t04_spiral_a08_spiral_weighted_fraction',
                                      't04_spiral_a09_no_spiral_weighted_fraction')
        task_dict['var_str'] = ('spiral','nospiral')
        task_dict['var_def'] = 'task04'
        task_dict['min_classifications'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['vf_prev'] = {'10':0.519,'20':0.715}
 
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
        task_dict['min_classifications'] = 10
        task_dict['ratio_type'] = 'log'
        task_dict['var_def'] = 'task05'
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['direct'][0,3]   = True
        task_dict['direct'][3,0]   = True
        task_dict['direct'][2,3]   = True
        task_dict['direct'][3,2]   = True
        task_dict['vf_prev'] = {'10':0.519,'20':0.715}

    if task is 'rounded':
        task_dict['task_name_count'] = 't07_rounded_total_weight'
        task_dict['task_names_wf'] = ('t07_rounded_a16_completely_round_weighted_fraction',
                                      't07_rounded_a17_in_between_weighted_fraction',
                                      't07_rounded_a18_cigar_shaped_weighted_fraction')
        task_dict['var_str'] = ('round','inbetween','cigar')
        task_dict['ratio_type'] = 'log'
        task_dict['var_def'] = 'task07'
        task_dict['min_classifications'] = 15
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a01_smooth_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['funcname']='sb'
        task_dict['direct'][0,1] = False
        task_dict['direct'][0,2] = False
        task_dict['direct'][1,2] = False
        task_dict['vf_prev'] = {'10':0.263,'20':0.469}

    if task is 'arms_winding':
        task_dict['task_name_count'] = 't10_arms_winding_total_weight'
        task_dict['task_names_wf'] = ('t10_arms_winding_a28_tight_weighted_fraction',
                                      't10_arms_winding_a29_medium_weighted_fraction',
                                      't10_arms_winding_a30_loose_weighted_fraction')
        task_dict['var_str'] = ('tight','medium','loose')
        task_dict['var_def'] = 'task10'
        task_dict['min_classifications'] = 5
        task_dict['ratio_type'] = 'log'
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a05_no_weighted_fraction','t04_spiral_a08_spiral_weighted_fraction')
        task_dict['direct']   = np.zeros((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['vf_prev'] = {'10':0.402,'20':0.619}

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
        task_dict['ratio_type'] = 'log'
        task_dict['min_classifications'] = 5
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
        task_dict['vf_prev'] = {'10':0.402,'20':0.619}

    if task is 'bulge_shape':
        task_dict['task_name_count'] = 't09_bulge_shape_total_weight'
        task_dict['task_names_wf'] = ('t09_bulge_shape_a25_rounded_weighted_fraction',
                                      't09_bulge_shape_a26_boxy_weighted_fraction',
                                      't09_bulge_shape_a27_no_bulge_weighted_fraction')
        task_dict['var_str'] = ('rounded','boxy','nobulge')
        task_dict['var_def'] = 'task09'
        task_dict['min_prob'] = 0.5
        task_dict['min_classifications'] = 10
        task_dict['dependent_tasks'] = ('t01_smooth_or_features_a02_features_or_disk_weighted_fraction','t02_edgeon_a04_yes_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['direct'][0,2] = False
        task_dict['direct'][2,0] = False
        task_dict['vf_prev'] = {'10':0.326,'20':0.602}

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
        task_dict['dependent_tasks'] = ('t06_odd_a14_yes_weighted_fraction')
        task_dict['direct']   = np.ones((len(task_dict['var_str']),len(task_dict['var_str'])),bool)
        task_dict['vf_prev'] = {'10':0.223,'20':0.420}

    return task_dict

def run_task(task_dict,
             nboot = 1, 
             fitbins = 1000,
             plot = True, 
             unset_vrange = True,
             stripe82 = False,
             depth='normal'
             ):
    
    """ Run all the steps in reduction pipeline
        for a single task. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    nboot : int
        Number of tries to permute weights. Higher numbers give
        better accuracy on the estimated uncertainty. 

    fitbins : int
        Number of bins along x,y axes on which to plot function

    plot : bool
        When `True`, create and save all plots in the reduction
        pipeline. Can be skipped to save time. 

    unset_vrange : bool
        When `True`, min and max of the colormap are set to the 
        range of the data in all plots. Only affects function
        if `plot` is `True`. 

    Returns
    -------
    None

    Raises
    -------
    UserWarning,RuntimeWarning

    Notes
    -------


    """

    import warnings

    tstart = time.time()

    warnings.filterwarnings('ignore', category=UserWarning,    append=True)
    warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

    assert type(task_dict) is dict, \
        "First argument in run_task must be a dictionary -- use galaxyzoo2.get_task_dict()"

    # Execute all the tasks to bin data, find baseline morphology, fit the baseline, 
    # and adjust the raw vote fractions

    bin_data_idl(task_dict,stripe82=stripe82,depth=depth)

    ntask = len(task_dict['task_names_wf'])

    # Determine baseline ratio and correction for all pairs of variables per task

    for idx1, vartop in enumerate(np.arange(ntask)):
        setdiff = np.concatenate((np.arange(vartop),np.arange(ntask - (vartop+1)) + (vartop+1) ))
        for idx2, varbot in enumerate(setdiff):

            rb = determine_ratio_baseline(task_dict,vartop,varbot,stripe82,depth)
            determine_ratio_baseline_sigma(task_dict,plot,vartop,varbot,stripe82,depth)
            # 10 bins assumed as a low minimum to achieve a fit
            if np.sum(rb > unknown_ratio) > 10 and task_dict['direct'][vartop,varbot] == False:
                fit_ratio_baseline(task_dict,nboot,plot,unset_vrange,vartop,varbot,stripe82,depth)
                if plot:
                    plot_ratio_baseline_fit(task_dict,fitbins,unset_vrange,vartop=vartop,varbot=varbot,stripe82=stripe82,depth=depth)
            elif task_dict['direct'][vartop,varbot] == False:
                    "No non-blank cells in baseline ratio for %s, responses %i and %i - setting task_dict['direct'] = True" % (task_dict['task_str'],vartop,varbot)
                    task_dict['direct'][vartop,varbot] = True
            determine_baseline_correction(task_dict,vartop,varbot,stripe82,depth)

            if plot:
                plot_ratio_baseline(task_dict,unset_vrange,vartop=vartop,varbot=varbot,stripe82=stripe82,depth=depth)
                plot_ratio_baseline_redshift(task_dict,vartop,varbot,stripe82=stripe82,depth=depth)
                plot_baseline_correction(task_dict,vartop=vartop,varbot=varbot,stripe82=stripe82,depth=depth)

    # Remaining tasks that do not require looping over variable pairs

    adjust_probabilities(task_dict,stripe82=stripe82,depth=depth)
    if plot:
        plot_galaxy_counts(task_dict,stripe82=stripe82,depth=depth)
        plot_type_fractions(task_dict,stripe82=stripe82,depth=depth)

    warnings.resetwarnings()

    tend = time.time()
    print 'Time elapsed for run_task on %s: %i seconds' % (task_dict['task_str'], tend - tstart)

    return None

def run_all_tasks(nboot = 1, 
                  fitbins = 1000,
                  plot = True, 
                  unset_vrange = True,
                  stripe82 = False,
                  depth='normal'
                  ):
    """ Runs full data pipeline for all tasks.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -------


    """


    tstart = time.time()

    tasklist = ['smooth','edgeon','bar','spiral','odd',
                'bulge','rounded','arms_winding','arms_number',
                'bulge_shape','odd_feature']

    for task in tasklist:
        td = get_task_dict(task)
        run_task(td,nboot,fitbins,plot,unset_vrange,stripe82,depth)

    tend = time.time()
    print 'Time elapsed for run_all_tasks: %i seconds' % (tend - tstart)

    return None

def determine_baseline_correction(task_dict, vartop = 0, varbot = 1, stripe82 = False, depth='normal'):

    """ Determine the correction ratio baseline for a single pair of
        responses to a GZ2 task. 

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
    correction : numpy.array
        3D array binned in (M,R,z); values are the correction
        C[vartop,varbot] determined from the difference between the
        baseline ratio and the data. 

    correction_masked : numpy.masked_array
        3D array binned in (M,R,z). Same values as `correction`, 
        except that bins without either well-sampled data or 
        a baseline measurement are masked. 

    ratio_masked : numpy.masked_array
        3D array binned in (M,R,z). Contains the morphology
        ratios for responses in each bin. Bins without well-sampled
        data or baseline measurements are masked. 

    """

    # Load best fit data

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

    var_def  = task_dict['var_def']
    var_str = task_dict['var_str']

    p_allvar = pyfits.open(fits_path_task+'%s_%sidlbinned.fits' % (var_def,file_str))
    d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
    d_var1 = np.squeeze(d_allvar[vartop,:,:,:])
    d_var2 = np.squeeze(d_allvar[varbot,:,:,:])

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict,stripe82,depth)

    p_allvar.close()

    ratio_baseline_masked = pickle.load(file(pkl_path+'%s_r%s%s_%slocal_ratio_baseline_masked.pkl' % (var_def,vartop,varbot,file_str),'rb')) 

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

    fit_filename = pkl_path+'%s_r%s%s_%sratio_baseline_fit.pkl' % (var_def,vartop,varbot,file_str)

    try:
        with open(fit_filename):
            if not task_dict['direct'][vartop,varbot]:
                pfit_params, pfit_params_err = pickle.load(file(pkl_path+'%s_r%s%s_%sratio_baseline_fit.pkl' % (var_def,vartop,varbot,file_str), 'r'))
                pfit = ratio_function(pfit_params, centers_mag, centers_size,task_dict['funcname'])
                correction = ratio - pfit
            else:
                correction = ratio - ratio_baseline_masked
    except IOError:
        correction = ratio - ratio_baseline_masked

    np.putmask(correction, mask, unknown_ratio)
    correction_masked = ma.masked_less_equal(correction, unknown_ratio)
    ratio_masked = ma.array(ratio, mask=correction_masked.mask)

    return correction, correction_masked, ratio_masked

def plot_baseline_correction(task_dict,ratio_vrange=(min_ratio,max_ratio),vartop=0,varbot=1,stripe82 = False,depth='normal'):

    """ Plot the derived correction for a single pair of responses to a
        GZ2 task. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    ratio_vrange: tuple
        Sets the min and max of the colormap. Only applies if
        `unset_vrange` = `False`.

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

    c,cmasked, ratio_masked = determine_baseline_correction(task_dict, vartop, varbot, stripe82, depth)

    if stripe82:
        file_str = s82_dict[depth]['str']
        mag_lim = s82_dict[depth]['maglim']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    p_allvar = pyfits.open(fits_path_task+'%s_%sidlbinned.fits' % (var_def,file_str))
    d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
    d_var1 = np.squeeze(d_allvar[vartop,:,:,:])
    d_var2 = np.squeeze(d_allvar[varbot,:,:,:])

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

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
                       vmin=ratio_vrange[0],vmax=ratio_vrange[1],
                       interpolation='nearest',origin='lower')
        ax.set_aspect('auto')

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(z)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(z))**2))
        absmag_lim = mag_lim - cosmology.dmod_flat(z)
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

    fig.savefig(plots_path+'%s_r%s%s_%sbaseline_correction.png' % (var_def,vartop,varbot,file_str), dpi=200)

    return None

def plot_baseline_correction_slice(task_dict, bslice=5,
                                   ratio_vrange = (min_ratio,max_ratio),
                                   vartop = 0, varbot=1, 
                                   smoothfunc = False,
                                   savefig=False,
                                   stripe82 = False, depth='normal'):

    """ Plot the data, fit, and correction for a single pair of responses to a
        GZ2 task at a particular redshift.  

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    bslice : int
        Index of the redshift slice in binned data cube to plot. 

    ratio_vrange: tuple
        Sets the min and max of the colormap. Only applies if
        `unset_vrange` = `False`.

    vartop: int
        Index for the variable in the numerator of the baseline
        morphology ratio. 

    varbot: int
        Index for the variable in the denominator of the baseline
        morphology ratio. 

    smoothfunc : bool
        When `True`, plot the analytic fit on a 1000x1000 grid. 

    savefig : bool
        When `True`, save the figure as a file. 

    Returns
    -------
    None

    """

    c,cmasked, ratio_masked = determine_baseline_correction(task_dict,vartop,varbot,stripe82,depth)

    if stripe82:
        file_str = s82_dict[depth]['str']
        mag_lim = s82_dict[depth]['maglim']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

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
                     vmin = ratio_vrange[0], vmax = ratio_vrange[1],
                     interpolation='nearest',origin='lower')
    ax1.set_title('Data ratio, %5.3f < z < %5.3f' % (centers_redshift[bslice]-zstep/2.,centers_redshift[bslice]+zstep/2.))
    ax1.set_aspect('auto')
    cb = fig.colorbar(im1)
    cb.set_label('ratio')

    SBlim_mag = (SBlim_app - cosmology.dmod_flat(0.07)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(0.07))**2))
    absmag_lim = mag_lim - cosmology.dmod_flat(0.07)
    size_1arcsec = cosmology.ang_scale_flat(0.07)
    ax1.autoscale(False)
    ax1.plot(SBlim_mag, SBlim_size,'w--')
    ax1.axhline(size_1arcsec, color='w', linestyle='dashed')
    ax1.axvline(absmag_lim,  color='w', linestyle='dashed')
    ax1.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
    ax1.set_xlabel(r'$M_r [mag]$',fontsize=22)
    plt.show()
    plt.draw()

    # Function

    pfit, pfit_err = pickle.load(file(pkl_path+'%s_r%s%s_%sratio_baseline_fit.pkl' % (var_def,vartop,varbot,file_str), 'r'))
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
                     vmin = ratio_vrange[0], vmax = ratio_vrange[1],
                     interpolation='nearest', 
                     origin='lower')
    cb = plt.colorbar(im2)
    ax2.set_title('Baseline fit')
    ax2.set_xlabel(r'$M_r [mag]$',fontsize=22)

    zslice_opaque = ma.copy(zslice)
    zslice_opaque.mask = np.logical_not(zslice_opaque.mask)
    zslice_opaque[np.logical_not(zslice_opaque.mask)] = 255.
    cmap_gray = cm.gray
    cmap_gray.set_bad(alpha=0.0)
    im2a = ax2.imshow(zslice_opaque.T, 
                      alpha=0.5,
                      extent = fit_extent,
                      vmin = ratio_vrange[0], vmax = ratio_vrange[1],
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
                     vmin = ratio_vrange[0], vmax = ratio_vrange[1],
                     interpolation='nearest',origin='lower')
    ax3.set_title('Correction')
    cb = fig.colorbar(im3)
    cb.set_label('baseline correction')

    SBlim_mag = (SBlim_app - cosmology.dmod_flat(0.07)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(0.07))**2))
    absmag_lim = mag_lim - cosmology.dmod_flat(0.07)
    size_1arcsec = cosmology.ang_scale_flat(0.07)
    ax3.plot(SBlim_mag, SBlim_size,'w--')
    ax3.axhline(size_1arcsec, color='w', linestyle='dashed')
    ax3.axvline(absmag_lim,  color='w', linestyle='dashed')
    ax3.set_xlabel(r'$M_r [mag]$',fontsize=22)
    ax3.set_xlim(fit_extent[:2])
    ax3.set_ylim(fit_extent[2:])
    ax3.autoscale(False)
    ax3.set_aspect('auto')

    if savefig:
        fig.savefig(plots_path+'%s_r%s%s_%sbaseline_correction_slice%03i.png' % (task_dict['var_def'],vartop,varbot,file_str,bslice), dpi=200)

    return None

def ratio_adj(p_el,p_sp,correction,ratio_type):
    
    """ Function to determine the adjusted elliptical/spiral
        vote ratio, based on raw data and a correction.

    Parameters
    ----------
    p_el : float
        vote fraction for elliptical galaxies. Assumed to be parameter
        that decreases when positive correction is applied. 

    p_sp : float
        vote fraction for spiral galaxies. Assumed to be parameter
        that increases when positive correction is applied. 

    correction : float
        Correction to apply to individual vote fractions based on
        difference in morphology ratio for this galaxy as function
        of (M,R,z) and a baseline relation. 

    ratio_type' : str
        Specifies form of the morphology ratio. 
        Options are either 'log' (default) or 'linear'.
      
    Returns
    -------
    ratio_adj : float
        elliptical/spiral ratio adjusted for classification bias

    Notes
    -------


    """

    assert ratio_type in ('linear','log'), \
        "Ratio type must be defined and either 'linear' or 'log'"

    if ratio_type is 'log':
        ratio_adj = (p_el / p_sp) / (10.**(correction))

    if ratio_type is 'linear':
        ratio_adj = (p_el / p_sp) / correction

    return ratio_adj

def p_x(p_el,p_sp):

    """ Gives the percentage of any votes not allocated to
        one of the two primary responses (el or sp).

    Parameters
    ----------
    p_el : float
        vote fraction for elliptical galaxies. Assumed to be parameter
        that decreases when positive correction is applied. 

    p_sp : float
        vote fraction for spiral galaxies. Assumed to be parameter
        that increases when positive correction is applied. 

    Returns
    -------
    p_x : float
        Fraction of votes not allocated to either of the primary responses.

    Notes
    -------


    """

    p_x = 1. - p_el - p_sp

    return p_x

def p_el_adj(p_el,p_sp,correction,ratio_type):

    """ Function to determine the adjusted elliptical vote fraction

    Parameters
    ----------
    p_el : float
        vote fraction for elliptical galaxies. Assumed to be parameter
        that decreases when positive correction is applied. 

    p_sp : float
        vote fraction for spiral galaxies. Assumed to be parameter
        that increases when positive correction is applied. 

    correction : float
        Correction to apply to individual vote fractions based on
        difference in morphology ratio for this galaxy as function
        of (M,R,z) and a baseline relation. 

    ratio_type' : str
        Specifies form of the morphology ratio. 
        Options are either 'log' (default) or 'linear'.
      
    Returns
    -------
    p_el_adj : float
        elliptical vote fraction adjusted for classification bias

    Notes
    -------


    """

    denom_el = 1./ratio_adj(p_el,p_sp,correction,ratio_type) + (p_x(p_el,p_sp) / p_el) + 1.
    p_el_adj = 1./denom_el

    return p_el_adj

def p_sp_adj(p_el,p_sp,correction,ratio_type):

    """ Function to determine the adjusted spiral vote fraction

    Parameters
    ----------
    p_el : float
        vote fraction for elliptical galaxies. Assumed to be parameter
        that decreases when positive correction is applied. 

    p_sp : float
        vote fraction for spiral galaxies. Assumed to be parameter
        that increases when positive correction is applied. 

    correction : float
        Correction to apply to individual vote fractions based on
        difference in morphology ratio for this galaxy as function
        of (M,R,z) and a baseline relation. 

    ratio_type' : str
        Specifies form of the morphology ratio. 
        Options are either 'log' (default) or 'linear'.
      
    Returns
    -------
    p_sp_adj : float
        spiral vote fraction adjusted for classification bias

    Notes
    -------


    """

    denom_sp = ratio_adj(p_el,p_sp,correction,ratio_type) + (p_x(p_el,p_sp) / p_sp) + 1.
    p_sp_adj = 1./denom_sp

    return p_sp_adj

def p_adj_new(f, corr):

    """ Function to determine the adjusted vote fraction
        for a given response, given all raw vote fractions
        and corrections between each pair. 

    Parameters
    ----------
    f : numpy.array
        vote fractions for all responses to the task. The first 
        entry is the vote fraction being corrected for.

    corr : numpy.array
        Corrections for each pair of responses in the task with respect
        to f[0]. The task order must be the same as that in f (once the
        response being corrected for has been removed). 

    Returns
    -------
    p_adj: float
        adjusted vote fraction

    Notes
    -------

    Example: there are three responses to a task (i,j,k). To determine
        the adjusted vote fraction for response j:

        >> f = numpy.array((f_j,f_i,f_k))
        >> corr = numpy.array((C_ji,C_jk))
        >> p_j = p_adj_new(f,corr)
        
    This would be improved if a dictionary were used instead of arrays.

    """

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

def adjust_probabilities(task_dict, stripe82=False, depth='normal', photoz=False):

    """ Adjust vote fractions for all responses to a task
        based on the derived corrections. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    stripe82 : bool
        When `True`, adjust probabilities for spectroscopic Stripe 82 data.

    photoz : bool
        When `True`, adjust probabilities for main sample data with
        photometric redshifts.

    Returns
    -------
    p_raw : numpy.array
        (N x M) array of raw vote fractions for task, where N is the number
        of galaxies in the sample and M the number of responses for the task.

    p_adj : numpy.array
        (N x M) array of adjusted vote fractions for task, where N is the number
        of galaxies in the sample and M the number of responses for the task.

    Notes
    -------
    Default settings run the analysis on the spectroscopic main-sample data. 

    """

    timestart2 = time.time()

    # Load in all raw probabilities

    if stripe82:
        file_str = s82_dict[depth]['str']
        p = pyfits.open(s82_dict[depth]['file'])
    elif photoz:
        p = pyfits.open(gz2_photoz_data_file)
        file_str = 'photoz_'
    else:
        p = pyfits.open(gz2_maglim_data_file)
        file_str = ''

    gzdata_all = p[1].data
    p.close()

    # Restrict to galaxies with redshifts

    if photoz:
        gzdata = gzdata_all
        redshift = gzdata['photoz']
        mr = (gzdata['PETROMAG_R'] - gzdata['EXTINCTION_R']) - cosmology.dmod_flat(redshift)
        r50_kpc = gzdata['PETROR50_R'] * cosmology.ang_scale_flat(redshift)
    else:
        gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]
        redshift = gzdata['REDSHIFT']
        mr = gzdata['PETROMAG_MR']
        r50_kpc = gzdata['PETROR50_R_KPC']

    ntask = len(task_dict['task_names_wf'])
    ngals = len(redshift)

    # Load in the bin sizes

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

    zstep = edges_redshift[1] - edges_redshift[0]
    magstep = edges_mag[1] - edges_mag[0]
    sizestep = edges_size[1] - edges_size[0]

    corr_test, corrmasked_test, ratiomasked_test = determine_baseline_correction(task_dict, 0, 1, stripe82, depth)
    cshape = corr_test.shape

    allcorr        = np.zeros((ntask,ntask-1,cshape[0],cshape[1],cshape[2]),dtype=float)
    allcorrmasked  = np.zeros_like(allcorr)
    allratiomasked = np.zeros_like(allcorr)

    # Loop over combinations of tasks and retrieve corrections

    for idx1, var1 in enumerate(np.arange(ntask)):
        setdiff = np.concatenate((np.arange(var1),np.arange(ntask - (var1+1)) + (var1+1) ))
        for idx2, var2 in enumerate(setdiff):
            corr, corrmasked, ratiomasked = determine_baseline_correction(task_dict, var1, var2, stripe82, depth)
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
          
        else: # if galaxy does not appear in the correction area. Should only happen if the galaxy doesn't have a redshift.
            outsidebinscount += 1
            p_adj[:,idx] = f

            '''
            corrarr[:,idx] = unknown_ratio
            corrbinarr[:,:,idx] = unknown_ratio
            '''

    pickle.dump(corrarr, open(pkl_path+'%s_%scorrarr.pkl' % (task_dict['var_def'], file_str),'wb')) 
    pickle.dump(corrbinarr, open(pkl_path+'%s_%scorrbinarr.pkl' % (task_dict['var_def'], file_str),'wb')) 
    pickle.dump(corrgoodarr, open(pkl_path+'%s_%scorrgoodarr.pkl' % (task_dict['var_def'], file_str),'wb')) 

    timeend2 = time.time()
    print ' '
    print '%7i galaxies had correction bins outside the binned volume' % outsidebinscount
    print '%7i galaxies had NaN corrections for all pairs of variables' % unknowncorrcount
    print '%7i galaxies had no non-zero corrections in any bin' % zerocorrcount
    print '%7i galaxies had finite corrections available' % corrcount
    
    if corrcount > 0:
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

    pickle.dump(p_raw, open(pkl_path+'%s_%sraw_probabilities.pkl' % (task_dict['var_def'], file_str),'wb')) 
    pickle.dump(p_adj, open(pkl_path+'%s_%sadj_probabilities.pkl' % (task_dict['var_def'], file_str),'wb')) 

    return p_raw,p_adj

def plot_galaxy_counts(task_dict,vmin=0,vmax=1000,stripe82=False,depth='normal'):

    """ Plot 2D histogram of galaxy counts for binned GZ2 data
        for each task. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    vmin : int or float
        minimum range of the colormap

    vmax : int or float
        maximum range of the colormap

    Returns
    -------
    None

    Notes
    -------


    """

    if stripe82:
        p = pyfits.open(s82_dict[depth]['file'])
        mag_lim = s82_dict[depth]['maglim']
        vote_threshold = s82_dict[depth]['vt']
        file_str = s82_dict[depth]['str']
    else:
        p = pyfits.open(gz2_maglim_data_file)
        mag_lim = appmag_lim_main
        vote_threshold = vt_mainsample
        file_str = ''

    gzall = p[1].data
    p.close()

    gzdata = gzall[(gzall[task_dict['task_name_count']] > vote_threshold) & np.isfinite(gzall['REDSHIFT'])]

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)
    zstep = centers_redshift[1] - centers_redshift[0]
    magmin,magmax = edges_mag.min(),edges_mag.max()
    sizemin,sizemax = edges_size.min(),edges_size.max()
    magstep = centers_mag[1] - centers_mag[0]
    sizestep = centers_size[1] - centers_size[0]

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
                           vmin=vmin,vmax=vmax,
                           interpolation='nearest',origin='lower')
        #ax.set_title('%s < z < %s' % (z-zstep/2.,z+zstep/2.))
        ax.set_aspect('auto')

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(z)- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(z))**2))
        absmag_lim = mag_lim - cosmology.dmod_flat(z)
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

    fig.savefig(plots_path+'%s_%sgalaxy_counts.png' % (task_dict['var_def'],file_str), dpi=200)

    return fig

def plot_type_fractions(task_dict, zlo = 0.01, zhi=0.085, zwidth=0.02, stripe82=False, depth='normal', ploterrors=False, photoz=False, plotadj=True):

    """ Plot the mean type fractions for all responses for each task
        as a function of redshift. Left plot shows all galaxies,
        right plot is a magnitude-limited sample. 

    Parameters
    ----------
    task_dict : dict
        Dictionary specifying parameters for the reduction of each
        GZ2 task. Called from `get_task_dict`.

    zlo : float
        minimum redshift range to plot

    zhi : float
        redshift at which to set the magnitude limit for the sample,
        based on the apparent magnitude sensitivity at that redshift

    zwidth : float
        width of the redshift bins

    stripe82 : bool
        When `True`, plot the results for the Stripe 82 spectroscopic sample

    ploterrors : bool
        When `True`, plot error bars on each type fraction


    Returns
    -------
    None

    Notes
    -------


    """

    # Load data 

    if stripe82:
        file_str = s82_dict[depth]['str']
        p = pyfits.open(s82_dict[depth]['file'])
        vote_threshold = s82_dict[depth]['vt']
        mag_lim = s82_dict[depth]['maglim']
    elif photoz:
        p = pyfits.open(gz2_photoz_data_file)
        file_str = 'photoz_'
        vote_threshold = vt_mainsample
        mag_lim = appmag_lim_main
    else:
        p = pyfits.open(gz2_maglim_data_file)
        file_str = ''
        vote_threshold = vt_mainsample
        mag_lim = appmag_lim_main

    gzdata_all = p[1].data
    p.close()

    if photoz:
        gzdata = gzdata_all
        redshift = gzdata['photoz']
        mr = (gzdata['PETROMAG_R'] - gzdata['EXTINCTION_R']) - cosmology.dmod_flat(redshift)
        r50_kpc = gzdata['PETROR50_R'] * cosmology.ang_scale_flat(redshift)
    else:
        gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]
        redshift = gzdata['REDSHIFT']
        mr = gzdata['PETROMAG_MR']
        r50_kpc = gzdata['PETROR50_R_KPC']

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

    task_counts = gzdata[task_dict['task_name_count']]

    zplotbins = np.arange(min(edges_redshift[0],zlo),edges_redshift[-1],zwidth)

    ntask = len(task_dict['var_str'])

    corrgoodarr = pickle.load(open(pkl_path+'%s_%scorrgoodarr.pkl' % (task_dict['var_def'],file_str),'rb')) 
    p_raw = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % (task_dict['var_def'],file_str),'rb')) 
    p_adj = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % (task_dict['var_def'],file_str),'rb')) 

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

    maglimval = mag_lim - cosmology.dmod_flat(zhi)

    # Create the two samples
    
    """
    # Old method
    n_all = np.sum((task_counts > task_dict['min_classifications']) & corrgoodarr)
    n_maglim = np.sum((mr < maglimval) & (task_counts > task_dict['min_classifications']) & corrgoodarr)
    """

    # Begin new method

    if task_dict['var_def'] not in ('task01','task06'):
        probarr_prevtask = gzdata[task_dict['dependent_tasks'][-1]] if isinstance(task_dict['dependent_tasks'],tuple) else gzdata[task_dict['dependent_tasks']]
        votearr_thistask = gzdata[task_dict['task_name_count']]
        taskprob_good = (probarr_prevtask >= task_dict['vf_prev'][str(vote_threshold)]) & (votearr_thistask >= vote_threshold)
    else:
        votearr_thistask = gzdata[task_dict['task_name_count']]
        taskprob_good = votearr_thistask >= vote_threshold

    """
    n_all = np.sum(taskprob_good & corrgoodarr)
    n_maglim = np.sum((mr < maglimval) & taskprob_good & corrgoodarr)
    """
    n_all = np.sum(taskprob_good)
    n_maglim = np.sum((mr < maglimval) & taskprob_good)

    # End new method

    # Loop over redshift bin and find the type fractions at each slice

    for idx,zbin in enumerate(zplotbins):

        gals_in_bin = (redshift >= zbin) & \
                      (redshift < (zbin+zwidth)) & \
                      taskprob_good
                      #corrgoodarr & \
                      #(task_counts > task_dict['min_classifications']) & \         Old method

        p_raw_typefrac[:,idx] = np.mean(p_raw[:,gals_in_bin],axis=1)
        p_adj_typefrac[:,idx] = np.mean(p_adj[:,gals_in_bin],axis=1)
        p_raw_typefrac_err[:,idx] = np.std(p_raw[:,gals_in_bin],axis=1)
        p_adj_typefrac_err[:,idx] = np.std(p_adj[:,gals_in_bin],axis=1)

        gals_in_bin_maglim = (redshift >= zbin) & \
                      (redshift < (zbin+zwidth)) & \
                      (mr < maglimval) & \
                      taskprob_good
                      #corrgoodarr & \
                      #(task_counts > task_dict['min_classifications']) & \         Old method

        p_raw_typefrac_maglim[:,idx] = np.mean(p_raw[:,gals_in_bin_maglim],axis=1)
        p_adj_typefrac_maglim[:,idx] = np.mean(p_adj[:,gals_in_bin_maglim],axis=1)
        p_raw_typefrac_err_maglim[:,idx] = np.std(p_raw[:,gals_in_bin_maglim],axis=1)
        p_adj_typefrac_err_maglim[:,idx] = np.std(p_adj[:,gals_in_bin_maglim],axis=1)

    pickle.dump(p_raw_typefrac,        open(pkl_path+'%s_%sraw_typefrac.pkl' % (task_dict['var_def'],file_str),'wb')) 
    pickle.dump(p_adj_typefrac,        open(pkl_path+'%s_%sadj_typefrac.pkl' % (task_dict['var_def'],file_str),'wb')) 
    pickle.dump(p_raw_typefrac_maglim, open(pkl_path+'%s_%sraw_typefrac_maglim.pkl' % (task_dict['var_def'],file_str),'wb')) 
    pickle.dump(p_adj_typefrac_maglim, open(pkl_path+'%s_%sadj_typefrac_maglim.pkl' % (task_dict['var_def'],file_str),'wb')) 

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

    if plotadj:
        colorarr = ['r','b','m','g','c','y','k'][::-1]
        for idx,task_str in enumerate(task_dict['var_str']):
            linecolor = colorarr.pop()
            ax1.plot(zplotbins, p_adj_typefrac[idx,:], color=linecolor, linestyle='--' ,linewidth=4)
            legend_adj_str.append('adj %s' % task_str)

    if ploterrors:
        ax1.errorbar(zplotbins, p_raw_typefrac[idx,:], yerr = p_raw_typefrac_err[idx,:], 
           color=linecolor, linestyle='-' ,linewidth=2,elinewidth=1)
        if plotadj:
            ax1.errorbar(zplotbins, p_adj_typefrac[idx,:], yerr = p_adj_typefrac_err[idx,:], 
               color=linecolor, linestyle='-' ,linewidth=4,elinewidth=1)


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

    if plotadj:
        colorarr = ['r','b','m','g','c','y','k'][::-1]
        for idx,task_str in enumerate(task_dict['var_str']):
            linecolor = colorarr.pop()
            ax2.plot(zplotbins, p_adj_typefrac_maglim[idx,:], color=linecolor, linestyle='--' ,linewidth=4)

    if ploterrors:
        ax2.errorbar(zplotbins, p_raw_typefrac_maglim[idx,:], yerr = p_raw_typefrac_err_maglim[idx,:], 
           color=linecolor, linestyle='-' ,linewidth=2,elinewidth=1)
        if plotadj:
            ax2.errorbar(zplotbins, p_adj_typefrac_maglim[idx,:], yerr = p_adj_typefrac_err_maglim[idx,:], 
               color=linecolor, linestyle='-' ,linewidth=4,elinewidth=1)

    ax2.axvline(zhi, color='k', linestyle='--')

    plt.legend(legend_raw_str, 'upper left', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.2)
    ax2.set_xlabel('redshift')
    ax2.set_ylabel('fraction')
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)
    ax2.text(0.15,1.0,'%i galaxies' % n_maglim)

    fig.savefig(plots_path+'%s_%stype_fractions_new.png' % (task_dict['var_def'],file_str), dpi=200)

    return None

def plot_all_baselines(paperplot=False,stripe82=False,depth='normal',plot_type='pdf'):

    """ Plot the baseline morphology ratios for all 11 tasks in GZ2

    Parameters
    ----------
    paperplot : bool
        When `True`, only plot a 2x2 version of the plot containing
        data for Tasks 01-04. 

    Returns
    -------
    None

    Notes
    -------


    """

    fig = plt.figure(13,figsize=(16,10))
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

    if stripe82:
        file_str = s82_dict[depth]['str']
        mag_lim = s82_dict[depth]['maglim']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    if paperplot:
        tasklist,titlenames = tasklist[:4],titlenames[:4]
        nx_plot, ny_plot = 2,2
        left_plot = (0,2)
        bottom_plot = (2,3)

    for idx,task_dict in enumerate(tasklist):

        var_def = task_dict['var_def']

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

        ratio_baseline_masked = pickle.load(open(pkl_path+'%s_r01_%slocal_ratio_baseline_masked.pkl' % (var_def,file_str),'rb'))

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
        ax.text(-16,14,titlenames[idx],horizontalalignment='right',color='w',fontsize=16)
        ax.text(-16,12.5,'Task '+'%02i'%(idx+1),horizontalalignment='right',color='w',fontsize=16)
        if idx in bottom_plot:
            ax.set_xlabel(r'$M_r [mag]$',fontsize=26,weight='bold')
        if idx in left_plot:
            ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=26,weight='bold')
        ax.set_xticks(np.arange(-24,-15,2))
        ax.set_yticks(np.arange(0,15,2))
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        """
        if paperplot or idx==10:
            cb = plt.colorbar(im,orientation='vertical')
            cb.set_label(r'$%s(N_{%s}/N_{%s})$' % (label_prefix,task_dict['var_str'][0],task_dict['var_str'][1]),fontsize=16)
        """

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
        absmag_lim = mag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
        absmag_lim_loz = mag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = mag_lim - cosmology.dmod_flat(0.25)
        size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--', lw=3)
        ax.axhline(size_1arcsec, color='w', linestyle='dashed', lw=3)


    fig.tight_layout()
    cax = fig.add_axes([0.77,0.15,0.20,0.05])
    cb = plt.colorbar(im,cax,orientation='horizontal')
    cb.set_label('log'+r'$_{10} (f_1/f_2)$' ,fontsize=16)
    cb.set_ticks(np.arange(-1.5,2.0,0.5))
    fig.savefig(paper_figures_path+'gz2_%sbaselines.%s' % (file_str,plot_type), dpi=200)

    return None

def plot_all_type_fractions(zlo = 0.01, zhi=0.085, stripe82=False, depth='normal', paperplot=False, axlabelsize=14, titlesize=16,
                            plot_type = 'pdf'):

    """ Plot the type fractions as function of redshift for all 11 tasks in GZ2

    Parameters
    ----------
    zlo : float
        minimum redshift range to plot

    zhi : float
        redshift at which to set the magnitude limit for the sample,
        based on the apparent magnitude sensitivity at that redshift

    stripe82 : bool
        When `True`, plot the results for the Stripe 82 spectroscopic sample

    paperplot : bool
        When `True`, only plot a 2x2 version of the plot containing
        data for Tasks 01-04. 

    Returns
    -------
    None

    Notes
    -------

    Legend does not play well with tight_fit. To make papers for figure:

        comment out legend line and run
        adjust plot window to appropriate size
        run again
        comment out tight_fit line
        uncomment legend line
        run for a third time


    """

    matplotlib.use('pdf')
    fig = plt.figure(14,figsize=(17,10))
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


    if stripe82:
        file_str = s82_dict[depth]['str']
        mag_lim = s82_dict[depth]['maglim']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    if paperplot:
        tasklist,titlenames = tasklist[:4],titlenames[:4]
        nx_plot, ny_plot = 2,2
        left_plot = (0,2)
        bottom_plot = (2,3)

    for idx,task_dict in enumerate(tasklist):

        var_def = task_dict['var_def']
        ntask = len(task_dict['var_str'])

        p_raw_typefrac = pickle.load(open(pkl_path+'%s_%sraw_typefrac.pkl' % (task_dict['var_def'],file_str),'rb')) 
        p_adj_typefrac = pickle.load(open(pkl_path+'%s_%sadj_typefrac.pkl' % (task_dict['var_def'],file_str),'rb')) 
        p_raw_typefrac_maglim = pickle.load(open(pkl_path+'%s_%sraw_typefrac_maglim.pkl' % (task_dict['var_def'],file_str),'rb')) 
        p_adj_typefrac_maglim = pickle.load(open(pkl_path+'%s_%sadj_typefrac_maglim.pkl' % (task_dict['var_def'],file_str),'rb')) 

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

        zwidth = 0.02
        zplotbins = np.arange(max(edges_redshift[0],zlo),edges_redshift[-1],zwidth)
        maglimval = mag_lim - cosmology.dmod_flat(zhi)

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

        font.set_size('small')
        plt.legend([s[4:] for s in legend_raw_str], 'upper right', shadow=True, fancybox=True, prop=font)

        ax1.set_xlim(-0.01,0.19)
        ax1.set_ylim(-0.01,1.01)
        if idx in bottom_plot:
            ax1.set_xlabel('redshift',fontsize=axlabelsize)
        if idx in left_plot:
            ax1.set_ylabel('fraction',fontsize=axlabelsize)
        ax1.set_aspect('auto')
        ax1.set_title(titlenames[idx],fontsize=titlesize)
        rc('xtick', labelsize=10)
        rc('ytick', labelsize=10)

    #fig.tight_layout()
    #fig.set_tight_layout(True)
    fig.savefig(paper_figures_path+'gz2_%stype_fractions.%s' % (file_str,plot_type), dpi=200)

    return None

def posterplot_debiasing(zlimval=0.085,stripe82=False, depth='normal'):

    """ Plot the type fractions for five binary tasks in GZ2. 

    Parameters
    ----------
    zlimval : float
        redshift at which to set the magnitude limit for the sample,
        based on the apparent magnitude sensitivity at that redshift

    stripe82 : bool
        When `True`, plot the results for the Stripe 82 spectroscopic sample

    Returns
    -------
    None

    Notes
    -------

    Uses the old debiasing technique, deprecated as of Jan 2013. Made plot on
    poster presented by KW at 2013 AAS meeting. 

    """

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

    if stripe82:
        mag_lim = s82_dict[depth]['maglim']
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''
        mag_lim = appmag_lim_main

    for idx,task_dict in enumerate(tasklist):

        var_def = task_dict['var_def']
        var1_str,var2_str = task_dict['var_str']

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

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
        ax.set_xlabel(r'$M_r [mag]$',fontsize=14)
        if idx is 0:
            ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        cb = plt.colorbar(im,orientation='vertical')
        cb.set_label(r'$%s(N_{%s}/N_{%s})$' % (label_prefix,var1_str,var2_str),fontsize=18)

        plt.xticks(np.arange(edges_mag[0],edges_mag[-1],2))

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
        absmag_lim = mag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
        absmag_lim_loz = mag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = mag_lim - cosmology.dmod_flat(0.25)
        size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        ax.axhline(size_1arcsec, color='w', linestyle='dashed')

        var_def = task_dict['var_def']
        var1_str = task_dict['var_str'][0]
        var2_str = task_dict['var_str'][1]

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

        zwidth = 0.02
        zplotbins = np.arange(edges_redshift[0],edges_redshift[-1],zwidth)
        maglimval = mag_lim - cosmology.dmod_flat(zlimval)

        p_el_raw_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,file_str),'rb')) 
        p_sp_raw_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,file_str),'rb')) 
        p_el_adj_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,file_str),'rb')) 
        p_sp_adj_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,file_str),'rb')) 

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

    fig.savefig(paper_figures_path+'gz2_%sposterplot_debiasing.png' % file_str, dpi=200)

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

def poster_table(stripe82=False, depth='normal'):

    """ Return combination of debiased and raw vote fractions for
        GZ2 data used on 2013 AAS poster (now deprecated).

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -------

    Used the old debiasing technique, deprecated as of Jan 2013. Made data for
    poster presented by KW at 2013 AAS meeting. 

    """

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

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
        
        p_el_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var1_str,file_str),'rb')) 
        p_sp_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var2_str,file_str),'rb')) 
        corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (var_def,file_str),'rb')) 
        task_counts = gzdata_withz[task_dict['task_name_count']]
        
        corrgood1 = corrarr > unknown_ratio
        corrgood = corrgood1

        redshift = gzdata_withz['REDSHIFT']
        mr = gzdata_withz['PETROMAG_MR']
        r50_kpc = gzdata_withz['PETROR50_R_KPC']
        
        if var_def is 'task02':
            corrarr2 = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % ('task01',file_str),'rb')) 
            corrgood2 = corrarr2 > unknown_ratio
            corrgood = corrgood1 | corrgood2
            dep_task_wf = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % ('task01','sp',file_str),'rb'))[corrgood]
        if var_def in ('task03','task04'):
            corrarr2 = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % ('task02',file_str),'rb')) 
            corrgood2 = corrarr2 > unknown_ratio
            corrgood = corrgood1 | corrgood2
            dep_task_wf = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % ('task02','notedgeon',file_str),'rb'))[corrgood]

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

def confidence_measures(task_dict,stripe82=False, depth='normal'):

    """ Compute the confidence measures for GZ2 from Lintott et al. (2011). Not functional.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -------

    To be finished - KW, 26 Feb 2013

    """

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

    var_def = task_dict['var_def']
    centers_redshift, centers_mag, centers_size, edges_redshift, edges_mag, edges_size = get_bins(task_dict, stripe82, depth)

    corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (var_def,file_str),'rb')) 
    corrbinarr      = pickle.load(open(pkl_path+'%s_%scorrbinarr.pkl' % (var_def,file_str),'rb')) 
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

def plot_confidence_measures(task_dict,stripe82=False,depth='normal'):

    """ Plot the confidence measures for GZ2 from Lintott et al. (2011). Not functional.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -------

    To be finished - KW, 26 Feb 2013

    """

    if stripe82:
        file_str = s82_dict[depth]['str']
    else:
        file_str = ''

    corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (task_dict['var_def'],file_str),'rb')) 
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

def gz1_comparison(plothist=True,plottrumpet=True,verbose=True, plot_type='pdf'):

    """ Compare the GZ1 and GZ2 results for consistency

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -------


    """

    # Load the GZ1 and GZ2 debiased data, already matched against each other in TOPCAT

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

    # Find the differences between GZ1/GZ2 debiased classifications

    data1_el = gz1_raw_el - gz2_raw_el
    data1_sp = gz1_raw_sp - gz2_raw_sp
    data1_el_clean = (gz1_raw_el >= 0.8) & (gz2_raw_el >= 0.8)
    data1_sp_clean = (gz1_raw_sp >= 0.8) & (gz2_raw_sp >= 0.8)

    if plothist:
        #fig = plt.figure(17,figsize=(5,8))
        fig = plt.figure(17,figsize=(5,4))
        fig.clf()

        """
        ax1 = fig.add_subplot(211)
        histML(data1_el, bins=25, ax=ax1, histtype='step', color='r',weights=np.zeros_like(data1_el) + 1./data1_el.size, range=(-1.,1.), linewidth=2, linestyle='dashed')
        histML(data1_sp, bins=25, ax=ax1, histtype='step', color='b',weights=np.zeros_like(data1_sp) + 1./data1_sp.size, range=(-1.,1.), linewidth=2, linestyle='dashed')
        histML(data1_el[data1_el_clean], bins=25, ax=ax1, histtype='step', color='r',weights=np.zeros_like(data1_el[data1_el_clean]) + 1./data1_el[data1_el_clean].size, range=(-1.,1.), linewidth=1, linestyle='solid')
        histML(data1_sp[data1_sp_clean], bins=25, ax=ax1, histtype='step', color='b',weights=np.zeros_like(data1_sp[data1_sp_clean]) + 1./data1_sp[data1_sp_clean].size, range=(-1.,1.), linewidth=1, linestyle='solid')
        ax1.set_xlabel(r'$f_{GZ1} - f_{GZ2}$')
        ax1.set_ylabel('fraction of total sample')
        ax1.set_xlim(-1.0,1.0)
        ax1.set_ylim(0,0.5)
        ax1.text(-0.9,0.45,'raw',fontsize='medium', weight='bold')
        #ax1.set_title('raw votes')

        plt.legend(('el','sp','el > 0.8','sp > 0.8'), 'upper right', shadow=True, fancybox=True, prop=legfont)
        """
        
        legfont = FontProperties()
        legfont.set_size('medium')

        #ax2 = fig.add_subplot(212)
        ax2 = fig.add_subplot(111)
        data2_el = gz1_adj_el - gz2_adj_el
        data2_sp = gz1_adj_sp - gz2_adj_sp
        data2_el_clean = (gz1_adj_el >= 0.8) & (gz2_adj_el >= 0.8)
        data2_sp_clean = (gz1_adj_sp >= 0.8) & (gz2_adj_sp >= 0.8)
        #histML(data2_el,                 bins=30, ax=ax2, histtype='step', alpha=1.0, color='r',weights=np.zeros_like(data2_el)                 + 1./data2_el.size, range=(-1.,1.), linewidth=4, linestyle='dashed')
        histML(data2_sp,                 bins=30, ax=ax2, histtype='step', alpha=1.0, color='b',weights=np.zeros_like(data2_sp)                 + 1./data2_sp.size, range=(-1.,1.), linewidth=4, linestyle='dashed')
        #histML(data2_el[data2_el_clean], bins=30, ax=ax2, histtype='stepfilled', lw=1, alpha=0.4, color='r',weights=np.zeros_like(data2_el[data2_el_clean]) + 1./data2_el[data2_el_clean].size, range=(-1.,1.), linestyle='solid')
        histML(data2_sp[data2_sp_clean], bins=30, ax=ax2, histtype='stepfilled', lw=1, alpha=0.4, color='b',weights=np.zeros_like(data2_sp[data2_sp_clean]) + 1./data2_sp[data2_sp_clean].size, range=(-1.,1.), linestyle='solid')
        ax2.set_xlabel(r'$p_{sp,GZ1} - p_{sp,GZ2}$')
        ax2.set_ylabel('fraction of total sample')
        ax2.set_xlim(-1.0,1.0)
        ax2.text(-0.9,0.45,'debiased',fontsize='medium', weight='bold')
        #ax2.set_title('debiased votes')
        
        print len(data2_sp),np.sum(data2_sp_clean)
        
        #plt.legend(('el','sp','el > 0.8','sp > 0.8'), 'upper right', shadow=True, fancybox=True, prop=legfont)
        plt.legend(('all galaxies',r'$p_{sp} > 0.8$'), 'upper right', shadow=True, fancybox=True, prop=legfont)

        fig.tight_layout()
        fig.savefig(paper_figures_path+'gz1_gz2.%s' % plot_type, dpi=200)

    # Try Steven's "trumpet-style" plot

    if plottrumpet:
        fig = plt.figure(18,(5,8))
        fig.clf()

        plotheight = 0.38

        #ax1 = fig.add_subplot(211)
        #ax1 = fig.add_axes([0.08,0.20,0.4,0.75])
        ax1 = fig.add_axes([0.15,0.59,0.8,plotheight])
        H, xedges, yedges = np.histogram2d(gz1_raw_sp - gz2_raw_sp, gz1_raw_sp, 
                                           range=[[-1., 1.], [0., 1.]], bins=(50, 50))
        xcens = xedges[:-1] + (xedges[1]-xedges[0])/2.
        ycens = yedges[:-1] + (yedges[1]-yedges[0])/2.
        CS = plt.contourf(xcens,ycens,np.log10(H).T,15,cmap=plt.cm.rainbow)
        #cb = plt.colorbar(CS,orientation='horizontal')
        #cb.set_label('log '+r'$N_{gal}$',fontsize=14)
        ax1.set_xlabel(r'$f_{GZ1} - f_{GZ2}$',fontsize=20)
        ax1.set_ylabel(r'$f_{sp,GZ1}$',fontsize=20)
        ax1.set_xlim(-1,1)
        ax1.set_ylim(0,1)
        ax1.text(-0.9,0.9,'raw',fontsize='medium', weight='bold')
        #ax1.set_title('raw votes')

        #ax3 = fig.add_subplot(212)
        #ax3 = fig.add_axes([0.58,0.20,0.4,0.75])
        ax3 = fig.add_axes([0.15,0.14,0.8,plotheight])
        H, xedges, yedges = np.histogram2d(gz1_adj_sp - gz2_adj_sp, gz1_adj_sp, 
                                           range=[[-1., 1.], [0., 1.]], bins=(50, 50))
        xcens = xedges[:-1] + (xedges[1]-xedges[0])/2.
        ycens = yedges[:-1] + (yedges[1]-yedges[0])/2.
        CS = plt.contourf(xcens,ycens,np.log10(H).T,15,cmap=plt.cm.rainbow)
        ax3.set_xlabel(r'$p_{GZ1} - p_{GZ2}$',fontsize=20)
        ax3.set_ylabel(r'$p_{sp,GZ1}$',fontsize=20)
        ax3.set_ylim(0,1)
        ax3.set_xlim(-1,1)
        ax3.text(-0.9,0.9,'debiased',fontsize='medium', weight='bold')
        #ax3.set_title('debiased votes')

        #cax = fig.add_axes([0.08,0.07,0.9,0.025])
        cax = fig.add_axes([0.08,0.06,0.9,0.015])
        cb = plt.colorbar(CS,cax,orientation='horizontal')
        cb.set_label('log '+r'$N_{gal}$',fontsize=14)
        #fig.tight_layout()
        fig.savefig(paper_figures_path+'gz1_gz2_trumpet.%s' % plot_type, dpi=200)

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


    if verbose:
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

def make_tables(makefits=True,stripe82=False, depth='normal', photoz=False,latex=False,imagelist=False):

    """ Make final GZ2 data products in FITS format. 

    Parameters
    ----------
    makefits : bool
        When `True`, write the FITS files to disk

    stripe82 : bool
        When `True`, adjust probabilities for spectroscopic Stripe 82 data.

    photoz : bool
        When `True`, adjust probabilities for main sample data with
        photometric redshifts.

    latex : bool
        When `True`, print a LaTeX-formatted sample of lines to the screen.
        Intended to be copied into the paper draft. 

    imagelist : bool
        When `True`, print a sample list of objid, RA, dec to the screen.
        Used in the SDSS Image List tool and for creating the figures
        in gz2_gallery.tex

    Returns
    -------
    None

    Notes
    -------

    Results are saved as binary FITS files. They can be converted to 
        CSV and VOTables in TOPCAT.

    """

    if stripe82:

        p = pyfits.open(s82_dict[depth]['file'])
        gzdata_all = p[1].data
        p.close()

        gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]

        file_str = s82_dict[depth]['str']
        vote_threshold = s82_dict[depth]['vt']

    else:
        if photoz:
            p = pyfits.open(gz2_photoz_data_file)
            gzdata = p[1].data
            p.close()

            file_str = 'photoz_'
            vote_threshold = vt_mainsample

        else:
            p = pyfits.open(gz2_maglim_data_file)
            gzdata_all = p[1].data
            p.close()
            
            gzdata = gzdata_all[np.isfinite(gzdata_all['REDSHIFT'])]

            file_str = ''
            vote_threshold = vt_mainsample
    
    
    raw_t01 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task01',file_str),'rb')) 
    adj_t01 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task01',file_str),'rb')) 
    raw_t02 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task02',file_str),'rb')) 
    adj_t02 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task02',file_str),'rb')) 
    raw_t03 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task03',file_str),'rb')) 
    adj_t03 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task03',file_str),'rb')) 
    raw_t04 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task04',file_str),'rb')) 
    adj_t04 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task04',file_str),'rb')) 
    raw_t05 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task05',file_str),'rb')) 
    adj_t05 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task05',file_str),'rb')) 
    raw_t06 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task06',file_str),'rb')) 
    adj_t06 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task06',file_str),'rb')) 
    raw_t07 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task07',file_str),'rb')) 
    adj_t07 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task07',file_str),'rb')) 
    raw_t08 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task08',file_str),'rb')) 
    adj_t08 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task08',file_str),'rb')) 
    raw_t09 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task09',file_str),'rb')) 
    adj_t09 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task09',file_str),'rb')) 
    raw_t10 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task10',file_str),'rb')) 
    adj_t10 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task10',file_str),'rb')) 
    raw_t11 = pickle.load(open(pkl_path+'%s_%sraw_probabilities.pkl' % ('task11',file_str),'rb')) 
    adj_t11 = pickle.load(open(pkl_path+'%s_%sadj_probabilities.pkl' % ('task11',file_str),'rb')) 

    assert len(gzdata) == raw_t01.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t02.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t03.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t04.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t05.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t06.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t07.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t08.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t09.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t10.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'
    assert len(gzdata) == raw_t11.shape[1], \
        'Length of original file must be the same as that of adjusted probabilities'

    raw_t01_a01 = raw_t01[0,:]
    raw_t01_a02 = raw_t01[1,:]
    raw_t01_a03 = raw_t01[2,:]
    adj_t01_a01 = adj_t01[0,:]
    adj_t01_a02 = adj_t01[1,:]
    adj_t01_a03 = adj_t01[2,:]

    raw_t02_a04 = raw_t02[0,:]
    raw_t02_a05 = raw_t02[1,:]
    adj_t02_a04 = adj_t02[0,:]
    adj_t02_a05 = adj_t02[1,:]

    raw_t03_a06 = raw_t03[0,:]
    raw_t03_a07 = raw_t03[1,:]
    adj_t03_a06 = adj_t03[0,:]
    adj_t03_a07 = adj_t03[1,:]

    raw_t04_a08 = raw_t04[0,:]
    raw_t04_a09 = raw_t04[1,:]
    adj_t04_a08 = adj_t04[0,:]
    adj_t04_a09 = adj_t04[1,:]

    raw_t05_a10 = raw_t05[0,:]
    raw_t05_a11 = raw_t05[1,:]
    raw_t05_a12 = raw_t05[2,:]
    raw_t05_a13 = raw_t05[3,:]
    adj_t05_a10 = adj_t05[0,:]
    adj_t05_a11 = adj_t05[1,:]
    adj_t05_a12 = adj_t05[2,:]
    adj_t05_a13 = adj_t05[3,:]

    raw_t06_a14 = raw_t06[0,:]
    raw_t06_a15 = raw_t06[1,:]
    adj_t06_a14 = adj_t06[0,:]
    adj_t06_a15 = adj_t06[1,:]

    raw_t07_a16 = raw_t07[0,:]
    raw_t07_a17 = raw_t07[1,:]
    raw_t07_a18 = raw_t07[2,:]
    adj_t07_a16 = adj_t07[0,:]
    adj_t07_a17 = adj_t07[1,:]
    adj_t07_a18 = adj_t07[2,:]

    raw_t08_a19 = raw_t08[0,:]
    raw_t08_a20 = raw_t08[1,:]
    raw_t08_a21 = raw_t08[2,:]
    raw_t08_a22 = raw_t08[3,:]
    raw_t08_a23 = raw_t08[4,:]
    raw_t08_a24 = raw_t08[5,:]
    raw_t08_a38 = raw_t08[6,:]
    adj_t08_a19 = adj_t08[0,:]
    adj_t08_a20 = adj_t08[1,:]
    adj_t08_a21 = adj_t08[2,:]
    adj_t08_a22 = adj_t08[3,:]
    adj_t08_a23 = adj_t08[4,:]
    adj_t08_a24 = adj_t08[5,:]
    adj_t08_a38 = adj_t08[6,:]

    raw_t09_a25 = raw_t09[0,:]
    raw_t09_a26 = raw_t09[1,:]
    raw_t09_a27 = raw_t09[2,:]
    adj_t09_a25 = adj_t09[0,:]
    adj_t09_a26 = adj_t09[1,:]
    adj_t09_a27 = adj_t09[2,:]

    raw_t10_a28 = raw_t10[0,:]
    raw_t10_a29 = raw_t10[1,:]
    raw_t10_a30 = raw_t10[2,:]
    adj_t10_a28 = adj_t10[0,:]
    adj_t10_a29 = adj_t10[1,:]
    adj_t10_a30 = adj_t10[2,:]

    raw_t11_a31 = raw_t11[0,:]
    raw_t11_a32 = raw_t11[1,:]
    raw_t11_a33 = raw_t11[2,:]
    raw_t11_a34 = raw_t11[3,:]
    raw_t11_a36 = raw_t11[4,:]
    raw_t11_a37 = raw_t11[5,:]
    adj_t11_a31 = adj_t11[0,:]
    adj_t11_a32 = adj_t11[1,:]
    adj_t11_a33 = adj_t11[2,:]
    adj_t11_a34 = adj_t11[3,:]
    adj_t11_a36 = adj_t11[4,:]
    adj_t11_a37 = adj_t11[5,:]

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

    """
    # Old method

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
    """

    # New method

    goodt01 = (gzdata['t01_smooth_or_features_total_weight'] >= vote_threshold)
    goodt02 = (gzdata['t01_smooth_or_features_total_weight'] >= vote_threshold) & (gzdata[edgeon['dependent_tasks']] >= edgeon['vf_prev'][str(vote_threshold)])
    goodt03 = (gzdata['t02_edgeon_total_weight'] >= bar['min_classifications']) & (gzdata[bar['dependent_tasks'][-1]] >= bar['vf_prev'][str(vote_threshold)])
    goodt04 = (gzdata['t02_edgeon_total_weight'] >= spiral['min_classifications']) & (gzdata[spiral['dependent_tasks'][-1]] >= spiral['vf_prev'][str(vote_threshold)])
    goodt05 = (gzdata['t02_edgeon_total_weight'] >= bulge['min_classifications']) & (gzdata[bulge['dependent_tasks'][-1]] >= bulge['vf_prev'][str(vote_threshold)])
    goodt06 = (gzdata['t06_odd_total_weight'] >= vote_threshold)
    goodt07 = (gzdata['t01_smooth_or_features_total_weight'] >= vote_threshold) & (gzdata[rounded['dependent_tasks']] >= rounded['vf_prev'][str(vote_threshold)])
    goodt08 = (gzdata['t06_odd_total_weight'] >= odd_feature['min_classifications']) & (gzdata[odd_feature['dependent_tasks']] >= odd_feature['vf_prev'][str(vote_threshold)])
    goodt09 = (gzdata['t02_edgeon_total_weight'] >= bulge_shape['min_classifications']) & (gzdata[bulge_shape['dependent_tasks'][-1]] >= bulge_shape['vf_prev'][str(vote_threshold)])
    goodt10 = (gzdata['t04_spiral_total_weight'] >= arms_winding['min_classifications']) & (gzdata[arms_winding['dependent_tasks'][-1]] >= arms_winding['vf_prev'][str(vote_threshold)])
    goodt11 = (gzdata['t04_spiral_total_weight'] >= arms_number['min_classifications']) & (gzdata[arms_number['dependent_tasks'][-1]] >= arms_number['vf_prev'][str(vote_threshold)])

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

    flag_t01_a01 = ((adj_t01_a01 >= cleanthresh_t01) & goodt01).astype(int)
    flag_t01_a02 = ((adj_t01_a02 >= cleanthresh_t01) & goodt01).astype(int)
    flag_t01_a03 = ((adj_t01_a03 >= cleanthresh_t01) & goodt01).astype(int)
    flag_t02_a04 = ((adj_t02_a04 >= cleanthresh_t02) & goodt02).astype(int)
    flag_t02_a05 = ((adj_t02_a05 >= cleanthresh_t02) & goodt02).astype(int)
    flag_t03_a06 = ((adj_t03_a06 >= cleanthresh_t03) & goodt03).astype(int)
    flag_t03_a07 = ((adj_t03_a07 >= cleanthresh_t03) & goodt03).astype(int)
    flag_t04_a08 = ((adj_t04_a08 >= cleanthresh_t04) & goodt04).astype(int)
    flag_t04_a09 = ((adj_t04_a09 >= cleanthresh_t04) & goodt04).astype(int)
    flag_t05_a10 = ((adj_t05_a10 >= cleanthresh_t05) & goodt05).astype(int)
    flag_t05_a11 = ((adj_t05_a11 >= cleanthresh_t05) & goodt05).astype(int)
    flag_t05_a12 = ((adj_t05_a12 >= cleanthresh_t05) & goodt05).astype(int)
    flag_t05_a13 = ((adj_t05_a13 >= cleanthresh_t05) & goodt05).astype(int)
    flag_t06_a14 = ((adj_t06_a14 >= cleanthresh_t06) & goodt06).astype(int)
    flag_t06_a15 = ((adj_t06_a15 >= cleanthresh_t06) & goodt06).astype(int)
    flag_t08_a19 = ((adj_t08_a19 >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a20 = ((adj_t08_a20 >= cleanthresh_t08) & goodt08).astype(int)
    flag_t07_a16 = ((adj_t07_a16 >= cleanthresh_t07) & goodt07).astype(int)
    flag_t07_a17 = ((adj_t07_a17 >= cleanthresh_t07) & goodt07).astype(int)
    flag_t07_a18 = ((adj_t07_a18 >= cleanthresh_t07) & goodt07).astype(int)
    flag_t08_a21 = ((adj_t08_a21 >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a22 = ((adj_t08_a22 >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a23 = ((adj_t08_a23 >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a24 = ((adj_t08_a24 >= cleanthresh_t08) & goodt08).astype(int)
    flag_t08_a38 = ((adj_t08_a38 >= cleanthresh_t08) & goodt08).astype(int)
    flag_t09_a25 = ((adj_t09_a25 >= cleanthresh_t09) & goodt09).astype(int)
    flag_t09_a26 = ((adj_t09_a26 >= cleanthresh_t09) & goodt09).astype(int)
    flag_t09_a27 = ((adj_t09_a27 >= cleanthresh_t09) & goodt09).astype(int)
    flag_t10_a28 = ((adj_t10_a28 >= cleanthresh_t10) & goodt10).astype(int)
    flag_t10_a29 = ((adj_t10_a29 >= cleanthresh_t10) & goodt10).astype(int)
    flag_t10_a30 = ((adj_t10_a30 >= cleanthresh_t10) & goodt10).astype(int)
    flag_t11_a31 = ((adj_t11_a31 >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a32 = ((adj_t11_a32 >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a33 = ((adj_t11_a33 >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a34 = ((adj_t11_a34 >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a36 = ((adj_t11_a36 >= cleanthresh_t11) & goodt11).astype(int)
    flag_t11_a37 = ((adj_t11_a37 >= cleanthresh_t11) & goodt11).astype(int)

    
    if makefits:
        col_objid = pyfits.Column(name = 'dr7objid', format='K', array=gzdata['objid'])
        col_ra = pyfits.Column(name = 'ra', format='E7.3', array=gzdata['RA'])
        col_dec = pyfits.Column(name = 'dec', format='E7.3', array=gzdata['DEC'])
        col_rastring = pyfits.Column(name = 'rastring', format='A11', array=[ra_deg_to_sex(ra) for ra in gzdata['RA']])
        col_decstring = pyfits.Column(name = 'decstring', format='A11', array=[dec_deg_to_sex(dec) for dec in gzdata['DEC']])
        col_sample = pyfits.Column(name = 'sample', format='A20', array=gzdata['sample'])
        col_n_class = pyfits.Column(name = 'total_classifications', format='I4', array=gzdata['t01_smooth_or_features_total_count'])
        col_n_votes = pyfits.Column(name = 'total_votes', format='I4', array=gzdata['total_count'])

        col_t01_a01_ct = pyfits.Column(name='t01_smooth_or_features_a01_smooth_count',                             format='I4', array = gzdata['t01_smooth_or_features_a01_smooth_count'])
        col_t01_a01_wc = pyfits.Column(name='t01_smooth_or_features_a01_smooth_weight',                            format='E5.3', array = gzdata['t01_smooth_or_features_a01_smooth_weight'])
        col_t01_a01_fr = pyfits.Column(name='t01_smooth_or_features_a01_smooth_fraction',                          format='E5.3', array = gzdata['t01_smooth_or_features_a01_smooth_fraction'])
        col_t01_a01_wf = pyfits.Column(name='t01_smooth_or_features_a01_smooth_weighted_fraction',                 format='E5.3', array = gzdata['t01_smooth_or_features_a01_smooth_weighted_fraction'])
        col_t01_a01_db = pyfits.Column(name='t01_smooth_or_features_a01_smooth_debiased',                          format='E5.3', array = adj_t01_a01)
        col_t01_a01_fl = pyfits.Column(name='t01_smooth_or_features_a01_smooth_flag',                              format='I4', array = flag_t01_a01)
        col_t01_a02_ct = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_count',                   format='I4', array = gzdata['t01_smooth_or_features_a02_features_or_disk_count'])
        col_t01_a02_wc = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_weight',                  format='E5.3', array = gzdata['t01_smooth_or_features_a02_features_or_disk_weight'])
        col_t01_a02_fr = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_fraction',                format='E5.3', array = gzdata['t01_smooth_or_features_a02_features_or_disk_fraction'])
        col_t01_a02_wf = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_weighted_fraction',       format='E5.3', array = gzdata['t01_smooth_or_features_a02_features_or_disk_weighted_fraction'])
        col_t01_a02_db = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_debiased',                format='E5.3', array = adj_t01_a02)
        col_t01_a02_fl = pyfits.Column(name='t01_smooth_or_features_a02_features_or_disk_flag',                    format='I4', array = flag_t01_a02)
        col_t01_a03_ct = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_count',                   format='I4', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_count'])
        col_t01_a03_wc = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_weight',                  format='E5.3', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_weight'])
        col_t01_a03_fr = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_fraction',                format='E5.3', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_fraction'])
        col_t01_a03_wf = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_weighted_fraction',       format='E5.3', array = gzdata['t01_smooth_or_features_a03_star_or_artifact_weighted_fraction'])
        col_t01_a03_db = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_debiased',                format='E5.3', array = adj_t01_a03)
        col_t01_a03_fl = pyfits.Column(name='t01_smooth_or_features_a03_star_or_artifact_flag',                    format='I4', array = flag_t01_a03)

        col_t02_a04_ct = pyfits.Column(name='t02_edgeon_a04_yes_count',                   format='I4', array = gzdata['t02_edgeon_a04_yes_count'])
        col_t02_a04_wc = pyfits.Column(name='t02_edgeon_a04_yes_weight',                  format='E5.3', array = gzdata['t02_edgeon_a04_yes_weight'])
        col_t02_a04_fr = pyfits.Column(name='t02_edgeon_a04_yes_fraction',                format='E5.3', array = gzdata['t02_edgeon_a04_yes_fraction'])
        col_t02_a04_wf = pyfits.Column(name='t02_edgeon_a04_yes_weighted_fraction',       format='E5.3', array = gzdata['t02_edgeon_a04_yes_weighted_fraction'])
        col_t02_a04_db = pyfits.Column(name='t02_edgeon_a04_yes_debiased',                format='E5.3', array = adj_t02_a04)
        col_t02_a04_fl = pyfits.Column(name='t02_edgeon_a04_yes_flag',                    format='I4', array = flag_t02_a04)
        col_t02_a05_ct = pyfits.Column(name='t02_edgeon_a05_no_count',                    format='I4', array = gzdata['t02_edgeon_a05_no_count'])
        col_t02_a05_wc = pyfits.Column(name='t02_edgeon_a05_no_weight',                   format='E5.3', array = gzdata['t02_edgeon_a05_no_weight'])
        col_t02_a05_fr = pyfits.Column(name='t02_edgeon_a05_no_fraction',                 format='E5.3', array = gzdata['t02_edgeon_a05_no_fraction'])
        col_t02_a05_wf = pyfits.Column(name='t02_edgeon_a05_no_weighted_fraction',        format='E5.3', array = gzdata['t02_edgeon_a05_no_weighted_fraction'])
        col_t02_a05_db = pyfits.Column(name='t02_edgeon_a05_no_debiased',                 format='E5.3', array = adj_t02_a05)
        col_t02_a05_fl = pyfits.Column(name='t02_edgeon_a05_no_flag',                     format='I4', array = flag_t02_a05)

        col_t03_a06_ct = pyfits.Column(name='t03_bar_a06_bar_count',                      format='I4', array = gzdata['t03_bar_a06_bar_count'])
        col_t03_a06_wc = pyfits.Column(name='t03_bar_a06_bar_weight',                     format='E5.3', array = gzdata['t03_bar_a06_bar_weight'])
        col_t03_a06_fr = pyfits.Column(name='t03_bar_a06_bar_fraction',                   format='E5.3', array = gzdata['t03_bar_a06_bar_fraction'])
        col_t03_a06_wf = pyfits.Column(name='t03_bar_a06_bar_weighted_fraction',          format='E5.3', array = gzdata['t03_bar_a06_bar_weighted_fraction'])
        col_t03_a06_db = pyfits.Column(name='t03_bar_a06_bar_debiased',                   format='E5.3', array = adj_t03_a06)
        col_t03_a06_fl = pyfits.Column(name='t03_bar_a06_bar_flag',                       format='I4', array = flag_t03_a06)
        col_t03_a07_ct = pyfits.Column(name='t03_bar_a07_no_bar_count',                   format='I4', array = gzdata['t03_bar_a07_no_bar_count'])
        col_t03_a07_wc = pyfits.Column(name='t03_bar_a07_no_bar_weight',                  format='E5.3', array = gzdata['t03_bar_a07_no_bar_weight'])
        col_t03_a07_fr = pyfits.Column(name='t03_bar_a07_no_bar_fraction',                format='E5.3', array = gzdata['t03_bar_a07_no_bar_fraction'])
        col_t03_a07_wf = pyfits.Column(name='t03_bar_a07_no_bar_weighted_fraction',       format='E5.3', array = gzdata['t03_bar_a07_no_bar_weighted_fraction'])
        col_t03_a07_db = pyfits.Column(name='t03_bar_a07_no_bar_debiased',                format='E5.3', array = adj_t03_a07)
        col_t03_a07_fl = pyfits.Column(name='t03_bar_a07_no_bar_flag',                    format='I4', array = flag_t03_a07)

        col_t04_a08_ct = pyfits.Column(name='t04_spiral_a08_spiral_count',                    format='I4', array = gzdata['t04_spiral_a08_spiral_count'])
        col_t04_a08_wc = pyfits.Column(name='t04_spiral_a08_spiral_weight',                   format='E5.3', array = gzdata['t04_spiral_a08_spiral_weight'])
        col_t04_a08_fr = pyfits.Column(name='t04_spiral_a08_spiral_fraction',                 format='E5.3', array = gzdata['t04_spiral_a08_spiral_fraction'])
        col_t04_a08_wf = pyfits.Column(name='t04_spiral_a08_spiral_weighted_fraction',        format='E5.3', array = gzdata['t04_spiral_a08_spiral_weighted_fraction'])
        col_t04_a08_db = pyfits.Column(name='t04_spiral_a08_spiral_debiased',                 format='E5.3', array =     adj_t04_a08)
        col_t04_a08_fl = pyfits.Column(name='t04_spiral_a08_spiral_flag',                     format='I4', array =    flag_t04_a08)
        col_t04_a09_ct = pyfits.Column(name='t04_spiral_a09_no_spiral_count',                 format='I4', array = gzdata['t04_spiral_a09_no_spiral_count'])
        col_t04_a09_wc = pyfits.Column(name='t04_spiral_a09_no_spiral_weight',                format='E5.3', array = gzdata['t04_spiral_a09_no_spiral_weight'])
        col_t04_a09_fr = pyfits.Column(name='t04_spiral_a09_no_spiral_fraction',              format='E5.3', array = gzdata['t04_spiral_a09_no_spiral_fraction'])
        col_t04_a09_wf = pyfits.Column(name='t04_spiral_a09_no_spiral_weighted_fraction',     format='E5.3', array = gzdata['t04_spiral_a09_no_spiral_weighted_fraction'])
        col_t04_a09_db = pyfits.Column(name='t04_spiral_a09_no_spiral_debiased',              format='E5.3', array =     adj_t04_a09)
        col_t04_a09_fl = pyfits.Column(name='t04_spiral_a09_no_spiral_flag',                  format='I4', array =    flag_t04_a09)

        col_t05_a10_ct = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_count',                    format='I4', array = gzdata['t05_bulge_prominence_a10_no_bulge_count'])
        col_t05_a10_wc = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_weight',                   format='E5.3', array = gzdata['t05_bulge_prominence_a10_no_bulge_weight'])
        col_t05_a10_fr = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_fraction',                 format='E5.3', array = gzdata['t05_bulge_prominence_a10_no_bulge_fraction'])
        col_t05_a10_wf = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_weighted_fraction',        format='E5.3', array = gzdata['t05_bulge_prominence_a10_no_bulge_weighted_fraction'])
        col_t05_a10_db = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_debiased',                 format='E5.3', array =     adj_t05_a10)
        col_t05_a10_fl = pyfits.Column(name='t05_bulge_prominence_a10_no_bulge_flag',                     format='I4', array =    flag_t05_a10)
        col_t05_a11_ct = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_count',             format='I4', array = gzdata['t05_bulge_prominence_a11_just_noticeable_count'])
        col_t05_a11_wc = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_weight',            format='E5.3', array = gzdata['t05_bulge_prominence_a11_just_noticeable_weight'])
        col_t05_a11_fr = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_fraction',          format='E5.3', array = gzdata['t05_bulge_prominence_a11_just_noticeable_fraction'])
        col_t05_a11_wf = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_weighted_fraction', format='E5.3', array = gzdata['t05_bulge_prominence_a11_just_noticeable_weighted_fraction'])
        col_t05_a11_db = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_debiased',          format='E5.3', array =     adj_t05_a11)
        col_t05_a11_fl = pyfits.Column(name='t05_bulge_prominence_a11_just_noticeable_flag',              format='I4', array =    flag_t05_a11)
        col_t05_a12_ct = pyfits.Column(name='t05_bulge_prominence_a12_obvious_count',                     format='I4', array = gzdata['t05_bulge_prominence_a12_obvious_count'])
        col_t05_a12_wc = pyfits.Column(name='t05_bulge_prominence_a12_obvious_weight',                    format='E5.3', array = gzdata['t05_bulge_prominence_a12_obvious_weight'])
        col_t05_a12_fr = pyfits.Column(name='t05_bulge_prominence_a12_obvious_fraction',                  format='E5.3', array = gzdata['t05_bulge_prominence_a12_obvious_fraction'])
        col_t05_a12_wf = pyfits.Column(name='t05_bulge_prominence_a12_obvious_weighted_fraction',         format='E5.3', array = gzdata['t05_bulge_prominence_a12_obvious_weighted_fraction'])
        col_t05_a12_db = pyfits.Column(name='t05_bulge_prominence_a12_obvious_debiased',                  format='E5.3', array =     adj_t05_a12)
        col_t05_a12_fl = pyfits.Column(name='t05_bulge_prominence_a12_obvious_flag',                      format='I4', array =    flag_t05_a12)
        col_t05_a13_ct = pyfits.Column(name='t05_bulge_prominence_a13_dominant_count',                    format='I4', array = gzdata['t05_bulge_prominence_a13_dominant_count'])
        col_t05_a13_wc = pyfits.Column(name='t05_bulge_prominence_a13_dominant_weight',                   format='E5.3', array = gzdata['t05_bulge_prominence_a13_dominant_weight'])
        col_t05_a13_fr = pyfits.Column(name='t05_bulge_prominence_a13_dominant_fraction',                 format='E5.3', array = gzdata['t05_bulge_prominence_a13_dominant_fraction'])
        col_t05_a13_wf = pyfits.Column(name='t05_bulge_prominence_a13_dominant_weighted_fraction',        format='E5.3', array = gzdata['t05_bulge_prominence_a13_dominant_weighted_fraction'])
        col_t05_a13_db = pyfits.Column(name='t05_bulge_prominence_a13_dominant_debiased',                 format='E5.3', array =     adj_t05_a13)
        col_t05_a13_fl = pyfits.Column(name='t05_bulge_prominence_a13_dominant_flag',                     format='I4', array =    flag_t05_a13)

        col_t06_a14_ct = pyfits.Column(name='t06_odd_a14_yes_count',                    format='I4', array = gzdata['t06_odd_a14_yes_count'])
        col_t06_a14_wc = pyfits.Column(name='t06_odd_a14_yes_weight',                   format='E5.3', array = gzdata['t06_odd_a14_yes_weight'])
        col_t06_a14_fr = pyfits.Column(name='t06_odd_a14_yes_fraction',                 format='E5.3', array = gzdata['t06_odd_a14_yes_fraction'])
        col_t06_a14_wf = pyfits.Column(name='t06_odd_a14_yes_weighted_fraction',        format='E5.3', array = gzdata['t06_odd_a14_yes_weighted_fraction'])
        col_t06_a14_db = pyfits.Column(name='t06_odd_a14_yes_debiased',                 format='E5.3', array =     adj_t06_a14)
        col_t06_a14_fl = pyfits.Column(name='t06_odd_a14_yes_flag',                     format='I4', array =    flag_t06_a14)
        col_t06_a15_ct = pyfits.Column(name='t06_odd_a15_no_count',                     format='I4', array = gzdata['t06_odd_a15_no_count'])
        col_t06_a15_wc = pyfits.Column(name='t06_odd_a15_no_weight',                    format='E5.3', array = gzdata['t06_odd_a15_no_weight'])
        col_t06_a15_fr = pyfits.Column(name='t06_odd_a15_no_fraction',                  format='E5.3', array = gzdata['t06_odd_a15_no_fraction'])
        col_t06_a15_wf = pyfits.Column(name='t06_odd_a15_no_weighted_fraction',         format='E5.3', array = gzdata['t06_odd_a15_no_weighted_fraction'])
        col_t06_a15_db = pyfits.Column(name='t06_odd_a15_no_debiased',                  format='E5.3', array =     adj_t06_a15)
        col_t06_a15_fl = pyfits.Column(name='t06_odd_a15_no_flag',                      format='I4', array =    flag_t06_a15)


        col_t07_a16_ct = pyfits.Column(name='t07_rounded_a16_completely_round_count',             format='I4', array = gzdata['t07_rounded_a16_completely_round_count'])
        col_t07_a16_wc = pyfits.Column(name='t07_rounded_a16_completely_round_weight',            format='E5.3', array = gzdata['t07_rounded_a16_completely_round_weight'])
        col_t07_a16_fr = pyfits.Column(name='t07_rounded_a16_completely_round_fraction',          format='E5.3', array = gzdata['t07_rounded_a16_completely_round_fraction'])
        col_t07_a16_wf = pyfits.Column(name='t07_rounded_a16_completely_round_weighted_fraction', format='E5.3', array = gzdata['t07_rounded_a16_completely_round_weighted_fraction'])
        col_t07_a16_db = pyfits.Column(name='t07_rounded_a16_completely_round_debiased',          format='E5.3', array = adj_t07_a16)
        col_t07_a16_fl = pyfits.Column(name='t07_rounded_a16_completely_round_flag',              format='I4', array = flag_t07_a16)
        col_t07_a17_ct = pyfits.Column(name='t07_rounded_a17_in_between_count',                   format='I4', array = gzdata['t07_rounded_a17_in_between_count'])
        col_t07_a17_wc = pyfits.Column(name='t07_rounded_a17_in_between_weight',                  format='E5.3', array = gzdata['t07_rounded_a17_in_between_weight'])
        col_t07_a17_fr = pyfits.Column(name='t07_rounded_a17_in_between_fraction',                format='E5.3', array = gzdata['t07_rounded_a17_in_between_fraction'])
        col_t07_a17_wf = pyfits.Column(name='t07_rounded_a17_in_between_weighted_fraction',       format='E5.3', array = gzdata['t07_rounded_a17_in_between_weighted_fraction'])
        col_t07_a17_db = pyfits.Column(name='t07_rounded_a17_in_between_debiased',                format='E5.3', array = adj_t07_a17)
        col_t07_a17_fl = pyfits.Column(name='t07_rounded_a17_in_between_flag',                    format='I4', array = flag_t07_a17)
        col_t07_a18_ct = pyfits.Column(name='t07_rounded_a18_cigar_shaped_count',                 format='I4', array = gzdata['t07_rounded_a18_cigar_shaped_count'])
        col_t07_a18_wc = pyfits.Column(name='t07_rounded_a18_cigar_shaped_weight',                format='E5.3', array = gzdata['t07_rounded_a18_cigar_shaped_weight'])
        col_t07_a18_fr = pyfits.Column(name='t07_rounded_a18_cigar_shaped_fraction',              format='E5.3', array = gzdata['t07_rounded_a18_cigar_shaped_fraction'])
        col_t07_a18_wf = pyfits.Column(name='t07_rounded_a18_cigar_shaped_weighted_fraction',     format='E5.3', array = gzdata['t07_rounded_a18_cigar_shaped_weighted_fraction'])
        col_t07_a18_db = pyfits.Column(name='t07_rounded_a18_cigar_shaped_debiased',              format='E5.3', array = adj_t07_a18)
        col_t07_a18_fl = pyfits.Column(name='t07_rounded_a18_cigar_shaped_flag',                  format='I4', array = flag_t07_a18)

        col_t08_a19_ct = pyfits.Column(name='t08_odd_feature_a19_ring_count',                             format='I4', array = gzdata['t08_odd_feature_a19_ring_count'])
        col_t08_a19_wc = pyfits.Column(name='t08_odd_feature_a19_ring_weight',                            format='E5.3', array = gzdata['t08_odd_feature_a19_ring_weight'])
        col_t08_a19_fr = pyfits.Column(name='t08_odd_feature_a19_ring_fraction',                          format='E5.3', array = gzdata['t08_odd_feature_a19_ring_fraction'])
        col_t08_a19_wf = pyfits.Column(name='t08_odd_feature_a19_ring_weighted_fraction',                 format='E5.3', array = gzdata['t08_odd_feature_a19_ring_weighted_fraction'])
        col_t08_a19_db = pyfits.Column(name='t08_odd_feature_a19_ring_debiased',                          format='E5.3', array = adj_t08_a19)    
        col_t08_a19_fl = pyfits.Column(name='t08_odd_feature_a19_ring_flag',                              format='I4', array = flag_t08_a19)  
        col_t08_a20_ct = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_count',                      format='I4', array = gzdata['t08_odd_feature_a20_lens_or_arc_count'])
        col_t08_a20_wc = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_weight',                     format='E5.3', array = gzdata['t08_odd_feature_a20_lens_or_arc_weight'])
        col_t08_a20_fr = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_fraction',                   format='E5.3', array = gzdata['t08_odd_feature_a20_lens_or_arc_fraction'])
        col_t08_a20_wf = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_weighted_fraction',          format='E5.3', array = gzdata['t08_odd_feature_a20_lens_or_arc_weighted_fraction'])
        col_t08_a20_db = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_debiased',                   format='E5.3', array = adj_t08_a20)  
        col_t08_a20_fl = pyfits.Column(name='t08_odd_feature_a20_lens_or_arc_flag',                       format='I4', array = flag_t08_a20)
        col_t08_a21_ct = pyfits.Column(name='t08_odd_feature_a21_disturbed_count',                        format='I4', array = gzdata['t08_odd_feature_a21_disturbed_count'])
        col_t08_a21_wc = pyfits.Column(name='t08_odd_feature_a21_disturbed_weight',                       format='E5.3', array = gzdata['t08_odd_feature_a21_disturbed_weight'])
        col_t08_a21_fr = pyfits.Column(name='t08_odd_feature_a21_disturbed_fraction',                     format='E5.3', array = gzdata['t08_odd_feature_a21_disturbed_fraction'])
        col_t08_a21_wf = pyfits.Column(name='t08_odd_feature_a21_disturbed_weighted_fraction',            format='E5.3', array = gzdata['t08_odd_feature_a21_disturbed_weighted_fraction'])
        col_t08_a21_db = pyfits.Column(name='t08_odd_feature_a21_disturbed_debiased',                     format='E5.3', array = adj_t08_a21)      
        col_t08_a21_fl = pyfits.Column(name='t08_odd_feature_a21_disturbed_flag',                         format='I4', array = flag_t08_a21)    
        col_t08_a22_ct = pyfits.Column(name='t08_odd_feature_a22_irregular_count',                        format='I4', array = gzdata['t08_odd_feature_a22_irregular_count'])
        col_t08_a22_wc = pyfits.Column(name='t08_odd_feature_a22_irregular_weight',                       format='E5.3', array = gzdata['t08_odd_feature_a22_irregular_weight'])
        col_t08_a22_fr = pyfits.Column(name='t08_odd_feature_a22_irregular_fraction',                     format='E5.3', array = gzdata['t08_odd_feature_a22_irregular_fraction'])
        col_t08_a22_wf = pyfits.Column(name='t08_odd_feature_a22_irregular_weighted_fraction',            format='E5.3', array = gzdata['t08_odd_feature_a22_irregular_weighted_fraction'])
        col_t08_a22_db = pyfits.Column(name='t08_odd_feature_a22_irregular_debiased',                     format='E5.3', array = adj_t08_a22)    
        col_t08_a22_fl = pyfits.Column(name='t08_odd_feature_a22_irregular_flag',                         format='I4', array = flag_t08_a22)  
        col_t08_a23_ct = pyfits.Column(name='t08_odd_feature_a23_other_count',                            format='I4', array = gzdata['t08_odd_feature_a23_other_count'])
        col_t08_a23_wc = pyfits.Column(name='t08_odd_feature_a23_other_weight',                           format='E5.3', array = gzdata['t08_odd_feature_a23_other_weight'])
        col_t08_a23_fr = pyfits.Column(name='t08_odd_feature_a23_other_fraction',                         format='E5.3', array = gzdata['t08_odd_feature_a23_other_fraction'])
        col_t08_a23_wf = pyfits.Column(name='t08_odd_feature_a23_other_weighted_fraction',                format='E5.3', array = gzdata['t08_odd_feature_a23_other_weighted_fraction'])
        col_t08_a23_db = pyfits.Column(name='t08_odd_feature_a23_other_debiased',                         format='E5.3', array = adj_t08_a23)       
        col_t08_a23_fl = pyfits.Column(name='t08_odd_feature_a23_other_flag',                             format='I4', array = flag_t08_a23)     
        col_t08_a24_ct = pyfits.Column(name='t08_odd_feature_a24_merger_count',                           format='I4', array = gzdata['t08_odd_feature_a24_merger_count'])
        col_t08_a24_wc = pyfits.Column(name='t08_odd_feature_a24_merger_weight',                          format='E5.3', array = gzdata['t08_odd_feature_a24_merger_weight'])
        col_t08_a24_fr = pyfits.Column(name='t08_odd_feature_a24_merger_fraction',                        format='E5.3', array = gzdata['t08_odd_feature_a24_merger_fraction'])
        col_t08_a24_wf = pyfits.Column(name='t08_odd_feature_a24_merger_weighted_fraction',               format='E5.3', array = gzdata['t08_odd_feature_a24_merger_weighted_fraction'])
        col_t08_a24_db = pyfits.Column(name='t08_odd_feature_a24_merger_debiased',                        format='E5.3', array = adj_t08_a24)     
        col_t08_a24_fl = pyfits.Column(name='t08_odd_feature_a24_merger_flag',                            format='I4', array = flag_t08_a24)   
        col_t08_a38_ct = pyfits.Column(name='t08_odd_feature_a38_dust_lane_count',                        format='I4', array = gzdata['t08_odd_feature_a38_dust_lane_count'])
        col_t08_a38_wc = pyfits.Column(name='t08_odd_feature_a38_dust_lane_weight',                       format='E5.3', array = gzdata['t08_odd_feature_a38_dust_lane_weight'])
        col_t08_a38_fr = pyfits.Column(name='t08_odd_feature_a38_dust_lane_fraction',                     format='E5.3', array = gzdata['t08_odd_feature_a38_dust_lane_fraction'])
        col_t08_a38_wf = pyfits.Column(name='t08_odd_feature_a38_dust_lane_weighted_fraction',            format='E5.3', array = gzdata['t08_odd_feature_a38_dust_lane_weighted_fraction'])
        col_t08_a38_db = pyfits.Column(name='t08_odd_feature_a38_dust_lane_debiased',                     format='E5.3', array = adj_t08_a38)   
        col_t08_a38_fl = pyfits.Column(name='t08_odd_feature_a38_dust_lane_flag',                         format='I4', array = flag_t08_a38) 

        col_t09_a25_ct = pyfits.Column(name='t09_bulge_shape_a25_rounded_count',                  format='I4', array = gzdata['t09_bulge_shape_a25_rounded_count'])
        col_t09_a25_wc = pyfits.Column(name='t09_bulge_shape_a25_rounded_weight',                 format='E5.3', array = gzdata['t09_bulge_shape_a25_rounded_weight'])
        col_t09_a25_fr = pyfits.Column(name='t09_bulge_shape_a25_rounded_fraction',               format='E5.3', array = gzdata['t09_bulge_shape_a25_rounded_fraction'])
        col_t09_a25_wf = pyfits.Column(name='t09_bulge_shape_a25_rounded_weighted_fraction',      format='E5.3', array = gzdata['t09_bulge_shape_a25_rounded_weighted_fraction'])
        col_t09_a25_db = pyfits.Column(name='t09_bulge_shape_a25_rounded_debiased',               format='E5.3', array = adj_t09_a25)
        col_t09_a25_fl = pyfits.Column(name='t09_bulge_shape_a25_rounded_flag',                   format='I4', array = flag_t09_a25)
        col_t09_a26_ct = pyfits.Column(name='t09_bulge_shape_a26_boxy_count',                     format='I4', array = gzdata['t09_bulge_shape_a26_boxy_count'])
        col_t09_a26_wc = pyfits.Column(name='t09_bulge_shape_a26_boxy_weight',                    format='E5.3', array = gzdata['t09_bulge_shape_a26_boxy_weight'])
        col_t09_a26_fr = pyfits.Column(name='t09_bulge_shape_a26_boxy_fraction',                  format='E5.3', array = gzdata['t09_bulge_shape_a26_boxy_fraction'])
        col_t09_a26_wf = pyfits.Column(name='t09_bulge_shape_a26_boxy_weighted_fraction',         format='E5.3', array = gzdata['t09_bulge_shape_a26_boxy_weighted_fraction'])
        col_t09_a26_db = pyfits.Column(name='t09_bulge_shape_a26_boxy_debiased',                  format='E5.3', array = adj_t09_a26)
        col_t09_a26_fl = pyfits.Column(name='t09_bulge_shape_a26_boxy_flag',                      format='I4', array = flag_t09_a26)
        col_t09_a27_ct = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_count',                 format='I4', array = gzdata['t09_bulge_shape_a27_no_bulge_count'])
        col_t09_a27_wc = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_weight',                format='E5.3', array = gzdata['t09_bulge_shape_a27_no_bulge_weight'])
        col_t09_a27_fr = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_fraction',              format='E5.3', array = gzdata['t09_bulge_shape_a27_no_bulge_fraction'])
        col_t09_a27_wf = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_weighted_fraction',     format='E5.3', array = gzdata['t09_bulge_shape_a27_no_bulge_weighted_fraction'])
        col_t09_a27_db = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_debiased',              format='E5.3', array = adj_t09_a27)
        col_t09_a27_fl = pyfits.Column(name='t09_bulge_shape_a27_no_bulge_flag',                  format='I4', array = flag_t09_a27)

        col_t10_a28_ct = pyfits.Column(name='t10_arms_winding_a28_tight_count',                   format='I4', array = gzdata['t10_arms_winding_a28_tight_count'])
        col_t10_a28_wc = pyfits.Column(name='t10_arms_winding_a28_tight_weight',                  format='E5.3', array = gzdata['t10_arms_winding_a28_tight_weight'])
        col_t10_a28_fr = pyfits.Column(name='t10_arms_winding_a28_tight_fraction',                format='E5.3', array = gzdata['t10_arms_winding_a28_tight_fraction'])
        col_t10_a28_wf = pyfits.Column(name='t10_arms_winding_a28_tight_weighted_fraction',       format='E5.3', array = gzdata['t10_arms_winding_a28_tight_weighted_fraction'])
        col_t10_a28_db = pyfits.Column(name='t10_arms_winding_a28_tight_debiased',                format='E5.3', array = adj_t10_a28)
        col_t10_a28_fl = pyfits.Column(name='t10_arms_winding_a28_tight_flag',                    format='I4', array = flag_t10_a28)
        col_t10_a29_ct = pyfits.Column(name='t10_arms_winding_a29_medium_count',                  format='I4', array = gzdata['t10_arms_winding_a29_medium_count'])
        col_t10_a29_wc = pyfits.Column(name='t10_arms_winding_a29_medium_weight',                 format='E5.3', array = gzdata['t10_arms_winding_a29_medium_weight'])
        col_t10_a29_fr = pyfits.Column(name='t10_arms_winding_a29_medium_fraction',               format='E5.3', array = gzdata['t10_arms_winding_a29_medium_fraction'])
        col_t10_a29_wf = pyfits.Column(name='t10_arms_winding_a29_medium_weighted_fraction',      format='E5.3', array = gzdata['t10_arms_winding_a29_medium_weighted_fraction'])
        col_t10_a29_db = pyfits.Column(name='t10_arms_winding_a29_medium_debiased',               format='E5.3', array = adj_t10_a29)
        col_t10_a29_fl = pyfits.Column(name='t10_arms_winding_a29_medium_flag',                   format='I4', array = flag_t10_a29)
        col_t10_a30_ct = pyfits.Column(name='t10_arms_winding_a30_loose_count',                   format='I4', array = gzdata['t10_arms_winding_a30_loose_count'])
        col_t10_a30_wc = pyfits.Column(name='t10_arms_winding_a30_loose_weight',                  format='E5.3', array = gzdata['t10_arms_winding_a30_loose_weight'])
        col_t10_a30_fr = pyfits.Column(name='t10_arms_winding_a30_loose_fraction',                format='E5.3', array = gzdata['t10_arms_winding_a30_loose_fraction'])
        col_t10_a30_wf = pyfits.Column(name='t10_arms_winding_a30_loose_weighted_fraction',       format='E5.3', array = gzdata['t10_arms_winding_a30_loose_weighted_fraction'])
        col_t10_a30_db = pyfits.Column(name='t10_arms_winding_a30_loose_debiased',                format='E5.3', array = adj_t10_a30)
        col_t10_a30_fl = pyfits.Column(name='t10_arms_winding_a30_loose_flag',                    format='I4', array = flag_t10_a30)

        col_t11_a31_ct = pyfits.Column(name='t11_arms_number_a31_1_count',                        format='I4', array = gzdata['t11_arms_number_a31_1_count'])
        col_t11_a31_wc = pyfits.Column(name='t11_arms_number_a31_1_weight',                       format='E5.3', array = gzdata['t11_arms_number_a31_1_weight'])
        col_t11_a31_fr = pyfits.Column(name='t11_arms_number_a31_1_fraction',                     format='E5.3', array = gzdata['t11_arms_number_a31_1_fraction'])
        col_t11_a31_wf = pyfits.Column(name='t11_arms_number_a31_1_weighted_fraction',            format='E5.3', array = gzdata['t11_arms_number_a31_1_weighted_fraction'])
        col_t11_a31_db = pyfits.Column(name='t11_arms_number_a31_1_debiased',                     format='E5.3', array = adj_t11_a31)    
        col_t11_a31_fl = pyfits.Column(name='t11_arms_number_a31_1_flag',                         format='I4', array = flag_t11_a31)  
        col_t11_a32_ct = pyfits.Column(name='t11_arms_number_a32_2_count',                        format='I4', array = gzdata['t11_arms_number_a32_2_count'])
        col_t11_a32_wc = pyfits.Column(name='t11_arms_number_a32_2_weight',                       format='E5.3', array = gzdata['t11_arms_number_a32_2_weight'])
        col_t11_a32_fr = pyfits.Column(name='t11_arms_number_a32_2_fraction',                     format='E5.3', array = gzdata['t11_arms_number_a32_2_fraction'])
        col_t11_a32_wf = pyfits.Column(name='t11_arms_number_a32_2_weighted_fraction',            format='E5.3', array = gzdata['t11_arms_number_a32_2_weighted_fraction'])
        col_t11_a32_db = pyfits.Column(name='t11_arms_number_a32_2_debiased',                     format='E5.3', array = adj_t11_a32)  
        col_t11_a32_fl = pyfits.Column(name='t11_arms_number_a32_2_flag',                         format='I4', array = flag_t11_a32)
        col_t11_a33_ct = pyfits.Column(name='t11_arms_number_a33_3_count',                        format='I4', array = gzdata['t11_arms_number_a33_3_count'])
        col_t11_a33_wc = pyfits.Column(name='t11_arms_number_a33_3_weight',                       format='E5.3', array = gzdata['t11_arms_number_a33_3_weight'])
        col_t11_a33_fr = pyfits.Column(name='t11_arms_number_a33_3_fraction',                     format='E5.3', array = gzdata['t11_arms_number_a33_3_fraction'])
        col_t11_a33_wf = pyfits.Column(name='t11_arms_number_a33_3_weighted_fraction',            format='E5.3', array = gzdata['t11_arms_number_a33_3_weighted_fraction'])
        col_t11_a33_db = pyfits.Column(name='t11_arms_number_a33_3_debiased',                     format='E5.3', array = adj_t11_a33)      
        col_t11_a33_fl = pyfits.Column(name='t11_arms_number_a33_3_flag',                         format='I4', array = flag_t11_a33)    
        col_t11_a34_ct = pyfits.Column(name='t11_arms_number_a34_4_count',                        format='I4', array = gzdata['t11_arms_number_a34_4_count'])
        col_t11_a34_wc = pyfits.Column(name='t11_arms_number_a34_4_weight',                       format='E5.3', array = gzdata['t11_arms_number_a34_4_weight'])
        col_t11_a34_fr = pyfits.Column(name='t11_arms_number_a34_4_fraction',                     format='E5.3', array = gzdata['t11_arms_number_a34_4_fraction'])
        col_t11_a34_wf = pyfits.Column(name='t11_arms_number_a34_4_weighted_fraction',            format='E5.3', array = gzdata['t11_arms_number_a34_4_weighted_fraction'])
        col_t11_a34_db = pyfits.Column(name='t11_arms_number_a34_4_debiased',                     format='E5.3', array = adj_t11_a34)    
        col_t11_a34_fl = pyfits.Column(name='t11_arms_number_a34_4_flag',                         format='I4', array = flag_t11_a34)  
        col_t11_a36_ct = pyfits.Column(name='t11_arms_number_a36_more_than_4_count',              format='I4', array = gzdata['t11_arms_number_a36_more_than_4_count'])
        col_t11_a36_wc = pyfits.Column(name='t11_arms_number_a36_more_than_4_weight',             format='E5.3', array = gzdata['t11_arms_number_a36_more_than_4_weight'])
        col_t11_a36_fr = pyfits.Column(name='t11_arms_number_a36_more_than_4_fraction',           format='E5.3', array = gzdata['t11_arms_number_a36_more_than_4_fraction'])
        col_t11_a36_wf = pyfits.Column(name='t11_arms_number_a36_more_than_4_weighted_fraction',  format='E5.3', array = gzdata['t11_arms_number_a36_more_than_4_weighted_fraction'])
        col_t11_a36_db = pyfits.Column(name='t11_arms_number_a36_more_than_4_debiased',           format='E5.3', array = adj_t11_a36)       
        col_t11_a36_fl = pyfits.Column(name='t11_arms_number_a36_more_than_4_flag',               format='I4', array = flag_t11_a36)     
        col_t11_a37_ct = pyfits.Column(name='t11_arms_number_a37_cant_tell_count',                format='I4', array = gzdata['t11_arms_number_a37_cant_tell_count'])
        col_t11_a37_wc = pyfits.Column(name='t11_arms_number_a37_cant_tell_weight',               format='E5.3', array = gzdata['t11_arms_number_a37_cant_tell_weight'])
        col_t11_a37_fr = pyfits.Column(name='t11_arms_number_a37_cant_tell_fraction',             format='E5.3', array = gzdata['t11_arms_number_a37_cant_tell_fraction'])
        col_t11_a37_wf = pyfits.Column(name='t11_arms_number_a37_cant_tell_weighted_fraction',    format='E5.3', array = gzdata['t11_arms_number_a37_cant_tell_weighted_fraction'])
        col_t11_a37_db = pyfits.Column(name='t11_arms_number_a37_cant_tell_debiased',             format='E5.3', array = adj_t11_a37)     
        col_t11_a37_fl = pyfits.Column(name='t11_arms_number_a37_cant_tell_flag',                 format='I4', array = flag_t11_a37)   


        primary_hdu = pyfits.PrimaryHDU()
        hdulist = pyfits.HDUList([primary_hdu])
        
        tb1_hdu = pyfits.new_table([\
                                    col_objid,  
                                    col_ra,
                                    col_dec,
                                    col_rastring,
                                    col_decstring,
                                    col_sample,
                                    col_n_class,
                                    col_n_votes,
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
        hdulist.writeto(fits_path_main+'dr10/gz2_%sdebiased.fits' % file_str,clobber=True)    

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

def get_cat(photoz=False,stripe82=False,depth='normal',final_tables=False,sample=False):


    """ Retrieve the saved GZ2 main sample catalog

    Arguments
    ----------

    Returns
    -------
    gzdata : astropy.io.fits.FITS
        FITS table with the GZ2 data


    Example
    -------

    >>> data = gz2cat()

    Notes
    -------

    Originally written by Kyle Willett, 2013-02-27

    """

    if sample:
        p = pyfits.open(fits_path_main+'debiased_main_specz_sample.fits')
    elif final_tables:
        if stripe82:
            p = pyfits.open(gz_path+'dr10/dr10_gz2_stripe82_%s.fits' % depth)
        elif photoz:
            p = pyfits.open(gz_path+'dr10/dr10_gz2_main_photoz.fits')
        else:
            p = pyfits.open(gz_path+'dr10/dr10_gz2_main_specz.fits')
    else:
        if stripe82:
            p = pyfits.open(fits_path_main+'dr10/gz2_s82_%s_debiased.fits' % depth)
        elif photoz:
            p = pyfits.open(fits_path_main+'dr10/gz2_photoz_debiased.fits')
        else:
            p = pyfits.open(fits_path_main+'dr10/gz2_debiased.fits')

    gzdata = p[1].data
    p.close()

    return gzdata

def which_flags(data,objid=588015508218577115):

    names = data.names

    flaglist = []
    for idx,name in enumerate(names):
        stub4 = name[-4:]
        if stub4 == 'flag':
            flaglist.append(name)

    if len(flaglist) > 0:
        ind = np.where(data['objid'] == objid)
        for flag in flaglist:
            val = data[flag][ind]
            if val == 1:
                print flag[:-5]

    return None

def objlist(gzdata,response):

    """ Print random list of objid, RA, dec to the screen 
        for use with SDSS Image List tool

    Parameters
    ----------
    gzdata : astropy.io.fits.FITS
        FITS file for GZ2 produced by `make_tables`

    response : str
        Response for which to retrieve examples of clean, debiased galaxies.
        Follows same format as gz2table headers, with suffix stripped.
        Example: 
        
            >> response = 't01_smooth_or_features_a01_smooth'
        

    Returns
    -------
    None

    Notes
    -------

    Results are saved as binary FITS files. They can be converted to 
        CSV and VOTables in TOPCAT.

    """

    ra = gzdata['RA'][gzdata[response+'_flag'].astype(bool)]
    dec = gzdata['DEC'][gzdata[response+'_flag'].astype(bool)]
    
    print ' '
    print 'ra dec'
    ll = len(ra) if len(ra) < 25 else 25
    for coo in sample(zip(ra,dec),ll):
        print "%6.3f, %6.3f" % (coo[0],coo[1])
    print ' '

    return None

def which_file_str(arg):

    d = {'normal':'s82_normal_',
         'coadd1':'s82_coadd1_',
         'coadd2':'s82_coadd2_',
         'photoz':'photoz_'}

    if arg in d.keys():
        return d['arg']
    else:
        return ''

def ra_deg_to_sex(ra):

    """
    ra_deg_to_sex, dec_deg_to_sex are adapted from https://code.google.com/p/agpy/source/browse/trunk/agpy/ratosexagesimal.py 
    """

    r = ra/15.
    deg = r - (r % 1)
    deg -= (deg % 1)
    min = (r % 1) * 60.
    sec = (min % 1) * 60
    min -= (min % 1)
    rs = "%02i:%02i:%05.2f" % (deg,min,sec)
    return rs

def dec_deg_to_sex(dec):
    deg = np.sign(dec) * (abs(dec) - (abs(dec) % 1))
    deg -= (deg%1)
    min = (abs(dec) % 1) * 60.
    sec = (min % 1) * 60
    min -= (min % 1)
    decs = "%+03i:%02i:%04.1f" % (deg,min,sec)
    return decs

def specobjid_cleanup():

    """
    Finds the best match for the multiple specobjIDs returned by a CasJobs search.
    Needed to do final match for creating FITS tables w/correct specobjID. 

    Shouldn't be any galaxies missing specobjIDs in the final table, since Steven's SQL
    code indicated that they all have spectra. Either track down bug in my code or
    ask him to supply them to me directly.
    """

    f1 = fits_path_main+'gz2_specobjid_z_duplicate.fits'
    
    p1 = pyfits.open(f1)
    zd = p1[1].data
    p1.close()

    p2 = pyfits.open(gz2sample_data_file)
    gz2data = p2[1].data
    p2.close()

    import collections

    zd_collection = collections.Counter(zd['dr7objid']).items()

    best_specobjid_list = []
    dr7objid_list = []

    # Horribly inefficient loop over unique dr7objIDs

    for (d,n) in zd_collection:
        mzind = (zd['dr7objid'] == d)
        if n > 1:
            multiple_z = zd['z'][mzind]
            one_z = gz2data['REDSHIFT'][gz2data['OBJID'] == d][0]
            diffz = (multiple_z - one_z)
            dd = diffz[np.isfinite(diffz)]
            if len(dd) > 0:
                best_specobjid_list.append(zd['specobjid'][mzind][diffz == dd.min()][0])
                dr7objid_list.append(d)
            else:
                best_specobjid_list.append(-9999)
                dr7objid_list.append(d)
        else:
            best_specobjid_list.append(zd['specobjid'][mzind][0])
            dr7objid_list.append(d)
        
    # Write out as FITS binary file
    col_specobj = pyfits.Column(name='best_specobjid',format='K', array = np.array(best_specobjid_list,dtype=long))   
    col_dr7obj = pyfits.Column(name='dr7objid',format='K', array = np.array(dr7objid_list,dtype=long))   

    primary_hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([primary_hdu])
    
    tb1_hdu = pyfits.new_table([col_specobj,col_dr7obj])
    tb1_hdu.name = 'BEST_SPECOBJ'
    
    hdulist.append(tb1_hdu)
    hdulist.writeto(fits_path_main+'gz2_best_specobj.fits',clobber=True)    

    return None

def classification_histogram():

    datafile_main = gz2_maglim_data_file
    datafile_s82_normal = s82_dict['normal']['file']
    datafile_s82_coadd1 = s82_dict['coadd1']['file']
    datafile_s82_coadd2 = s82_dict['coadd2']['file']

    p = pyfits.open(datafile_main)
    main = p[1].data
    p.close()

    p = pyfits.open(datafile_s82_normal)
    normal = p[1].data
    p.close()

    p = pyfits.open(datafile_s82_coadd1)
    coadd1 = p[1].data
    p.close()

    p = pyfits.open(datafile_s82_coadd2)
    coadd2 = p[1].data
    p.close()

    data_main = main['t01_smooth_or_features_total_count']
    data_normal = normal['t01_smooth_or_features_total_count']
    data_coadd1 = coadd1['t01_smooth_or_features_total_count']
    data_coadd2 = coadd2['t01_smooth_or_features_total_count']

    fig = plt.figure(19)
    fig.clf()
    ax = fig.add_subplot(111)
    bintype = 40

    color1 = (228/255., 26/255., 28/255.) 
    color2 = (55/255., 126/255., 184/255.)
    color3 = (77/255., 175/255., 74/255.) 
    color4 = (152/255., 78/255., 163/255.)

    histML(data_main,   bins=bintype, ax=ax, alpha = 0.6, histtype='stepfilled', lw=1, color=color1, weights=np.zeros_like(data_main)   + 1./data_main.size,   range=(0,80))
    histML(data_normal, bins=bintype, ax=ax, alpha = 0.6, histtype='stepfilled', lw=1, color=color2, weights=np.zeros_like(data_normal) + 1./data_normal.size, range=(0,80))
    histML(data_coadd1, bins=bintype, ax=ax, alpha = 0.6, histtype='stepfilled', lw=1, color=color3, weights=np.zeros_like(data_coadd1) + 1./data_coadd1.size, range=(0,80))
    histML(data_coadd2, bins=bintype, ax=ax, alpha = 0.6, histtype='stepfilled', lw=1, color=color4, weights=np.zeros_like(data_coadd2) + 1./data_coadd2.size, range=(0,80))
    ax.set_xlabel('classification count', fontsize=20)
    ax.set_ylabel('fraction', fontsize=20)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)

    legfont = FontProperties()
    legfont.set_size('large')
    plt.legend(('main','stripe82_normal','stripe82_coadd1','stripe82_coadd2'), 'upper right', shadow=True, fancybox=True, prop=legfont)
    
    fig.tight_layout()
    fig.savefig(paper_figures_path+'classification_histogram.pdf', dpi=200)

    return None

def correlation_matrix():

    gz2data = get_cat()

    dbs_names = gz2data.names[12::6]
    data = np.zeros((gz2data.shape[0]),dtype=float)

    for n in dbs_names:
        data = np.vstack((data,gz2data[n]))

    data = data[1:,:]

    # Set up a mask for data insufficiently answered from previous questions

    vote_threshold = vt_mainsample

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

    goodt01 = ((gz2data['t01_smooth_or_features_a01_smooth_weight'] + gz2data['t01_smooth_or_features_a02_features_or_disk_weight'] + gz2data['t01_smooth_or_features_a03_star_or_artifact_weight']) >= vote_threshold)
    goodt02 = ((gz2data['t01_smooth_or_features_a01_smooth_weight'] + gz2data['t01_smooth_or_features_a02_features_or_disk_weight'] + gz2data['t01_smooth_or_features_a03_star_or_artifact_weight']) >= vote_threshold) & (gz2data[edgeon['dependent_tasks']] >= edgeon['vf_prev'][str(vote_threshold)])
    goodt03 = ((gz2data['t02_edgeon_a04_yes_weight'] + gz2data['t02_edgeon_a05_no_weight']) >= bar['min_classifications']) & (gz2data[bar['dependent_tasks'][-1]] >= bar['vf_prev'][str(vote_threshold)])
    goodt04 = ((gz2data['t02_edgeon_a04_yes_weight'] + gz2data['t02_edgeon_a05_no_weight']) >= spiral['min_classifications']) & (gz2data[spiral['dependent_tasks'][-1]] >= spiral['vf_prev'][str(vote_threshold)])
    goodt05 = ((gz2data['t02_edgeon_a04_yes_weight'] + gz2data['t02_edgeon_a05_no_weight']) >= bulge['min_classifications']) & (gz2data[bulge['dependent_tasks'][-1]] >= bulge['vf_prev'][str(vote_threshold)])
    goodt06 = ((gz2data['t06_odd_a14_yes_weight'] + gz2data['t06_odd_a15_no_weight']) >= vote_threshold)
    goodt07 = ((gz2data['t01_smooth_or_features_a01_smooth_weight'] + gz2data['t01_smooth_or_features_a02_features_or_disk_weight'] + gz2data['t01_smooth_or_features_a03_star_or_artifact_weight']) >= vote_threshold) & (gz2data[rounded['dependent_tasks']] >= rounded['vf_prev'][str(vote_threshold)])
    goodt08 = ((gz2data['t06_odd_a14_yes_weight'] + gz2data['t06_odd_a15_no_weight']) >= odd_feature['min_classifications']) & (gz2data[odd_feature['dependent_tasks']] >= odd_feature['vf_prev'][str(vote_threshold)])
    goodt09 = ((gz2data['t02_edgeon_a04_yes_weight'] + gz2data['t02_edgeon_a05_no_weight']) >= bar['min_classifications']) & (gz2data[bulge_shape['dependent_tasks'][-1]] >= bulge_shape['vf_prev'][str(vote_threshold)])
    goodt10 = ((gz2data['t04_spiral_a08_spiral_weight'] + gz2data['t04_spiral_a09_no_spiral_weight']) >= arms_winding['min_classifications']) & (gz2data[arms_winding['dependent_tasks'][-1]] >= arms_winding['vf_prev'][str(vote_threshold)])
    goodt11 = ((gz2data['t04_spiral_a08_spiral_weight'] + gz2data['t04_spiral_a09_no_spiral_weight']) >= arms_number['min_classifications']) & (gz2data[arms_number['dependent_tasks'][-1]] >= arms_number['vf_prev'][str(vote_threshold)])

    # Create a mask for data considered to be well-answered

    qmask = np.vstack((goodt01,goodt01,goodt01,
                      goodt02,goodt02,
                      goodt03,goodt03,
                      goodt04,goodt04,
                      goodt05,goodt05,goodt05,goodt05,
                      goodt06,goodt06,
                      goodt07,goodt07,goodt07,
                      goodt08,goodt08,goodt08,goodt08,goodt08,goodt08,goodt08,
                      goodt09,goodt09,goodt09,
                      goodt10,goodt10,goodt10,
                      goodt11,goodt11,goodt11,goodt11,goodt11,goodt11))
    
    # Apply mask to data
    dm = ma.array(data,mask=qmask)

    # Create correlation coefficient arrays for both unmasked and masked data

    R_all = ma.corrcoef(dm)
    R_masked = np.corrcoef(data)

    # Save to pkl files for plotting 

    pickle.dump(R_all, open(pkl_path+'correlation_matrix.pkl','wb')) 
    pickle.dump(R_masked, open(pkl_path+'correlation_matrix_masked.pkl','wb')) 

    return None

def correlation_matrix_plot(mask=True):

    if mask:
        R = pickle.load(open(pkl_path+'correlation_matrix_masked.pkl','rb')) 
        corrmatrix = np.ma.masked_equal(np.tril(R),0.).T[::-1]
    else:
        R = pickle.load(open(pkl_path+'correlation_matrix.pkl','rb')) 
        corrmatrix = np.ma.masked_equal(np.tril(R),0.).T[::-1]

    # Plotting the correlation matrix

    fig = plt.figure(21,figsize=(10,8))
    fig.clf()
    ax = fig.add_subplot(111)

    cm.BrBG.set_bad(color='grey')
    pc = ax.pcolormesh(corrmatrix,cmap=cm.BrBG,vmin=-1.,vmax=1.)
    ax.set_position((0.15,0.02,0.85,0.85))
    cb = plt.colorbar(pc,orientation='vertical')
    cb.set_label('linear correlation coefficient',fontsize=20)
    ax.set_xlim(0,37)
    ax.set_ylim(0,37)
    ax.axis('off')

    ax.text( 1,38,'smooth or features',    rotation=45.,fontsize='small',va='bottom')
    ax.text( 4,38,'edge-on',               rotation=45.,fontsize='small',va='bottom')
    ax.text( 6,38,'bar',                   rotation=45.,fontsize='small',va='bottom')
    ax.text( 8,38,'spiral',                rotation=45.,fontsize='small',va='bottom')
    ax.text(10,38,'bulge prominence',      rotation=45.,fontsize='small',va='bottom')
    ax.text(13,38,'anything odd',          rotation=45.,fontsize='small',va='bottom')
    ax.text(16,38,'rounded',               rotation=45.,fontsize='small',va='bottom')
    ax.text(21,38,'odd feature',           rotation=45.,fontsize='small',va='bottom')
    ax.text(26,38,'bulge shape',           rotation=45.,fontsize='small',va='bottom')
    ax.text(29,38,'arms winding',          rotation=45.,fontsize='small',va='bottom')
    ax.text(33,38,'arms number',           rotation=45.,fontsize='small',va='bottom')

    ax.text( 0 + 0.5,37.3,'s',fontsize='x-small',ha='center')
    ax.text( 1 + 0.5,37.3,'f',fontsize='x-small',ha='center')
    ax.text( 2 + 0.5,37.3,'a',fontsize='x-small',ha='center')
    ax.text( 3 + 0.5,37.3,'y',fontsize='x-small',ha='center')
    ax.text( 4 + 0.5,37.3,'n',fontsize='x-small',ha='center')
    ax.text( 5 + 0.5,37.3,'y',fontsize='x-small',ha='center')
    ax.text( 6 + 0.5,37.3,'n',fontsize='x-small',ha='center')
    ax.text( 7 + 0.5,37.3,'y',fontsize='x-small',ha='center')
    ax.text( 8 + 0.5,37.3,'n',fontsize='x-small',ha='center')
    ax.text( 9 + 0.5,37.3,'n',fontsize='x-small',ha='center')
    ax.text(10 + 0.5,37.3,'j',fontsize='x-small',ha='center')
    ax.text(11 + 0.5,37.3,'o',fontsize='x-small',ha='center')
    ax.text(12 + 0.5,37.3,'d',fontsize='x-small',ha='center')
    ax.text(13 + 0.5,37.3,'y',fontsize='x-small',ha='center')
    ax.text(14 + 0.5,37.3,'n',fontsize='x-small',ha='center')
    ax.text(15 + 0.5,37.3,'r',fontsize='x-small',ha='center')
    ax.text(16 + 0.5,37.3,'i',fontsize='x-small',ha='center')
    ax.text(17 + 0.5,37.3,'c',fontsize='x-small',ha='center')
    ax.text(18 + 0.5,37.3,'r',fontsize='x-small',ha='center')
    ax.text(19 + 0.5,37.3,'l',fontsize='x-small',ha='center')
    ax.text(20 + 0.5,37.3,'d',fontsize='x-small',ha='center')
    ax.text(21 + 0.5,37.3,'i',fontsize='x-small',ha='center')
    ax.text(22 + 0.5,37.3,'o',fontsize='x-small',ha='center')
    ax.text(23 + 0.5,37.3,'m',fontsize='x-small',ha='center')
    ax.text(24 + 0.5,37.3,'dl',fontsize='x-small',ha='center')
    ax.text(25 + 0.5,37.3,'r',fontsize='x-small',ha='center')
    ax.text(26 + 0.5,37.3,'b',fontsize='x-small',ha='center')
    ax.text(27 + 0.5,37.3,'n',fontsize='x-small',ha='center')
    ax.text(28 + 0.5,37.3,'t',fontsize='x-small',ha='center')
    ax.text(29 + 0.5,37.3,'m',fontsize='x-small',ha='center')
    ax.text(30 + 0.5,37.3,'l',fontsize='x-small',ha='center')
    ax.text(31 + 0.5,37.3,'1',fontsize='x-small',ha='center')
    ax.text(32 + 0.5,37.3,'2',fontsize='x-small',ha='center')
    ax.text(33 + 0.5,37.3,'3',fontsize='x-small',ha='center')
    ax.text(34 + 0.5,37.3,'4',fontsize='x-small',ha='center')
    ax.text(35 + 0.5,37.3,'+',fontsize='x-small',ha='center')
    ax.text(36 + 0.5,37.3,'?',fontsize='x-small',ha='center')

    ax.text(-0.3,37 -  0 - 0.5,'s', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  1 - 0.5,'f', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  2 - 0.5,'a', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  3 - 0.5,'y', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  4 - 0.5,'n', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  5 - 0.5,'y', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  6 - 0.5,'n', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  7 - 0.5,'y', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  8 - 0.5,'n', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 -  9 - 0.5,'n', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 10 - 0.5,'j', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 11 - 0.5,'o', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 12 - 0.5,'d', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 13 - 0.5,'y', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 14 - 0.5,'n', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 15 - 0.5,'r', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 16 - 0.5,'i', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 17 - 0.5,'c', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 18 - 0.5,'r', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 19 - 0.5,'l', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 20 - 0.5,'d', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 21 - 0.5,'i', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 22 - 0.5,'o', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 23 - 0.5,'m', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 24 - 0.5,'dl',fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 25 - 0.5,'r', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 26 - 0.5,'b', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 27 - 0.5,'n', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 28 - 0.5,'t', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 29 - 0.5,'m', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 30 - 0.5,'l', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 31 - 0.5,'1', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 32 - 0.5,'2', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 33 - 0.5,'3', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 34 - 0.5,'4', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 35 - 0.5,'+', fontsize='x-small',va='center',ha='right')
    ax.text(-0.3,37 - 36 - 0.5,'?', fontsize='x-small',va='center',ha='right')

    ax.text(-1,37 -  2,'smooth or features',    rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 -  4,'edge-on',               rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 -  6,'bar',                   rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 -  8,'spiral',                rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 - 11,'bulge prominence',      rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 - 14,'anything odd',          rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 - 17,'rounded',               rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 - 22,'odd feature',           rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 - 26,'bulge shape',           rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 - 29,'arms winding',          rotation=0.,fontsize='small',ha='right')
    ax.text(-1,37 - 34,'arms number',           rotation=0.,fontsize='small',ha='right')

    ax.hlines(37. -  3.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. -  5.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. -  7.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. -  9.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. - 13.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. - 15.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. - 18.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. - 25.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. - 28.,0.,37.,color='k',linestyles='dashed')
    ax.hlines(37. - 31.,0.,37.,color='k',linestyles='dashed')

    ax.vlines( 3.,0.,37.,color='k',linestyles='dashed')
    ax.vlines( 5.,0.,37.,color='k',linestyles='dashed')
    ax.vlines( 7.,0.,37.,color='k',linestyles='dashed')
    ax.vlines( 9.,0.,37.,color='k',linestyles='dashed')
    ax.vlines(13.,0.,37.,color='k',linestyles='dashed')
    ax.vlines(15.,0.,37.,color='k',linestyles='dashed')
    ax.vlines(18.,0.,37.,color='k',linestyles='dashed')
    ax.vlines(25.,0.,37.,color='k',linestyles='dashed')
    ax.vlines(28.,0.,37.,color='k',linestyles='dashed')
    ax.vlines(31.,0.,37.,color='k',linestyles='dashed')

    fig.savefig(paper_figures_path+'correlation_matrix.pdf', dpi=200)

    return None

def simulation_example_gals(data,verbose=False):

    # Pull images of example morphologies for a range of sizes, luminosities, and redshifts in GZ2

    el = (data['t01_smooth_or_features_a01_smooth_flag'] == 1)
    sp_bar_medium = (data['t03_bar_a06_bar_flag'] == 1) & (data['t10_arms_winding_a29_medium_flag'] == 1)
    sp_nobar_medium = (data['t03_bar_a07_no_bar_flag'] == 1) & (data['t10_arms_winding_a29_medium_flag'] == 1)
    edgeon = (data['t02_edgeon_a04_yes_flag'] == 1)
    merger = (data['t08_odd_feature_a24_merger_flag'] == 1)

    # Slice on redshift first

    dz = 0.02
    zarr = np.arange(0.02,0.25,dz)

    simfile_path = gz_path+'simulations/GZ2_randomish_examples/'

    for z in zarr:

        z_range = (data['REDSHIFT'] > (z - dz/2.)) & (data['REDSHIFT'] < (z + dz/2.))

        dl = 1
        lumarr = np.arange(-23,-19,dl)
        ds = 3
        sizearr = np.arange(2,15,ds)

        for size in sizearr:

            #lum_range = (data['PETROMAG_MR'] > (lum - dl/2.)) & (data['PETROMAG_MR'] < (lum + dl/2.)) 
            size_range = (data['PETROR50_R_KPC'] > (size - ds/2.)) & (data['PETROR50_R_KPC'] < (size + ds/2.)) 

            ind_el             = (el & z_range & size_range)
            ind_sp_bar_medium   = (sp_bar_medium & z_range & size_range)
            ind_sp_nobar_medium = (sp_nobar_medium & z_range & size_range)
            ind_edgeon         = (edgeon & z_range & size_range)
            ind_merger         = (merger & z_range & size_range)

            inds = (ind_el, ind_sp_bar_medium, ind_sp_nobar_medium, ind_edgeon, ind_merger)
            ind_names = ('Elliptical','Spiral barred medium','Spiral no bar medium','Edge-on','Merger')
            ind_shorts = ('ell','bar','nob','edg','mer')

            for ind_set,ind_name,ind_short in zip(inds,ind_names,ind_shorts):
                if np.sum(ind_set) >= 1:
                    ra,dec,location = sample(zip(data['ra'][ind_set],data['dec'][ind_set],data['location'][ind_set]),1)[0]
                    if verbose:
                        print "%s, z = %4.2f, size = %2i kpc, %6.3f, %6.3f" % (ind_name,z,size,ra,dec)
                    else:
                        print "%6.3f, %6.3f" % (ra,dec)

                    # Download files to folder
                    locstring = location.split('/')[-1]
                    urllib.urlretrieve(location, '%sgz2_%4.2fz_%02ikpc_%s_%s' % (simfile_path,z,size,ind_short,locstring))

                else:
                    if verbose:
                        print "No examples found in GZ2 for %s, z = %4.2f, size = %2i kpc" % (ind_name,z,size)
                    else:
                        print '0.0, 0.0'


    return None    

def make_movie_frames(data):

    # Pull images of example morphologies for a range of sizes, luminosities, and redshifts in GZ2

    edgeon       = (data['t02_edgeon_a04_yes_flag'] == 1)
    el_cigar     = (data['t07_rounded_a18_cigar_shaped_flag'] == 1)
    el_between   = (data['t07_rounded_a17_in_between_flag'] == 1)
    el_round     = (data['t07_rounded_a16_completely_round_flag'] == 1)
    sp_tight     = (data['t10_arms_winding_a28_tight_flag'] == 1)
    sp_medium    = (data['t10_arms_winding_a29_medium_flag'] == 1)
    sp_loose     = (data['t10_arms_winding_a30_loose_flag'] == 1)
    sp_bar       = (data['t03_bar_a06_bar_flag'] == 1)
    sp_nobar     = (data['t03_bar_a07_no_bar_flag'] == 1)
    sp_bulge_no  = (data['t02_edgeon_a05_no_flag'] == 1) & (data['t05_bulge_prominence_a10_no_bulge_debiased'] >= 0.5)
    sp_bulge_ju  = (data['t02_edgeon_a05_no_flag'] == 1) & (data['t05_bulge_prominence_a11_just_noticeable_debiased'] >= 0.5)
    sp_bulge_ob  = (data['t02_edgeon_a05_no_flag'] == 1) & (data['t05_bulge_prominence_a12_obvious_debiased'] >= 0.5)
    sp_bulge_do  = (data['t02_edgeon_a05_no_flag'] == 1) & (data['t05_bulge_prominence_a13_dominant_debiased'] >= 0.4)
    ring         = (data['t08_odd_feature_a19_ring_debiased'] >= 0.5) & (data['t06_odd_a14_yes_debiased'] >= 0.5)
    lens         = (data['t08_odd_feature_a20_lens_or_arc_debiased'] >= 0.3) & (data['t06_odd_a14_yes_debiased'] >= 0.5)
    disturbed    = (data['t08_odd_feature_a21_disturbed_debiased'] >= 0.5) & (data['t06_odd_a14_yes_debiased'] >= 0.5)
    irregular    = (data['t08_odd_feature_a22_irregular_debiased'] >= 0.5) & (data['t06_odd_a14_yes_debiased'] >= 0.5)
    other        = (data['t08_odd_feature_a23_other_debiased'] >= 0.5) & (data['t06_odd_a14_yes_debiased'] >= 0.5)
    dustlane     = (data['t08_odd_feature_a38_dust_lane_debiased'] >= 0.5) & (data['t06_odd_a14_yes_debiased'] >= 0.5)
    merger       = (data['t08_odd_feature_a24_merger_debiased'] >= 0.5) & (data['t06_odd_a14_yes_debiased'] >= 0.5)

    morphs = (edgeon, el_cigar, el_between, el_round, sp_tight, sp_medium, sp_loose, sp_bar, sp_nobar, sp_bulge_no, sp_bulge_ju, sp_bulge_ob, sp_bulge_do, ring, lens, disturbed, irregular, other, dustlane, merger)

    moviepath = gz_path + 'movie/sorted1000/'

    npercat = 50
    for midx,m in enumerate(morphs):
        loc = sample(data['location'][m],npercat)
        for lidx,location in enumerate(loc):
            urllib.urlretrieve(location, '%sgz2_2000_%04i.jpg' % (moviepath,(midx*npercat + lidx)))

    return None    
