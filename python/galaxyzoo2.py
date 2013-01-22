import gc
import numpy as np
import numpy.ma as ma
import cPickle as pickle
import time
import cosmology
import subprocess
import os
import warnings

from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
from scipy.stats import scoreatpercentile
from scipy.stats import tsem
from scipy.optimize import fmin_powell
from scipy.signal import convolve2d
from astropy.io import fits as pyfits
from astroML.plotting import hist as histML

"""

Code is hosted on github at willettk/galaxyzoo2.

To do:

 - find the proper bin size for each task
 - find the proper limits on probability, number of classifications per galaxy, number of galaxies per bin for each task
 - input different probabilities for different variables
 - Odd features task is not yet quantified. Sheer number of unrelated classifications make it difficult to establish a local baseline.

 - compute the correction
 - apply the correction to Task 01
 - compare results to Steven's catalog

"""

gz_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/'
fits_path = gz_path+'fits/'
pkl_path = gz_path+'pickle/'
plots_path = gz_path+'plots/'
dropbox_figs_path = gz_path+'gz2dropbox/figures/'
gz2table_data_file = fits_path+'gz2table.fits'
gz2_both_data_file = fits_path+'gz2main_table_sample.fits'
gz2_stripe82_data_file = fits_path+'gz2_stripe82_normal.fits'
gz2_full_data_file = fits_path+'gz2_original_extra_s82norm_table_sample.fits'

min_ratio = -2.0
max_ratio = 2.0
range_ratio = max_ratio - min_ratio
unknown_ratio = -999
no_correction_applied = unknown_ratio + 10
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

SBlim_app = 23.0
appmag_lim = 17.0
SBlim_size = np.arange(200) / 10.0

gc.collect()
plt.ion()

def plot_redshifts():
    
    gz2main_fitsfile = '/Users/willettk/Astronomy/Research/GalaxyZoo/fits/gz2main_table_sample.fits'
    hdulist = pyfits.open(gz2main_fitsfile)
    gz2main = hdulist[1].data
    hdulist.close()

    redshift_all = gz2main['redshift']
    redshift_finite = redshift_all[np.isfinite(redshift_all)]

    fig = plt.figure(11)
    ax = fig.add_subplot(111)
#   histML(redshift_finite, bins='blocks', ax=ax, histtype='stepfilled', color='r', ec='r', normed=True)
    histML(redshift_finite, bins='scotts', ax=ax, histtype='step', color='b',normed=True)
    histML(redshift_finite, bins='freedman', ax=ax, histtype='step', color='y',normed=True)
    histML(redshift_finite, bins='knuth', ax=ax, histtype='step', color='y',normed=True)
    ax.set_xlabel('Redshift')
    ax.set_ylabel('Frequency')
    
    #plt.show()

    return None

def weighted_parameter(task_dict, data):

    wfnames = task_dict['task_names_wf']
    n_responses = len(wfnames)
    temparr = np.zeros((n_responses,len(data)))
    for index,responses in enumerate(wfnames):
        temparr[index,:] = data[responses]

    w = np.arange(n_responses)/np.float(n_responses - 1)
    weights = np.resize(w,temparr.T.shape).T

    wbp_val = np.sum(weights * temparr,axis=0)

    return wbp_val

def bin_data_idl(task_dict,
                           zmin = 0.00, zmax = 0.26, zstep = 0.01,
                           magmin = -24, magmax = -16, magstep = 0.25,
                           sizemin = 0, sizemax = 15, sizestep = 0.5
                           ):

    tstart = time.time()
    p = pyfits.open(gz2_full_data_file)
    gzdata = p[1].data
    p.close()

    gzdata_hasredshift = gzdata[np.isfinite(gzdata['REDSHIFT'])]

    redshift = gzdata_hasredshift['REDSHIFT']
    mr = gzdata_hasredshift['PETROMAG_MR']
    r50_kpc = gzdata_hasredshift['PETROR50_R_KPC']
    r90 = gzdata_hasredshift['PETROR90_R']

    # Add masks for magnitude, size, and surface brightness limits

    dmod = cosmology.dmod_flat(redshift)
    ang_scale = cosmology.ang_scale_flat(redshift)
    absmag_lim = appmag_lim - dmod
    r90_kpc = r90 * ang_scale

    absmag_padding = 0.0
    sb_padding = 0.0

    taskcount_mask = gzdata_hasredshift[task_dict['task_name_count']] < task_dict['min_classifications']
    magnitude_mask = mr > (absmag_lim - absmag_padding)
    size_mask = r90_kpc < (3.0 * ang_scale)
    """
    Removed the surface brightness requirement - no selection algorithms for the GZ2 sample should depend
    on this, and I'm worried that the signs in the equation below are not correct. 

    surfacebrightness_mask = (mr + dmod + \
        2.5*np.log10(6.283185*(r50_kpc/ang_scale)**2)) > (SBlim_app - sb_padding)
    """

    totalmask = magnitude_mask | size_mask | taskcount_mask #| surfacebrightness_mask

    print ' '
    print '%7i galaxies removed from sample due to magnitude cutoff'%(np.sum(magnitude_mask.astype(int)))
    print '%7i galaxies removed from sample due to angular size cutoff'%(np.sum(size_mask.astype(int)))
    print '%7i galaxies removed from sample due to minimum number of classifications'%(np.sum(taskcount_mask.astype(int)))
    #print '%7i galaxies removed from sample due to surface brightness cutoff'%(np.sum(surfacebrightness_mask.astype(int)))
    print '%6.2f percent of the total is kept for %s' % ((1. - np.sum(totalmask.astype(float))/len(gzdata_hasredshift)) * 100., task_dict['var_def'])
    print ' '

    gzgood = gzdata_hasredshift[np.logical_not(totalmask)]

    # Write the masked GZ2 data to a FITS file to be binned by IDL script

    pyfits.writeto(fits_path+'%s_data_for_idl.fits' % task_dict['var_def'],
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
        idlfilecheck = not os.path.isfile(fits_path+'idlfilecreated')
        time.sleep(0.5)

    os.remove(fits_path+'idlfilecreated')

    # End of timer

    tend = time.time()
    print 'Time elapsed to bin GZ2 likelihood data for %s: %i seconds' % (task_dict['var_def'], tend - tstart)

    return None

def get_bins(task_dict):

    idl_binned = pyfits.open(fits_path+'%s_idlbinned.fits' % (task_dict['var_def']))
    bindata = idl_binned[1].data

    edges_redshift = np.squeeze(bindata['EDGES_REDSHIFT'])
    edges_mag = np.squeeze(bindata['EDGES_MAG'])
    edges_size = np.squeeze(bindata['EDGES_SIZE'])
    centers_redshift = np.squeeze(bindata['CENTERS_REDSHIFT'])
    centers_mag = np.squeeze(bindata['CENTERS_MAG'])
    centers_size = np.squeeze(bindata['CENTERS_SIZE'])

    return centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size 

def determine_ratio_baseline(task_dict, vartop = 0, varbot = 1):

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load the GZ2 sample data

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[vartop],bintype))
        p_var2 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[varbot],bintype))
        d_var1 = p_var1[0].data.astype(int)                         # Data is pre-binned
        d_var2 = p_var2[0].data.astype(int)

        zbins = p_var1['REDSHIFT_BIN_CENTERS'].data['centers']
        magbins = p_var1['MR_BIN_CENTERS'].data['centers']
        sizebins = p_var1['R50_KPC_BIN_CENTERS'].data['centers']
        
        p_var1.close()
        p_var2.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
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

    pyfits.writeto(fits_path+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot),
                   ratio_baseline, clobber=True)    
    pyfits.writeto(fits_path+'%s_r%s%s_local_counts_baseline.fits' % (var_def,vartop,varbot),
                   counts_baseline, clobber=True)    
    pyfits.writeto(fits_path+'%s_r%s%s_local_redshift_baseline.fits' % (var_def,vartop,varbot),
                   redshift_baseline, clobber=True)    

    pickle.dump(ratio_baseline_masked, open(pkl_path+'%s_r%s%s_local_ratio_baseline_masked.pkl' % (var_def,vartop,varbot),'wb')) 
    pickle.dump(counts_baseline_masked, open(pkl_path+'%s_r%s%s_local_counts_baseline_masked.pkl' % (var_def,vartop,varbot),'wb')) 
    pickle.dump(redshift_baseline_masked, open(pkl_path+'%s_r%s%s_local_redshift_baseline_masked.pkl' % (var_def,vartop,varbot),'wb')) 

    return None

def determine_ratio_baseline_sigma(task_dict, plot=False, vartop=0, varbot=1):

    var_def = task_dict['var_def']

    data = pyfits.getdata(fits_path+'%s_r%s%s_local_counts_baseline.fits' % (var_def,vartop,varbot))
    var1 = data[:,:,0].astype(np.float)
    var2 = data[:,:,1].astype(np.float)
    mask_var1 = var1 < 1
    mask_var2 = var2 < 1
    mask_all = np.logical_and(mask_var1, mask_var2)
    mask = np.logical_or(mask_var1, mask_var2)
    ratio = pyfits.getdata(fits_path+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot))
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

    pyfits.writeto(fits_path+'%s_r%s%s_local_ratio_baseline_sigma.fits' % (var_def,vartop,varbot),
           sigma, clobber=True)    

    return None

def plot_ratio_baseline(task_dict,
                        unset_vrange = False,
                        ratio_vrange = (min_ratio,max_ratio),
                        plot_mag_lims=False,
                        vartop=0, varbot=1
                        ):

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    f_ratio = fits_path+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot)
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

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load the GZ2 sample data

    """
    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var1_str,bintype))
        p_var2 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var2_str,bintype))
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
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
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

def ratio_function_flip(p, x, y):
    a, b, c, d, e, f, g, h, i = p
    z = np.zeros((len(x), len(y)), np.float)
    x0 = b**-(a + h*y**i) + c
    x1 = d + e*(x0-c)
    for i in range(len(x)):
        z[i] = (-1. * f / (1.0 + np.exp((x0 - x[i])/x1))) + g
    return z

def ratio_minfunc_flip(p, x, y, z, w):
    f = ratio_function_flip(p, x, y)
    r2 = (z - f)**2                                     # Difference between model and data
    df = (w > 0).ravel().astype(np.int).sum() - len(p)  # Degrees of freedom: N_cells - N_parameters
    s = (w * r2).sum() / df                             # Chi-squared statistic
    return s

def fit_ratio_baseline(task_dict, nboot=50, plot=False, unset_vrange=False, vartop=0, varbot=1):

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    sigma = pyfits.getdata(fits_path+'%s_r%s%s_local_ratio_baseline_sigma.fits' % (var_def,vartop,varbot))
    ratio = pyfits.getdata(fits_path+'%s_r%s%s_local_ratio_baseline.fits' % (var_def,vartop,varbot))
    
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
                             centers_mag, centers_size, ratio, weights, nboot)
    chi2r = ratio_minfunc(p, centers_mag, centers_size, ratio, weights)

    print ' '
    print 'Number of bootstrap iterations: %i' % nboot
    print 'param = %s' % p
    print 'param errors = %s' % perr
    print 'param errors (percent) = %s' % (perr/p)
    print 'chi2r = %f' % chi2r
    print ' '

    ratio_baseline_fit = ratio_function(p, centers_mag, centers_size)
    pickle.dump((p, perr), file(pkl_path+'%s_r%s%s_ratio_baseline_fit.pkl' % (var_def,vartop,varbot), 'w'))
    res = ratio - ratio_baseline_fit
    normres = res/sigma
    np.putmask(normres, weights < 0.000001, unknown_ratio)

    if plot:
        plot_ratio_function(p,centers_mag,centers_size, unset_vrange)

    pyfits.writeto(fits_path+'%s_r%s%s_local_ratio_baseline_fit.fits' % (var_def,vartop,varbot),
                   ratio_baseline_fit, clobber=True)    
    pyfits.writeto(fits_path+'%s_r%s%s_local_ratio_baseline_fitres.fits' % (var_def,vartop,varbot),
                   res, clobber=True) 
    pyfits.writeto(fits_path+'%s_r%s%s_local_ratio_baseline_fitnormres.fits' % (var_def,vartop,varbot),
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

    ratio_baseline_fit = pyfits.getdata(fits_path+'%s_r%s%s_local_ratio_baseline_fit.fits' % (var_def,vartop,varbot))

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

    fitarray_fn = ratio_function(pfit, magbinsplot, sizebinsplot)
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
                        fignum=4,
                        flip=False,
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

    if flip:
        fit = ratio_function_flip(p,x,y)
    else:
        fit = ratio_function(p,x,y)

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
                 'reverse' : False
                 }

    # Change entries for various GZ2 tasks

    if task is 'smooth':
        pass
 
    if task is 'edgeon':
        task_dict['task_name_count'] = 't02_edgeon_total_weight'
        task_dict['task_names_wf'] = ('t02_edgeon_a05_no_weighted_fraction',
                                      't02_edgeon_a04_yes_weighted_fraction')
        task_dict['var_str'] = ('notedgeon','edgeon')
        task_dict['var_def'] = 'task02'
        task_dict['min_classifications'] = 10
        task_dict['min_galperbin'] = 10
        task_dict['direct'] = True
 
    if task is 'bar':
        task_dict['task_name_count'] = 't03_bar_total_weight'
        task_dict['task_names_wf'] = ('t03_bar_a07_no_bar_weighted_fraction',
                                      't03_bar_a06_bar_weighted_fraction')
        task_dict['var_str'] = ('nobar','bar')
        task_dict['var_def'] = 'task03'
        task_dict['min_classifications'] = 10
        task_dict['min_galperbin'] = 10
        task_dict['direct'] = True
 
    if task is 'spiral':
        task_dict['task_name_count'] = 't04_spiral_total_weight'
        task_dict['task_names_wf'] = ('t04_spiral_a09_no_spiral_weighted_fraction',
                                      't04_spiral_a08_spiral_weighted_fraction')
        task_dict['var_str'] = ('nospiral','spiral')
        #task_dict['task_names_wf'] = ('t04_spiral_a08_spiral_weighted_fraction',
        #                              't04_spiral_a09_no_spiral_weighted_fraction')
        #task_dict['var_str'] = ('spiral','nospiral')
        task_dict['var_def'] = 'task04'
        task_dict['min_classifications'] = 10
        task_dict['min_galperbin'] = 10
        #task_dict['direct'] = True
 
    if task is 'odd':
        task_dict['task_name_count'] = 't06_odd_total_weight'
        task_dict['task_names_wf'] = ('t06_odd_a15_no_weighted_fraction',
                                      't06_odd_a14_yes_weighted_fraction')
        task_dict['var_str'] = ('notodd','odd')
        task_dict['var_def'] = 'task06'
        #task_dict['ratio_type'] = 'linear'
        task_dict['ratio_type'] = 'log'
        task_dict['min_classifications'] = 30
        task_dict['direct'] = True
                     
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

    if task is 'rounded':
        task_dict['task_name_count'] = 't07_rounded_total_weight'
        task_dict['task_names_wf'] = ('t07_rounded_a16_completely_round_weighted_fraction',
                                      't07_rounded_a17_in_between_weighted_fraction',
                                      't07_rounded_a18_cigar_shaped_weighted_fraction')
        task_dict['var_str'] = ('round','inbetween','cigar')
        task_dict['wp_lim'] = 0.3
        task_dict['wp_type'] = 'weighted'
        task_dict['reverse'] = False
        task_dict['ratio_type'] = 'log'
        task_dict['var_def'] = 'task07'
        task_dict['min_classifications'] = 20

    if task is 'arms_winding':
        task_dict['task_name_count'] = 't10_arms_winding_total_weight'
        task_dict['task_names_wf'] = ('t10_arms_winding_a28_tight_weighted_fraction',
                                      't10_arms_winding_a29_medium_weighted_fraction',
                                      't10_arms_winding_a30_loose_weighted_fraction')
        task_dict['var_str'] = ('tight','medium','loose')
        task_dict['var_def'] = 'task10'
        task_dict['wp_lim'] = 0.4
        task_dict['wp_type'] = 'weighted'
        task_dict['reverse'] = True
        task_dict['ratio_type'] = 'log'

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
        task_dict['reverse'] = True
        task_dict['ratio_type'] = 'log'

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

    if task is 'odd_feature':
        task_dict['task_name_count'] = 't08_odd_feature_total_weight'
        task_dict['task_names_wf'] = ('t08_bulge_shape_a19_ring_weighted_fraction',
                                      't08_bulge_shape_a20_lens_or_arc_weighted_fraction',
                                      't08_bulge_shape_a21_disturbed_weighted_fraction',
                                      't08_bulge_shape_a22_irregular_weighted_fraction',
                                      't08_bulge_shape_a23_other_weighted_fraction',
                                      't08_bulge_shape_a24_merger_weighted_fraction',
                                      't08_bulge_shape_a38_dust_lane_weighted_fraction')
        task_dict['var_str'] = ('ring','lens','disturbed','irregular','other','merger','dustlane')
        task_dict['var_def'] = 'task08'
        task_dict['min_prob'] = 0.5
        task_dict['min_classifications'] = 10

    return task_dict

def run_task(task_dict,
             nboot = 1, 
             fitbins = 1000,
             plot = True, 
             unset_vrange = True,
             vartop=0,
             varbot=1
             ):
    
    tstart = time.time()

    warnings.filterwarnings('ignore', category=UserWarning,    append=True)
    warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

    assert type(task_dict) is dict, \
        "First argument in run_task must be a dictionary -- use galaxyzoo2.get_task_dict()"

    # Execute all the fitting tasks and generate plots

    """
    if task_dict['bintype'] is 'counts':
        create_binned_data(task_dict,plot)
    """
    if task_dict['bintype'] is 'rawlikelihood':
        #create_binned_data_rawlikelihood(task_dict,plot)
        bin_data_idl(task_dict)

    determine_ratio_baseline(task_dict,vartop,varbot)
    determine_ratio_baseline_sigma(task_dict,plot,vartop,varbot)
    fit_ratio_baseline(task_dict,nboot,plot,unset_vrange,vartop,varbot)
    determine_baseline_correction(task_dict,vartop,varbot)

    if plot:
        plot_ratio_baseline(task_dict,unset_vrange,vartop=vartop,varbot=varbot)
        plot_ratio_baseline_fit(task_dict,fitbins,unset_vrange,vartop=vartop,varbot=varbot)
        plot_ratio_baseline_redshift(task_dict,vartop,varbot)
        plot_baseline_correction(task_dict,vartop=vartop,varbot=varbot)
        plot_galaxy_counts(task_dict)

    if len(task_dict['task_names_wf']) == 2:
        adjust_probabilities(task_dict)
        plot_type_fractions(task_dict)

    if len(task_dict['task_names_wf']) == 3:
        adjust_probabilities_n3(task_dict)
        plot_type_fractions_n3(task_dict)

    if len(task_dict['task_names_wf']) == 4:
        adjust_probabilities_n4(task_dict)
        plot_type_fractions_n4(task_dict)

    if len(task_dict['task_names_wf']) == 5:
        adjust_probabilities_n5(task_dict)
        plot_type_fractions_n5(task_dict)

    if len(task_dict['task_names_wf']) == 6:
        adjust_probabilities_n6(task_dict)
        plot_type_fractions_n6(task_dict)

    warnings.resetwarnings()

    tend = time.time()
    print 'Time elapsed for run_task on %s: %i seconds' % (task_dict['task_str'], tend - tstart)

    return None

def run_all_tasks(bintype='rawlikelihood'):

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
                'bulge_shape']

    for task in tasklist:
        td = get_task_dict(task)
        td['bintype']=bintype
        run_task(td)

    plot_all_baselines()
    plot_all_type_fractions()
    posterplot_debiasing()

    tend = time.time()
    print 'Time elapsed for run_all_tasks: %i seconds' % (tend - tstart)

    return None

def determine_baseline_correction(task_dict, vartop = 0, varbot = 1, spb=False):

    # Load best fit data

    var_def  = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[vartop],bintype))
        p_var2 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[varbot],bintype))
        d_var1 = p_var1[0].data.astype(int)                         # Data is pre-binned
        d_var2 = p_var2[0].data.astype(int)

        zbins = p_var1['REDSHIFT_BIN_CENTERS'].data['centers']
        magbins = p_var1['MR_BIN_CENTERS'].data['centers']
        sizebins = p_var1['R50_KPC_BIN_CENTERS'].data['centers']
        
        p_var1.close()
        p_var2.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        d_var1 = np.squeeze(d_allvar[vartop,:,:,:])
        d_var2 = np.squeeze(d_allvar[varbot,:,:,:])

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    # Steven's correction from GZ1 for Task 01
    if spb:
        pfit_params = np.array([ -3.76894097,1.8810826,-22.66409396,0.16529829,0.28365128,-2.15835938,1.04379026,1.22666975,0.85754666])

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

    if task_dict['direct']:
        correction = ratio - ratio_baseline_masked
    else:
        pfit_params, pfit_params_err = pickle.load(file(pkl_path+'%s_r%s%s_ratio_baseline_fit.pkl' % (var_def,vartop,varbot), 'r'))
        pfit = ratio_function(pfit_params, centers_mag, centers_size)
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
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[vartop],bintype))
        p_var2 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[varbot],bintype))
        d_var1 = p_var1[0].data.astype(int)                         # Data is pre-binned
        d_var2 = p_var2[0].data.astype(int)

        zbins = p_var1['REDSHIFT_BIN_CENTERS'].data['centers']
        edges_mag = p_var1['MR_BIN_EDGES'].data['edges']
        edges_size = p_var1['R50_KPC_BIN_EDGES'].data['edges']

        p_var1.close()
        p_var2.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
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
        fit = ratio_function(pfit, np.linspace(edges_mag[0],edges_mag[-1],1000),np.linspace(edges_size[0],edges_size[-1],1000))
    else:
        fit = ratio_function(pfit, centers_mag, centers_size)
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

def p_x(p_el,p_sp):

    p_x = 1. - p_el - p_sp

    return p_x

def ratio_adj(p_el,p_sp,correction,ratio_type):
    
    assert ratio_type in ('linear','log'), \
        "Ratio type must be defined and either 'linear' or 'log'"

    if ratio_type is 'log':
        ratio_adj = (p_el / p_sp) / (10.**(correction))

    if ratio_type is 'linear':
        ratio_adj = (p_el / p_sp) / correction

    return ratio_adj

def p_el_adj(p_el,p_sp,correction,ratio_type):

    denom_el = 1./ratio_adj(p_el,p_sp,correction,ratio_type) + (p_x(p_el,p_sp) / p_el) + 1.
    p_el_adj = 1./denom_el

    return p_el_adj

def p_sp_adj(p_el,p_sp,correction,ratio_type):

    denom_sp = ratio_adj(p_el,p_sp,correction,ratio_type) + (p_x(p_el,p_sp) / p_sp) + 1.
    p_sp_adj = 1./denom_sp

    return p_sp_adj

def p_n6_adj(f_i, f_j, f_k, f_l, f_m, f_n, corr_ji, corr_ki, corr_li, corr_mi, corr_ni, ratio_type):

    K_ji = 10.**(corr_ji) 
    term_ji = K_ji * (f_j / f_i)
    K_ki = 10.**(corr_ki) 
    term_ki = K_ki * (f_k / f_i)
    K_li = 10.**(corr_li) 
    term_li = K_li * (f_l / f_i)
    K_mi = 10.**(corr_mi) 
    term_mi = K_mi * (f_m / f_i)
    K_ni = 10.**(corr_ni) 
    term_ni = K_ni * (f_n / f_i)

    denom = term_ji + term_ki + term_li + term_mi + term_ni + 1.

    p_adj = 1. / denom

    return p_adj

def p_n5_adj(f_i, f_j, f_k, f_l, f_m, corr_ji, corr_ki, corr_li, corr_mi, ratio_type):

    K_ji = 10.**(corr_ji) 
    term_ji = K_ji * (f_j / f_i)
    K_ki = 10.**(corr_ki) 
    term_ki = K_ki * (f_k / f_i)
    K_li = 10.**(corr_li) 
    term_li = K_li * (f_l / f_i)
    K_mi = 10.**(corr_mi) 
    term_mi = K_mi * (f_m / f_i)

    denom = term_ji + term_ki + term_li + term_mi + 1.

    p_adj = 1. / denom

    return p_adj

def p_n4_adj(f_i, f_j, f_k, f_l, corr_ji, corr_ki, corr_li, ratio_type):

    K_ji = 10.**(corr_ji) 
    term_ji = K_ji * (f_j / f_i)
    K_ki = 10.**(corr_ki) 
    term_ki = K_ki * (f_k / f_i)
    K_li = 10.**(corr_li) 
    term_li = K_li * (f_l / f_i)

    denom = term_ji + term_ki + term_li + 1.

    p_adj = 1. / denom

    return p_adj

def p_n3_adj(f_i, f_j, f_k, corr_ji, corr_ki, ratio_type):

    K_ji = 10.**(corr_ji) 
    term_ji = K_ji * (f_j / f_i)
    K_ki = 10.**(corr_ki) 
    term_ki = K_ki * (f_k / f_i)

    denom = term_ji + term_ki + 1.

    p_adj = 1. / denom

    return p_adj

def adjust_probabilities_n6(task_dict, plot=False, stripe82=False):

    # Load in the raw probabilities from the GZ2 catalog

    if stripe82:
        p = pyfits.open(gz2_stripe82_data_file)
        s82_str = 'stripe82_'
    else:
        p = pyfits.open(gz2_full_data_file)
        s82_str = ''

    gzdata = p[1].data
    p.close()

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load in the bin sizes from the fits

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[0],bintype))

        edges_redshift = p_var1['REDSHIFT_BIN_EDGES'].data['edges']
        edges_mag = p_var1['MR_BIN_EDGES'].data['edges']
        edges_size = p_var1['R50_KPC_BIN_EDGES'].data['edges']

        p_var1.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    """
    Take each galaxy, find which bin it corresponds to, adjust probability via Bamford+09 method
    """

    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']

    zstep = edges_redshift[1] - edges_redshift[0]
    magstep = edges_mag[1] - edges_mag[0]
    sizestep = edges_size[1] - edges_size[0]

    # Loop over combinations of tasks

    corr_12, corrmasked_12, ratiomasked_12 = determine_baseline_correction(task_dict, 0, 1)
    corr_13, corrmasked_13, ratiomasked_13 = determine_baseline_correction(task_dict, 0, 2)
    corr_14, corrmasked_14, ratiomasked_14 = determine_baseline_correction(task_dict, 0, 3)
    corr_15, corrmasked_15, ratiomasked_15 = determine_baseline_correction(task_dict, 0, 4)
    corr_16, corrmasked_16, ratiomasked_16 = determine_baseline_correction(task_dict, 0, 5)
    corr_21, corrmasked_21, ratiomasked_21 = determine_baseline_correction(task_dict, 1, 0)
    corr_23, corrmasked_23, ratiomasked_23 = determine_baseline_correction(task_dict, 1, 2)
    corr_24, corrmasked_24, ratiomasked_24 = determine_baseline_correction(task_dict, 1, 3)
    corr_25, corrmasked_25, ratiomasked_25 = determine_baseline_correction(task_dict, 1, 4)
    corr_26, corrmasked_26, ratiomasked_26 = determine_baseline_correction(task_dict, 1, 5)
    corr_32, corrmasked_32, ratiomasked_32 = determine_baseline_correction(task_dict, 2, 1)
    corr_31, corrmasked_31, ratiomasked_31 = determine_baseline_correction(task_dict, 2, 0)
    corr_34, corrmasked_34, ratiomasked_34 = determine_baseline_correction(task_dict, 2, 3)
    corr_35, corrmasked_35, ratiomasked_35 = determine_baseline_correction(task_dict, 2, 4)
    corr_36, corrmasked_36, ratiomasked_36 = determine_baseline_correction(task_dict, 2, 5)
    corr_42, corrmasked_42, ratiomasked_42 = determine_baseline_correction(task_dict, 3, 1)
    corr_43, corrmasked_43, ratiomasked_43 = determine_baseline_correction(task_dict, 3, 2)
    corr_41, corrmasked_41, ratiomasked_41 = determine_baseline_correction(task_dict, 3, 0)
    corr_45, corrmasked_45, ratiomasked_45 = determine_baseline_correction(task_dict, 3, 4)
    corr_46, corrmasked_46, ratiomasked_46 = determine_baseline_correction(task_dict, 3, 5)
    corr_52, corrmasked_52, ratiomasked_52 = determine_baseline_correction(task_dict, 4, 1)
    corr_53, corrmasked_53, ratiomasked_53 = determine_baseline_correction(task_dict, 4, 2)
    corr_54, corrmasked_54, ratiomasked_54 = determine_baseline_correction(task_dict, 4, 3)
    corr_51, corrmasked_51, ratiomasked_51 = determine_baseline_correction(task_dict, 4, 0)
    corr_56, corrmasked_56, ratiomasked_56 = determine_baseline_correction(task_dict, 4, 5)
    corr_62, corrmasked_62, ratiomasked_62 = determine_baseline_correction(task_dict, 5, 1)
    corr_63, corrmasked_63, ratiomasked_63 = determine_baseline_correction(task_dict, 5, 2)
    corr_64, corrmasked_64, ratiomasked_64 = determine_baseline_correction(task_dict, 5, 3)
    corr_65, corrmasked_65, ratiomasked_65 = determine_baseline_correction(task_dict, 5, 4)
    corr_61, corrmasked_61, ratiomasked_61 = determine_baseline_correction(task_dict, 5, 0)

    p_el_raw = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw = gzdata_withz[task_dict['task_names_wf'][2]]      # Designate the third variable (for now) as 'xx'
    p_yy_raw = gzdata_withz[task_dict['task_names_wf'][3]]      
    p_zz_raw = gzdata_withz[task_dict['task_names_wf'][4]]     
    p_ww_raw = gzdata_withz[task_dict['task_names_wf'][5]]     
    p_task_counts = gzdata_withz[task_dict['task_name_count']]

    # Start adjusting

    timestart2 = time.time()

    p_el_adj_arr2 = np.zeros_like(p_el_raw)
    p_sp_adj_arr2 = np.zeros_like(p_sp_raw)
    p_xx_adj_arr2 = np.zeros_like(p_xx_raw)
    p_yy_adj_arr2 = np.zeros_like(p_yy_raw)
    p_zz_adj_arr2 = np.zeros_like(p_zz_raw)
    p_ww_adj_arr2 = np.zeros_like(p_ww_raw)
    corrarr = np.zeros_like(p_sp_raw)
    corrbinarr = np.zeros((3,len(p_el_raw)))                    # Three bins give the (z,M,r) coordinates of the matching bin

    corrshape = corr_12.shape                                   # Should probably assert that corrections are all the same shape.
    allvotecount = 0
    outbincount = 0
    corrcount = 0
    unknowncount = 0

    # Loop over all galaxies in the GZ catalog with redshifts

    for idx, (el, sp, xx, yy, zz, ww) in enumerate(zip(p_el_raw, p_sp_raw, p_xx_raw, p_yy_raw, p_zz_raw, p_ww_raw)):

        z2 = redshift[idx]
        m2 = mr[idx]
        s2 = r50_kpc[idx]

        zbin = np.abs(z2 >= edges_redshift).argmin() - 1
        mbin = np.abs(m2 >= edges_mag).argmin() - 1
        sbin = np.abs(s2 >= edges_size).argmin() - 1

        if (zbin < corrshape[0]) and (mbin < corrshape[1]) and (sbin < corrshape[2]):
            applycorr_12 = 0. if corr_12[zbin,mbin,sbin] == unknown_ratio else corr_12[zbin,mbin,sbin]
            applycorr_13 = 0. if corr_13[zbin,mbin,sbin] == unknown_ratio else corr_13[zbin,mbin,sbin]
            applycorr_14 = 0. if corr_14[zbin,mbin,sbin] == unknown_ratio else corr_14[zbin,mbin,sbin]
            applycorr_15 = 0. if corr_15[zbin,mbin,sbin] == unknown_ratio else corr_15[zbin,mbin,sbin]
            applycorr_16 = 0. if corr_16[zbin,mbin,sbin] == unknown_ratio else corr_16[zbin,mbin,sbin]
            applycorr_21 = 0. if corr_21[zbin,mbin,sbin] == unknown_ratio else corr_21[zbin,mbin,sbin]
            applycorr_23 = 0. if corr_23[zbin,mbin,sbin] == unknown_ratio else corr_23[zbin,mbin,sbin]
            applycorr_24 = 0. if corr_24[zbin,mbin,sbin] == unknown_ratio else corr_24[zbin,mbin,sbin]
            applycorr_25 = 0. if corr_25[zbin,mbin,sbin] == unknown_ratio else corr_25[zbin,mbin,sbin]
            applycorr_26 = 0. if corr_26[zbin,mbin,sbin] == unknown_ratio else corr_26[zbin,mbin,sbin]
            applycorr_32 = 0. if corr_32[zbin,mbin,sbin] == unknown_ratio else corr_32[zbin,mbin,sbin]
            applycorr_31 = 0. if corr_31[zbin,mbin,sbin] == unknown_ratio else corr_31[zbin,mbin,sbin]
            applycorr_34 = 0. if corr_34[zbin,mbin,sbin] == unknown_ratio else corr_34[zbin,mbin,sbin]
            applycorr_35 = 0. if corr_35[zbin,mbin,sbin] == unknown_ratio else corr_35[zbin,mbin,sbin]
            applycorr_36 = 0. if corr_36[zbin,mbin,sbin] == unknown_ratio else corr_36[zbin,mbin,sbin]
            applycorr_42 = 0. if corr_42[zbin,mbin,sbin] == unknown_ratio else corr_42[zbin,mbin,sbin]
            applycorr_43 = 0. if corr_43[zbin,mbin,sbin] == unknown_ratio else corr_43[zbin,mbin,sbin]
            applycorr_41 = 0. if corr_41[zbin,mbin,sbin] == unknown_ratio else corr_41[zbin,mbin,sbin]
            applycorr_45 = 0. if corr_45[zbin,mbin,sbin] == unknown_ratio else corr_45[zbin,mbin,sbin]
            applycorr_46 = 0. if corr_46[zbin,mbin,sbin] == unknown_ratio else corr_46[zbin,mbin,sbin]
            applycorr_52 = 0. if corr_52[zbin,mbin,sbin] == unknown_ratio else corr_52[zbin,mbin,sbin]
            applycorr_53 = 0. if corr_53[zbin,mbin,sbin] == unknown_ratio else corr_53[zbin,mbin,sbin]
            applycorr_54 = 0. if corr_54[zbin,mbin,sbin] == unknown_ratio else corr_54[zbin,mbin,sbin]
            applycorr_51 = 0. if corr_51[zbin,mbin,sbin] == unknown_ratio else corr_51[zbin,mbin,sbin]
            applycorr_56 = 0. if corr_56[zbin,mbin,sbin] == unknown_ratio else corr_56[zbin,mbin,sbin]
            applycorr_62 = 0. if corr_62[zbin,mbin,sbin] == unknown_ratio else corr_62[zbin,mbin,sbin]
            applycorr_63 = 0. if corr_63[zbin,mbin,sbin] == unknown_ratio else corr_63[zbin,mbin,sbin]
            applycorr_64 = 0. if corr_64[zbin,mbin,sbin] == unknown_ratio else corr_64[zbin,mbin,sbin]
            applycorr_65 = 0. if corr_65[zbin,mbin,sbin] == unknown_ratio else corr_65[zbin,mbin,sbin]
            applycorr_61 = 0. if corr_61[zbin,mbin,sbin] == unknown_ratio else corr_61[zbin,mbin,sbin]

            allcorrs = np.array((applycorr_12,
                                 applycorr_13,
                                 applycorr_14,
                                 applycorr_15,
                                 applycorr_16,
                                 applycorr_21,
                                 applycorr_23,
                                 applycorr_24,
                                 applycorr_25,
                                 applycorr_26,
                                 applycorr_32,
                                 applycorr_31,
                                 applycorr_34,
                                 applycorr_35,
                                 applycorr_36,
                                 applycorr_42,
                                 applycorr_43,
                                 applycorr_41,
                                 applycorr_45,
                                 applycorr_46,
                                 applycorr_52,
                                 applycorr_53,
                                 applycorr_54,
                                 applycorr_51,
                                 applycorr_56,
                                 applycorr_62,
                                 applycorr_63,
                                 applycorr_64,
                                 applycorr_65,
                                 applycorr_61
                                )
                               )

            if np.sum(allcorrs) == 0.:
                unknowncount +=1
                p_el_new = el
                p_sp_new = sp
                p_xx_new = xx
                p_yy_new = yy
                p_zz_new = zz
                p_ww_new = ww
                corrarr[idx] = unknown_ratio
                corrbinarr[:,idx] = unknown_ratio
            else:
                corrcount += 1
                p_el_new = p_n6_adj(el,sp,xx,yy,zz,ww,
                                    applycorr_12,applycorr_13,applycorr_14,applycorr_15,applycorr_16,
                                    task_dict['ratio_type']) if el > 0. else el
                p_sp_new = p_n6_adj(sp,el,xx,yy,zz,ww,
                                    applycorr_21,applycorr_23,applycorr_24,applycorr_25,applycorr_26,
                                    task_dict['ratio_type']) if sp > 0. else sp
                p_xx_new = p_n6_adj(xx,el,sp,yy,zz,ww,
                                    applycorr_31,applycorr_32,applycorr_34,applycorr_35,applycorr_36,
                                    task_dict['ratio_type']) if xx > 0. else xx
                p_yy_new = p_n6_adj(yy,el,sp,xx,zz,ww,
                                    applycorr_41,applycorr_42,applycorr_43,applycorr_45,applycorr_46,
                                    task_dict['ratio_type']) if yy > 0. else yy
                p_zz_new = p_n6_adj(zz,el,sp,xx,yy,ww,
                                    applycorr_51,applycorr_52,applycorr_53,applycorr_54,applycorr_56,
                                    task_dict['ratio_type']) if zz > 0. else zz
                p_ww_new = p_n6_adj(ww,el,sp,xx,yy,zz,
                                    applycorr_61,applycorr_62,applycorr_63,applycorr_64,applycorr_65,
                                    task_dict['ratio_type']) if ww > 0. else ww
                corrbinarr[:,idx] = zbin,mbin,sbin
                corrarr[idx] = applycorr_12

            p_el_adj_arr2[idx] = p_el_new
            p_sp_adj_arr2[idx] = p_sp_new
            p_xx_adj_arr2[idx] = p_xx_new
            p_yy_adj_arr2[idx] = p_yy_new
            p_zz_adj_arr2[idx] = p_zz_new
            p_ww_adj_arr2[idx] = p_ww_new
          
        else: # if galaxy does not appear in the correction area
            outbincount += 1
            p_el_adj_arr2[idx] = el
            p_sp_adj_arr2[idx] = sp
            p_xx_adj_arr2[idx] = xx
            p_yy_adj_arr2[idx] = yy
            p_zz_adj_arr2[idx] = zz
            p_ww_adj_arr2[idx] = ww
            corrbinarr[:,idx] = unknown_ratio
            corrarr[idx] = unknown_ratio


    timeend2 = time.time()
    print ' '
    #print '%7i galaxies had 100 percent of the vote for %s or %s' % (allvotecount,var_str[vartop],var_str[varbot])
    print '%7i galaxies had correction bins larger than volume' % outbincount
    print '%7i galaxies corrected' % corrcount
    print 'Time elapsed to adjust probabilities: %i seconds' % (timeend2 - timestart2)
    print ' '

    p.close()

    pickle.dump(p_el_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][0], s82_str),'wb')) 
    pickle.dump(p_sp_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][1], s82_str),'wb')) 
    pickle.dump(p_xx_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][2], s82_str),'wb')) 
    pickle.dump(p_yy_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][3], s82_str),'wb')) 
    pickle.dump(p_zz_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][4], s82_str),'wb')) 
    pickle.dump(p_ww_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][5], s82_str),'wb')) 
    pickle.dump(corrarr, open(pkl_path+'%s_%scorrarr.pkl' % (var_def, s82_str),'wb')) 
    pickle.dump(corrbinarr, open(pkl_path+'%s_%scorrbinarr.pkl' % (var_def, s82_str),'wb')) 

    return p_el_adj_arr2, p_sp_adj_arr2, p_xx_adj_arr2, p_yy_adj_arr2, p_zz_adj_arr2, p_ww_adj_arr2

def adjust_probabilities_n5(task_dict, plot=False, stripe82=False):

    # Load in the raw probabilities from the GZ2 catalog

    if stripe82:
        p = pyfits.open(gz2_stripe82_data_file)
        s82_str = 'stripe82_'
    else:
        p = pyfits.open(gz2_full_data_file)
        s82_str = ''

    gzdata = p[1].data
    p.close()

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load in the bin sizes from the fits

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[0],bintype))

        edges_redshift = p_var1['REDSHIFT_BIN_EDGES'].data['edges']
        edges_mag = p_var1['MR_BIN_EDGES'].data['edges']
        edges_size = p_var1['R50_KPC_BIN_EDGES'].data['edges']

        p_var1.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    """
    Take each galaxy, find which bin it corresponds to, adjust probability via Bamford+09 method
    """

    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']

    zstep = edges_redshift[1] - edges_redshift[0]
    magstep = edges_mag[1] - edges_mag[0]
    sizestep = edges_size[1] - edges_size[0]

    # Loop over combinations of tasks

    corr_12, corrmasked_12, ratiomasked_12 = determine_baseline_correction(task_dict, 0, 1)
    corr_13, corrmasked_13, ratiomasked_13 = determine_baseline_correction(task_dict, 0, 2)
    corr_14, corrmasked_14, ratiomasked_14 = determine_baseline_correction(task_dict, 0, 3)
    corr_15, corrmasked_15, ratiomasked_15 = determine_baseline_correction(task_dict, 0, 4)
    corr_21, corrmasked_21, ratiomasked_21 = determine_baseline_correction(task_dict, 1, 0)
    corr_23, corrmasked_23, ratiomasked_23 = determine_baseline_correction(task_dict, 1, 2)
    corr_24, corrmasked_24, ratiomasked_24 = determine_baseline_correction(task_dict, 1, 3)
    corr_25, corrmasked_25, ratiomasked_25 = determine_baseline_correction(task_dict, 1, 4)
    corr_32, corrmasked_32, ratiomasked_32 = determine_baseline_correction(task_dict, 2, 1)
    corr_31, corrmasked_31, ratiomasked_31 = determine_baseline_correction(task_dict, 2, 0)
    corr_34, corrmasked_34, ratiomasked_34 = determine_baseline_correction(task_dict, 2, 3)
    corr_35, corrmasked_35, ratiomasked_35 = determine_baseline_correction(task_dict, 2, 4)
    corr_42, corrmasked_42, ratiomasked_42 = determine_baseline_correction(task_dict, 3, 1)
    corr_43, corrmasked_43, ratiomasked_43 = determine_baseline_correction(task_dict, 3, 2)
    corr_41, corrmasked_41, ratiomasked_41 = determine_baseline_correction(task_dict, 3, 0)
    corr_45, corrmasked_45, ratiomasked_45 = determine_baseline_correction(task_dict, 3, 4)
    corr_52, corrmasked_52, ratiomasked_52 = determine_baseline_correction(task_dict, 4, 1)
    corr_53, corrmasked_53, ratiomasked_53 = determine_baseline_correction(task_dict, 4, 2)
    corr_54, corrmasked_54, ratiomasked_54 = determine_baseline_correction(task_dict, 4, 3)
    corr_51, corrmasked_51, ratiomasked_51 = determine_baseline_correction(task_dict, 4, 0)

    p_el_raw = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw = gzdata_withz[task_dict['task_names_wf'][2]]      # Designate the third variable (for now) as 'xx'
    p_yy_raw = gzdata_withz[task_dict['task_names_wf'][3]]      
    p_zz_raw = gzdata_withz[task_dict['task_names_wf'][4]]     
    p_task_counts = gzdata_withz[task_dict['task_name_count']]

    # Start adjusting

    timestart2 = time.time()

    p_el_adj_arr2 = np.zeros_like(p_el_raw)
    p_sp_adj_arr2 = np.zeros_like(p_sp_raw)
    p_xx_adj_arr2 = np.zeros_like(p_xx_raw)
    p_yy_adj_arr2 = np.zeros_like(p_yy_raw)
    p_zz_adj_arr2 = np.zeros_like(p_zz_raw)
    corrarr = np.zeros_like(p_sp_raw)
    corrbinarr = np.zeros((3,len(p_el_raw)))                    # Three bins give the (z,M,r) coordinates of the matching bin

    corrshape = corr_12.shape                                   # Should probably assert that corrections are all the same shape.
    allvotecount = 0
    outbincount = 0
    corrcount = 0
    unknowncount = 0

    # Loop over all galaxies in the GZ catalog with redshifts

    for idx, (el, sp, xx, yy, zz) in enumerate(zip(p_el_raw, p_sp_raw, p_xx_raw, p_yy_raw, p_zz_raw)):

        z2 = redshift[idx]
        m2 = mr[idx]
        s2 = r50_kpc[idx]

        zbin = np.abs(z2 >= edges_redshift).argmin() - 1
        mbin = np.abs(m2 >= edges_mag).argmin() - 1
        sbin = np.abs(s2 >= edges_size).argmin() - 1

        if (zbin < corrshape[0]) and (mbin < corrshape[1]) and (sbin < corrshape[2]):
            applycorr_12 = 0. if corr_12[zbin,mbin,sbin] == unknown_ratio else corr_12[zbin,mbin,sbin]
            applycorr_13 = 0. if corr_13[zbin,mbin,sbin] == unknown_ratio else corr_13[zbin,mbin,sbin]
            applycorr_14 = 0. if corr_14[zbin,mbin,sbin] == unknown_ratio else corr_14[zbin,mbin,sbin]
            applycorr_15 = 0. if corr_15[zbin,mbin,sbin] == unknown_ratio else corr_15[zbin,mbin,sbin]
            applycorr_21 = 0. if corr_21[zbin,mbin,sbin] == unknown_ratio else corr_21[zbin,mbin,sbin]
            applycorr_23 = 0. if corr_23[zbin,mbin,sbin] == unknown_ratio else corr_23[zbin,mbin,sbin]
            applycorr_24 = 0. if corr_24[zbin,mbin,sbin] == unknown_ratio else corr_24[zbin,mbin,sbin]
            applycorr_25 = 0. if corr_25[zbin,mbin,sbin] == unknown_ratio else corr_25[zbin,mbin,sbin]
            applycorr_32 = 0. if corr_32[zbin,mbin,sbin] == unknown_ratio else corr_32[zbin,mbin,sbin]
            applycorr_31 = 0. if corr_31[zbin,mbin,sbin] == unknown_ratio else corr_31[zbin,mbin,sbin]
            applycorr_34 = 0. if corr_34[zbin,mbin,sbin] == unknown_ratio else corr_34[zbin,mbin,sbin]
            applycorr_35 = 0. if corr_35[zbin,mbin,sbin] == unknown_ratio else corr_35[zbin,mbin,sbin]
            applycorr_42 = 0. if corr_42[zbin,mbin,sbin] == unknown_ratio else corr_42[zbin,mbin,sbin]
            applycorr_43 = 0. if corr_43[zbin,mbin,sbin] == unknown_ratio else corr_43[zbin,mbin,sbin]
            applycorr_41 = 0. if corr_41[zbin,mbin,sbin] == unknown_ratio else corr_41[zbin,mbin,sbin]
            applycorr_45 = 0. if corr_45[zbin,mbin,sbin] == unknown_ratio else corr_45[zbin,mbin,sbin]
            applycorr_52 = 0. if corr_52[zbin,mbin,sbin] == unknown_ratio else corr_52[zbin,mbin,sbin]
            applycorr_53 = 0. if corr_53[zbin,mbin,sbin] == unknown_ratio else corr_53[zbin,mbin,sbin]
            applycorr_54 = 0. if corr_54[zbin,mbin,sbin] == unknown_ratio else corr_54[zbin,mbin,sbin]
            applycorr_51 = 0. if corr_51[zbin,mbin,sbin] == unknown_ratio else corr_51[zbin,mbin,sbin]

            allcorrs = np.array((applycorr_12,
                                 applycorr_13,
                                 applycorr_14,
                                 applycorr_15,
                                 applycorr_21,
                                 applycorr_23,
                                 applycorr_24,
                                 applycorr_25,
                                 applycorr_32,
                                 applycorr_31,
                                 applycorr_34,
                                 applycorr_35,
                                 applycorr_52,
                                 applycorr_53,
                                 applycorr_54,
                                 applycorr_51
                                )
                               )

            if np.sum(allcorrs) == 0.:
                unknowncount +=1
                p_el_new = el
                p_sp_new = sp
                p_xx_new = xx
                p_yy_new = yy
                p_zz_new = zz
                corrarr[idx] = unknown_ratio
                corrbinarr[:,idx] = unknown_ratio
            else:
                corrcount += 1
                p_el_new = p_n5_adj(el,sp,xx,yy,zz,
                                    applycorr_12,applycorr_13,applycorr_14,applycorr_15,
                                    task_dict['ratio_type']) if el > 0. else el
                p_sp_new = p_n5_adj(sp,el,xx,yy,zz,
                                    applycorr_21,applycorr_23,applycorr_24,applycorr_25,
                                    task_dict['ratio_type']) if sp > 0. else sp
                p_xx_new = p_n5_adj(xx,el,sp,yy,zz,
                                    applycorr_31,applycorr_32,applycorr_34,applycorr_35,
                                    task_dict['ratio_type']) if xx > 0. else xx
                p_yy_new = p_n5_adj(yy,el,sp,xx,zz,
                                    applycorr_41,applycorr_42,applycorr_43,applycorr_45,
                                    task_dict['ratio_type']) if yy > 0. else yy
                p_zz_new = p_n5_adj(zz,el,sp,xx,yy,
                                    applycorr_51,applycorr_52,applycorr_53,applycorr_54,
                                    task_dict['ratio_type']) if zz > 0. else zz
                corrbinarr[:,idx] = zbin,mbin,sbin
                corrarr[idx] = applycorr_12

            p_el_adj_arr2[idx] = p_el_new
            p_sp_adj_arr2[idx] = p_sp_new
            p_xx_adj_arr2[idx] = p_xx_new
            p_yy_adj_arr2[idx] = p_yy_new
            p_zz_adj_arr2[idx] = p_zz_new
          
        else: # if galaxy does not appear in the correction area
            outbincount += 1
            p_el_adj_arr2[idx] = el
            p_sp_adj_arr2[idx] = sp
            p_xx_adj_arr2[idx] = xx
            p_yy_adj_arr2[idx] = yy
            p_zz_adj_arr2[idx] = zz
            corrbinarr[:,idx] = unknown_ratio
            corrarr[idx] = unknown_ratio


    timeend2 = time.time()
    print ' '
    #print '%7i galaxies had 100 percent of the vote for %s or %s' % (allvotecount,var_str[vartop],var_str[varbot])
    print '%7i galaxies had correction bins larger than volume' % outbincount
    print '%7i galaxies corrected' % corrcount
    print 'Time elapsed to adjust probabilities: %i seconds' % (timeend2 - timestart2)
    print ' '

    p.close()

    pickle.dump(p_el_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][0], s82_str),'wb')) 
    pickle.dump(p_sp_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][1], s82_str),'wb')) 
    pickle.dump(p_xx_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][2], s82_str),'wb')) 
    pickle.dump(p_yy_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][3], s82_str),'wb')) 
    pickle.dump(p_zz_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][4], s82_str),'wb')) 
    pickle.dump(corrarr, open(pkl_path+'%s_%scorrarr.pkl' % (var_def, s82_str),'wb')) 
    pickle.dump(corrbinarr, open(pkl_path+'%s_%scorrbinarr.pkl' % (var_def, s82_str),'wb')) 

    return p_el_adj_arr2, p_sp_adj_arr2, p_xx_adj_arr2, p_yy_adj_arr2, p_zz_adj_arr2

def adjust_probabilities_n4(task_dict, plot=False, stripe82=False):

    # Load in the raw probabilities from the GZ2 catalog

    if stripe82:
        p = pyfits.open(gz2_stripe82_data_file)
        s82_str = 'stripe82_'
    else:
        p = pyfits.open(gz2_full_data_file)
        s82_str = ''

    gzdata = p[1].data
    p.close()

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load in the bin sizes from the fits

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[0],bintype))

        edges_redshift = p_var1['REDSHIFT_BIN_EDGES'].data['edges']
        edges_mag = p_var1['MR_BIN_EDGES'].data['edges']
        edges_size = p_var1['R50_KPC_BIN_EDGES'].data['edges']

        p_var1.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    """
    Take each galaxy, find which bin it corresponds to, adjust probability via Bamford+09 method
    """

    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']

    zstep = edges_redshift[1] - edges_redshift[0]
    magstep = edges_mag[1] - edges_mag[0]
    sizestep = edges_size[1] - edges_size[0]

    # Loop over combinations of tasks

    corr_12, corrmasked_12, ratiomasked_12 = determine_baseline_correction(task_dict, 0, 1)
    corr_13, corrmasked_13, ratiomasked_13 = determine_baseline_correction(task_dict, 0, 2)
    corr_14, corrmasked_14, ratiomasked_14 = determine_baseline_correction(task_dict, 0, 3)
    corr_21, corrmasked_21, ratiomasked_21 = determine_baseline_correction(task_dict, 1, 0)
    corr_23, corrmasked_23, ratiomasked_23 = determine_baseline_correction(task_dict, 1, 2)
    corr_24, corrmasked_24, ratiomasked_24 = determine_baseline_correction(task_dict, 1, 3)
    corr_32, corrmasked_32, ratiomasked_32 = determine_baseline_correction(task_dict, 2, 1)
    corr_31, corrmasked_31, ratiomasked_31 = determine_baseline_correction(task_dict, 2, 0)
    corr_34, corrmasked_34, ratiomasked_34 = determine_baseline_correction(task_dict, 2, 3)
    corr_42, corrmasked_42, ratiomasked_42 = determine_baseline_correction(task_dict, 3, 1)
    corr_43, corrmasked_43, ratiomasked_43 = determine_baseline_correction(task_dict, 3, 2)
    corr_41, corrmasked_41, ratiomasked_41 = determine_baseline_correction(task_dict, 3, 0)

    p_el_raw = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw = gzdata_withz[task_dict['task_names_wf'][2]]      # Designate the third variable (for now) as 'xx'
    p_yy_raw = gzdata_withz[task_dict['task_names_wf'][3]]      # Designate the third variable (for now) as 'xx'
    p_task_counts = gzdata_withz[task_dict['task_name_count']]

    # Start adjusting

    timestart2 = time.time()

    p_el_adj_arr2 = np.zeros_like(p_el_raw)
    p_sp_adj_arr2 = np.zeros_like(p_sp_raw)
    p_xx_adj_arr2 = np.zeros_like(p_xx_raw)
    p_yy_adj_arr2 = np.zeros_like(p_yy_raw)
    corrarr = np.zeros_like(p_sp_raw)
    corrbinarr = np.zeros((3,len(p_el_raw)))                    # Three bins give the (z,M,r) coordinates of the matching bin

    corrshape = corr_12.shape                                   # Should probably assert that corrections are all the same shape.
    allvotecount = 0
    outbincount = 0
    corrcount = 0
    unknowncount = 0

    # Loop over all galaxies in the GZ catalog with redshifts

    for idx, (el, sp, xx, yy) in enumerate(zip(p_el_raw, p_sp_raw, p_xx_raw, p_yy_raw)):

        z2 = redshift[idx]
        m2 = mr[idx]
        s2 = r50_kpc[idx]

        zbin = np.abs(z2 >= edges_redshift).argmin() - 1
        mbin = np.abs(m2 >= edges_mag).argmin() - 1
        sbin = np.abs(s2 >= edges_size).argmin() - 1

        if (zbin < corrshape[0]) and (mbin < corrshape[1]) and (sbin < corrshape[2]):
            applycorr_12 = 0. if corr_12[zbin,mbin,sbin] == unknown_ratio else corr_12[zbin,mbin,sbin]
            applycorr_13 = 0. if corr_13[zbin,mbin,sbin] == unknown_ratio else corr_13[zbin,mbin,sbin]
            applycorr_14 = 0. if corr_14[zbin,mbin,sbin] == unknown_ratio else corr_14[zbin,mbin,sbin]
            applycorr_21 = 0. if corr_21[zbin,mbin,sbin] == unknown_ratio else corr_21[zbin,mbin,sbin]
            applycorr_23 = 0. if corr_23[zbin,mbin,sbin] == unknown_ratio else corr_23[zbin,mbin,sbin]
            applycorr_24 = 0. if corr_24[zbin,mbin,sbin] == unknown_ratio else corr_24[zbin,mbin,sbin]
            applycorr_32 = 0. if corr_32[zbin,mbin,sbin] == unknown_ratio else corr_32[zbin,mbin,sbin]
            applycorr_31 = 0. if corr_31[zbin,mbin,sbin] == unknown_ratio else corr_31[zbin,mbin,sbin]
            applycorr_34 = 0. if corr_34[zbin,mbin,sbin] == unknown_ratio else corr_34[zbin,mbin,sbin]
            applycorr_42 = 0. if corr_42[zbin,mbin,sbin] == unknown_ratio else corr_42[zbin,mbin,sbin]
            applycorr_43 = 0. if corr_43[zbin,mbin,sbin] == unknown_ratio else corr_43[zbin,mbin,sbin]
            applycorr_41 = 0. if corr_41[zbin,mbin,sbin] == unknown_ratio else corr_41[zbin,mbin,sbin]

            allcorrs = np.array((applycorr_12,
                                 applycorr_13,
                                 applycorr_14,
                                 applycorr_21,
                                 applycorr_23,
                                 applycorr_24,
                                 applycorr_32,
                                 applycorr_31,
                                 applycorr_34,
                                 applycorr_42,
                                 applycorr_43,
                                 applycorr_41,
                                )
                               )

            if np.sum(allcorrs) == 0.:
                unknowncount +=1
                p_el_new = el
                p_sp_new = sp
                p_xx_new = xx
                p_yy_new = yy
                corrarr[idx] = unknown_ratio
                corrbinarr[:,idx] = unknown_ratio
            else:
                corrcount += 1
                p_el_new = p_n4_adj(el,sp,xx,yy,
                                    applycorr_12,applycorr_13,applycorr_14,
                                    task_dict['ratio_type']) if el > 0. else el
                p_sp_new = p_n4_adj(sp,el,xx,yy,
                                    applycorr_21,applycorr_23,applycorr_24,
                                    task_dict['ratio_type']) if sp > 0. else sp
                p_xx_new = p_n4_adj(xx,el,sp,yy,
                                    applycorr_31,applycorr_32,applycorr_34,
                                    task_dict['ratio_type']) if xx > 0. else xx
                p_yy_new = p_n4_adj(yy,el,sp,xx,
                                    applycorr_41,applycorr_42,applycorr_43,
                                    task_dict['ratio_type']) if yy > 0. else yy
                corrbinarr[:,idx] = zbin,mbin,sbin
                corrarr[idx] = applycorr_12

            p_el_adj_arr2[idx] = p_el_new
            p_sp_adj_arr2[idx] = p_sp_new
            p_xx_adj_arr2[idx] = p_xx_new
            p_yy_adj_arr2[idx] = p_yy_new
          
        else: # if galaxy does not appear in the correction area
            outbincount += 1
            p_el_adj_arr2[idx] = el
            p_sp_adj_arr2[idx] = sp
            p_xx_adj_arr2[idx] = xx
            p_yy_adj_arr2[idx] = yy
            corrbinarr[:,idx] = unknown_ratio
            corrarr[idx] = unknown_ratio


    timeend2 = time.time()
    print ' '
    #print '%7i galaxies had 100 percent of the vote for %s or %s' % (allvotecount,var_str[vartop],var_str[varbot])
    print '%7i galaxies had correction bins larger than volume' % outbincount
    print '%7i galaxies corrected' % corrcount
    print 'Time elapsed to adjust probabilities: %i seconds' % (timeend2 - timestart2)
    print ' '

    p.close()

    pickle.dump(p_el_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][0], s82_str),'wb')) 
    pickle.dump(p_sp_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][1], s82_str),'wb')) 
    pickle.dump(p_xx_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][2], s82_str),'wb')) 
    pickle.dump(p_yy_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][3], s82_str),'wb')) 
    pickle.dump(corrarr, open(pkl_path+'%s_%scorrarr.pkl' % (var_def, s82_str),'wb')) 
    pickle.dump(corrbinarr, open(pkl_path+'%s_%scorrbinarr.pkl' % (var_def, s82_str),'wb')) 

    return p_el_adj_arr2, p_sp_adj_arr2, p_xx_adj_arr2, p_yy_adj_arr2

def adjust_probabilities_n3(task_dict, plot=False, stripe82=False):

    # Load in the raw probabilities from the GZ2 catalog

    if stripe82:
        p = pyfits.open(gz2_stripe82_data_file)
        s82_str = 'stripe82_'
    else:
        p = pyfits.open(gz2_full_data_file)
        s82_str = ''

    gzdata = p[1].data
    p.close()

    var_def = task_dict['var_def']
    var_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load in the bin sizes from the fits

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var_str[0],bintype))

        edges_redshift = p_var1['REDSHIFT_BIN_EDGES'].data['edges']
        edges_mag = p_var1['MR_BIN_EDGES'].data['edges']
        edges_size = p_var1['R50_KPC_BIN_EDGES'].data['edges']

        p_var1.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    """
    Take each galaxy, find which bin it corresponds to, adjust probability via Bamford+09 method
    """

    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']

    zstep = edges_redshift[1] - edges_redshift[0]
    magstep = edges_mag[1] - edges_mag[0]
    sizestep = edges_size[1] - edges_size[0]

    # Need six different corrections: K 12, 23, 13, 21, 32, 31

    varsarr = ((1,2),(2,3),(1,3),(2,1),(3,2),(3,1))

    # Loop over combinations of tasks

    corr_12, corrmasked_12, ratiomasked_12 = determine_baseline_correction(task_dict, 0, 1)
    corr_13, corrmasked_13, ratiomasked_13 = determine_baseline_correction(task_dict, 0, 2)
    corr_21, corrmasked_21, ratiomasked_21 = determine_baseline_correction(task_dict, 1, 0)
    corr_23, corrmasked_23, ratiomasked_23 = determine_baseline_correction(task_dict, 1, 2)
    corr_31, corrmasked_31, ratiomasked_31 = determine_baseline_correction(task_dict, 2, 0)
    corr_32, corrmasked_32, ratiomasked_32 = determine_baseline_correction(task_dict, 2, 1)

    p_el_raw = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw = gzdata_withz[task_dict['task_names_wf'][2]]      # Designate the third variable (for now) as 'xx'
    p_task_counts = gzdata_withz[task_dict['task_name_count']]

    # Start adjusting

    timestart2 = time.time()

    p_el_adj_arr2 = np.zeros_like(p_el_raw)
    p_sp_adj_arr2 = np.zeros_like(p_sp_raw)
    p_xx_adj_arr2 = np.zeros_like(p_xx_raw)
    corrarr = np.zeros_like(p_sp_raw)
    corrbinarr = np.zeros((3,len(p_el_raw)))                    # Three bins give the (z,M,r) coordinates of the matching bin

    corrshape = corr_12.shape                                   # Should probably assert that corrections are all the same shape.
    allvotecount = 0
    outbincount = 0
    corrcount = 0
    unknowncount = 0

    # Loop over all galaxies in the GZ catalog with redshifts

    for idx, (el, sp, xx) in enumerate(zip(p_el_raw, p_sp_raw, p_xx_raw)):

        z2 = redshift[idx]
        m2 = mr[idx]
        s2 = r50_kpc[idx]

        zbin = np.abs(z2 >= edges_redshift).argmin() - 1
        mbin = np.abs(m2 >= edges_mag).argmin() - 1
        sbin = np.abs(s2 >= edges_size).argmin() - 1

        if (zbin < corrshape[0]) and (mbin < corrshape[1]) and (sbin < corrshape[2]):
            applycorr_12 = 0. if corr_12[zbin,mbin,sbin] == unknown_ratio else corr_12[zbin,mbin,sbin]
            applycorr_13 = 0. if corr_13[zbin,mbin,sbin] == unknown_ratio else corr_13[zbin,mbin,sbin]
            applycorr_21 = 0. if corr_21[zbin,mbin,sbin] == unknown_ratio else corr_21[zbin,mbin,sbin]
            applycorr_23 = 0. if corr_23[zbin,mbin,sbin] == unknown_ratio else corr_23[zbin,mbin,sbin]
            applycorr_31 = 0. if corr_31[zbin,mbin,sbin] == unknown_ratio else corr_31[zbin,mbin,sbin]
            applycorr_32 = 0. if corr_32[zbin,mbin,sbin] == unknown_ratio else corr_32[zbin,mbin,sbin]
            allcorrs = np.array((applycorr_12,
                                 applycorr_13,
                                 applycorr_21,
                                 applycorr_23,
                                 applycorr_31,
                                 applycorr_32
                                )
                               )

            if np.sum(allcorrs) == 0.:
                unknowncount +=1
                p_el_new = el
                p_sp_new = sp
                p_xx_new = xx
                corrarr[idx] = unknown_ratio
                corrbinarr[:,idx] = unknown_ratio
            else:
                corrcount += 1
                p_el_new = p_n3_adj(el,sp,xx,applycorr_12,applycorr_13,task_dict['ratio_type']) if el > 0. else el
                p_sp_new = p_n3_adj(sp,el,xx,applycorr_21,applycorr_23,task_dict['ratio_type']) if sp > 0. else sp
                p_xx_new = p_n3_adj(xx,el,sp,applycorr_31,applycorr_32,task_dict['ratio_type']) if xx > 0. else xx
                corrbinarr[:,idx] = zbin,mbin,sbin
                corrarr[idx] = applycorr_12

            p_el_adj_arr2[idx] = p_el_new
            p_sp_adj_arr2[idx] = p_sp_new
            p_xx_adj_arr2[idx] = p_xx_new
          
        else: # if galaxy does not appear in the correction area
            outbincount += 1
            p_el_adj_arr2[idx] = el
            p_sp_adj_arr2[idx] = sp
            p_xx_adj_arr2[idx] = xx
            corrbinarr[:,idx] = unknown_ratio
            corrarr[idx] = unknown_ratio


    timeend2 = time.time()
    print ' '
    #print '%7i galaxies had 100 percent of the vote for %s or %s' % (allvotecount,var_str[vartop],var_str[varbot])
    print '%7i galaxies had correction bins larger than volume' % outbincount
    print '%7i galaxies corrected' % corrcount
    print 'Time elapsed to adjust probabilities: %i seconds' % (timeend2 - timestart2)
    print ' '

    if plot:
        
        # Downsample the galaxies to reduce plotting time

        nsamp = 1000
        goodind = np.arange(len(corrarr))[np.arange(len(corrarr))[corrarr > unknown_ratio]]
        randind = goodind[np.random.random_integers(0,len(goodind),nsamp)]
        corrhi = 1.0
        corrlo = -0.5
        rand_correction = (corrhi-corrlo) * np.random.random(nsamp) + corrlo
        constant_correction = np.zeros(nsamp) + 0.5

        # Set plot parameters

        fig = plt.figure(8)
        fig.clf()

        ax_el = fig.add_subplot(221)
        ax_el.scatter(p_el_raw[randind],p_el_adj_arr2[randind], marker='.')
        ax_el.set_xlabel(r'$p_{%s}$ raw' % var_str[0])
        ax_el.set_ylabel(r'$p_{%s}$ adjusted' % var_str[0])

        ax_sp = fig.add_subplot(222)
        ax_sp.scatter(p_sp_raw[randind],p_sp_adj_arr2[randind], marker='.')
        ax_sp.set_xlabel(r'$p_{%s}$ raw' % var_str[varbot])
        ax_sp.set_ylabel(r'$p_{%s}$ adjusted' % var_str[varbot])

        ax_el_sim = fig.add_subplot(223)
        ax_el_sim.scatter(p_el_raw[randind],p_el_adj(p_el_raw[randind],p_sp_raw[randind],rand_correction,task_dict['ratio_type']), 
                          marker='.',color='g')
        ax_el_sim.scatter(p_el_raw[randind],p_el_adj(p_el_raw[randind],p_sp_raw[randind],constant_correction,task_dict['ratio_type']), 
                          marker='.',color='r')
        ax_el_sim.set_xlabel(r'$p_{%s}$ raw' % var_str[0])
        ax_el_sim.set_ylabel(r'$p_{%s}$ adjusted' % var_str[0])

        plt.legend(('Random correction', 'Constant '+r'$C = 0.5$'), 'upper left', shadow=True, fancybox=True)

        ax_sp_sim = fig.add_subplot(224)
        ax_sp_sim.scatter(p_sp_raw[randind],p_sp_adj(p_el_raw[randind],p_sp_raw[randind],rand_correction,task_dict['ratio_type']), 
                          marker='.',color='g')
        ax_sp_sim.scatter(p_sp_raw[randind],p_sp_adj(p_el_raw[randind],p_sp_raw[randind],constant_correction,task_dict['ratio_type']), 
                          marker='.',color='r')
        ax_sp_sim.set_xlabel(r'$p_{%s}$ raw' % var_str[varbot])
        ax_sp_sim.set_ylabel(r'$p_{%s}$ adjusted' % var_str[varbot])

        plt.legend(('Random correction', 'Constant '+r'$C = 0.5$'), 'upper left', shadow=True, fancybox=True)

    p.close()

    pickle.dump(p_el_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][0], s82_str),'wb')) 
    pickle.dump(p_sp_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][1], s82_str),'wb')) 
    pickle.dump(p_xx_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,task_dict['var_str'][2], s82_str),'wb')) 
    pickle.dump(corrarr, open(pkl_path+'%s_%scorrarr.pkl' % (var_def, s82_str),'wb')) 
    pickle.dump(corrbinarr, open(pkl_path+'%s_%scorrbinarr.pkl' % (var_def, s82_str),'wb')) 

    return p_el_adj_arr2, p_sp_adj_arr2, p_xx_adj_arr2

def adjust_probabilities(task_dict, plot=False, stripe82=False):

    # Load in the raw probabilities from the GZ2 catalog

    if stripe82:
        p = pyfits.open(gz2_stripe82_data_file)
        s82_str = 'stripe82_'
    else:
        p = pyfits.open(gz2_full_data_file)
        s82_str = ''

    gzdata = p[1].data
    p.close()

    correction, correction_masked, ratio_masked = determine_baseline_correction(task_dict)

    var_def = task_dict['var_def']
    var1_str,var2_str = task_dict['var_str']
    bintype = task_dict['bintype']

    # Load in the bin sizes from the fits

    if bintype is 'counts':
        p_var1 = pyfits.open(fits_path+'%s_%s_binned_%s.fits' % (var_def,var1_str,bintype))

        edges_redshift = p_var1['REDSHIFT_BIN_EDGES'].data['edges']
        edges_mag = p_var1['MR_BIN_EDGES'].data['edges']
        edges_size = p_var1['R50_KPC_BIN_EDGES'].data['edges']

        p_var1.close()

    if bintype is 'rawlikelihood':
        #p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        #d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        p_allvar = pyfits.open(fits_path+'%s_idlbinned.fits' % var_def)
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        p_allvar.close()

    """
    Take each galaxy, find which bin it corresponds to, adjust probability via Bamford+09 method
    """

    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']

    zstep = edges_redshift[1] - edges_redshift[0]
    magstep = edges_mag[1] - edges_mag[0]
    sizestep = edges_size[1] - edges_size[0]

    p_el_raw = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw = gzdata_withz[task_dict['task_names_wf'][1]]
    p_task_counts = gzdata_withz[task_dict['task_name_count']]

    # Start adjusting

    timestart2 = time.time()

    p_el_adj_arr2 = np.zeros_like(p_el_raw)
    p_sp_adj_arr2 = np.zeros_like(p_sp_raw)
    corrarr = np.zeros_like(p_sp_raw)
    corrbinarr = np.zeros((3,len(p_el_raw)))

    corrshape = correction.shape
    allvotecount = 0
    outbincount = 0
    corrcount = 0
    unknowncount = 0

    # Loop over all galaxies in the GZ catalog with redshifts

    for idx, (el, sp) in enumerate(zip(p_el_raw, p_sp_raw)):

        z2 = redshift[idx]
        m2 = mr[idx]
        s2 = r50_kpc[idx]

        zbin = np.abs(z2 >= edges_redshift).argmin() - 1
        mbin = np.abs(m2 >= edges_mag).argmin() - 1
        sbin = np.abs(s2 >= edges_size).argmin() - 1

        if (zbin < corrshape[0]) and (mbin < corrshape[1]) and (sbin < corrshape[2]):
            corrbin2 = correction[zbin,mbin,sbin]

            if el == 0. or sp == 0.:
                allvotecount += 1
                p_sp_new = sp
                p_el_new = el
                corrarr[idx] = 0.
                corrbinarr[:,idx] = no_correction_applied
            elif corrbin2 == unknown_ratio:
                unknowncount +=1
                p_el_new = el
                p_sp_new = sp
                corrarr[idx] = unknown_ratio
                corrbinarr[:,idx] = unknown_ratio
            else:
                corrcount += 1
                p_el_new = p_el_adj(el,sp,corrbin2,task_dict['ratio_type'])
                p_sp_new = p_sp_adj(el,sp,corrbin2,task_dict['ratio_type'])
                corrbinarr[:,idx] = zbin,mbin,sbin
                corrarr[idx] = corrbin2

            p_el_adj_arr2[idx] = p_el_new
            p_sp_adj_arr2[idx] = p_sp_new
          
        else:
            outbincount += 1
            p_el_adj_arr2[idx] = el
            p_sp_adj_arr2[idx] = sp
            corrbinarr[:,idx] = unknown_ratio
            corrarr[idx] = unknown_ratio


    timeend2 = time.time()
    print ' '
    print '%7i galaxies had 100 percent of the vote for %s or %s' % (allvotecount,var1_str,var2_str)
    print '%7i galaxies had correction bins larger than volume' % outbincount
    print '%7i galaxies corrected' % corrcount
    print 'Time elapsed to adjust probabilities: %i seconds' % (timeend2 - timestart2)
    print ' '

    if plot:
        
        # Downsample the galaxies to reduce plotting time

        nsamp = 1000
        goodind = np.arange(len(corrarr))[np.arange(len(corrarr))[corrarr > unknown_ratio]]
        randind = goodind[np.random.random_integers(0,len(goodind),nsamp)]
        corrhi = 1.0
        corrlo = -0.5
        rand_correction = (corrhi-corrlo) * np.random.random(nsamp) + corrlo
        constant_correction = np.zeros(nsamp) + 0.5

        # Set plot parameters

        fig = plt.figure(8)
        fig.clf()

        ax_el = fig.add_subplot(221)
        ax_el.scatter(p_el_raw[randind],p_el_adj_arr2[randind], marker='.')
        ax_el.set_xlabel(r'$p_{%s}$ raw' % var1_str)
        ax_el.set_ylabel(r'$p_{%s}$ adjusted' % var1_str)

        ax_sp = fig.add_subplot(222)
        ax_sp.scatter(p_sp_raw[randind],p_sp_adj_arr2[randind], marker='.')
        ax_sp.set_xlabel(r'$p_{%s}$ raw' % var2_str)
        ax_sp.set_ylabel(r'$p_{%s}$ adjusted' % var2_str)

        ax_el_sim = fig.add_subplot(223)
        ax_el_sim.scatter(p_el_raw[randind],p_el_adj(p_el_raw[randind],p_sp_raw[randind],rand_correction,task_dict['ratio_type']), 
                          marker='.',color='g')
        ax_el_sim.scatter(p_el_raw[randind],p_el_adj(p_el_raw[randind],p_sp_raw[randind],constant_correction,task_dict['ratio_type']), 
                          marker='.',color='r')
        ax_el_sim.set_xlabel(r'$p_{%s}$ raw' % var1_str)
        ax_el_sim.set_ylabel(r'$p_{%s}$ adjusted' % var1_str)

        plt.legend(('Random correction', 'Constant '+r'$C = 0.5$'), 'upper left', shadow=True, fancybox=True)

        ax_sp_sim = fig.add_subplot(224)
        ax_sp_sim.scatter(p_sp_raw[randind],p_sp_adj(p_el_raw[randind],p_sp_raw[randind],rand_correction,task_dict['ratio_type']), 
                          marker='.',color='g')
        ax_sp_sim.scatter(p_sp_raw[randind],p_sp_adj(p_el_raw[randind],p_sp_raw[randind],constant_correction,task_dict['ratio_type']), 
                          marker='.',color='r')
        ax_sp_sim.set_xlabel(r'$p_{%s}$ raw' % var2_str)
        ax_sp_sim.set_ylabel(r'$p_{%s}$ adjusted' % var2_str)

        plt.legend(('Random correction', 'Constant '+r'$C = 0.5$'), 'upper left', shadow=True, fancybox=True)

    p.close()

    pickle.dump(p_el_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var1_str, s82_str),'wb')) 
    pickle.dump(p_sp_adj_arr2, open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var2_str, s82_str),'wb')) 
    pickle.dump(corrarr, open(pkl_path+'%s_%scorrarr.pkl' % (var_def, s82_str),'wb')) 
    pickle.dump(corrbinarr, open(pkl_path+'%s_%scorrbinarr.pkl' % (var_def, s82_str),'wb')) 

    #pickle.dump(p_el_raw, open(pkl_path+'%s_%s_%sraw.pkl' % (var_def,var1_str,s82_str),'wb')) 
    #pickle.dump(p_sp_raw, open(pkl_path+'%s_%s_%sraw.pkl' % (var_def,var2_str,s82_str),'wb')) 
    #pickle.dump(p_task_counts, open(pkl_path+'%s_%stask_counts.pkl' % (var_def,s82_str),'wb')) 
    #pickle.dump(redshift, open(pkl_path+'%s_%sredshift.pkl' % (var_def,s82_str),'wb')) 

    return None

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

def plot_type_fractions_n6(task_dict, zlo = 0.01, zhi=0.085, stripe82=False):

    # Load data 

    if stripe82:
        s82_str = 'stripe82_'
        p = pyfits.open(gz2_stripe82_data_file)
    else:
        s82_str = ''
        p = pyfits.open(gz2_full_data_file)

    gzdata = p[1].data
    p.close()
    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']
    p_el_raw_values = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw_values = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw_values = gzdata_withz[task_dict['task_names_wf'][2]]
    p_yy_raw_values = gzdata_withz[task_dict['task_names_wf'][3]]
    p_zz_raw_values = gzdata_withz[task_dict['task_names_wf'][4]]
    p_ww_raw_values = gzdata_withz[task_dict['task_names_wf'][5]]
    task_counts = gzdata_withz[task_dict['task_name_count']]

    zwidth = 0.02
    zplotbins = np.arange(min(edges_redshift[0],zlo),edges_redshift[-1],zwidth)

    n = len(task_dict['var_str'])
    var_def = task_dict['var_def']
    var1_str,var2_str,var3_str,var4_str,var5_str,var6_str = task_dict['var_str']

    p_el_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var1_str,s82_str),'rb')) 
    p_sp_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var2_str,s82_str),'rb')) 
    p_xx_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var3_str,s82_str),'rb')) 
    p_yy_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var4_str,s82_str),'rb')) 
    p_zz_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var5_str,s82_str),'rb')) 
    p_ww_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var6_str,s82_str),'rb')) 
    corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (var_def,s82_str),'rb')) 

    # Only plot galaxies that had a correction available

    corrgood = corrarr > unknown_ratio
    elgood = np.isfinite(p_el_adj_values)
    spgood = np.isfinite(p_sp_adj_values)
    xxgood = np.isfinite(p_xx_adj_values)
    yygood = np.isfinite(p_yy_adj_values)
    zzgood = np.isfinite(p_zz_adj_values)
    wwgood = np.isfinite(p_ww_adj_values)
    goodinds = corrgood & elgood & spgood & xxgood & yygood & zzgood & wwgood
    p_el_adj_values = p_el_adj_values[goodinds]
    p_sp_adj_values = p_sp_adj_values[goodinds]
    p_xx_adj_values = p_xx_adj_values[goodinds]
    p_yy_adj_values = p_yy_adj_values[goodinds]
    p_zz_adj_values = p_zz_adj_values[goodinds]
    p_ww_adj_values = p_ww_adj_values[goodinds]
    task_counts_adj = task_counts[goodinds]

    # Empty arrays for the binned type fractions

    p_el_raw_typefrac = np.zeros_like(zplotbins)
    p_sp_raw_typefrac = np.zeros_like(zplotbins)
    p_el_adj_typefrac = np.zeros_like(zplotbins)
    p_sp_adj_typefrac = np.zeros_like(zplotbins)
    p_el_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_el_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_raw_typefrac = np.zeros_like(zplotbins)
    p_xx_adj_typefrac = np.zeros_like(zplotbins)
    p_xx_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_yy_raw_typefrac = np.zeros_like(zplotbins)
    p_yy_adj_typefrac = np.zeros_like(zplotbins)
    p_yy_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_yy_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_zz_raw_typefrac = np.zeros_like(zplotbins)
    p_zz_adj_typefrac = np.zeros_like(zplotbins)
    p_zz_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_zz_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_ww_raw_typefrac = np.zeros_like(zplotbins)
    p_ww_adj_typefrac = np.zeros_like(zplotbins)
    p_ww_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_ww_adj_typefrac_maglim = np.zeros_like(zplotbins)

    """
    Method 1: take the mean of the raw likelihoods for each galaxy per redshift bin
    Method 2: take the mean of each likelihood weighted by the total number of votes

    Results are within 5% of each other (fraction of < 0.01) at all redshifts.
    """

    # Create a second, magnitude-limited sample
    
    maglimval = appmag_lim - cosmology.dmod_flat(zhi)
    
    mr_raw = mr
    redshift_raw = redshift
    maglim_raw = mr_raw < maglimval
    task_counts_maglim_raw = task_counts[maglim_raw]
    p_el_raw_values_maglim = p_el_raw_values[maglim_raw]
    p_sp_raw_values_maglim = p_sp_raw_values[maglim_raw]
    p_xx_raw_values_maglim = p_xx_raw_values[maglim_raw]
    p_yy_raw_values_maglim = p_yy_raw_values[maglim_raw]
    p_zz_raw_values_maglim = p_zz_raw_values[maglim_raw]
    p_ww_raw_values_maglim = p_ww_raw_values[maglim_raw]

    mr_adj = mr[goodinds]
    redshift_adj = redshift[goodinds]
    maglim_adj = mr_adj < maglimval
    task_counts_maglim_adj = task_counts_adj[maglim_adj]
    p_el_adj_values_maglim = p_el_adj_values[maglim_adj]
    p_sp_adj_values_maglim = p_sp_adj_values[maglim_adj]
    p_xx_adj_values_maglim = p_xx_adj_values[maglim_adj]
    p_yy_adj_values_maglim = p_yy_adj_values[maglim_adj]
    p_zz_adj_values_maglim = p_zz_adj_values[maglim_adj]
    p_ww_adj_values_maglim = p_ww_adj_values[maglim_adj]

    # Loop over redshift bin and find the type fractions at each slice

    for idx,zbin in enumerate(zplotbins):
        gals_in_bin_raw = (redshift_raw >= zbin) & (redshift_raw < (zbin+zwidth))
        weights_raw = task_counts[gals_in_bin_raw]
        p_el_raw_typefrac[idx] = np.mean(p_el_raw_values[gals_in_bin_raw])
        p_sp_raw_typefrac[idx] = np.mean(p_sp_raw_values[gals_in_bin_raw])
        p_xx_raw_typefrac[idx] = np.mean(p_xx_raw_values[gals_in_bin_raw])
        p_yy_raw_typefrac[idx] = np.mean(p_yy_raw_values[gals_in_bin_raw])
        p_zz_raw_typefrac[idx] = np.mean(p_zz_raw_values[gals_in_bin_raw])
        p_ww_raw_typefrac[idx] = np.mean(p_ww_raw_values[gals_in_bin_raw])
        p_el_raw_typefrac_err[idx] = tsem(p_el_raw_values[gals_in_bin_raw])
        p_sp_raw_typefrac_err[idx] = tsem(p_sp_raw_values[gals_in_bin_raw])
        p_xx_raw_typefrac_err[idx] = tsem(p_xx_raw_values[gals_in_bin_raw])
        p_yy_raw_typefrac_err[idx] = tsem(p_yy_raw_values[gals_in_bin_raw])
        p_zz_raw_typefrac_err[idx] = tsem(p_zz_raw_values[gals_in_bin_raw])
        p_ww_raw_typefrac_err[idx] = tsem(p_ww_raw_values[gals_in_bin_raw])

        gals_in_bin_adj = (redshift_adj >= zbin) & (redshift_adj < (zbin+zwidth))
        weights_adj = task_counts_adj[gals_in_bin_adj]
        p_el_adj_typefrac[idx] = np.mean(p_el_adj_values[gals_in_bin_adj])
        p_sp_adj_typefrac[idx] = np.mean(p_sp_adj_values[gals_in_bin_adj])
        p_xx_adj_typefrac[idx] = np.mean(p_xx_adj_values[gals_in_bin_adj])
        p_yy_adj_typefrac[idx] = np.mean(p_yy_adj_values[gals_in_bin_adj])
        p_zz_adj_typefrac[idx] = np.mean(p_zz_adj_values[gals_in_bin_adj])
        p_ww_adj_typefrac[idx] = np.mean(p_ww_adj_values[gals_in_bin_adj])
        p_el_adj_typefrac_err[idx] = tsem(p_el_adj_values[gals_in_bin_adj])
        p_sp_adj_typefrac_err[idx] = tsem(p_sp_adj_values[gals_in_bin_adj])
        p_xx_adj_typefrac_err[idx] = tsem(p_xx_adj_values[gals_in_bin_adj])
        p_yy_adj_typefrac_err[idx] = tsem(p_yy_adj_values[gals_in_bin_adj])
        p_zz_adj_typefrac_err[idx] = tsem(p_zz_adj_values[gals_in_bin_adj])
        p_ww_adj_typefrac_err[idx] = tsem(p_ww_adj_values[gals_in_bin_adj])

        gals_in_bin_maglim_raw = (redshift_raw[maglim_raw] >= zbin) & (redshift_raw[maglim_raw] < (zbin+zwidth))
        weights_maglim_raw = task_counts_maglim_raw[gals_in_bin_maglim_raw]
        p_el_raw_typefrac_maglim[idx] = np.sum(p_el_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_sp_raw_typefrac_maglim[idx] = np.sum(p_sp_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_xx_raw_typefrac_maglim[idx] = np.sum(p_xx_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_yy_raw_typefrac_maglim[idx] = np.sum(p_yy_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_zz_raw_typefrac_maglim[idx] = np.sum(p_zz_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_ww_raw_typefrac_maglim[idx] = np.sum(p_ww_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)

        gals_in_bin_maglim_adj = (redshift_adj[maglim_adj] >= zbin) & (redshift_adj[maglim_adj] < (zbin+zwidth))
        weights_maglim_adj = task_counts_maglim_adj[gals_in_bin_maglim_adj]
        p_el_adj_typefrac_maglim[idx] = np.sum(p_el_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_sp_adj_typefrac_maglim[idx] = np.sum(p_sp_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_xx_adj_typefrac_maglim[idx] = np.sum(p_xx_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_yy_adj_typefrac_maglim[idx] = np.sum(p_yy_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_zz_adj_typefrac_maglim[idx] = np.sum(p_zz_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_ww_adj_typefrac_maglim[idx] = np.sum(p_ww_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)

    pickle.dump(p_el_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 

    pickle.dump(p_xx_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_yy_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_zz_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var5_str,s82_str),'wb')) 
    pickle.dump(p_zz_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var5_str,s82_str),'wb')) 
    pickle.dump(p_zz_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var5_str,s82_str),'wb')) 
    pickle.dump(p_zz_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var5_str,s82_str),'wb')) 

    pickle.dump(p_ww_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var6_str,s82_str),'wb')) 
    pickle.dump(p_ww_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var6_str,s82_str),'wb')) 
    pickle.dump(p_ww_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var6_str,s82_str),'wb')) 
    pickle.dump(p_ww_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var6_str,s82_str),'wb')) 

    # Plot results

    fig = plt.figure(12, (12,7))
    fig.clf()
    ax1 = fig.add_subplot(121)
    font = FontProperties()

    ax1.plot(zplotbins, p_el_raw_typefrac, color='r', linestyle='-' ,linewidth=2)
    ax1.plot(zplotbins, p_sp_raw_typefrac, color='b', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_xx_raw_typefrac, color='m', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_yy_raw_typefrac, color='g', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_zz_raw_typefrac, color='c', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_ww_raw_typefrac, color='y', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_el_adj_typefrac, color='r', linestyle='-' ,linewidth=4)
    ax1.plot(zplotbins, p_sp_adj_typefrac, color='b', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_xx_adj_typefrac, color='m', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_yy_adj_typefrac, color='g', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_zz_adj_typefrac, color='c', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_ww_adj_typefrac, color='y', linestyle='--',linewidth=4)
    ax1.errorbar(zplotbins, p_el_raw_typefrac, yerr = p_el_raw_typefrac_err, color='r', linestyle='-' ,linewidth=2)
    ax1.errorbar(zplotbins, p_sp_raw_typefrac, yerr = p_sp_raw_typefrac_err, color='b', linestyle='--',linewidth=2)
    ax1.errorbar(zplotbins, p_xx_raw_typefrac, yerr = p_xx_raw_typefrac_err, color='m', linestyle='--',linewidth=2)
    ax1.errorbar(zplotbins, p_yy_raw_typefrac, yerr = p_yy_raw_typefrac_err, color='g', linestyle='--',linewidth=2)
    ax1.errorbar(zplotbins, p_zz_raw_typefrac, yerr = p_zz_raw_typefrac_err, color='c', linestyle='--',linewidth=2)
    ax1.errorbar(zplotbins, p_ww_raw_typefrac, yerr = p_ww_raw_typefrac_err, color='y', linestyle='--',linewidth=2)
    ax1.errorbar(zplotbins, p_el_adj_typefrac, yerr = p_el_adj_typefrac_err, color='r', linestyle='-' ,linewidth=4)
    ax1.errorbar(zplotbins, p_sp_adj_typefrac, yerr = p_sp_adj_typefrac_err, color='b', linestyle='--',linewidth=4)
    ax1.errorbar(zplotbins, p_xx_adj_typefrac, yerr = p_xx_adj_typefrac_err, color='m', linestyle='--',linewidth=4)
    ax1.errorbar(zplotbins, p_yy_adj_typefrac, yerr = p_yy_adj_typefrac_err, color='g', linestyle='--',linewidth=4)
    ax1.errorbar(zplotbins, p_zz_adj_typefrac, yerr = p_zz_adj_typefrac_err, color='c', linestyle='--',linewidth=4)
    ax1.errorbar(zplotbins, p_ww_adj_typefrac, yerr = p_ww_adj_typefrac_err, color='y', linestyle='--',linewidth=4)

    font.set_size(10)
    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'raw %s' % var4_str, 'raw %s' % var5_str, 'raw %s' % var6_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str, 'debiased %s' % var4_str, 'debiased %s' % var5_str, 'debiased %s' % var6_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax1.set_xlim(-0.01,0.25)
    ax1.set_ylim(-0.01,1.2)
    ax1.set_xlabel('redshift')
    ax1.set_ylabel('fraction')
    ax1.set_title('All GZ2 galaxies with redshifts, %s' % var_def)
    ax1.text(0.15,1.0,str(len(p_el_raw_values))+' galaxies')

    ax2 = fig.add_subplot(122)

    ax2.plot(zplotbins, p_el_raw_typefrac_maglim, color='r', linestyle='-' ,linewidth=2)
    ax2.plot(zplotbins, p_sp_raw_typefrac_maglim, color='b', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_xx_raw_typefrac_maglim, color='m', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_yy_raw_typefrac_maglim, color='g', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_zz_raw_typefrac_maglim, color='c', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_ww_raw_typefrac_maglim, color='y', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_el_adj_typefrac_maglim, color='r', linestyle='-' ,linewidth=4)
    ax2.plot(zplotbins, p_sp_adj_typefrac_maglim, color='b', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_xx_adj_typefrac_maglim, color='m', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_yy_adj_typefrac_maglim, color='g', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_zz_adj_typefrac_maglim, color='c', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_ww_adj_typefrac_maglim, color='y', linestyle='--',linewidth=4)

    ax2.axvline(zlo, color='k', linestyle='--')
    ax2.axvline(zhi, color='k', linestyle='--')

    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'raw %s' % var4_str, 'raw %s' % var5_str, 'raw %s' % var6_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str, 'debiased %s' % var4_str, 'debiased %s' % var5_str, 'debiased %s' % var6_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.2)
    ax2.set_xlabel('redshift')
    ax2.set_ylabel('fraction')
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)
    ax2.text(0.15,1.0,(str(len(p_el_adj_values))+' galaxies'))

    fig.savefig(plots_path+'%s_%stype_fractions_n5.png' % (task_dict['var_def'],s82_str), dpi=200)

    return None

def plot_type_fractions_n5(task_dict, zlo = 0.01, zhi=0.085, stripe82=False):

    # Load data 

    if stripe82:
        s82_str = 'stripe82_'
        p = pyfits.open(gz2_stripe82_data_file)
    else:
        s82_str = ''
        p = pyfits.open(gz2_full_data_file)

    gzdata = p[1].data
    p.close()
    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']
    p_el_raw_values = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw_values = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw_values = gzdata_withz[task_dict['task_names_wf'][2]]
    p_yy_raw_values = gzdata_withz[task_dict['task_names_wf'][3]]
    p_zz_raw_values = gzdata_withz[task_dict['task_names_wf'][4]]
    task_counts = gzdata_withz[task_dict['task_name_count']]

    zwidth = 0.02
    zplotbins = np.arange(min(edges_redshift[0],zlo),edges_redshift[-1],zwidth)

    n = len(task_dict['var_str'])
    var_def = task_dict['var_def']
    var1_str,var2_str,var3_str,var4_str,var5_str = task_dict['var_str']

    p_el_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var1_str,s82_str),'rb')) 
    p_sp_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var2_str,s82_str),'rb')) 
    p_xx_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var3_str,s82_str),'rb')) 
    p_yy_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var4_str,s82_str),'rb')) 
    p_zz_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var5_str,s82_str),'rb')) 
    corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (var_def,s82_str),'rb')) 

    # Only plot galaxies that had a correction available

    corrgood = corrarr > unknown_ratio
    elgood = np.isfinite(p_el_adj_values)
    spgood = np.isfinite(p_sp_adj_values)
    xxgood = np.isfinite(p_xx_adj_values)
    yygood = np.isfinite(p_yy_adj_values)
    zzgood = np.isfinite(p_zz_adj_values)
    goodinds = corrgood & elgood & spgood & xxgood & yygood & zzgood
    p_el_adj_values = p_el_adj_values[goodinds]
    p_sp_adj_values = p_sp_adj_values[goodinds]
    p_xx_adj_values = p_xx_adj_values[goodinds]
    p_yy_adj_values = p_yy_adj_values[goodinds]
    p_zz_adj_values = p_zz_adj_values[goodinds]
    task_counts_adj = task_counts[goodinds]

    # Empty arrays for the binned type fractions

    p_el_raw_typefrac = np.zeros_like(zplotbins)
    p_sp_raw_typefrac = np.zeros_like(zplotbins)
    p_el_adj_typefrac = np.zeros_like(zplotbins)
    p_sp_adj_typefrac = np.zeros_like(zplotbins)
    p_el_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_el_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_raw_typefrac = np.zeros_like(zplotbins)
    p_xx_adj_typefrac = np.zeros_like(zplotbins)
    p_xx_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_yy_raw_typefrac = np.zeros_like(zplotbins)
    p_yy_adj_typefrac = np.zeros_like(zplotbins)
    p_yy_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_yy_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_zz_raw_typefrac = np.zeros_like(zplotbins)
    p_zz_adj_typefrac = np.zeros_like(zplotbins)
    p_zz_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_zz_adj_typefrac_maglim = np.zeros_like(zplotbins)

    """
    Method 1: take the mean of the raw likelihoods for each galaxy per redshift bin
    Method 2: take the mean of each likelihood weighted by the total number of votes

    Results are within 5% of each other (fraction of < 0.01) at all redshifts.
    """

    # Create a second, magnitude-limited sample
    
    maglimval = appmag_lim - cosmology.dmod_flat(zhi)
    
    mr_raw = mr
    redshift_raw = redshift
    maglim_raw = mr_raw < maglimval
    task_counts_maglim_raw = task_counts[maglim_raw]
    p_el_raw_values_maglim = p_el_raw_values[maglim_raw]
    p_sp_raw_values_maglim = p_sp_raw_values[maglim_raw]
    p_xx_raw_values_maglim = p_xx_raw_values[maglim_raw]
    p_yy_raw_values_maglim = p_yy_raw_values[maglim_raw]
    p_zz_raw_values_maglim = p_zz_raw_values[maglim_raw]

    mr_adj = mr[goodinds]
    redshift_adj = redshift[goodinds]
    maglim_adj = mr_adj < maglimval
    task_counts_maglim_adj = task_counts_adj[maglim_adj]
    p_el_adj_values_maglim = p_el_adj_values[maglim_adj]
    p_sp_adj_values_maglim = p_sp_adj_values[maglim_adj]
    p_xx_adj_values_maglim = p_xx_adj_values[maglim_adj]
    p_yy_adj_values_maglim = p_yy_adj_values[maglim_adj]
    p_zz_adj_values_maglim = p_zz_adj_values[maglim_adj]

    # Loop over redshift bin and find the type fractions at each slice

    for idx,zbin in enumerate(zplotbins):
        gals_in_bin_raw = (redshift_raw >= zbin) & (redshift_raw < (zbin+zwidth))
        weights_raw = task_counts[gals_in_bin_raw]
        p_el_raw_typefrac[idx] = np.mean(p_el_raw_values[gals_in_bin_raw])
        p_sp_raw_typefrac[idx] = np.mean(p_sp_raw_values[gals_in_bin_raw])
        p_xx_raw_typefrac[idx] = np.mean(p_xx_raw_values[gals_in_bin_raw])
        p_yy_raw_typefrac[idx] = np.mean(p_yy_raw_values[gals_in_bin_raw])
        p_zz_raw_typefrac[idx] = np.mean(p_zz_raw_values[gals_in_bin_raw])

        gals_in_bin_adj = (redshift_adj >= zbin) & (redshift_adj < (zbin+zwidth))
        weights_adj = task_counts_adj[gals_in_bin_adj]
        p_el_adj_typefrac[idx] = np.mean(p_el_adj_values[gals_in_bin_adj])
        p_sp_adj_typefrac[idx] = np.mean(p_sp_adj_values[gals_in_bin_adj])
        p_xx_adj_typefrac[idx] = np.mean(p_xx_adj_values[gals_in_bin_adj])
        p_yy_adj_typefrac[idx] = np.mean(p_yy_adj_values[gals_in_bin_adj])
        p_zz_adj_typefrac[idx] = np.mean(p_zz_adj_values[gals_in_bin_adj])

        gals_in_bin_maglim_raw = (redshift_raw[maglim_raw] >= zbin) & (redshift_raw[maglim_raw] < (zbin+zwidth))
        weights_maglim_raw = task_counts_maglim_raw[gals_in_bin_maglim_raw]
        p_el_raw_typefrac_maglim[idx] = np.sum(p_el_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_sp_raw_typefrac_maglim[idx] = np.sum(p_sp_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_xx_raw_typefrac_maglim[idx] = np.sum(p_xx_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_yy_raw_typefrac_maglim[idx] = np.sum(p_yy_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_zz_raw_typefrac_maglim[idx] = np.sum(p_zz_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)

        gals_in_bin_maglim_adj = (redshift_adj[maglim_adj] >= zbin) & (redshift_adj[maglim_adj] < (zbin+zwidth))
        weights_maglim_adj = task_counts_maglim_adj[gals_in_bin_maglim_adj]
        p_el_adj_typefrac_maglim[idx] = np.sum(p_el_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_sp_adj_typefrac_maglim[idx] = np.sum(p_sp_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_xx_adj_typefrac_maglim[idx] = np.sum(p_xx_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_yy_adj_typefrac_maglim[idx] = np.sum(p_yy_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_zz_adj_typefrac_maglim[idx] = np.sum(p_zz_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)

    pickle.dump(p_el_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 

    pickle.dump(p_xx_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_yy_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_zz_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var5_str,s82_str),'wb')) 
    pickle.dump(p_zz_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var5_str,s82_str),'wb')) 
    pickle.dump(p_zz_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var5_str,s82_str),'wb')) 
    pickle.dump(p_zz_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var5_str,s82_str),'wb')) 

    # Plot results

    fig = plt.figure(12, (12,7))
    fig.clf()
    ax1 = fig.add_subplot(121)
    font = FontProperties()

    ax1.plot(zplotbins, p_el_raw_typefrac, color='r', linestyle='-' ,linewidth=2)
    ax1.plot(zplotbins, p_sp_raw_typefrac, color='b', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_xx_raw_typefrac, color='m', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_yy_raw_typefrac, color='g', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_zz_raw_typefrac, color='c', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_el_adj_typefrac, color='r', linestyle='-' ,linewidth=4)
    ax1.plot(zplotbins, p_sp_adj_typefrac, color='b', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_xx_adj_typefrac, color='m', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_yy_adj_typefrac, color='g', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_zz_adj_typefrac, color='c', linestyle='--',linewidth=4)

    font.set_size(10)
    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'raw %s' % var4_str, 'raw %s' % var5_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str, 'debiased %s' % var4_str, 'debiased %s' % var5_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax1.set_xlim(-0.01,0.25)
    ax1.set_ylim(-0.01,1.2)
    ax1.set_xlabel('redshift')
    ax1.set_ylabel('fraction')
    ax1.set_title('All GZ2 galaxies with redshifts, %s' % var_def)
    ax1.text(0.15,1.0,str(len(p_el_raw_values))+' galaxies')

    ax2 = fig.add_subplot(122)

    ax2.plot(zplotbins, p_el_raw_typefrac_maglim, color='r', linestyle='-' ,linewidth=2)
    ax2.plot(zplotbins, p_sp_raw_typefrac_maglim, color='b', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_xx_raw_typefrac_maglim, color='m', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_yy_raw_typefrac_maglim, color='g', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_zz_raw_typefrac_maglim, color='c', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_el_adj_typefrac_maglim, color='r', linestyle='-' ,linewidth=4)
    ax2.plot(zplotbins, p_sp_adj_typefrac_maglim, color='b', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_xx_adj_typefrac_maglim, color='m', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_yy_adj_typefrac_maglim, color='g', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_zz_adj_typefrac_maglim, color='c', linestyle='--',linewidth=4)

    ax2.axvline(zlo, color='k', linestyle='--')
    ax2.axvline(zhi, color='k', linestyle='--')

    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'raw %s' % var4_str, 'raw %s' % var5_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str, 'debiased %s' % var4_str, 'debiased %s' % var5_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.2)
    ax2.set_xlabel('redshift')
    ax2.set_ylabel('fraction')
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)
    ax2.text(0.15,1.0,(str(len(p_el_adj_values))+' galaxies'))

    fig.savefig(plots_path+'%s_%stype_fractions_n5.png' % (task_dict['var_def'],s82_str), dpi=200)

    return None

def plot_type_fractions_n4(task_dict, zlo = 0.01, zhi=0.085, stripe82=False):

    # Load data 

    if stripe82:
        s82_str = 'stripe82_'
        p = pyfits.open(gz2_stripe82_data_file)
    else:
        s82_str = ''
        p = pyfits.open(gz2_full_data_file)

    gzdata = p[1].data
    p.close()
    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']
    p_el_raw_values = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw_values = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw_values = gzdata_withz[task_dict['task_names_wf'][2]]
    p_yy_raw_values = gzdata_withz[task_dict['task_names_wf'][3]]
    task_counts = gzdata_withz[task_dict['task_name_count']]

    zwidth = 0.02
    zplotbins = np.arange(min(edges_redshift[0],zlo),edges_redshift[-1],zwidth)

    n = len(task_dict['var_str'])
    var_def = task_dict['var_def']
    var1_str,var2_str,var3_str,var4_str = task_dict['var_str']

    p_el_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var1_str,s82_str),'rb')) 
    p_sp_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var2_str,s82_str),'rb')) 
    p_xx_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var3_str,s82_str),'rb')) 
    p_yy_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var4_str,s82_str),'rb')) 
    corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (var_def,s82_str),'rb')) 

    # Only plot galaxies that had a correction available

    corrgood = corrarr > unknown_ratio
    elgood = np.isfinite(p_el_adj_values)
    spgood = np.isfinite(p_sp_adj_values)
    xxgood = np.isfinite(p_xx_adj_values)
    yygood = np.isfinite(p_yy_adj_values)
    goodinds = corrgood & elgood & spgood & xxgood & yygood
    p_el_adj_values = p_el_adj_values[goodinds]
    p_sp_adj_values = p_sp_adj_values[goodinds]
    p_xx_adj_values = p_xx_adj_values[goodinds]
    p_yy_adj_values = p_yy_adj_values[goodinds]
    task_counts_adj = task_counts[goodinds]

    # Empty arrays for the binned type fractions

    p_el_raw_typefrac = np.zeros_like(zplotbins)
    p_sp_raw_typefrac = np.zeros_like(zplotbins)
    p_el_adj_typefrac = np.zeros_like(zplotbins)
    p_sp_adj_typefrac = np.zeros_like(zplotbins)
    p_el_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_el_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_raw_typefrac = np.zeros_like(zplotbins)
    p_xx_adj_typefrac = np.zeros_like(zplotbins)
    p_xx_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_yy_raw_typefrac = np.zeros_like(zplotbins)
    p_yy_adj_typefrac = np.zeros_like(zplotbins)
    p_yy_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_yy_adj_typefrac_maglim = np.zeros_like(zplotbins)

    """
    Method 1: take the mean of the raw likelihoods for each galaxy per redshift bin
    Method 2: take the mean of each likelihood weighted by the total number of votes

    Results are within 5% of each other (fraction of < 0.01) at all redshifts.
    """

    # Create a second, magnitude-limited sample
    
    maglimval = appmag_lim - cosmology.dmod_flat(zhi)
    
    mr_raw = mr
    redshift_raw = redshift
    maglim_raw = mr_raw < maglimval
    task_counts_maglim_raw = task_counts[maglim_raw]
    p_el_raw_values_maglim = p_el_raw_values[maglim_raw]
    p_sp_raw_values_maglim = p_sp_raw_values[maglim_raw]
    p_xx_raw_values_maglim = p_xx_raw_values[maglim_raw]
    p_yy_raw_values_maglim = p_yy_raw_values[maglim_raw]

    mr_adj = mr[goodinds]
    redshift_adj = redshift[goodinds]
    maglim_adj = mr_adj < maglimval
    task_counts_maglim_adj = task_counts_adj[maglim_adj]
    p_el_adj_values_maglim = p_el_adj_values[maglim_adj]
    p_sp_adj_values_maglim = p_sp_adj_values[maglim_adj]
    p_xx_adj_values_maglim = p_xx_adj_values[maglim_adj]
    p_yy_adj_values_maglim = p_yy_adj_values[maglim_adj]

    # Loop over redshift bin and find the type fractions at each slice

    for idx,zbin in enumerate(zplotbins):
        gals_in_bin_raw = (redshift_raw >= zbin) & (redshift_raw < (zbin+zwidth))
        weights_raw = task_counts[gals_in_bin_raw]
        p_el_raw_typefrac[idx] = np.mean(p_el_raw_values[gals_in_bin_raw])
        p_sp_raw_typefrac[idx] = np.mean(p_sp_raw_values[gals_in_bin_raw])
        p_xx_raw_typefrac[idx] = np.mean(p_xx_raw_values[gals_in_bin_raw])
        p_yy_raw_typefrac[idx] = np.mean(p_yy_raw_values[gals_in_bin_raw])

        gals_in_bin_adj = (redshift_adj >= zbin) & (redshift_adj < (zbin+zwidth))
        weights_adj = task_counts_adj[gals_in_bin_adj]
        p_el_adj_typefrac[idx] = np.mean(p_el_adj_values[gals_in_bin_adj])
        p_sp_adj_typefrac[idx] = np.mean(p_sp_adj_values[gals_in_bin_adj])
        p_xx_adj_typefrac[idx] = np.mean(p_xx_adj_values[gals_in_bin_adj])
        p_yy_adj_typefrac[idx] = np.mean(p_yy_adj_values[gals_in_bin_adj])

        gals_in_bin_maglim_raw = (redshift_raw[maglim_raw] >= zbin) & (redshift_raw[maglim_raw] < (zbin+zwidth))
        weights_maglim_raw = task_counts_maglim_raw[gals_in_bin_maglim_raw]
        p_el_raw_typefrac_maglim[idx] = np.sum(p_el_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_sp_raw_typefrac_maglim[idx] = np.sum(p_sp_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_xx_raw_typefrac_maglim[idx] = np.sum(p_xx_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_yy_raw_typefrac_maglim[idx] = np.sum(p_yy_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)

        gals_in_bin_maglim_adj = (redshift_adj[maglim_adj] >= zbin) & (redshift_adj[maglim_adj] < (zbin+zwidth))
        weights_maglim_adj = task_counts_maglim_adj[gals_in_bin_maglim_adj]
        p_el_adj_typefrac_maglim[idx] = np.sum(p_el_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_sp_adj_typefrac_maglim[idx] = np.sum(p_sp_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_xx_adj_typefrac_maglim[idx] = np.sum(p_xx_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_yy_adj_typefrac_maglim[idx] = np.sum(p_yy_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)

    pickle.dump(p_el_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 

    pickle.dump(p_xx_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_yy_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var4_str,s82_str),'wb')) 
    pickle.dump(p_yy_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var4_str,s82_str),'wb')) 

    # Plot results

    fig = plt.figure(12, (12,7))
    fig.clf()
    ax1 = fig.add_subplot(121)
    font = FontProperties()

    ax1.plot(zplotbins, p_el_raw_typefrac, color='r', linestyle='-' ,linewidth=2)
    ax1.plot(zplotbins, p_sp_raw_typefrac, color='b', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_xx_raw_typefrac, color='m', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_yy_raw_typefrac, color='g', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_el_adj_typefrac, color='r', linestyle='-' ,linewidth=4)
    ax1.plot(zplotbins, p_sp_adj_typefrac, color='b', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_xx_adj_typefrac, color='m', linestyle='--',linewidth=4)
    ax1.plot(zplotbins, p_yy_adj_typefrac, color='g', linestyle='--',linewidth=4)

    font.set_size(10)
    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'raw %s' % var4_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str, 'debiased %s' % var4_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax1.set_xlim(-0.01,0.25)
    ax1.set_ylim(-0.01,1.2)
    ax1.set_xlabel('redshift')
    ax1.set_ylabel('fraction')
    ax1.set_title('All GZ2 galaxies with redshifts, %s' % var_def)
    ax1.text(0.15,1.0,str(len(p_el_raw_values))+' galaxies')

    ax2 = fig.add_subplot(122)

    ax2.plot(zplotbins, p_el_raw_typefrac_maglim, color='r', linestyle='-' ,linewidth=2)
    ax2.plot(zplotbins, p_sp_raw_typefrac_maglim, color='b', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_xx_raw_typefrac_maglim, color='m', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_yy_raw_typefrac_maglim, color='g', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_el_adj_typefrac_maglim, color='r', linestyle='-' ,linewidth=4)
    ax2.plot(zplotbins, p_sp_adj_typefrac_maglim, color='b', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_xx_adj_typefrac_maglim, color='m', linestyle='--',linewidth=4)
    ax2.plot(zplotbins, p_yy_adj_typefrac_maglim, color='g', linestyle='--',linewidth=4)

    ax2.axvline(zlo, color='k', linestyle='--')
    ax2.axvline(zhi, color='k', linestyle='--')

    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'raw %s' % var4_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str, 'debiased %s' % var4_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.2)
    ax2.set_xlabel('redshift')
    ax2.set_ylabel('fraction')
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)
    ax2.text(0.15,1.0,(str(len(p_el_adj_values))+' galaxies'))

    fig.savefig(plots_path+'%s_%stype_fractions_n4.png' % (task_dict['var_def'],s82_str), dpi=200)

    return None

def plot_type_fractions_n3(task_dict, zlo = 0.01, zhi=0.085, stripe82=False):

    # Load data 

    if stripe82:
        s82_str = 'stripe82_'
        p = pyfits.open(gz2_stripe82_data_file)
    else:
        s82_str = ''
        p = pyfits.open(gz2_full_data_file)

    gzdata = p[1].data
    p.close()
    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']
    p_el_raw_values = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw_values = gzdata_withz[task_dict['task_names_wf'][1]]
    p_xx_raw_values = gzdata_withz[task_dict['task_names_wf'][2]]
    task_counts = gzdata_withz[task_dict['task_name_count']]

    zwidth = 0.02
    zplotbins = np.arange(min(edges_redshift[0],zlo),edges_redshift[-1],zwidth)

    n = len(task_dict['var_str'])
    var_def = task_dict['var_def']
    var1_str,var2_str,var3_str = task_dict['var_str']

    p_el_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var1_str,s82_str),'rb')) 
    p_sp_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var2_str,s82_str),'rb')) 
    p_xx_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var3_str,s82_str),'rb')) 
    corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (var_def,s82_str),'rb')) 

    # Only plot galaxies that had a correction available

    corrgood = corrarr > unknown_ratio
    elgood = np.isfinite(p_el_adj_values)
    spgood = np.isfinite(p_sp_adj_values)
    xxgood = np.isfinite(p_xx_adj_values)
    task_lim = task_counts > task_dict['min_classifications']

    goodinds = corrgood & elgood & spgood & xxgood & task_lim

    redshift = redshift[goodinds]
    r50_kpc = r50_kpc[goodinds]
    mr = mr[goodinds]

    p_el_raw_values = p_el_raw_values[goodinds]
    p_sp_raw_values = p_sp_raw_values[goodinds]
    p_xx_raw_values = p_xx_raw_values[goodinds]

    p_el_adj_values = p_el_adj_values[goodinds]
    p_sp_adj_values = p_sp_adj_values[goodinds]
    p_xx_adj_values = p_xx_adj_values[goodinds]

    corrarr = corrarr[goodinds]

    task_counts = task_counts[goodinds]

    # Empty arrays for the binned type fractions

    p_el_raw_typefrac = np.zeros_like(zplotbins)
    p_sp_raw_typefrac = np.zeros_like(zplotbins)
    p_el_adj_typefrac = np.zeros_like(zplotbins)
    p_sp_adj_typefrac = np.zeros_like(zplotbins)
    p_xx_raw_typefrac = np.zeros_like(zplotbins)
    p_xx_adj_typefrac = np.zeros_like(zplotbins)
    p_el_raw_typefrac_err = np.zeros_like(zplotbins)
    p_sp_raw_typefrac_err = np.zeros_like(zplotbins)
    p_el_adj_typefrac_err = np.zeros_like(zplotbins)
    p_sp_adj_typefrac_err = np.zeros_like(zplotbins)
    p_xx_raw_typefrac_err = np.zeros_like(zplotbins)
    p_xx_adj_typefrac_err = np.zeros_like(zplotbins)
    p_el_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_el_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_xx_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_el_raw_typefrac_maglim_err = np.zeros_like(zplotbins)
    p_sp_raw_typefrac_maglim_err = np.zeros_like(zplotbins)
    p_el_adj_typefrac_maglim_err = np.zeros_like(zplotbins)
    p_sp_adj_typefrac_maglim_err = np.zeros_like(zplotbins)
    p_xx_raw_typefrac_maglim_err = np.zeros_like(zplotbins)
    p_xx_adj_typefrac_maglim_err = np.zeros_like(zplotbins)

    """
    Method 1: take the mean of the raw likelihoods for each galaxy per redshift bin
    Method 2: take the mean of each likelihood weighted by the total number of votes

    Results are within 5% of each other (fraction of < 0.01) at all redshifts.
    """

    # Create a second, magnitude-limited sample
    
    maglimval = appmag_lim - cosmology.dmod_flat(zhi)
    
    maglim = mr < maglimval
    task_counts_maglim = task_counts[maglim]
    p_el_raw_values_maglim = p_el_raw_values[maglim]
    p_sp_raw_values_maglim = p_sp_raw_values[maglim]
    p_xx_raw_values_maglim = p_xx_raw_values[maglim]

    p_el_adj_values_maglim = p_el_adj_values[maglim]
    p_sp_adj_values_maglim = p_sp_adj_values[maglim]
    p_xx_adj_values_maglim = p_xx_adj_values[maglim]

    # Loop over redshift bin and find the type fractions at each slice

    for idx,zbin in enumerate(zplotbins):
        gals_in_bin_raw = (redshift_raw >= zbin) & (redshift_raw < (zbin+zwidth))
        weights_raw = task_counts[gals_in_bin_raw]
        p_el_raw_typefrac[idx] = np.mean(p_el_raw_values[gals_in_bin_raw])
        p_sp_raw_typefrac[idx] = np.mean(p_sp_raw_values[gals_in_bin_raw])
        p_xx_raw_typefrac[idx] = np.mean(p_xx_raw_values[gals_in_bin_raw])
        p_el_raw_typefrac_err[idx] = tsem(p_el_raw_values[gals_in_bin_raw])
        p_sp_raw_typefrac_err[idx] = tsem(p_sp_raw_values[gals_in_bin_raw])
        p_xx_raw_typefrac_err[idx] = tsem(p_xx_raw_values[gals_in_bin_raw])

        gals_in_bin_adj = (redshift_adj >= zbin) & (redshift_adj < (zbin+zwidth))
        weights_adj = task_counts[gals_in_bin_adj]
        p_el_adj_typefrac[idx] = np.mean(p_el_adj_values[gals_in_bin_adj])
        p_sp_adj_typefrac[idx] = np.mean(p_sp_adj_values[gals_in_bin_adj])
        p_xx_adj_typefrac[idx] = np.mean(p_xx_adj_values[gals_in_bin_adj])
        p_el_adj_typefrac_err[idx] = tsem(p_el_adj_values[gals_in_bin_adj])
        p_sp_adj_typefrac_err[idx] = tsem(p_sp_adj_values[gals_in_bin_adj])
        p_xx_adj_typefrac_err[idx] = tsem(p_xx_adj_values[gals_in_bin_adj])

        gals_in_bin_maglim_raw = (redshift_raw[maglim_raw] >= zbin) & (redshift_raw[maglim_raw] < (zbin+zwidth))
        weights_maglim_raw = task_counts_maglim[gals_in_bin_maglim_raw]
        """
        p_el_raw_typefrac_maglim[idx] = np.sum(p_el_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_sp_raw_typefrac_maglim[idx] = np.sum(p_sp_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_xx_raw_typefrac_maglim[idx] = np.sum(p_xx_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        """
        p_el_raw_typefrac_maglim[idx] = np.mean(p_el_raw_values_maglim[gals_in_bin_maglim_raw])
        p_sp_raw_typefrac_maglim[idx] = np.mean(p_sp_raw_values_maglim[gals_in_bin_maglim_raw])
        p_xx_raw_typefrac_maglim[idx] = np.mean(p_xx_raw_values_maglim[gals_in_bin_maglim_raw])
        p_el_raw_typefrac_maglim_err[idx] = tsem(p_el_raw_values_maglim[gals_in_bin_maglim_raw])
        p_sp_raw_typefrac_maglim_err[idx] = tsem(p_sp_raw_values_maglim[gals_in_bin_maglim_raw])
        p_xx_raw_typefrac_maglim_err[idx] = tsem(p_xx_raw_values_maglim[gals_in_bin_maglim_raw])

        gals_in_bin_maglim_adj = (redshift_adj[maglim_adj] >= zbin) & (redshift_adj[maglim_adj] < (zbin+zwidth))
        weights_maglim_adj = task_counts_maglim[gals_in_bin_maglim_adj]

        p_el_adj_typefrac_maglim[idx] =  np.mean(p_el_adj_values_maglim[gals_in_bin_maglim_adj])
        p_sp_adj_typefrac_maglim[idx] =  np.mean(p_sp_adj_values_maglim[gals_in_bin_maglim_adj])
        p_xx_adj_typefrac_maglim[idx] =  np.mean(p_xx_adj_values_maglim[gals_in_bin_maglim_adj])
        p_el_adj_typefrac_maglim_err[idx] = tsem(p_el_adj_values_maglim[gals_in_bin_maglim_adj])
        p_sp_adj_typefrac_maglim_err[idx] = tsem(p_sp_adj_values_maglim[gals_in_bin_maglim_adj])
        p_xx_adj_typefrac_maglim_err[idx] = tsem(p_xx_adj_values_maglim[gals_in_bin_maglim_adj])

    pickle.dump(p_el_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 

    pickle.dump(p_xx_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 
    pickle.dump(p_xx_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var3_str,s82_str),'wb')) 

    # Plot results

    fig = plt.figure(12, (12,7))
    fig.clf()
    ax1 = fig.add_subplot(121)
    font = FontProperties()

    ax1.errorbar(zplotbins, p_el_raw_typefrac, yerr=p_el_raw_typefrac_err, color='r', linestyle='-' ,linewidth=2, capsize=0, elinewidth=1)
    ax1.errorbar(zplotbins, p_sp_raw_typefrac, yerr=p_sp_raw_typefrac_err, color='b', linestyle='--',linewidth=2, capsize=0, elinewidth=1)
    ax1.errorbar(zplotbins, p_xx_raw_typefrac, yerr=p_xx_raw_typefrac_err, color='m', linestyle='--',linewidth=2, capsize=0, elinewidth=1)
    ax1.errorbar(zplotbins, p_el_adj_typefrac, yerr=p_el_adj_typefrac_err, color='r', linestyle='-' ,linewidth=4, capsize=0, elinewidth=2)
    ax1.errorbar(zplotbins, p_sp_adj_typefrac, yerr=p_sp_adj_typefrac_err, color='b', linestyle='--',linewidth=4, capsize=0, elinewidth=2)
    ax1.errorbar(zplotbins, p_xx_adj_typefrac, yerr=p_xx_adj_typefrac_err, color='m', linestyle='--',linewidth=4, capsize=0, elinewidth=2)

    font.set_size(10)
    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax1.set_xlim(-0.01,0.25)
    ax1.set_ylim(-0.01,1.2)
    ax1.set_xlabel('redshift')
    ax1.set_ylabel('mean vote fraction')
    ax1.set_title('All GZ2 galaxies with redshifts, %s' % var_def)
    ax1.text(0.15,1.0,str(len(p_el_raw_values))+' galaxies')

    ax2 = fig.add_subplot(122)

    ax2.errorbar(zplotbins, p_el_raw_typefrac_maglim, yerr=p_el_raw_typefrac_maglim_err, color='r', linestyle='-' ,linewidth=2, capsize=0, elinewidth=1)
    ax2.errorbar(zplotbins, p_sp_raw_typefrac_maglim, yerr=p_sp_raw_typefrac_maglim_err, color='b', linestyle='--',linewidth=2, capsize=0, elinewidth=1)
    ax2.errorbar(zplotbins, p_xx_raw_typefrac_maglim, yerr=p_xx_raw_typefrac_maglim_err, color='m', linestyle='--',linewidth=2, capsize=0, elinewidth=1)
    ax2.errorbar(zplotbins, p_el_adj_typefrac_maglim, yerr=p_el_adj_typefrac_maglim_err, color='r', linestyle='-' ,linewidth=4, capsize=0, elinewidth=2)
    ax2.errorbar(zplotbins, p_sp_adj_typefrac_maglim, yerr=p_sp_adj_typefrac_maglim_err, color='b', linestyle='--',linewidth=4, capsize=0, elinewidth=2)
    ax2.errorbar(zplotbins, p_xx_adj_typefrac_maglim, yerr=p_xx_adj_typefrac_maglim_err, color='m', linestyle='--',linewidth=4, capsize=0, elinewidth=2)

    ax2.axvline(zlo, color='k', linestyle='--')
    ax2.axvline(zhi, color='k', linestyle='--')

    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'raw %s' % var3_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str, 'debiased %s' % var3_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.2)
    ax2.set_xlabel('redshift')
    ax2.set_ylabel('mean vote fraction')
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)
    ax2.text(0.15,1.0,(str(len(p_el_adj_values))+' galaxies'))

    fig.savefig(plots_path+'%s_%stype_fractions_n3.png' % (task_dict['var_def'],s82_str), dpi=200)

    return None

def plot_type_fractions(task_dict, zlo = 0.01, zhi=0.085, stripe82=False):

    # Load data 

    if stripe82:
        s82_str = 'stripe82_'
        p = pyfits.open(gz2_stripe82_data_file)
    else:
        s82_str = ''
        p = pyfits.open(gz2_full_data_file)

    gzdata = p[1].data
    p.close()
    gzdata_withz = gzdata[np.isfinite(gzdata['redshift'])]

    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

    redshift = gzdata_withz['REDSHIFT']
    mr = gzdata_withz['PETROMAG_MR']
    r50_kpc = gzdata_withz['PETROR50_R_KPC']
    p_el_raw_values = gzdata_withz[task_dict['task_names_wf'][0]]
    p_sp_raw_values = gzdata_withz[task_dict['task_names_wf'][1]]
    task_counts = gzdata_withz[task_dict['task_name_count']]

    zwidth = 0.02
    zplotbins = np.arange(max(edges_redshift[0],zlo),edges_redshift[-1],zwidth)

    n = len(task_dict['var_str'])
    var_def = task_dict['var_def']
    var1_str,var2_str = task_dict['var_str']

    p_el_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var1_str,s82_str),'rb')) 
    p_sp_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var2_str,s82_str),'rb')) 
    if n == 3:
        var3_str = task_dict['var_str'][2]
        p_xx_adj_values = pickle.load(open(pkl_path+'%s_%s_%sadj.pkl' % (var_def,var3_str,s82_str),'rb')) 
    corrarr         = pickle.load(open(pkl_path+'%s_%scorrarr.pkl' % (var_def,s82_str),'rb')) 

    # Only plot galaxies that had a correction available

    corrgood = corrarr > unknown_ratio
    p_el_adj_values = p_el_adj_values[corrgood]
    p_sp_adj_values = p_sp_adj_values[corrgood]
    task_counts_adj = task_counts[corrgood]

    # Empty arrays for the binned type fractions

    p_el_raw_typefrac = np.zeros_like(zplotbins)
    p_sp_raw_typefrac = np.zeros_like(zplotbins)
    p_el_adj_typefrac = np.zeros_like(zplotbins)
    p_sp_adj_typefrac = np.zeros_like(zplotbins)
    p_el_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_raw_typefrac_maglim = np.zeros_like(zplotbins)
    p_el_adj_typefrac_maglim = np.zeros_like(zplotbins)
    p_sp_adj_typefrac_maglim = np.zeros_like(zplotbins)

    """
    Method 1: take the mean of the raw likelihoods for each galaxy per redshift bin
    Method 2: take the mean of each likelihood weighted by the total number of votes

    Results are within 5% of each other (fraction of < 0.01) at all redshifts.
    """

    # Create a second, magnitude-limited sample
    
    maglimval = appmag_lim - cosmology.dmod_flat(zhi)
    
    mr_raw = mr
    redshift_raw = redshift
    maglim_raw = mr_raw < maglimval
    task_counts_maglim_raw = task_counts[maglim_raw]
    p_el_raw_values_maglim = p_el_raw_values[maglim_raw]
    p_sp_raw_values_maglim = p_sp_raw_values[maglim_raw]

    mr_adj = mr[corrgood]
    redshift_adj = redshift[corrgood]
    maglim_adj = mr_adj < maglimval
    task_counts_maglim_adj = task_counts_adj[maglim_adj]
    p_el_adj_values_maglim = p_el_adj_values[maglim_adj]
    p_sp_adj_values_maglim = p_sp_adj_values[maglim_adj]

    # Loop over redshift bin and find the type fractions at each slice

    for idx,zbin in enumerate(zplotbins):
        gals_in_bin_raw = (redshift_raw >= zbin) & (redshift_raw < (zbin+zwidth))
        weights_raw = task_counts[gals_in_bin_raw]
        p_el_raw_typefrac[idx] = np.mean(p_el_raw_values[gals_in_bin_raw])
        p_sp_raw_typefrac[idx] = np.mean(p_sp_raw_values[gals_in_bin_raw])

        gals_in_bin_adj = (redshift_adj >= zbin) & (redshift_adj < (zbin+zwidth))
        weights_adj = task_counts_adj[gals_in_bin_adj]
        p_el_adj_typefrac[idx] = np.mean(p_el_adj_values[gals_in_bin_adj])
        p_sp_adj_typefrac[idx] = np.mean(p_sp_adj_values[gals_in_bin_adj])

        gals_in_bin_maglim_raw = (redshift_raw[maglim_raw] >= zbin) & (redshift_raw[maglim_raw] < (zbin+zwidth))
        weights_maglim_raw = task_counts_maglim_raw[gals_in_bin_maglim_raw]
        p_el_raw_typefrac_maglim[idx] = np.sum(p_el_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)
        p_sp_raw_typefrac_maglim[idx] = np.sum(p_sp_raw_values_maglim[gals_in_bin_maglim_raw] * weights_maglim_raw) / np.sum(weights_maglim_raw)

        gals_in_bin_maglim_adj = (redshift_adj[maglim_adj] >= zbin) & (redshift_adj[maglim_adj] < (zbin+zwidth))
        weights_maglim_adj = task_counts_maglim_adj[gals_in_bin_maglim_adj]
        p_el_adj_typefrac_maglim[idx] = np.sum(p_el_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)
        p_sp_adj_typefrac_maglim[idx] = np.sum(p_sp_adj_values_maglim[gals_in_bin_maglim_adj] * weights_maglim_adj) / np.sum(weights_maglim_adj)

    pickle.dump(p_el_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac,        open(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac,        open(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_raw_typefrac_maglim, open(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 
    pickle.dump(p_el_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'wb')) 
    pickle.dump(p_sp_adj_typefrac_maglim, open(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'wb')) 

    # Plot results

    fig = plt.figure(12, (12,7))
    fig.clf()
    ax1 = fig.add_subplot(121)
    font = FontProperties()

    ax1.plot(zplotbins, p_el_raw_typefrac, color='r', linestyle='-' ,linewidth=2)
    ax1.plot(zplotbins, p_sp_raw_typefrac, color='b', linestyle='--',linewidth=2)
    ax1.plot(zplotbins, p_el_adj_typefrac, color='r', linestyle='-' ,linewidth=4)
    ax1.plot(zplotbins, p_sp_adj_typefrac, color='b', linestyle='--',linewidth=4)

    font.set_size(10)
    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax1.set_xlim(-0.01,0.25)
    ax1.set_ylim(-0.01,1.2)
    ax1.set_xlabel('redshift')
    ax1.set_ylabel('fraction')
    ax1.set_title('All GZ2 galaxies with redshifts, %s' % var_def)
    ax1.text(0.15,1.0,str(len(p_el_raw_values))+' galaxies')

    ax2 = fig.add_subplot(122)

    ax2.plot(zplotbins, p_el_raw_typefrac_maglim, color='r', linestyle='-' ,linewidth=2)
    ax2.plot(zplotbins, p_sp_raw_typefrac_maglim, color='b', linestyle='--',linewidth=2)
    ax2.plot(zplotbins, p_el_adj_typefrac_maglim, color='r', linestyle='-' ,linewidth=4)
    ax2.plot(zplotbins, p_sp_adj_typefrac_maglim, color='b', linestyle='--',linewidth=4)

    ax2.axvline(zlo, color='k', linestyle='--')
    ax2.axvline(zhi, color='k', linestyle='--')

    plt.legend(('raw %s' % var1_str, 'raw %s' % var2_str, 'debiased %s' % var1_str, 'debiased %s' % var2_str), 'upper left', shadow=True, fancybox=True, prop=font)

    ax2.set_xlim(-0.01,0.25)
    ax2.set_ylim(-0.01,1.2)
    ax2.set_xlabel('redshift')
    ax2.set_ylabel('fraction')
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)
    ax2.text(0.15,1.0,(str(len(p_el_adj_values))+' galaxies'))

    fig.savefig(plots_path+'%s_%stype_fractions.png' % (task_dict['var_def'],s82_str), dpi=200)

    return None

def plot_all_baselines():

    fig = plt.figure(13, (12,7))
    fig.clf()

    taskstrings=('smooth','edgeon','bar','spiral','odd')
    titlenames = ('Smooth or features','Edge-on','Bar', 'Spiral structure','Anything odd?')
    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')
    odd = get_task_dict('odd')

    tasklist = (smooth,edgeon,bar,spiral,odd)

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

        ax = fig.add_subplot(2,3,idx+1,aspect=1)

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
        if idx in (3,4):
            ax.set_xlabel(r'$M_R [mag]$',fontsize=16)
        if idx in (0,3):
            ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
        ax.set_aspect('auto')
        rc(('xtick','ytick'), labelsize=12)
        cb = plt.colorbar(im,orientation='vertical')
        cb.set_label(r'$%s(N_{%s}/N_{%s})$' % (label_prefix,var1_str,var2_str),fontsize=16)

        SBlim_mag = (SBlim_app - cosmology.dmod_flat(np.mean(centers_redshift))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(centers_redshift)))**2))
        absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(centers_redshift))
        absmag_lim_loz = appmag_lim - cosmology.dmod_flat(0.0005)
        absmag_lim_hiz = appmag_lim - cosmology.dmod_flat(0.25)
        size_1arcsec = cosmology.ang_scale_flat(np.mean(centers_redshift))
        ax.autoscale(False)
        ax.plot(SBlim_mag, SBlim_size,'w--')
        ax.axhline(size_1arcsec, color='w', linestyle='dashed')

    fig.savefig(dropbox_figs_path+'gz2_baselines.eps', dpi=200)

    return None

def plot_all_type_fractions(zlo = 0.01, zhi=0.085, stripe82=False):

    fig = plt.figure(14, (12,7))
    fig.clf()

    taskstrings=('smooth','edgeon','bar','spiral','odd')
    titlenames = ('Smooth or features','Edge-on','Bar', 'Spiral structure','Anything odd?')
    smooth = get_task_dict('smooth')
    edgeon = get_task_dict('edgeon')
    bar = get_task_dict('bar')
    spiral = get_task_dict('spiral')
    odd = get_task_dict('odd')

    tasklist = (smooth,edgeon,bar,spiral,odd)

    for idx,task_dict in enumerate(tasklist):

        var_def = task_dict['var_def']
        var1_str,var2_str = task_dict['var_str']

        centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = get_bins(task_dict)

        zwidth = 0.02
        zplotbins = np.arange(max(edges_redshift[0],zlo),edges_redshift[-1],zwidth)
        maglimval = appmag_lim - cosmology.dmod_flat(zhi)

        s82_str = ''
        if stripe82:
            s82_str = 'stripe82_'

        p_el_raw_typefrac        = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var1_str,s82_str),'rb')) 
        p_sp_raw_typefrac        = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac.pkl' % (var_def,var2_str,s82_str),'rb')) 
        p_el_adj_typefrac        = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var1_str,s82_str),'rb')) 
        p_sp_adj_typefrac        = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac.pkl' % (var_def,var2_str,s82_str),'rb')) 
        p_el_raw_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'rb')) 
        p_sp_raw_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sraw_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'rb')) 
        p_el_adj_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var1_str,s82_str),'rb')) 
        p_sp_adj_typefrac_maglim = pickle.load(file(pkl_path+'%s_%s_%sadj_typefrac_maglim.pkl' % (var_def,var2_str,s82_str),'rb')) 

        ax2 = fig.add_subplot(2,3,idx+1)
        font = FontProperties()

        ax2.plot(zplotbins, p_el_raw_typefrac_maglim, color='r', linestyle='-' ,linewidth=2)
        ax2.plot(zplotbins, p_sp_raw_typefrac_maglim, color='b', linestyle='--',linewidth=2)
        ax2.plot(zplotbins, p_el_adj_typefrac_maglim, color='r', linestyle='-' ,linewidth=4)
        ax2.plot(zplotbins, p_sp_adj_typefrac_maglim, color='b', linestyle='--',linewidth=4)

        ax2.axvline(zlo, color='k', linestyle='--')
        ax2.axvline(zhi, color='k', linestyle='--')

        plt.legend((var1_str,var2_str), 'center right', shadow=True, fancybox=True, prop=font)

        ax2.set_xlim(-0.01,0.25)
        ax2.set_ylim(-0.01,1.01)
        if idx in (3,4):
            ax2.set_xlabel('redshift')
        if idx in (0,3):
            ax2.set_ylabel('fraction')
        ax2.set_title(titlenames[idx])

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
    ax2.set_title(r'$M_r < %4.2f$' % maglimval)

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

    p = pyfits.open(gz2_full_data_file)
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
