"""

Modules for attempting to determine the 'overall' ratio of tasks, instead of one vs. the other. 

For example this uses the ratio p_el / (p_el + p_sp), rather than just p_el / p_sp

Moved here temporarily to avoid clutter in galaxyzoo2.py

"""
def fit_overall_ratio_baseline(task_dict, nboot=50, plot=False, unset_vrange=False):

    var_def = task_dict['var_def']
    var1_str = task_dict['var_str'][0]
    var2_str = task_dict['var_str'][1]
    bintype = task_dict['bintype']

    if bintype is 'counts':
        f_bins = fits_path+'%s_%s_binned_%s.fits' % (var_def,var1_str,bintype)
        p_bins = pyfits.open(f_bins)

    if bintype is 'rawlikelihood':
        f_bins = fits_path+'%s_binned_%s.fits' % (var_def,bintype)
        p_bins = pyfits.open(f_bins)

    zbins = p_bins['REDSHIFT_BIN_CENTERS'].data['centers']
    magbins = p_bins['MR_BIN_CENTERS'].data['centers']
    sizebins = p_bins['R50_KPC_BIN_CENTERS'].data['centers']

    zedges = p_bins['REDSHIFT_BIN_EDGES'].data['edges']
    magedges = p_bins['MR_BIN_EDGES'].data['edges']
    sizeedges = p_bins['R50_KPC_BIN_EDGES'].data['edges']

    p_bins.close()

    for r_idx, ratiovar_str in enumerate(task_dict['var_str']):

        sigma = pyfits.getdata(fits_path+'%s_%s_local_overall_ratio_baseline_sigma.fits' % (var_def,ratiovar_str))
        ratio = pyfits.getdata(fits_path+'%s_%s_local_overall_ratio_baseline.fits' % (var_def,ratiovar_str))
        
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

        # Normalize the weights

        weights /= weights.max()
        
        # Fit baseline ratio with function

        if r_idx==0:
            p, perr = bootstrap_fmin(ratio_minfunc, pinit,
                                     magbins, sizebins, ratio, weights, nboot)
            chi2r = ratio_minfunc(p, magbins, sizebins, ratio, weights)

            ratio_baseline_fit = ratio_function(p, magbins, sizebins)
            pickle.dump((p, perr), file(pkl_path+'%s_ratio_baseline_fit.pkl' % var_def, 'w'))
            res = ratio - ratio_baseline_fit
            normres = res/sigma
            np.putmask(normres, weights < 0.000001, unknown_ratio)

            if plot:
                plot_ratio_function(p,magbins,sizebins, unset_vrange,fignum=r_idx+1)

        if r_idx==1:
            p, perr = bootstrap_fmin(ratio_minfunc_flip, pinit,
                                     magbins, sizebins, ratio, weights, nboot)
            chi2r = ratio_minfunc_flip(p, magbins, sizebins, ratio, weights)

            ratio_baseline_fit = ratio_function_flip(p, magbins, sizebins)
            pickle.dump((p, perr), file(pkl_path+'%s_ratio_baseline_fit.pkl' % var_def, 'w'))
            res = ratio - ratio_baseline_fit
            normres = res/sigma
            np.putmask(normres, weights < 0.000001, unknown_ratio)

            if plot:
                plot_ratio_function(p,magbins,sizebins, unset_vrange,fignum=r_idx+1)

        print 'param =', p
        print 'param errors =', perr
        print 'param errors (%) =', perr/p
        print 'chi2r =', chi2r

        pyfits.writeto(fits_path+'%s_%s_local_overall_ratio_baseline_fit.fits' % (var_def,ratiovar_str),
                       ratio_baseline_fit, clobber=True)    
        pyfits.writeto(fits_path+'%s_%s_local_overall_ratio_baseline_fitres.fits' % (var_def,ratiovar_str),
                       res, clobber=True)    
        pyfits.writeto(fits_path+'%s_%s_local_overall_ratio_baseline_fitnormres.fits' % (var_def,ratiovar_str),
                       normres, clobber=True)    

    return None

def determine_overall_ratio_baseline(task_dict):

    var_def = task_dict['var_def']
    var1_str = task_dict['var_str'][0]
    var2_str = task_dict['var_str'][1]
    bintype = task_dict['bintype']

    # Load the GZ2 sample data

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

    if bintype is 'rawlikelihood':
        p_allvar = pyfits.open(fits_path+'%s_binned_%s.fits' % (var_def,bintype))
        d_allvar = p_allvar[0].data.astype(float)                         # Data is pre-binned
        d_var1 = np.squeeze(d_allvar[0,:,:,:])
        d_var2 = np.squeeze(d_allvar[1,:,:,:])

        zbins = p_allvar['REDSHIFT_BIN_CENTERS'].data['centers']
        magbins = p_allvar['MR_BIN_CENTERS'].data['centers']
        sizebins = p_allvar['R50_KPC_BIN_CENTERS'].data['centers']

        p_allvar.close()

    # Bin sizes are set when FITS file is created

    n_magbin = len(magbins)
    n_sizebin = len(sizebins)

    # Trim the data so it only includes counts in the relevant magnitude and size bins

    ratiovars = (d_var1,d_var2)
    ratiovar_strings = (var1_str,var2_str)

    # Loop over all variables

    for ratiovar, ratiovar_str in zip(ratiovars,ratiovar_strings):

        # Empty 2-D arrays to store the ratio and number of galaxies for chosen variables

        ratio_baseline = np.zeros((n_magbin, n_sizebin), np.float) + unknown_ratio
        counts_baseline = np.zeros((n_magbin, n_sizebin, 2), np.int) + unknown_ratio
        redshift_baseline = np.zeros((n_magbin, n_sizebin), np.float) + unknown_ratio

        # Loop over each slice in redshift space. Cells without entries will be filled in the first available loop.

        for z_index, zsel in enumerate(zbins[1:]):
            rv_slice = ratiovar[z_index,:,:]
            var1_slice = d_var1[z_index, :, :]
            var2_slice = d_var2[z_index, :, :]
         
            # Create mask for cells with low counts
            mask = (var1_slice + var2_slice) < task_dict['min_galperbin']
         
            # Compute the task ratio for the entire array

            if task_dict['ratio_type'] is 'linear':
                ratio = rv_slice.astype(float)/(var1_slice + var2_slice)
            else:
                ratio = np.log10(rv_slice.astype(float)/(var1_slice + var2_slice))

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

            # Populate the 2-D empty arrays with the ratio for all non-masked cells

            ratio_baseline[select] = ratio[select]
            counts_baseline[select] = np.transpose([var1_slice[select], var2_slice[select]])
            redshift_baseline[select_redshift] = redshift_this_slice[select_redshift]

        ratio_baseline_masked = ma.masked_less_equal(ratio_baseline, unknown_ratio)
        counts_baseline_masked = ma.masked_less_equal(counts_baseline, unknown_ratio)
        redshift_baseline_masked = ma.masked_less_equal(redshift_baseline, unknown_ratio)

        # Write the results to FITS files

        pyfits.writeto(fits_path+'%s_%s_local_overall_ratio_baseline.fits' % (var_def,ratiovar_str),
                       ratio_baseline, clobber=True)    
        pyfits.writeto(fits_path+'%s_%s_local_overall_counts_baseline.fits' % (var_def,ratiovar_str),
                       counts_baseline, clobber=True)    
        pyfits.writeto(fits_path+'%s_%s_local_overall_redshift_baseline.fits' % (var_def,ratiovar_str),
                       redshift_baseline, clobber=True)    

        pickle.dump(ratio_baseline_masked, open(pkl_path+'%s_%s_local_overall_ratio_baseline_masked.pkl' % (var_def,ratiovar_str),'wb')) 
        pickle.dump(counts_baseline_masked, open(pkl_path+'%s_%s_local_overall_counts_baseline_masked.pkl' % (var_def,ratiovar_str),'wb')) 
        pickle.dump(redshift_baseline_masked, open(pkl_path+'%s_%s_local_overall_redshift_baseline_masked.pkl' % (var_def,ratiovar_str),'wb')) 

    return None

def determine_overall_ratio_baseline_sigma(task_dict, plot=False):

    var_def = task_dict['var_def']
    for r_idx, ratiovar_str in enumerate(task_dict['var_str']):

        data  = pyfits.getdata(fits_path+'%s_%s_local_overall_counts_baseline.fits' % (var_def,ratiovar_str))
        ratio = pyfits.getdata(fits_path+'%s_%s_local_overall_ratio_baseline.fits' % (var_def,ratiovar_str))
        el = data[:,:,0].astype(np.float)
        sp = data[:,:,1].astype(np.float)
        mask_el = el < 1
        mask_sp = sp < 1
        mask_all = np.logical_and(mask_el, mask_sp)
        mask = np.logical_or(mask_el, mask_sp)
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
            im = ax.imshow(sigma_masked.T, interpolation='nearest', origin='lower')
            cb = plt.colorbar(im)
            ax.set_aspect('auto')
            ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
            ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)

        pyfits.writeto(fits_path+'%s_%s_local_overall_ratio_baseline_sigma.fits' % (var_def,ratiovar_str),
               sigma, clobber=True)    

    return None

def plot_overall_ratio_baseline(task_dict,
                        unset_vrange = False,
                        ratio_vrange = (min_ratio,max_ratio),
                        plot_mag_lims=False,
                        ):

    var_def = task_dict['var_def']
    var1_str = task_dict['var_str'][0]
    var2_str = task_dict['var_str'][1]
    bintype = task_dict['bintype']

    if task_dict['ratio_type'] is 'linear':
        label_prefix = ''
    else:
        label_prefix = 'log_{10}'

    titlenames = ('counts','ratio')

    fig = plt.figure(3,(14,8))
    fig.clf()

    # Loop over all variables

    ratiovar_strings = (var1_str,var2_str)
    for r_idx, ratiovar_str in enumerate(ratiovar_strings):

        f_ratio = fits_path+'%s_%s_local_overall_ratio_baseline.fits' % (var_def,ratiovar_str)
        p_ratio = pyfits.open(f_ratio)

        if bintype is 'counts':
            f_bins = fits_path+'%s_%s_binned_%s.fits' % (var_def,var1_str,bintype)
            p_bins = pyfits.open(f_bins)

        if bintype is 'rawlikelihood':
            f_bins = fits_path+'%s_binned_%s.fits' % (var_def,bintype)
            p_bins = pyfits.open(f_bins)

        zbins = p_bins['REDSHIFT_BIN_CENTERS'].data['centers']
        magbins = p_bins['MR_BIN_CENTERS'].data['centers']
        sizebins = p_bins['R50_KPC_BIN_CENTERS'].data['centers']

        zedges = p_bins['REDSHIFT_BIN_EDGES'].data['edges']
        magedges = p_bins['MR_BIN_EDGES'].data['edges']
        sizeedges = p_bins['R50_KPC_BIN_EDGES'].data['edges']

        ratio_baseline_masked = pickle.load(open(pkl_path+'%s_%s_local_overall_ratio_baseline_masked.pkl' % (var_def,ratiovar_str),'rb'))
        counts_baseline_masked = pickle.load(open(pkl_path+'%s_%s_local_overall_counts_baseline_masked.pkl' % (var_def,ratiovar_str),'rb'))
        sumcounts = np.sum(counts_baseline_masked, axis=2)

        if unset_vrange:
            ratio_vrange = (np.min(ratio_baseline_masked),np.max(ratio_baseline_masked))

        labelnames = (r'$N_{%s} + N_{%s}$' % (var1_str,var2_str),r'$%s(N_{%s}/(N_{%s}+N_{%s}))$' % (label_prefix,ratiovar_str,var1_str,var2_str))

        for index, data in enumerate((sumcounts, ratio_baseline_masked)):
            ax = fig.add_subplot(len(ratiovar_strings),2,index+1+(r_idx*len(ratiovar_strings)),aspect=1)
            cmap = cm.jet
            cmap.set_bad('k')
            vrange = [(0,np.max(sumcounts)),ratio_vrange]
            imextent=(magedges[0],magedges[-1],sizeedges[0],sizeedges[-1])
            im = ax.imshow(data.T, 
                           vmin = vrange[index][0],
                           vmax = vrange[index][1],
                           extent=imextent,
                           interpolation='nearest',
                           origin='lower'
                           )
            ax.set_title('%s/%s ' % (var1_str,var2_str)+titlenames[index])
            ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
            ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)
            ax.set_aspect('auto')
            rc(('xtick','ytick'), labelsize=12)
            cb = plt.colorbar(im,orientation='vertical')
            cb.set_label(labelnames[index],fontsize=16)

            SBlim_mag = (SBapp - cosmology.dmod_flat(np.mean(zbins))- 2.5*np.log10(6.283185*(SBlim_size/cosmology.ang_scale_flat(np.mean(zbins)))**2))
            absmag_lim = appmag_lim - cosmology.dmod_flat(np.mean(zbins))
            absmag_lim_loz = appmag_lim - cosmology.dmod_flat(0.0005)
            absmag_lim_hiz = appmag_lim - cosmology.dmod_flat(0.25)
            r50_lim = cosmology.ang_scale_flat(np.mean(zbins))
            ax.autoscale(False)
            ax.plot(SBlim_mag, SBlim_size,'w--')
            if plot_mag_lims:
                ax.axvline(absmag_lim,color='w',linestyle='dashed')
                ax.axvline(absmag_lim_loz,color='g',linestyle='dashed')
                ax.axvline(absmag_lim_hiz,color='g',linestyle='dashed')

            ax.axhline(r50_lim, color='w', linestyle='dashed')
            #plt.show()
            #plt.draw()

        fig.savefig(plots_path+'%s_overall_ratio_baseline.png' % var_def, dpi=200)
        p_ratio.close()
        p_bins.close()

    return None

