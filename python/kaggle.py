import numpy as np
import random
import math
from astropy.io import fits as pyfits
from astropy.io import ascii
from astropy.table import Table

def split_datasets():

    # Create the test vs. training datasets for the Kaggle GZ2 competition, based
    # on random selection of GZ2 galaxies according to Steven Bamford's weighted Voronoi tessellation
    # by redshift, size, and luminosity. 

    p = pyfits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/kaggle/kaggle_gz2_wvt.fits')
    wvt = p[1].data
    p.close()

    # Proceed individually in each redshift bin

    z_good = (wvt['REDSHIFT_SIMPLE_BIN'] >= 0)
    wvt_good = (wvt['WVT_BIN'] > 0)
    z_wvt = z_good & wvt_good

    wvt_reduced = wvt[z_wvt]

    unique_redshift_bins = np.unique(wvt_reduced['REDSHIFT_SIMPLE_BIN'])

    # Goal size for test sample
    n_test = 80000
    ngal = len(wvt_reduced)

    # Test sample should be 75% private, 25% public
    f_test = n_test / float(ngal)
    f_training = 1. - f_test
    f_test_private = 0.75
    f_test_public = 0.25

    private = []
    public = []
    training = []

    # Loop over redshift bins, since WVT bin naming is not unique

    for zb in unique_redshift_bins:
        wvt_zbin = wvt_reduced[wvt_reduced['REDSHIFT_SIMPLE_BIN'] == zb]
        wvt_bins = np.unique(wvt_zbin['WVT_BIN'])

        # Loop over WVT bin

        for wb in wvt_bins:
            wb_temp = list((wvt_zbin[wvt_zbin['WVT_BIN'] == wb])['dr7objid'])
            lt = len(wb_temp)
            random.shuffle(wb_temp)

            slice_training = int(math.floor(lt*f_training))
            slice_private  = int(math.floor((lt - slice_training) * f_test_private) + slice_training)

            training.extend(wb_temp[:slice_training])
            private.extend(wb_temp[slice_training:slice_private])
            public.extend(wb_temp[slice_private:])


    # Save data as CSV file

    d_all = training + private + public
    s_all = ['training']*len(training)+['test_private']*len(private)+['test_public']*len(public)
    data = Table({'dr7objid': d_all,
                  'kaggle_set': s_all},
                 names=['dr7objid', 'kaggle_set'])
    ascii.write(data,'/Users/willettk/Astronomy/Research/GalaxyZoo/kaggle/kaggle_gz2_testtraining.csv',delimiter=',')

    return None
