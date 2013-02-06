from astropy.io import fits as pyfits
import numpy as np
from galaxyzoo2 import fits_path

gz2_full_data_file = fits_path+'gz2_original_extra_s82norm_table_sample.fits'
p = pyfits.open(gz2_full_data_file)
gzdata = p[1].data
p.close()

pdisk = gzdata['t01_smooth_or_features_a02_features_or_disk_weighted_fraction']
pnotedge = gzdata['t02_edgeon_a05_no_weighted_fraction']
nbar = gzdata['t03_bar_a06_bar_weight']
pbar = gzdata['t03_bar_a06_bar_weighted_fraction']
redshift = gzdata['REDSHIFT']

klm = ((pdisk*pnotedge) >= 0.25) & (nbar >= 10) & (pbar >= 0.5) & (redshift <= 0.01)


col1 = pyfits.Column(name = 'objid', format='A18', array=gzdata['objid_1'][klm])
col2 = pyfits.Column(name = 'p_bar', format='E', array=pbar[klm])
col3 = pyfits.Column(name = 'petror50_r', format='E', array=gzdata['PETROR50_R'][klm])
col4 = pyfits.Column(name = 'petror90_r', format='E', array=gzdata['PETROR90_R'][klm])
col5 = pyfits.Column(name = 'p_disk', format='E', array=pdisk[klm])
col6 = pyfits.Column(name = 'p_notedge', format='E', array=pnotedge[klm])
col7 = pyfits.Column(name = 'n_bar', format='E', array=nbar[klm])
col8 = pyfits.Column(name = 'redshift', format='E', array=redshift[klm])

primary_hdu = pyfits.PrimaryHDU(np.array((10,10)))
hdulist = pyfits.HDUList([primary_hdu])

tb1_hdu = pyfits.new_table(pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8]))
tb1_hdu.name = 'GZDATA'

hdulist.append(tb1_hdu)

hdulist.writeto(fits_path+'klm_evla_bars_gas.fits',clobber=True)    

"""
Remaining steps to send data to Karen:

In TOPCAT, convert the `objid` column to a long integer by using parseLong($1)
Export file from TOPCAT as VOTable
Upload VOTable to CasJobs SDSS-II site
Cross-match data using SQL query:

    SELECT KLM.*, G.modelmag_r, log10(G.expAB_r) as log_expAB_r, G.ra, G.dec
    into mydb.klm_all 
    from mydb.klm2 as KLM
    JOIN GALAXY as G on G.objid = KLM.objid

Download data as FITS or CSV file.
"""
