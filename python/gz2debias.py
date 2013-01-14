from math import *
import numpy as np
import pyfits
import cosmology
import sys
from matplotlib.pyplot import *

import gc
gc.collect()

data_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/fits/'
gz2table_file = 'gz2table.fits'
gz2sample_file = 'gz2sample.fits'
gz2ts = 'gz2_table_sample_match_coadd2.fits'

def plot_a2():

	#gz2table = pyfits.open(data_path + gz2table_file)[1]
	#gz2sample = pyfits.open(data_path + gz2sample_file)[1]
	g = pyfits.open(data_path + gz2ts)
	gz2 = g[1]
	p_el = gz2.data['t01_smooth_or_features_a01_smooth_weighted_fraction']
	#p_sp = gz2.data['t01_smooth_or_features_a02_features_or_disk_weighted_fraction']
	redshift = gz2.data['redshift']

	plot(redshift,p_el,'r.')

	#gz2table.close()
	#gz2sample.close()
	g.close()

