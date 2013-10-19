import imp
import numpy as np
import datetime

gzdir = '/Users/willettk/Astronomy/Research/GalaxyZoo/'
nedpy = imp.load_source('','/Users/willettk/python/nedpy/nedpy.py')
gz2 = imp.load_source('',gzdir+'python/galaxyzoo2.py')

photoz = gz2.get_cat(photoz=True)
ra = photoz['ra']
dec = photoz['dec']

chunks = range(3000,len(ra),1000)
chunks.append(len(ra))

now = datetime.datetime.now()

for firstind,lastind in zip(chunks[:-1],chunks[1:]):
    zlist = []

    for r,d in zip(ra[firstind:lastind],dec[firstind:lastind]):
        try:
            ret = nedpy.query_ned_nearpos(ra=r,dec=d,sr=2.0)
            z = ret['Redshift']
            goodz = z[np.isfinite(z)]
            if len(goodz) > 0:
                zlist.append(goodz[0])
            else:
                zlist.append(np.nan)
        except TypeError,ValueError:
            zlist.append(np.nan)
    
        if len(zlist) % 100 == 0:
            print 'Set %5i, %5i galaxies in %3i seconds' % (firstind,len(zlist),(datetime.datetime.now() - now).seconds)
            now = datetime.datetime.now()
    
    np.save('%snedcheck/ned_photoz_%i.npy' % (gzdir,firstind),np.array(zlist))
