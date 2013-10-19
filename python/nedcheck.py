from random import sample 
from galaxyzoo2 import get_cat

# Select 3000 random galaxies from each sample

specz = get_cat()
photoz = get_cat(photoz=True)
stripe82 = get_cat(stripe82=True,depth='normal') 
coadd1 = get_cat(stripe82=True,depth='coadd1') 
coadd2 = get_cat(stripe82=True,depth='coadd2') 

samples = (specz,photoz,stripe82,coadd1,coadd2)
samplestrings = ('specz','photoz','stripe82','coadd1','coadd2')

for samp,sampstr in zip(samples,samplestrings):

    batchtemplatefile = '/Users/willettk/Astronomy/general/batch_sim.txt'
    b1 = open(batchtemplatefile,'rb')
    b = b1.readlines()
    b1.close()
    
    writedir = '/Users/willettk/Astronomy/Research/GalaxyZoo/nedcheck/'
    obj = open('%swillettk_gz2_%s.txt' % (writedir,sampstr), 'wb')
    
    for line in b[:30]:
        obj.write(line)
    
    obj.write('OUTPUT_FILENAME    willettk_gz2_%s\n' % sampstr)
    
    for line in b[32:124]:
        obj.write(line)
    
    cootemp = sample(zip(samp['rastring'],samp['decstring']),3000)
    for c in cootemp:
        obj.write('SDSS J%s%s\n' % (c[0].translate(None,':'),c[1].translate(None,':')))
        
    for line in b[125:]:
        obj.write(line)
    
    obj.close()


