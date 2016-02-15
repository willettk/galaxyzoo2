from astropy.io import fits as pyfits
import numpy as N
import os.path

def combine_fits_tables(tables, outfile, idfield='OBJID'):
    cols = []
    matches = {}
    lengths = []
    for t in tables:
        columns = t[1].columns
        data = t[1].data
        lengths.append(len(data))
        uniqidx, uniqid = N.unique1d(data.field(idfield), True)
        repeatedidx = N.ones(len(data), N.bool)
        repeatedidx[uniqidx] = False
        repeatedid = data.field(idfield)[repeatedidx]
        if len(repeatedid) > 0:
            print 'Repeated ids:', repeatedid
            print 'The repeats have been discarded'
            data = data[uniqidx]
        order = N.argsort(data.field(idfield))
        for c in columns:
            name = c.name.upper()
            if name not in [ci.name for ci in cols]:
                cols.append(pyfits.Column(name=name, format=c.format,
                                          array=data.field(name)[order]))
            else:
                matches[name] = data.field(name)[order]
    print 'Table lengths, with repeats:',
    print lengths
    tbhdu=pyfits.new_table(cols)
    data = tbhdu.data
    # check that repeated columns are identical
    print 'Repeated columns:',
    print matches.keys()
    for name in matches.keys():
        if N.any(data.field(name) != matches[name]):
            print 'Warning!  Repeated columns (%s) differ.'%name
            print data.field(name)[:5]
            print matches[name][:5]
    file_exists = os.path.isfile(outfile)
    if file_exists:
        os.remove(outfile)
    tbhdu.writeto(outfile)

def concatenate_fits_tables(tables, outfile, uniquecol=None):
    tmaster = tables[0]
    nhdu = len(tmaster)
    hdus = [pyfits.PrimaryHDU()]
    for t in tables[1:]:
        nhdut = len(t)
        if nhdu != nhdut:
            raise ValueError('files do not all have same number of hdus')
    for i in range(nhdu):
        nrows1 = len(tmaster[i].data)
        nrows = nrows1
        ncols1 = len(tmaster[i].columns.names)
        print 'ncols1:', ncols1
        names = tmaster[i].columns.names
        for t in tables[1:]:
            nrowst = len(t[i].data)
            ncolst = len(t[i].columns.names)
            print 'ncolst:', ncolst
            namest = t[i].columns.names
            newnames = []
            for n in names:
                if n in namest:
                    newnames.append(n)
            names = newnames
            #if ncols1 != ncolst:
            #    raise ValueError('tables do not all have same number of columns')
            nrows += nrowst
        cols = []
        for j in range(len(tmaster[i].columns)):
            if tmaster[i].columns.names[j] in names:
                cols.append(tmaster[i].columns[j])
        hdu = pyfits.new_table(cols, nrows=nrows)
        print names
        startrow = nrows1
        for t in tables[1:]:
            nrowst = len(t[i].data)
            for name in names:
                hdu.data.field(name)[startrow:startrow+nrowst]=t[i].data.field(name) 
            startrow += nrowst
        if uniquecol is not None:
            ar, idx = N.unique(hdu.data.field(uniquecol), return_index=True)
            hdu.data = hdu.data[idx]
        hdus.append(hdu)
    hdulist = pyfits.HDUList(hdus)
    hdulist.info()
    hdulist.writeto(outfile)
