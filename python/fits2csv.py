import pyfits

# This has now been profiled and made around 40 times faster than before!

def fits2csv(infilename, outfilename, sep=',', linesep='\n', maxlength=32,
             selectednames=None, limit=None, gzipped=False, keeporder=True,
             blankfill=None):
    d = pyfits.getdata(infilename)
    rec2csv(d, outfilename, sep, linesep, maxlength,
            selectednames, limit, gzipped, keeporder, blankfill)

def rec2csv(d, outfilename, sep=',', linesep='\n', maxlength=32,
            selectednames=None, limit=None, gzipped=False, keeporder=True,
            blankfill=None):
    if limit is not None:
        d = d[eval(limit)]
    if selectednames is None:
        names = d.names
    else:
        if keeporder:
            names = [n for n in d.names if n in selectednames]
        else:
            names = selectednames
    columns = [(d.field(n)).astype('S%i'%maxlength) for n in names]
    if gzipped:
        import gzip
        fout = gzip.open(outfilename, 'wb')
    else:
        fout = file(outfilename, 'w')
    fout.write(sep.join(names)+linesep)
    for row in range(len(d)):
	s = [c[row] for c in columns]
        if blankfill is not None:
            for i, si in enumerate(s):
                if si.strip() == '':
                    s[i] = blankfill
	fout.write(sep.join(s)+linesep)
    fout.close()

def fits2csv_round(infilename, outfilename, sep=',', linesep='\n', maxlength=32,
                   selectednames=None, limit=None, round=4, gzipped=False):
    d = pyfits.getdata(infilename)
    if limit is not None:
        d = d[eval(limit)]
    if selectednames is None:
        selectednames = d.names
    floatnames = [n[0] for n in d.dtype.descr
                  if (n[0] in selectednames) and ('f' in n[1])]
    othernames = [n[0] for n in d.dtype.descr
                  if (n[0] in selectednames) and ('f' not in n[1])]
    floatcolumns = dict([(c, d.field(c)) for c in floatnames])
    othercolumns = dict([(c, d.field(c)) for c in othernames])
    n = len(d)
    del d
    if gzipped:
        import gzip
        fout = gzip.open(outfilename, 'wb')
    else:
        fout = file(outfilename, 'w')
    fout.write(sep.join(selectednames)+linesep)
    roundformat = '%.' + '%if'%round
    for row in range(n):
        s = []
        for sn in selectednames:
            if sn in floatnames:
                s.append(roundformat%(floatcolumns[sn][row]))
            else:
                s.append('%s'%(othercolumns[sn][row]))
	fout.write(sep.join(s)+linesep)
    fout.close()

def fits2latex(infilename, outfilename, maxlength=32,
               selectednames=None, limit=None, round=3):
    fits2csv_round(infilename, outfilename,
                   sep=' & ', linesep='\\\\\n', maxlength=maxlength,
                   selectednames=selectednames, limit=limit, round=round)
