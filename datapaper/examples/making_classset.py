In [848]: specz = gz2.get_cat()

In [855]: classset = set(specz['gz2class'])

In [862]: ll = []

In [864]: for item in classset:
    ll.append((np.sum(specz['gz2class'].count(item)),item))
   .....:     

In [882]: mostcommon = []

In [883]: for l in ll:
    if l[0] > 5000:
        mostcommon.append(l[1])
   .....:         


In [890]: for mc in mostcommon:
    print mc,'http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%6.3f&dec=%7.4f&scale=0.30&opt=&width=256&height=256' % (specz['ra'][(specz['gz2class'] == mc) & (specz['redshift'] > 0.05) & (specz['redshift'] < 0.055)][0],specz['dec'][(specz['gz2class'] == mc) & (specz['redshift'] > 0.05) & (specz['redshift'] < 0.055)][0])
   .....:     
SbB http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=163.883&dec= 8.8470&scale=0.30&opt=&width=256&height=256
Sc http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=179.787&dec=45.9913&scale=0.30&opt=&width=256&height=256
Sb http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=131.207&dec=43.4505&scale=0.30&opt=&width=256&height=256
Sd http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=173.222&dec=57.1641&scale=0.30&opt=&width=256&height=256
Ei http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=116.578&dec=18.3687&scale=0.30&opt=&width=256&height=256
Ec http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=226.239&dec=16.3398&scale=0.30&opt=&width=256&height=256
Er http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=175.197&dec=46.5405&scale=0.30&opt=&width=256&height=256
Sc2m http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=231.997&dec=51.6779&scale=0.30&opt=&width=256&height=256
Ser http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=129.596&dec=28.6489&scale=0.30&opt=&width=256&height=256
Sen http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=198.768&dec= 8.4868&scale=0.30&opt=&width=256&height=256
Sc?t http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=323.174&dec=-0.1278&scale=0.30&opt=&width=256&height=256
ScB http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=218.730&dec= 8.8535&scale=0.30&opt=&width=256&height=256

