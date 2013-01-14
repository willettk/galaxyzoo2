# astrotools.py

# Useful astronomy utility functions

# Steven Bamford
# Created September 2005

from math import *
import numpy as N
#from ppgplot_spb import *
from checkarray import checkarray
from cosmology import *
import string
import sys
import copy

# convert strings of RA (hh:mm:ss.s) and Dec (dd:mm:ss.s) in to
# decimal degrees (d.ddd)
def radec_to_deg(ra, dec):
    ra_h = float(ra[: ra.find(':')])
    ra_m = float(ra[ra.find(':')+1 : ra.rfind(':')])
    ra_s = float(ra[ra.rfind(':')+1 :])
    dec_d = float(dec[: dec.find(':')])
    dec_m = float(dec[dec.find(':')+1 : dec.rfind(':')])
    dec_s = float(dec[dec.rfind(':')+1 :])
    ra_deg, dec_deg = rahms_decdms_to_deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)
    return ra_deg, dec_deg


def radec_deg_to_hms_dms(ra_deg, dec_deg, sec_precision=3, safe=True,
                         asstring=True):
    if safe:
        ra_deg = checkarray(ra_deg, double=True)
        dec_deg = checkarray(dec_deg, double=True)
    ra_abs = N.absolute(ra_deg)
    ra_signval = ra_deg/ra_abs
    ra_sign = N.where(ra_signval > 0, '', '-')
    ra_hfull = ra_abs / 15.0
    ra_h = ra_hfull.astype(N.int)
    ra_mfull = (ra_hfull - ra_h) * 60.0
    ra_m = ra_mfull.astype(N.int)
    ra_s = (ra_mfull - ra_m) * 60.0
    dec_abs = N.absolute(dec_deg) 
    dec_signval = dec_deg/dec_abs
    dec_sign = N.where(dec_signval > 0, '+', '-')
    dec_d = dec_abs.astype(N.int)
    dec_mfull = (dec_abs - dec_d) * 60.0
    dec_m = dec_mfull.astype(N.int)
    dec_s = (dec_mfull - dec_m) * 60.0
    ra_sec_precision = max(0,sec_precision-1)
    dec_sec_precision = max(0,sec_precision)
    ra_format = '%s%02i:%02i:%'+'0%i.%i'%(3+ra_sec_precision, ra_sec_precision)+'f'
    dec_format = '%s%02i:%02i:%'+'0%i.%i'%(3+dec_sec_precision, dec_sec_precision)+'f'
    n = len(ra_deg)
    if asstring:
        ra = N.zeros(n, dtype='S%i'%(10+sec_precision))
        dec = N.zeros(n, dtype='S%i'%(10+sec_precision))
        for i in range(n):
            ra[i] = ra_format%(ra_sign[i], ra_h[i], ra_m[i], ra_s[i])
            dec[i] = dec_format%(dec_sign[i], dec_d[i], dec_m[i], dec_s[i])
            if n == 1:
                ra = ra[0]
                dec = dec[0]
        return ra, dec
    else:
        return ra_signval*ra_h, ra_m, ra_s, dec_signval*dec_d, dec_m, dec_s

def rahms_decdms_to_deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, safe=True):
    if safe:
        ra_h, ra_m, ra_s = [checkarray(i, double=True)
                            for i in (ra_h, ra_m, ra_s)]
        dec_d, dec_m, dec_s = [checkarray(i, double=True)
                               for i in (dec_d, dec_m, dec_s)]
    ra_sign = (ra_h+0.5) / abs(ra_h+0.5)
    ra_h = N.absolute(ra_h)
    ra_hours = ra_sign * (ra_h + ra_m / 60.0 + ra_s / 3600.0)
    ra_deg = ra_hours * 15.0
    dec_sign = (dec_d+0.5) / abs(dec_d+0.5)
    dec_d = N.absolute(dec_d)
    dec_deg = dec_sign * (dec_d + dec_m / 60.0 + dec_s / 3600.0)
    # put onto ranges RA: 0..360 deg, Dec: -90..90
    while (dec_deg >= 360.0).any():
        dec_deg[decdeg >= 360.0] -= 360.0
    while (dec_deg < 0.0).any():
        dec_deg[decdeg < 0.0] += 360.0
    select = (dec_deg < 270.0)&(dec_deg > 90.0)
    ra_deg[select] += 180.0
    dec_deg[select] *= -1
    dec_deg[select] += 180
    select = (dec_deg > -270.0)&(dec_deg < -90.0)
    ra_deg[select] += 180.0
    dec_deg[select] *= -1
    dec_deg[select] -= 180
    select = (dec_deg >= 270.0)
    dec_deg[select] -= 360
    select = (dec_deg <= -270.0)
    dec_deg[select] += 360
    while (ra_deg >= 360.0).any():
        ra_deg[radeg >= 360.0] -= 360.0
    while (ra_deg < 0.0).any():
        ra_deg[radeg < 0.0] += 360.0
    if len(ra_deg) == 1:
        ra_deg = ra_deg[0]
        dec_deg = dec_deg[0]
    return ra_deg, dec_deg


# calculate angle between two co-planar angles in decimal degrees
def calc_ang_diff(a1, a2):
    # Put angles on range 0 -- 360
    if a1 >= 0.0:
        while a1 >= 360.0:  a1 = a1 - 360.0
    else:
        while a1 <= -360.0:  a1 = a1 + 360.0
    if a2 >= 0.0:
        while a2 >= 360.0:  a2 = a2 - 360.0
    else:
        while a2 <= -360.0:  a2 = a2 + 360.0
    # numerical difference
    delta = a2 - a1
    # put difference on range -180 -- 180
    if delta > 180.0:
        delta = delta - 360.0
    elif delta <= -180.0:
        delta = 360.0 + delta
    return delta


def calc_ang_dist(ra1, dec1, ra2, dec2, units='degrees',
                  safe=True, keeparray=False):
    DPIBY2 = 0.5 * pi
    
    if safe:
        ra1, dec1, ra2, dec2 = [checkarray(i, double=True)
                                for i in (ra1, dec1, ra2, dec2)]
    
    if units == "hrdeg":
        convRA = pi / 12.0
        convDEC = pi / 180.0
    elif units == "radians":
        convRA = 1.0
        convDEC = 1.0
    else: # degrees
        convRA = pi / 180.0
        convDEC = pi / 180.0
        
    theta1 = dec1*convDEC + DPIBY2
    theta2 = dec2*convDEC + DPIBY2
    cosgamma = (N.sin(theta1) * N.sin(theta2) * N.cos((ra1-ra2)*convRA) + 
               N.cos(theta1) * N.cos(theta2))
    adist = 0.0 * cosgamma
    ivalid = (cosgamma < 1.0).nonzero()[0]
    if len(ivalid) > 0:
        adist[ivalid] = N.arccos(cosgamma[ivalid]) / convDEC
    if (not keeparray) and (adist.shape == (1,)):
        return adist[0]
    else:
        return adist


def calc_ang_dist_old(ra1, dec1, ra2, dec2):
    # approx
    print 'APPROX - get djs_diff_angle from idlutils for correct calculation'
    delta_dec = dec2 - dec1
    delta_ra = (ra2 - ra1) * N.cos((pi/180)*(dec2+dec1)/2.0)
    delta = N.sqrt(delta_dec**2 + delta_ra**2)
    delta = delta * 3600.0
    return delta

def angles_to_xyz(r,phi,theta):
# Convert spherical coords (r,phi,theta) into Cartesion coords (x,y,z).

#  The angles must be in the following ranges:
#    0 <= phi < 360
#    0 <= theta <= 180
#  where theta=0 corresponds to the N pole, and theta=180 is the S pole.
#  If you want to convert from RA and DEC, pass the following
#  arguments (in degrees):  RA, 90-DEC
#  From IDLUTILS.
 
   dradeg = 180.0/pi
   stheta = N.sin(theta / dradeg)
   x = r * N.cos(phi / dradeg) * stheta
   y = r * N.sin(phi / dradeg) * stheta
   z = r * N.cos(theta / dradeg)
 
   return x, y, z

def radec_to_xy(ra, dec, racentre=None, deccentre=None, rotate=None):
    dradeg = 180.0/pi
    cdelt = N.ones(2, N.float)/60.0
    if racentre is None:
        racentre = N.mean(ra)
    if deccentre is None:
        deccentre = N.mean(dec)
    racentre /= dradeg
    deccentre /= dradeg
    radif = ra/dradeg - racentre
    dec = dec / dradeg
    h = N.sin(dec)*N.sin(deccentre) + N.cos(dec)*N.cos(deccentre)*N.cos(radif)
    xsi = N.cos(dec)*N.sin(radif)/h
    eta = (N.sin(dec)*N.cos(deccentre) - N.cos(dec)*N.sin(deccentre)*N.cos(radif))/h
    xsi = xsi*dradeg
    eta = eta*dradeg
    xsi = xsi/cdelt[0]
    eta = eta/cdelt[1]

    if rotate is not None:
        rotate /= dradeg
        cd = N.array([[N.cos(rotate), -N.sin(rotate)],
                      [N.sin(rotate), N.cos(rotate)]])
        x = cd[0,0]*xsi + cd[0,1]*eta
        y = cd[1,0]*xsi + cd[1,1]*eta
    else:
        x = xsi
        y = eta
    
    return x, y

def distance(x1, y1, z1, x2, y2, z2):
    return N.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def _normdistance(gra, gdec, gredshift, cra, cdec, credshift, rvir, deltaz):
    # This function is not used - it has been included directly in the
    # dr6_sample.py code to speed up.
    angscale = ang_scale_flat(credshift) # kpc/arcsec
    angscale *= 3.6 # Mpc/deg
    angdist = calc_ang_dist(gra, gdec, cra, cdec) # degrees
    angdist *= angscale # Mpc
    angdist /= rvir # virial radii
    losdist = dC_flat(gredshift) - dC_flat(credshift) # Mpc
    # adjust losdist to be smaller when redshifts are similar,
    # but unchanged otherwise.
    # Fix losdist=rvir when |gredshift-credshift| < deltaz
    # and ajustment < 1 per cent for |gredshift-credshift| > 2*deltaz
    A = -N.log(1 - rvir/(dC_flat(credshift+deltaz/2.0) - dC_flat(credshift-deltaz/2.0)))
    B = N.log(N.log(100)/A)/N.log(2.0)
    C = 1 - N.exp(-A * N.abs(gredshift - credshift)**B / (deltaz**B))
    losdist *= C
    losdist /= rvir # virial radii
    dist = N.sqrt(angdist**2 + losdist**2)
    return dist

def Tfunction(p, x):
    p0, p1, q0, q1, q2 = p
    return p0 + p1*x + q0*N.tanh((x-q1)/q2)

def Tfunction04(p, x):
    p0, p1, q0, q1, q2 = p
    return p0 + p1*(x+20) + q0*N.tanh((x-q1)/q2)


def lumfn_blue_baldry04(Mr):
    # Schechter parameters for blue sequence
    #-20.90  0.00254  -0.25  0.00157  -1.54
    #  0.08  0.00032   0.23  0.00043   0.08
    pass

def lumfn_red_baldry04(Mr):
    # Schechter parameters for red sequence
    #-21.34  0.00304  -0.69  0.00000  -1.00
    #0.02  0.00008   0.02 -9.99999  -9.99
    pass

def Cur_divide_baldry06(Mstar):
    p = (2.18, 0, 0.38, 10.26, 0.85)
    return Tfunction(p, Mstar)

def Cur_divide_baldry04(Mr):
    # updated to DR4
    # http://www.astro.ljmu.ac.uk/~ikb/research/bimodality-paperI.html
    p = (2.192, 0, -0.354, -20.58, 2.05)
    return Tfunction(p, Mr)

def Cur_red_baldry04(Mr):
    # updated to DR4
    # http://www.astro.ljmu.ac.uk/~ikb/research/bimodality-paperI.html
    p = (2.386, -0.079, -0.076, -19.75, 1.04)
    # old, from paper
    #p = (2.279, -0.037, -0.108, -19.81, 0.96)
    return Tfunction04(p, Mr)

def Cur_red_baldry04_Mstar(Mstar):
    Cur = Cur_divide_baldry06(Mstar)
    deltaCur = 1.0
    oldCur = Cur 
    while deltaCur > 0.001:
        Mr = Mr_from_logMstar_baldry06(Mstar, Cur)
        Cur = Cur_red_baldry04(Mr)
        deltaCur = (N.absolute(Cur - oldCur)).max()
        oldCur = Cur
    return Mr, Cur

def Cur_blue_baldry04_Mstar(Mstar):
    Cur = Cur_divide_baldry06(Mstar)
    deltaCur = 1.0
    oldCur = Cur 
    while deltaCur > 0.001:
        Mr = Mr_from_logMstar_baldry06(Mstar, Cur)
        Cur = Cur_blue_baldry04(Mr)
        deltaCur = (N.absolute(Cur - oldCur)).max()
        oldCur = Cur
    return Mr, Cur

def Cur_sigma_red_baldry04(Mr):
    # updated to DR4
    # http://www.astro.ljmu.ac.uk/~ikb/research/bimodality-paperI.html
    p = (0.140, 0.000, 0.049, -21.23, 1.13)
    # old, from paper
    #p = (0.152, 0.008, 0.044, -19.91, 0.94)
    return Tfunction04(p, Mr)

def Cur_blue_baldry04(Mr):
    # updated to DR4
    # http://www.astro.ljmu.ac.uk/~ikb/research/bimodality-paperI.html
    p = (1.866, -0.062, -0.361, -21.29, 1.54)
    # old, from paper
    #p = (1.790, -0.053, -0.363, -20.75, 1.12)
    return Tfunction04(p, Mr)

def Cur_sigma_blue_baldry04(Mr):
    # updated to DR4
    # http://www.astro.ljmu.ac.uk/~ikb/research/bimodality-paperI.html
    p = (0.276, 0.000, -0.038, -19.88, 0.43)
    # old, from paper
    #p = (0.298, 0.014, -0.067, -19.90, 0.58)
    return Tfunction04(p, Mr)

def check_Cur_divide_baldry06():
    pgaqt()
    pgsetup()
    pgenv(8.85, 11.65, 0.0, 3.5)
    pglab('log(\(2563)/\(2563)\d\(2281)\u)', 'restframe (u-r)\dmodel\u', '')
    Mstar = N.arange(8.8, 11.8, 0.1)
    Cur = Cur_divide_baldry06(Mstar)
    pgline(Mstar, Cur)
    pgend()

def logML_baldry06(Cur):
    MLa = -0.95 + 0.56*Cur
    MLb = -0.16 + 0.18*Cur
    ML = N.where(MLa < MLb, MLa, MLb)
    return ML

def check_logML_baldry06():
    pgaqt()
    pgsetup()
    pgenv(0.5, 3.5, -1.0, 1.0)
    pglab('restframe (u-r)\dmodel\u', 'log(\(2563)/L\dr\u) (solar units)', '')
    Cur = N.arange(0.5, 3.6, 0.1)
    ML = logML_baldry06(Cur)
    pgline(Cur, ML)
    pgend()

def logMstar_baldry06(Mr, Cur):
    Mr_solar = 4.62
    logML = logML_baldry06(Cur)
    logM = (Mr_solar - Mr)/2.5 + logML
    return logM

def Mr_from_logMstar_baldry06(Mstar, Cur):
    Mr_solar = 4.62
    logML = logML_baldry06(Cur)
    Mr = Mr_solar - 2.5*(Mstar - logML)
    return Mr

# Utilties for matching catalogs by RA & Dec or ID
# H. Ferguson 11/22/05

def matchsorted(ra,dec,ra1,dec1,tol):
    """ Find closest ra,dec within tol to a target in an ra-sorted list of ra,dec.
        Arguments:
          ra - Right Ascension decimal degrees (numpy sorted in ascending order)
          dec - Declination decimal degrees (numpy array)
          ra1 - RA to match (scalar, decimal degrees)
          ra1 - Dec to match (scalar, decimal degrees)
          tol - Matching tolerance in decimal degrees. 
        Returns:
          ibest - index of the best match within tol; -1 if no match within tol
          sep - separation (defaults to tol if no match within tol)
    """
    i1 = N.searchsorted(ra,ra1-tol)-1
    i2 = N.searchsorted(ra,ra1+tol)+1
    if i1 < 0:
        i1 = 0
    sep = calc_ang_dist(ra[i1:i2],dec[i1:i2],ra1,dec1, units='degrees',
                        safe=False, keeparray=True)
    indices = N.argsort(sep)
    if sep[indices[0]] > tol:
        return -1, tol
    ibest = indices[0] + i1
    return ibest, sep[indices[0]]
    
def matchpos(ra1,dec1,ra2,dec2,tol):
    """ Match two sets of ra,dec within a tolerance.
        Longer catalog should go first
        Arguments:
          ra1 - Right Ascension decimal degrees (numpy array)
          dec1 - Declination decimal degrees (numpy array)
          ra2 - Right Ascension decimal degrees (numpy array)
          dec2 - Declination decimal degrees (numpy array)
          tol - Matching tolerance in degrees. 
        Returns:
          ibest - indices of the best matches within tol; -1 if no match within tol
          sep - separations (defaults to tol if no match within tol)
    """
    ra1, dec1, ra2, dec2 = [checkarray(i, double=True)
                            for i in (ra1, dec1, ra2, dec2)]
    indices = N.argsort(ra1)
    rasorted = ra1[indices]
    decsorted = dec1[indices]
    ibest = []
    sep = []
    m = len(ra1)
    n = len(ra2)
    print '%i rows in catalog 1,  %i rows in catalog 2'%(m, n)
    if n > 1000:
        npc = n//100
        print '%:',
    for i in range(n):
        if n > 1000:
            if i%npc == 0:
                print i//npc,
                sys.stdout.flush()
        j, s = matchsorted(rasorted,decsorted,ra2[i],dec2[i],tol)
        if j < 0:
            ibest += [j]
        else:
            ibest += [indices[j]]
        sep += [s]
    if n > 1000: print
    return N.array(ibest), N.array(sep)

def matchjoin(si1,si2,ibest,sep=[],dict1={},dict2={}):
    """ Keep only elements that match in both catalogs. 
        Arguments:
          si1 -- object with data arrays as attributes (e.g. si1.mag_auto, si1.id)
          si2 -- object with data arrays as attributes 
          ibest -- indices of si1 that match si2 (in order of si2)
        Keyword Arguments:
          sep -- array of separations, will be added as an attribute to second object
          dict1 -- dictionary of attributes for si1 (e.g. {'mag_auto':0,'id':1} )
          dict2 -- dictionary of attributes for si2 
        The objects si1 and si2 would normally be sextractor objects returned
        by sextutils.sextractor(). In that case, the column names are stored in
        si1._d and si2._d. If the objects are not sextractor objects, the user can
        provide separate dictionaries whose keys are the object attributes 
	that correspond to the numpy data arrays. 
       
        Returns:
          s1, s2 -- object arrays that include only the matching objects contained
              in both s1 and s2.
    """
    if len(dict1) == 0:
        dict1 = si1._d
    if len(dict2) == 0:
        dict2 = si2._d
    indices = N.compress(ibest > 0,ibest)
    # print indices
    s1 = copy.deepcopy(si1)
    s2 = copy.deepcopy(si2)
    for k in dict1.keys():
        if type(s1.__dict__[k]) == type([]):
             s1.__dict__[k] = char.array(s1.__dict__[k])
        s1.__dict__[k] = s1.__dict__[k][indices]
    flag = ibest > 0
    for k in dict2.keys():
        if type(s2.__dict__[k]) == type([]):
             s2.__dict__[k] = char.array(s2.__dict__[k])
        s2.__dict__[k] = N.compress(flag,s2.__dict__[k])
    if len(sep) > 0:
        s2.separation = N.compress(flag,sep)
    return s1,s2

def matchids(id1,id2):
    """ Match two sets of ids. 
        Returns: 
          ibest -- array of indices of i1 that match i2; -1 if no match
    """
    indices = N.argsort(id1)
    idsorted = id1[indices]
    ibest = []
    print len(id2)
    for i in range(len(id2)):
        print i
        j = matchidsorted(idsorted,id2[i])
        if j < 0:
            ibest += [j]
        else:
            ibest += [indices[j]]
    return N.array(ibest)

def matchidsorted(ids,targetid):
    """ Find id matches, return index in i1 that matches targetid; -1 if no match. """
    i1 = N.searchsorted(ids,targetid)
    if targetid == ids[i1]:
        ibest = i1
    else:
        ibest = -1 
    return ibest
