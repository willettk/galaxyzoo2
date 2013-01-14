# rostat.py

# median = MDIAN1, biwt = XBIWT, mad = XMAD
# converted to python by hand from rostat.f
# by Steven Bamford, 03 June 2004

# confidence interval output added to biwt from
# subroutine CONFIDENCE in rostat.f

# t-tests adapted for biweight from Numerical
# Recipes functions with dof_biwt = 0.7*dof

# fit by minimising MAD adapted from
# Numerical Recipes

# ks-test adapted from Numerical Recipes

import numpy as N
import nr
from math import sqrt, exp
import pickle
import sys
from os.path import join as pjoin

def median(x, sorted=0):
    if not sorted: x = N.sort(x,axis=0)
    n = len(x)
    n2 = n/2
    if 2*n2 == n:
        xmed = 0.5*(x[n2-1]+x[n2])
    else:
        xmed = x[n2]
    return xmed    

def quartiles(x, sorted=0):
    if not sorted: x = N.sort(x)
    n = len(x)
    n2 = n/2
    n4 = n/4
    if 2*n2 == n:
        xlq = 0.75*x[n4] + 0.25*x[n4+1]
        xhq = 0.25*x[3*n4] + 0.75*x[3*n4+1]
    else:
        xlq = x[n4]
        xhq = x[3*n4]
    return (xlq, xhq)


def biwt(xdata):
    # The subroutine XBIWT provides an estimator of the location and
    #    scale of the data set XDATA.  The scale uses the Biweight function
    #    in the general formula of "A-estimators." This formula is given
    #    on page of 416 in UREDA (formula 4). The BIWEIGHT scale estimate
    #    is returned as the value XSBIWT. The BIWEIGHT function is given
    #    by:
    #
    #                                  u((1-u*u)**2)     abs(u) <= 1
    #                         f(u) =
    #                                  0                 abs(u) >  1
    #
    #    where u is defined by
    #
    #                         u = (XDATA[i] - M) / c*MAD  .
    #
    #    M, MAD, and c are the median, the median absolute deviation from
    #    the median, and the tuning constant respectively. The tuning
    #    constant is a parameter which is chosen depending on the sample
    #    size and the specific function being used for the scale estimate.
    #    (See page 417 in UREDA).  Here we take c = 9.0.
    #
    # The biweght location is found using the formula:
    #
    #                         T = M + (sums)
    #
    #                         where M is the sample median and sums are
    #                         as given on page 421 in UREDA
    #
    #                         the tuning constant c is set to 6.0 for calculation
    #                         of the location as reccommended by Tukey ()
    #
    # NOTE that the biweight is meant to be an iterated estimator, but one
    #    commonly only takes the first step of the iteration.  Here we report
    #    both the one-step estimators (XLBIWT1, XSBIWT1) and the preferred
    #    fully iterated versions (XLBIWT, XSBIWT).
    n = len(xdata)
    small = 0.0001
    #    sort the data and find the median
    xm = median(xdata)
    #    call xmad to find the median absolute deviation
    xmadm = mad(xdata, xm)
    #    must choose value of the tuning constant "c"
    #       here c = 6.0 for the location estimator and
    #       9.0 for the scale estimator
    c1 = 6.0
    c2 = 9.0
    if (xmadm <= small):
        xlb = xm
        xsb = xmadm
    else:
        xlb = xm
        xsb = xmadm
        u1 = N.zeros(xdata.shape, N.float32)
        u2 = N.zeros(xdata.shape, N.float32)
        xlb_change = xsb_change = 1.0
        while xlb_change > small and xsb_change > small:
            for i in range(n):
                u1[i] = (xdata[i] - xlb)/(c1*xmadm)
                u2[i] = (xdata[i] - xlb)/(c2*xmadm)
            s1 = s2 = s3 = s4 = 0.0
            for i in range(n):
                if (abs(u2[i]) < 1.0):
                    s1 = s1 + (xdata[i] - xlb)**2 * (1.0 - u2[i]**2)**4
                    s2 = s2 + (1.0 - u2[i]**2) * (1.0 - 5.0*u2[i]**2)
                if (abs(u1[i]) < 1.0):
                    s3 = s3 + (xdata[i] - xlb) * (1.0 - u1[i]**2)**2
                    s4 = s4 + (1.0 - u1[i]**2)**2
            xlb_new = xlb + s3/s4
            xsb_new = (float(n)/float(n-1)**0.5) * (s1**0.5 / abs(s2))
            xlb_change = abs((xlb_new - xlb) / xlb_new)
            xsb_change = abs((xsb_new - xsb) / xsb_new)
            xlb = xlb_new
            xsb = xsb_new

    # load t-table created by make_t_table.py
    t_table_file = file('/tmp/t_table.pickle')
    t_table = pickle.load(t_table_file)
    t_table_file.close()

    # calculate confidence intervals
    ndof = n-1
    idof = int(0.7*ndof)
    tsigma = [1.00, 1.64, 1.96, 2.58]
    tprob = [0.68, 0.90, 0.95, 0.99]
    xtsc = []
    for tp in tprob:
        tpt = 0.5 + tp/2.0
        tvalue = t_table.t(idof,tpt)
        xtsc.append(tvalue*xsb/n**0.5)
    ci = N.array([tsigma, tprob, xtsc])
    
    return xlb, xsb, ci

 
def mad(xdata, xmed):
    # The XMAD subroutine calculates the Median Absolute Deviation from
    #    the sample median. The median, M , is subtracted from each
    #    ORDERED statistic and then the absolute value is taken. This new
    #    set of of statistics is then resorted so that they are ORDERED
    #    statistics. The MAD is then defined to be the median of this
    #    new set of statistics and is returned as XMADM. The MAD can
    #    be defined:
    #
    #                   XMADM = median{ abs(x[i] - M) }
    #
    #    where the x[i] are the values passed in the array XDATA, and
    #    the median, M, is passed in the array XLETTER. The set of stats
    #    in the brackets is assumed to be resorted. For more information
    #    see page 408 in UREDA.
    n = len(xdata)
    dhalf, n1, n2 = (0.5, 1, 2)
    xdata2 = N.absolute(xdata - xmed)
    xdata2 = N.sort(xdata2,0)
    if (float(n)/float(n2) - int(n/n2) == 0):
        i1 = n/n2
        i2 = n/n2 - n1
        xmadm = dhalf*(xdata2[i1] + xdata2[i2])
    else:
        i1 = n/n2
        xmadm = xdata2[i1]
    return xmadm

def biwt_ttest(n1, ave1, var1, n2, ave2, var2):
    df = 0.7 * (n1+n2-2)
    svar = ((n1-1)*var1+(n2-1)*var2)/df
    t = (ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2))
    prob = nr.betai(0.5*df,0.5,df/(df+t**2))
    return t, prob


def biwt_tutest(n1, ave1, var1, n2, ave2, var2):
    t = (ave1-ave2)/sqrt(var1/n1+var2/n2)
    df = (var1/n1+var2/n2)**2/((var1/n1)**2/(n1-1)+(var2/n2)**2/(n2-1))
    df = 0.7 * df
    prob = nr.betai(0.5*df,0.5,df/(df+t**2))
    return t, prob


def med_func(x, y, sig, b):
    small = 1.0e-8
    aa = median(y - b*x)
    d = (y - aa - b*x)
    mad = median(N.absolute(d))
    s = mad  / 0.6745
    d /= sig
    sign = N.compress(N.absolute(d) > small, d)
    sign = sign / N.absolute(sign)
    x = N.compress(N.absolute(d) > small, x)
    sum = N.sum(sign * x)
    return sum, s, aa


def biwt_func(x, y, sig, b): # Problems?!?
    aa = median(y - b*x)
    d = (y - aa - b*x)
    mad = median(N.absolute(d))
    s = mad  / 0.6745
    d /= sig
    # biweight
    c = 6.0
    f = d*(1-d**2/c**2)**2
    sum = N.sum(N.compress(N.absolute(d) <= c, x*f))
    # lorentzian
    #f = d/(1+0.5*d**2)
    #sum = N.sum(x*f)
    # MAD
    #small = 1.0e-8
    #sign = N.compress(N.absolute(d) > small, d)
    #sign = sign / N.absolute(sign)
    #sum = N.sum(N.compress(N.absolute(d) > small, x)*sign)
    return sum, s, aa


def fit(x, y, sig=1.0, func=med_func):
    aa, bb, siga, sigb, chi2 = nr.fit(x, y)
    b1 = bb
    f1, sigma, aa = func(x, y, sig, b1)
    if abs(f1) > 1.0e-8:
	b2 = bb + 3.0*sigb * f1/abs(f1)
    else:
	b2 = bb + 3.0*sigb
    f2, sigma, aa = func(x, y, sig, b2)
    #print 'Initial:', aa, b1, b2, f1, f2
    while f1*f2 > 0.0:
        bb = 2.0*b2-b1
        b1 = b2
        f1 = f2
        b2 = bb
        f2, sigma, aa = func(x, y, sig, b2)
        #print 'Bracketing:', aa, b1, b2, f1, f2
    sigb = 0.01*sigb
    while abs(b2-b1) > sigb:
        bb = 0.5 * (b1+b2)
        if (bb == b1 or bb == b2): break
        f, sigma, aa = func(x, y, sig, bb)
        if f*f1 >= 0.0:
            f1 = f
            b1 = bb
        else:
            f2 = f
            b2 = bb
        #print 'Bisecting:', aa, b1, b2, f1, f2
    return aa, bb, sigma


def fit_intercept(x, y, b):
    sum, s, a = med_func(x, y, 1.0, b)
    return a, s


def ks_test(x, y):
    x = N.sort(x)
    y = N.sort(y)
    nx = len(x)
    ny = len(y)
    jx = jy = 0
    fnx = fny = 0.0
    fnx_arr = N.zeros(nx+1, N.float32)
    fny_arr = N.zeros(ny+1, N.float32)
    binx = N.zeros(nx+1, N.float32)
    biny = N.zeros(ny+1, N.float32)
    d = 0.0
    while (jx < nx or jy < ny):
        if jx < nx:
            dx = x[jx]
        else:
            dx = x[-1]
        if jy < ny:
            dy = y[jy]
        else:
            dy = y[-1]
        if (dx <= dy or jy >= ny-1) and jx < nx:
            fnx = (jx+1)/float(nx)
            binx[jx+1] = dx
            fnx_arr[jx+1] = fnx
            jx += 1
        if (dy <= dx or jx >= nx-1) and jy < ny:
            fny = (jy+1)/float(ny)
            biny[jy+1] = dy
            fny_arr[jy+1] = fny
            jy += 1
        dt = abs(fny-fnx)
        if dt > d:
            d = dt
        #print '%5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f'%(dx, dy, fnx, fny, dt, d)
    n = sqrt(nx*ny/float(nx+ny))
    prob = probks((n+0.12+0.11/n)*d)
    binx[0] = binx[1] - (binx[-1] - binx[1])
    biny[0] = biny[1] - (biny[-1] - biny[1])
    return d, prob, [binx, fnx_arr], [biny, fny_arr]


def probks(alam):
    EPS1 = 0.001
    EPS2 = 1.0e-8
    fac=2.0
    sum=0.0
    termbf=0.0
    a2 = -2.0*alam*alam
    for j in range(1, 100):
        term = fac*exp(a2*j**2)
        sum += term
        if (abs(term) <= EPS1*termbf or abs(term) <= EPS2*sum):
            return sum
        fac = -fac
        termbf=abs(term)
    print 'Warning: KS-test probability function did not converge'
    return 1.0


def kuiper_test(x, y):
    x = N.sort(x)
    y = N.sort(y)
    nx = len(x)
    ny = len(y)
    jx = jy = 0
    fnx = fny = 0.0
    fnx_arr = N.zeros(nx+1, N.float32)
    fny_arr = N.zeros(ny+1, N.float32)
    binx = N.zeros(nx+1, N.float32)
    biny = N.zeros(ny+1, N.float32)
    d1 = d2 = 0.0
    while (jx < nx or jy < ny):
        if jx < nx:
            dx = x[jx]
        else:
            dx = x[-1]
        if jy < ny:
            dy = y[jy]
        else:
            dy = y[-1]
        if (dx <= dy or jy >= ny-1) and jx < nx:
            fnx = (jx+1)/float(nx)
            binx[jx+1] = dx
            fnx_arr[jx+1] = fnx
            jx += 1
        if (dy <= dx or jx >= nx-1) and jy < ny:
            fny = (jy+1)/float(ny)
            biny[jy+1] = dy
            fny_arr[jy+1] = fny
            jy += 1
        dt1 = (fny-fnx)
        dt2 = (fnx-fny)
	d1 = max(d1, dt1)
	d2 = max(d2, dt2)
    v = d1 + d2
    n = sqrt(nx*ny/float(nx+ny))
    prob = probkuiper((n+0.155+0.24/n)*v)
    binx[0] = binx[1] - (binx[-1] - binx[1])
    biny[0] = biny[1] - (biny[-1] - biny[1])
    return v, prob, [binx, fnx_arr], [biny, fny_arr]


def probkuiper(alam):
    if alam < 0.4:
	return 1.0
    else:
	EPS1 = 0.001
	EPS2 = 1.0e-8
	sum=0.0
	termbf=0.0
	alam2 = alam**2
	for j in range(1, 100):
	    j2 = j**2
	    term = 2.0 * (4.0*j2*alam2 - 1.0) * exp(-2.0*j2*alam2)
	    sum += term
	    if (abs(term) <= EPS1*termbf or abs(term) <= EPS2*sum):
		return sum
	    termbf=abs(term)
	print 'Warning: Kuiper test probability function did not converge'
	return 1.0
