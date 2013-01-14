# cosmology.py
# (formerly distance_modulus.py)

from math import log, sqrt, pi, sin, cos, exp
from nr import erffc
#from ppgplot_spb import *
from checkarray import checkarray
import scipy.integrate
import numpy as N

# WMAP 1-year results
h0_WMAP1 = 0.71
H0_WMAP1 = h0_WMAP1 * 100.0
omega_m0_WMAP1 = 0.135 / h0_WMAP1**2

# WMAP 3-year results
h0_WMAP3 = 0.73
H0_WMAP3 = h0_WMAP3 * 100.0
omega_m0_WMAP3 = 0.127 / h0_WMAP3**2

# WMAP 7-year results
h0_WMAP7 = 0.704
H0_WMAP7 = h0_WMAP7 * 100.0
omega_m0_WMAP7 = 0.2726

# WMAP 9-year results
h0_WMAP9 = 0.7181
H0_WMAP9 = h0_WMAP9 * 100.0
omega_m0_WMAP9 = 0.2726

"""
H0_std = 70.0
omega_m0_std = 0.30
omega_lambda0_std = 1.0 - omega_m0_std
"""

H0_std = H0_WMAP9
omega_m0_std = omega_m0_WMAP9
omega_lambda0_std = 1.0 - omega_m0_std

H0_classical = 75.0
q0_classical = 0.05

c0 = 299792.458  # km/s

# mass of Sun in kg
M_sun_kg = 1.9889e30  # kg

# Mpc in kilometres
Mpc_km = 3.0857e19  # km

# Newton's gravitational constant
G_N = 6.673e-11  # m**3 kg**(-1) s**(-2)

# Number of square degrees over full sky
sky_sq_deg = 4*pi * (180/pi)**2

# The distance modulus for a flat universe with
# omega_matter + omega_lambda = 1
def dmod_flat(z, H0=H0_std, omega_m0=omega_m0_std):
    dL = dL_flat(z, H0, omega_m0)
    mu = 5.0 * N.log10(dL*1.E6) - 5.0
    return mu

# Luminosity distance for a flat universe with
# omega_matter + omega_lambda = 1
# in Mpc
def dL_flat(z, H0=H0_std, omega_m0=omega_m0_std):
    dL = (1.0+z) * (c0/H0) * dc(omega_m0, z)
    return dL

# Comoving distance for a flat universe with
# omega_matter + omega_lambda = 1
# in Mpc
def dC_flat(z, H0=H0_std, omega_m0=omega_m0_std):
    dC =  (c0/H0) * dc(omega_m0, z)
    return dC

# Angular diameter distance for a flat universe with
# omega_matter + omega_lambda = 1
def dA_flat(z, H0=H0_std, omega_m0=omega_m0_std):
    dA = dL_flat(z, H0, omega_m0) / (1.0+z)**2
    return dA

# Angular scale distance for a flat universe
arcsec_in_rad = 180*60*60/pi
kpc_in_Mpc = 10**3
ang_scale_conversion = kpc_in_Mpc / arcsec_in_rad

def ang_scale_flat(z, H0=H0_std, omega_m0=omega_m0_std):
    # kpc/arcsec
    dA = dA_flat(z, H0, omega_m0)
    ang_scale = dA * ang_scale_conversion
    return ang_scale

# The distance modulus for a classical universe with
# omega_lambda = 0
def dmod_classical(z, H0=H0_classical, q0=q0_classical):
    dL = dL_classical(z, H0, q0)
    mu = 5.0 * N.log10(dL*1.E6) - 5.0
    return mu

# Luminosity distance for a classical universe with
# omega_lambda = 0
# in Mpc
def dL_classical(z, H0=H0_classical, q0=q0_classical):
    dL = c0/(H0*q0*q0) * (q0*z + (q0-1.0) * (sqrt(1.0 + 2.0*q0*z)-1.))
    return dL

# Angular diameter distance for a classical universe with
# omega_lambda = 0
def dA_classical(z, H0=H0_classical, q0=q0_classical):
    dA = dL_classical(z, H0, q0) / (1.0+z)**2
    return dA

# calculate comoving distance (to a factor) for a flat lambda universe
# old integration method
def dc_old(omega_m0, z, dz=None):
    # probably not very efficient:
    # takes 0.25 cpu seconds for z=1.5
    #       0.08 cpu seconds for z=0.5
    # should use clever integration technique instead
    if dz is None: dz = 0.00001
    d = 0.0
    rom0 = 1.0/omega_m0
    rsom0 = 1.0/sqrt(omega_m0)
    for z1 in (N.arange(int(z/dz)+1) * dz):
        d = d + rsom0 * dz/sqrt((1.0+z1)**3.0 - 1.0 + rom0)
    return d

# calculate comoving distance (to a factor) for a flat lambda universe
# improved integration method
def dc_no_array(omega_m0, z):
    # function to integrate
    rom0 = 1.0/omega_m0
    rsom0 = 1.0/sqrt(omega_m0)
    def rEz(z1):
        return rsom0 / sqrt((1.0+z1)**3.0 - 1.0 + rom0)
    d, derr = scipy.integrate.quad(rEz, 0.0, z)
    return d


# calculate comoving distance (to a factor) for a flat lambda universe
# improved integration method, added support for arrays
def dc(omega_m0, z):
    tz = str(type(z))
    if 'float' in tz or 'int' in tz:
	z = [z]
    z = N.asarray(z)
    d = N.zeros(z.shape)
    for i, zi in enumerate(z):
	# function to integrate
	rom0 = 1.0/omega_m0
	rsom0 = 1.0/sqrt(omega_m0)
	def rEz(z1):
	    return rsom0 / sqrt((1.0+z1)**3.0 - 1.0 + rom0)
	di, dierr = scipy.integrate.quad(rEz, 0.0, zi, limit=100)
	d[i] = di
    if len(d) == 1:
	d = d[0]
    return d

# calculate the look-back time for a flat lambda universe
def lt_flat_old(z, H0=H0_std, omega_m0=omega_m0_std, dz=None):
    # probably not very efficient
    # should use clever integration technique instead
    # dz=0.0001 gives an accuracy of 0.01 Gyr for all z
    # dz=0.00001 gives an accuracy of 0.001 Gyr for all z
    if dz is None: dz = 0.00001
    t = 0.0
    omega_lambda = 1.0 - omega_m0
    for z1 in (N.arange(int(z/dz)+1) * dz):
        zfactor = 1.0 + z1
        t = t + dz / ( zfactor * sqrt( zfactor**2 * (1+omega_m0*z1) - z1*(2+z1)*omega_lambda ) )
    # t currently a fraction of the Hubble time
    # convert into Gyr
    mperpc = 3.085677581e16
    secperday = 31557600
    H0persec = H0 * 1.0e-3 / mperpc
    H0perGyr = H0persec * secperday * 1.0e9
    t = t / H0perGyr
    return t

# calculate the look-back time for a flat lambda universe
# improved integration method
def lt_flat(z, H0=H0_std, omega_m0=omega_m0_std):
    omega_lambda = 1.0 - omega_m0
    z = N.asarray(z)
    d = N.zeros(z.shape)
    for i, zi in enumerate(z):
        # function to integrate
        def intfn(z1):
            zfactor = 1.0 + z1
            return 1.0 / (zfactor * sqrt(zfactor**2 * (1+omega_m0*z1) - 
                                         z1*(2+z1)*omega_lambda))
        t, terr = scipy.integrate.quad(intfn, 0.0, zi)
        # t currently a fraction of the Hubble time
        # convert into Gyr
        mperpc = 3.085677581e16
        secperday = 31557600
        H0persec = H0 * 1.0e-3 / mperpc
        H0perGyr = H0persec * secperday * 1.0e9
        d[i] = t / H0perGyr
    if len(d) == 1:
	d = d[0]
    return d

    return t

# calculate age of a flat universe
def age_flat(z, H0=H0_std, omega_m0=omega_m0_std):
    mperpc = 3.085677581e16
    secperday = 31557600
    H0persec = H0 * 1.0e-3 / mperpc
    H0perGyr = H0persec * secperday * 1.0e9
    soom0 = sqrt(1.0-omega_m0)
    age_now = 2.0/(3.0*H0perGyr) * 1.0/(2*soom0) * log((1+soom0)/(1-soom0))
    lbt = lt_flat(z, H0, omega_m0)
    age = age_now - lbt
    return age

# calculate ratio of comoving volume elements at two redshifts
# deduced from MIT astronomy course notes by Edmund Bertschinger, 1999
#http://ocw.mit.edu/NR/rdonlyres/Physics/8-942Fall2001/2F658E61-68A8-40F4-9168-B7AD0E23CA49/0/cosmog.pdf
def vol_ratio(z1, z2, H0=H0_std,
              omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    # ratio of comoving volume element at z1 to that at z2
    Hz1 = H(z1, H0, omega_m0, omega_lambda0)
    Hz2 = H(z2, H0, omega_m0, omega_lambda0)
    dc1 = dc(omega_m0, z1)
    dc2 = dc(omega_m0, z2)
    return Hz2/Hz1 * (dc1/dc2)**2


# calculate the comoving volume enclosed between two redshifts
# deduced from eqn. 9 in Shen et al, 2003, MNRAS, 343, 978
# All-sky, in cubic Mpc
def vol(zmin, zmax, H0=H0_std,
        omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    zmax = checkarray(zmax)
    zmin = checkarray(zmin)
    v = N.zeros(zmax.shape)
    # function to integrate
    def intfn(z):
        return (4*pi * dA_flat(z, H0, omega_m0)**2 * (1+z)**2 * c0 / 
                H(z, H0, omega_m0, omega_lambda0))
    for i in range(len(zmax)):
        zmini = zmin[i]
        zmaxi = zmax[i]
        vi, vierr = scipy.integrate.quadrature(intfn, zmini, zmaxi, tol=1.0e-3, maxiter=100)
        # this tol is easily sufficient for any reasonable zmin, zmax
        v[i] = vi
    if len(v) == 1:
	v = v[0]
    return v

# calculate the maximum redshift at which an object with absolute
# magnitude mabs could be included in a survey limited to apparent
# magnitude mapplim.  Optimised to do many objects in one go.
def zmax(mabs, mapplim, zlow, zhigh, nz=None):
    mabs = checkarray(mabs)
    if nz is None:
        deltaz = 0.001
        nz = int((zhigh-zlow)/deltaz)
    deltaz = (zhigh-zlow)/float(nz)
    zlist = N.arange(zhigh+deltaz/10, zlow, -deltaz)
    dmodlist = dmod_flat(zlist)
    mabslist = mapplim - dmodlist
    i = N.searchsorted(mabslist, mabs, side='left')
    ihigh = i == nz+1
    N.putmask(i, ihigh, nz)
    ilow = i == 0
    N.putmask(i, ilow, 1)
    z1 = zlist[i-1]
    z2 = zlist[i]
    m1 = mabslist[i-1]
    m2 = mabslist[i]
    s = (z2-z1)/(m2-m1)
    zmax = z1 + s*(mabs-m1)
    N.putmask(zmax, ilow, zhigh)
    N.putmask(zmax, ihigh, zlow)
    if zmax.shape == (1,):
       zmax = zmax[0]
    return zmax

# calculate the maximum redshift at which an object with absolute
# size rkpc could be included in a survey limited to apparent
# size raslim.  Optimised to do many objects in one go.
def zmax_size(rkpc, raslim, zlow, zhigh, nz=None):
    rkpc = checkarray(rkpc)
    if nz is None:
        deltaz = 0.001
        nz = int((zhigh-zlow)/deltaz)
    deltaz = (zhigh-zlow)/float(nz)
    zlist = N.arange(zlow, zhigh+deltaz/10, deltaz)
    angscalelist = ang_scale_flat(zlist)
    rkpclist = raslim * angscalelist
    i = N.searchsorted(rkpclist, rkpc)
    ihigh = i == nz+1
    N.putmask(i, ihigh, nz)
    ilow = i == 0
    N.putmask(i, ilow, 1)
    z1 = zlist[i-1]
    z2 = zlist[i]
    r1 = rkpclist[i-1]
    r2 = rkpclist[i]
    s = (z2-z1)/(r2-r1)
    zmax = z1 + s*(rkpc-r1)
    N.putmask(zmax, ilow, zlow)
    N.putmask(zmax, ihigh, zhigh)
    if zmax.shape == (1,):
       zmax = zmax[0]
    return zmax

# calculate the maximum redshift at which an object with absolute
# surface brightness sbabs (could be included in a survey limited to
# apparent surface brightness sbapplim.
# Optimised to do many objects in one go.
def zmax_sb(sbabs, sbapplim, zlow, zhigh, nz=None):
    sbabs = checkarray(sbabs)
    if nz is None:
        deltaz = 0.001
        nz = int((zhigh-zlow)/deltaz)
    deltaz = (zhigh-zlow)/float(nz)
    zlist = N.arange(zhigh+deltaz/10, zlow, -deltaz)
    dmodlist = dmod_flat(zlist)
    angscalelist = ang_scale_flat(zlist)
    sbabslist = sbapplim - dmodlist + 2.5*N.log10(angscalelist**2)
    i = N.searchsorted(sbabslist, sbabs)
    ihigh = i == nz+1
    N.putmask(i, ihigh, nz)
    ilow = i == 0
    N.putmask(i, ilow, 1)
    z1 = zlist[i-1]
    z2 = zlist[i]
    sb1 = sbabslist[i-1]
    sb2 = sbabslist[i]
    s = (z2-z1)/(sb2-sb1)
    zmax = z1 + s*(sbabs-sb1)
    N.putmask(zmax, ilow, zhigh)
    N.putmask(zmax, ihigh, zlow)
    if zmax.shape == (1,):
       zmax = zmax[0]
    return zmax


#----------------------------------------------------------------------
# Calculate some cosmological structure formation information
# from the equations in Mo & White, 2002, MNRAS, 336, 112.

# calculate the number of haloes per unit comoving volume at
# redshift z with mass in the interval [M, M+dM]
# using the Sheth & Tormen (1999, 2001, 2002) function
def n(M, z, dM,
      omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    A = 0.322
    a = 0.707
    q = 0.3
    rho_m0 = omega_m0_std * rho_c(H(z))
    nu0 = nu_from_M(M, z,  omega_m0, omega_lambda0)
    nu1 = nu_from_M(M+dM, z,  omega_m0, omega_lambda0)
    nu0_primed = sqrt(a) * nu0
    print 'nu0_primed:', nu0_primed
    nu1_primed = sqrt(a) * nu1
    dnu_primed_by_dM = (nu1_primed - nu0_primed) / dM
    print 'dnu_primed_by_dM:', dnu_primed_by_dM
    term1 = 1.0 + 1.0/(nu0_primed**(2*q))
    term2 = sqrt(2.0/pi) * (rho_m0 / M) * dnu_primed_by_dM
    term3 = exp(-(nu0_primed**2) / 2.0)
    print 'term1:', term1
    print 'term2:', term2
    print 'term3:', term3
    return A * term1 * term2 * term3

# calculate the fraction of all mass in haloes with mass exceeding M
# (eqn. 15)
def F(nu):
    return 0.4 * (1.0 + 0.4/(nu**0.4)) * erffc(0.85*nu/sqrt(2))

# calculate the sigmaM corresponding to a given nu,
# the sigma of the halo in terms of a density contrast
def sigmaM_from_nu(nu, z, omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    delta_c = 1.69
    return delta_c / (nu * D(z, omega_m0, omega_lambda0))

# calculate nu, the sigma of the halo in terms of a density
# contrast, given sigmaM
def nu_from_M(M, z, H0=H0_std, omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    delta_c = 1.69
    R = R_Lagrangian(M, H0, omega_m0)
    varM = var(R, H0, omega_m0)
    sigmaM = sqrt(varM)
    Dz = D(z, omega_m0, omega_lambda0)
    print 'Dz:', Dz
    print 'sigmaM:', sigmaM
    return delta_c / (sigmaM * Dz)

# calculate the growth factor for linear perturbations
# following Carroll, Press & Turner Edwin (1992), as given by
# Mo & White (2002)
def D(z, omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    omega_mz = omega_m(z, omega_m0, omega_lambda0)
    omega_lambdaz = omega_lambda(z, omega_m0, omega_lambda0)
    return g(z, omega_mz, omega_lambdaz) / (g(0, omega_m0, omega_lambda0) * (1.0+z))

# calculate function used in D (eqn. 10)
def g(z, omega_m, omega_lambda):
    term1 = omega_m**(4.0/7.0) - omega_lambda
    term2 = (1 + omega_m/2.0)*(1+omega_lambda/70.0)
    return (5.0/2.0) * omega_m / (term1 + term2)

# calcuate omega_matter at any redshift (eqn. 2)
def omega_m(z, omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    return omega_m0 * (1.0+z)**3 / (E(z, omega_m0, omega_lambda0)**2)

# calcuate omega_lambda at any redshift (eqn. 11)
def omega_lambda(z, omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    return omega_lambda0 / (E(z, omega_m0, omega_lambda0)**2)

# calculate the factor giving the evolution of the hubble constant
# as a function of redshift (eqn. 3)
def E(z, omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    omega_total = omega_m0 + omega_lambda0
    return N.sqrt(omega_lambda0 + (1.0 - omega_total) * (1.0+z)**2 +
                  omega_m0 * (1.0+z)**3)

# calculate the Hubble constant at any redshift (eqn. 2)
def H(z, H0=H0_std, omega_m0=omega_m0_std, omega_lambda0=omega_lambda0_std):
    return H0 * E(z, omega_m0, omega_lambda0)

# calculate Lagrangian radius of a halo of mass M
def R_Lagrangian(M, H0=H0_std, omega_m0=omega_m0_std):
    return (3.0 * M / (4 * pi * omega_m0 * rho_c(H0)))**(1.0/3.0)

# calculate critical density of universe given Hubble constant
def rho_c(H):
    # H in km s**(-1) Mpc**(-1)
    H = H / Mpc_km# s**(-1)
    rho = 3.0 * H**2 / (8.0 * pi * G_N)  # kg m**(-3)
    return rho / M_sun_kg  * (Mpc_km * 1000.0)**3  # M_sun Mpc**(-3)

# transfer function representing differential growth since early times
# from Bardeen et al. (1986) as given in Mo & White (2002) eqn. 8.
def T(k, H0=H0_std, omega_m0=omega_m0_std):
    q = k / (omega_m0 * (H0/100.0)**2)
    if abs(q) < 1.0e-6:
        term1 = 1.0
    else:
        term1 = log(1 + 2.34*q) / (2.34*q)
    term2 = 1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4
    return term1 * term2**(-1.0/4.0)

# approximation of CDM power spectrum (eqn. 7)
def P(k, H0=H0_std, omega_m0=omega_m0_std):
    return k * T(k, H0, omega_m0)**2

# calculate normalisation of P
def intP(H0=H0_std, omega_m0=omega_m0_std, dk=None):
    if dk is None:
            dk = 0.0001
    #pgopen('P_func.ps/ps')
    k_array = []
    term_array = []
    sum = 0.0
    i = 0
    kmax = 0
    while 1:
        k = (i+0.0)*dk 
        term = P(k, H0, omega_m0)
        sum = sum + term * dk
        k_array.append(k)
        term_array.append(term)
        if i > 0:
            if kmax == 0:
                if (term / sum) < 0.01:
                    kmax = i*2
            else:
                if i >= kmax:
                    #print 'integration took %i steps'%i
                    break
            #print k, sum, term, term*dk / sum
        else:
            pass
            #print k, sum, term
        i = i+1
    #k_array = N.array(k_array)
    #term_array = N.array(term_array)
    #pgenv(0.0, 1.0, 0.0, 12.0)
    #pgline(k_array, term_array/sum)
    #pgclos()
    return sum
    

# variance of Gaussian density field at a given radius
def var(R, H0=H0_std, omega_m0=omega_m0_std, dk=None):
    # this is integrated using a very simple method
    # the function consists of a number of similarly spaced,
    # decreasing peaks, with only the first three contributing
    # significantly to the integral (first is much bigger)
    if dk is None:
            dk = 0.01 / R  # seems to give about 1000 steps and
                           # decent accuracy (<1%) for 0.03 < R < 100
    #f = file('var_func', 'w')
    #pgopen('var_func.ps/ps')
    #k_array = []
    #term_array = []
    sum = 0.0
    i = 0
    kmax = None
    while 1:
        k = i*dk 
        term = k**2 * P(k, H0, omega_m0) * W(k*R)**2
        sum = sum + term * dk
        #f.write('%16.6g%16.6g%16.6g\n'%(k, term, sum))
        #k_array.append(k)
        #term_array.append(term)
        if i > 0:
            if kmax is None:
                if (term / sum) < 0.01:
                    # given good sampling provided by choice of dk
                    # above, this occurs at end of first peak, so
                    # setting kmax to three times this point gets
                    # all three significant peaks.
                    kmax = i*3
            else:
                if i >= kmax:
                    #print 'integration took %i steps'%i
                    break
            #print k, sum, term, term*dk / sum
        else:
            pass
            #print k, sum, term
        i = i+1
    #f.close()
    #k_array = N.array(k_array)
    #term_array = N.array(term_array)
    #pgenv(0.0, k, 0.0, 2.0e-5)
    #pgline(k_array, term_array)
    #pgclos()
    arb_scaling = 1.0e18
    return sum / (2.0 * pi**2) / intP() * arb_scaling


# Fourier transform of a top-hat filter with radius R
# for use by var (x=kR)
def W(x):
    if x < 1e-6:
        return 1.0
    else:
        return 3.0 * (sin(x) - x*cos(x)) / x**3
