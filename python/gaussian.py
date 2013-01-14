# functions for plotting and fitting gaussians
from math import sqrt, pi, log10
import numpy as N

def g_residuals(p, x, y):
    ampl, centre, sigma = p
    g = gaussian(x, ampl, centre, sigma)
    return y - g

def multi_g_residuals(p, n, x, y):
    ampl = p[0::3]
    centre = p[1::3]
    sigma = p[2::3]
    g = multi_gaussian(x, ampl, centre, sigma)
    return y - g    

def g2d_residuals(p, shape, z):
    ampl, xcentre, ycentre, xsigma, ysigma = p
    g = gaussian2d(shape, ampl, xcentre, ycentre, xsigma, ysigma)
    r = N.ravel(z - g)
    return r

def gaussian(x, ampl, centre, sigma):
    # with unit area
    return (ampl/(sigma*sqrt(2.0*pi)) *
	    N.exp(-(x-centre)**2/(2.0*sigma**2)))

def gaussian_unit_max(x, ampl, centre, sigma):
    # with unit maximum
    return ampl * N.exp(-(x-centre)**2/(2.0*sigma**2))

def gaussian2dx(x, y, ampl, xcentre, ycentre, xsigma, ysigma):
    return (ampl/(xsigma*ysigma*2*pi) *
	    N.exp(-(x-xcentre)**2/(2.0*xsigma**2)) *
	    N.exp(-(y-ycentre)**2/(2.0*ysigma**2)))

def multi_gaussian_derivatives(x, ampl, centre, sigma, shift):
    dg_dampl = []
    dg_dcentre = []
    dg_dsigma = []
    for i in range(len(ampl)):
	dg_dampl.append(1.0/(sigma[i]*sqrt(2*pi)) *
			N.exp(-(x-centre[i]-shift)**2/(2.0*sigma[i]**2)))
	dg_dcentre.append(ampl[i]/(sigma[i]*sqrt(2*pi)) *
			   N.exp(-(x-centre[i]-shift)**2/(2.0*sigma[i]**2)) *
			   ((x-centre[i]-shift)/sigma[i]**2))
	dg_dsigma.append(-ampl[i]/sqrt(2*pi) *
			  N.exp(-(x-centre[i]-shift)**2/(2.0*sigma[i]**2)) *
			  (1.0/sigma[i]**2 - (x-centre[i]-shift)**2/sigma[i]**4))
    if len(dg_dcentre) > 1:
	dg_dshift = N.sum(dg_dcentre, 0)
    else:
	dg_dshift = dg_dcentre
    return dg_dampl, dg_dcentre, dg_dsigma, dg_dshift

def gaussian_derivatives(x, ampl, centre, sigma):
    dg_dampl = (1.0/(sigma*sqrt(2*pi)) *
		N.exp(-(x-centre)**2/(2.0*sigma**2)))
    dg_dcentre = -(ampl/(sigma*sqrt(2*pi)) *
		  N.exp(-(x-centre)**2/(2.0*sigma**2)) *
		  (-(x-centre)/sigma**2))
    dg_dsigma = -(ampl/sqrt(2*pi) *
		  N.exp(-(x-centre)**2/(2.0*sigma**2)) *
		  (1.0/sigma**2 - (x-centre)**2/sigma**4))
    return dg_dampl, dg_dcentre, dg_dsigma

def multi_gaussian_derivatives_numerical(x, ampl, centre, sigma):
    eps = 1.0e-8
    g = multi_gaussian(x, ampl, centre, sigma)
    dg_dampl = []
    for i in range(len(ampl)):
	a = N.array(ampl)
	a[i] += eps
	ge = multi_gaussian(x, a, centre, sigma)
	dg_dampl.append((ge - g)/eps)
    dg_dcentre = []
    for i in range(len(centre)):
	c = N.array(centre)
	c[i] += eps
	ge = multi_gaussian(x, ampl, c, sigma)
	dg_dcentre.append((ge - g)/eps)
    dg_dsigma = []
    for i in range(len(sigma)):
	s = N.array(sigma)
	s[i] += eps
	ge = multi_gaussian(x, ampl, centre, s)
	dg_dsigma.append((ge - g)/eps)
    return dg_dampl, dg_dcentre, dg_dsigma

def gaussian_derivatives_numerical(x, ampl, centre, sigma):
    eps = 1.0e-8
    g = gaussian(x, ampl, centre, sigma)
    ge = gaussian(x, ampl+eps, centre, sigma)
    dg_dampl = (ge - g)/eps
    ge = gaussian(x, ampl, centre+eps, sigma)
    dg_dcentre = (ge - g)/eps
    ge = gaussian(x, ampl, centre, sigma+eps)
    dg_dsigma = (ge - g)/eps
    return dg_dampl, dg_dcentre, dg_dsigma


def gaussian2d(shape, ampl, xcentre, ycentre, xsigma, ysigma):
    def make_dx(j, i):
	return (i-float(xcentre))
    def make_dy(j, i):
	return (j-float(ycentre))
    dx = N.fromfunction(make_dx, shape)
    dy = N.fromfunction(make_dy, shape)
    dx2 = dx**2
    dy2 = dy**2
    return (ampl/(xsigma*ysigma*2*pi) *
	    N.exp(-dx2/(2.0*xsigma**2)) *
	    N.exp(-dy2/(2.0*ysigma**2)))

def multi_gaussian(x, ampl, centre, sigma, shift=0.0):
    if type(x) != type(1.0):
	g = N.zeros(x.shape, N.float64)
    else:
	g = 0
    for i in range(len(ampl)):
	g = g + gaussian(x, ampl[i], centre[i]+shift, sigma[i])
    return g
