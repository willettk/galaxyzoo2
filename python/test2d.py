import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import cPickle as pickle
import galaxyzoo2 as gz2
from scipy.optimize import fmin_powell

x = np.arange(-24.,-16.,0.1)
y = np.arange(0.,15.,0.1)
gz_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/'
pkl_path = gz_path+'pickle/'

def tilt_function(p,x,y):

    z = np.zeros((len(x),len(y)),np.float)
    
    for i in np.arange(len(x)):
        z[i] = p[0]*(x[i] - p[1]) + p[2]*(y - p[3]) + p[4]
    return z

    
def ratio_minfunc(p, x, y, z, w):
    f = tilt_function(p, x, y)
    r2 = (z - f)**2                                     # Difference between model and data
    df = (w > 0).ravel().astype(np.int).sum() - len(p)  # Degrees of freedom: N_cells - N_parameters
    s = (w * r2).sum() / df                             # Chi-squared statistic
    return s

def bootstrap_fmin(f, p0, x, y, z, mask, nboot):
    plist = np.zeros((nboot, len(p0)), np.float)                                # Zero array for number of tries, number of parameters
    ndata = len(mask.ravel().nonzero()[0])                                      # number of non-masked cells in the data
    for i in range(nboot):
        bootmask = np.zeros(mask.shape)                                         # Grid of zeros same size as mask, data
        while bootmask.sum() < ndata:                                           # Loop until grid is filled at data locations
            rx = int(np.random.uniform(0, x.shape))                             # Pick random indices from the bin arrays
            ry = int(np.random.uniform(0, y.shape))                             # 
            if mask[rx,ry] > 0.000001:                                          # If existing data mask in that cell is non-zero
                bootmask[rx,ry] += 1                                            #   then put a 1 in bootmask cell
                                                                                # Preserves total weight, but assigns random weights
                                                                                #   to each data cell. Perturbs to get multiple guesses. 
        bootmask *= mask                                                        # Multiply new weights by original mask
        plist[i] = fmin_powell(f, p0, (x, y, z, bootmask),disp=0)               # Fit parameters based on new mask passed to function
    p = fmin_powell(f, p0, (x, y, z, mask))                                     # Fit parameters based on original (unweighted) mask
    pmed = np.median(plist, axis=0)                                             # Find the median best fit parameters for each variation in the list
    perr = np.array([(gz2.scoreatpercentile(plist[:,k], 84.0) -                     # The error is the mean of the 84th percentile parameter
            gz2.scoreatpercentile(plist[:,k], 16.0))/2.0 for k in range(len(p))])   #   and the 16th percentile parameter. Ah. Looping is only to estimate error.
    return p, perr

def plot_func(nboot=2):

    pinit = np.array([0.5,-20.0,0.5,7.0,0.1])

    chi2r = 0.0
    niter = 0
    sys = 0.0

    task_dict = gz2.get_task_dict('arms_number')
    data = pickle.load(open(pkl_path+'%s_r01_local_ratio_baseline_masked.pkl' % task_dict['var_def'],'rb'))
    centers_redshift, centers_mag, centers_size,edges_redshift, edges_mag, edges_size = gz2.get_bins(task_dict)

    weights = np.ones_like(data)
    p, perr = bootstrap_fmin(ratio_minfunc, pinit,
                             centers_mag, centers_size, data, weights, nboot)
    chi2r = ratio_minfunc(p, centers_mag, centers_size, data, weights)

    print ' '
    print 'Number of bootstrap iterations: %i' % nboot
    print 'param = %s' % p
    print 'param errors = %s' % perr
    print 'param errors (percent) = %s' % (perr/p)
    print 'chi2r = %f' % chi2r
    print ' '

    ratio_baseline_fit = tilt_function(p, x, y)

    z = ratio_baseline_fit

    fig = plt.figure(88)
    fig.clf()
    ax = fig.add_subplot(111)
    im = ax.imshow(z.T, 
        extent = (min(x),max(x),min(y),max(y)),
        interpolation='nearest', 
        origin='lower')
    ax.set_aspect('auto')
    ax.set_xlabel(r'$M_R [mag]$',fontsize=22)
    ax.set_ylabel(r'$R_{50} [kpc]$',fontsize=22)

    imf = ax.imshow(z.T,
                   alpha=1.0,
                   extent = (min(x),max(x),min(y),max(y)),
                   vmin=-1.5,vmax=1.5,
                   interpolation='nearest', 
                   origin='lower')
    cb = plt.colorbar(imf)

    # Plot the local masked baseline relation as semi-transparent layer

    cmap_gray = cm.gray
    cmap_gray.set_bad(alpha=0.0)
    im = ax.imshow(data.T, 
                   alpha=0.5,
                   extent = (min(x),max(x),min(y),max(y)),
                   cmap = cmap_gray,
                   vmin=-1.5,vmax=1.5,
                   interpolation='nearest',
                   origin='lower')

    ax.autoscale(False)
    ax.set_aspect('auto')
