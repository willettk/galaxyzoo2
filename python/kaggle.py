import numpy as np
import math
import imp
import copy
from numpy import random
from pprint import pprint
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

from os import sys
sys.path.append('/Users/willettk/Astronomy/Research/GalaxyZoo/python')
from galaxyzoo2 import get_task_dict

kaggle_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/kaggle'

tasklist = ['smooth','edgeon','bar','spiral','bulge',\
    'odd','rounded','odd_feature','bulge_shape','arms_winding','arms_number']

def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def sample_color_morph():

    # Create the test vs. training datasets for the Kaggle GZ2 competition.
    # Now base the selection on a balance of color, f_spiral, and f_merger

    with fits.open('%s/gz2_morph_colors.fits' % kaggle_path) as f:
        data = f[1].data
    
    spiral = data['t01_smooth_or_features_a02_features_or_disk_debiased']
    merger = data['t08_odd_feature_a24_merger_debiased']
    odd = data['t06_odd_a14_yes_debiased']
    colour = data['ur_color']

    # Mask the central spiral bin?

    most_spirals = (colour > 1.5) & (colour < 2.0) & (spiral > 0.8)
    all_ind = np.arange(len(most_spirals))
    spiral_ind = all_ind[most_spirals]
    not_spiral_ind = all_ind[~most_spirals]
    random.shuffle(spiral_ind)
    combined_ind = np.union1d(not_spiral_ind,spiral_ind[::2])
    some_spirals = np.zeros(len(data)).astype(bool)
    for s in combined_ind:
        some_spirals[s] = True

    spiral_somespirals = data[some_spirals]['t01_smooth_or_features_a02_features_or_disk_debiased']
    merger_somespirals = data[some_spirals]['t08_odd_feature_a24_merger_debiased']
    colour_somespirals = data[some_spirals]['ur_color']

# Try this: for each f_spiral bin as fn. of color, fit Gaussian to distribution. Keep everything redder than 1-2 sigma; randomly select from sample for -1 to +1 sigma, though. Do the same with ellipticals and the blue edge. Play with the sampling ratio as appropriate.

    xarr = np.linspace(0,4,100)
    total_left = len(spiral)
    decrease_factor = 2.0

    print ''
    print 'Starting with %6i galaxies' % total_left
    print ''

    fig = plt.figure(4)
    fig.clf()
    thresholds = (0.0,0.1,0.2,0.3,0.6,0.7,0.8,0.9)
    types = ['elliptical']*4 + ['spiral']*4 

    removed_indices = np.array((),dtype=int)
    allind = np.arange(len(spiral))

    for idx,(s,t) in enumerate(zip(thresholds,types)):
        spiralbin = (spiral > s) & (spiral <= s + 0.1) & (colour > 0.) & (colour < 4.)
        h, bin_edges = np.histogram(colour[spiralbin],bins=20)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
        coeff, var_matrix = curve_fit(gauss, bin_centres, h, p0=[5000.,1.5,1.])
        amp,mu,sigma = coeff

        if t is 'spiral':
            linepos,textpos = mu+sigma,mu-sigma
            halign = 'right'
            arrowdir = -sigma
            histcolor,gausscolor = 'blue','red'
            colorcut = colour < linepos
        else:
            linepos,textpos = mu-sigma,mu+sigma
            halign = 'left'
            arrowdir = sigma
            histcolor,gausscolor = 'red','blue'
            colorcut = colour > linepos

        ax = fig.add_subplot(4,2,idx+1)
        ax.hist(colour[spiralbin],bins=20,range=(0,4),color=histcolor)
        ax.plot(xarr,gauss(xarr,*coeff),color=gausscolor,linewidth=2)

        yl = ax.get_ylim()
        ax.set_ylim(bottom=0)
        ax.set_xlim((0,4))

        ax.vlines(linepos,yl[0],yl[1],linewidth=2,color='green')
        ax.arrow(linepos, 0.9*yl[1], arrowdir, 0., fc="g", ec="g", head_width=100., head_length=0.2 )
        ax.text(textpos, 0.9*yl[1], 'Decrease',color='g',ha=halign)

        ax.set_xlabel(r'$(g-r)$')
        ax.set_ylabel('Count')
        ax.set_title(r'%3.1f $< f_{spiral} < $ %3.1f' % (s,s+0.1))

        # Remove middle +- 1 sigma of spirals by some factor

        cut = spiralbin & colorcut
        things_removed = np.sum(cut) * (1. - 1./decrease_factor)
        total_left -= int(things_removed)
        print 'Keeping %2i percent from the color/morphology cut removes %5i/%5i galaxies (%6i total left) with %3.1f < f_spiral < %3.1f' % \
            ((1./decrease_factor)*100,things_removed,np.sum(cut),total_left,s,s+0.1)

        # plot new bins

        potentially_removed_ind = allind[cut]
        rp = random.permutation(potentially_removed_ind)
        new_ind = rp[things_removed:]
        removed = rp[:things_removed]
        old_ind = allind[spiralbin & ~colorcut]

        removed_indices = np.append(removed_indices,removed)

        ax.hist(colour[np.union1d(old_ind,new_ind)],bins=20,range=(0,4),color='k',histtype='step',linewidth=3)

        #fig.tight_layout()

    fig.savefig('%s/hist_spiral_elliptical.png' % kaggle_path)

    # Plot the remaining intermediate galaxies

    fig = plt.figure(5)
    fig.clf()
    thresholds = (0.4,0.5)
    df_middle = 2.

    print ''

    for idx,s in enumerate(thresholds):
        spiralbin = (spiral > s) & (spiral <= s + 0.1) & (colour > 0.) & (colour < 4.)
        h, bin_edges = np.histogram(colour[spiralbin],bins=20)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
        coeff, var_matrix = curve_fit(gauss, bin_centres, h, p0=[5000.,1.5,1.])
        amp,mu,sigma = coeff

        linepos,textpos = mu-sigma,mu+sigma
        halign = 'left'
        arrowdir = +sigma
        histcolor,gausscolor = 'blue','red'
        colorcut = colour > linepos

        ax = fig.add_subplot(1,2,idx+1)
        ax.hist(colour[spiralbin],bins=20,range=(0,4),color=histcolor)
        ax.plot(xarr,gauss(xarr,*coeff),color=gausscolor,linewidth=2)

        yl = ax.get_ylim()
        ax.set_ylim(bottom=0)

        ax.vlines(linepos,yl[0],yl[1],linewidth=2,color='green')
        ax.arrow(linepos, 0.9*yl[1], arrowdir, 0., fc="g", ec="g", head_width=100., head_length=0.2 )
        ax.text(textpos, 0.9*yl[1], 'Decrease',color='g',ha=halign)

        ax.set_xlabel(r'$(g-r)$')
        ax.set_ylabel('Count')
        ax.set_title(r'%3.1f $< f_{spiral} < $ %3.1f' % (s,s+0.1))

        # Remove middle +- 1 sigma of spirals by some factor

        cut = spiralbin & colorcut
        things_removed = np.sum(cut) * (1. - 1./df_middle)
        total_left -= int(things_removed)
        print 'Keeping %2i percent from the color/morphology cut removes %5i/%5i galaxies (%6i total left) with %3.1f < f_spiral < %3.1f' % \
            ((1./df_middle)*100,things_removed,np.sum(cut),total_left,s,s+0.1)

        # plot new bins

        potentially_removed_ind = allind[cut]
        rp = random.permutation(potentially_removed_ind)
        new_ind = rp[things_removed:]
        removed = rp[:things_removed]
        old_ind = allind[spiralbin & ~colorcut]

        removed_indices = np.append(removed_indices,removed)

        ax.hist(colour[np.union1d(old_ind,new_ind)],bins=20,range=(0,4),color='k',histtype='step',linewidth=3)

        #fig.tight_layout()

    fig.savefig('%s/hist_intermediate.png' % kaggle_path)

    removed_indices.sort()


# Number of intermediate galaxies is roughly the same; don't remove them
# Put mergers back in 

    fig = plt.figure(6)
    fig.clf()

    t_odd = 0.430
    t_merger = 0.5

    removed_mergers = [merger[r] if odd[r] > t_odd else 0. for r in removed_indices]
    n_removed_mergers = len(removed_mergers)

    ax = fig.add_subplot(111)
    ax.hist(removed_mergers,bins=20,range=(0,1),color='g')

    yl = ax.get_ylim()
    ax.set_ylim(bottom=0)

    ax.set_xlabel(r'$f_{merger}$')
    ax.set_ylabel('Count')
    ax.set_title('Galaxies removed from color/spiral cuts')

    fig.savefig('%s/hist_mergers.png' % kaggle_path)

    #fig.tight_layout()

    confirmed_merger = (odd > t_odd) & (merger > t_merger)
    removed_mergers = np.intersect1d(allind[confirmed_merger],removed_indices).astype(int)
    print ''
    print '%i galaxies removed from sample with f_odd > %5.3f and f_merger > %5.3f' % (len(removed_mergers),t_odd,t_merger)

    # Put the removed mergers back in

    still_good = np.union1d(np.setdiff1d(allind,removed_indices).astype(int),removed_mergers).astype(int)
    
    culled_data = data[still_good]

    print ''
    print 'Culled data length is %i' % len(culled_data)


    # plot

    H1,edges1 = np.histogramdd(np.array((colour_somespirals,spiral_somespirals)).T,bins=(20,20),range=((0,5),(0,1)))
    H2,edges2 = np.histogramdd(np.array((colour_somespirals,merger_somespirals)).T,bins=(20,20),range=((0,5),(0,1)))
    H3,edges3 = np.histogramdd(np.array((spiral_somespirals,merger_somespirals)).T,bins=(20,20),range=((0,1),(0,1)))

    fig = plt.figure(2)
    fig.clf()
    ax1 = fig.add_subplot(3,1,1)
    ax2 = fig.add_subplot(3,1,2)
    ax3 = fig.add_subplot(3,1,3)
    
    img1 = ax1.imshow(H1.T,origin='lower',interpolation='nearest',extent=[edges1[0][0], edges1[0][-1], edges1[1][0], edges1[1][-1]],aspect='auto')
    ax1.set_xlabel(r'$(u-r)$',fontsize=20)
    ax1.set_ylabel(r'$f_{spiral}$',fontsize=20)
    ax1.set_title('Observed GZ2 distribution',fontsize=25)
    cb1 = plt.colorbar(img1,ax=ax1)
    cb1.set_label(r'$N_{galaxies}$',fontsize=20)

    img2 = ax2.imshow(np.log(H2.T+1),origin='lower',interpolation='nearest',extent=[edges2[0][0], edges2[0][-1], edges2[1][0], edges2[1][-1]],aspect='auto')
    ax2.set_xlabel(r'$(u-r)$',fontsize=20)
    ax2.set_ylabel(r'$f_{merger}$',fontsize=20)
    cb2 = plt.colorbar(img2,ax=ax2)
    cb1.set_label(r'log $N_{galaxies}$',fontsize=20)
    cb2.set_label(r'log $N_{galaxies}$',fontsize=20)

    img3 = ax3.imshow(np.log(H3.T+1),origin='lower',interpolation='nearest',extent=[edges3[0][0], edges3[0][-1], edges3[1][0], edges3[1][-1]],aspect='auto')
    ax3.set_xlabel(r'$f_{spiral}$',fontsize=20)
    ax3.set_ylabel(r'$f_{merger}$',fontsize=20)
    cb3 = plt.colorbar(img3,ax=ax3)
    cb1.set_label(r'log $N_{galaxies}$',fontsize=20)
    cb3.set_label(r'log $N_{galaxies}$',fontsize=20)

    plt.show()

    fig.savefig('%s/color_spiral_merger.png' % kaggle_path)

    # Write the data to a FITS file?

    return culled_data

def split_data(data):

    # Randomly select from split data

    # Goal size for test sample
    n_test = 80000
    ngal = len(data)

    # Test sample should be 75% private, 25% public
    f_test = n_test / float(ngal)
    f_training = 1. - f_test
    f_test_private = 0.75
    f_test_public = 0.25

    private = []
    public = []
    training = []

    random_order = random.permutation(np.arange(len(data)))
    data_randomized = data[random_order]

    data_test_private = data_randomized[:int(n_test * 0.75)]
    data_test_public = data_randomized[int(n_test * 0.75):n_test]
    data_training = data_randomized[n_test:]
    
    # Write the data to a csv file

    ascii.write(data_test_private,'%s/kaggle_gz2_test_private.csv' % kaggle_path,delimiter=',')
    ascii.write(data_test_public,'%s/kaggle_gz2_test_public.csv' % kaggle_path,delimiter=',')
    ascii.write(data_training,'%s/kaggle_gz2_training.csv' % kaggle_path,delimiter=',')

    print 'Test set (private) has %5i galaxies' % len(data_test_private)
    print 'Test set (public)  has %5i galaxies' % len(data_test_public)
    print 'Training set       has %5i galaxies' % len(data_training)

    return None

def normalize_debiased_vote_fractions(dataset = 'test_public'):

    '''
    Datasets = 'test_public','test_private', and 'training'
    '''

    '''
    with fits.open('%s/kaggle_gz2_%s.fits' % (kaggle_path,dataset)) as p:
        d = p[1].data
    '''
    d = ascii.read('%s/kaggle_gz2_%s.csv' % (kaggle_path,dataset))

    # Read in each line from the kaggle_gz2_test file

    task01 = []
    task02 = []
    task03 = []
    task04 = []
    task05 = []
    task06 = []
    task07 = []
    task08 = []
    task09 = []
    task10 = []
    task11 = []

    taskarrays = (task01, task02, task03, task04, task05, task06, task07, task08, task09, task10, task11)

    for line in d:

    # Find the matching vote fractions for each task in GZ2

    # If N_votes < 0. and sum of debiased votes != 1., then normalize to 1

        for task,taskarr in zip(tasklist,taskarrays):
            td = get_task_dict(task)
            task_names_weights = [wf[:-17]+'weight' for wf in td['task_names_wf']]
            task_names_debiased = [wf[:-17]+'debiased' for wf in td['task_names_wf']]
            weights = [line[weight_name] for weight_name in task_names_weights]
            weight_sum = np.sum(weights)
            debiased = [line[debiased_name] for debiased_name in task_names_debiased]
            debiased_sum = np.sum(debiased)

            # Normalize each task sum to 1.

            if (weight_sum > 0.) and (debiased_sum != 1.) and (debiased_sum > 0.):
                debiased_normalized = [debiased_old/debiased_sum for debiased_old in debiased]
                taskarr.extend(debiased_normalized)
            else:
                taskarr.extend(debiased)

    # Convert the task floats to strings (necessary?)

    task01_str = ['%9.6f' % t for t in task01]
    task02_str = ['%9.6f' % t for t in task02]
    task03_str = ['%9.6f' % t for t in task03]
    task04_str = ['%9.6f' % t for t in task04]
    task05_str = ['%9.6f' % t for t in task05]
    task06_str = ['%9.6f' % t for t in task06]
    task07_str = ['%9.6f' % t for t in task07]
    task08_str = ['%9.6f' % t for t in task08]
    task09_str = ['%9.6f' % t for t in task09]
    task10_str = ['%9.6f' % t for t in task10]
    task11_str = ['%9.6f' % t for t in task11]

    #weights = np.zeros(len(d)).astype(int)

    # Rewrite the data solutions file as a csv

    newdata = Table({'dr7objid': d['dr7objid'],
                   'Class1.1': task01_str[0::3],
                   'Class1.2': task01_str[1::3],
                   'Class1.3': task01_str[2::3],
                   'Class2.1': task02_str[0::2],
                   'Class2.2': task02_str[1::2],
                   'Class3.1': task03_str[0::2],
                   'Class3.2': task03_str[1::2],
                   'Class4.1': task04_str[0::2],
                   'Class4.2': task04_str[1::2],
                   'Class5.1': task05_str[0::4],
                   'Class5.2': task05_str[1::4],
                   'Class5.3': task05_str[2::4],
                   'Class5.4': task05_str[3::4],
                   'Class6.1': task06_str[0::2],
                   'Class6.2': task06_str[1::2],
                   'Class7.1': task07_str[0::3],
                   'Class7.2': task07_str[1::3],
                   'Class7.3': task07_str[2::3],
                   'Class8.1': task08_str[0::7],
                   'Class8.2': task08_str[1::7],
                   'Class8.3': task08_str[2::7],
                   'Class8.4': task08_str[3::7],
                   'Class8.5': task08_str[4::7],
                   'Class8.6': task08_str[5::7],
                   'Class8.7': task08_str[6::7],
                   'Class9.1': task09_str[0::3],
                   'Class9.2': task09_str[1::3],
                   'Class9.3': task09_str[2::3],
                  'Class10.1': task10_str[0::3],
                  'Class10.2': task10_str[1::3],
                  'Class10.3': task10_str[2::3],
                  'Class11.1': task11_str[0::6],
                  'Class11.2': task11_str[1::6],
                  'Class11.3': task11_str[2::6],
                  'Class11.4': task11_str[3::6],
                  'Class11.5': task11_str[4::6],
                  'Class11.6': task11_str[5::6], \
                 #'Weight': weights \
                  'Usage': [dataset] * len(d)
                  },
                  names =['dr7objid',
                          'Class1.1',
                          'Class1.2',
                          'Class1.3',
                          'Class2.1',
                          'Class2.2',
                          'Class3.1',
                          'Class3.2',
                          'Class4.1',
                          'Class4.2',
                          'Class5.1',
                          'Class5.2',
                          'Class5.3',
                          'Class5.4',
                          'Class6.1',
                          'Class6.2',
                          'Class7.1',
                          'Class7.2',
                          'Class7.3',
                          'Class8.1',
                          'Class8.2',
                          'Class8.3',
                          'Class8.4',
                          'Class8.5',
                          'Class8.6',
                          'Class8.7',
                          'Class9.1',
                          'Class9.2',
                          'Class9.3',
                          'Class10.1',
                          'Class10.2',
                          'Class10.3',
                          'Class11.1',
                          'Class11.2',
                          'Class11.3',
                          'Class11.4',
                          'Class11.5',
                          'Class11.6', \
                          #'Weight' \
                          #'Usage'
                          ]
                           )
    ascii.write(newdata,'%s/kaggle_gz2_%s_normalized.csv' % (kaggle_path,dataset),delimiter=',')

    # What is it I need to rewrite?
        # Solutions file for the training, test public, and test private sets  (like TestSolutions.csv)
        # Make sure the dr7objid corresponds to David's GalaxyID in the solutions file (match this later, either
        #   by letting David do it himself or getting some agreement and doing it in topcat. Suppose it doesn't
        #   really matter; could renumber the galaxies at this stage?

    return None

def cumulative_probability(dataset='test_public'):

    task_class = {
        't01_smooth_or_features_a01_smooth'           :  'Class1.1',
        't01_smooth_or_features_a02_features_or_disk' :  'Class1.2',
        't01_smooth_or_features_a03_star_or_artifact' :  'Class1.3',
        't02_edgeon_a04_yes'                          :  'Class2.1',
        't02_edgeon_a05_no'                           :  'Class2.2',
        't03_bar_a06_bar'                             :  'Class3.1',
        't03_bar_a07_no_bar'                          :  'Class3.2',
        't04_spiral_a08_spiral'                       :  'Class4.1',
        't04_spiral_a09_no_spiral'                    :  'Class4.2',
        't05_bulge_prominence_a10_no_bulge'           :  'Class5.1',
        't05_bulge_prominence_a11_just_noticeable'    :  'Class5.2',
        't05_bulge_prominence_a12_obvious'            :  'Class5.3',
        't05_bulge_prominence_a13_dominant'           :  'Class5.4',
        't06_odd_a14_yes'                             :  'Class6.1',
        't06_odd_a15_no'                              :  'Class6.2',
        't07_rounded_a16_completely_round'            :  'Class7.1',
        't07_rounded_a17_in_between'                  :  'Class7.2',
        't07_rounded_a18_cigar_shaped'                :  'Class7.3',
        't08_odd_feature_a19_ring'                    :  'Class8.1',
        't08_odd_feature_a20_lens_or_arc'             :  'Class8.2',
        't08_odd_feature_a21_disturbed'               :  'Class8.3',
        't08_odd_feature_a22_irregular'               :  'Class8.4',
        't08_odd_feature_a23_other'                   :  'Class8.5',
        't08_odd_feature_a24_merger'                  :  'Class8.6',
        't08_odd_feature_a38_dust_lane'               :  'Class8.7',
        't09_bulge_shape_a25_rounded'                 :  'Class9.1',
        't09_bulge_shape_a26_boxy'                    :  'Class9.2',
        't09_bulge_shape_a27_no_bulge'                :  'Class9.3',
        't10_arms_winding_a28_tight'                  : 'Class10.1',
        't10_arms_winding_a29_medium'                 : 'Class10.2',
        't10_arms_winding_a30_loose'                  : 'Class10.3',
        't11_arms_number_a31_1'                       : 'Class11.1',
        't11_arms_number_a32_2'                       : 'Class11.2',
        't11_arms_number_a33_3'                       : 'Class11.3',
        't11_arms_number_a34_4'                       : 'Class11.4',
        't11_arms_number_a36_more_than_4'             : 'Class11.5',
        't11_arms_number_a37_cant_tell'               : 'Class11.6'
    }

    task_class_inv = {v:k for k, v in task_class.items()}

    normdata = ascii.read('%s/kaggle_gz2_%s_normalized.csv' % (kaggle_path,dataset))
    cumudata = Table(normdata,copy=True)        # Copying is necessary; otherwise it overwrites the same cell

    for n,c in zip(normdata,cumudata):
        for t in tasklist:
            td = get_task_dict(t)

            # Select tasks that have branching nodes above them in the GZ2 tree

            # This likely has an error - Task 06 really does have nodes above it, but are not
            # treated as such since they don't have Task 01 as a dependent task in galaxyzoo2.py
            # Can be fixed, but would need to rerun entire data set. - KW, 28 Dec 2013

            if 'dependent_tasks' in td.keys():
                dependent_tasks = td['dependent_tasks']
                cf = 1.

                # Compute cumulative probability of reaching that node

                if type(dependent_tasks) is tuple:
                    for deptask in dependent_tasks:
                        cf *= n[task_class[deptask.rsplit('_weighted_fraction')[0]]]
                        #print 'Corrected data for task %s, line %i; CF = %5.3f' % (t,idx,cf)
                else:
                    cf *= n[task_class[dependent_tasks.rsplit('_weighted_fraction')[0]]]
                    #print 'Corrected data for task %s, line %i; CF = %5.3f' % (t,idx,cf)

                # Now multiply the relevant normalized probabilities by the new factor

                tnum = td['var_def'][-2:]
                if cf > 0. and cf < 1.:
                    for k in task_class.keys():
                        if k[1:3] == tnum:
                            taskkey = task_class[k]
                            c[taskkey] = (n[taskkey] * cf)

    # Write the results

    ascii.write(cumudata,'%s/kaggle_gz2_%s_norm_cumu.csv' % (kaggle_path,dataset),delimiter=',')

    return None

if __name__ == '__main__':

    resampled_data = sample_color_morph()
    split_data(resampled_data)
    for dset in ('test_public','test_private','training'):
        normalize_debiased_vote_fractions(dataset=dset)
        cumulative_probability(dataset=dset)

