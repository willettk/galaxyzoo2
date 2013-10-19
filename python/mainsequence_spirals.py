from astropy.io import fits as pyfits
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

def get_mpadata():

    p = pyfits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/fits/mpajhu_gz2.fits')
    mpajhu = p[1].data
    p.close()

    obliques = (mpajhu['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.430) & \
               (mpajhu['t02_edgeon_a05_no_debiased'] >= 0.715) & \
               (mpajhu['t04_spiral_a08_spiral_debiased'] >= 0.619) & \
               (mpajhu['t04_spiral_a08_spiral_weight'] >= 20)

    '''
    Should check later as to how including low S/N galaxies affects results
    '''
    starforming = (mpajhu['bpt'] == 1)# | (mpajhu['bpt'] == 2)

    return mpajhu, obliques, starforming

def most_massive_spiral(n=1, printtable=False):

    mpajhu, obliques, starforming = get_mpadata()

    odata = mpajhu[obliques]

    avg_ind = np.argmax(odata['AVG_MASS'])
    median_ind = np.argmax(odata['MEDIAN_MASS'])

    avg_inds = odata['AVG_MASS'].argsort()[-n:][::-1] 

    '''
    if avg_ind != median_ind:
        print 'Average =  %6.1f, Median = %6.1f', (odata[avg_ind]['AVG_MASS'],odata[median_ind]['MEDIAN_MASS'])
    '''

    data = odata[avg_inds]

    from galaxyzoo2 import ra_deg_to_sex, dec_deg_to_sex

    if printtable:
        for im,mm in enumerate(data):
            print 'MS%02i  & J%s%s  & %5.2f & %5.2f & %5.2f & %5.2f & %f & %5.2f & %5.2f & %6.3f  \\\\' % (im+1,ra_deg_to_sex(mm['RA']).translate(None,':'),dec_deg_to_sex(mm['DEC']).translate(None,':'),mm['PETROMAG_R'], mm['PETROMAG_MR'],mm['PETROMAG_U'] - mm['PETROMAG_R'],mm['PETROMAG_G'] - mm['PETROMAG_R'],mm['REDSHIFT'],mm['PETROR50_R'],mm['PETROR50_R_KPC'],mm['AVG_MASS']) 

    return data

def fit_lines(savefig=False):

    mpajhu, obliques, starforming = get_mpadata()

    ms = mpajhu[obliques & starforming]

    ms_mass = ms['AVG_MASS']
    ms_sfr = ms['AVG_SFR']

    arms_tight = (ms['t10_arms_winding_a28_tight_debiased'] >= 0.5)
    arms_medium = (ms['t10_arms_winding_a29_medium_debiased'] >= 0.5)
    arms_loose = (ms['t10_arms_winding_a30_loose_debiased'] >= 0.5)

    err = np.mean((np.abs(ms['P16_SFR'] - ms['AVG_SFR']),np.abs(ms['P84_SFR'] - ms['AVG_SFR'])))
    weights = 1./err**2

    xarr = np.linspace(6,12,100)

    (a,b),res,_,_,_ = np.polyfit(ms_mass,ms_sfr,1,full=True)
    yarr = np.polyval([a,b],xarr)
    (a_tight,b_tight),res_tight,_,_,_ = np.polyfit(ms_mass[arms_tight],ms_sfr[arms_tight],1,full=True)
    yarr_tight = np.polyval([a_tight,b_tight],xarr)
    (a_medium,b_medium),res_medium,_,_,_ = np.polyfit(ms_mass[arms_medium],ms_sfr[arms_medium],1,full=True)
    yarr_medium = np.polyval([a_medium,b_medium],xarr)
    (a_loose,b_loose),res_loose,_,_,_ = np.polyfit(ms_mass[arms_loose],ms_sfr[arms_loose],1,full=True)
    yarr_loose = np.polyval([a_loose,b_loose],xarr)

    print 'Spirals, %f, %f, %f' % (a,b,res/len(ms_mass))
    print 'Tight, %f, %f, %f' % (a_tight,b_tight,res_tight/len(arms_tight))
    print 'Medium, %f, %f, %f' % (a_medium,b_medium,res_medium/len(arms_medium))
    print 'Loose, %f, %f, %f' % (a_loose,b_loose,res_loose/len(arms_loose))

    plt.ion()

    fig = plt.figure(23,figsize=(9,8))
    fig.clf()

    # Plot 1

    ax = fig.add_subplot(221)
    ax.plot(xarr,yarr,color='grey',linestyle='--',lw=3)
    ax.plot(xarr,yarr_tight,color='red',linestyle='--',lw=3)
    ax.scatter(ms_mass,ms_sfr,s=1,color='grey')
    ax.scatter(ms_mass[arms_tight],ms_sfr[arms_tight],s=2,color='red')
    ax.set_xlabel('log (Stellar mass ['+r'$M_{sun}$'+'])',fontsize=20)
    ax.set_ylabel('log (SFR [('+r'$M_{sun}$'+'/yr])',fontsize=20)
    ax.set_xlim(6,12)
    ax.set_ylim(-4,2)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(20)

    legfont = FontProperties()
    legfont.set_size('medium')
    plt.legend(('Arms tight','All spirals'), 'upper left', scatterpoints=1, shadow=True, fancybox=True, prop=legfont)
    
    # Plot 2

    ax = fig.add_subplot(222)
    ax.plot(xarr,yarr,color='grey',linestyle='--',lw=3)
    ax.plot(xarr,yarr_medium,color='blue',linestyle='--',lw=3)
    ax.scatter(ms_mass,ms_sfr,s=1,color='grey')
    ax.scatter(ms_mass[arms_medium],ms_sfr[arms_medium],s=2,color='blue')
    ax.set_xlabel('log (Stellar mass ['+r'$M_{sun}$'+'])',fontsize=20)
    ax.set_ylabel('log (SFR [('+r'$M_{sun}$'+'/yr])',fontsize=20)
    ax.set_xlim(6,12)
    ax.set_ylim(-4,2)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(20)

    legfont = FontProperties()
    legfont.set_size('medium')
    plt.legend(('Arms medium','All spirals'), 'upper left', scatterpoints=1, shadow=True, fancybox=True, prop=legfont)
    
    # Plot 3

    ax = fig.add_subplot(223)
    ax.plot(xarr,yarr,color='grey',linestyle='--',lw=3)
    ax.plot(xarr,yarr_loose,color='green',linestyle='--',lw=3)
    ax.scatter(ms_mass,ms_sfr,s=1,color='grey')
    ax.scatter(ms_mass[arms_loose],ms_sfr[arms_loose],s=2,color='green')
    ax.set_xlabel('log (Stellar mass ['+r'$M_{sun}$'+'])',fontsize=20)
    ax.set_ylabel('log (SFR [('+r'$M_{sun}$'+'/yr])',fontsize=20)
    ax.set_xlim(6,12)
    ax.set_ylim(-4,2)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(20)

    legfont = FontProperties()
    legfont.set_size('medium')
    plt.legend(('Arms loose','All spirals'), 'upper left', scatterpoints=1, shadow=True, fancybox=True, prop=legfont)
    
    # Plot 4

    ax = fig.add_subplot(224)
    ax.plot(xarr,yarr,color='black',linestyle='-',lw=1)
    ax.plot(xarr,yarr_tight,color='red',linestyle='-',lw=1)
    ax.plot(xarr,yarr_medium,color='blue',linestyle='-',lw=1)
    ax.plot(xarr,yarr_loose,color='green',linestyle='-',lw=1)
    ax.set_xlabel('log (Stellar mass ['+r'$M_{sun}$'+'])',fontsize=20)
    ax.set_ylabel('log (SFR [('+r'$M_{sun}$'+'/yr])',fontsize=20)
    ax.set_xlim(6,12)
    ax.set_ylim(-4,2)

    ax.text(9.0,-1.5,r'$\alpha,\beta$=%6.3f,%6.3f' % (a,b),color='black')
    ax.text(9.0,-2.0,r'$\alpha,\beta$=%6.3f,%6.3f' % (a_tight,b_tight),color='red')
    ax.text(9.0,-2.5,r'$\alpha,\beta$=%6.3f,%6.3f' % (a_medium,b_medium),color='blue')
    ax.text(9.0,-3.0,r'$\alpha,\beta$=%6.3f,%6.3f' % (a_loose,b_loose),color='green')

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(20)

    legfont = FontProperties()
    legfont.set_size('medium')
    plt.legend(('All spirals','Tight','Medium','Loose'), 'upper left', scatterpoints=1, shadow=True, fancybox=True, prop=legfont)

    plt.draw()
    plt.show()
    
    fig.tight_layout()
    if savefig:
        fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/zurich/mainsequence_spiraltightness.pdf', dpi=200)

    return None

def get_lackner():


    p = pyfits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/fits/lackner_mpa_gz2.fits')
    lackner = p[1].data
    p.close()

    obliques = (lackner['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.430) & \
               (lackner['t02_edgeon_a05_no_debiased'] >= 0.715) & \
               (lackner['t04_spiral_a08_spiral_debiased'] >= 0.619) & \
               (lackner['t04_spiral_a08_spiral_weight'] >= 20)

    starforming = (lackner['bpt'] == 1)# | (lackner['bpt'] == 2)

    return lackner, obliques

def plot_bt_dispersion(savefig=False):

    lackner, obliques = get_lackner()

    p_medium = lackner['t10_arms_winding_a29_medium_debiased']
    p_2armed = lackner['t11_arms_number_a32_2_debiased']
    p_bar = lackner['t03_bar_a06_bar_debiased']
    p_features = lackner['t01_smooth_or_features_a02_features_or_disk_debiased']

    # Split by mass range

    mass_lo = lackner['MEDIAN_MASS'] < 9.5
    mass_md = (lackner['MEDIAN_MASS'] >= 9.5) & (lackner['MEDIAN_MASS'] <= 10.5)
    mass_hi = lackner['MEDIAN_MASS'] > 10.5

    # Split by B/T ratio

    bt_q1 = lackner['bulge_to_tot_r'] <= 0.25
    bt_q2 = (lackner['bulge_to_tot_r'] > 0.25) & (lackner['bulge_to_tot_r'] <= 0.50)
    bt_q3 = (lackner['bulge_to_tot_r'] > 0.50) & (lackner['bulge_to_tot_r'] <= 0.75)
    bt_q4 = lackner['bulge_to_tot_r'] > 0.75

    plt.ion()
    fig = plt.figure(2,figsize=(9,8))
    fig.clf()

    fs = 10

    # Plots

    mass_arr = (mass_lo,mass_md,mass_hi)
    bt_arr = (bt_q1,bt_q2,bt_q3,bt_q4)
    mass_arr_text = ('Low','Medium','High')
    bt_arr_text = ('1st','2nd','3rd','4th')

    for im,m in enumerate(mass_arr):
        for ib,b in enumerate(bt_arr):

            ax = fig.add_subplot(3,4,im*4 + ib + 1)
            #ax.hist((p_medium[m & b])[p_medium[m & b] > 0.],bins=10,histtype='stepfilled',alpha=0.6,color='green')
            #ax.hist((p_2armed[m & b])[p_2armed[m & b] > 0.],bins=10,histtype='stepfilled',alpha=0.6,color='red')
            #ax.hist((p_bar[m & b])[p_bar[m & b] > 0.],bins=10,histtype='stepfilled',alpha=0.6,color='blue')
            ax.hist((p_features[m & b]),bins=10,histtype='stepfilled',alpha=0.6,color='purple')
            ax.set_title('%s mass, B/T %s quartile' % (mass_arr_text[im],bt_arr_text[ib]),fontsize=fs*0.8)
            ax.set_xlabel(r'$p_{GZ2}$',fontsize=fs)
            ax.set_xlim(-0.05,1.05)
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(fs)

    plt.draw()
    plt.show()
    
    fig.tight_layout()
    if savefig:
        fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/zurich/bulgeratio_morph.pdf', dpi=200)

    return None


