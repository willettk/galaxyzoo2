# coding: utf-8

# 24 May 2013
# Scripts to compute the level of agreement between GZ2 and various other catalogs; makes Table 3 in paper


p = gz2.pyfits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/fits/hc_gz2_debiased.fits')
efigi = p[1].data
p.close()

gz2merger = (efigi['t06_odd_a14_yes_debiased'] >= 0.42) & (efigi['t06_odd_a14_yes_weight'] >= 20) & (efigi['t08_odd_feature_a24_merger_debiased'] >= 0.8)
efigi_gz2_merger = efigi[gz2merger]
multiplicity = (efigi_gz2_merger['Multiplicity'] >= 0.25) & np.logical_not((efigi_gz2_merger['Multiplicity_inf'] == 0.) & (efigi_gz2_merger['Multiplicity_sup'] == 1.))
contamination = (efigi_gz2_merger['Contamination'] >= 0.50) & np.logical_not((efigi_gz2_merger['Contamination_inf'] == 0.) & (efigi_gz2_merger['Contamination_sup'] == 1.))
perturbation = (efigi_gz2_merger['Perturbation'] >= 0.50) & np.logical_not((efigi_gz2_merger['Perturbation_inf'] == 0.) & (efigi_gz2_merger['Perturbation_sup'] == 1.))
allmergers = multiplicity | contamination | perturbation
print np.sum(allmergers),np.sum(gz2merger),np.sum(allmergers).astype(float)/np.sum(gz2merger)


gz2early = (hc['t01_smooth_or_features_a01_smooth_debiased'] >= 0.5) & (hc['total_classifications'] >= 20)
hcearly = (hc['pE/S0'][gz2early] >= 0.5)
print np.sum(hcearly),np.sum(gz2early),np.sum(hcearly).astype(float)/np.sum(gz2early)

gz2early_clean = (hc['t01_smooth_or_features_a01_smooth_debiased'] >= 0.8) & (hc['total_classifications'] >= 20)
hcearly_clean = (hc['pE/S0'][gz2early_clean] >= 0.5)
print np.sum(hcearly_clean),np.sum(gz2early_clean),np.sum(hcearly_clean).astype(float)/np.sum(gz2early_clean)

gz2late = (hc['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5) & (hc['total_classifications'] >= 20)
hclate = (hc['pE/S0'][gz2late] <= 0.5)
print np.sum(hclate),np.sum(gz2late),np.sum(hclate).astype(float)/np.sum(gz2late)

gz2late_clean = (hc['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.8) & (hc['total_classifications'] >= 20)
hclate_clean = (hc['pE/S0'][gz2late_clean] <= 0.5)
print np.sum(hclate_clean),np.sum(gz2late_clean),np.sum(hclate_clean).astype(float)/np.sum(gz2late_clean)


import galaxyzoo2 as gz2
import numpy as np
file='/Users/willettk/Astronomy/Research/GalaxyZoo/fits/gz1_gz2_debiased.fits'
p = gz2.pyfits.open(file)
gz1gz2 = p[1].data
p.close()
gz1gz2.names
gz1gz2['t01_smooth_or_features_a02_features_or_disk'] >= 0.5
gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5
np.sum(gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5)
gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5
gz2late = gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5
gz1gz2[gz2late]
gz1gz2['P_CS'][gz2late]
np.mean(gz1gz2['P_CS'][gz2late])
gz1gz2['P_CS'][gz2late] >= 0.5
np.sum(gz1gz2['P_CS'][gz2late] >= 0.5)
np.sum(gz1gz2['P_CS'][gz2late] >= 0.5) / 138022.
gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.8
gz2late_clean = gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.8
np.sum(gz1gz2['P_CS'][gz2late_clean] >= 0.5)
np.sum(gz2late_clean)
np.sum(gz1gz2['P_CS'][gz2late_clean] >= 0.5) / np.sum(gz2late_clean)
np.sum(gz1gz2['P_CS'][gz2late_clean] >= 0.5) / np.sum(gz2late_clean).astype(float)
gz2early = gz1gz2['t01_smooth_or_features_a01_smooth_debiased'] >= 0.5
gz2early_clean = gz1gz2['t01_smooth_or_features_a01_smooth_debiased'] >= 0.8
np.sum(gz2early),np.sum(gz2early_clean)
np.sum(gz1gz2['P_EL'][gz2early] >= 0.5)
np.sum(gz1gz2['P_EL'][gz2early_clean] >= 0.5)
91934./101153
25772./26314
gz2early = (gz1gz2['t06_odd_a14_debiased'] >= 0.42) & (gz1gz2['t08_odd_feature_a24_merger_debiased'] >= 0.5)
gz2early = (gz1gz2['t06_odd_a14_yes_debiased'] >= 0.42) & (gz1gz2['t08_odd_feature_a24_merger_debiased'] >= 0.5)
gz2merger = (gz1gz2['t06_odd_a14_yes_debiased'] >= 0.42) & (gz1gz2['t08_odd_feature_a24_merger_debiased'] >= 0.5)
gz2early = gz1gz2['t01_smooth_or_features_a01_smooth_debiased'] >= 0.5
gz2merger_clean = (gz1gz2['t06_odd_a14_yes_debiased'] >= 0.42) & (gz1gz2['t08_odd_feature_a24_merger_debiased'] >= 0.8)
np.sum(gz2merger),np.sum(gz2merger_clean)
np.sum(gz1gz2['P_MG'][gz2early_clean] >= 0.5)
np.sum(gz1gz2['P_MG'][gz2merger] >= 0.5)
np.sum(gz1gz2['P_MG'][gz2merger_clean] >= 0.5)
2000./6657
333./526
gz2merger = (gz1gz2['t06_odd_a14_yes_debiased'] >= 0.42) & (gz1gz2['t08_odd_feature_a24_merger_debiased'] >= 0.5) & (gz1gz2['t06_odd_a14_yes_weight'] >= 20)
np.sum(gz2merger),np.sum(gz2merger_clean)
gz2merger_clean = (gz1gz2['t06_odd_a14_yes_debiased'] >= 0.42) & (gz1gz2['t08_odd_feature_a24_merger_debiased'] >= 0.8) & (gz1gz2['t06_odd_a14_yes_weight'] >= 20)
np.sum(gz2merger),np.sum(gz2merger_clean)
np.sum(gz1gz2['P_MG'][gz2merger] >= 0.5)
np.sum(gz1gz2['P_MG'][gz2merger_clean] >= 0.5)
330./487
1989./5331
gz2late = (gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5) & (gz1gz2['total_classifications'] >= 20)
gz2late_clean = (gz1gz2['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.8) & (gz1gz2['total_classifications'] >= 20)
np.sum(gz2late),np.sum(gz2late_clean)
np.sum(gz1gz2['P_CS'][gz2late] >= 0.5)
np.sum(gz1gz2['P_CS'][gz2late_clean] >= 0.5)
np.sum(gz1gz2['P_CS'][gz2late] >= 0.5).astype(float) / 138009.
68241./79203
gz2early = (gz1gz2['t01_smooth_or_features_a01_smooth_debiased'] >= 0.5) & (gz1gz2['total_classifications'] >= 20)
gz2early_clean = (gz1gz2['t01_smooth_or_features_a01_smooth_debiased'] >= 0.8) & (gz1gz2['total_classifications'] >= 20)
np.sum(gz2early),np.sum(gz1gz2['P_CS'][gz2late] >= 0.5).astype(float),np.sum(gz1gz2['P_CS'][gz2late] >= 0.5).astype(float)/np.sum(gz2early)
np.sum(gz2early),np.sum(gz1gz2['P_CS'][gz2late] >= 0.5),np.sum(gz1gz2['P_CS'][gz2late] >= 0.5).astype(float)/np.sum(gz2early)
np.sum(gz2early),np.sum(gz1gz2['P_EL'][gz2early] >= 0.5).astype(float),np.sum(gz1gz2['P_EL'][gz2early] >= 0.5).astype(float)/np.sum(gz2early)
np.sum(gz2early),np.sum(gz1gz2['P_EL'][gz2early] >= 0.5),np.sum(gz1gz2['P_EL'][gz2early] >= 0.5).astype(float)/np.sum(gz2early)
np.sum(gz2early_clean),np.sum(gz1gz2['P_EL'][gz2early_clean] >= 0.5),np.sum(gz1gz2['P_EL'][gz2early_clean] >= 0.5).astype(float)/np.sum(gz2early_clean)
p = gz2.pyfits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/cjl/gz2_nair.fits')
gz2na10 = p[1].data
p.close()
gz2early = (gz2na10['t01_smooth_or_features_a01_smooth_debiased_1'] >= 0.5) & (gz1gz2['total_classifications_1'] >= 20)
gz2early = (gz2na10['t01_smooth_or_features_a01_smooth_debiased'] >= 0.5) & (gz2na10['total_classifications'] >= 20)
gz2early_clean = (gz2na10['t01_smooth_or_features_a01_smooth_debiased'] >= 0.8) & (gz2na10['total_classifications'] >= 20)
np.sum(gz2early),np.sum(gz2early_clean)
np.sum(gz2na10['TType'][gz2early] < 0)
np.sum(gz2na10['TType'][gz2early_clean] < 0)
np.sum(gz2na10['TType'][gz2early_clean] <= 0)
np.sum(gz2na10['TType'][gz2early] <= 0)
1930./1995
4354./4999
gz2late = (gz2na10['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5) & (gz2na10['total_classifications'] >= 20)
gz2late_clean = (gz2na10['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.8) & (gz2na10['total_classifications'] >= 20)
gz2late.sum
gz2late.sum()
gz2late.sum(),gz2late_clean.sum()
np.sum((gz2na10['TType'][gz2late] >= 1) & (gz2na10['TType'][gz2late] <=10))
np.sum((gz2na10['TType'][gz2late_clean] >= 1) & (gz2na10['TType'][gz2late_clean] <=10))
np.sum((gz2na10['TType'][gz2late_clean] >= 1))
np.sum((gz2na10['TType'][gz2late] >= 1))
6624/7515.
5203/5481.
gz2bar = (gz2na10['t03_bar_a06_bar_debiased'] >= 0.5) & (gz2na10['t02_edgeon_a05_no_weight'] >= 20) & (gz2na10['t01_smooth_or_features_a02_features_or_disk_debiased'] > 0.43) & (gz2na10['t02_edgeon_a05_no_debiased'] > 0.715)
gz2bar = (gz2na10['t03_bar_a06_bar_debiased'] >= 0.5) & (gz2na10['t02_edgeon_a05_no_weight_1'] >= 20) & (gz2na10['t01_smooth_or_features_a02_features_or_disk_debiased'] > 0.43) & (gz2na10['t02_edgeon_a05_no_debiased'] > 0.715)
gz2bar_clean = (gz2na10['t03_bar_a06_bar_debiased'] >= 0.8) & (gz2na10['t02_edgeon_a05_no_weight_1'] >= 20) & (gz2na10['t01_smooth_or_features_a02_features_or_disk_debiased'] > 0.43) & (gz2na10['t02_edgeon_a05_no_debiased'] > 0.715)
gz2bar.sum(),gz2bar_clean.sum
gz2bar.sum(),gz2bar_clean.sum()
np.sum((gz2na10['Bar'][gz2bar] > 0))
np.sum((gz2na10['Bar'][gz2bar_clean] > 0))
1235/1446.
618/651.
gz2ring = (gz2na10['t06_odd_a14_yes_debiased'] >= 0.42) & (gz2na10['t06_odd_a14_yes_weight_1'] >= 20) & (gz2na10['t08_odd_feature_a19_ring_debiased'] >= 0.5)
gz2ring_clean = (gz2na10['t06_odd_a14_yes_debiased'] >= 0.42) & (gz2na10['t06_odd_a14_yes_weight_1'] >= 20) & (gz2na10['t08_odd_feature_a19_ring_debiased'] >= 0.8)
gz2ring.sum(),gz2ring_clean.sum()
np.sum((gz2na10['Ring'][gz2ring] > 0))
np.sum((gz2na10['Ring'][gz2ring_clean] > 0))
401./438
487/549.
gz2ring = (gz2na10['t06_odd_a14_yes_debiased'] >= 0.42) & (gz2na10['t06_odd_a14_yes_weight_1'] >= 20) & (gz2na10['t08_odd_feature_a24_merger_debiased'] >= 0.5)
gz2merger = (gz2na10['t06_odd_a14_yes_debiased'] >= 0.42) & (gz2na10['t06_odd_a14_yes_weight_1'] >= 20) & (gz2na10['t08_odd_feature_a24_merger_debiased'] >= 0.5)
gz2merger_clean = (gz2na10['t06_odd_a14_yes_debiased'] >= 0.42) & (gz2na10['t06_odd_a14_yes_weight_1'] >= 20) & (gz2na10['t08_odd_feature_a24_merger_debiased'] >= 0.8)
gz2merger.sum(),gz2merger_clean.sum()
np.sum((gz2na10['Pair'][gz2ring] > 0) | (gz2na10['dist'][gz2ring] > 0)
)
np.sum((gz2na10['Pair'][gz2ring] > 0) | (gz2na10['dist'][gz2ring] > 0))
np.sum((gz2na10['Pair'][gz2ring_clean] > 0) | (gz2na10['dist'][gz2ring_clean] > 0))
np.sum((gz2na10['Pair'][gz2merger] > 0) | (gz2na10['dist'][gz2merger] > 0))
np.sum((gz2na10['Pair'][gz2merger_clean] > 0) | (gz2na10['dist'][gz2merger_clean] > 0))
242/253.
p = gz2.pyfits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/cjl/gz2_efigi_debiased.fits')
efigi = p[1].data
p.close()
gz2early = (efigi['t01_smooth_or_features_a01_smooth_debiased'] >= 0.5) & (efigi['total_classifications'] >= 20)
gz2early_clean = (efigi['t01_smooth_or_features_a01_smooth_debiased'] >= 0.8) & (efigi['total_classifications'] >= 20)
gz2early.sum(),gz2early_clean.sum()
np.sum((efigi['T'][g2early] > 0))
np.sum((efigi['T'][gx2early] > 0))
np.sum((efigi['T'][gz2early] > 0))
np.sum((efigi['T'][gz2early_clean] > 0))
np.sum((efigi['T'][gz2early] < 0))
np.sum((efigi['T'][gz2early_clean] < 0))
400/519.
181/214.
gz2late = (efigi['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.5) & (efigi['total_classifications'] >= 20)
gz2late_clean = (efigi['t01_smooth_or_features_a02_features_or_disk_debiased'] >= 0.8) & (efigi['total_classifications'] >= 20)
np.sum((efigi['T'][gz2late] >= 0))
np.sum((efigi['T'][gz2late_clean] >= 0))
gz2late.sum(),gz2late_clean.sum()
2073/2160.
1645/1675.
gz2bar = (efigi['t03_bar_a06_bar_debiased'] >= 0.5) & (efigi['t02_edgeon_a05_no_weight_1'] >= 20) & (efigi['t01_smooth_or_features_a02_features_or_disk_debiased'] > 0.43) & (efigi['t02_edgeon_a05_no_debiased'] > 0.715)
gz2bar = (efigi['t03_bar_a06_bar_debiased'] >= 0.5) & (efigi['t02_edgeon_a05_no_weight'] >= 20) & (efigi['t01_smooth_or_features_a02_features_or_disk_debiased'] > 0.43) & (efigi['t02_edgeon_a05_no_debiased'] > 0.715)
gz2bar_clean = (efigi['t03_bar_a06_bar_debiased'] >= 0.8) & (efigi['t02_edgeon_a05_no_weight'] >= 20) & (efigi['t01_smooth_or_features_a02_features_or_disk_debiased'] > 0.43) & (efigi['t02_edgeon_a05_no_debiased'] > 0.715)
efigibar = ((efigi['Bar_Length'][gz2late] >= 0.5) & np.logical_not((efigi['Bar_Length_inf'][gz2late] == 0.) & (efigi['Bar_Length_sup'][gz2late] == 1.))
)
efigibar = ((efigi['Bar_Length'][gz2late] >= 0.5) & np.logical_not((efigi['Bar_Length_inf'][gz2late] == 0.) & (efigi['Bar_Length_sup'][gz2late] == 1.)))
np.sum(efigibar)
np.sum(efigi['Bar_Length'][gz2late] >= 0.5)
efigibar_clean = ((efigi['Bar_Length'][gz2late_clean] >= 0.5) & np.logical_not((efigi['Bar_Length_inf'][gz2late_clean] == 0.) & (efigi['Bar_Length_sup'][gz2late_clean] == 1.)))
np.sum(efigibar),np.sum(gz2bar),np.sum(efigibar)/np.sum(gz2bar).astype(float)
efigibar = ((efigi['Bar_Length'][gz2bar] >= 0.5) & np.logical_not((efigi['Bar_Length_inf'][gz2bar] == 0.) & (efigi['Bar_Length_sup'][gz2bar] == 1.)))
efigibar_clean = ((efigi['Bar_Length'][gz2bar_clean] >= 0.5) & np.logical_not((efigi['Bar_Length_inf'][gz2bar_clean] == 0.) & (efigi['Bar_Length_sup'][gz2bar_clean] == 1.)))
np.sum(efigibar),np.sum(gz2bar),np.sum(efigibar)/np.sum(gz2bar).astype(float)
np.sum(efigibar_clean),np.sum(gz2bar_clean),np.sum(efigibar_clean)/np.sum(gz2bar_clean).astype(float)
gz2ring = (efigi['t06_odd_a14_yes_debiased'] >= 0.42) & (efigi['t06_odd_a14_yes_weight'] >= 20) & (efigi['t08_odd_feature_a24_merger_debiased'] >= 0.5)
gz2ring_clean = (efigi['t06_odd_a14_yes_debiased'] >= 0.42) & (efigi['t06_odd_a14_yes_weight'] >= 20) & (efigi['t08_odd_feature_a19_ring_debiased'] >= 0.8)
gz2ring = (efigi['t06_odd_a14_yes_debiased'] >= 0.42) & (efigi['t06_odd_a14_yes_weight'] >= 20) & (efigi['t08_odd_feature_a19_ring_debiased'] >= 0.5)
efigi_gz2_ring = efigi[gz2ring]
efigi_gz2_ring_clean = efigi[gz2ring_clean]
iring = (efigi_gz2_ring['Inner_Ring'] >= 0.5) & np.logical_not((efigi_gz2_ring['Inner_Ring_inf'] == 0.) & (efigi_gz2_ring['Inner_Ring_sup'] == 1.))
oring = (efigi_gz2_ring['Outer_Ring'] >= 0.5) & np.logical_not((efigi_gz2_ring['Outer_Ring_inf'] == 0.) & (efigi_gz2_ring['Outer_Ring_sup'] == 1.))
pring = (efigi_gz2_ring['Pseudo_Ring'] >= 0.5) & np.logical_not((efigi_gz2_ring['Pseudo_Ring_inf'] == 0.) & (efigi_gz2_ring['Pseudo_Ring_sup'] == 1.))
allrings = iring | oring | pring
np.sum(allrings),np.sum(gz2ring),np.sum(allrings).astype(float)/np.sum(gz2ring)
iring_clean = (efigi_gz2_ring_clean['Inner_Ring'] >= 0.5) & np.logical_not((efigi_gz2_ring_clean['Inner_Ring_inf'] == 0.) & (efigi_gz2_ring_clean['Inner_Ring_sup'] == 1.))
oring_clean = (efigi_gz2_ring_clean['Outer_Ring'] >= 0.5) & np.logical_not((efigi_gz2_ring_clean['Outer_Ring_inf'] == 0.) & (efigi_gz2_ring_clean['Outer_Ring_sup'] == 1.))
pring_clean = (efigi_gz2_ring_clean['Pseudo_Ring'] >= 0.5) & np.logical_not((efigi_gz2_ring_clean['Pseudo_Ring_inf'] == 0.) & (efigi_gz2_ring_clean['Pseudo_Ring_sup'] == 1.))
allrings_clean = iring_clean | oring_clean | pring_clean
np.sum(allrings_clean),np.sum(gz2ring_clean),np.sum(allrings_clean).astype(float)/np.sum(gz2ring_clean)
get_ipython().magic(u'paste')
get_ipython().magic(u'paste')
len(efigi)
len(efigi) * 0.02
p = gz2.pyfits.open('/Users/willettk/Astronomy/Research/GalaxyZoo/fits/hc_gz2_debiased.fits')
hc = p[1].data
p.close()
get_ipython().magic(u'paste')
get_ipython().magic(u'paste')
get_ipython().magic(u'save /Users/willettk/Astronomy/Research/GalaxyZoo/python/comparing_catalogs_history.py')
