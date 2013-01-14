
;+
; NAME:
;       
;	MORAN_DEMO
;
; PURPOSE:
;
;	Plot some example GZ2 (not debiased) graphs to show Ed Moran
;
; INPUTS:
;
;	
;
; OUTPUTS:
;
;	
;
; KEYWORDS:
;
;	
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Oct 12
;-

device,retain=2

;gz2main = mrdfits('~/Astronomy/Research/GalaxyZoo/fits/gz2main_table_sample.fits',1,/silent)
;stripe82 = mrdfits('~/Astronomy/Research/GalaxyZoo/fits/coadd1_coadd2_sample.fits',1,/silent)

mpc80 = 0.01866d
mpc50 = 0.01174d
mpc20 = 0.00473d

lowz1 = where(gz2main.redshift ge mpc20 and gz2main.redshift lt mpc50)
lowz2 = where(gz2main.redshift ge mpc50 and gz2main.redshift lt mpc80)

bs = 0.05

ps_start, filename='~/Astronomy/Research/GalaxyZoo/dwarfagn/moran_demo.eps', /color, /quiet, xs=6, ys=8, /encap
!p.multi=[0,1,3]

tcount = 20

cs = 2
ls = cs * 0.5

task01 = where(gz2main.t01_smooth_or_features_total_weight ge tcount)
d1 = gz2main[setintersection(lowz1,task01)].t01_smooth_or_features_a01_smooth_weighted_fraction
d2 = gz2main[setintersection(lowz2,task01)].t01_smooth_or_features_a01_smooth_weighted_fraction

cghistoplot, d1, $
	title='Smooth (Task 01)', $
	xtitle='Weighted vote fraction', $
	charsize = cs, $
	datacolor='red', $
	/oprob, $
	;yr=[0,1], $
	probcolor='red', $
	thick=3, $
	/outline, $
	/freq, $
	binsize=bs

cghistoplot, d2, $
	datacolor='blue', $
	/oprob, $
	probcolor='blue', $
	/oplot, $
	/outline, $
	/freq, $
	binsize=bs

al_legend, pos=[0.3,0.95], /normal, ['20 - 50 Mpc', '50 - 80 Mpc'], color=['red','blue'], linestyle=0, thick=[3,3], charsize=ls

kstwo, d1, d2, dks, probks
print,'Task 01 KS prob: ', probks

task03 = where(gz2main.t03_bar_total_weight ge tcount)

d1 = gz2main[setintersection(lowz1,task03)].t03_bar_a06_bar_weighted_fraction
d2 = gz2main[setintersection(lowz2,task03)].t03_bar_a06_bar_weighted_fraction

cghistoplot, d1, $
	title='Barred (Task 03)', $
	xtitle='Weighted vote fraction', $
	charsize = cs, $
	datacolor='red', $
	/oprob, $
	;yr=[0,1], $
	probcolor='red', $
	thick=3, $
	/outline, $
	/freq, $
	binsize=bs

cghistoplot, d2, $
	datacolor='blue', $
	/oprob, $
	probcolor='blue', $
	/oplot, $
	/outline, $
	/freq, $
	binsize=bs

;al_legend, /top, /right, ['20 - 50 Mpc', '50 - 80 Mpc'], color=['red','blue'], linestyle=0, thick=[3,3], charsize=ls

kstwo, d1, d2, dks, probks
print,'Task 03 KS prob: ', probks

task05 = where(gz2main.t05_bulge_prominence_total_weight ge tcount)

d1 = gz2main[setintersection(lowz1,task05)].t05_bulge_prominence_a12_obvious_weighted_fraction + $
	gz2main[setintersection(lowz1,task05)].t05_bulge_prominence_a13_dominant_weighted_fraction
d2 = gz2main[setintersection(lowz2,task05)].t05_bulge_prominence_a12_obvious_weighted_fraction + $
	gz2main[setintersection(lowz2,task05)].t05_bulge_prominence_a13_dominant_weighted_fraction

cghistoplot, d1, $
	title='Obvious or dominant bulge (Task 05)', $
	xtitle='Weighted vote fraction', $
	charsize = cs, $
	datacolor='red', $
	/oprob, $
;	yr=[0,1], $
	probcolor='red', $
	thick=3, $
	/outline, $
	/freq, $
	binsize=bs

cghistoplot, d2, $
	datacolor='blue', $
	/oplot, $
	/oprob, $
	probcolor='blue', $
	/outline, $
	/freq, $
	binsize=bs

;al_legend, /top, /right, ['20 - 50 Mpc', '50 - 80 Mpc'], color=['red','blue'], linestyle=0, thick=[3,3], charsize=ls

kstwo, d1, d2, dks, probks
print,'Task 05 KS prob: ', probks

ps_end

end
