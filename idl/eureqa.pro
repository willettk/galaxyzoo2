
;+
; NAME:
;       
;	EUREQA
;
; PURPOSE:
;
;	
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
;       Written by K. Willett                5 Mar 13
;-

pro eureqa, stop=stop, ps=ps

fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
csvdir = '~/Astronomy/Research/GalaxyZoo/csv/'
figdir = '~/Astronomy/Research/GalaxyZoo/datapaper/figures/'

; NA10 data

nagz2=mrdfits(fitsdir+'na_gz2_willettk.fit',1,/silent)

na_gz_features_or_disk_wfraction = nagz2.T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION
na_gz_not_edgeon_wfraction = nagz2.T02_EDGEON_A05_NO_WEIGHTED_FRACTION
na_gz_spiral_wfraction = nagz2.T04_SPIRAL_A08_SPIRAL_WEIGHTED_FRACTION
na_gz_tight_wfraction = nagz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION
na_gz_medium_wfraction = nagz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION
na_gz_loose_wfraction = nagz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION
na_gz_no_bulge_wfraction = nagz2.T05_BULGE_PROMINENCE_A10_NO_BULGE_WEIGHTED_FRACTION
na_gz_just_noticeable_wfraction = nagz2.T05_BULGE_PROMINENCE_A11_JUST_NOTICEABLE_WEIGHTED_FRACTION
na_gz_obvious_wfraction = nagz2.T05_BULGE_PROMINENCE_A12_OBVIOUS_WEIGHTED_FRACTION
na_gz_dominant_wfraction = nagz2.T05_BULGE_PROMINENCE_A13_DOMINANT_WEIGHTED_FRACTION
na_ttype = nagz2.ttype
na_totalbc = nagz2.T05_BULGE_PROMINENCE_TOTAL_COUNT
na_totalsc = nagz2.T10_ARMS_WINDING_TOTAL_COUNT

; EFIGI data

;efigi=mrdfits(fitsdir+'efigi_gz2_willettk.fit',1,/silent)
;
;efigi_gz_tight_wfraction = efigi.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION
;efigi_gz_medium_wfraction = efigi.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION
;efigi_gz_loose_wfraction = efigi.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION
;efigi_gz_no_bulge_wfraction = efigi.T05_BULGE_PROMINENCE_A10_NO_BULGE_WEIGHTED_FRACTION
;efigi_gz_just_noticeable_wfraction = efigi.T05_BULGE_PROMINENCE_A11_JUST_NOTICEABLE_WEIGHTED_FRACTION
;efigi_gz_obvious_wfraction = efigi.T05_BULGE_PROMINENCE_A12_OBVIOUS_WEIGHTED_FRACTION
;efigi_gz_dominant_wfraction = efigi.T05_BULGE_PROMINENCE_A13_DOMINANT_WEIGHTED_FRACTION
;efigi_ttype = efigi.t
;efigi_totalbc = efigi.T05_BULGE_PROMINENCE_TOTAL_COUNT
;efigi_totalsc = efigi.T10_ARMS_WINDING_TOTAL_COUNT

bulgecount = 10
spiralcount = 10

; Goal: see if there is a linear combination of parameters for the GZ2 data that can predict T-type.

; Try it on only the spiral sequence (S0 - Sd)

; Only include galaxies verified as features or disk: f_features > 0.5, f_notedgeon > 0.5, n_notedgeon > 10

; Three sets of parameters:
; bulge only
; spiral only
; bulge and spiral

; Construct parameter for bulge-only, then histogram by true (NA10) T-Types

; First run: include E and S0 (T-Type -5 to 0)

; Stopped run after 6m25s
; Best fit for linear combination

!p.multi=[0,1,2]

if keyword_set(ps) then begin
	ps_start, filename=figdir+'eureqa.eps',/color, /quiet, /encap, xsize=4,ysize=6
	cs=1.2
	th=3
	thickline=5
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

a10 = na_gz_no_bulge_wfraction
a11 = na_gz_just_noticeable_wfraction 
a12 = na_gz_obvious_wfraction
a13 = na_gz_dominant_wfraction
b_c = -2.44
b_a10 = 10.7
b_a11 = 7.14
b_a12 = 4.3
b_a13 = -1.

b_fit = b_c + b_a10 * a10 + b_a11 * a11 + b_a12 * a12 + b_a13 * a13

e_ind  = where(na_gz_features_or_disk_wfraction ge 0.5 and na_gz_not_edgeon_wfraction ge 0.5 and na_totalbc ge 10 and na_ttype eq -5)
s0_ind = where(na_gz_features_or_disk_wfraction ge 0.5 and na_gz_not_edgeon_wfraction ge 0.5 and na_totalbc ge 10 and na_ttype ge -3 and na_ttype le 0)
sa_ind = where(na_gz_features_or_disk_wfraction ge 0.5 and na_gz_not_edgeon_wfraction ge 0.5 and na_totalbc ge 10 and na_ttype ge  1 and na_ttype le 2)
sb_ind = where(na_gz_features_or_disk_wfraction ge 0.5 and na_gz_not_edgeon_wfraction ge 0.5 and na_totalbc ge 10 and na_ttype ge  3 and na_ttype le 4)
sc_ind = where(na_gz_features_or_disk_wfraction ge 0.5 and na_gz_not_edgeon_wfraction ge 0.5 and na_totalbc ge 10 and na_ttype ge  5 and na_ttype le 6)
sd_ind = where(na_gz_features_or_disk_wfraction ge 0.5 and na_gz_not_edgeon_wfraction ge 0.5 and na_totalbc ge 10 and na_ttype ge  7 and na_ttype le 8)

;b_binsize = 0.5
;cghistoplot, b_fit[e_ind], $
;	binsize = b_binsize, $
;	/outline, $
;    xr=[-3,10], $
;    yr=[0,700], $
;	xtitle='GZ2 combined bulge prominence', $
;	datacolor='black'
;
;cghistoplot, b_fit[s0_ind], $
;	binsize = b_binsize, $
;	/outline, $
;    /oplot, $
;	datacolor='tomato'
;
;cghistoplot, b_fit[sa_ind], $
;	binsize = b_binsize, $
;	/outline, $
;    /oplot, $
;	datacolor='yellow'
;
;cghistoplot, b_fit[sb_ind], $
;	binsize = b_binsize, $
;	/outline, $
;    /oplot, $
;	datacolor='green'
;
;cghistoplot, b_fit[sc_ind], $
;	binsize = b_binsize, $
;	/outline, $
;    /oplot, $
;	datacolor='dark green'
;
;cghistoplot, b_fit[sd_ind], $
;	binsize = b_binsize, $
;	/outline, $
;    /oplot, $
;	datacolor='blue'
;
;cgplot, b_fit[e_ind], na_ttype[e_ind], $
;    xr=[-3,10], $
;    yr=[-3,10], $
;    xtitle='Predicted T-type', $
;    ytitle='Actual NA10 T-type', $
;    psym = 16, $
;    color='black'
;
;cgplot, b_fit[s0_ind], na_ttype[s0_ind], $
;    /overplot, $
;    psym = 16, $
;    color='tomato'
;
;cgplot, b_fit[sa_ind], na_ttype[sa_ind], $
;    /overplot, $
;    psym = 16, $
;    color='yellow'
;
;cgplot, b_fit[sb_ind], na_ttype[sb_ind], $
;    /overplot, $
;    psym = 16, $
;    color='green'
;
;cgplot, b_fit[sc_ind], na_ttype[sc_ind], $
;    /overplot, $
;    psym = 16, $
;    color='dark green'
;
;cgplot, b_fit[sd_ind], na_ttype[sd_ind], $
;    /overplot, $
;    psym = 16, $
;    color='blue'

; Not a very good separation at all. Formulize also crashed when trying to save report. Damn.

; Try adding spiral prominence; smaller sample of galaxies, since adding "yes" to spiral question
; After 10m and 92300 generations (4.3e9 evaluations), best fit is:
; T = 4.63 + 4.17*a10 - 2.27*a12 - 8.38*a13
; No significant improvements in adding any of the spiral tightness terms.

bs_c = 4.63
bs_a10 = 4.17
bs_a12 = -2.27
bs_a13 = -8.38

bs_fit = bs_c + bs_a10 * a10 + bs_a12 * a12 + bs_a13 * a13

bs_binsize = 0.5
cghistoplot, bs_fit[e_ind], $
	binsize = bs_binsize, $
	/outline, $
    thick = th, $
    charsize=cs, $
    xr=[-6,10], /xstyle, $
    yr=[0,700], $
	xtitle='Predicted T-type', $
	ytitle='Number', $
	datacolor='black'

cghistoplot, bs_fit[s0_ind], $
	binsize = bs_binsize, $
    thick = th, $
	/outline, $
    /oplot, $
	datacolor='tomato'

cghistoplot, bs_fit[sa_ind], $
	binsize = bs_binsize, $
    thick = th, $
	/outline, $
    /oplot, $
	datacolor='yellow'

cghistoplot, bs_fit[sb_ind], $
	binsize = bs_binsize, $
    thick = th, $
	/outline, $
    /oplot, $
	datacolor='green'

cghistoplot, bs_fit[sc_ind], $
	binsize = bs_binsize, $
    thick = th, $
	/outline, $
    /oplot, $
	datacolor='dark green'

cghistoplot, bs_fit[sd_ind], $
	binsize = bs_binsize, $
    thick = th, $
	/outline, $
    /oplot, $
	datacolor='blue'

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

cgplot, bs_fit[e_ind], na_ttype[e_ind], $
    xr=[-6,10], /xstyle, $
    yr=[-6,10], /ystyle, $
    thick = th, $
    charsize = cs, $
    xtitle='Predicted T-type', $
    ytitle='NA10 T-type', $
    psym = 16, $
    color='black'

cgplot, bs_fit[s0_ind], na_ttype[s0_ind], $
    /overplot, $
    thick = th, $
    psym = 16, $
    color='tomato'

cgplot, bs_fit[sa_ind], na_ttype[sa_ind], $
    /overplot, $
    thick = th, $
    psym = 16, $
    color='yellow'

cgplot, bs_fit[sb_ind], na_ttype[sb_ind], $
    /overplot, $
    psym = 16, $
    thick = th, $
    color='green'

cgplot, bs_fit[sc_ind], na_ttype[sc_ind], $
    /overplot, $
    psym = 16, $
    thick = th, $
    color='dark green'

cgplot, bs_fit[sd_ind], na_ttype[sd_ind], $
    /overplot, $
    psym = 16, $
    thick = th, $
    color='blue'

cgplots, !x.crange, !y.crange, $
    linestyle=2, $
    thick = th, $
    color='black'

if keyword_set(ps) then ps_end

stop

; Instead of fitting the GZ2 fractions for T-types, try instead looking at average distributions of
; GZ2 fractions for a given T-type

!p.multi=[0,1,6]

bsize = 0.05

cghistoplot, a10[e_ind], title='E galaxies', binsize = bsize, /outline, /freq, xr=[0,1], yr=[0,1], xtitle='GZ2 vote fraction', datacolor='red',thick=th, ytickformat='(f3.1)', charsize=2
cghistoplot, a11[e_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='blue'
cghistoplot, a12[e_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='orange'
cghistoplot, a13[e_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='green'

cghistoplot, a10[s0_ind], title='S0 galaxies', binsize = bsize, /outline, /freq, xr=[0,1], yr=[0,1], xtitle='GZ2 vote fraction', datacolor='red',thick=th, ytickformat='(f3.1)', charsize=2
cghistoplot, a11[s0_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='blue'
cghistoplot, a12[s0_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='orange'
cghistoplot, a13[s0_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='green'

cghistoplot, a10[sa_ind], title='Sa galaxies', binsize = bsize, /outline, /freq, xr=[0,1], yr=[0,1], xtitle='GZ2 vote fraction', datacolor='red',thick=th, ytickformat='(f3.1)', charsize=2
cghistoplot, a11[sa_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='blue'
cghistoplot, a12[sa_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='orange'
cghistoplot, a13[sa_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='green'

cghistoplot, a10[sb_ind], title='Sb galaxies', binsize = bsize, /outline, /freq, xr=[0,1], yr=[0,1], xtitle='GZ2 vote fraction', datacolor='red',thick=th, ytickformat='(f3.1)', charsize=2
cghistoplot, a11[sb_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='blue'
cghistoplot, a12[sb_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='orange'
cghistoplot, a13[sb_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='green'

cghistoplot, a10[sc_ind], title='Sc galaxies', binsize = bsize, /outline, /freq, xr=[0,1], yr=[0,1], xtitle='GZ2 vote fraction', datacolor='red',thick=th, ytickformat='(f3.1)', charsize=2
cghistoplot, a11[sc_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='blue'
cghistoplot, a12[sc_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='orange'
cghistoplot, a13[sc_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='green'

cghistoplot, a10[sd_ind], title='Sd galaxies', binsize = bsize, /outline, /freq, xr=[0,1], yr=[0,1], xtitle='GZ2 vote fraction', datacolor='red',thick=th, ytickformat='(f3.1)', charsize=2
cghistoplot, a11[sd_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='blue'
cghistoplot, a12[sd_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='orange'
cghistoplot, a13[sd_ind], binsize=bsize, /outline, /freq, /oplot, thick=th, datacolor='green'

print,mean(a11[e_ind]/a12[e_ind],/nan),stddev(a11[e_ind]/a12[e_ind],/nan)
print,mean(a11[s0_ind]/a12[s0_ind],/nan),stddev(a11[s0_ind]/a12[s0_ind],/nan)
print,mean(a11[sa_ind]/a12[sa_ind],/nan),stddev(a11[sa_ind]/a12[sa_ind],/nan)
print,mean(a11[sb_ind]/a12[sb_ind],/nan),stddev(a11[sb_ind]/a12[sb_ind],/nan)
print,mean(a11[sc_ind]/a12[sc_ind],/nan),stddev(a11[sc_ind]/a12[sc_ind],/nan)
print,mean(a11[sd_ind]/a12[sd_ind],/nan),stddev(a11[sd_ind]/a12[sd_ind],/nan)

if keyword_set(stop) then stop

end
