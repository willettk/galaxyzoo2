
;+
; NAME:
;       
;	EFIGI_GZ2
;
; PURPOSE:
;
;	Compare the EFIGI and GZ2 catalogs
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
;       Written by K. Willett                May 2012
;-

pro efigi_gz2, ps=ps, stop=stop, count=count, gz2armcount=gz2armcount

; Bars

file='~/Astronomy/Research/GalaxyZoo/csv/bars_efigi_gz2_axial_willettk.csv'
readcol, file, $
	objid,$
	efigi_bar,efigi_bar_inf,efigi_bar_sup,$
	gz2_bar_wfraction,gz2_bar_count,$
	gz2_nobar_wfraction,gz2_nobar_count,$
	gz2_edgeon_wfraction,gz2_edgeon_count, $
	expAB_r,expAB_g,$
	format='a,f,f,f,f,i,f,i,f,i,f,f', $
	skip=1,/silent

	; Remove galaxies with no bar attribute defined

	badefigi = where(efigi_bar eq 0.5 and efigi_bar_inf eq 0.0 and efigi_bar_sup eq 1.0, bad_efcount)
	if bad_efcount gt 0 then begin
		good_efind = set_difference(indgen(n_elements(efigi_bar)),badefigi)
		objid = objid[good_efind]
		efigi_bar = efigi_bar[good_efind]
		efigi_bar_inf = efigi_bar_inf[good_efind]
		efigi_bar_sup = efigi_bar_sup[good_efind]
		gz2_bar_wfraction = gz2_bar_wfraction[good_efind]
		gz2_bar_count = gz2_bar_count[good_efind]
		gz2_nobar_wfraction = gz2_nobar_wfraction[good_efind]
		gz2_nobar_count = gz2_nobar_count[good_efind]
		gz2_edgeon_wfraction = gz2_edgeon_wfraction[good_efind]
		gz2_edgeon_count = gz2_edgeon_count[good_efind]
		expAB_r = expAB_r[good_efind]
		expAB_g = expAB_g[good_efind]
		print,strtrim(bad_efcount,2)+' galaxies removed with no EFIGI bar attribute defined.'
	endif

	if n_elements(count) gt 0 then begin
		countind = where(gz2_bar_count + gz2_nobar_count ge count)
		efigi_bar = efigi_bar[countind]
		efigi_bar_sup = efigi_bar_sup[countind]
		efigi_bar_inf = efigi_bar_inf[countind]
		gz2_bar_wfraction = gz2_bar_wfraction[countind]
		counttitle='_'+strtrim(count,2)
	endif else counttitle=''

	title=counttitle

	print,'Total number of overlap galaxies between samples: ',strtrim(n_elements(objid))

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/datapaper/figures/efigi_bars'+counttitle+'.eps',$
		/color, /quiet, /encap, xs=20, ys=7
	cs=1.5
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

!p.multi=[0,2,1]

; Plot 1

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.07,0.1,0.42,0.9], $
	title='All galaxies', $
	xtitle='EFIGI bar length', $
	ytitle='GZ2 bar weighted fraction'

efbin = 0.25
gz2bin = 0.1

barcell=fltarr(5,10)
bccount_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((gz2_bar_wfraction gt gz2val) $
			and (gz2_bar_wfraction le gz2val+gz2bin) $
			and (efigi_bar eq efval),bccount)
		bccount_arr[i,j] = bccount
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_arr[i,j])/max(bccount_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, !x.crange, !y.crange, linestyle=2, thick=th

cgcolorbar, /vertical, /right, position=[0.1,0.43,0.9,0.45],range=[0,max(bccount_arr)]

corr_na_gz = correlate(gz2_bar_wfraction,efigi_bar)
print,''
print,'Correlation (bars, all galaxies): ',corr_na_gz,n_elements(efigi_bar)

	axialcut = where(alog10(1./expab_r) lt 0.3,axcount)

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.55,0.1,0.90,0.9], $
	title='Disks, not edge-on', $
	xtitle='EFIGI bar length', $
	ytitle='GZ2 bar weighted fraction'

efbin = 0.25
gz2bin = 0.1

barcell=fltarr(5,10)
bccount_ax_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((gz2_bar_wfraction gt gz2val) $
			and (gz2_bar_wfraction le gz2val+gz2bin) $
			and (efigi_bar eq efval)$
			and alog10(1./expab_r) lt 0.3, $
			bccount_ax)
		bccount_ax_arr[i,j] = bccount_ax
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_ax_arr[i,j])/max(bccount_ax_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, !x.crange, !y.crange, linestyle=2, thick=th

cgcolorbar, /vertical, /right, position=[0.1,0.93,0.9,0.95],range=[0,max(bccount_ax_arr)]

corr_na_gz = correlate(gz2_bar_wfraction[axialcut],efigi_bar[axialcut])
print,''
print,'Correlation (bars, face-on galaxies): ',corr_na_gz,axcount


	; Weighted mean of the sample

	efigi_weight = ((efigi_bar_sup - efigi_bar_inf) + 0.25)^(-2)
	efigi_err = (efigi_bar_sup - efigi_bar_inf) + 0.25
	print, wmean(efigi_bar,efigi_err)
	print, wmean(efigi_bar[where(gz2_bar_wfraction ge 0.5 and alog10(1./expab_r) lt 0.3)],efigi_err[where(gz2_bar_wfraction ge 0.5 and alog10(1./expab_r) lt 0.3)])

if keyword_set(ps) then ps_end

; Contamination

;    e_pgc,$
;    e_RAdeg,$
;    e_DEdeg,$
;    e_T,$
;    e_T_inf,$
;    e_T_sup,$
;    e_Bulge_to_Total,$
;    e_Bulge_to_Total_inf,$
;    e_Bulge_to_Total_sup,$
;    e_Arm_Strength,$
;    e_Arm_Strength_inf,$
;    e_Arm_Strength_sup,$
;    e_Arm_Curvature,$
;    e_Arm_Curvature_inf,$
;    e_Arm_Curvature_sup,$
;    e_Arm_Rotation,$
;    e_Arm_Rotation_inf,$
;    e_Arm_Rotation_sup,$
;    e_Bar_Length,$
;    e_Bar_Length_inf,$
;    e_Bar_Length_sup,$
;    e_Inner_Ring,$
;    e_Inner_Ring_inf,$
;    e_Inner_Ring_sup,$
;    e_Outer_Ring	,$
;    e_Outer_Ring_inf,$
;    e_Outer_Ring_sup,$
;    e_Pseudo_Ring,$
;    e_Pseudo_Ring_inf,$
;    e_Pseudo_Ring_sup,$
;    e_Perturbation,$
;    e_Perturbation_inf,$
;    e_Perturbation_sup,$
;    e_Visible_Dust,$
;    e_Visible_Dust_inf,$
;    e_Visible_Dust_sup,$
;    e_Dust_Dispersion,$
;    e_Dust_Dispersion_inf,$
;    e_Dust_Dispersion_sup,$
;    e_Flocculence,$
;    e_Flocculence_inf,$
;    e_Flocculence_sup,$
;    e_Hot_Spots,$
;    e_Hot_Spots_inf,$
;    e_Hot_Spots_sup,$
;    e_Inclination,$
;    e_Inclination_inf,$
;    e_Inclination_sup,$
;    e_Contamination,$
;    e_Contamination_inf,$
;    e_Contamination_sup,$
;    e_Multiplicity,$
;    e_Multiplicity_inf,$
;    e_Multiplicity_sup,$

restore, '~/Astronomy/Research/GalaxyZoo/EFIGI/efigi_all.sav'

file='~/Astronomy/Research/GalaxyZoo/csv/contamination_efigi_gz2_willettk.csv'
readcol, file, $
	c_objid,c_pgc,$
	gz2_artifact_wfraction,$
	gz2_artifact_count, $
	format='a,f,f,i', $
	skip=1,/silent

match, e_pgc, c_pgc, eind, cind

!p.multi=[0,1,1]
;cgplot, e_contamination[eind], gz2_artifact_wfraction[cind], $
;	psym=5, $
;	charsize=cs, $
;	xtitle='EFIGI contamination', $
;	ytitle='GZ2 artifact weighted fraction'
;
;cgplot, indgen(1), $
;	charsize=cs, $
;	xr=[-0.125,1.125],/xstyle, $
;	yr=[0,1],/ystyle, $
;	xtickv = [0,0.25,0.50,0.75,1.0], $
;	xticks=4, $
;	position=[0.07,0.1,0.42,0.9], $
;	title='All galaxies', $
;	xtitle='EFIGI contamination', $
;	ytitle='GZ2 artifact weighted fraction'
;
;efbin = 0.25
;gz2bin = 0.1
;
;gz2_artifact = gz2_artifact_wfraction[cind]
;efigi_artifact = e_contamination[eind]
;
;contcell=fltarr(5,10)
;ccount_arr = fltarr(5,10)
;for i=0,4 do begin
;	for j=0,9 do begin
;		efval = efbin*i
;		gz2val = gz2bin * j
;		contcell = where((gz2_artifact gt gz2val) $
;			and (gz2_artifact le gz2val+gz2bin) $
;			and (efigi_artifact eq efval),ccount)
;		ccount_arr[i,j] = ccount
;	endfor
;endfor
;
;cgloadct, 3,/reverse
;
;for i=0,4 do begin
;	for j=0,9 do begin
;		barcolor=fix(float(ccount_arr[i,j])/max(ccount_arr)*255.)
;		efval = efbin*i
;		gz2val = gz2bin * j
;		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
;	endfor
;endfor
;
;cgplots, !x.crange, !y.crange, linestyle=2, thick=th
;
;cgcolorbar, /vertical, /right, position=[0.1,0.43,0.9,0.45],range=[0,max(ccount_arr)]



; Arm properties

!p.multi=[0,2,3]

efigi = mrdfits('~/Astronomy/Research/GalaxyZoo/fits/efigi_gz2_willettk.fit',1)

if n_elements(gz2armcount) eq 0 then gz2armcount = 0

; Axial cut match for arm curvature

match, objid[where(alog10(1./expab_r) lt 0.3,axcount)], efigi[where(efigi.arm_curvature ge 0 and efigi.t10_arms_winding_total_count ge gz2armcount)].objid, junkind1, axial_armcurve_ind

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/datapaper/figures/efigi_arm_curvature'+counttitle+'.eps',$
		/color, /quiet, /encap, xs=8, ys=7
	cs=2
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

; Tight wind 

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.07,0.66,0.42,0.93], $
	title='All galaxies', $
;	xtitle='EFIGI arm curvature', $
	ytitle='GZ2 tight fraction'

efbin = 0.25
gz2bin = 0.1

goodarmind = where(efigi.arm_curvature ge 0 and efigi.t10_arms_winding_total_count ge gz2armcount)

barcell=fltarr(5,10)
bccount_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((efigi.t10_arms_winding_a28_tight_weighted_fraction gt gz2val) $
			and (efigi.t10_arms_winding_a28_tight_weighted_fraction le gz2val+gz2bin) $
			and (efigi.t10_arms_winding_total_count ge gz2armcount) $
			and ((efigi.arm_curvature_sup - efigi.arm_curvature_inf) ne 1) $
			and (efigi.arm_curvature eq efval),bccount)
		bccount_arr[i,j] = bccount
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_arr[i,j])/max(bccount_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, !x.crange, !y.crange, linestyle=2, thick=th

cgcolorbar, /vertical, /right, position=[0.66,0.43,0.93,0.45],range=[0,max(bccount_arr)]

corr_na_gz = correlate(efigi[goodarmind].t10_arms_winding_a28_tight_weighted_fraction,efigi[goodarmind].arm_curvature)
print,''
print,'Correlation (tight spirals, all galaxies): ',corr_na_gz,n_elements(efigi[goodarmind].arm_curvature)

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.55,0.66,0.90,0.93], $
	title='Disks, not edge-on', $
;	xtitle='EFIGI arm curvature', $
	ytitle='GZ2 tight fraction'

efbin = 0.25
gz2bin = 0.1

barcell=fltarr(5,10)
bccount_ax_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((efigi.t10_arms_winding_a28_tight_weighted_fraction gt gz2val) $
			and (efigi.t10_arms_winding_a28_tight_weighted_fraction le gz2val+gz2bin) $
			and ((efigi.arm_curvature_sup - efigi.arm_curvature_inf) ne 1) $
			and (efigi.t10_arms_winding_total_count ge gz2armcount) $
			and (efigi.arm_curvature eq efval)$
			and alog10(1./expab_r) lt 0.3, $
			bccount_ax)
		bccount_ax_arr[i,j] = bccount_ax
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_ax_arr[i,j])/max(bccount_ax_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, !x.crange, !y.crange, linestyle=2, thick=th

cgcolorbar, /vertical, /right, position=[0.66,0.93,0.93,0.95],range=[0,max(bccount_ax_arr)]

corr_na_gz = correlate(efigi[axial_armcurve_ind].t10_arms_winding_a28_tight_weighted_fraction,efigi[axial_armcurve_ind].arm_curvature)
print,''
print,'Correlation (tight spirals, face-on galaxies): ',corr_na_gz,n_elements(axial_armcurve_ind)


	; Weighted mean of the sample

	efigi_weight = ((efigi_bar_sup - efigi_bar_inf) + 0.25)^(-2)
	efigi_err = (efigi_bar_sup - efigi_bar_inf) + 0.25
	;print, wmean(efigi.arm_curvature,efigi_err)
	;print, wmean(efigi[where(efigi.t10_arms_winding_a28_tight_weighted_fraction ge 0.5 and alog10(1./expab_r) lt 0.3)].arm_curvature,efigi_err[where(efigi.t10_arms_winding_a28_tight_weighted_fraction ge 0.5 and alog10(1./expab_r) lt 0.3)])

; Medium wind

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.07,0.35,0.42,0.60], $
;	xtitle='EFIGI arm curvature', $
	ytitle='GZ2 medium fraction'

efbin = 0.25
gz2bin = 0.1

barcell=fltarr(5,10)
bccount_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((efigi.t10_arms_winding_a29_medium_weighted_fraction gt gz2val) $
			and (efigi.t10_arms_winding_a29_medium_weighted_fraction le gz2val+gz2bin) $
			and (efigi.t10_arms_winding_total_count ge gz2armcount) $
			and ((efigi.arm_curvature_sup - efigi.arm_curvature_inf) ne 1) $
			and (efigi.arm_curvature eq efval),bccount)
		bccount_arr[i,j] = bccount
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_arr[i,j])/max(bccount_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, !x.crange, [0.5,0.5], linestyle=2, thick=th
xarr=fillarr(0.01,-2,2)
;cgoplot, xarr, exp(-1 * (xarr-0.5)^2/(2*(0.25)^2)),linestyle=2,thick=th,color='black'

cgcolorbar, /vertical, /right, position=[0.35,0.43,0.60,0.45],range=[0,max(bccount_arr)]

corr_na_gz = correlate(efigi[goodarmind].t10_arms_winding_a29_medium_weighted_fraction,efigi[goodarmind].arm_curvature)
print,''
print,'Correlation (medium spirals, all galaxies): ',corr_na_gz,n_elements(efigi[goodarmind].arm_curvature)

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.55,0.35,0.90,0.60], $
;	xtitle='EFIGI arm curvature', $
	ytitle='GZ2 medium fraction'

efbin = 0.25
gz2bin = 0.1

barcell=fltarr(5,10)
bccount_ax_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((efigi.t10_arms_winding_a29_medium_weighted_fraction gt gz2val) $
			and (efigi.t10_arms_winding_a29_medium_weighted_fraction le gz2val+gz2bin) $
			and ((efigi.arm_curvature_sup - efigi.arm_curvature_inf) ne 1) $
			and (efigi.t10_arms_winding_total_count ge gz2armcount) $
			and (efigi.arm_curvature eq efval)$
			and alog10(1./expab_r) lt 0.3, $
			bccount_ax)
		bccount_ax_arr[i,j] = bccount_ax
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_ax_arr[i,j])/max(bccount_ax_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, !x.crange, [0.5,0.5], linestyle=2, thick=th
;cgoplot, xarr, exp(-1 * (xarr-0.5)^2/(2*(0.25)^2)),linestyle=2,thick=th,color='black'

cgcolorbar, /vertical, /right, position=[0.35,0.93,0.60,0.95],range=[0,max(bccount_ax_arr)]

corr_na_gz = correlate(efigi[axial_armcurve_ind].t10_arms_winding_a29_medium_weighted_fraction,efigi[axial_armcurve_ind].arm_curvature)
print,''
print,'Correlation (medium spirals, face-on galaxies): ',corr_na_gz,n_elements(axial_armcurve_ind)


	; Weighted mean of the sample

	efigi_weight = ((efigi_bar_sup - efigi_bar_inf) + 0.25)^(-2)
	efigi_err = (efigi_bar_sup - efigi_bar_inf) + 0.25
	;print, wmean(efigi.arm_curvature,efigi_err)
	;print, wmean(efigi[where(efigi.t10_arms_winding_a29_medium_weighted_fraction ge 0.5 and alog10(1./expab_r) lt 0.3)].arm_curvature,efigi_err[where(efigi.t10_arms_winding_a29_medium_weighted_fraction ge 0.5 and alog10(1./expab_r) lt 0.3)])

; Loose wind

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.07,0.10,0.42,0.30], $
	xtitle='EFIGI arm curvature', $
	ytitle='GZ2 loose fraction'

efbin = 0.25
gz2bin = 0.1

barcell=fltarr(5,10)
bccount_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((efigi.t10_arms_winding_a30_loose_weighted_fraction gt gz2val) $
			and (efigi.t10_arms_winding_a30_loose_weighted_fraction le gz2val+gz2bin) $
			and (efigi.t10_arms_winding_total_count ge gz2armcount) $
			and ((efigi.arm_curvature_sup - efigi.arm_curvature_inf) ne 1) $
			and (efigi.arm_curvature eq efval),bccount)
		bccount_arr[i,j] = bccount
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_arr[i,j])/max(bccount_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, reverse(!x.crange), !y.crange, linestyle=2, thick=th

cgcolorbar, /vertical, /right, position=[0.10,0.43,0.30,0.45],range=[0,max(bccount_arr)]

corr_na_gz = correlate(efigi[goodarmind].t10_arms_winding_a30_loose_weighted_fraction,efigi[goodarmind].arm_curvature)
print,''
print,'Correlation (loose spirals, all galaxies): ',corr_na_gz,n_elements(efigi[goodarmind].arm_curvature)

cgplot, indgen(1), $
	charsize=cs, $
	xr=[-0.125,1.125],/xstyle, $
	yr=[0,1],/ystyle, $
	xtickv = [0,0.25,0.50,0.75,1.0], $
	xticks=4, $
	position=[0.55,0.10,0.90,0.30], $
	xtitle='EFIGI arm curvature', $
	ytitle='GZ2 loose fraction'

efbin = 0.25
gz2bin = 0.1

barcell=fltarr(5,10)
bccount_ax_arr = fltarr(5,10)
for i=0,4 do begin
	for j=0,9 do begin
		efval = efbin*i
		gz2val = gz2bin * j
		barcell = where((efigi.t10_arms_winding_a30_loose_weighted_fraction gt gz2val) $
			and (efigi.t10_arms_winding_a30_loose_weighted_fraction le gz2val+gz2bin) $
			and (efigi.t10_arms_winding_total_count ge gz2armcount) $
			and ((efigi.arm_curvature_sup - efigi.arm_curvature_inf) ne 1) $
			and (efigi.arm_curvature eq efval)$
			and alog10(1./expab_r) lt 0.3, $
			bccount_ax)
		bccount_ax_arr[i,j] = bccount_ax
	endfor
endfor

cgloadct, 3,/reverse

for i=0,4 do begin
	for j=0,9 do begin
		barcolor=fix(float(bccount_ax_arr[i,j])/max(bccount_ax_arr)*255.)
		efval = efbin*i
		gz2val = gz2bin * j
		cgcolorfill, [efval,efval,efval+efbin,efval+efbin,efval]-efbin/2., [gz2val,gz2val+gz2bin,gz2val+gz2bin,gz2val,gz2val],color=barcolor, /data
	endfor
endfor

;cgplots, reverse(!x.crange), !y.crange, linestyle=2, thick=th

cgcolorbar, /vertical, /right, position=[0.10,0.93,0.30,0.95],range=[0,max(bccount_ax_arr)]

corr_na_gz = correlate(efigi[axial_armcurve_ind].t10_arms_winding_a30_loose_weighted_fraction,efigi[axial_armcurve_ind].arm_curvature)
print,''
print,'Correlation (loose spirals, face-on galaxies): ',corr_na_gz,n_elements(axial_armcurve_ind)


	; Weighted mean of the sample

	efigi_weight = ((efigi_bar_sup - efigi_bar_inf) + 0.25)^(-2)
	efigi_err = (efigi_bar_sup - efigi_bar_inf) + 0.25
	;print, wmean(efigi.arm_curvature,efigi_err)
	;print, wmean(efigi[where(efigi.t10_arms_winding_a30_loose_weighted_fraction ge 0.5 and alog10(1./expab_r) lt 0.3)].arm_curvature,efigi_err[where(efigi.t10_arms_winding_a30_loose_weighted_fraction ge 0.5 and alog10(1./expab_r) lt 0.3)])

if keyword_set(ps) then ps_end


if keyword_set(stop) then stop

end
