;+
; NAME:
;       
;	GZ_BULGE_ROWS
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
;       Written by K. Willett                Jun 12
;-

pro gz2_bulge_rows, ps=ps, stop=stop

fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
csvdir = '~/Astronomy/Research/GalaxyZoo/csv/'
figdir = '~/Astronomy/Research/GalaxyZoo/datapaper/figures/'

; NA10 data

nagz2=mrdfits(fitsdir+'na_gz2_willettk.fit',1,/silent)

na_gz_features_or_disk_wfraction = nagz2.T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION
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

efigi=mrdfits(fitsdir+'efigi_gz2_willettk.fit',1,/silent)

efigi_gz_tight_wfraction = efigi.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION
efigi_gz_medium_wfraction = efigi.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION
efigi_gz_loose_wfraction = efigi.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION
efigi_gz_no_bulge_wfraction = efigi.T05_BULGE_PROMINENCE_A10_NO_BULGE_WEIGHTED_FRACTION
efigi_gz_just_noticeable_wfraction = efigi.T05_BULGE_PROMINENCE_A11_JUST_NOTICEABLE_WEIGHTED_FRACTION
efigi_gz_obvious_wfraction = efigi.T05_BULGE_PROMINENCE_A12_OBVIOUS_WEIGHTED_FRACTION
efigi_gz_dominant_wfraction = efigi.T05_BULGE_PROMINENCE_A13_DOMINANT_WEIGHTED_FRACTION
efigi_ttype = efigi.t
efigi_totalbc = efigi.T05_BULGE_PROMINENCE_TOTAL_COUNT
efigi_totalsc = efigi.T10_ARMS_WINDING_TOTAL_COUNT

bulgecount = 10
spiralcount = 10

; Colorful histogram for spiral tightness

!p.multi=[0,2,3]

if keyword_set(ps) then begin
	ps_start, filename=figdir+'spiraltightness_color.ps',/color, /quiet
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

bs = 0.1

; NA tight spirals

tight_countlim = na_gz_tight_wfraction[where(na_totalsc ge spiralcount)]
ttype_countlim = na_ttype[where(na_totalsc ge spiralcount)]
print,'NA10 | GZ2 spiral galaxies', n_elements(ttype_countlim)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='tight spiral weighted fraction (GZ2)', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	junk =  where(tight_countlim gt x0 and tight_countlim le x1,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; EFIGI tight spirals

tight_countlim_efigi = efigi_gz_tight_wfraction[where(efigi_totalsc ge spiralcount)]
ttype_countlim_efigi = efigi_ttype[where(efigi_totalsc ge spiralcount)]
print,'EFIGI | GZ2 spiral galaxies', n_elements(ttype_countlim_efigi)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='tight spiral weighted fraction (GZ2)', $
	ytitle='EFIGI fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_countlim_efigi ge 7 and ttype_countlim_efigi le 8,csd)
	z_sc = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_countlim_efigi ge 5 and ttype_countlim_efigi le 6,csc)
	z_sb = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_countlim_efigi ge 3 and ttype_countlim_efigi le 4,csb)
	z_sa = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_countlim_efigi ge 1 and ttype_countlim_efigi le 2,csa)
	z_s0 = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_countlim_efigi ge -3 and ttype_countlim_efigi le 0,cs0)
	z_es = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_countlim_efigi ge -6 and ttype_countlim_efigi le -4,ces)

	junk =  where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_countlim_efigi le 8,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; NA medium spirals

medium_countlim = na_gz_medium_wfraction[where(na_totalsc ge spiralcount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='medium spiral weighted fraction (GZ2)', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	junk =  where(medium_countlim gt x0 and medium_countlim le x1,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; EFIGI medium spirals

medium_countlim_efigi = efigi_gz_medium_wfraction[where(efigi_totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='medium spiral weighted fraction (GZ2)', $
	ytitle='EFIGI fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_countlim_efigi ge 7 and ttype_countlim_efigi le 8,csd)
	z_sc = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_countlim_efigi ge 5 and ttype_countlim_efigi le 6,csc)
	z_sb = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_countlim_efigi ge 3 and ttype_countlim_efigi le 4,csb)
	z_sa = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_countlim_efigi ge 1 and ttype_countlim_efigi le 2,csa)
	z_s0 = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_countlim_efigi ge -3 and ttype_countlim_efigi le 0,cs0)
	z_es = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_countlim_efigi ge -6 and ttype_countlim_efigi le -4,ces)

	junk =  where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_countlim_efigi le 8,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; NA loose spirals

loose_countlim = na_gz_loose_wfraction[where(na_totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='loose spiral weighted fraction (GZ2)', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	junk =  where(loose_countlim gt x0 and loose_countlim le x1,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; EFIGI loose spirals

loose_countlim_efigi = efigi_gz_loose_wfraction[where(efigi_totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='loose spiral weighted fraction (GZ2)', $
	ytitle='EFIGI fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_countlim_efigi ge 7 and ttype_countlim_efigi le 8,csd)
	z_sc = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_countlim_efigi ge 5 and ttype_countlim_efigi le 6,csc)
	z_sb = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_countlim_efigi ge 3 and ttype_countlim_efigi le 4,csb)
	z_sa = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_countlim_efigi ge 1 and ttype_countlim_efigi le 2,csa)
	z_s0 = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_countlim_efigi ge -3 and ttype_countlim_efigi le 0,cs0)
	z_es = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_countlim_efigi ge -6 and ttype_countlim_efigi le -4,ces)

	junk =  where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_countlim_efigi le 8,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

if keyword_set(ps) then ps_end



; Colorful histogram for bulge prominence

!p.multi=[0,2,4]

if keyword_set(ps) then begin
	ps_start, filename=figdir+'bulgeprominence_color.ps',/color, /quiet
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

sd = where(na_ttype ge 7  and na_ttype le 10)
sc = where(na_ttype ge 5  and na_ttype le 6)
sb = where(na_ttype ge 3  and na_ttype le 4)
sa = where(na_ttype ge 1  and na_ttype le 2)
s0 = where(na_ttype ge -3 and na_ttype le 0)
es = where(na_ttype eq -5)

bs = 0.1

; No bulges

no_bulge_countlim = na_gz_no_bulge_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]
print,'NA10 | GZ2 bulge galaxies', n_elements(ttype_bc_countlim)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='no_bulge weighted fraction (GZ2)', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(no_bulge_countlim gt x0 and no_bulge_countlim le x1 and ttype_bc_countlim ge 7,csd)
	z_sc = where(no_bulge_countlim gt x0 and no_bulge_countlim le x1 and ttype_bc_countlim ge 5 and ttype_bc_countlim le 6,csc)
	z_sb = where(no_bulge_countlim gt x0 and no_bulge_countlim le x1 and ttype_bc_countlim ge 3 and ttype_bc_countlim le 4,csb)
	z_sa = where(no_bulge_countlim gt x0 and no_bulge_countlim le x1 and ttype_bc_countlim ge 1 and ttype_bc_countlim le 2,csa)
	z_s0 = where(no_bulge_countlim gt x0 and no_bulge_countlim le x1 and ttype_bc_countlim ge -3 and ttype_bc_countlim le 0,cs0)
	z_es = where(no_bulge_countlim gt x0 and no_bulge_countlim le x1 and ttype_bc_countlim ge -5 and ttype_bc_countlim le -4,ces)

	junk =  where(no_bulge_countlim gt x0 and no_bulge_countlim le x1,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; EFIGI No bulge

no_bulge_countlim_efigi = efigi_gz_no_bulge_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]
print,'EFIGI | GZ2 bulge galaxies', n_elements(ttype_bc_countlim_efigi)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='no_bulge weighted fraction (GZ2)', $
	ytitle='EFIGI fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(no_bulge_countlim_efigi gt x0 and no_bulge_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 7 and ttype_bc_countlim_efigi le 8,csd)
	z_sc = where(no_bulge_countlim_efigi gt x0 and no_bulge_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 5 and ttype_bc_countlim_efigi le 6,csc)
	z_sb = where(no_bulge_countlim_efigi gt x0 and no_bulge_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 3 and ttype_bc_countlim_efigi le 4,csb)
	z_sa = where(no_bulge_countlim_efigi gt x0 and no_bulge_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 1 and ttype_bc_countlim_efigi le 2,csa)
	z_s0 = where(no_bulge_countlim_efigi gt x0 and no_bulge_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -3 and ttype_bc_countlim_efigi le 0,cs0)
	z_es = where(no_bulge_countlim_efigi gt x0 and no_bulge_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -6 and ttype_bc_countlim_efigi le -4,ces)

	junk =  where(no_bulge_countlim_efigi gt x0 and no_bulge_countlim_efigi le x1 and ttype_bc_countlim_efigi le 8,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Just noticeable bulges

just_noticeable_countlim = na_gz_just_noticeable_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='just_noticeable weighted fraction (GZ2)', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(just_noticeable_countlim gt x0 and just_noticeable_countlim le x1 and ttype_bc_countlim ge 7,csd)
	z_sc = where(just_noticeable_countlim gt x0 and just_noticeable_countlim le x1 and ttype_bc_countlim ge 5 and ttype_bc_countlim le 6,csc)
	z_sb = where(just_noticeable_countlim gt x0 and just_noticeable_countlim le x1 and ttype_bc_countlim ge 3 and ttype_bc_countlim le 4,csb)
	z_sa = where(just_noticeable_countlim gt x0 and just_noticeable_countlim le x1 and ttype_bc_countlim ge 1 and ttype_bc_countlim le 2,csa)
	z_s0 = where(just_noticeable_countlim gt x0 and just_noticeable_countlim le x1 and ttype_bc_countlim ge -3 and ttype_bc_countlim le 0,cs0)
	z_es = where(just_noticeable_countlim gt x0 and just_noticeable_countlim le x1 and ttype_bc_countlim ge -5 and ttype_bc_countlim le -4,ces)

	junk =  where(just_noticeable_countlim gt x0 and just_noticeable_countlim le x1,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Just noticeable

just_noticeable_countlim_efigi = efigi_gz_just_noticeable_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='just_noticeable weighted fraction (GZ2)', $
	ytitle='EFIGI fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(just_noticeable_countlim_efigi gt x0 and just_noticeable_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 7 and ttype_bc_countlim_efigi le 8,csd)
	z_sc = where(just_noticeable_countlim_efigi gt x0 and just_noticeable_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 5 and ttype_bc_countlim_efigi le 6,csc)
	z_sb = where(just_noticeable_countlim_efigi gt x0 and just_noticeable_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 3 and ttype_bc_countlim_efigi le 4,csb)
	z_sa = where(just_noticeable_countlim_efigi gt x0 and just_noticeable_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 1 and ttype_bc_countlim_efigi le 2,csa)
	z_s0 = where(just_noticeable_countlim_efigi gt x0 and just_noticeable_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -3 and ttype_bc_countlim_efigi le 0,cs0)
	z_es = where(just_noticeable_countlim_efigi gt x0 and just_noticeable_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -6 and ttype_bc_countlim_efigi le -4,ces)

	junk =  where(just_noticeable_countlim_efigi gt x0 and just_noticeable_countlim_efigi le x1 and ttype_bc_countlim_efigi le 8,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Obvious bulges

obvious_countlim = na_gz_obvious_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='obvious weighted fraction (GZ2)', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(obvious_countlim gt x0 and obvious_countlim le x1 and ttype_bc_countlim ge 7,csd)
	z_sc = where(obvious_countlim gt x0 and obvious_countlim le x1 and ttype_bc_countlim ge 5 and ttype_bc_countlim le 6,csc)
	z_sb = where(obvious_countlim gt x0 and obvious_countlim le x1 and ttype_bc_countlim ge 3 and ttype_bc_countlim le 4,csb)
	z_sa = where(obvious_countlim gt x0 and obvious_countlim le x1 and ttype_bc_countlim ge 1 and ttype_bc_countlim le 2,csa)
	z_s0 = where(obvious_countlim gt x0 and obvious_countlim le x1 and ttype_bc_countlim ge -3 and ttype_bc_countlim le 0,cs0)
	z_es = where(obvious_countlim gt x0 and obvious_countlim le x1 and ttype_bc_countlim ge -5 and ttype_bc_countlim le -4,ces)

	junk =  where(obvious_countlim gt x0 and obvious_countlim le x1,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Obvious

obvious_countlim_efigi = efigi_gz_obvious_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='obvious weighted fraction (GZ2)', $
	ytitle='EFIGI fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(obvious_countlim_efigi gt x0 and obvious_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 7 and ttype_bc_countlim_efigi le 8,csd)
	z_sc = where(obvious_countlim_efigi gt x0 and obvious_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 5 and ttype_bc_countlim_efigi le 6,csc)
	z_sb = where(obvious_countlim_efigi gt x0 and obvious_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 3 and ttype_bc_countlim_efigi le 4,csb)
	z_sa = where(obvious_countlim_efigi gt x0 and obvious_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 1 and ttype_bc_countlim_efigi le 2,csa)
	z_s0 = where(obvious_countlim_efigi gt x0 and obvious_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -3 and ttype_bc_countlim_efigi le 0,cs0)
	z_es = where(obvious_countlim_efigi gt x0 and obvious_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -6 and ttype_bc_countlim_efigi le -4,ces)

	junk =  where(obvious_countlim_efigi gt x0 and obvious_countlim_efigi le x1 and ttype_bc_countlim_efigi le 8,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Dominant bulge

dominant_countlim = na_gz_dominant_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='dominant weighted fraction (GZ2)', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(dominant_countlim gt x0 and dominant_countlim le x1 and ttype_bc_countlim ge 7,csd)
	z_sc = where(dominant_countlim gt x0 and dominant_countlim le x1 and ttype_bc_countlim ge 5 and ttype_bc_countlim le 6,csc)
	z_sb = where(dominant_countlim gt x0 and dominant_countlim le x1 and ttype_bc_countlim ge 3 and ttype_bc_countlim le 4,csb)
	z_sa = where(dominant_countlim gt x0 and dominant_countlim le x1 and ttype_bc_countlim ge 1 and ttype_bc_countlim le 2,csa)
	z_s0 = where(dominant_countlim gt x0 and dominant_countlim le x1 and ttype_bc_countlim ge -3 and ttype_bc_countlim le 0,cs0)
	z_es = where(dominant_countlim gt x0 and dominant_countlim le x1 and ttype_bc_countlim ge -5 and ttype_bc_countlim le -4,ces)

	junk =  where(dominant_countlim gt x0 and dominant_countlim le x1,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Dominant

dominant_countlim_efigi = efigi_gz_dominant_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='dominant weighted fraction (GZ2)', $
	ytitle='EFIGI fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(dominant_countlim_efigi gt x0 and dominant_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 7 and ttype_bc_countlim_efigi le 8,csd)
	z_sc = where(dominant_countlim_efigi gt x0 and dominant_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 5 and ttype_bc_countlim_efigi le 6,csc)
	z_sb = where(dominant_countlim_efigi gt x0 and dominant_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 3 and ttype_bc_countlim_efigi le 4,csb)
	z_sa = where(dominant_countlim_efigi gt x0 and dominant_countlim_efigi le x1 and ttype_bc_countlim_efigi ge 1 and ttype_bc_countlim_efigi le 2,csa)
	z_s0 = where(dominant_countlim_efigi gt x0 and dominant_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -3 and ttype_bc_countlim_efigi le 0,cs0)
	z_es = where(dominant_countlim_efigi gt x0 and dominant_countlim_efigi le x1 and ttype_bc_countlim_efigi ge -6 and ttype_bc_countlim_efigi le -4,ces)

	junk =  where(dominant_countlim_efigi gt x0 and dominant_countlim_efigi le x1 and ttype_bc_countlim_efigi le 8,c_all)

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

if keyword_set(ps) then ps_end





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Weighted bulge dominance vs. weighted spiral tightness, color-coded by morphology
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





efigi_wbulge = (0. * efigi_gz_no_bulge_wfraction + 0.33 * efigi_gz_just_noticeable_wfraction + 0.66 * efigi_gz_obvious_wfraction + 1.0 * efigi_gz_dominant_wfraction)
efigi_wspiral = (0. * efigi_gz_tight_wfraction + 0.5 * efigi_gz_medium_wfraction + 1.0 * efigi_gz_loose_wfraction)

na_wbulge = (0. * na_gz_no_bulge_wfraction + 0.33 * na_gz_just_noticeable_wfraction + 0.66 * na_gz_obvious_wfraction + 1.0 * na_gz_dominant_wfraction)
na_wspiral = (0. * na_gz_tight_wfraction + 0.5 * na_gz_medium_wfraction + 1.0 * na_gz_loose_wfraction)

 
if keyword_set(ps) then begin
	ps_start, filename=figdir+'wbulge_wspiral.eps',/color, /quiet, /encap, xs=7, ys=10
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

!p.multi=[0,1,2]

ss = 0.8

; NA10

cgplot, na_wbulge, na_wspiral, $
	xtitle='GZ2 combined bulge prominence', $
	ytitle='GZ2 combined spiral prominence', $
	title='NA10', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.1,1.1], /ystyle, $
	/nodata

cgplot, na_wbulge[where(na_ttype ge 7 and na_ttype le 8 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $ 
	na_wspiral[where(na_ttype ge 7 and na_ttype le 8 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $
	/overplot, $
	psym = 15, $
	symsize=ss, $
	color='blue'

cgplot, na_wbulge[where(na_ttype ge 5 and na_ttype le 6 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $ 
	na_wspiral[where(na_ttype ge 5 and na_ttype le 6 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $
	/overplot, $
	psym = 15, $
	symsize=ss, $
	color='dark green'

cgplot, na_wbulge[where(na_ttype ge 3 and na_ttype le 4 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $ 
	na_wspiral[where(na_ttype ge 3 and na_ttype le 4 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $
	/overplot, $
	psym = 15, $
	symsize=ss, $
	color='green'

cgplot, na_wbulge[where(na_ttype ge 1 and na_ttype le 2 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $ 
	na_wspiral[where(na_ttype ge 1 and na_ttype le 2 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $
	/overplot, $
	psym = 15, $
	symsize=ss, $
	color='yellow'

cgplot, na_wbulge[where(na_ttype ge -3 and na_ttype le 0 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $ 
	na_wspiral[where(na_ttype ge -3 and na_ttype le 0 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $
	/overplot, $
	psym = 15, $
	symsize=ss, $
	color='tomato'

cgplot, na_wbulge[where(na_ttype ge -6 and na_ttype le -4 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $ 
	na_wspiral[where(na_ttype ge -6 and na_ttype le -4 and na_totalbc ge bulgecount and na_totalsc ge spiralcount)], $
	/overplot, $
	psym = 15, $
	symsize=ss, $
	color='black'

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right


cgplot, efigi_wbulge, efigi_wspiral, $
	xtitle='GZ2 combined bulge prominence', $
	ytitle='GZ2 combined spiral prominence', $
	title='EFIGI', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.1,1.1], /ystyle, $
	/nodata

cgplot, efigi_wbulge[where(efigi_ttype ge 7 and efigi_ttype le 8 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $ 
	efigi_wspiral[where(efigi_ttype ge 7 and efigi_ttype le 8 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $
	/overplot, $
	psym = 16, $
	symsize=ss, $
	color='blue'

cgplot, efigi_wbulge[where(efigi_ttype ge 5 and efigi_ttype le 6 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $ 
	efigi_wspiral[where(efigi_ttype ge 5 and efigi_ttype le 6 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $
	/overplot, $
	psym = 16, $
	symsize=ss, $
	color='dark green'

cgplot, efigi_wbulge[where(efigi_ttype ge 3 and efigi_ttype le 4 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $ 
	efigi_wspiral[where(efigi_ttype ge 3 and efigi_ttype le 4 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $
	/overplot, $
	psym = 16, $
	symsize=ss, $
	color='green'

cgplot, efigi_wbulge[where(efigi_ttype ge 1 and efigi_ttype le 2 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $ 
	efigi_wspiral[where(efigi_ttype ge 1 and efigi_ttype le 2 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $
	/overplot, $
	psym = 16, $
	symsize=ss, $
	color='yellow'

cgplot, efigi_wbulge[where(efigi_ttype ge -3 and efigi_ttype le 0 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $ 
	efigi_wspiral[where(efigi_ttype ge -3 and efigi_ttype le 0 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $
	/overplot, $
	psym = 16, $
	symsize=ss, $
	color='tomato'

cgplot, efigi_wbulge[where(efigi_ttype ge -6 and efigi_ttype le -4 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $ 
	efigi_wspiral[where(efigi_ttype ge -6 and efigi_ttype le -4 and efigi_totalbc ge bulgecount and efigi_totalsc ge spiralcount)], $
	/overplot, $
	psym = 16, $
	symsize=ss, $
	color='black'

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right


if keyword_set(ps) then ps_end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot spiral tightness vs. pitch angle
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


gz2 = mrdfits(fitsdir+'gz2table.fits',1,/silent)
pa_gz2 = mrdfits(fitsdir+'pitchangle_gz2.fits',1,/silent)

; need to match gz2 structure against the pitch angle catalog; use SQL or TOPCAT

gz_tight_wfraction = pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION
gz_medium_wfraction = pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION
gz_loose_wfraction = pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION
gz_no_bulge_wfraction = pa_gz2.T05_BULGE_PROMINENCE_A10_NO_BULGE_WEIGHTED_FRACTION
gz_just_noticeable_wfraction = pa_gz2.T05_BULGE_PROMINENCE_A11_JUST_NOTICEABLE_WEIGHTED_FRACTION
gz_obvious_wfraction = pa_gz2.T05_BULGE_PROMINENCE_A12_OBVIOUS_WEIGHTED_FRACTION
gz_dominant_wfraction = pa_gz2.T05_BULGE_PROMINENCE_A13_DOMINANT_WEIGHTED_FRACTION

gz2_spiral_ind = where(pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount)

wbulge = (0. * gz_no_bulge_wfraction + 0.33 * gz_just_noticeable_wfraction + 0.66 * gz_obvious_wfraction + 1.0 * gz_dominant_wfraction)
wspiral = (0. * gz_tight_wfraction + 0.5 * gz_medium_wfraction + 1.0 * gz_loose_wfraction)

;pitchanglefile = csvdir+'pitch_angle.csv'
;readcol, pitchanglefile, $
;	name, fit_state, chirality_maj, chirality_alenWtd, chirality_wtdPangSum, chirality_longestArc, axisRatio, minAxsLen, majAxsLen, majAxsAngle, contourBrtRatio, chirality_votes_maj, chirality_votes_alenWtd, alenWtdPangSum, chirality_longest_arc, top2_chirality_agreement, pa_longest, pa_avg, pa_avg_abs, pa_avg_domChiralityOnly, pa_alenWtd_avg, pa_alenWtd_avg_abs, pa_alenWtd_avg_domChiralityOnly, pa_totBrtWtd, pa_avgBrtWtd, $
;	delimiter=',', $
;	skipline=1, $
;	/quiet, $
; 	format='a, a, a, a, a, a, f, f, f, f, f, a, a, f, i, a, f, f, f, f, f, f, f, f, f'
;
;pitchangle_fitsname= fitsdir+'pitch_angle.fits'
;
;pa_sing = {name:name[0], fit_state:fit_state[0], chirality_maj:chirality_maj[0], chirality_alenWtd:chirality_alenWtd[0], chirality_wtdPangSum:chirality_wtdPangSum[0], chirality_longestArc:chirality_longestArc[0], axisRatio:axisRatio[0], minAxsLen:minAxsLen[0], majAxsLen:majAxsLen[0], majAxsAngle:majAxsAngle[0], contourBrtRatio:contourBrtRatio[0], chirality_votes_maj:chirality_votes_maj[0], chirality_votes_alenWtd:chirality_votes_alenWtd[0], alenWtdPangSum:alenWtdPangSum[0], chirality_longest_arc:chirality_longest_arc[0], top2_chirality_agreement:top2_chirality_agreement[0], pa_longest:pa_longest[0], pa_avg:pa_avg[0], pa_avg_abs:pa_avg_abs[0], pa_avg_domChiralityOnly:pa_avg_domChiralityOnly[0], pa_alenWtd_avg:pa_alenWtd_avg[0], pa_alenWtd_avg_abs:pa_alenWtd_avg_abs[0], pa_alenWtd_avg_domChiralityOnly:pa_alenWtd_avg_domChiralityOnly[0], pa_totBrtWtd:pa_totBrtWtd[0], pa_avgBrtWtd:pa_avgBrtWtd[0]}
;pa = replicate(pa_sing,n_elements(name))
;
; pa.(0)=name
; pa.(1)=fit_state
; pa.(2)=chirality_maj
; pa.(3)=chirality_alenWtd
; pa.(4)=chirality_wtdPangSum
; pa.(5)=chirality_longestArc
; pa.(6)=axisRatio
; pa.(7)=minAxsLen
; pa.(8)=majAxsLen
; pa.(9)=majAxsAngle
; pa.(10)=contourBrtRatio
; pa.(11)=chirality_votes_maj
; pa.(12)=chirality_votes_alenWtd
; pa.(13)=alenWtdPangSum
; pa.(14)=chirality_longest_arc
; pa.(15)=top2_chirality_agreement
; pa.(16)=pa_longest
; pa.(17)=pa_avg
; pa.(18)=pa_avg_abs
; pa.(19)=pa_avg_domChiralityOnly
; pa.(20)=pa_alenWtd_avg
; pa.(21)=pa_alenWtd_avg_abs
; pa.(22)=pa_alenWtd_avg_domChiralityOnly
; pa.(23)=pa_totBrtWtd
; pa.(24)=pa_avgBrtWtd
;
;mwrfits, pa, pitchangle_fitsname

pa = mrdfits(fitsdir+'pitch_angle.fits',1,/silent)

;!p.multi=[0,1,2]
;
;cgplot, abs(pa_alenWtd_avg_domChiralityOnly[gz2_spiral_ind]), wspiral[gz2_spiral_ind], $
;	psym = 16, $
;	xtitle='Pitch angle', $
;	ytitle='GZ2 spiral prominence'
;
;cgplot, abs(pa_alenWtd_avg_domChiralityOnly[gz2_spiral_ind]), wbulge[gz2_spiral_ind], $
;	psym = 16, $
;	xtitle='Pitch angle', $
;	ytitle='GZ2 bulge prominence'

; See whether plurality vote affects pitch angle

!p.multi=[0,1,2]

pa_gz2_tight_ind = where($
	pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount and $
	pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION and $
	pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION)

pa_gz2_medium_ind = where($
	pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount and $
	pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION and $
	pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION)

pa_gz2_loose_ind = where($
	pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount and $
	pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION and $
	pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION)


pa_binsize = 2

cghistoplot, abs(pa_gz2[pa_gz2_tight_ind].pa_alenWtd_avg_domChiralityOnly), $
	binsize = pa_binsize, $
	/outline, $
	xr=[0,50], $
	yr=[0,2000], $
	xtitle='Pitch angle [deg]', $
	datacolor='blue'

cghistoplot, abs(pa_gz2[pa_gz2_medium_ind].pa_alenWtd_avg_domChiralityOnly), $
	binsize = pa_binsize, $
	/outline, $
	/oplot, $
	datacolor='green'

cghistoplot, abs(pa_gz2[pa_gz2_loose_ind].pa_alenWtd_avg_domChiralityOnly), $
	binsize = pa_binsize, $
	/outline, $
	/oplot, $
	datacolor='red'

kstwo, abs(pa_gz2[pa_gz2_tight_ind].pa_alenWtd_avg_domChiralityOnly), abs(pa_gz2[pa_gz2_medium_ind].pa_alenWtd_avg_domChiralityOnly), d,p & print,'KS2: tight and medium', d,p
print,mean(abs(pa_gz2[pa_gz2_tight_ind].pa_alenWtd_avg_domChiralityOnly)), stddev(abs(pa_gz2[pa_gz2_tight_ind].pa_alenWtd_avg_domChiralityOnly))
kstwo, abs(pa_gz2[pa_gz2_tight_ind].pa_alenWtd_avg_domChiralityOnly), abs(pa_gz2[pa_gz2_loose_ind].pa_alenWtd_avg_domChiralityOnly), d,p & print,'KS2: tight and loose', d,p
print,mean(abs(pa_gz2[pa_gz2_medium_ind].pa_alenWtd_avg_domChiralityOnly)), stddev(abs(pa_gz2[pa_gz2_medium_ind].pa_alenWtd_avg_domChiralityOnly))
kstwo, abs(pa_gz2[pa_gz2_loose_ind].pa_alenWtd_avg_domChiralityOnly), abs(pa_gz2[pa_gz2_medium_ind].pa_alenWtd_avg_domChiralityOnly), d,p & print,'KS2: medium and loose', d,p
print,mean(abs(pa_gz2[pa_gz2_loose_ind].pa_alenWtd_avg_domChiralityOnly)), stddev(abs(pa_gz2[pa_gz2_loose_ind].pa_alenWtd_avg_domChiralityOnly))

; Reproduce Figure 2 from Davis and Hayes

anglebin = fillarr(pa_binsize,0,50)
nang = n_elements(anglebin)

tightproportion = fltarr(nang)
mediumproportion = fltarr(nang)
looseproportion = fltarr(nang)

for i=0,nang-1 do begin
	tight = where($
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) ge anglebin[i] and $
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) lt anglebin[i]+pa_binsize and $
		pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount and $
		pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION and $
		pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION, $
		ntight)

	medium = where($
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) ge anglebin[i] and $
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) lt anglebin[i]+pa_binsize and $
		pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount and $
		pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION and $
		pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION, $
		nmedium)

	tight = where($
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) ge anglebin[i] and $
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) lt anglebin[i]+pa_binsize and $
		pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount and $
		pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION and $
		pa_gz2.T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION gt pa_gz2.T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION, $
		nloose)

	all = where($
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) ge anglebin[i] and $
		abs(pa_gz2.pa_alenWtd_avg_domChiralityOnly) lt anglebin[i]+pa_binsize and $
		pa_gz2.T10_ARMS_WINDING_TOTAL_COUNT ge spiralcount, $
		nall)

	tightproportion[i] = float(ntight) / float(nall)
	mediumproportion[i] = float(nmedium) / float(nall)
	looseproportion[i] = float(nloose) / float(nall)
endfor


cgplot, anglebin, tightproportion, color='blue', xtitle='Pitch angle [deg]', ytitle='GZ2 vote proportion',yr=[0,1]
cgplot, anglebin, mediumproportion, color='green', /overplot
cgplot, anglebin, looseproportion, color='red', /overplot

if keyword_set(stop) then stop
end
