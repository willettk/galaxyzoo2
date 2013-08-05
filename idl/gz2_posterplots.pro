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

pro gz2_posterplots, ps=ps, stop=stop

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

; Colorful histogram for bulge prominence

!p.multi=[0,2,4]

if keyword_set(ps) then begin
	ps_start, filename=figdir+'poster_na_efigi_bulge_ttype.eps',/color, /quiet, xs=14, ys=16, /encap
	cs=4.0
    csnum = 1.0
    cslegend = 1.3
    cthick = 2.0
	th=3
	thickline=5
	thinline=1
endif else begin
	cs=2
    csnum = 1
    cslegend = 1.0
    cthick = 1.0
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
xticks = 14
xticknames = [' ','0.0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1.0',' ',' ',' ']

; No bulges

no_bulge_countlim = na_gz_no_bulge_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]
print,'NA10 | GZ2 bulge galaxies', n_elements(ttype_bc_countlim)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
;    xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='no bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
    title='Nair & Abraham', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; EFIGI No bulge

no_bulge_countlim_efigi = efigi_gz_no_bulge_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]
print,'EFIGI | GZ2 bulge galaxies', n_elements(ttype_bc_countlim_efigi)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='no bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
    title='EFIGI', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Just noticeable bulges

just_noticeable_countlim = na_gz_just_noticeable_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='just noticeable bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Just noticeable

just_noticeable_countlim_efigi = efigi_gz_just_noticeable_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='just noticeable bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Obvious bulges

obvious_countlim = na_gz_obvious_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='obvious bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Obvious

obvious_countlim_efigi = efigi_gz_obvious_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='obvious bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Dominant bulge

dominant_countlim = na_gz_dominant_wfraction[where(na_totalbc ge bulgecount)]
ttype_bc_countlim = na_ttype[where(na_totalbc ge bulgecount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='dominant bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; Dominant

dominant_countlim_efigi = efigi_gz_dominant_wfraction[where(efigi_totalbc ge bulgecount)]
ttype_bc_countlim_efigi = efigi_ttype[where(efigi_totalbc ge bulgecount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='dominant bulge vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

if keyword_set(ps) then ps_end

; Colorful histogram for spiral arm tightness

!p.multi=[0,2,3]

if keyword_set(ps) then begin
	ps_start, filename=figdir+'poster_na_efigi_spiral_ttype.eps',/color, /quiet, xs=14, ys=16, /encap
	cs=4.0
    csnum = 1.0
    cslegend = 1.3
    cthick = 2.0
	th=3
	thickline=5
	thinline=1
endif else begin
	cs=2
    csnum = 1
    cslegend = 1.0
    cthick = 1.0
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
xticks = 14
xticknames = [' ','0.0',' ','0.2',' ','0.4',' ','0.6',' ','0.8',' ','1.0',' ',' ',' ']

; Loose

loose_countlim = na_gz_loose_wfraction[where(na_totalsc ge spiralcount)]
ttype_sc_countlim = na_ttype[where(na_totalsc ge spiralcount)]
print,'NA10 | GZ2 bulge galaxies', n_elements(ttype_sc_countlim)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
;    xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='loose spiral vote fraction (GZ2)', $
	ytitle='T-type fraction', $
    title='Nair & Abraham', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_sc_countlim ge 7,csd)
	z_sc = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_sc_countlim ge 5 and ttype_sc_countlim le 6,csc)
	z_sb = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_sc_countlim ge 3 and ttype_sc_countlim le 4,csb)
	z_sa = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_sc_countlim ge 1 and ttype_sc_countlim le 2,csa)
	z_s0 = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_sc_countlim ge -3 and ttype_sc_countlim le 0,cs0)
	z_es = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_sc_countlim ge -5 and ttype_sc_countlim le -4,ces)

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

; EFIGI loose

loose_countlim_efigi = efigi_gz_loose_wfraction[where(efigi_totalsc ge spiralcount)]
ttype_sc_countlim_efigi = efigi_ttype[where(efigi_totalsc ge spiralcount)]
print,'EFIGI | GZ2 bulge galaxies', n_elements(ttype_sc_countlim_efigi)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='loose spiral vote fraction (GZ2)', $
	ytitle='T-type fraction', $
    title='EFIGI', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 7 and ttype_sc_countlim_efigi le 8,csd)
	z_sc = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 5 and ttype_sc_countlim_efigi le 6,csc)
	z_sb = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 3 and ttype_sc_countlim_efigi le 4,csb)
	z_sa = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 1 and ttype_sc_countlim_efigi le 2,csa)
	z_s0 = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_sc_countlim_efigi ge -3 and ttype_sc_countlim_efigi le 0,cs0)
	z_es = where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_sc_countlim_efigi ge -6 and ttype_sc_countlim_efigi le -4,ces)

	junk =  where(loose_countlim_efigi gt x0 and loose_countlim_efigi le x1 and ttype_sc_countlim_efigi le 8,c_all)

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

; NA10 medium

medium_countlim = na_gz_medium_wfraction[where(na_totalsc ge spiralcount)]
ttype_sc_countlim = na_ttype[where(na_totalsc ge spiralcount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='medium spiral vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_sc_countlim ge 7,csd)
	z_sc = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_sc_countlim ge 5 and ttype_sc_countlim le 6,csc)
	z_sb = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_sc_countlim ge 3 and ttype_sc_countlim le 4,csb)
	z_sa = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_sc_countlim ge 1 and ttype_sc_countlim le 2,csa)
	z_s0 = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_sc_countlim ge -3 and ttype_sc_countlim le 0,cs0)
	z_es = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_sc_countlim ge -5 and ttype_sc_countlim le -4,ces)

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; EFIGI medium

medium_countlim_efigi = efigi_gz_medium_wfraction[where(efigi_totalsc ge spiralcount)]
ttype_sc_countlim_efigi = efigi_ttype[where(efigi_totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='medium spiral vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 7 and ttype_sc_countlim_efigi le 8,csd)
	z_sc = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 5 and ttype_sc_countlim_efigi le 6,csc)
	z_sb = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 3 and ttype_sc_countlim_efigi le 4,csb)
	z_sa = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 1 and ttype_sc_countlim_efigi le 2,csa)
	z_s0 = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_sc_countlim_efigi ge -3 and ttype_sc_countlim_efigi le 0,cs0)
	z_es = where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_sc_countlim_efigi ge -6 and ttype_sc_countlim_efigi le -4,ces)

	junk =  where(medium_countlim_efigi gt x0 and medium_countlim_efigi le x1 and ttype_sc_countlim_efigi le 8,c_all)

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; NA10 tight

tight_countlim = na_gz_tight_wfraction[where(na_totalsc ge spiralcount)]
ttype_sc_countlim = na_ttype[where(na_totalsc ge spiralcount)]


cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='tight spiral vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_sc_countlim ge 7,csd)
	z_sc = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_sc_countlim ge 5 and ttype_sc_countlim le 6,csc)
	z_sb = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_sc_countlim ge 3 and ttype_sc_countlim le 4,csb)
	z_sa = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_sc_countlim ge 1 and ttype_sc_countlim le 2,csa)
	z_s0 = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_sc_countlim ge -3 and ttype_sc_countlim le 0,cs0)
	z_es = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_sc_countlim ge -5 and ttype_sc_countlim le -4,ces)

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

; EFIGI tight

tight_countlim_efigi = efigi_gz_tight_wfraction[where(efigi_totalsc ge spiralcount)]
ttype_sc_countlim_efigi = efigi_ttype[where(efigi_totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
    ;xticks = xticks, xtickname=xticknames, $
	charsize=cs, $
    charthick = cthick, $
	xtitle='tight spiral vote fraction (GZ2)', $
	ytitle='T-type fraction', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.19,1.19],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 7 and ttype_sc_countlim_efigi le 8,csd)
	z_sc = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 5 and ttype_sc_countlim_efigi le 6,csc)
	z_sb = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 3 and ttype_sc_countlim_efigi le 4,csb)
	z_sa = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_sc_countlim_efigi ge 1 and ttype_sc_countlim_efigi le 2,csa)
	z_s0 = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_sc_countlim_efigi ge -3 and ttype_sc_countlim_efigi le 0,cs0)
	z_es = where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_sc_countlim_efigi ge -6 and ttype_sc_countlim_efigi le -4,ces)

	junk =  where(tight_countlim_efigi gt x0 and tight_countlim_efigi le x1 and ttype_sc_countlim_efigi le 8,c_all)

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=csnum


endfor

;al_legend, charsize=cslegend, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

if keyword_set(ps) then ps_end





if keyword_set(stop) then stop
end

