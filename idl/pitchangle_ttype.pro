
;+
; NAME:
;       
;	PITCHANGLE_TTYPE
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
;       Written by K. Willett                Jun 2012
;-


; NA10 data

nagz2=mrdfits('~/Astronomy/Research/GalaxyZoo/na_gz2_willettk.fit',1,/silent)

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

efigi=mrdfits('~/Astronomy/Research/GalaxyZoo/efigi_gz2_willettk.fit',1,/silent)

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

gz2 = mrdfits('~/Astronomy/Research/GalaxyZoo/gz2table.fits',1,/silent)
pa_gz2 = mrdfits('~/Astronomy/Research/GalaxyZoo/pitchangle_gz2.fits',1,/silent)

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

pa = mrdfits('~/Astronomy/Research/GalaxyZoo/pitch_angle.fits',1)


; Pitch angle as function of morphology (NA10)

ndeg = 50
bs = 2
nbins=ceil(float(ndeg)/bs)

na_pitchangle_countlim = na_gz_tight_wfraction[where(na_totalsc ge spiralcount)]
na_ttype_countlim = na_ttype[where(na_totalsc ge spiralcount)]
print,'NA10 | GZ2 spiral galaxies', n_elements(na_ttype_countlim)

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='Pitch angle [deg]', $
	ytitle='NA10 fraction per bin', $
	xr=[-0.1,1.3], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,nbins-1 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(na_tight_countlim gt x0 and na_tight_countlim le x1 and na_ttype_countlim ge 7,csd)
	z_sc = where(na_tight_countlim gt x0 and na_tight_countlim le x1 and na_ttype_countlim ge 5 and na_ttype_countlim le 6,csc)
	z_sb = where(na_tight_countlim gt x0 and na_tight_countlim le x1 and na_ttype_countlim ge 3 and na_ttype_countlim le 4,csb)
	z_sa = where(na_tight_countlim gt x0 and na_tight_countlim le x1 and na_ttype_countlim ge 1 and na_ttype_countlim le 2,csa)
	z_s0 = where(na_tight_countlim gt x0 and na_tight_countlim le x1 and na_ttype_countlim ge -3 and na_ttype_countlim le 0,cs0)
	z_es = where(na_tight_countlim gt x0 and na_tight_countlim le x1 and na_ttype_countlim ge -5 and na_ttype_countlim le -4,ces)

	junk =  where(na_tight_countlim gt x0 and na_tight_countlim le x1,c_all)

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



end
