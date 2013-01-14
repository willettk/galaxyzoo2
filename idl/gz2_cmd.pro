
;+
; NAME:
;       
;	GZ2_CMD
;
; PURPOSE:
;
;	Plot the morphological demographics of the CMD
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
;       Written by K. Willett            Aug 12    
;-

;restore, '~/Astronomy/Research/GalaxyZoo/gz2.sav'
;gz2main = gz2[where(strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra')]

count = 10

features = where(gz2main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge 0.8)
smooth = where(gz2main.t01_smooth_or_features_a01_smooth_weighted_fraction ge 0.8)
bar = where(gz2main.t03_bar_a06_bar_weighted_fraction ge 0.8)
nobar = where(gz2main.t03_bar_a07_no_bar_weighted_fraction ge 0.8)
nobulge = where(gz2main.T05_BULGE_PROMINENCE_A10_NO_BULGE_WEIGHTED_FRACTION ge 0.8 and gz2main.T05_BULGE_PROMINENCE_TOTAL_WEIGHT gt count)
justnoticeable = where(gz2main.T05_BULGE_PROMINENCE_A11_JUST_NOTICEABLE_WEIGHTED_FRACTION ge 0.8 and gz2main.T05_BULGE_PROMINENCE_TOTAL_WEIGHT gt count)
obvious = where(gz2main.T05_BULGE_PROMINENCE_A12_OBVIOUS_WEIGHTED_FRACTION ge 0.8 and gz2main.T05_BULGE_PROMINENCE_TOTAL_WEIGHT gt count)
dominant = where(gz2main.T05_BULGE_PROMINENCE_A13_DOMINANT_WEIGHTED_FRACTION ge 0.8 and gz2main.T05_BULGE_PROMINENCE_TOTAL_WEIGHT gt count)
rounded = where(gz2main.T07_ROUNDED_A16_COMPLETELY_ROUND_WEIGHTED_FRACTION ge 0.8)
inbetween = where(gz2main.T07_ROUNDED_A17_IN_BETWEEN_WEIGHTED_FRACTION ge 0.8)
cigar = where(gz2main.T07_ROUNDED_A18_CIGAR_SHAPED_WEIGHTED_FRACTION ge 0.8)

!p.multi=[0,2,2]

rmag = gz2main.petromag_mr
color = gz2main.petromag_mu - gz2main.petromag_mr

rmagarr = fillarr(0.2,-24,-18+0.1)
colorarr = fillarr(0.1,0.5,3.0+0.1)
smoothind = where(gz2main.T01_SMOOTH_OR_FEATURES_TOTAL_WEIGHT gt count)
barind = where(gz2main.T03_BAR_TOTAL_WEIGHT gt count)
bulgeind = where(gz2main.T05_BULGE_PROMINENCE_TOTAL_WEIGHT gt count)
roundedind = where(gz2main.T07_ROUNDED_TOTAL_WEIGHT gt count)

smoothhist = hist_2d(rmag[smoothind],color[smoothind], min1 = -24, max1=-18, min2 = 0.5, max2 = 3.0, bin1=0.2,bin2=0.1)
barhist = hist_2d(rmag[barind],color[barind], min1 = -24, max1=-18, min2 = 0.5, max2 = 3.0, bin1=0.2,bin2=0.1)
bulgehist = hist_2d(rmag[bulgeind],color[bulgeind], min1 = -24, max1=-18, min2 = 0.5, max2 = 3.0, bin1=0.2,bin2=0.1)
roundedhist = hist_2d(rmag[roundedind],color[roundedind], min1 = -24, max1=-18, min2 = 0.5, max2 = 3.0, bin1=0.2,bin2=0.1)

;; Spirals vs. ellipticals
;
;cgplot, rmag[smooth], color[smooth], $
;	xr=[-18,-24], /xstyle, $
;	yr = [0.5,3], $
;	xtitle='M!Ir!N', $
;	ytitle='u-r', $
;	title='Task 01 - smooth vs. features', $
;	color='red', $
;	charsize=2, $
;	psym = 3
;cgplot, rmag[features], color[features], psym=3, /overplot, color='black'
;
;cgcontour, /over, smoothhist, rmagarr, colorarr, thick=3, nlevels=10, label=0
;
;; Bar vs. no bar
;
;cgplot, rmag[nobar], color[nobar], $
;	xr=[-18,-24], /xstyle, $
;	yr = [0.5,3], $
;	xtitle='M!Ir!N', $
;	ytitle='u-r', $
;	title='Task 03 - bar vs. no bar', $
;	color='black', $
;	charsize=2, $
;	psym = 3
;cgplot, rmag[bar], color[bar], psym=16, symsize=0.2, /overplot, color='red'
;
;cgcontour, /over, barhist, rmagarr, colorarr, thick=3, nlevels=10, label=0
;
;; Bulge prominence
;
;cgplot, rmag[justnoticeable], color[justnoticeable], $
;	xr=[-18,-24], /xstyle, $
;	yr = [0.5,3], $
;	xtitle='M!Ir!N', $
;	ytitle='u-r', $
;	title='Task 05 - bulge prominence', $
;	color='green', $
;	charsize=2, $
;	symsize = 0.2, $
;	psym = 16
;cgplot, rmag[obvious], color[obvious], psym=16, symsize=0.2, /overplot, color='blue'
;cgplot, rmag[nobulge], color[nobulge], psym=16, symsize=0.2, /overplot, color='red'
;cgplot, rmag[dominant], color[dominant], psym=16, symsize=0.2, /overplot, color='black'
;
;cgcontour, /over, bulgehist, rmagarr, colorarr, thick=3, nlevels=10, label=0
;
;; Roundedness
;
;cgplot, rmag[rounded], color[rounded], $
;	xr=[-18,-24], /xstyle, $
;	yr = [0.5,3], $
;	xtitle='M!Ir!N', $
;	ytitle='u-r', $
;	title='Task 07 - roundedness', $
;	color='red', $
;	charsize=2, $
;	symsize=0.2, $ 
;	psym = 16
;
;cgplot, rmag[inbetween], color[inbetween], psym=16, symsize=0.2, /overplot, color='green'
;cgplot, rmag[cigar], color[cigar], psym=16, symsize=0.2, /overplot, color='blue'
;
;cgcontour, /over, roundedhist, rmagarr, colorarr, thick=3, nlevels=10, label=0

; Bulge prominence

cgplot, rmag[justnoticeable], color[justnoticeable], $
	xr=[-18,-24], /xstyle, $
	yr = [0.0,3], $
	xtitle='M!Ir!N', $
	ytitle='u-r', $
	title='No bulge', $
	color='dark green', $
	charsize=2, $
	symsize = 0.4, $
	psym = 16

cgcontour, /over, barhist, rmagarr, colorarr, thick=2, nlevels=10, label=0

cgplot, rmag[obvious], color[obvious], $
	xr=[-18,-24], /xstyle, $
	yr = [0.0,3], $
	xtitle='M!Ir!N', $
	ytitle='u-r', $
	title='Just noticeable', $
	color='blue', $
	charsize=2, $
	symsize = 0.6, $
	psym = 16

cgcontour, /over, barhist, rmagarr, colorarr, thick=2, nlevels=10, label=0

cgplot, rmag[nobulge], color[nobulge], $
	xr=[-18,-24], /xstyle, $
	yr = [0.0,3], $
	xtitle='M!Ir!N', $
	ytitle='u-r', $
	title='Obvious bulge', $
	color='red', $
	charsize=2, $
	symsize = 0.6, $
	psym = 16

cgcontour, /over, barhist, rmagarr, colorarr, thick=2, nlevels=10, label=0

cgplot, rmag[dominant], color[dominant],$
	xr=[-18,-24], /xstyle, $
	yr = [0.0,3], $
	xtitle='M!Ir!N', $
	ytitle='u-r', $
	title='Dominant bulge', $
	color='purple', $
	charsize=2, $
	symsize = 1, $
	psym = 16

cgcontour, /over, barhist, rmagarr, colorarr, thick=2, nlevels=10, label=0

stop

end
