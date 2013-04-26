
;+
; NAME:
;       
;	HC_BULGE
;
; PURPOSE:
;
;	Plot Huertas-Company classifications vs. GZ2 vote fractions for bulge prominence
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
;       Written by K. Willett                Apr 13
;-

pro hc_bulge, ps=ps, stop=stop

gz2dir  = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = gz2dir+'fits/'
figdir  = gz2dir+'datapaper/figures/'
savdir  = gz2dir+'sav/'

gz2hcfile = fitsdir+'gz2_sample_table_hc.fits'

restore,savdir+'gz2hc.sav'

tagnames=tag_names(gz2hc)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]-2	

gz2main = gz2hc[where((strtrim(gz2hc.sample,2) eq 'original' or strtrim(gz2hc.sample,2) eq 'extra'),maincount)]

; Plot the HC probabilities as function of bulge dominance

mintask = 10
minprob = 0.5
taskind = where(gz2main.T02_EDGEON_A05_NO_WEIGHT ge 30,nt)
taskind = where(gz2main.T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION ge 0.430 and gz2main.T02_EDGEON_A05_NO_WEIGHTED_FRACTION ge 0.715,nt)
print,nt

h1 = double(hist_2d(gz2main[taskind].(wfind[10]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))
h2 = double(hist_2d(gz2main[taskind].(wfind[11]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))
h3 = double(hist_2d(gz2main[taskind].(wfind[12]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))
h4 = double(hist_2d(gz2main[taskind].(wfind[13]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))

y = fillarr(0.05, 0.0, 1.05)
x = fillarr(0.05, 0.0, 1.05)
whiteback

if keyword_set(ps) then begin
	ps_start, filename=figdir+'hc_gz2_bulge_contour.eps',/color, /quiet,/encap, xsize=8, ysize=7
	cs=1.5
	cbsize = 0.9
	th=3
	thickline=5
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
	cbsize = 1
endelse

!p.multi=[0,2,2]

ncolors = 50
cgloadct, 33, ncolors=ncolors, bottom=1
levels = fillarr(0.07,0,3.5)
;levels = fillarr(10,1,500)

cgcontour, alog10(h1+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='No bulge'
cgcontour, alog10(h2+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='Just noticeable'
cgcontour, alog10(h3+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='Obvious bulge'
cgcontour, alog10(h4+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='Dominant bulge'

cgcolorbar,position=[0.46,0.50,0.60,0.54],title='log(N!Igal!N + 1)',ncolors=ncolors,charsize=cbsize

if keyword_set(ps) then ps_end

if keyword_set(stop) then stop

end
