
;+
; NAME:
;       
;	EDGEON_CIGAR
;
; PURPOSE:
;
;	Examine the overlap between edge-on disks and cigar-shaped smooth galaxies in GZ2
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

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
figdir = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'

;gz2 = mrdfits(fitsdir+'gz2main_table_sample.fits',1,/silent)

tagnames=tag_names(gz2)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]-2	

x_taskweight = gz2.T07_ROUNDED_TOTAL_WEIGHT
y_taskweight = gz2.T02_EDGEON_TOTAL_WEIGHT

count=20
ind = where(x_taskweight gt count and y_taskweight ge count,nind)

x=gz2[ind].T07_ROUNDED_A18_CIGAR_SHAPED_WEIGHTED_FRACTION 
y=gz2[ind].T02_EDGEON_A04_YES_WEIGHTED_FRACTION                                     

binsize=0.05
xarr = fillarr(binsize,0.,1.)
yarr = fillarr(binsize,0.,1.)

h = hist_2d(x,y,bin1=binsize,bin2=binsize,min1=0,min2=0,max1=1,max2=1)

cgloadct,3

!p.multi=[0,1,2]

cgcontour,alog10(h+1), xarr, yarr, /fill, position=[0.1,0.55,0.80,0.95], xtitle='Cigar-shaped fraction', ytitle='Edge-on fraction'
cgimage,alog10(h),xarr,yarr,/scale,/axes,axkeywords={xrange:[0,1],yrange:[0,1]}, position=[0.1,0.14,0.80,0.9], xtitle='Cigar-shaped fraction', ytitle='Edge-on fraction'
cgcolorbar,/vert, range=[alog10(min(h)+1),alog10(max(h)+1)], title='log(N!Igal!N + 1)'

cgtext,0.85,0.94,'N = '+strtrim(nind,2), /normal, charsize=2

; How many galaxies have large numbers of votes for both edge-on and cigar-shaped?
; What fraction of edge-on galaxies extend 


end
