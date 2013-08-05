
;+
; NAME:
;       
;	HC_BAR
;
; PURPOSE:
;
;	Plot Huertas-Company classifications vs. GZ2 bar fractions
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

pro hc_bar, ps=ps, stop=stop

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

;taskind = where(gz2main.T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION ge 0.430 and gz2main.T02_EDGEON_A05_NO_WEIGHTED_FRACTION ge 0.715 and gz2main.T03_BAR_TOTAL_COUNT ge 20,nt)
taskind = where(gz2main.T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION ge 0.430 and gz2main.T02_EDGEON_A05_NO_WEIGHTED_FRACTION ge 0.715,nt)
print,nt

if keyword_set(ps) then begin
	ps_start, filename=figdir+'hc_gz2_bar.eps',/color, /quiet,/encap, xsize=4, ysize=4
	cs=1.2
	lsize = 0.9
	th=3
	thickline=5
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
	lsize = 1
endelse

!p.multi=[0,1,1]

bs = 0.1
binarr = fillarr(bs,0,1.0)
val_el = binarr*0d
err_el = binarr*0d
val_s0 = binarr*0d
err_s0 = binarr*0d
val_ab = binarr*0d
err_ab = binarr*0d
val_cd = binarr*0d
err_cd = binarr*0d

disks = gz2main[taskind]

for i=0,n_elements(binarr)-1 do begin
    val_el[i] =   mean(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].pe)
    err_el[i] = stddev(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].pe)
    val_s0[i] =   mean(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].ps0)
    err_s0[i] = stddev(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].ps0)
    val_ab[i] =   mean(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].psab)
    err_ab[i] = stddev(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].psab)
    val_cd[i] =   mean(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].pscd)
    err_cd[i] = stddev(disks[where(disks.t03_bar_a06_bar_weighted_fraction ge binarr[i] and disks.t03_bar_a06_bar_weighted_fraction lt binarr[i]+bs)].pscd)
endfor

ploterror, binarr+bs/2., val_cd, err_cd, $
    psym = -4, $
    color='purple', $
    errcolor='purple', $
    charsize=cs, $
    thick = th, $
    errthick = th-2, $
    xthick = th, $
    ythick = th, $
    xr=[0,1], /xstyle, $
    yr=[0,1], /ystyle, $
    position=[0.17,0.15,0.95,0.95], $
    xtitle='GZ2 bar vote fraction', $
    ytitle='HC probability'

oploterror, binarr+bs/2., val_el, err_el, psym=-5, thick = th, errthick = th-2, color='red'
oploterror, binarr+bs/2., val_s0, err_s0, psym=-6, thick = th, errthick = th-2, color='blue'
oploterror, binarr+bs/2., val_ab, err_ab, psym=-7, thick = th, errthick = th-2, color='dark green'

al_legend, /top, /right, ['E','S0','Sab','Scd'], color=['red','blue','dark green','purple'], charsize=lsize, thick=th, psym=[-5,-6,-7,-4]

if keyword_set(ps) then ps_end,/pdf

if keyword_set(stop) then stop

end

