
;+
; NAME:
;       
;	HC2
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
;       Written by K. Willett                
;-

pro hc2, ps=ps, stop=stop

gz2dir  = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = gz2dir+'fits/'
figdir  = gz2dir+'datapaper/figures/'
savdir  = gz2dir+'sav/'

restore,savdir+'gz2hc.sav'

tagnames=tag_names(gz2hc)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]	

gz2main = gz2hc[where((strtrim(gz2hc.sample,2) eq 'original' or strtrim(gz2hc.sample,2) eq 'extra'),maincount)]

gz2_el = where(gz2main.(wfind[1]) ge 0.8)
gz2_sp = where(gz2main.(wfind[2]) ge 0.8)
gz2_unc = where(gz2main.(wfind[1]) le 0.8 and gz2main.(wfind[2]) le 0.8 and gz2main.(wfind[3]) le 0.8)

; Plot the HC probabilities as function of bulge dominance

h1 = double(hist_2d(gz2main.pe_s0,gz2main.(wfind[10]),bin1 = 0.05, bin2 = 0.05))
h2 = double(hist_2d(gz2main.pe_s0,gz2main.(wfind[11]),bin1 = 0.05, bin2 = 0.05))
h3 = double(hist_2d(gz2main.pe_s0,gz2main.(wfind[12]),bin1 = 0.05, bin2 = 0.05))
h4 = double(hist_2d(gz2main.pe_s0,gz2main.(wfind[13]),bin1 = 0.05, bin2 = 0.05))

x = fillarr(0.05,0,1.00)
y = fillarr(0.05,0,1.05)

whiteback
!p.multi=[0,2,2]

ncolors = 50
cgloadct, 33, ncolors=ncolors, bottom=1

cgcontour, alog10(h1+1), x, y, /fill, nlevels=ncolors, c_colors=indgen(ncolors)+1, olevels=o1, ytitle='GZ2 bulge w.f.', title='No bulge'
cgcontour, alog10(h2+1), x, y, /fill, nlevels=ncolors, c_colors=indgen(ncolors)+1, olevels=o2, ytitle='GZ2 bulge w.f.', title='Just noticeable'
cgcontour, alog10(h3+1), x, y, /fill, nlevels=ncolors, c_colors=indgen(ncolors)+1, olevels=o3, ytitle='GZ2 bulge w.f.', title='Obvious'
cgcontour, alog10(h4+1), x, y, /fill, nlevels=ncolors, c_colors=indgen(ncolors)+1, olevels=o4, ytitle='GZ2 bulge w.f.', title='Dominant'

cgcolorbar, pos=[0.1,0.51,0.472,0.52], ncolors=ncolors, bottom=1, charsize=1.0, range=[min(alog10(h1+1)),max(alog10(h1+1))]
cgcolorbar, pos=[0.6,0.51,0.970,0.52], ncolors=ncolors, bottom=1, charsize=1.0, range=[min(alog10(h2+1)),max(alog10(h2+1))]
cgcolorbar, pos=[0.1,0.02,0.472,0.03], ncolors=ncolors, bottom=1, charsize=1.0, range=[min(alog10(h3+1)),max(alog10(h3+1))]
cgcolorbar, pos=[0.6,0.02,0.970,0.03], ncolors=ncolors, bottom=1, charsize=1.0, range=[min(alog10(h4+1)),max(alog10(h4+1))]

if keyword_set(stop) then stop

end

