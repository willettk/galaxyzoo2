
;+
; NAME:
;       
;	NAIRABRAHAM
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
;	Numbers of barred, lensed, ringed galaxies are slightly different than reported in NA10. 
;	Color coding of ringed galaxies is incorrect (Fig 22)
;	Panel e) of Fig. 22 does not match the published data
;
; REVISION HISTORY
;       Written by K. Willett                Apr 12
;-

na = mrdfits('~/Astronomy/Research/GalaxyZoo/nair_simbad.fit',1,hdr,/silent)

;na = na[where(na.tt lt 99)]


; Figure 19: Histogrammed T-type as function of redshift

!p.multi=[0,1,2]
cs=2.0

sd = where(na.tt ge 7 and na.tt le 10)
sc = where(na.tt ge 5 and na.tt le 6)
sb = where(na.tt ge 3 and na.tt le 4)
sa = where(na.tt ge 1 and na.tt le 2)
s0 = where(na.tt ge -3 and na.tt le 0)
es = where(na.tt eq -5)

bs = 0.01

cgplot, indgen(10), $
	/nodata, $
	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='z', $
	ytitle='Number', $
	title='T-type', $
	xr=[0,0.10], $
	yr=[-500,4000],/ystyle

for i=0,8 do begin
	x0 = i*0.01
	x1 = (i+1)*0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.tt ge 7,csd)
	z_sc = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.tt ge 5 and na.tt le 6,csc)
	z_sb = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.tt ge 3 and na.tt le 4,csb)
	z_sa = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.tt ge 1 and na.tt le 2,csa)
	z_s0 = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.tt ge -3 and na.tt le 0,cs0)
	z_es = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.tt ge -5 and na.tt le -4,ces)

	;print,csd,csc,csb,csa,cs0,ces

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif


endfor

cgplot, indgen(10), $
	/nodata, $
	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='z', $
	ytitle='Fraction', $
	xr=[0,0.10], $
	yr=[-0.2,1.2],/ystyle

for i=0,8 do begin
	x0 = i*0.01
	x1 = (i+1)*0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(na.zs gt x0 and na.zs le x1 and na.tt ge 7,csd)
	z_sc = where(na.zs gt x0 and na.zs le x1 and na.tt ge 5 and na.tt le 6,csc)
	z_sb = where(na.zs gt x0 and na.zs le x1 and na.tt ge 3 and na.tt le 4,csb)
	z_sa = where(na.zs gt x0 and na.zs le x1 and na.tt ge 1 and na.tt le 2,csa)
	z_s0 = where(na.zs gt x0 and na.zs le x1 and na.tt ge -3 and na.tt le 0,cs0)
	z_es = where(na.zs gt x0 and na.zs le x1 and na.tt ge -5 and na.tt le -4,ces)

	junk =  where(na.zs gt i*0.01 and na.zs le (i+1)*0.01,c_all)

	;print,csd,csc,csb,csa,cs0,ces

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

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=1.5


endfor



stop


; Figure 20: Histogrammed stellar mass as function of redshift

!p.multi=[0,1,2]
cs=2.0

bs = 0.01

cgplot, indgen(10), $
	/nodata, $
	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='z', $
	ytitle='Number', $
	title='Stellar mass', $
	xr=[0,0.10], $
	yr=[-500,4000],/ystyle

for i=0,8 do begin
	x0 = i*0.01
	x1 = (i+1)*0.01
	y0 = 0.
	y1 = 0.
	z_sm_blue = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 8. and na.log_m_ le 9.5,c1)
	z_sm_darkgreen = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 9.5 and na.log_m_ le 10.2,c2)
	z_sm_green = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 10.2 and na.log_m_ le 10.5,c3)
	z_sm_yellow = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 10.5 and na.log_m_ le 10.8,c4)
	z_sm_tomato = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 10.8 and na.log_m_ le 11.2,c5)
	z_sm_black = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 11.2,c6)

	;print,csd,csc,csb,csa,cs0,ces

	if c6 gt 0 then begin
		y0 = y1
		y1 = y0 + c6
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if c5 gt 0 then begin
		y0 = y1
		y1 = y0 + c5
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if c4 gt 0 then begin
		y0 = y1
		y1 = y0 + c4
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if c3 gt 0 then begin
		y0 = y1
		y1 = y0 + c3
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if c2 gt 0 then begin
		y0 = y1
		y1 = y0 + c2
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if c1 gt 0 then begin
		y0 = y1
		y1 = y0 + c1
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif


endfor

cgplot, indgen(10), $
	/nodata, $
	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='z', $
	ytitle='Fraction', $
	xr=[0,0.10], $
	yr=[-0.2,1.2],/ystyle

for i=0,8 do begin
	x0 = i*0.01
	x1 = (i+1)*0.01
	y0 = 0.
	y1 = 0.
	z_sm_blue = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 8. and na.log_m_ le 9.5,c1)
	z_sm_darkgreen = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 9.5 and na.log_m_ le 10.2,c2)
	z_sm_green = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 10.2 and na.log_m_ le 10.5,c3)
	z_sm_yellow = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 10.5 and na.log_m_ le 10.8,c4)
	z_sm_tomato = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 10.8 and na.log_m_ le 11.2,c5)
	z_sm_black = where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 11.2,c6)

	junk =  where(na.zs gt i*0.01 and na.zs le (i+1)*0.01 and na.log_m_ gt 8.,c_all)

	;print,csd,csc,csb,csa,cs0,ces

	if c6 gt 0 then begin
		y0 = y1
		y1 = y0 + c6/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if c5 gt 0 then begin
		y0 = y1
		y1 = y0 + c5/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if c4 gt 0 then begin
		y0 = y1
		y1 = y0 + c4/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if c3 gt 0 then begin
		y0 = y1
		y1 = y0 + c3/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if c2 gt 0 then begin
		y0 = y1
		y1 = y0 + c2/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if c1 gt 0 then begin
		y0 = y1
		y1 = y0 + c1/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=1.5


endfor










; SEEING

bar = where(na.bar gt 1)
strongbar = where(na.bar eq 2)
medbar = where(na.bar eq 4)
weakbar = where(na.bar eq 8)

!p.multi=[0,3,2]

; Nbar vs. seeing

bs = 0.1
cghistoplot, na[bar].seeing, $
	binsize=bs, $
	xr=[0,3], $
	yr=[0,500], $
	title='Bars', $
	/outline,$
	datacolor='black'

cghistoplot, na[medbar].seeing, $
	binsize=bs, $
	/outline, $
	/fillpolygon, $
	polycolor='purple', $
	/oplot

cghistoplot, na[strongbar].seeing, $
	binsize=bs, $
	/fillpolygon, $
	/outline, $
	polycolor='red', $
	/oplot

cghistoplot, na[weakbar].seeing, $
	binsize=bs, $
	/outline, $
	/oplot, $
	datacolor='cyan'

print,''
print,'Bars: ',n_elements(bar)
print,'Weak bar: ',n_elements(weakbar)
print,'Medium bar: ',n_elements(medbar)
print,'Strong bar: ',n_elements(strongbar)

; Nring vs. seeing

ring = where(na.ring gt 1)
innerring = where(na.ring eq 4)
outerring = where(na.ring eq 8)
ioring = where(na.ring eq 12)

cghistoplot, na[ring].seeing, $
	binsize=bs, $
	xr=[0,3], $
	yr=[0,500], $
	title='Rings', $
	/outline,$
	datacolor='black'

cghistoplot, na[innerring].seeing, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='blue', $
	/outline, $
	/oplot

cghistoplot, na[outerring].seeing, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='purple', $
	/oplot

cghistoplot, na[ioring].seeing, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='red', $
	/oplot

print,''
print,'Rings: ',n_elements(ring)
print,'Inner ring: ',n_elements(innerring)
print,'Outer ring: ',n_elements(outerring)
print,'I+O ring: ',n_elements(ioring)

; Nlens vs. seeing

lens = where(na.lens gt 1)
innerlens = where(na.lens eq 2)
outerlens = where(na.lens eq 4)

cghistoplot, na[lens].seeing, $
	binsize=bs, $
	xr=[0,3], $
	yr=[0,120], $
	title='Lenses', $
	/outline,$
	datacolor='black'

cghistoplot, na[innerlens].seeing, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='blue', $
	/outline, $
	/oplot

cghistoplot, na[outerlens].seeing, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='purple', $
	datacolor='purple', $
	/oplot

cghistoplot, na[innerlens].seeing, $
	binsize=bs, $
	/outline, $
	datacolor='blue',$
	/oplot

print,''
print,'Lenses: ',n_elements(lens)
print,'Inner lens: ',n_elements(innerlens)
print,'Outer lens: ',n_elements(outerlens)

; Fractional bar vs. seeing

barhist = histogram(na.seeing,binsize=0.1,min=1.0,max=2.2,omin=omin,reverse_indices=ri)
nb = n_elements(barhist)

barfrac = fltarr(nb)
barfrac_weak = fltarr(nb)
barfrac_med = fltarr(nb)
barfrac_strong = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nbar = where(na[ind].bar gt 1,barcount)
		nbar_weak = where(na[ind].bar eq 8,barcount_weak)
		nbar_med = where(na[ind].bar eq 4,barcount_med)
		nbar_strong = where(na[ind].bar eq 2,barcount_strong)
		barfrac[i] = float(barcount) / n_elements(ind)
		barfrac_weak[i] = float(barcount_weak) / n_elements(ind)
		barfrac_med[i] = float(barcount_med) / n_elements(ind)
		barfrac_strong[i] = float(barcount_strong) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,barfrac, $
	psym=10, $
	xr=[0,3], $
	yr=[0,0.4], $
	ytitle='f!Ibar!N'
cgplots,[omin,omin],[barfrac[0],0]

cgplot, xarr, barfrac_weak, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[barfrac_weak[0],0],color='blue'

cgplot, xarr, barfrac_med, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[barfrac_med[0],0],color='purple'

cgplot, xarr, barfrac_strong, $
	/overplot, $
	psym=10, $
	color='red'
cgplots,[omin,omin],[barfrac_strong[0],0],color='red'

; Fractional ring vs. seeing

ringhist = histogram(na.seeing,binsize=0.1,min=1.0,max=2.4,omin=omin,reverse_indices=ri)
nb = n_elements(ringhist)

ringfrac = fltarr(nb)
ringfrac_inner = fltarr(nb)
ringfrac_outer = fltarr(nb)
ringfrac_io = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nring = where(na[ind].ring gt 1,ringcount)
		nring_inner = where(na[ind].ring eq 2,ringcount_inner)
		nring_outer = where(na[ind].ring eq 4,ringcount_outer)
		nring_io = where(na[ind].ring eq 6,ringcount_io)
		ringfrac[i] = float(ringcount) / n_elements(ind)
		ringfrac_inner[i] = float(ringcount_inner) / n_elements(ind)
		ringfrac_outer[i] = float(ringcount_outer) / n_elements(ind)
		ringfrac_io[i] = float(ringcount_io) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,ringfrac, $
	psym=10, $
	xr=[0,3], $
	yr=[0,0.4], $
	ytitle='f!Iring!N'
cgplots,[omin,omin],[ringfrac[0],0]

cgplot, xarr, ringfrac_inner, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[ringfrac_inner[0],0],color='blue'

cgplot, xarr, ringfrac_outer, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[ringfrac_outer[0],0],color='purple'

cgplot, xarr, ringfrac_io, $
	/overplot, $
	psym=10, $
	color='red'
cgplots,[omin,omin],[ringfrac_io[0],0],color='red'

; Fractional lens vs. seeing

lenshist = histogram(na.seeing,binsize=0.1,min=1.0,max=2.0,omin=omin,reverse_indices=ri)
nb = n_elements(lenshist)

lensfrac = fltarr(nb)
lensfrac_inner = fltarr(nb)
lensfrac_outer = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nlens = where(na[ind].lens gt 1,lenscount)
		nlens_inner = where(na[ind].lens eq 2,lenscount_inner)
		nlens_outer = where(na[ind].lens eq 4,lenscount_outer)
		lensfrac[i] = float(lenscount) / n_elements(ind)
		lensfrac_inner[i] = float(lenscount_inner) / n_elements(ind)
		lensfrac_outer[i] = float(lenscount_outer) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,lensfrac, $
	psym=10, $
	xr=[0,3], $
	yr=[0,0.1], $
	ytitle='f!Ilens!N'
cgplots,[omin,omin],[lensfrac[0],0]

cgplot, xarr, lensfrac_inner, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[lensfrac_inner[0],0],color='blue'

cgplot, xarr, lensfrac_outer, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[lensfrac_outer[0],0],color='purple'





; B/A





!p.multi=[0,3,2]

; Nbar vs. b_a

bs = 0.05
cghistoplot, na[bar].b_a, $
	binsize=bs, $
	xr=[0,1], $
	yr=[0,300], $
	title='Bars', $
	xtitle='b/a', $
	/outline,$
	datacolor='black'

cghistoplot, na[medbar].b_a, $
	binsize=bs, $
	/outline, $
	/fillpolygon, $
	polycolor='purple', $
	/oplot

cghistoplot, na[strongbar].b_a, $
	binsize=bs, $
	/fillpolygon, $
	/outline, $
	polycolor='red', $
	/oplot

cghistoplot, na[weakbar].b_a, $
	binsize=bs, $
	/outline, $
	/oplot, $
	datacolor='cyan'

print,''
print,'Bars: ',n_elements(bar)
print,'Weak bar: ',n_elements(weakbar)
print,'Medium bar: ',n_elements(medbar)
print,'Strong bar: ',n_elements(strongbar)

; Nring vs. b_a

ring = where(na.ring gt 1)
innerring = where(na.ring eq 4)
outerring = where(na.ring eq 8)
ioring = where(na.ring eq 12)

cghistoplot, na[ring].b_a, $
	binsize=bs, $
	xr=[0,1], $
	yr=[0,300], $
	title='Rings', $
	xtitle='b/a', $
	/outline,$
	datacolor='black'

cghistoplot, na[innerring].b_a, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='blue', $
	/outline, $
	/oplot

cghistoplot, na[outerring].b_a, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='purple', $
	/oplot

cghistoplot, na[ioring].b_a, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='red', $
	/oplot

print,''
print,'Rings: ',n_elements(ring)
print,'Inner ring: ',n_elements(innerring)
print,'Outer ring: ',n_elements(outerring)
print,'I+O ring: ',n_elements(ioring)

; Nlens vs. b_a

lens = where(na.lens gt 1)
innerlens = where(na.lens eq 2)
outerlens = where(na.lens eq 4)

cghistoplot, na[lens].b_a, $
	binsize=bs, $
	xr=[0,1], $
	yr=[0,70], $
	title='Lenses', $
	xtitle='b/a', $
	/outline,$
	datacolor='black'

cghistoplot, na[innerlens].b_a, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='blue', $
	/outline, $
	/oplot

cghistoplot, na[outerlens].b_a, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='purple', $
	datacolor='purple', $
	/oplot

cghistoplot, na[innerlens].b_a, $
	binsize=bs, $
	/outline, $
	datacolor='blue',$
	/oplot

print,''
print,'Lenses: ',n_elements(lens)
print,'Inner lens: ',n_elements(innerlens)
print,'Outer lens: ',n_elements(outerlens)

; Fractional bar vs. b_a

barhist = histogram(na.b_a,binsize=0.05,min=0.15,max=1.0,omin=omin,reverse_indices=ri)
nb = n_elements(barhist)

barfrac = fltarr(nb)
barfrac_weak = fltarr(nb)
barfrac_med = fltarr(nb)
barfrac_strong = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nbar = where(na[ind].bar gt 1,barcount)
		nbar_weak = where(na[ind].bar eq 8,barcount_weak)
		nbar_med = where(na[ind].bar eq 4,barcount_med)
		nbar_strong = where(na[ind].bar eq 2,barcount_strong)
		barfrac[i] = float(barcount) / n_elements(ind)
		barfrac_weak[i] = float(barcount_weak) / n_elements(ind)
		barfrac_med[i] = float(barcount_med) / n_elements(ind)
		barfrac_strong[i] = float(barcount_strong) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,barfrac, $
	psym=10, $
	xr=[0,1], $
	yr=[0,0.4], $
	xtitle='b/a', $
	ytitle='f!Ibar!N'
cgplots,[omin,omin],[barfrac[0],0]

cgplot, xarr, barfrac_weak, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[barfrac_weak[0],0],color='blue'

cgplot, xarr, barfrac_med, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[barfrac_med[0],0],color='purple'

cgplot, xarr, barfrac_strong, $
	/overplot, $
	psym=10, $
	color='red'
cgplots,[omin,omin],[barfrac_strong[0],0],color='red'

; Fractional ring vs. b_a

ringhist = histogram(na.b_a,binsize=0.05,min=0.15,max=1.0,omin=omin,reverse_indices=ri)
nb = n_elements(ringhist)

ringfrac = fltarr(nb)
ringfrac_inner = fltarr(nb)
ringfrac_outer = fltarr(nb)
ringfrac_io = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nring = where(na[ind].ring gt 1,ringcount)
		nring_inner = where(na[ind].ring eq 2,ringcount_inner)
		nring_outer = where(na[ind].ring eq 4,ringcount_outer)
		nring_io = where(na[ind].ring eq 6,ringcount_io)
		ringfrac[i] = float(ringcount) / n_elements(ind)
		ringfrac_inner[i] = float(ringcount_inner) / n_elements(ind)
		ringfrac_outer[i] = float(ringcount_outer) / n_elements(ind)
		ringfrac_io[i] = float(ringcount_io) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,ringfrac, $
	psym=10, $
	xr=[0,1], $
	yr=[0,0.4], $
	xtitle='b/a', $
	ytitle='f!Iring!N'
cgplots,[omin,omin],[ringfrac[0],0]

cgplot, xarr, ringfrac_inner, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[ringfrac_inner[0],0],color='blue'

cgplot, xarr, ringfrac_outer, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[ringfrac_outer[0],0],color='purple'

cgplot, xarr, ringfrac_io, $
	/overplot, $
	psym=10, $
	color='red'
cgplots,[omin,omin],[ringfrac_io[0],0],color='red'

; Fractional lens vs. b_a

lenshist = histogram(na.b_a,binsize=0.05,min=0.15,max=1.0,omin=omin,reverse_indices=ri)
nb = n_elements(lenshist)

lensfrac = fltarr(nb)
lensfrac_inner = fltarr(nb)
lensfrac_outer = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nlens = where(na[ind].lens gt 1,lenscount)
		nlens_inner = where(na[ind].lens eq 2,lenscount_inner)
		nlens_outer = where(na[ind].lens eq 4,lenscount_outer)
		lensfrac[i] = float(lenscount) / n_elements(ind)
		lensfrac_inner[i] = float(lenscount_inner) / n_elements(ind)
		lensfrac_outer[i] = float(lenscount_outer) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,lensfrac, $
	psym=10, $
	xr=[0,1], $
	yr=[0,0.1], $
	xtitle='b/a', $
	ytitle='f!Ilens!N'
cgplots,[omin,omin],[lensfrac[0],0]

cgplot, xarr, lensfrac_inner, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[lensfrac_inner[0],0],color='blue'

cgplot, xarr, lensfrac_outer, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[lensfrac_outer[0],0],color='purple'






; REDSHIFT






!p.multi=[0,3,2]

; Nbar vs. zs

bs = 0.005
xmin = 0.0
xmax = 0.10

cghistoplot, na[bar].zs, $
	binsize=bs, $
	xr = [xmin,xmax], $
	yr=[0,500], $
	title='Bars', $
	xtitle='z', $
	/outline,$
	datacolor='black'

cghistoplot, na[medbar].zs, $
	binsize=bs, $
	/outline, $
	/fillpolygon, $
	polycolor='purple', $
	/oplot

cghistoplot, na[strongbar].zs, $
	binsize=bs, $
	/fillpolygon, $
	/outline, $
	polycolor='red', $
	/oplot

cghistoplot, na[weakbar].zs, $
	binsize=bs, $
	/outline, $
	/oplot, $
	datacolor='cyan'

print,''
print,'Bars: ',n_elements(bar)
print,'Weak bar: ',n_elements(weakbar)
print,'Medium bar: ',n_elements(medbar)
print,'Strong bar: ',n_elements(strongbar)

; Nring vs. zs

ring = where(na.ring gt 1)
innerring = where(na.ring eq 4)
outerring = where(na.ring eq 8)
ioring = where(na.ring eq 12)

cghistoplot, na[ring].zs, $
	binsize=bs, $
	xr = [xmin,xmax], $
	yr=[0,500], $
	title='Rings', $
	xtitle='z', $
	/outline,$
	datacolor='black'

cghistoplot, na[innerring].zs, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='blue', $
	/outline, $
	/oplot

cghistoplot, na[outerring].zs, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='purple', $
	/oplot

cghistoplot, na[ioring].zs, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='red', $
	/oplot

print,''
print,'Rings: ',n_elements(ring)
print,'Inner ring: ',n_elements(innerring)
print,'Outer ring: ',n_elements(outerring)
print,'I+O ring: ',n_elements(ioring)

; Nlens vs. zs

lens = where(na.lens gt 1)
innerlens = where(na.lens eq 2)
outerlens = where(na.lens eq 4)

cghistoplot, na[lens].zs, $
	binsize=bs, $
	xr = [xmin,xmax], $
	yr=[0,100], $
	title='Lenses', $
	xtitle='z', $
	/outline,$
	datacolor='black'

cghistoplot, na[innerlens].zs, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='blue', $
	/outline, $
	/oplot

cghistoplot, na[outerlens].zs, $
	/outline, $
	binsize=bs, $
	/fillpolygon, $
	polycolor='purple', $
	datacolor='purple', $
	/oplot

cghistoplot, na[innerlens].zs, $
	binsize=bs, $
	/outline, $
	datacolor='blue',$
	/oplot

print,''
print,'Lenses: ',n_elements(lens)
print,'Inner lens: ',n_elements(innerlens)
print,'Outer lens: ',n_elements(outerlens)

; Fractional bar vs. zs

hbs = 0.005
hmin = 0.01
hmax = 0.090

barhist = histogram(na.zs,binsize=hbs,min=hmin,max=hmax,omin=omin,reverse_indices=ri)
nb = n_elements(barhist)

barfrac = fltarr(nb)
barfrac_weak = fltarr(nb)
barfrac_med = fltarr(nb)
barfrac_strong = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nbar = where(na[ind].bar gt 1,barcount)
		nbar_weak = where(na[ind].bar eq 8,barcount_weak)
		nbar_med = where(na[ind].bar eq 4,barcount_med)
		nbar_strong = where(na[ind].bar eq 2,barcount_strong)
		barfrac[i] = float(barcount) / n_elements(ind)
		barfrac_weak[i] = float(barcount_weak) / n_elements(ind)
		barfrac_med[i] = float(barcount_med) / n_elements(ind)
		barfrac_strong[i] = float(barcount_strong) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,barfrac, $
	psym=10, $
	xr = [xmin,xmax], $
	yr=[0,0.4], $
	xtitle='z', $
	ytitle='f!Ibar!N'
cgplots,[omin,omin],[barfrac[0],0]

cgplot, xarr, barfrac_weak, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[barfrac_weak[0],0],color='blue'

cgplot, xarr, barfrac_med, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[barfrac_med[0],0],color='purple'

cgplot, xarr, barfrac_strong, $
	/overplot, $
	psym=10, $
	color='red'
cgplots,[omin,omin],[barfrac_strong[0],0],color='red'

; Fractional ring vs. zs

ringhist = histogram(na.zs,binsize=hbs,min=hmin,max=hmax,omin=omin,reverse_indices=ri)
nb = n_elements(ringhist)

ringfrac = fltarr(nb)
ringfrac_inner = fltarr(nb)
ringfrac_outer = fltarr(nb)
ringfrac_io = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nring = where(na[ind].ring gt 1,ringcount)
		nring_inner = where(na[ind].ring eq 2,ringcount_inner)
		nring_outer = where(na[ind].ring eq 4,ringcount_outer)
		nring_io = where(na[ind].ring eq 6,ringcount_io)
		ringfrac[i] = float(ringcount) / n_elements(ind)
		ringfrac_inner[i] = float(ringcount_inner) / n_elements(ind)
		ringfrac_outer[i] = float(ringcount_outer) / n_elements(ind)
		ringfrac_io[i] = float(ringcount_io) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,ringfrac, $
	psym=10, $
	xr = [xmin,xmax], $
	yr=[0,0.4], $
	xtitle='z', $
	ytitle='f!Iring!N'
cgplots,[omin,omin],[ringfrac[0],0]

cgplot, xarr, ringfrac_inner, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[ringfrac_inner[0],0],color='blue'

cgplot, xarr, ringfrac_outer, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[ringfrac_outer[0],0],color='purple'

cgplot, xarr, ringfrac_io, $
	/overplot, $
	psym=10, $
	color='red'
cgplots,[omin,omin],[ringfrac_io[0],0],color='red'

; Fractional lens vs. zs

lenshist = histogram(na.zs,binsize=hbs,min=hmin,max=hmax,omin=omin,reverse_indices=ri)
nb = n_elements(lenshist)

lensfrac = fltarr(nb)
lensfrac_inner = fltarr(nb)
lensfrac_outer = fltarr(nb)

for i=0,nb-1 do begin
	ind1 = ri[i]
	ind2 = ri[i+1]
	if ind1 ne ind2 then begin
		ind = ri[ri[i]:ri[i+1]-1]
		nlens = where(na[ind].lens gt 1,lenscount)
		nlens_inner = where(na[ind].lens eq 2,lenscount_inner)
		nlens_outer = where(na[ind].lens eq 4,lenscount_outer)
		lensfrac[i] = float(lenscount) / n_elements(ind)
		lensfrac_inner[i] = float(lenscount_inner) / n_elements(ind)
		lensfrac_outer[i] = float(lenscount_outer) / n_elements(ind)
	endif
endfor

xarr = indgen(nb)*bs + omin

cgplot, xarr,lensfrac, $
	psym=10, $
	xr = [xmin,xmax], $
	yr=[0,0.1], $
	xtitle='z', $
	ytitle='f!Ilens!N'
cgplots,[omin,omin],[lensfrac[0],0]

cgplot, xarr, lensfrac_inner, $
	/overplot, $
	psym=10, $
	color='blue'
cgplots,[omin,omin],[lensfrac_inner[0],0],color='blue'

cgplot, xarr, lensfrac_outer, $
	/overplot, $
	psym=10, $
	color='purple'
cgplots,[omin,omin],[lensfrac_outer[0],0],color='purple'



; Figure 24 - histograms of T-type by fine structure, active class

!p.multi=[0,2,1]

cghistoplot, na.tt, $
	binsize=1, $
	datacolor='black', $
	xtitle='T-type', $
	/outline, $
	yr=[0,3000], $
	xr=[-5,11], $
	xstyle=9, $
	charsize=cs, $
	xticks=7, xtickv=[-4.5,-2,1,3,5,7,8.5,9.5],xtickname=['E','S0','Sa','Sb','Sc','Sd','Sm','Im']

cghistoplot, na[bars].tt, $
	/oplot, $
	datacolor='blue', $
	/outline, $
	/fillpoly,polycolor='blue'

cghistoplot, na[ring].tt, $
	/oplot, $
	datacolor='purple', $
	/outline, $
	/fillpoly,polycolor='purple'

cghistoplot, na[bars].tt, $
	/oplot, $
	datacolor='blue', $
	/outline

cghistoplot, na[lens].tt, $
	/oplot, $
	datacolor='red', $
	/outline, $
	/fillpoly,polycolor='red'

al_legend, psym=28, colors=['black','blue','purple','red'], ['All','Bar','Ring','Lens'], charsize=cs, /top, /right

; By AGN type (Kauffmann criterion)

agn = where(na.agn1 gt 0)
liner = where(na.agn1 eq 3)
seyfert = where(na.agn1 eq 2)

cghistoplot, na.tt, $
	binsize=1, $
	datacolor='black', $
	xtitle='T-type', $
	/outline, $
	yr=[0,3000], $
	xr=[-5,11], $
	xstyle=9, $
	charsize=cs, $
	xticks=7, xtickv=[-4.5,-2,1,3,5,7,8.5,9.5],xtickname=['E','S0','Sa','Sb','Sc','Sd','Sm','Im']

cghistoplot, na[agn].tt, $
	/oplot, $
	datacolor='blue', $
	/outline, $
	/fillpoly,polycolor='blue'

cghistoplot, na[liner].tt, $
	/oplot, $
	datacolor='purple', $
	/outline, $
	/fillpoly,polycolor='purple'

cghistoplot, na[seyfert].tt, $
	/oplot, $
	datacolor='red', $
	/outline, $
	/fillpoly,polycolor='red'

al_legend, psym=28, colors=['black','blue','purple','red'], ['All','AGN','LINER','Seyfert'], charsize=cs, /top, /right


; Table 3

e = where(na.tt eq -5)
es0 = where(na.tt eq -3)


;print, 'E    ', n_elements(where(na.tt eq -5)), 


stop

end
