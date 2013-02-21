
;+
; NAME:
;       
;	HUERTASCOMPANY
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

pro huertascompany, ps=ps, stop=stop

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
figdir = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'
gz2hcfile = fitsdir+'gz2_sample_table_hc.fits'

restore,gz2dir+'gz2hc.sav'

tagnames=tag_names(gz2hc)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]-2	

gz2main = gz2hc[where((strtrim(gz2hc.sample,2) eq 'original' or strtrim(gz2hc.sample,2) eq 'extra'),maincount)]

; Plot thresholded ellipticals and spirals against HC classifications

gz2_el = where(gz2main.(wfind[1]) ge 0.8)
gz2_sp = where(gz2main.(wfind[2]) ge 0.8)
gz2_unc = where(gz2main.(wfind[1]) le 0.8 and gz2main.(wfind[2]) le 0.8 and gz2main.(wfind[3]) le 0.8)

;gz2_el = where(gz2main.(wfind[1]) ge 0.8 and gz2main.petromag_r lt 16)			; Apparent magnitude cut
;gz2_sp = where(gz2main.(wfind[2]) ge 0.8 and gz2main.petromag_r lt 16)
; gz2_el = where(gz2main.(wfind[1]) ge 0.8 and gz2main.redshift lt 0.05)		; Low-redshift galaxies
; gz2_sp = where(gz2main.(wfind[2]) ge 0.8 and gz2main.redshift lt 0.05)
; gz2_el = where(gz2main.(wfind[1]) ge 0.8 and gz2main.(wfind[17]) lt 0.1 and gz2main.(wfind[18]) lt 0.1)	; Galaxies that are smooth and round/
; gz2_sp = where(gz2main.(wfind[2]) ge 0.8 and gz2main.(wfind[17]) lt 0.1 and gz2main.(wfind[18]) lt 0.1)
; gz2_el = where(gz2main.(wfind[1]) ge 0.8 and gz2main.(wfind[14]) lt 0.25)		; Only look at galaxies that are not odd. 
; gz2_sp = where(gz2main.(wfind[2]) ge 0.8 and gz2main.(wfind[14]) lt 0.25)

!p.multi=[0,3,1]

bs = 0.04

if keyword_set(ps) then begin
	ps_start, filename=figdir+'hc_histogram.eps',/color, /quiet, /encap, xsize=7, ysize=2.5
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

cghistoplot, gz2main[gz2_el].pe_s0, $
	xr=[-0.05,1.05], /xstyle, $
	yr=[-0.01,0.30], /ystyle, $
	thick=th, $
	/freq, $
	/outline, $
	binsize = bs, $
	datacolor='red', $
	title='GZ2 smooth', $
	xtitle='HC early-type probability'

;cghistoplot, 1. - gz2main[gz2_el].pe_s0, $
;	thick=th, $
;	/freq, $
;	/outline, $
;	binsize = bs, $
;	/oplot, $
;	datacolor='blue'

cghistoplot, 1. - gz2main[gz2_sp].pe_s0, $
	xr=[-0.05,1.05], /xstyle, $
	yr=[-0.01,0.30], /ystyle, $
	thick=th, $
	/freq, $
	/outline, $
	binsize = bs, $
	datacolor='blue', $
	title='GZ2 features/disk', $
	xtitle='HC late-type probability'

;cghistoplot, 1. - gz2main[gz2_sp].pe_s0, $
;	thick=th, $
;	/freq, $
;	/outline, $
;	binsize = bs, $
;	/oplot, $
;	datacolor='blue'

cghistoplot, 1d - gz2main[gz2_unc].pe_s0, $
	xr=[-0.05,1.05], /xstyle, $
	yr=[-0.01,0.30], /ystyle, $
	thick=th, $
	/freq, $
	/outline, $
	binsize = bs, $
	datacolor='black', $
	title='GZ2 uncertain', $
	xtitle='HC late-type probability'

if keyword_set(ps) then ps_end

; Interesting. Changing the threshold doesn't change the distribution of HC probabilities for late-type galaxies, but it DOES increase the fraction of low-prob HC early-types that GZ2 users classified as being smooth. Why would HC find a population of galaxies they classify as late-type, but GZ2 users insist are smooth?

	; How many galaxies are in each peak?

	lo_hc_smooth = where(gz2main.(wfind[1]) ge 0.8 and gz2main.pe_s0 le 0.2, locount)	; 17,599 galaxies
	hi_hc_smooth = where(gz2main.(wfind[1]) ge 0.8 and gz2main.pe_s0 ge 0.8, hicount)	; 57,696 galaxies

	; Average weighted fraction for Answer 01

	lo_avg_a01 = mean(gz2main[lo_hc_smooth].(wfind[1]))	; 0.868
	hi_avg_a01 = mean(gz2main[hi_hc_smooth].(wfind[1]))	; 0.886		Votes are slightly higher for the high HC confidence early types; only <2% difference, though

	; Mean number of votes per answer for the two populations

	lo_nvotes_a01 = mean(gz2main[lo_hc_smooth].(wfind[1]-2))	; 36.9
	hi_nvotes_a01 = mean(gz2main[hi_hc_smooth].(wfind[1]-2))	; 38.6		
									; High HC early types have a few more votes, but not significantly different. 
									; Mean for all strong ellipticals is 38.2 votes/galaxy
									; Mean for total GZ2 sample is 27.4 votes/galaxy

	; Are they classified as being odd? If so, what is the category breakdown?

		; Set the 'odd' classification as at least 25% of users clicking the "odd" button

	lo_odd = where(gz2main.(wfind[1]) ge 0.8 and gz2main.pe_s0 le 0.2 and gz2main.(wfind[14]) gt 0.25, lo_oddcount)	; 1,345 / 17,599 =  7.6%
	hi_odd = where(gz2main.(wfind[1]) ge 0.8 and gz2main.pe_s0 ge 0.8 and gz2main.(wfind[14]) gt 0.25, hi_oddcount)	; 7,099 / 57,696 = 12.3%

		; The high-confidence HC early types are somewhat more likely to be classified as odd. Opposite of what I'd expect, given that odd galaxies
		; in the Zoo2 tree might fool the HC algorithm, but be picked up by the users as not genuine ellipticals. 

		; Is the distribution the same if I take out the odd galaxies?
		; 	Yes. Still a significant peak at ~0.1 HC early-type probability. Late-type distribution is also unchanged. 

		; Does it depend on galaxy shape at all?

		lo_roundbulge = where((gz2main[lo_hc_smooth].(wfind[16]) gt gz2main[lo_hc_smooth].(wfind[17])) and (gz2main[lo_hc_smooth].(wfind[16]) gt gz2main[lo_hc_smooth].(wfind[18])), lo_roundcount)
		lo_betweenbulge = where((gz2main[lo_hc_smooth].(wfind[17]) gt gz2main[lo_hc_smooth].(wfind[16])) and (gz2main[lo_hc_smooth].(wfind[17]) gt gz2main[lo_hc_smooth].(wfind[18])), lo_betweencount) 
		lo_cigarbulge = where((gz2main[lo_hc_smooth].(wfind[18]) gt gz2main[lo_hc_smooth].(wfind[17])) and (gz2main[lo_hc_smooth].(wfind[18]) gt gz2main[lo_hc_smooth].(wfind[16])), lo_cigarcount)
		hi_roundbulge = where((gz2main[hi_hc_smooth].(wfind[16]) gt gz2main[hi_hc_smooth].(wfind[17])) and (gz2main[hi_hc_smooth].(wfind[16]) gt gz2main[hi_hc_smooth].(wfind[18])), hi_roundcount) 
		hi_betweenbulge = where((gz2main[hi_hc_smooth].(wfind[17]) gt gz2main[hi_hc_smooth].(wfind[16])) and (gz2main[hi_hc_smooth].(wfind[17]) gt gz2main[hi_hc_smooth].(wfind[18])), hi_betweencount) 
		hi_cigarbulge = where((gz2main[hi_hc_smooth].(wfind[18]) gt gz2main[hi_hc_smooth].(wfind[17])) and (gz2main[hi_hc_smooth].(wfind[18]) gt gz2main[hi_hc_smooth].(wfind[16])), hi_cigarcount) 

		;print, [lo_roundcount,lo_betweencount,lo_cigarcount]/float(locount)
		;print, 1. - total([lo_roundcount,lo_betweencount,lo_cigarcount]/float(locount))
		;print, [hi_roundcount,hi_betweencount,hi_cigarcount]/float(hicount)
		;print, 1. - total([hi_roundcount,hi_betweencount,hi_cigarcount]/float(hicount))

			; HC lo-prob galaxies: 18% round, 58% between,   23% cigar, 0.7% no majority
			; HC hi-prob galaxies: 59% round, 41% between,  0.1% cigar, 0.6% no majority

			; Interesting. Almost none of the high-prob HC early types include the Zoo2 smooth galaxies that are cigar shaped. Assuming
			; some typical tri-axial distribution, there should be a limit on the true inclusion of smooth galaxies; too high of an 
			; axial ratio likely means that it's an inclined disk, but with no visible bulge or features (possibly a lens, via the old
			; Buta & Combes definition). HC might be trained against galaxies with too high of an axial ratio. 

			; Even if limiting the sample to only galaxies that Zoo2 users identified as being round, there's still a small bump at the
			; lower end of the sample. The peak does significantly decrease, going from about 6% to 3% when only round galaxies are considered.

			; So we're still disagreeing on _something_. 

	; Redshift cut?

	lo_redshift = mean(gz2main[lo_hc_smooth].redshift,/nan)		; Mean redshift z=0.068
	hi_redshift = mean(gz2main[hi_hc_smooth].redshift,/nan)		; Mean redshift z=0.100

		; Higher redshifts for the high-confidence galaxies. This could definitely also be a contributing factor; known bias means
		; more Zoo2 votes will go toward smoothness rather than disk. Maybe this will go away with debiasing of the sample, then. 

		; Looking only at low-redshift galaxies INCREASES the low peak. Weird. 
	
	; Colors of the two galaxies

	lo_color = mean(gz2main[lo_hc_smooth].petromag_g - gz2main[lo_hc_smooth].petromag_r,/nan)		; Mean g-r color: 0.67
	hi_color = mean(gz2main[hi_hc_smooth].petromag_g - gz2main[hi_hc_smooth].petromag_r,/nan)		; Mean g-r color: 0.97

		; The low-probability galaxies are significantly BLUER than those in the higher peak, with mean colors approaching that in the green valley. 
		; This is an indication that they might be at least partially composed of late-type galaxies, even if the GZ users disagree with them. 
		; The SVM does use color as one of the metadata parameters in classification: g-r and r-i. They state that this affects ~10% of galaxies classified.

	; Apparent magnitude

	lo_magr = mean(gz2main[lo_hc_smooth].petromag_r,/nan)		; Mean r-mag: 16.7
	hi_magr = mean(gz2main[hi_hc_smooth].petromag_r,/nan)		; Mean r-mag: 16.3	(slightly brighter)

		 ; Magnitude cut of m_r < 16 DOES essentially eliminate the low-HC early types. 

	; What are the T-Types of these galaxies as classified by NA10 or EFIGI?

		; NA10

		na_gz2_hc_file = fitsdir+'na_gz2_hc.fits'
		ngh = mrdfits(na_gz2_hc_file,1,/silent)
		ngh_tagnames=tag_names(ngh)
		nt = n_elements(ngh_tagnames)
		
		temparr=intarr(nt) 
		for i=0,nt-1 do if (strlen(ngh_tagnames[i]) ge 18) and (strmid(ngh_tagnames[i],strlen(ngh_tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 
		
		ngh_wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
		ngh_wcind = [0,where(temparr)]	
		
		nghmain = ngh[where((strtrim(ngh.sample,2) eq 'original' or strtrim(ngh.sample,2) eq 'extra'),maincount)]

		lo_ngh_smooth = where(nghmain.(ngh_wfind[1]) ge 0.8 and nghmain.pe_s0 le 0.2, locount)	; 67 galaxies
		hi_ngh_smooth = where(nghmain.(ngh_wfind[1]) ge 0.8 and nghmain.pe_s0 ge 0.8, hicount)	; 2408 galaxies

		; Locount is low because NA10 catalog only goes down to a depth of 16 mag. Since I showed that the apparent magnitude cut removes almost
		; all galaxies on the spurious tail of the Zoo2 early types, this makes sense. 

		lo_ttype = mode(nghmain[lo_ngh_smooth].ttype)		; 99	uncertain (17/67 galaxies). 7/67 are ellipticals, and 23/67 are various late-type spirals
		hi_ttype = mode(nghmain[hi_ngh_smooth].ttype)		; -5	elliptical (1632/2408). Only 14 are late-type spirals or unknown. Others are S0. 

		; EFIGI

		efigi_gz2_hc_file = fitsdir+'efigi_gz2_hc.fits'
		egh = mrdfits(efigi_gz2_hc_file,1,/silent)
		egh_tagnames=tag_names(egh)
		nt = n_elements(egh_tagnames)
		
		temparr=intarr(nt) 
		for i=0,nt-1 do if (strlen(egh_tagnames[i]) ge 18) and (strmid(egh_tagnames[i],strlen(egh_tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 
		
		egh_wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
		egh_wcind = [0,where(temparr)]	
		
		eghmain = egh[where((strtrim(egh.sample,2) eq 'original' or strtrim(egh.sample,2) eq 'extra'),maincount)]

		lo_egh_smooth = where(eghmain.(egh_wfind[1]) ge 0.8 and eghmain.pe_s0 le 0.2, locount)	; 20 galaxies
		hi_egh_smooth = where(eghmain.(egh_wfind[1]) ge 0.8 and eghmain.pe_s0 ge 0.8, hicount)	; 153 galaxies

		; Very low overlap on the EFIGI/GZ2/HC elliptical samples. 

		lo_ttype = mode(eghmain[lo_egh_smooth].t)		; 16/20 galaxies are dwarf spheroid ellipticals
		hi_ttype = mode(eghmain[hi_egh_smooth].t)		; -5	elliptical (82/153). All 153 galaxies are between -6 and -1 (E and S0). Cool.

	; What are the splits among E and S0 galaxies from the HC probabilities?

		lo_prob_e = gz2main[lo_hc_smooth].pe
		lo_prob_s0 = gz2main[lo_hc_smooth].ps0
		hi_prob_e = gz2main[hi_hc_smooth].pe
		hi_prob_s0 = gz2main[hi_hc_smooth].ps0

		lo_e = where(lo_prob_e gt lo_prob_s0, lell) 	;  0.1%
		lo_s0 = where(lo_prob_s0 gt lo_prob_e, ls0) 	; 99.9%
		hi_e = where(hi_prob_e gt hi_prob_s0, hell) 	; 74%
		hi_s0 = where(hi_prob_s0 gt hi_prob_e, hs0) 	; 26%

		; Amongst the low ones, there's also a much higher likelihood that the low probability ones are S0. That's basically a tautulogy; 
		; anything above 0% is likely to have only votes coming from there. 

	; Plot the overall agreement for Zoo2 vs. HC probabilities

	hcbinsize = 0.05
	gzbinsize = 0.05
	hcbins = fillarr(hcbinsize,0,1.)
	gzbins = fillarr(gzbinsize,0,1.)
	ni = n_elements(hcbins)
	nj = n_elements(gzbins)

	restore, gz2dir+'hc_grids.sav'

	if keyword_set(ps) then begin
		ps_start, filename=figdir+'hc_gz2.eps',/color, /quiet,/encap, xsize=7, ysize=3
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

		;grid_el = fltarr(ni,nj)
		;grid_sp = fltarr(ni,nj)

		;for i = 0,ni-1 do begin & for j=0,nj-1 do begin & grid_el[i,j] = n_elements(where((gz2main.pe_s0 gt hcbins[i]) and (gz2main.pe_s0 le (hcbins[i]+hcbinsize)) and (gz2main.(wfind[1]) gt gzbins[j]) and (gz2main.(wfind[1]) le (gzbins[j]+gzbinsize)))) & endfor & endfor
		;for i = 0,ni-1 do begin & for j=0,nj-1 do begin & grid_sp[i,j] = n_elements(where((1.-gz2main.pe_s0 gt hcbins[i]) and (1.-gz2main.pe_s0 le (hcbins[i]+hcbinsize)) and (gz2main.(wfind[2]) gt gzbins[j]) and (gz2main.(wfind[2]) le (gzbins[j]+gzbinsize)))) & endfor & endfor

		;grid_el_frac = grid_el/max(grid_el)
		;grid_el_frac_log = alog10(grid_el)/max(alog10(grid_el))
		;grid_sp_frac = grid_sp/max(grid_sp)
		;grid_sp_frac_log = alog10(grid_sp)/max(alog10(grid_sp))

		cgloadct,3
		!p.multi=[0,2,1]
		cgplot,indgen(10),/nodata,xr=[0,1],yr=[0,1], xtitle='HC E-S0 probability',ytitle='GZ2 smooth w.f.', charsize = cs
		for i=0,ni-1 do begin & for j=0,nj-1 do begin & cgcolorfill,replicate(hcbins[i],5)+[0,0,1,1,0]*hcbinsize,replicate(gzbins[j],5)+[0,1,1,0,0]*gzbinsize,color=fix(grid_el_frac_log[i,j]*255) & endfor & endfor
		
		cgloadct,1
		cgplot,indgen(10),/nodata,xr=[0,1],yr=[0,1], xtitle='HC Sab-Scd probability',ytitle='GZ2 features/disk w.f.', charsize = cs
		for i=0,ni-1 do begin & for j=0,nj-1 do begin & cgcolorfill,replicate(hcbins[i],5)+[0,0,1,1,0]*hcbinsize,replicate(gzbins[j],5)+[0,1,1,0,0]*gzbinsize,color=fix(grid_sp_frac_log[i,j]*255) & endfor & endfor
	
	if keyword_set(ps) then ps_end

	; Plot same results for subclasses of HC morphology

	if keyword_set(ps) then begin
		ps_start, filename=figdir+'hc_gz2_subclass.eps',/color, /quiet,/encap, xsize=7, ysize=7
		cs=0.8
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

		;grid_pe = fltarr(ni,nj)
		;grid_ps0 = fltarr(ni,nj)
		;grid_psab = fltarr(ni,nj)
		;grid_pscd = fltarr(ni,nj)

		;for i = 0,ni-1 do begin 
		;	for j=0,nj-1 do begin 
		;		grid_pe[i,j] = n_elements(where((1.-gz2main.pe gt hcbins[i]) and (1.-gz2main.pe le (hcbins[i]+hcbinsize)) and (gz2main.(wfind[1]) gt gzbins[j]) and (gz2main.(wfind[1]) le (gzbins[j]+gzbinsize)))) 
		;		grid_ps0[i,j] = n_elements(where((1.-gz2main.ps0 gt hcbins[i]) and (1.-gz2main.ps0 le (hcbins[i]+hcbinsize)) and (gz2main.(wfind[1]) gt gzbins[j]) and (gz2main.(wfind[1]) le (gzbins[j]+gzbinsize)))) 
		;		grid_psab[i,j] = n_elements(where((1.-gz2main.psab gt hcbins[i]) and (1.-gz2main.psab le (hcbins[i]+hcbinsize)) and (gz2main.(wfind[1]) gt gzbins[j]) and (gz2main.(wfind[1]) le (gzbins[j]+gzbinsize)))) 
		;		grid_pscd[i,j] = n_elements(where((1.-gz2main.pscd gt hcbins[i]) and (1.-gz2main.pscd le (hcbins[i]+hcbinsize)) and (gz2main.(wfind[1]) gt gzbins[j]) and (gz2main.(wfind[1]) le (gzbins[j]+gzbinsize)))) 
		;	endfor 
		;endfor

		;grid_pe_frac = grid_pe/max(grid_pe)
		;grid_pe_frac_log = alog10(grid_pe)/max(alog10(grid_pe))
		;grid_ps0_frac = grid_ps0/max(grid_ps0)
		;grid_ps0_frac_log = alog10(grid_ps0)/max(alog10(grid_ps0))
		;grid_psab_frac = grid_psab/max(grid_psab)
		;grid_psab_frac_log = alog10(grid_psab)/max(alog10(grid_psab))
		;grid_pscd_frac = grid_pscd/max(grid_pscd)
		;grid_pscd_frac_log = alog10(grid_pscd)/max(alog10(grid_pscd))

		;save, grid_el_frac, grid_el_frac_log, grid_sp_frac, grid_sp_frac_log, grid_pe_frac, grid_pe_frac_log, grid_ps0_frac, grid_ps0_frac_log, grid_psab_frac, grid_psab_frac_log, grid_pscd_frac, grid_pscd_frac_log, filename=gz2dir+'hc_grids.sav'

		; Find the median line (weight the mean by fraction in bin)

		avg_pe = fltarr(nj)
		avg_ps0 = fltarr(nj)
		avg_psab = fltarr(nj)
		avg_pscd = fltarr(nj)

		for i=0, ni-1 do begin
			avg_pe[i] = total(gzbins * grid_pe_frac_log[i,*]) / total(gzbins)
			avg_ps0[i] = total(gzbins * grid_ps0_frac_log[i,*]) / total(gzbins)
			avg_psab[i] = total(gzbins * grid_psab_frac_log[i,*]) / total(gzbins)
			avg_pscd[i] = total(gzbins * grid_pscd_frac_log[i,*]) / total(gzbins)
		endfor

		!p.multi=[0,2,2]

		cgloadct,0
		cgplot,indgen(10),/nodata,xr=[0,1],yr=[0,1], xtitle='HC P(E)',ytitle='GZ2 smooth w.f.', charsize = cs

		for i=0,ni-1 do begin 
			for j=0,nj-1 do begin 
				cgcolorfill,replicate(hcbins[i],5)+[0,0,1,1,0]*hcbinsize,replicate(gzbins[j],5)+[0,1,1,0,0]*gzbinsize,color=fix(grid_pe_frac_log[i,j]*255)
			endfor
		endfor
		
		cgplot, /overplot, hcbins + hcbinsize, avg_pe, thick = 3, color='black'

		cgplot,indgen(10),/nodata,xr=[0,1],yr=[0,1], xtitle='HC P(S0)',ytitle='GZ2 smooth w.f.', charsize = cs

		for i=0,ni-1 do begin 
			for j=0,nj-1 do begin 
				cgcolorfill,replicate(hcbins[i],5)+[0,0,1,1,0]*hcbinsize,replicate(gzbins[j],5)+[0,1,1,0,0]*gzbinsize,color=fix(grid_ps0_frac_log[i,j]*255)
			endfor
		endfor
		
		cgplot, /overplot, hcbins + hcbinsize, avg_ps0, thick = 3, color='black'

		cgplot,indgen(10),/nodata,xr=[0,1],yr=[0,1], xtitle='HC P(Sab)',ytitle='GZ2 smooth w.f.', charsize = cs

		for i=0,ni-1 do begin 
			for j=0,nj-1 do begin 
				cgcolorfill,replicate(hcbins[i],5)+[0,0,1,1,0]*hcbinsize,replicate(gzbins[j],5)+[0,1,1,0,0]*gzbinsize,color=fix(grid_psab_frac_log[i,j]*255)
			endfor
		endfor
		
		cgplot, /overplot, hcbins + hcbinsize, avg_psab, thick = 3, color='black'

		cgplot,indgen(10),/nodata,xr=[0,1],yr=[0,1], xtitle='HC P(Scd)',ytitle='GZ2 smooth w.f.', charsize = cs

		for i=0,ni-1 do begin 
			for j=0,nj-1 do begin 
				cgcolorfill,replicate(hcbins[i],5)+[0,0,1,1,0]*hcbinsize,replicate(gzbins[j],5)+[0,1,1,0,0]*gzbinsize,color=fix(grid_pscd_frac_log[i,j]*255)
			endfor
		endfor
		
		cgplot, /overplot, hcbins + hcbinsize, avg_pscd, thick = 3, color='black'

	if keyword_set(ps) then ps_end

; Plot the HC probabilities as function of bulge dominance

mintask = 10
minprob = 0.5
taskind = where(gz2main.T02_EDGEON_A05_NO_WEIGHT ge 30,nt)
print,nt

h1 = double(hist_2d(gz2main[taskind].(wfind[10]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))
h2 = double(hist_2d(gz2main[taskind].(wfind[11]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))
h3 = double(hist_2d(gz2main[taskind].(wfind[12]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))
h4 = double(hist_2d(gz2main[taskind].(wfind[13]),     gz2main[taskind].pscd,bin1 = 0.05, bin2 = 0.05, min1 = 0.00, max1 = 1.00, min2 = 0.00, max2 = 1.00))

y = fillarr(0.05, 0.0, 1.05)
x = fillarr(0.05, 0.0, 1.05)
whiteback

if keyword_set(ps) then begin
	ps_start, filename=figdir+'hc_gz2_bulge_contour.eps',/color, /quiet,/encap, xsize=15, ysize=15
	cs=2
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

!p.multi=[0,2,2]

ncolors = 50
cgloadct, 33, ncolors=ncolors, bottom=1
levels = fillarr(0.07,0,3.5)
;levels = fillarr(10,1,500)

cgcontour, alog10(h1+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='No bulge'
cgcontour, alog10(h2+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='Just noticeable'
cgcontour, alog10(h3+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='Obvious bulge'
cgcontour, alog10(h4+1), x, y, /fill, levels=levels, c_colors=indgen(ncolors)+1, ytitle='HC Scd probability', xtitle='GZ2 vote fraction', charsize = cs, title='Dominant bulge'

cgcolorbar,position=[0.45,0.52,0.60,0.54],title='log(N!Igal!N + 1)',ncolors=ncolors

if keyword_set(ps) then ps_end

if keyword_set(stop) then stop

end
