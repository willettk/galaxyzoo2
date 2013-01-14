
;+
; NAME:
;       
;	GZ2_BIAS_DEMO
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

pro gz2_bias_demo_rawlikelihood, ps=ps, volumelimited=volumelimited, plotone = plotone

device,retain=2

timestart = systime(1)

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
figsdir = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'
gz2tablesamplefile = fitsdir+'gz2_table_sample_match_coadd2.fits'

	;gz2 = mrdfits(gz2tablesamplefile, 1, /silent)
	;save,gz2,filename=fitsdir+'gz2.sav'

restore,gz2dir+'gz2.sav'

tagnames=tag_names(gz2)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]	

gz2main = gz2[where((strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra'),maincount)]
coadd2 = gz2[where(strtrim(gz2.sample,2) eq 'stripe82_coadd_2' and (gz2.petromag_r - gz2.extinction_r) le 17.0 and gz2.petror90_r ge 3.,s82c_count)]
normal = gz2[where(strtrim(gz2.sample,2) eq 'stripe82'         and (gz2.petromag_r - gz2.extinction_r) le 17.0 and gz2.petror90_r ge 3.,s82n_count)]

; Volume-limiting

	absr = gz2main.petromag_mr
	z = gz2main.redshift
	
	delz_vol=0.001
	zarr_vol = fillarr(delz_vol, delz_vol, 0.30)
	nz_vol = n_elements(zarr_vol)
	maxmag = fltarr(nz_vol)
	for i=0, nz_vol-1 do begin & maxmag[i] = max(absr[where(z gt zarr_vol[i] and z le zarr_vol[i]+delz_vol)]) & endfor
	
	expr='p[0] + p[1] * alog10(x/p[2])'                                              
	start=[-21., -6., 0.2]
	p = mpfitexpr(expr,zarr_vol, maxmag, maxmag*0.1,start,/quiet)
	
	if keyword_set(ps) then begin
		ps_start, filename=figsdir+'gz2_volumelimited.ps', /color, /quiet
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

	;!p.multi=[0,1,2]
	;cgplot, z, absr, psym=3, xtitle='Redshift', ytitle='M!Ir!N', charsize=2
	;cgplot, /over, coadd2.redshift, coadd2.petromag_mr, psym=3, color='green'
	;cgplot, /over, zarr_vol, p[0] + p[1] * alog10(zarr_vol/p[2]), color='red'
	;
	vollim = lonarr(nz_vol)
	for i=0, nz_vol - 1 do begin & vollim[i] = n_elements(where(z le zarr_vol[i] and absr le p[0] + p[1] * alog10(zarr_vol[i]/p[2]))) & endfor
	;
	;cgplot, zarr_vol, vollim, xtitle='Redshift', ytitle='Size of volume-limited sample', charsize=2
	;cgplot, /over, replicate(zarr_vol[where(vollim eq max(vollim))],2), !y.crange, color='red'
	;cgtext, 0.7, 0.3, string(zarr_vol[where(vollim eq max(vollim))],format='(f6.3)'), /normal, charsize=1.5
	
	if keyword_set(ps) then ps_end

	zlim = (zarr_vol[where(vollim eq max(vollim))])[0]
	absmag_lim = p[0] + p[1] * alog10(zlim/p[2])

	vl = where(gz2main.redshift lt zlim and gz2main.petromag_mr lt absmag_lim)
	vl_s82 = where(coadd2.redshift lt zlim and coadd2.petromag_mr lt absmag_lim)

	if keyword_set(volumelimited) then begin
		gz2main = gz2main[vl]
		coadd2 = coadd2[vl_s82]
		vltitle='_volumelimited'
		zmax = zlim
	endif else begin
		vltitle=''
		;zmax = 0.25
		zmax = 0.20
	endelse

; Two plots: thresholded likelihoods (try f > 0.8), and raw vote summation for all answers as a function of redshift. Repeat as function of surface brightness, physical size, and apparent magnitude

yrange = [0,1.2]

delz=0.02
zarr = fillarr(delz,0,zmax)
nz = n_elements(zarr)
xrange = [-0.02,0.30]
xvar_main = gz2main.redshift
xvar_coadd2 = coadd2.redshift
xvar_normal = normal.redshift
xlabel = 'Redshift'

delr50 = 2
r50arr = fillarr(delr50,0,15)
nr50 = n_elements(r50arr)

delrmag = 0.5
rmagarr = fillarr(delrmag,13,17)
nrmag = n_elements(rmagarr)

	xrange=[min(zarr),max(zarr)*1.1]

if keyword_set(ps) then begin

	if keyword_set(plotone) then begin
		plotname='gz2_bias_demo_rawlikelihood_task01'+vltitle+'.ps' 
		encap=1
	endif
	plotname='gz2_bias_demo_rawlikelihood_alltasks'+vltitle+'.ps'
	;plotname='gz2_bias_demo_alltasks_r50'+vltitle+'.ps'
	;plotname='gz2_bias_demo_alltasks_rmag'+vltitle+'.ps'

	ps_start, filename=figsdir+plotname, /color, /quiet, encap = encap
	cs=1.2
	legendcs = 0.9
	th=3
	thickline=10
	thinline=1
endif else begin
	cs=2
	legendcs = cs
	th=1
	th=3
	thickline=3
	thinline=1
endelse

	; Task 01 - smooth, features, or artifact

	t01_votelim = 10
	
	a01_thresh = fltarr(nz)
	a02_thresh = fltarr(nz)
	a03_thresh = fltarr(nz)
	t01_unc_thresh = fltarr(nz)
	
	a01_wmean = fltarr(nz)
	a02_wmean = fltarr(nz)
	a03_wmean = fltarr(nz)
	
	a01_sumvf = fltarr(nz)
	a02_sumvf = fltarr(nz)
	a03_sumvf = fltarr(nz)
	
		s82c_a01_thresh = fltarr(nz)
		s82c_a02_thresh = fltarr(nz)
		s82c_a03_thresh = fltarr(nz)
		s82c_t01_unc_thresh = fltarr(nz)

		s82c_a01_raw = fltarr(nz)
		s82c_a02_raw = fltarr(nz)
		s82c_a03_raw = fltarr(nz)
	
		s82n_a01_thresh = fltarr(nz)
		s82n_a02_thresh = fltarr(nz)
		s82n_a03_thresh = fltarr(nz)
		s82n_t01_unc_thresh = fltarr(nz)

		s82n_a01_raw = fltarr(nz)
		s82n_a02_raw = fltarr(nz)
		s82n_a03_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz
		a01_thresh_ind = where(gz2main.(wfind[1]) gt 0.8 and gz2main.(wfind[3]+2) ge t01_votelim and xvar_main ge zlo and xvar_main lt zhi, a01_count)
		a02_thresh_ind = where(gz2main.(wfind[2]) gt 0.8 and gz2main.(wfind[3]+2) ge t01_votelim and xvar_main ge zlo and xvar_main lt zhi, a02_count)
		a03_thresh_ind = where(gz2main.(wfind[3]) gt 0.8 and gz2main.(wfind[3]+2) ge t01_votelim and xvar_main ge zlo and xvar_main lt zhi, a03_count)
		t01_ind        = where(                              gz2main.(wfind[3]+2) ge t01_votelim and xvar_main ge zlo and xvar_main lt zhi, t01_count)
	
		a01_thresh[i] = float(a01_count)/float(t01_count)
		a02_thresh[i] = float(a02_count)/float(t01_count)
		a03_thresh[i] = float(a03_count)/float(t01_count)
		t01_unc_thresh[i] = 1. - float(a01_count + a02_count + a03_count)/float(t01_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t01_ind].(wfind[3]+2)
		a01_wmean[i] = total(gz2main[t01_ind].(wfind[1]) * weights) / total(weights)
		a02_wmean[i] = total(gz2main[t01_ind].(wfind[2]) * weights) / total(weights)
		a03_wmean[i] = total(gz2main[t01_ind].(wfind[3]) * weights) / total(weights)

		; Sum the vote fractions

		a01_sumvf[i] = total(gz2main[t01_ind].(wfind[1])) / t01_count
		a02_sumvf[i] = total(gz2main[t01_ind].(wfind[2])) / t01_count
		a03_sumvf[i] = total(gz2main[t01_ind].(wfind[3])) / t01_count

			; Stripe 82 - coadd2

			s82c_a01_thresh_ind = where(coadd2.(wfind[1]) gt 0.8 and coadd2.(wfind[3]+2) ge t01_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a01_count)
			s82c_a02_thresh_ind = where(coadd2.(wfind[2]) gt 0.8 and coadd2.(wfind[3]+2) ge t01_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a02_count)
			s82c_a03_thresh_ind = where(coadd2.(wfind[3]) gt 0.8 and coadd2.(wfind[3]+2) ge t01_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a03_count)
			s82c_t01_ind        = where(                             coadd2.(wfind[3]+2) ge t01_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t01_count)
	
			s82c_a01_thresh[i] = float(s82c_a01_count)/float(s82c_t01_count)
			s82c_a02_thresh[i] = float(s82c_a02_count)/float(s82c_t01_count)
			s82c_a03_thresh[i] = float(s82c_a03_count)/float(s82c_t01_count)
			s82c_t01_unc_thresh[i] = 1. - float(s82c_a01_count + s82c_a02_count + s82c_a03_count)/float(s82c_t01_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82c_weights = coadd2[s82c_t01_ind].(wfind[3]+2)
			s82c_a01_raw[i] = total(coadd2[s82c_t01_ind].(wfind[1]) * s82c_weights) / total(s82c_weights)
			s82c_a02_raw[i] = total(coadd2[s82c_t01_ind].(wfind[2]) * s82c_weights) / total(s82c_weights)
			s82c_a03_raw[i] = total(coadd2[s82c_t01_ind].(wfind[3]) * s82c_weights) / total(s82c_weights)

			; normaldepth

			s82n_a01_thresh_ind = where(normal.(wfind[1]) gt 0.8 and normal.(wfind[3]+2) ge t01_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a01_count)
			s82n_a02_thresh_ind = where(normal.(wfind[2]) gt 0.8 and normal.(wfind[3]+2) ge t01_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a02_count)
			s82n_a03_thresh_ind = where(normal.(wfind[3]) gt 0.8 and normal.(wfind[3]+2) ge t01_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a03_count)
			s82n_t01_ind        = where(                             normal.(wfind[3]+2) ge t01_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t01_count)
	
			s82n_a01_thresh[i] = float(s82n_a01_count)/float(s82n_t01_count)
			s82n_a02_thresh[i] = float(s82n_a02_count)/float(s82n_t01_count)
			s82n_a03_thresh[i] = float(s82n_a03_count)/float(s82n_t01_count)
			s82n_t01_unc_thresh[i] = 1. - float(s82n_a01_count + s82n_a02_count + s82n_a03_count)/float(s82n_t01_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82n_weights = normal[s82n_t01_ind].(wfind[3]+2)
			s82n_a01_raw[i] = total(normal[s82n_t01_ind].(wfind[1]) * s82n_weights) / total(s82n_weights)
			s82n_a02_raw[i] = total(normal[s82n_t01_ind].(wfind[2]) * s82n_weights) / total(s82n_weights)
			s82n_a03_raw[i] = total(normal[s82n_t01_ind].(wfind[3]) * s82n_weights) / total(s82n_weights)
	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a01_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 01', charsize = cs
	cgplot, /over, zarr, a02_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a03_thresh, color='green', thick=thickline
	cgplot, /over, zarr, t01_unc_thresh, color='orange', thick=thickline
	
		;cgplot, /over, zarr, s82c_a01_thresh, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a02_thresh, color='blue', thick=thinline
		;cgplot, /over, zarr, s82c_a03_thresh, color='green', thick=thinline
		;cgplot, /over, zarr, s82c_t01_unc_thresh, color='orange', thick=thinline
	
		;cgplot, /over, zarr, s82n_a01_thresh, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a02_thresh, color='blue', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a03_thresh, color='green', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_t01_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['Smooth','Features','Artifact','Unclassified'],color=['red','blue','green','orange'],psym=28, charsize = legendcs
	
	cgplot,        zarr, a01_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 01', charsize = cs
	cgplot, /over, zarr, a02_wmean, color='blue', thick=thickline
	cgplot, /over, zarr, a03_wmean, color='green', thick=thickline
	
		cgplot, /over, zarr, a01_sumvf, color='red', thick=thinline
		cgplot, /over, zarr, a02_sumvf, color='blue', thick=thinline
		cgplot, /over, zarr, a03_sumvf, color='green', thick=thinline

		;cgplot, /over, zarr, s82c_a01_raw, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a02_raw, color='blue', thick=thinline
		;cgplot, /over, zarr, s82c_a03_raw, color='green', thick=thinline

		;cgplot, /over, zarr, s82n_a01_raw, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a02_raw, color='blue', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a03_raw, color='green', thick=thinline, linestyle=2
	
	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize=legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	if ~keyword_set(plotone) then begin

	; Task 02 - edge-on or not?

	t02_votelim = 10
	
	a04_thresh = fltarr(nz)
	a05_thresh = fltarr(nz)
	t02_unc_thresh = fltarr(nz)
	
	a04_wmean = fltarr(nz)
	a05_wmean = fltarr(nz)
	
	a04_sumvf = fltarr(nz)
	a05_sumvf = fltarr(nz)
	
		s82c_a04_thresh = fltarr(nz)
		s82c_a05_thresh = fltarr(nz)
		s82c_t02_unc_thresh = fltarr(nz)

		s82c_a04_raw = fltarr(nz)
		s82c_a05_raw = fltarr(nz)
	
		s82n_a04_thresh = fltarr(nz)
		s82n_a05_thresh = fltarr(nz)
		s82n_t02_unc_thresh = fltarr(nz)

		s82n_a04_raw = fltarr(nz)
		s82n_a05_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz

		a04_thresh_ind = where(gz2main.(wfind[4]) gt 0.8 and gz2main.(wfind[5]+2) ge t02_votelim and xvar_main ge zlo and xvar_main lt zhi, a04_count)
		a05_thresh_ind = where(gz2main.(wfind[5]) gt 0.8 and gz2main.(wfind[5]+2) ge t02_votelim and xvar_main ge zlo and xvar_main lt zhi, a05_count)
		t02_ind     = where(                              gz2main.(wfind[5]+2) ge t02_votelim and xvar_main ge zlo and xvar_main lt zhi, t02_count)
	
		a04_thresh[i] = float(a04_count)/float(t02_count)
		a05_thresh[i] = float(a05_count)/float(t02_count)
		t02_unc_thresh[i] = 1. - float(a04_count + a05_count)/float(t02_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t02_ind].(wfind[5]+2)
		a04_wmean[i] = total(gz2main[t02_ind].(wfind[4]) * weights) / total(weights)
		a05_wmean[i] = total(gz2main[t02_ind].(wfind[5]) * weights) / total(weights)

		; Sum the vote fractions
	
		a04_sumvf[i] = total(gz2main[t02_ind].(wfind[4])) / t02_count
		a05_sumvf[i] = total(gz2main[t02_ind].(wfind[5])) / t02_count

			; Stripe 82 - coadd2

			s82c_a04_thresh_ind = where(coadd2.(wfind[4]) gt 0.8 and coadd2.(wfind[5]+2) ge t02_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a04_count)
			s82c_a05_thresh_ind = where(coadd2.(wfind[5]) gt 0.8 and coadd2.(wfind[5]+2) ge t02_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a05_count)
			s82c_t02_ind     = where(                              coadd2.(wfind[5]+2) ge t02_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t02_count)
	
			s82c_a04_thresh[i] = float(s82c_a04_count)/float(s82c_t02_count)
			s82c_a05_thresh[i] = float(s82c_a05_count)/float(s82c_t02_count)
			s82c_t02_unc_thresh[i] = 1. - float(s82c_a04_count + s82c_a05_count)/float(s82c_t02_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82c_weights = coadd2[s82c_t02_ind].(wfind[5]+2)
			s82c_a04_raw[i] = total(coadd2[s82c_t02_ind].(wfind[4]) * s82c_weights) / total(s82c_weights)
			s82c_a05_raw[i] = total(coadd2[s82c_t02_ind].(wfind[5]) * s82c_weights) / total(s82c_weights)

			; Stripe 82 - normal

			s82n_a04_thresh_ind = where(normal.(wfind[4]) gt 0.8 and normal.(wfind[5]+2) ge t02_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a04_count)
			s82n_a05_thresh_ind = where(normal.(wfind[5]) gt 0.8 and normal.(wfind[5]+2) ge t02_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a05_count)
			s82n_t02_ind     = where(                              normal.(wfind[5]+2) ge t02_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t02_count)
	
			s82n_a04_thresh[i] = float(s82n_a04_count)/float(s82n_t02_count)
			s82n_a05_thresh[i] = float(s82n_a05_count)/float(s82n_t02_count)
			s82n_t02_unc_thresh[i] = 1. - float(s82n_a04_count + s82n_a05_count)/float(s82n_t02_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82n_weights = normal[s82n_t02_ind].(wfind[5]+2)
			s82n_a04_raw[i] = total(normal[s82n_t02_ind].(wfind[4]) * s82n_weights) / total(s82n_weights)
			s82n_a05_raw[i] = total(normal[s82n_t02_ind].(wfind[5]) * s82n_weights) / total(s82n_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a04_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 02', charsize = cs
	cgplot, /over, zarr, a05_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, t02_unc_thresh, color='orange', thick=thickline
	
		;cgplot, /over, zarr, s82c_a04_thresh, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a05_thresh, color='blue', thick=thinline
		;cgplot, /over, zarr, s82c_t02_unc_thresh, color='orange', thick=thinline
	
		;cgplot, /over, zarr, s82n_a04_thresh, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a05_thresh, color='blue', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_t02_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['Edge-on','Not edge-on','Unclassified'],color=['red','blue','orange'],psym=28, charsize = legendcs
	
	cgplot,        zarr, a04_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 02', charsize = cs
	cgplot, /over, zarr, a05_wmean, color='blue', thick=thickline

		cgplot, /over, zarr, a04_sumvf, color='red', thick=thinline
		cgplot, /over, zarr, a05_sumvf, color='blue', thick=thinline

		;cgplot, /over, zarr, s82c_a04_raw, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a05_raw, color='blue', thick=thinline

		;cgplot, /over, zarr, s82n_a04_raw, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a05_raw, color='blue', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 03 - bar or not?

	t03_votelim = 10
	
	a06_thresh = fltarr(nz)
	a07_thresh = fltarr(nz)
	t03_unc_thresh = fltarr(nz)
	
	a06_wmean = fltarr(nz)
	a07_wmean = fltarr(nz)
	
	a06_sumvf = fltarr(nz)
	a07_sumvf = fltarr(nz)
	
		s82c_a06_thresh = fltarr(nz)
		s82c_a07_thresh = fltarr(nz)
		s82c_t03_unc_thresh = fltarr(nz)

		s82c_a06_raw = fltarr(nz)
		s82c_a07_raw = fltarr(nz)
	
		s82n_a06_thresh = fltarr(nz)
		s82n_a07_thresh = fltarr(nz)
		s82n_t03_unc_thresh = fltarr(nz)

		s82n_a06_raw = fltarr(nz)
		s82n_a07_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz
		a06_thresh_ind = where(gz2main.(wfind[6]) gt 0.8 and gz2main.(wfind[7]+2) ge t03_votelim and xvar_main ge zlo and xvar_main lt zhi, a06_count)
		a07_thresh_ind = where(gz2main.(wfind[7]) gt 0.8 and gz2main.(wfind[7]+2) ge t03_votelim and xvar_main ge zlo and xvar_main lt zhi, a07_count)
		t03_ind        = where(                              gz2main.(wfind[7]+2) ge t03_votelim and xvar_main ge zlo and xvar_main lt zhi, t03_count)
	
		a06_thresh[i] = float(a06_count)/float(t03_count)
		a07_thresh[i] = float(a07_count)/float(t03_count)
		t03_unc_thresh[i] = 1. - float(a06_count + a07_count)/float(t03_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t03_ind].(wfind[7]+2)
		a06_wmean[i] = total(gz2main[t03_ind].(wfind[6]) * weights) / total(weights)
		a07_wmean[i] = total(gz2main[t03_ind].(wfind[7]) * weights) / total(weights)
	
		; Sum the vote fractions
	
		a06_sumvf[i] = total(gz2main[t03_ind].(wfind[6])) / t03_count
		a07_sumvf[i] = total(gz2main[t03_ind].(wfind[7])) / t03_count
	
			; Stripe 82 - coadd2

			s82c_a06_thresh_ind = where(coadd2.(wfind[6]) gt 0.8 and coadd2.(wfind[7]+2) ge t03_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a06_count)
			s82c_a07_thresh_ind = where(coadd2.(wfind[7]) gt 0.8 and coadd2.(wfind[7]+2) ge t03_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a07_count)
			s82c_t03_ind     = where(                                coadd2.(wfind[7]+2) ge t03_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t03_count)
	
			s82c_a06_thresh[i] = float(s82c_a06_count)/float(s82c_t03_count)
			s82c_a07_thresh[i] = float(s82c_a07_count)/float(s82c_t03_count)
			s82c_t03_unc_thresh[i] = 1. - float(s82c_a06_count + s82c_a07_count)/float(s82c_t03_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82c_weights = coadd2[s82c_t03_ind].(wfind[7]+2)
			s82c_a06_raw[i] = total(coadd2[s82c_t03_ind].(wfind[6]) * s82c_weights) / total(s82c_weights)
			s82c_a07_raw[i] = total(coadd2[s82c_t03_ind].(wfind[7]) * s82c_weights) / total(s82c_weights)

			; Stripe 82 - normal

			s82n_a06_thresh_ind = where(normal.(wfind[6]) gt 0.8 and normal.(wfind[7]+2) ge t03_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a06_count)
			s82n_a07_thresh_ind = where(normal.(wfind[7]) gt 0.8 and normal.(wfind[7]+2) ge t03_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a07_count)
			s82n_t03_ind     = where(                                normal.(wfind[7]+2) ge t03_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t03_count)
	
			s82n_a06_thresh[i] = float(s82n_a06_count)/float(s82n_t03_count)
			s82n_a07_thresh[i] = float(s82n_a07_count)/float(s82n_t03_count)
			s82n_t03_unc_thresh[i] = 1. - float(s82n_a06_count + s82n_a07_count)/float(s82n_t03_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82n_weights = normal[s82n_t03_ind].(wfind[7]+2)
			s82n_a06_raw[i] = total(normal[s82n_t03_ind].(wfind[6]) * s82n_weights) / total(s82n_weights)
			s82n_a07_raw[i] = total(normal[s82n_t03_ind].(wfind[7]) * s82n_weights) / total(s82n_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a06_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 03', charsize = cs
	cgplot, /over, zarr, a07_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, t03_unc_thresh, color='orange', thick=thickline
	
		;cgplot, /over, zarr, s82c_a06_thresh, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a07_thresh, color='blue', thick=thinline
		;cgplot, /over, zarr, s82c_t03_unc_thresh, color='orange', thick=thinline
	
		;cgplot, /over, zarr, s82n_a06_thresh, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a07_thresh, color='blue', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_t03_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['Bar','No bar','Unclassified'],color=['red','blue','orange'],psym=28, charsize = legendcs
	
	cgplot,        zarr, a06_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 03', charsize = cs
	cgplot, /over, zarr, a07_wmean, color='blue', thick=thickline

		cgplot, /over, zarr, a06_sumvf, color='red', thick=thinline
		cgplot, /over, zarr, a07_sumvf, color='blue', thick=thinline

		;cgplot, /over, zarr, s82c_a06_raw, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a07_raw, color='blue', thick=thinline

		;cgplot, /over, zarr, s82n_a06_raw, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a07_raw, color='blue', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 04 - spiral or not?

	t04_votelim = 10
	
	a08_thresh = fltarr(nz)
	a09_thresh = fltarr(nz)
	t04_unc_thresh = fltarr(nz)
	
	a08_wmean = fltarr(nz)
	a09_wmean = fltarr(nz)
	
	a08_sumvf = fltarr(nz)
	a09_sumvf = fltarr(nz)
	
		s82c_a08_thresh = fltarr(nz)
		s82c_a09_thresh = fltarr(nz)
		s82c_t04_unc_thresh = fltarr(nz)

		s82c_a08_raw = fltarr(nz)
		s82c_a09_raw = fltarr(nz)
	
		s82n_a08_thresh = fltarr(nz)
		s82n_a09_thresh = fltarr(nz)
		s82n_t04_unc_thresh = fltarr(nz)

		s82n_a08_raw = fltarr(nz)
		s82n_a09_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz
		a08_thresh_ind = where(gz2main.(wfind[8]) gt 0.8 and gz2main.(wfind[9]+2) ge t04_votelim and xvar_main ge zlo and xvar_main lt zhi, a08_count)
		a09_thresh_ind = where(gz2main.(wfind[9]) gt 0.8 and gz2main.(wfind[9]+2) ge t04_votelim and xvar_main ge zlo and xvar_main lt zhi, a09_count)
		t04_ind        = where(                              gz2main.(wfind[9]+2) ge t04_votelim and xvar_main ge zlo and xvar_main lt zhi, t04_count)
	
		a08_thresh[i] = float(a08_count)/float(t04_count)
		a09_thresh[i] = float(a09_count)/float(t04_count)
		t04_unc_thresh[i] = 1. - float(a08_count + a09_count)/float(t04_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t04_ind].(wfind[9]+2)
		a08_wmean[i] = total(gz2main[t04_ind].(wfind[8]) * weights) / total(weights)
		a09_wmean[i] = total(gz2main[t04_ind].(wfind[9]) * weights) / total(weights)
	
		; Sum the vote fractions

		a08_sumvf[i] = total(gz2main[t04_ind].(wfind[8])) / t04_count
		a09_sumvf[i] = total(gz2main[t04_ind].(wfind[9])) / t04_count

			; Stripe 82 - coadd2

			s82c_a08_thresh_ind = where(coadd2.(wfind[8]) gt 0.8 and coadd2.(wfind[9]+2) ge t04_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a08_count)
			s82c_a09_thresh_ind = where(coadd2.(wfind[9]) gt 0.8 and coadd2.(wfind[9]+2) ge t04_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a09_count)
			s82c_t04_ind        = where(                             coadd2.(wfind[9]+2) ge t04_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t04_count)
	
			s82c_a08_thresh[i] = float(s82c_a08_count)/float(s82c_t04_count)
			s82c_a09_thresh[i] = float(s82c_a09_count)/float(s82c_t04_count)
			s82c_t04_unc_thresh[i] = 1. - float(s82c_a08_count + s82c_a09_count)/float(s82c_t04_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82c_weights = coadd2[s82c_t04_ind].(wfind[9]+2)
			s82c_a08_raw[i] = total(coadd2[s82c_t04_ind].(wfind[8]) * s82c_weights) / total(s82c_weights)
			s82c_a09_raw[i] = total(coadd2[s82c_t04_ind].(wfind[9]) * s82c_weights) / total(s82c_weights)

			; Stripe 82 - normal

			s82n_a08_thresh_ind = where(normal.(wfind[8]) gt 0.8 and normal.(wfind[9]+2) ge t04_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a08_count)
			s82n_a09_thresh_ind = where(normal.(wfind[9]) gt 0.8 and normal.(wfind[9]+2) ge t04_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a09_count)
			s82n_t04_ind        = where(                             normal.(wfind[9]+2) ge t04_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t04_count)
	
			s82n_a08_thresh[i] = float(s82n_a08_count)/float(s82n_t04_count)
			s82n_a09_thresh[i] = float(s82n_a09_count)/float(s82n_t04_count)
			s82n_t04_unc_thresh[i] = 1. - float(s82n_a08_count + s82n_a09_count)/float(s82n_t04_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82n_weights = normal[s82n_t04_ind].(wfind[9]+2)
			s82n_a08_raw[i] = total(normal[s82n_t04_ind].(wfind[8]) * s82n_weights) / total(s82n_weights)
			s82n_a09_raw[i] = total(normal[s82n_t04_ind].(wfind[9]) * s82n_weights) / total(s82n_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a08_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 04', charsize = cs
	cgplot, /over, zarr, a09_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, t04_unc_thresh, color='orange', thick=thickline
	
		;cgplot, /over, zarr, s82c_a08_thresh, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a09_thresh, color='blue', thick=thinline
		;cgplot, /over, zarr, s82c_t04_unc_thresh, color='orange', thick=thinline
	
		;cgplot, /over, zarr, s82n_a08_thresh, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a09_thresh, color='blue', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_t04_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['Spiral','No spiral','Unclassified'],color=['red','blue','orange'],psym=28, charsize = legendcs
	
	cgplot,        zarr, a08_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 04', charsize = cs
	cgplot, /over, zarr, a09_wmean, color='blue', thick=thickline

		cgplot, /over, zarr, a08_sumvf, color='red', thick=thinline
		cgplot, /over, zarr, a09_sumvf, color='blue', thick=thinline

		;cgplot, /over, zarr, s82c_a08_raw, color='red', thick=thinline
		;cgplot, /over, zarr, s82c_a09_raw, color='blue', thick=thinline

		;cgplot, /over, zarr, s82n_a08_raw, color='red', thick=thinline, linestyle=2
		;cgplot, /over, zarr, s82n_a09_raw, color='blue', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs
	
	; Task 05 - bulge prominence

	t05_votelim = 10
	
	a10_thresh = fltarr(nz)
	a11_thresh = fltarr(nz)
	a12_thresh = fltarr(nz)
	a13_thresh = fltarr(nz)
	t05_unc_thresh = fltarr(nz)
	
	a10_wmean = fltarr(nz)
	a11_wmean = fltarr(nz)
	a12_wmean = fltarr(nz)
	a13_wmean = fltarr(nz)
	
	a10_sumvf = fltarr(nz)
	a11_sumvf = fltarr(nz)
	a12_sumvf = fltarr(nz)
	a13_sumvf = fltarr(nz)
	
		s82c_a10_thresh = fltarr(nz)
		s82c_a11_thresh = fltarr(nz)
		s82c_a12_thresh = fltarr(nz)
		s82c_a13_thresh = fltarr(nz)
		s82c_t05_unc_thresh = fltarr(nz)

		s82c_a10_raw = fltarr(nz)
		s82c_a11_raw = fltarr(nz)
		s82c_a12_raw = fltarr(nz)
		s82c_a13_raw = fltarr(nz)
	
		s82n_a10_thresh = fltarr(nz)
		s82n_a11_thresh = fltarr(nz)
		s82n_a12_thresh = fltarr(nz)
		s82n_a13_thresh = fltarr(nz)
		s82n_t05_unc_thresh = fltarr(nz)

		s82n_a10_raw = fltarr(nz)
		s82n_a11_raw = fltarr(nz)
		s82n_a12_raw = fltarr(nz)
		s82n_a13_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz
		a10_thresh_ind = where(gz2main.(wfind[10]) gt 0.8 and gz2main.(wfind[13]+2) ge t05_votelim and xvar_main ge zlo and xvar_main lt zhi, a10_count)
		a11_thresh_ind = where(gz2main.(wfind[11]) gt 0.8 and gz2main.(wfind[13]+2) ge t05_votelim and xvar_main ge zlo and xvar_main lt zhi, a11_count)
		a12_thresh_ind = where(gz2main.(wfind[12]) gt 0.8 and gz2main.(wfind[13]+2) ge t05_votelim and xvar_main ge zlo and xvar_main lt zhi, a12_count)
		a13_thresh_ind = where(gz2main.(wfind[13]) gt 0.8 and gz2main.(wfind[13]+2) ge t05_votelim and xvar_main ge zlo and xvar_main lt zhi, a13_count)
		t05_ind        = where(                               gz2main.(wfind[13]+2) ge t05_votelim and xvar_main ge zlo and xvar_main lt zhi, t05_count)
	
		a10_thresh[i] = float(a10_count)/float(t05_count)
		a11_thresh[i] = float(a11_count)/float(t05_count)
		a12_thresh[i] = float(a12_count)/float(t05_count)
		a13_thresh[i] = float(a13_count)/float(t05_count)
		t05_unc_thresh[i] = 1. - float(a10_count + a11_count + a12_count + a13_count)/float(t05_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t05_ind].(wfind[13]+2)
		a10_wmean[i] = total(gz2main[t05_ind].(wfind[10]) * weights) / total(weights)
		a11_wmean[i] = total(gz2main[t05_ind].(wfind[11]) * weights) / total(weights)
		a12_wmean[i] = total(gz2main[t05_ind].(wfind[12]) * weights) / total(weights)
		a13_wmean[i] = total(gz2main[t05_ind].(wfind[13]) * weights) / total(weights)

		; Sum the vote fractions
	
		a10_sumvf[i] = total(gz2main[t05_ind].(wfind[10])) / t05_count
		a11_sumvf[i] = total(gz2main[t05_ind].(wfind[11])) / t05_count
		a12_sumvf[i] = total(gz2main[t05_ind].(wfind[12])) / t05_count
		a13_sumvf[i] = total(gz2main[t05_ind].(wfind[13])) / t05_count

			; Stripe 82 - coadd2

			s82c_a10_thresh_ind = where(coadd2.(wfind[10]) gt 0.8 and coadd2.(wfind[13]+2) ge t05_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a10_count)
			s82c_a11_thresh_ind = where(coadd2.(wfind[11]) gt 0.8 and coadd2.(wfind[13]+2) ge t05_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a11_count)
			s82c_a12_thresh_ind = where(coadd2.(wfind[12]) gt 0.8 and coadd2.(wfind[13]+2) ge t05_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a12_count)
			s82c_a13_thresh_ind = where(coadd2.(wfind[13]) gt 0.8 and coadd2.(wfind[13]+2) ge t05_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a13_count)
			s82c_t05_ind        = where(                              coadd2.(wfind[13]+2) ge t05_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t05_count)
	
			s82c_a10_thresh[i] = float(s82c_a10_count)/float(s82c_t05_count)
			s82c_a11_thresh[i] = float(s82c_a11_count)/float(s82c_t05_count)
			s82c_a12_thresh[i] = float(s82c_a12_count)/float(s82c_t05_count)
			s82c_a13_thresh[i] = float(s82c_a13_count)/float(s82c_t05_count)
			s82c_t05_unc_thresh[i] = 1. - float(s82c_a10_count + s82c_a11_count + s82c_a12_count + s82c_a13_count)/float(s82c_t05_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82c_weights = coadd2[s82c_t05_ind].(wfind[13]+2)
			s82c_a10_raw[i] = total(coadd2[s82c_t05_ind].(wfind[10]) * s82c_weights) / total(s82c_weights)
			s82c_a11_raw[i] = total(coadd2[s82c_t05_ind].(wfind[11]) * s82c_weights) / total(s82c_weights)
			s82c_a12_raw[i] = total(coadd2[s82c_t05_ind].(wfind[12]) * s82c_weights) / total(s82c_weights)
			s82c_a13_raw[i] = total(coadd2[s82c_t05_ind].(wfind[13]) * s82c_weights) / total(s82c_weights)

			; Stripe 82 - normal

			s82n_a10_thresh_ind = where(normal.(wfind[10]) gt 0.8 and normal.(wfind[13]+2) ge t05_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a10_count)
			s82n_a11_thresh_ind = where(normal.(wfind[11]) gt 0.8 and normal.(wfind[13]+2) ge t05_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a11_count)
			s82n_a12_thresh_ind = where(normal.(wfind[12]) gt 0.8 and normal.(wfind[13]+2) ge t05_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a12_count)
			s82n_a13_thresh_ind = where(normal.(wfind[13]) gt 0.8 and normal.(wfind[13]+2) ge t05_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a13_count)
			s82n_t05_ind        = where(                              normal.(wfind[13]+2) ge t05_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t05_count)
	
			s82n_a10_thresh[i] = float(s82n_a10_count)/float(s82n_t05_count)
			s82n_a11_thresh[i] = float(s82n_a11_count)/float(s82n_t05_count)
			s82n_a12_thresh[i] = float(s82n_a12_count)/float(s82n_t05_count)
			s82n_a13_thresh[i] = float(s82n_a13_count)/float(s82n_t05_count)
			s82n_t05_unc_thresh[i] = 1. - float(s82n_a10_count + s82n_a11_count + s82n_a12_count + s82n_a13_count)/float(s82n_t05_count)
	
			; Mean weighted by total number of votes for each galaxy on that task
	
			s82n_weights = normal[s82n_t05_ind].(wfind[13]+2)
			s82n_a10_raw[i] = total(normal[s82n_t05_ind].(wfind[10]) * s82n_weights) / total(s82n_weights)
			s82n_a11_raw[i] = total(normal[s82n_t05_ind].(wfind[11]) * s82n_weights) / total(s82n_weights)
			s82n_a12_raw[i] = total(normal[s82n_t05_ind].(wfind[12]) * s82n_weights) / total(s82n_weights)
			s82n_a13_raw[i] = total(normal[s82n_t05_ind].(wfind[13]) * s82n_weights) / total(s82n_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a10_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 05', charsize = cs
	cgplot, /over, zarr, a11_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a12_thresh, color='green', thick=thickline
	cgplot, /over, zarr, a13_thresh, color='purple', thick=thickline
	cgplot, /over, zarr, t05_unc_thresh, color='orange', thick=thickline
	
	    ;cgplot, /over, zarr, s82c_a10_thresh, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a11_thresh, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a12_thresh, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_a13_thresh, color='purple', thick=thinline
	    ;cgplot, /over, zarr, s82c_t05_unc_thresh, color='orange', thick=thinline
	
	    ;cgplot, /over, zarr, s82n_a10_thresh, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a11_thresh, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a12_thresh, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a13_thresh, color='purple', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_t05_unc_thresh, color='orange', thick=thinline, linestyle=2

	al_legend,/top,/left,['No bulge','Just noticeable','Obvious','Dominant','Unclassified'],color=['red','blue','green','purple','orange'],psym=28, charsize = legendcs
       
	cgplot,        zarr, a10_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 05', charsize = cs
	cgplot, /over, zarr, a11_wmean, color='blue', thick=thickline
	cgplot, /over, zarr, a12_wmean, color='green', thick=thickline
	cgplot, /over, zarr, a13_wmean, color='purple', thick=thickline

	    cgplot, /over, zarr, a10_sumvf, color='red', thick=thinline
	    cgplot, /over, zarr, a11_sumvf, color='blue', thick=thinline
	    cgplot, /over, zarr, a12_sumvf, color='green', thick=thinline
	    cgplot, /over, zarr, a13_sumvf, color='purple', thick=thinline

	    ;cgplot, /over, zarr, s82c_a10_raw, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a11_raw, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a12_raw, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_a13_raw, color='purple', thick=thinline

	    ;cgplot, /over, zarr, s82n_a10_raw, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a11_raw, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a12_raw, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a13_raw, color='purple', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 06 - anything odd?

	t06_votelim = 10
	
	a14_thresh = fltarr(nz)
	a15_thresh = fltarr(nz)
	t06_unc_thresh = fltarr(nz)
	
	a14_wmean = fltarr(nz)
	a15_wmean = fltarr(nz)
	
	a14_sumvf = fltarr(nz)
	a15_sumvf = fltarr(nz)
	
	    s82n_a14_thresh = fltarr(nz)
	    s82n_a15_thresh = fltarr(nz)
	    s82n_t06_unc_thresh = fltarr(nz)

	    s82n_a14_raw = fltarr(nz)
	    s82n_a15_raw = fltarr(nz)
	    
	    s82c_a14_thresh = fltarr(nz)
	    s82c_a15_thresh = fltarr(nz)
	    s82c_t06_unc_thresh = fltarr(nz)

	    s82c_a14_raw = fltarr(nz)
	    s82c_a15_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz

		a14_thresh_ind = where(gz2main.(wfind[14]) gt 0.8 and gz2main.(wfind[15]+2) ge t06_votelim and xvar_main ge zlo and xvar_main lt zhi, a14_count)
		a15_thresh_ind = where(gz2main.(wfind[15]) gt 0.8 and gz2main.(wfind[15]+2) ge t06_votelim and xvar_main ge zlo and xvar_main lt zhi, a15_count)
		t06_ind     = where(                              gz2main.(wfind[15]+2) ge t06_votelim and xvar_main ge zlo and xvar_main lt zhi, t06_count)
	
		a14_thresh[i] = float(a14_count)/float(t06_count)
		a15_thresh[i] = float(a15_count)/float(t06_count)
		t06_unc_thresh[i] = 1. - float(a14_count + a15_count)/float(t06_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t06_ind].(wfind[15]+2)
		a14_wmean[i] = total(gz2main[t06_ind].(wfind[14]) * weights) / total(weights)
		a15_wmean[i] = total(gz2main[t06_ind].(wfind[15]) * weights) / total(weights)

		; Sum the vote fractions
	
		a14_sumvf[i] = total(gz2main[t06_ind].(wfind[14])) / t06_count
		a15_sumvf[i] = total(gz2main[t06_ind].(wfind[15])) / t06_count

		    ; Stripe 82 - normal

		    s82n_a14_thresh_ind = where(normal.(wfind[14]) gt 0.8 and normal.(wfind[15]+2) ge t06_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a14_count)
		    s82n_a15_thresh_ind = where(normal.(wfind[15]) gt 0.8 and normal.(wfind[15]+2) ge t06_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a15_count)
		    s82n_t06_ind        = where(                              normal.(wfind[15]+2) ge t06_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t06_count)
	
		    s82n_a14_thresh[i] = float(s82n_a14_count)/float(s82n_t06_count)
		    s82n_a15_thresh[i] = float(s82n_a15_count)/float(s82n_t06_count)
		    s82n_t06_unc_thresh[i] = 1. - float(s82n_a14_count + s82n_a15_count)/float(s82n_t06_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82n_weights = normal[s82n_t06_ind].(wfind[15]+2)
		    s82n_a14_raw[i] = total(normal[s82n_t06_ind].(wfind[14]) * s82n_weights) / total(s82n_weights)
		    s82n_a15_raw[i] = total(normal[s82n_t06_ind].(wfind[15]) * s82n_weights) / total(s82n_weights)

		    ; Stripe 82 - coadd2

		    s82c_a14_thresh_ind = where(coadd2.(wfind[14]) gt 0.8 and coadd2.(wfind[15]+2) ge t06_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a14_count)
		    s82c_a15_thresh_ind = where(coadd2.(wfind[15]) gt 0.8 and coadd2.(wfind[15]+2) ge t06_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a15_count)
		    s82c_t06_ind        = where(                              coadd2.(wfind[15]+2) ge t06_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t06_count)
	
		    s82c_a14_thresh[i] = float(s82c_a14_count)/float(s82c_t06_count)
		    s82c_a15_thresh[i] = float(s82c_a15_count)/float(s82c_t06_count)
		    s82c_t06_unc_thresh[i] = 1. - float(s82c_a14_count + s82c_a15_count)/float(s82c_t06_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82c_weights = normal[s82c_t06_ind].(wfind[15]+2)
		    s82c_a14_raw[i] = total(normal[s82c_t06_ind].(wfind[14]) * s82c_weights) / total(s82c_weights)
		    s82c_a15_raw[i] = total(normal[s82c_t06_ind].(wfind[15]) * s82c_weights) / total(s82c_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a14_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 06', charsize = cs
	cgplot, /over, zarr, a15_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, t06_unc_thresh, color='orange', thick=thickline
	
	    ;cgplot, /over, zarr, s82c_a14_thresh, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a15_thresh, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_t06_unc_thresh, color='orange', thick=thinline
	
	    ;cgplot, /over, zarr, s82n_a14_thresh, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a15_thresh, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_t06_unc_thresh, color='orange', thick=thinline, linestyle=2
	    
	al_legend,/top,/left,['Something odd','Nothing odd','Unclassified'],color=['red','blue','orange'],psym=28, charsize = legendcs
	
	cgplot,        zarr, a14_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 06', charsize = cs
	cgplot, /over, zarr, a15_wmean, color='blue', thick=thickline

	    cgplot, /over, zarr, a14_sumvf, color='red', thick=thinline
	    cgplot, /over, zarr, a15_sumvf, color='blue', thick=thinline

	    ;cgplot, /over, zarr, s82c_a14_raw, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a15_raw, color='blue', thick=thinline

	    ;cgplot, /over, zarr, s82n_a14_raw, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a15_raw, color='blue', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 07 - roundness of the smooth galaxy

	t07_votelim = 10
	
	a16_thresh = fltarr(nz)
	a17_thresh = fltarr(nz)
	a18_thresh = fltarr(nz)
	t07_unc_thresh = fltarr(nz)
	
	a16_wmean = fltarr(nz)
	a17_wmean = fltarr(nz)
	a18_wmean = fltarr(nz)
	
	a16_sumvf = fltarr(nz)
	a17_sumvf = fltarr(nz)
	a18_sumvf = fltarr(nz)
	
	    s82n_a16_thresh = fltarr(nz)
	    s82n_a17_thresh = fltarr(nz)
	    s82n_a18_thresh = fltarr(nz)
	    s82n_t07_unc_thresh = fltarr(nz)

	    s82n_a16_raw = fltarr(nz)
	    s82n_a17_raw = fltarr(nz)
	    s82n_a18_raw = fltarr(nz)
	    
	    s82c_a16_thresh = fltarr(nz)
	    s82c_a17_thresh = fltarr(nz)
	    s82c_a18_thresh = fltarr(nz)
	    s82c_t07_unc_thresh = fltarr(nz)

	    s82c_a16_raw = fltarr(nz)
	    s82c_a17_raw = fltarr(nz)
	    s82c_a18_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz

		; Main sample

		a16_thresh_ind = where(gz2main.(wfind[16]) gt 0.8 and gz2main.(wfind[18]+2) ge t07_votelim and xvar_main ge zlo and xvar_main lt zhi, a16_count)
		a17_thresh_ind = where(gz2main.(wfind[17]) gt 0.8 and gz2main.(wfind[18]+2) ge t07_votelim and xvar_main ge zlo and xvar_main lt zhi, a17_count)
		a18_thresh_ind = where(gz2main.(wfind[18]) gt 0.8 and gz2main.(wfind[18]+2) ge t07_votelim and xvar_main ge zlo and xvar_main lt zhi, a18_count)
		t07_ind     = where(                              gz2main.(wfind[18]+2) ge t07_votelim and xvar_main ge zlo and xvar_main lt zhi, t07_count)
	
		a16_thresh[i] = float(a16_count)/float(t07_count)
		a17_thresh[i] = float(a17_count)/float(t07_count)
		a18_thresh[i] = float(a18_count)/float(t07_count)
		t07_unc_thresh[i] = 1. - float(a16_count + a17_count + a18_count)/float(t07_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t07_ind].(wfind[18]+2)
		a16_wmean[i] = total(gz2main[t07_ind].(wfind[16]) * weights) / total(weights)
		a17_wmean[i] = total(gz2main[t07_ind].(wfind[17]) * weights) / total(weights)
		a18_wmean[i] = total(gz2main[t07_ind].(wfind[18]) * weights) / total(weights)

		; Sum the vote fractions
	
		a16_sumvf[i] = total(gz2main[t07_ind].(wfind[16])) / t07_count
		a17_sumvf[i] = total(gz2main[t07_ind].(wfind[17])) / t07_count
		a18_sumvf[i] = total(gz2main[t07_ind].(wfind[18])) / t07_count

		    ; Stripe 82 - normal

		    s82n_a16_thresh_ind = where(normal.(wfind[16]) gt 0.8 and normal.(wfind[18]+2) ge t07_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a16_count)
		    s82n_a17_thresh_ind = where(normal.(wfind[17]) gt 0.8 and normal.(wfind[18]+2) ge t07_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a17_count)
		    s82n_a18_thresh_ind = where(normal.(wfind[18]) gt 0.8 and normal.(wfind[18]+2) ge t07_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a18_count)
		    s82n_t07_ind        = where(                              normal.(wfind[18]+2) ge t07_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t07_count)
	
		    s82n_a16_thresh[i] = float(s82n_a16_count)/float(s82n_t07_count)
		    s82n_a17_thresh[i] = float(s82n_a17_count)/float(s82n_t07_count)
		    s82n_a18_thresh[i] = float(s82n_a18_count)/float(s82n_t07_count)
		    s82n_t07_unc_thresh[i] = 1. - float(s82n_a16_count + s82n_a17_count + s82n_a18_count)/float(s82n_t07_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82n_weights = normal[s82n_t07_ind].(wfind[18]+2)
		    s82n_a16_raw[i] = total(normal[s82n_t07_ind].(wfind[16]) * s82n_weights) / total(s82n_weights)
		    s82n_a17_raw[i] = total(normal[s82n_t07_ind].(wfind[17]) * s82n_weights) / total(s82n_weights)
		    s82n_a18_raw[i] = total(normal[s82n_t07_ind].(wfind[18]) * s82n_weights) / total(s82n_weights)

		    ; Stripe 82 - coadd2

		    s82c_a16_thresh_ind = where(coadd2.(wfind[16]) gt 0.8 and coadd2.(wfind[18]+2) ge t07_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a16_count)
		    s82c_a17_thresh_ind = where(coadd2.(wfind[17]) gt 0.8 and coadd2.(wfind[18]+2) ge t07_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a17_count)
		    s82c_a18_thresh_ind = where(coadd2.(wfind[18]) gt 0.8 and coadd2.(wfind[18]+2) ge t07_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a18_count)
		    s82c_t07_ind        = where(                              coadd2.(wfind[18]+2) ge t07_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t07_count)
	
		    s82c_a16_thresh[i] = float(s82c_a16_count)/float(s82c_t07_count)
		    s82c_a17_thresh[i] = float(s82c_a17_count)/float(s82c_t07_count)
		    s82c_a18_thresh[i] = float(s82c_a18_count)/float(s82c_t07_count)
		    s82c_t07_unc_thresh[i] = 1. - float(s82c_a16_count + s82c_a17_count + s82c_a18_count)/float(s82c_t07_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82c_weights = coadd2[s82c_t07_ind].(wfind[18]+2)
		    s82c_a16_raw[i] = total(coadd2[s82c_t07_ind].(wfind[16]) * s82c_weights) / total(s82c_weights)
		    s82c_a17_raw[i] = total(coadd2[s82c_t07_ind].(wfind[17]) * s82c_weights) / total(s82c_weights)
		    s82c_a18_raw[i] = total(coadd2[s82c_t07_ind].(wfind[18]) * s82c_weights) / total(s82c_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a16_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 07', charsize = cs
	cgplot, /over, zarr, a17_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a18_thresh, color='green', thick=thickline
	cgplot, /over, zarr, t07_unc_thresh, color='orange', thick=thickline
	
	    ;cgplot, /over, zarr, s82c_a16_thresh, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a17_thresh, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a18_thresh, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_t07_unc_thresh, color='orange', thick=thinline
	
	    ;cgplot, /over, zarr, s82n_a16_thresh, color='red', thick=thinline , linestyle=2
	    ;cgplot, /over, zarr, s82n_a17_thresh, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a18_thresh, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_t07_unc_thresh, color='orange', thick=thinline, linestyle=2
	    
	al_legend,/top,/left,['Round galaxy','In between','Cigar-shaped','Unclassified'],color=['red','blue','green','orange'],psym=28, charsize = legendcs
	
	cgplot,        zarr, a16_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 07', charsize = cs
	cgplot, /over, zarr, a17_wmean, color='blue', thick=thickline
	cgplot, /over, zarr, a18_wmean, color='green', thick=thickline

	    cgplot, /over, zarr, a16_sumvf, color='red', thick=thinline
	    cgplot, /over, zarr, a17_sumvf, color='blue', thick=thinline
	    cgplot, /over, zarr, a18_sumvf, color='green', thick=thinline

	    ;cgplot, /over, zarr, s82c_a16_raw, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a17_raw, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a18_raw, color='green', thick=thinline

	    ;cgplot, /over, zarr, s82n_a16_raw, color='red', thick=thinline,  linestyle=2
	    ;cgplot, /over, zarr, s82n_a17_raw, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a18_raw, color='green', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 08 - odd features

	t08_votelim = 10
	
	a19_thresh = fltarr(nz)
	a20_thresh = fltarr(nz)
	a21_thresh = fltarr(nz)
	a22_thresh = fltarr(nz)
	a23_thresh = fltarr(nz)
	a24_thresh = fltarr(nz)
	a38_thresh = fltarr(nz)
	t08_unc_thresh = fltarr(nz)
	
	a19_wmean = fltarr(nz)
	a20_wmean = fltarr(nz)
	a21_wmean = fltarr(nz)
	a22_wmean = fltarr(nz)
	a23_wmean = fltarr(nz)
	a24_wmean = fltarr(nz)
	a38_wmean = fltarr(nz)
	
	a19_sumvf = fltarr(nz)
	a20_sumvf = fltarr(nz)
	a21_sumvf = fltarr(nz)
	a22_sumvf = fltarr(nz)
	a23_sumvf = fltarr(nz)
	a24_sumvf = fltarr(nz)
	a38_sumvf = fltarr(nz)
	
	    s82n_a19_thresh = fltarr(nz)
	    s82n_a20_thresh = fltarr(nz)
	    s82n_a21_thresh = fltarr(nz)
	    s82n_a22_thresh = fltarr(nz)
	    s82n_a23_thresh = fltarr(nz)
	    s82n_a24_thresh = fltarr(nz)
	    s82n_a38_thresh = fltarr(nz)
	    s82n_t08_unc_thresh = fltarr(nz)
	    s82n_a19_raw = fltarr(nz)
	    s82n_a20_raw = fltarr(nz)
	    s82n_a21_raw = fltarr(nz)
	    s82n_a22_raw = fltarr(nz)
	    s82n_a23_raw = fltarr(nz)
	    s82n_a24_raw = fltarr(nz)
	    s82n_a38_raw = fltarr(nz)
	    
	    s82c_a19_thresh = fltarr(nz)
	    s82c_a20_thresh = fltarr(nz)
	    s82c_a21_thresh = fltarr(nz)
	    s82c_a22_thresh = fltarr(nz)
	    s82c_a23_thresh = fltarr(nz)
	    s82c_a24_thresh = fltarr(nz)
	    s82c_a38_thresh = fltarr(nz)
	    s82c_t08_unc_thresh = fltarr(nz)
	    s82c_a19_raw = fltarr(nz)
	    s82c_a20_raw = fltarr(nz)
	    s82c_a21_raw = fltarr(nz)
	    s82c_a22_raw = fltarr(nz)
	    s82c_a23_raw = fltarr(nz)
	    s82c_a24_raw = fltarr(nz)
	    s82c_a38_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz

		; Main sample

		a19_thresh_ind = where(gz2main.(wfind[19]) gt 0.8 and gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, a19_count)
		a20_thresh_ind = where(gz2main.(wfind[20]) gt 0.8 and gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, a20_count)
		a21_thresh_ind = where(gz2main.(wfind[21]) gt 0.8 and gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, a21_count)
		a22_thresh_ind = where(gz2main.(wfind[22]) gt 0.8 and gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, a22_count)
		a23_thresh_ind = where(gz2main.(wfind[23]) gt 0.8 and gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, a23_count)
		a24_thresh_ind = where(gz2main.(wfind[24]) gt 0.8 and gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, a24_count)
		a38_thresh_ind = where(gz2main.(wfind[25]) gt 0.8 and gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, a38_count)
		t08_ind        = where(                               gz2main.(wfind[25]+2) ge t08_votelim and xvar_main ge zlo and xvar_main lt zhi, t08_count)
	
		a19_thresh[i] = float(a19_count)/float(t08_count)
		a20_thresh[i] = float(a20_count)/float(t08_count)
		a21_thresh[i] = float(a21_count)/float(t08_count)
		a22_thresh[i] = float(a22_count)/float(t08_count)
		a23_thresh[i] = float(a23_count)/float(t08_count)
		a24_thresh[i] = float(a24_count)/float(t08_count)
		a38_thresh[i] = float(a38_count)/float(t08_count)
		t08_unc_thresh[i] = 1. - float(a19_count + a20_count + a21_count + a22_count + a23_count + a24_count + a38_count)/float(t08_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t08_ind].(wfind[25]+2)
		a19_wmean[i] = total(gz2main[t08_ind].(wfind[19]) * weights) / total(weights)
		a20_wmean[i] = total(gz2main[t08_ind].(wfind[20]) * weights) / total(weights)
		a21_wmean[i] = total(gz2main[t08_ind].(wfind[21]) * weights) / total(weights)
		a22_wmean[i] = total(gz2main[t08_ind].(wfind[22]) * weights) / total(weights)
		a23_wmean[i] = total(gz2main[t08_ind].(wfind[23]) * weights) / total(weights)
		a24_wmean[i] = total(gz2main[t08_ind].(wfind[24]) * weights) / total(weights)
		a38_wmean[i] = total(gz2main[t08_ind].(wfind[25]) * weights) / total(weights)

		; Sum the vote fractions
	
		a19_sumvf[i] = total(gz2main[t08_ind].(wfind[19])) / t08_count
		a20_sumvf[i] = total(gz2main[t08_ind].(wfind[20])) / t08_count
		a21_sumvf[i] = total(gz2main[t08_ind].(wfind[21])) / t08_count
		a22_sumvf[i] = total(gz2main[t08_ind].(wfind[22])) / t08_count
		a23_sumvf[i] = total(gz2main[t08_ind].(wfind[23])) / t08_count
		a24_sumvf[i] = total(gz2main[t08_ind].(wfind[24])) / t08_count
		a38_sumvf[i] = total(gz2main[t08_ind].(wfind[25])) / t08_count

		    ; Stripe 82 - normal

		    s82n_a19_thresh_ind = where(normal.(wfind[19]) gt 0.8 and normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a19_count)
		    s82n_a20_thresh_ind = where(normal.(wfind[20]) gt 0.8 and normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a20_count)
		    s82n_a21_thresh_ind = where(normal.(wfind[21]) gt 0.8 and normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a21_count)
		    s82n_a22_thresh_ind = where(normal.(wfind[22]) gt 0.8 and normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a22_count)
		    s82n_a23_thresh_ind = where(normal.(wfind[23]) gt 0.8 and normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a23_count)
		    s82n_a24_thresh_ind = where(normal.(wfind[24]) gt 0.8 and normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a24_count)
		    s82n_a38_thresh_ind = where(normal.(wfind[25]) gt 0.8 and normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a38_count)
		    s82n_t08_ind        = where(                              normal.(wfind[25]+2) ge t08_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t08_count)
	
		    s82n_a19_thresh[i] = float(s82n_a19_count)/float(s82n_t08_count)
		    s82n_a20_thresh[i] = float(s82n_a20_count)/float(s82n_t08_count)
		    s82n_a21_thresh[i] = float(s82n_a21_count)/float(s82n_t08_count)
		    s82n_a22_thresh[i] = float(s82n_a22_count)/float(s82n_t08_count)
		    s82n_a23_thresh[i] = float(s82n_a23_count)/float(s82n_t08_count)
		    s82n_a24_thresh[i] = float(s82n_a24_count)/float(s82n_t08_count)
		    s82n_a38_thresh[i] = float(s82n_a38_count)/float(s82n_t08_count)
		    s82n_t08_unc_thresh[i] = 1. - float(s82n_a19_count + s82n_a20_count + s82n_a21_count + s82n_a22_count + s82n_a23_count + s82n_a24_count + s82n_a38_count)/float(s82n_t08_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82n_weights = normal[s82n_t08_ind].(wfind[25]+2)
		    s82n_a19_raw[i] = total(normal[s82n_t08_ind].(wfind[19]) * s82n_weights) / total(s82n_weights)
		    s82n_a20_raw[i] = total(normal[s82n_t08_ind].(wfind[20]) * s82n_weights) / total(s82n_weights)
		    s82n_a21_raw[i] = total(normal[s82n_t08_ind].(wfind[21]) * s82n_weights) / total(s82n_weights)
		    s82n_a22_raw[i] = total(normal[s82n_t08_ind].(wfind[22]) * s82n_weights) / total(s82n_weights)
		    s82n_a23_raw[i] = total(normal[s82n_t08_ind].(wfind[23]) * s82n_weights) / total(s82n_weights)
		    s82n_a24_raw[i] = total(normal[s82n_t08_ind].(wfind[24]) * s82n_weights) / total(s82n_weights)
		    s82n_a38_raw[i] = total(normal[s82n_t08_ind].(wfind[25]) * s82n_weights) / total(s82n_weights)

		    ; Stripe 82 - coadd

		    s82c_a19_thresh_ind = where(coadd2.(wfind[19]) gt 0.8 and coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a19_count)
		    s82c_a20_thresh_ind = where(coadd2.(wfind[20]) gt 0.8 and coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a20_count)
		    s82c_a21_thresh_ind = where(coadd2.(wfind[21]) gt 0.8 and coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a21_count)
		    s82c_a22_thresh_ind = where(coadd2.(wfind[22]) gt 0.8 and coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a22_count)
		    s82c_a23_thresh_ind = where(coadd2.(wfind[23]) gt 0.8 and coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a23_count)
		    s82c_a24_thresh_ind = where(coadd2.(wfind[24]) gt 0.8 and coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a24_count)
		    s82c_a38_thresh_ind = where(coadd2.(wfind[25]) gt 0.8 and coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a38_count)
		    s82c_t08_ind        = where(                              coadd2.(wfind[25]+2) ge t08_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t08_count)
	
		    s82c_a19_thresh[i] = float(s82c_a19_count)/float(s82c_t08_count)
		    s82c_a20_thresh[i] = float(s82c_a20_count)/float(s82c_t08_count)
		    s82c_a21_thresh[i] = float(s82c_a21_count)/float(s82c_t08_count)
		    s82c_a22_thresh[i] = float(s82c_a22_count)/float(s82c_t08_count)
		    s82c_a23_thresh[i] = float(s82c_a23_count)/float(s82c_t08_count)
		    s82c_a24_thresh[i] = float(s82c_a24_count)/float(s82c_t08_count)
		    s82c_a38_thresh[i] = float(s82c_a38_count)/float(s82c_t08_count)
		    s82c_t08_unc_thresh[i] = 1. - float(s82c_a19_count + s82c_a20_count + s82c_a21_count + s82c_a22_count + s82c_a23_count + s82c_a24_count + s82c_a38_count)/float(s82c_t08_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82c_weights = coadd2[s82c_t08_ind].(wfind[25]+2)
		    s82c_a19_raw[i] = total(coadd2[s82c_t08_ind].(wfind[19]) * s82c_weights) / total(s82c_weights)
		    s82c_a20_raw[i] = total(coadd2[s82c_t08_ind].(wfind[20]) * s82c_weights) / total(s82c_weights)
		    s82c_a21_raw[i] = total(coadd2[s82c_t08_ind].(wfind[21]) * s82c_weights) / total(s82c_weights)
		    s82c_a22_raw[i] = total(coadd2[s82c_t08_ind].(wfind[22]) * s82c_weights) / total(s82c_weights)
		    s82c_a23_raw[i] = total(coadd2[s82c_t08_ind].(wfind[23]) * s82c_weights) / total(s82c_weights)
		    s82c_a24_raw[i] = total(coadd2[s82c_t08_ind].(wfind[24]) * s82c_weights) / total(s82c_weights)
		    s82c_a38_raw[i] = total(coadd2[s82c_t08_ind].(wfind[25]) * s82c_weights) / total(s82c_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a19_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 08', charsize = cs
	cgplot, /over, zarr, a20_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a21_thresh, color='green', thick=thickline
	cgplot, /over, zarr, a22_thresh, color='purple', thick=thickline
	cgplot, /over, zarr, a23_thresh, color='yellow', thick=thickline
	cgplot, /over, zarr, a24_thresh, color='grey', thick=thickline
	cgplot, /over, zarr, a38_thresh, color='brown', thick=thickline
	cgplot, /over, zarr, t08_unc_thresh, color='orange', thick=thickline
	
	    ;cgplot, /over, zarr, s82c_a19_thresh, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a20_thresh, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a21_thresh, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_a22_thresh, color='purple', thick=thinline
	    ;cgplot, /over, zarr, s82c_a23_thresh, color='yellow', thick=thinline
	    ;cgplot, /over, zarr, s82c_a24_thresh, color='grey', thick=thinline
	    ;cgplot, /over, zarr, s82c_a38_thresh, color='brown', thick=thinline
	    ;cgplot, /over, zarr, s82c_t08_unc_thresh, color='orange', thick=thinline
	    ;
	    ;cgplot, /over, zarr, s82n_a19_thresh, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a20_thresh, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a21_thresh, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a22_thresh, color='purple', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a23_thresh, color='yellow', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a24_thresh, color='grey', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a38_thresh, color='brown', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_t08_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['Ring','Lens or arc','Disturbed','Irregular','Other','Merger','Dust lane','Unclassified'],color=['red','blue','green','purple','yellow','grey','brown','orange'],psym=28, charsize = legendcs

	cgplot,        zarr, a19_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 08', charsize = cs
	cgplot, /over, zarr, a20_wmean, color='blue', thick=thickline
	cgplot, /over, zarr, a21_wmean, color='green', thick=thickline
	cgplot, /over, zarr, a22_wmean, color='purple', thick=thickline
	cgplot, /over, zarr, a23_wmean, color='yellow', thick=thickline
	cgplot, /over, zarr, a24_wmean, color='grey', thick=thickline
	cgplot, /over, zarr, a38_wmean, color='brown', thick=thickline

	    cgplot, /over, zarr, a19_sumvf, color='red', thick=thinline
	    cgplot, /over, zarr, a20_sumvf, color='blue', thick=thinline
	    cgplot, /over, zarr, a21_sumvf, color='green', thick=thinline
	    cgplot, /over, zarr, a22_sumvf, color='purple', thick=thinline
	    cgplot, /over, zarr, a23_sumvf, color='yellow', thick=thinline
	    cgplot, /over, zarr, a24_sumvf, color='grey', thick=thinline
	    cgplot, /over, zarr, a38_sumvf, color='brown', thick=thinline

	    ;cgplot, /over, zarr, s82c_a19_raw, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a20_raw, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a21_raw, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_a22_raw, color='purple', thick=thinline
	    ;cgplot, /over, zarr, s82c_a23_raw, color='yellow', thick=thinline
	    ;cgplot, /over, zarr, s82c_a24_raw, color='grey', thick=thinline
	    ;cgplot, /over, zarr, s82c_a38_raw, color='brown', thick=thinline

	    ;cgplot, /over, zarr, s82n_a19_raw, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a20_raw, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a21_raw, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a22_raw, color='purple', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a23_raw, color='yellow', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a24_raw, color='grey', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a38_raw, color='brown', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 09 - bulge shape

	t09_votelim = 10
	
	a25_thresh = fltarr(nz)
	a26_thresh = fltarr(nz)
	a27_thresh = fltarr(nz)
	t09_unc_thresh = fltarr(nz)

	a25_wmean = fltarr(nz)
	a26_wmean = fltarr(nz)
	a27_wmean = fltarr(nz)
	
	a25_sumvf = fltarr(nz)
	a26_sumvf = fltarr(nz)
	a27_sumvf = fltarr(nz)
	
	    s82n_a25_thresh = fltarr(nz)
	    s82n_a26_thresh = fltarr(nz)
	    s82n_a27_thresh = fltarr(nz)
	    s82n_t09_unc_thresh = fltarr(nz)
	    s82n_a25_raw = fltarr(nz)
	    s82n_a26_raw = fltarr(nz)
	    s82n_a27_raw = fltarr(nz)
	    
	    s82c_a25_thresh = fltarr(nz)
	    s82c_a26_thresh = fltarr(nz)
	    s82c_a27_thresh = fltarr(nz)
	    s82c_t09_unc_thresh = fltarr(nz)
	    s82c_a25_raw = fltarr(nz)
	    s82c_a26_raw = fltarr(nz)
	    s82c_a27_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz

		; Main sample

		a25_thresh_ind = where(gz2main.(wfind[26]) gt 0.8 and gz2main.(wfind[28]+2) ge t09_votelim and xvar_main ge zlo and xvar_main lt zhi, a25_count)
		a26_thresh_ind = where(gz2main.(wfind[27]) gt 0.8 and gz2main.(wfind[28]+2) ge t09_votelim and xvar_main ge zlo and xvar_main lt zhi, a26_count)
		a27_thresh_ind = where(gz2main.(wfind[28]) gt 0.8 and gz2main.(wfind[28]+2) ge t09_votelim and xvar_main ge zlo and xvar_main lt zhi, a27_count)
		t09_ind     = where(                              gz2main.(wfind[28]+2) ge t09_votelim and xvar_main ge zlo and xvar_main lt zhi, t09_count)
	
		a25_thresh[i] = float(a25_count)/float(t09_count)
		a26_thresh[i] = float(a26_count)/float(t09_count)
		a27_thresh[i] = float(a27_count)/float(t09_count)
		t09_unc_thresh[i] = 1. - float(a25_count + a26_count + a27_count)/float(t09_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t09_ind].(wfind[28]+2)
		a25_wmean[i] = total(gz2main[t09_ind].(wfind[26]) * weights) / total(weights)
		a26_wmean[i] = total(gz2main[t09_ind].(wfind[27]) * weights) / total(weights)
		a27_wmean[i] = total(gz2main[t09_ind].(wfind[28]) * weights) / total(weights)

		; Sum the vote fractions
	
		a25_sumvf[i] = total(gz2main[t09_ind].(wfind[26])) / t09_count
		a26_sumvf[i] = total(gz2main[t09_ind].(wfind[27])) / t09_count
		a27_sumvf[i] = total(gz2main[t09_ind].(wfind[28])) / t09_count

		    ; Stripe 82 - normal

		    s82n_a25_thresh_ind = where(normal.(wfind[26]) gt 0.8 and normal.(wfind[28]+2) ge t09_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a25_count)
		    s82n_a26_thresh_ind = where(normal.(wfind[27]) gt 0.8 and normal.(wfind[28]+2) ge t09_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a26_count)
		    s82n_a27_thresh_ind = where(normal.(wfind[28]) gt 0.8 and normal.(wfind[28]+2) ge t09_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a27_count)
		    s82n_t09_ind        = where(                              normal.(wfind[28]+2) ge t09_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t09_count)
	
		    s82n_a25_thresh[i] = float(s82n_a25_count)/float(s82n_t09_count)
		    s82n_a26_thresh[i] = float(s82n_a26_count)/float(s82n_t09_count)
		    s82n_a27_thresh[i] = float(s82n_a27_count)/float(s82n_t09_count)
		    s82n_t09_unc_thresh[i] = 1. - float(s82n_a25_count + s82n_a26_count + s82n_a27_count)/float(s82n_t09_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82n_weights = normal[s82n_t09_ind].(wfind[28]+2)
		    s82n_a25_raw[i] = total(normal[s82n_t09_ind].(wfind[26]) * s82n_weights) / total(s82n_weights)
		    s82n_a26_raw[i] = total(normal[s82n_t09_ind].(wfind[27]) * s82n_weights) / total(s82n_weights)
		    s82n_a27_raw[i] = total(normal[s82n_t09_ind].(wfind[28]) * s82n_weights) / total(s82n_weights)

		    ; Stripe 82 - coadd2

		    s82c_a25_thresh_ind = where(coadd2.(wfind[26]) gt 0.8 and coadd2.(wfind[28]+2) ge t09_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a25_count)
		    s82c_a26_thresh_ind = where(coadd2.(wfind[27]) gt 0.8 and coadd2.(wfind[28]+2) ge t09_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a26_count)
		    s82c_a27_thresh_ind = where(coadd2.(wfind[28]) gt 0.8 and coadd2.(wfind[28]+2) ge t09_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a27_count)
		    s82c_t09_ind        = where(                              coadd2.(wfind[28]+2) ge t09_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t09_count)
	
		    s82c_a25_thresh[i] = float(s82c_a25_count)/float(s82c_t09_count)
		    s82c_a26_thresh[i] = float(s82c_a26_count)/float(s82c_t09_count)
		    s82c_a27_thresh[i] = float(s82c_a27_count)/float(s82c_t09_count)
		    s82c_t09_unc_thresh[i] = 1. - float(s82c_a25_count + s82c_a26_count + s82c_a27_count)/float(s82c_t09_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82c_weights = coadd2[s82c_t09_ind].(wfind[28]+2)
		    s82c_a25_raw[i] = total(coadd2[s82c_t09_ind].(wfind[26]) * s82c_weights) / total(s82c_weights)
		    s82c_a26_raw[i] = total(coadd2[s82c_t09_ind].(wfind[27]) * s82c_weights) / total(s82c_weights)
		    s82c_a27_raw[i] = total(coadd2[s82c_t09_ind].(wfind[28]) * s82c_weights) / total(s82c_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a25_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 09', charsize = cs
	cgplot, /over, zarr, a26_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a27_thresh, color='green', thick=thickline
	cgplot, /over, zarr, t09_unc_thresh, color='orange', thick=thickline
	
	    ;cgplot, /over, zarr, s82c_a25_thresh, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a26_thresh, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a27_thresh, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_t09_unc_thresh, color='orange', thick=thinline
	    ;
	    ;cgplot, /over, zarr, s82n_a25_thresh, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a26_thresh, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a27_thresh, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_t09_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['Rounded bulge','Boxy bulge','No bulge','Unclassified'],color=['red','blue','green','orange'],psym=28, charsize = legendcs

	cgplot,        zarr, a25_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 09', charsize = cs
	cgplot, /over, zarr, a26_wmean, color='blue', thick=thickline
	cgplot, /over, zarr, a27_wmean, color='green', thick=thickline

	    cgplot, /over, zarr, a25_sumvf, color='red', thick=thinline
	    cgplot, /over, zarr, a26_sumvf, color='blue', thick=thinline
	    cgplot, /over, zarr, a27_sumvf, color='green', thick=thinline

	    ;cgplot, /over, zarr, s82c_a25_raw, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a26_raw, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a27_raw, color='green', thick=thinline

	    ;cgplot, /over, zarr, s82n_a25_raw, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a26_raw, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a27_raw, color='green', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 10 - arms winding

	t10_votelim = 10
	
	a28_thresh = fltarr(nz)
	a29_thresh = fltarr(nz)
	a30_thresh = fltarr(nz)
	t10_unc_thresh = fltarr(nz)

	a28_wmean = fltarr(nz)
	a29_wmean = fltarr(nz)
	a30_wmean = fltarr(nz)
	
	a28_sumvf = fltarr(nz)
	a29_sumvf = fltarr(nz)
	a30_sumvf = fltarr(nz)
	
	    s82n_a28_thresh = fltarr(nz)
	    s82n_a29_thresh = fltarr(nz)
	    s82n_a30_thresh = fltarr(nz)
	    s82n_t10_unc_thresh = fltarr(nz)
	    s82n_a28_raw = fltarr(nz)
	    s82n_a29_raw = fltarr(nz)
	    s82n_a30_raw = fltarr(nz)
	    
	    s82c_a28_thresh = fltarr(nz)
	    s82c_a29_thresh = fltarr(nz)
	    s82c_a30_thresh = fltarr(nz)
	    s82c_t10_unc_thresh = fltarr(nz)
	    s82c_a28_raw = fltarr(nz)
	    s82c_a29_raw = fltarr(nz)
	    s82c_a30_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz

		; Main sample

		a28_thresh_ind = where(gz2main.(wfind[29]) gt 0.8 and gz2main.(wfind[31]+2) ge t10_votelim and xvar_main ge zlo and xvar_main lt zhi, a28_count)
		a29_thresh_ind = where(gz2main.(wfind[30]) gt 0.8 and gz2main.(wfind[31]+2) ge t10_votelim and xvar_main ge zlo and xvar_main lt zhi, a29_count)
		a30_thresh_ind = where(gz2main.(wfind[31]) gt 0.8 and gz2main.(wfind[31]+2) ge t10_votelim and xvar_main ge zlo and xvar_main lt zhi, a30_count)
		t10_ind     = where(                              gz2main.(wfind[31]+2) ge t10_votelim and xvar_main ge zlo and xvar_main lt zhi, t10_count)
	
		a28_thresh[i] = float(a28_count)/float(t10_count)
		a29_thresh[i] = float(a29_count)/float(t10_count)
		a30_thresh[i] = float(a30_count)/float(t10_count)
		t10_unc_thresh[i] = 1. - float(a28_count + a29_count + a30_count)/float(t10_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t10_ind].(wfind[13]+2)
		a28_wmean[i] = total(gz2main[t10_ind].(wfind[29]) * weights) / total(weights)
		a29_wmean[i] = total(gz2main[t10_ind].(wfind[30]) * weights) / total(weights)
		a30_wmean[i] = total(gz2main[t10_ind].(wfind[31]) * weights) / total(weights)

		; Sum the vote fractions
	
		a28_sumvf[i] = total(gz2main[t10_ind].(wfind[29])) / t10_count
		a29_sumvf[i] = total(gz2main[t10_ind].(wfind[30])) / t10_count
		a30_sumvf[i] = total(gz2main[t10_ind].(wfind[31])) / t10_count

		    ; Stripe 82 - normal

		    s82n_a28_thresh_ind = where(normal.(wfind[29]) gt 0.8 and normal.(wfind[31]+2) ge t10_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a28_count)
		    s82n_a29_thresh_ind = where(normal.(wfind[30]) gt 0.8 and normal.(wfind[31]+2) ge t10_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a29_count)
		    s82n_a30_thresh_ind = where(normal.(wfind[31]) gt 0.8 and normal.(wfind[31]+2) ge t10_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a30_count)
		    s82n_t10_ind        = where(                              normal.(wfind[31]+2) ge t10_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t10_count)
	
		    s82n_a28_thresh[i] = float(s82n_a28_count)/float(s82n_t10_count)
		    s82n_a29_thresh[i] = float(s82n_a29_count)/float(s82n_t10_count)
		    s82n_a30_thresh[i] = float(s82n_a30_count)/float(s82n_t10_count)
		    s82n_t10_unc_thresh[i] = 1. - float(s82n_a28_count + s82n_a29_count + s82n_a30_count)/float(s82n_t10_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82n_weights = normal[s82n_t10_ind].(wfind[13]+2)
		    s82n_a28_raw[i] = total(normal[s82n_t10_ind].(wfind[29]) * s82n_weights) / total(s82n_weights)
		    s82n_a29_raw[i] = total(normal[s82n_t10_ind].(wfind[30]) * s82n_weights) / total(s82n_weights)
		    s82n_a30_raw[i] = total(normal[s82n_t10_ind].(wfind[31]) * s82n_weights) / total(s82n_weights)

		    ; Stripe 82 - coadd2

		    s82c_a28_thresh_ind = where(coadd2.(wfind[29]) gt 0.8 and coadd2.(wfind[31]+2) ge t10_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a28_count)
		    s82c_a29_thresh_ind = where(coadd2.(wfind[30]) gt 0.8 and coadd2.(wfind[31]+2) ge t10_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a29_count)
		    s82c_a30_thresh_ind = where(coadd2.(wfind[31]) gt 0.8 and coadd2.(wfind[31]+2) ge t10_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a30_count)
		    s82c_t10_ind        = where(                              coadd2.(wfind[31]+2) ge t10_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t10_count)
	
		    s82c_a28_thresh[i] = float(s82c_a28_count)/float(s82c_t10_count)
		    s82c_a29_thresh[i] = float(s82c_a29_count)/float(s82c_t10_count)
		    s82c_a30_thresh[i] = float(s82c_a30_count)/float(s82c_t10_count)
		    s82c_t10_unc_thresh[i] = 1. - float(s82c_a28_count + s82c_a29_count + s82c_a30_count)/float(s82c_t10_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82c_weights = coadd2[s82c_t10_ind].(wfind[13]+2)
		    s82c_a28_raw[i] = total(coadd2[s82c_t10_ind].(wfind[29]) * s82c_weights) / total(s82c_weights)
		    s82c_a29_raw[i] = total(coadd2[s82c_t10_ind].(wfind[30]) * s82c_weights) / total(s82c_weights)
		    s82c_a30_raw[i] = total(coadd2[s82c_t10_ind].(wfind[31]) * s82c_weights) / total(s82c_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a28_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 10', charsize = cs
	cgplot, /over, zarr, a29_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a30_thresh, color='green', thick=thickline
	cgplot, /over, zarr, t10_unc_thresh, color='orange', thick=thickline
	
	    ;cgplot, /over, zarr, s82n_a28_thresh, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a29_thresh, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a30_thresh, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_t10_unc_thresh, color='orange', thick=thinline, linestyle=2
	    ;
	    ;cgplot, /over, zarr, s82c_a28_thresh, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a29_thresh, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a30_thresh, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_t10_unc_thresh, color='orange', thick=thinline
	
	al_legend,/top,/left,['Tight spiral','Medium spiral','Loose spiral','Unclassified'],color=['red','blue','green','orange'],psym=28, charsize = legendcs

	cgplot,        zarr, a28_wmean, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 10', charsize = cs
	cgplot, /over, zarr, a29_wmean, color='blue', thick=thickline
	cgplot, /over, zarr, a30_wmean, color='green', thick=thickline

	    cgplot, /over, zarr, a28_sumvf, color='red', thick=thinline
	    cgplot, /over, zarr, a29_sumvf, color='blue', thick=thinline
	    cgplot, /over, zarr, a30_sumvf, color='green', thick=thinline

	    ;cgplot, /over, zarr, s82c_a28_raw, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a29_raw, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a30_raw, color='green', thick=thinline

	    ;cgplot, /over, zarr, s82n_a28_raw, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a29_raw, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a30_raw, color='green', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	; Task 11 - arms number

	t11_votelim = 10
	
	a31_thresh = fltarr(nz)
	a32_thresh = fltarr(nz)
	a33_thresh = fltarr(nz)
	a34_thresh = fltarr(nz)
	a36_thresh = fltarr(nz)
	a37_thresh = fltarr(nz)
	t11_unc_thresh = fltarr(nz)

	a31_wmean = fltarr(nz)
	a32_wmean = fltarr(nz)
	a33_wmean = fltarr(nz)
	a34_wmean = fltarr(nz)
	a36_wmean = fltarr(nz)
	a37_wmean = fltarr(nz)
	
	a31_sumvf = fltarr(nz)
	a32_sumvf = fltarr(nz)
	a33_sumvf = fltarr(nz)
	a34_sumvf = fltarr(nz)
	a36_sumvf = fltarr(nz)
	a37_sumvf = fltarr(nz)
	
	    s82n_a31_thresh = fltarr(nz)
	    s82n_a32_thresh = fltarr(nz)
	    s82n_a33_thresh = fltarr(nz)
	    s82n_a34_thresh = fltarr(nz)
	    s82n_a36_thresh = fltarr(nz)
	    s82n_a37_thresh = fltarr(nz)
	    s82n_t11_unc_thresh = fltarr(nz)
	    s82n_a31_raw = fltarr(nz)
	    s82n_a32_raw = fltarr(nz)
	    s82n_a33_raw = fltarr(nz)
	    s82n_a34_raw = fltarr(nz)
	    s82n_a36_raw = fltarr(nz)
	    s82n_a37_raw = fltarr(nz)
	    
	    s82c_a31_thresh = fltarr(nz)
	    s82c_a32_thresh = fltarr(nz)
	    s82c_a33_thresh = fltarr(nz)
	    s82c_a34_thresh = fltarr(nz)
	    s82c_a36_thresh = fltarr(nz)
	    s82c_a37_thresh = fltarr(nz)
	    s82c_t11_unc_thresh = fltarr(nz)
	    s82c_a31_raw = fltarr(nz)
	    s82c_a32_raw = fltarr(nz)
	    s82c_a33_raw = fltarr(nz)
	    s82c_a34_raw = fltarr(nz)
	    s82c_a36_raw = fltarr(nz)
	    s82c_a37_raw = fltarr(nz)
	
	for i=0,nz-1 do begin
		zlo = zarr[i]
		zhi = zarr[i]+delz

		; Main sample

		a31_thresh_ind = where(gz2main.(wfind[32]) gt 0.8 and gz2main.(wfind[37]+2) ge t11_votelim and xvar_main ge zlo and xvar_main lt zhi, a31_count)
		a32_thresh_ind = where(gz2main.(wfind[33]) gt 0.8 and gz2main.(wfind[37]+2) ge t11_votelim and xvar_main ge zlo and xvar_main lt zhi, a32_count)
		a33_thresh_ind = where(gz2main.(wfind[34]) gt 0.8 and gz2main.(wfind[37]+2) ge t11_votelim and xvar_main ge zlo and xvar_main lt zhi, a33_count)
		a34_thresh_ind = where(gz2main.(wfind[35]) gt 0.8 and gz2main.(wfind[37]+2) ge t11_votelim and xvar_main ge zlo and xvar_main lt zhi, a34_count)
		a36_thresh_ind = where(gz2main.(wfind[36]) gt 0.8 and gz2main.(wfind[37]+2) ge t11_votelim and xvar_main ge zlo and xvar_main lt zhi, a36_count)
		a37_thresh_ind = where(gz2main.(wfind[37]) gt 0.8 and gz2main.(wfind[37]+2) ge t11_votelim and xvar_main ge zlo and xvar_main lt zhi, a37_count)
		t11_ind        = where(                               gz2main.(wfind[37]+2) ge t11_votelim and xvar_main ge zlo and xvar_main lt zhi, t11_count)
	
		a31_thresh[i] = float(a31_count)/float(t11_count)
		a32_thresh[i] = float(a32_count)/float(t11_count)
		a33_thresh[i] = float(a33_count)/float(t11_count)
		a34_thresh[i] = float(a34_count)/float(t11_count)
		a36_thresh[i] = float(a36_count)/float(t11_count)
		a37_thresh[i] = float(a37_count)/float(t11_count)
		t11_unc_thresh[i] = 1. - float(a31_count + a32_count + a33_count + a34_count + a36_count + a37_count)/float(t11_count)
	
		; Mean weighted by total number of votes for each galaxy on that task
	
		weights = gz2main[t11_ind].(wfind[37]+2)
		a31_wmean[i] = total(gz2main[t11_ind].(wfind[32]) * weights) / total(weights)
		a32_wmean[i] = total(gz2main[t11_ind].(wfind[33]) * weights) / total(weights)
		a33_wmean[i] = total(gz2main[t11_ind].(wfind[34]) * weights) / total(weights)
		a34_wmean[i] = total(gz2main[t11_ind].(wfind[35]) * weights) / total(weights)
		a36_wmean[i] = total(gz2main[t11_ind].(wfind[36]) * weights) / total(weights)
		a37_wmean[i] = total(gz2main[t11_ind].(wfind[37]) * weights) / total(weights)

		; Mean weighted by total number of votes for each galaxy on that task
	
		a31_sumvf[i] = total(gz2main[t11_ind].(wfind[32])) / t11_count
		a32_sumvf[i] = total(gz2main[t11_ind].(wfind[33])) / t11_count
		a33_sumvf[i] = total(gz2main[t11_ind].(wfind[34])) / t11_count
		a34_sumvf[i] = total(gz2main[t11_ind].(wfind[35])) / t11_count
		a36_sumvf[i] = total(gz2main[t11_ind].(wfind[36])) / t11_count
		a37_sumvf[i] = total(gz2main[t11_ind].(wfind[37])) / t11_count

		    ; Stripe 82 - normal

		    s82n_a31_thresh_ind = where(normal.(wfind[32]) gt 0.8 and normal.(wfind[37]+2) ge t11_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a31_count)
		    s82n_a32_thresh_ind = where(normal.(wfind[33]) gt 0.8 and normal.(wfind[37]+2) ge t11_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a32_count)
		    s82n_a33_thresh_ind = where(normal.(wfind[34]) gt 0.8 and normal.(wfind[37]+2) ge t11_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a33_count)
		    s82n_a34_thresh_ind = where(normal.(wfind[35]) gt 0.8 and normal.(wfind[37]+2) ge t11_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a34_count)
		    s82n_a36_thresh_ind = where(normal.(wfind[36]) gt 0.8 and normal.(wfind[37]+2) ge t11_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a36_count)
		    s82n_a37_thresh_ind = where(normal.(wfind[37]) gt 0.8 and normal.(wfind[37]+2) ge t11_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_a37_count)
		    s82n_t11_ind        = where(                              normal.(wfind[37]+2) ge t11_votelim and xvar_normal ge zlo and xvar_normal lt zhi, s82n_t11_count)
	
		    s82n_a31_thresh[i] = float(s82n_a31_count)/float(s82n_t11_count)
		    s82n_a32_thresh[i] = float(s82n_a32_count)/float(s82n_t11_count)
		    s82n_a33_thresh[i] = float(s82n_a33_count)/float(s82n_t11_count)
		    s82n_a34_thresh[i] = float(s82n_a34_count)/float(s82n_t11_count)
		    s82n_a36_thresh[i] = float(s82n_a36_count)/float(s82n_t11_count)
		    s82n_a37_thresh[i] = float(s82n_a37_count)/float(s82n_t11_count)
		    s82n_t11_unc_thresh[i] = 1. - float(s82n_a31_count + s82n_a32_count + s82n_a33_count + s82n_a34_count + s82n_a36_count + s82n_a37_count)/float(s82n_t11_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82n_weights = normal[s82n_t11_ind].(wfind[37]+2)
		    s82n_a31_raw[i] = total(normal[s82n_t11_ind].(wfind[32]) * s82n_weights) / total(s82n_weights)
		    s82n_a32_raw[i] = total(normal[s82n_t11_ind].(wfind[33]) * s82n_weights) / total(s82n_weights)
		    s82n_a33_raw[i] = total(normal[s82n_t11_ind].(wfind[34]) * s82n_weights) / total(s82n_weights)
		    s82n_a34_raw[i] = total(normal[s82n_t11_ind].(wfind[35]) * s82n_weights) / total(s82n_weights)
		    s82n_a36_raw[i] = total(normal[s82n_t11_ind].(wfind[36]) * s82n_weights) / total(s82n_weights)
		    s82n_a37_raw[i] = total(normal[s82n_t11_ind].(wfind[37]) * s82n_weights) / total(s82n_weights)

		    ; Stripe 82 - coadd2

		    s82c_a31_thresh_ind = where(coadd2.(wfind[32]) gt 0.8 and coadd2.(wfind[37]+2) ge t11_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a31_count)
		    s82c_a32_thresh_ind = where(coadd2.(wfind[33]) gt 0.8 and coadd2.(wfind[37]+2) ge t11_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a32_count)
		    s82c_a33_thresh_ind = where(coadd2.(wfind[34]) gt 0.8 and coadd2.(wfind[37]+2) ge t11_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a33_count)
		    s82c_a34_thresh_ind = where(coadd2.(wfind[35]) gt 0.8 and coadd2.(wfind[37]+2) ge t11_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a34_count)
		    s82c_a36_thresh_ind = where(coadd2.(wfind[36]) gt 0.8 and coadd2.(wfind[37]+2) ge t11_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a36_count)
		    s82c_a37_thresh_ind = where(coadd2.(wfind[37]) gt 0.8 and coadd2.(wfind[37]+2) ge t11_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_a37_count)
		    s82c_t11_ind        = where(                              coadd2.(wfind[37]+2) ge t11_votelim and xvar_coadd2 ge zlo and xvar_coadd2 lt zhi, s82c_t11_count)
	
		    s82c_a31_thresh[i] = float(s82c_a31_count)/float(s82c_t11_count)
		    s82c_a32_thresh[i] = float(s82c_a32_count)/float(s82c_t11_count)
		    s82c_a33_thresh[i] = float(s82c_a33_count)/float(s82c_t11_count)
		    s82c_a34_thresh[i] = float(s82c_a34_count)/float(s82c_t11_count)
		    s82c_a36_thresh[i] = float(s82c_a36_count)/float(s82c_t11_count)
		    s82c_a37_thresh[i] = float(s82c_a37_count)/float(s82c_t11_count)
		    s82c_t11_unc_thresh[i] = 1. - float(s82c_a31_count + s82c_a32_count + s82c_a33_count + s82c_a34_count + s82c_a36_count + s82c_a37_count)/float(s82c_t11_count)
	
		    ; Mean weighted by total number of votes for each galaxy on that task
	
		    s82c_weights = coadd2[s82c_t11_ind].(wfind[37]+2)
		    s82c_a31_raw[i] = total(coadd2[s82c_t11_ind].(wfind[32]) * s82c_weights) / total(s82c_weights)
		    s82c_a32_raw[i] = total(coadd2[s82c_t11_ind].(wfind[33]) * s82c_weights) / total(s82c_weights)
		    s82c_a33_raw[i] = total(coadd2[s82c_t11_ind].(wfind[34]) * s82c_weights) / total(s82c_weights)
		    s82c_a34_raw[i] = total(coadd2[s82c_t11_ind].(wfind[35]) * s82c_weights) / total(s82c_weights)
		    s82c_a36_raw[i] = total(coadd2[s82c_t11_ind].(wfind[36]) * s82c_weights) / total(s82c_weights)
		    s82c_a37_raw[i] = total(coadd2[s82c_t11_ind].(wfind[37]) * s82c_weights) / total(s82c_weights)

	endfor
	
	!p.multi=[0,2,1]
	cgplot,        zarr, a31_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 11', charsize = cs
	cgplot, /over, zarr, a32_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a33_thresh, color='green', thick=thickline
	cgplot, /over, zarr, a34_thresh, color='purple', thick=thickline
	cgplot, /over, zarr, a36_thresh, color='grey', thick=thickline
	cgplot, /over, zarr, a37_thresh, color='brown', thick=thickline
	cgplot, /over, zarr, t11_unc_thresh, color='orange', thick=thickline
	
	    ;cgplot, /over, zarr, s82c_a31_thresh, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a32_thresh, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a33_thresh, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_a34_thresh, color='purple', thick=thinline
	    ;cgplot, /over, zarr, s82c_a36_thresh, color='grey', thick=thinline
	    ;cgplot, /over, zarr, s82c_a37_thresh, color='brown', thick=thinline
	    ;cgplot, /over, zarr, s82c_t11_unc_thresh, color='orange', thick=thinline
	    ;
	    ;cgplot, /over, zarr, s82n_a31_thresh, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a32_thresh, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a33_thresh, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a34_thresh, color='purple', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a36_thresh, color='grey', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a37_thresh, color='brown', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_t11_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['1 arm','2 arms','3 arms','4 arms','4+ arms',"Can't tell",'Unclassified'],color=['red','blue','green','purple','grey','brown','orange'],psym=28, charsize = legendcs

	cgplot,        zarr, a31_sumvf, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, title='GZ2 Task 11', charsize = cs
	cgplot, /over, zarr, a32_sumvf, color='blue', thick=thickline
	cgplot, /over, zarr, a33_sumvf, color='green', thick=thickline
	cgplot, /over, zarr, a34_sumvf, color='purple', thick=thickline
	cgplot, /over, zarr, a36_sumvf, color='grey', thick=thickline
	cgplot, /over, zarr, a37_sumvf, color='brown', thick=thickline

	    cgplot, /over, zarr, a31_sumvf, color='red', thick=thinline
	    cgplot, /over, zarr, a32_sumvf, color='blue', thick=thinline
	    cgplot, /over, zarr, a33_sumvf, color='green', thick=thinline
	    cgplot, /over, zarr, a34_sumvf, color='purple', thick=thinline
	    cgplot, /over, zarr, a36_sumvf, color='grey', thick=thinline
	    cgplot, /over, zarr, a37_sumvf, color='brown', thick=thinline

	    ;cgplot, /over, zarr, s82c_a31_raw, color='red', thick=thinline
	    ;cgplot, /over, zarr, s82c_a32_raw, color='blue', thick=thinline
	    ;cgplot, /over, zarr, s82c_a33_raw, color='green', thick=thinline
	    ;cgplot, /over, zarr, s82c_a34_raw, color='purple', thick=thinline
	    ;cgplot, /over, zarr, s82c_a36_raw, color='grey', thick=thinline
	    ;cgplot, /over, zarr, s82c_a37_raw, color='brown', thick=thinline

	    ;cgplot, /over, zarr, s82n_a31_raw, color='red', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a32_raw, color='blue', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a33_raw, color='green', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a34_raw, color='purple', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a36_raw, color='grey', thick=thinline, linestyle=2
	    ;cgplot, /over, zarr, s82n_a37_raw, color='brown', thick=thinline, linestyle=2

	;al_legend,/top,/left,['Main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize = legendcs
	al_legend,/top,/left,['W. mean','Sum'],color='black',linestyle=0, thick=[thickline,thinline], linsize=0.4, charsize=legendcs

	endif	; plotone

if keyword_set(ps) then ps_end

; End time

timeend = systime(1)

print,''
print,'Time elapsed: '+string(timeend-timestart,format='(f5.1)')+' sec.'
print,''

;stop

end
