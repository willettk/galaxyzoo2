;+
; NAME:
;       
;	GZ2_BIAS_DEMO_TASK01
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
; Changed yrange to [0,1.0] - 2013-03-15
;-

pro gz2_bias_demo_task01, ps=ps

timestart = systime(1)

savdir = '~/Astronomy/Research/GalaxyZoo/sav/'
fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
figsdir = '~/Astronomy/Research/GalaxyZoo/datapaper/figures/'
gz2tablesamplefile = fitsdir+'gz2_table_sample_match_coadd2.fits'

	;gz2 = mrdfits(gz2tablesamplefile, 1, /silent)
	;save,gz2,filename=fitsdir+'gz2.sav'

restore,savdir+'gz2.sav'

tagnames=tag_names(gz2)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]	

gz2main = gz2[where((strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra'),maincount)]
coadd2 = gz2[where(strtrim(gz2.sample,2) eq 'stripe82_coadd_2' and (gz2.petromag_r - gz2.extinction_r) le 17.0 and gz2.petror90_r ge 3.,s82c_count)]
normal = gz2[where(strtrim(gz2.sample,2) eq 'stripe82'         and (gz2.petromag_r - gz2.extinction_r) le 17.0 and gz2.petror90_r ge 3.,s82n_count)]

; Two plots: thresholded likelihoods (try f > 0.8), and raw vote summation for all answers as a function of redshift. Repeat as function of surface brightness, physical size, and apparent magnitude

yrange = [0,1.0]

delz=0.02
xlabel = 'redshift'
zarr = fillarr(delz,0,0.20)
nz = n_elements(zarr)

xvar_main = gz2main.redshift
xvar_coadd2 = coadd2.redshift
xvar_normal = normal.redshift

xrange=[min(zarr),max(zarr)*1.1]

if keyword_set(ps) then begin

	ps_start, filename=figsdir+'gz2_bias_demo_task01.ps', /color, /quiet, xsize=7, ysize=3.5
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

	t01_votelim = 10
	
	a01_thresh = fltarr(nz)
	a02_thresh = fltarr(nz)
	a03_thresh = fltarr(nz)
	t01_unc_thresh = fltarr(nz)
	
	a01_raw = fltarr(nz)
	a02_raw = fltarr(nz)
	a03_raw = fltarr(nz)
	
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
		a01_raw[i] = total(gz2main[t01_ind].(wfind[1]) * weights) / total(weights)
		a02_raw[i] = total(gz2main[t01_ind].(wfind[2]) * weights) / total(weights)
		a03_raw[i] = total(gz2main[t01_ind].(wfind[3]) * weights) / total(weights)

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
	cgplot,        zarr, a01_thresh, color='red', thick=thickline, xtitle=xlabel,ytitle='fraction', xr=xrange, /xstyle, yr=yrange, /ystyle, charsize = cs, $
        position=[0.10,0.15,0.50,0.95]
	cgplot, /over, zarr, a02_thresh, color='blue', thick=thickline
	cgplot, /over, zarr, a03_thresh, color='green', thick=thickline
	cgplot, /over, zarr, t01_unc_thresh, color='orange', thick=thickline
	
		cgplot, /over, zarr, s82c_a01_thresh, color='red', thick=thinline
		cgplot, /over, zarr, s82c_a02_thresh, color='blue', thick=thinline
		cgplot, /over, zarr, s82c_a03_thresh, color='green', thick=thinline
		cgplot, /over, zarr, s82c_t01_unc_thresh, color='orange', thick=thinline
	
		cgplot, /over, zarr, s82n_a01_thresh, color='red', thick=thinline, linestyle=2
		cgplot, /over, zarr, s82n_a02_thresh, color='blue', thick=thinline, linestyle=2
		cgplot, /over, zarr, s82n_a03_thresh, color='green', thick=thinline, linestyle=2
		cgplot, /over, zarr, s82n_t01_unc_thresh, color='orange', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['smooth','features','artifact','unclassified'],color=['red','blue','green','orange'],psym=28, charsize = legendcs
	
	cgplot,        zarr, a01_raw, color='red', thick=thickline, xtitle=xlabel,                   xr=xrange, /xstyle, yr=yrange, /ystyle, charsize = cs, $
        position=[0.60,0.15,1.0,0.95]
	cgplot, /over, zarr, a02_raw, color='blue', thick=thickline
	cgplot, /over, zarr, a03_raw, color='green', thick=thickline
	
		cgplot, /over, zarr, s82c_a01_raw, color='red', thick=thinline
		cgplot, /over, zarr, s82c_a02_raw, color='blue', thick=thinline
		cgplot, /over, zarr, s82c_a03_raw, color='green', thick=thinline

		cgplot, /over, zarr, s82n_a01_raw, color='red', thick=thinline, linestyle=2
		cgplot, /over, zarr, s82n_a02_raw, color='blue', thick=thinline, linestyle=2
		cgplot, /over, zarr, s82n_a03_raw, color='green', thick=thinline, linestyle=2
	
	al_legend,/top,/left,['main','S82 normal','S82 coadd2'],color='black',linestyle=[0,2,0], thick=[thickline,thinline,thinline], linsize=0.4, charsize=legendcs

if keyword_set(ps) then ps_end

; End time

timeend = systime(1)

print,''
print,'Time elapsed: '+string(timeend-timestart,format='(f5.1)')+' sec.'
print,''

end

