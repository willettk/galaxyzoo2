
;+
; NAME:
;       
;	PETROR_TESTS.pro
;
; PURPOSE:
;
;	Find out what the conversion between angular and physical size in the GZ2 tables is. 
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


device, retain=2

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
figsdir = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'

red, /default
	
;if n_elements(lowz) gt 0 then begin
;	coadd_all = mrdfits(gz2dir+'fits/coadd1_coadd2_sample.fits',1,/silent)
;	coadd = coadd_all[where(coadd_all.redshift gt 0. and coadd_all.redshift le lowz)]
;endif else coadd = mrdfits(gz2dir+'fits/coadd1_coadd2.fits',1,/silent)

coadd = mrdfits(gz2dir+'fits/coadd1_coadd2_sample.fits',1,/silent)

tagnames=tag_names(coadd)
nt = n_elements(tagnames)

wfind_temp=intarr(nt) 
wcind_temp=intarr(nt) 

for i=0,nt-1 do begin
	if (strlen(tagnames[i]) ge 20) and (strmid(tagnames[i],strlen(tagnames[i])-19,19) eq 'WEIGHTED_FRACTION_2') then wfind_temp[i]=1 
	if (strlen(tagnames[i]) ge 15) and (strmid(tagnames[i],strlen(tagnames[i])-14,14) eq 'TOTAL_WEIGHT_2') then wcind_temp[i]=1 
endfor

wfind = where(wfind_temp, nw)	; Weighted fraction
wcind = where(wcind_temp, nc)	; Total weight (count)

; Select only galaxies with redshifts

goodz = where(finite(coadd.redshift))

coaddz = coadd[goodz]

; Angular size from the table

tablesize = coaddz.petror50_r_kpc

; Angular size from my conversion calculation

rad2arcsec = 206265.0d
mysize = coaddz.petror50_r * dangular(coaddz.redshift,/kpc)/rad2arcsec

red, omega0=0.3, omegalambda=0.7, h100=0.70
mysize_def = coaddz.petror50_r * dangular(coaddz.redshift,/kpc)/rad2arcsec

cgplot, tablesize, mysize_def - tablesize, psym=3, color='cyan'

stop

end
