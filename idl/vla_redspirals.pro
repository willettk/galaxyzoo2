
;+
; NAME:
;       
;	
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


!p.multi=[0,2,2]

zarr = fillarr(1d-4,0.01,0.10)
cgplot, zarr, snu(1d21,zarr), xtitle='Redshift', ytitle='Expected flux density [mJy]', thick=2, yr=[0,1], title='Sensitivity to radio sources'
cgplot, zarr, snu(10^18.0,zarr), /over, linestyle=2
cgplot, zarr, snu(10^18.5,zarr), /over, linestyle=2
cgplot, zarr, snu(10^19.0,zarr), /over, linestyle=2
cgplot, zarr, snu(10^19.5,zarr), /over, linestyle=2
cgplot, zarr, snu(10^20.0,zarr), /over, linestyle=2
cgplot, zarr, snu(10^20.5,zarr), /over, linestyle=2
cgplot, zarr, snu(10^21.5,zarr), /over, linestyle=2
cgplot, zarr, snu(10^22.0,zarr), /over, linestyle=2
cgplots, !x.crange, [0.1,0.1], color='blue'
cgplots, !x.crange, [0.020,0.020], color='blue'

file='~/Astronomy/proposals/redspirals_chandra/masters2010/mnr_16503_sm_TableA1.txt'

readcol, file, $
	objid, ra, dec, redshift, $
	absmag_r, absmag_r_err, color_gr, color_gr_err, redness, ellipticity, fracdeV, $
	format='a,f,f,f,f,f,f,f,f,f', $
	comment='%', $
	/silent, $
	skip=4
	
g=[11.57, 12.56, 12.46, 13.45, 13.23, 11.77, 12.25, 12.58, 17.17, 12.47]

r=[10.68, 11.40, 11.43, 11.90, 12.44, 10.92, 11.38, 11.73, 15.16, 11.70]

dm=[31.02, 31.56, 31.86, 31.15, 32.86, 31.03, 31.04, 32.30, 31.75, 32.20]

magarr = fillarr(1d-2,-23,-16)
cgplot, r - dm, g-r, psym = 16, xtitle='M!Ir!N',ytitle='(g-r)', xr=[-16,-23], /xstyle, title='CMD of late-type LINERs from Nagar'
cgplot, /over, magarr, 0.63 - 0.02*(magarr + 20)

	; What is the total time it would take to observe all the LINERs in the sample at 20 uJy sensitivity?

	tobs = 4.	; minutes
	; Even at z=0.08, a 10^21 galaxy is still a 56 mJy source
	; print,n_elements(redshift) * tobs / 60. 
	; 293 galaxies, even at 4m apiece, would take 19 hours. Could ask for half of the closer sample.

a = mrdfits('~/Astronomy/proposals/redspirals_vla/nagar2005_table1.fits',1,/silent)

liner = where(strmid(a.atype,0,1) eq 'L')
liner_e = where(strmid(a.atype,0,1) eq 'L' and a.ttype le 0 and a.l_lp15ghz ne '<' and a.lp15ghz ne 0)
liner_s = where(strmid(a.atype,0,1) eq 'L' and a.ttype gt 0 and a.l_lp15ghz ne '<' and a.lp15ghz ne 0)
linerlim_e = where(strmid(a.atype,0,1) eq 'L' and a.ttype le 0 and a.l_lp15ghz eq '<' and a.lp15ghz ne 0)
linerlim_s = where(strmid(a.atype,0,1) eq 'L' and a.ttype gt 0 and a.l_lp15ghz eq '<' and a.lp15ghz ne 0)
cghistoplot,a[liner_e].lp15ghz, datacolor='blue', bin=1, thick=3, yr=[0,15], xr=[17,25],/outline, title='Nagar LINER radio luminosities'
cghistoplot,a[liner_s].lp15ghz, /oplot, datacolor='red', bin=1, thick=3                ,/outline
cghistoplot,a[linerlim_e].lp15ghz, /oplot, datacolor='blue', bin=1, thick = 1          ,/outline
cghistoplot,a[linerlim_s].lp15ghz, /oplot, datacolor='red', bin=1, thick = 1           ,/outline
cgplots, replicate(mean(a[liner_e].lp15ghz),2), !y.crange, color='red',linestyle=2, thick=2
cgplots, replicate(mean(a[liner_s].lp15ghz),2), !y.crange, color='blue',linestyle=2, thick=2

end
