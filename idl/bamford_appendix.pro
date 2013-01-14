
;+
; NAME:
;       
;	BAMFORD_APPENDIX
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

pro bamford_appendix, ps=ps 

timestart = systime(1)

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
figsdir = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'

gz2main = mrdfits(fitsdir+'gz2main_table_sample.fits',1,/silent)

tagnames=tag_names(gz2main)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]	

absr = gz2main.petromag_mr
z = gz2main.redshift
	
;; Figure A2 (Bamford et al. 2009)
;
;magbin1 = where(absr gt -23.00 and absr lt -22.75,bc1)
;magbin2 = where(absr gt -21.25 and absr lt -21.00,bc2)
;magbin3 = where(absr gt -19.75 and absr lt -19.50,bc3)
;
;	t01_votelim = 10
;	t01_galaxy_lim = 10
;	
;	a01_mag1_raw = fltarr(nz)
;	a02_mag1_raw = fltarr(nz)
;	a03_mag1_raw = fltarr(nz)
;	
;	a01_mag2_raw = fltarr(nz)
;	a02_mag2_raw = fltarr(nz)
;	a03_mag2_raw = fltarr(nz)
;	
;	a01_mag3_raw = fltarr(nz)
;	a02_mag3_raw = fltarr(nz)
;	a03_mag3_raw = fltarr(nz)
;	
;	for i=0,nz-1 do begin
;		zlo = zarr[i]
;		zhi = zarr[i]+delz
;		t01_magbin1_ind = where(absr gt -23.00 and absr lt -22.75 and gz2main.(wfind[3]+2) ge t01_votelim and gz2main.redshift ge zlo and gz2main.redshift lt zhi, t01_magbin1_count)
;		t01_magbin2_ind = where(absr gt -21.25 and absr lt -21.00 and gz2main.(wfind[3]+2) ge t01_votelim and gz2main.redshift ge zlo and gz2main.redshift lt zhi, t01_magbin2_count)
;		t01_magbin3_ind = where(absr gt -19.75 and absr lt -19.50 and gz2main.(wfind[3]+2) ge t01_votelim and gz2main.redshift ge zlo and gz2main.redshift lt zhi, t01_magbin3_count)
;	
;		; Mean weighted by total number of votes for each galaxy on that task
;	
;		weights_mag1 = gz2main[t01_magbin1_ind].(wfind[3]+2)
;		a01_mag1_raw[i] = total(gz2main[t01_magbin1_ind].(wfind[1]) * weights_mag1) / total(weights_mag1)
;		a02_mag1_raw[i] = total(gz2main[t01_magbin1_ind].(wfind[2]) * weights_mag1) / total(weights_mag1)
;		a03_mag1_raw[i] = total(gz2main[t01_magbin1_ind].(wfind[3]) * weights_mag1) / total(weights_mag1)
;
;		weights_mag2 = gz2main[t01_magbin2_ind].(wfind[3]+2)
;		a01_mag2_raw[i] = total(gz2main[t01_magbin2_ind].(wfind[1]) * weights_mag2) / total(weights_mag2)
;		a02_mag2_raw[i] = total(gz2main[t01_magbin2_ind].(wfind[2]) * weights_mag2) / total(weights_mag2)
;		a03_mag2_raw[i] = total(gz2main[t01_magbin2_ind].(wfind[3]) * weights_mag2) / total(weights_mag2)
;
;		weights_mag3 = gz2main[t01_magbin3_ind].(wfind[3]+2)
;		a01_mag3_raw[i] = total(gz2main[t01_magbin3_ind].(wfind[1]) * weights_mag3) / total(weights_mag3)
;		a02_mag3_raw[i] = total(gz2main[t01_magbin3_ind].(wfind[2]) * weights_mag3) / total(weights_mag3)
;		a03_mag3_raw[i] = total(gz2main[t01_magbin3_ind].(wfind[3]) * weights_mag3) / total(weights_mag3)
;
;		if t01_magbin1_count lt t01_galaxy_lim then begin
;			a01_mag1_raw[i] = !values.f_nan
;			a02_mag1_raw[i] = !values.f_nan
;			a03_mag1_raw[i] = !values.f_nan
;		endif
;
;		if t01_magbin2_count lt t01_galaxy_lim then begin
;			a01_mag2_raw[i] = !values.f_nan
;			a02_mag2_raw[i] = !values.f_nan
;			a03_mag2_raw[i] = !values.f_nan
;		endif
;
;		if t01_magbin3_count lt t01_galaxy_lim then begin
;			a01_mag3_raw[i] = !values.f_nan
;			a02_mag3_raw[i] = !values.f_nan
;			a03_mag3_raw[i] = !values.f_nan
;		endif
;
;
;	endfor
;
;!p.multi=[0,3,1]
;
;cgplot,        zarr, a01_mag1_raw, color='red', thick=thickline, xtitle='Redshift',ytitle='fraction',xr=[-0.02,0.30], /xstyle, yr=[0,1], /ystyle, title='GZ2 Task 01'
;cgplot, /over, zarr, a02_mag1_raw, color='blue', thick=thickline
;cgplot, /over, zarr, a03_mag1_raw, color='green', thick=thickline
;
;cgtext, /data, 0.15, 0.90, 'N!Igal!N ='+string(n_elements(where(absr gt -23.00 and absr lt -22.75 and gz2main.(wfind[3]+2) ge t01_votelim)),format='(i7)')
;
;al_legend,/top,/left,['smooth','features/disk','star/artifact'],color=['red','blue','green'],linestyle=0, thick=thickline, linsize=0.4
;
;cgplot,        zarr, a01_mag2_raw, color='red', thick=thickline, xtitle='Redshift',ytitle='fraction',xr=[-0.02,0.30], /xstyle, yr=[0,1], /ystyle, title='GZ2 Task 01'
;cgplot, /over, zarr, a02_mag2_raw, color='blue', thick=thickline
;cgplot, /over, zarr, a03_mag2_raw, color='green', thick=thickline
;
;cgtext, /data, 0.15, 0.90, 'N!Igal!N ='+string(n_elements(where(absr gt -21.25 and absr lt -21.00 and gz2main.(wfind[3]+2) ge t01_votelim)),format='(i7)')
;
;al_legend,/top,/left,['smooth','features/disk','star/artifact'],color=['red','blue','green'],linestyle=0, thick=thickline, linsize=0.4
;
;cgplot,        zarr, a01_mag3_raw, color='red', thick=thickline, xtitle='Redshift',ytitle='fraction',xr=[-0.02,0.30], /xstyle, yr=[0,1], /ystyle, title='GZ2 Task 01'
;cgplot, /over, zarr, a02_mag3_raw, color='blue', thick=thickline
;cgplot, /over, zarr, a03_mag3_raw, color='green', thick=thickline
;
;cgtext, /data, 0.15, 0.90, 'N!Igal!N ='+string(n_elements(where(absr gt -19.75 and absr lt -19.50 and gz2main.(wfind[3]+2) ge t01_votelim)),format='(i7)')
;
;al_legend,/top,/left,['smooth','features/disk','star/artifact'],color=['red','blue','green'],linestyle=0, thick=thickline, linsize=0.4

; Figure A3

z1 = [0.03, 0.04]
z2 = [0.08, 0.09]
z3 = [0.14, 0.15]
zbins =[[z1],[z2],[z3]]

zbins = [transpose(fillarr(0.01,0,0.24)),transpose(fillarr(0.01,0.01,0.25))]
nplots = (size(zbins))[2]

mag_r50_binlim = 30
mag_r50_votelim = 10

ratio_color_min = -2.2
ratio_color_max = 2.8

delmag = 0.2d
delr50 = 0.5d
minmag = -24.5d
maxmag = -17d
minr50 = 0d
maxr50 = 11d

magarr = fillarr(delmag, minmag, maxmag)
r50arr = fillarr(delr50, minr50, maxr50)
nm = n_elements(magarr)
nr = n_elements(r50arr)
sp_el_ratio = fltarr(nm,nr,nplots)
n_el = fltarr(nm,nr,nplots)
n_sp = fltarr(nm,nr,nplots)

!p.multi=[0,1,1]

ncolors=180
cgloadct, 5, clip=[0,ncolors-1]

for k=0,nplots-1 do begin

	zmean=mean(zbins[*,k])
	dl = lumdist(zmean, /wmap7, /silent)
	dmod = 5*alog10(dl*1e6) - 5.
	ang_scale = 1./zang(1.,zmean,/wmap7,/silent)
	SBapp = 23.0
	SBlim_size = indgen(200) / 10.0
	SBlim_mag = (SBapp - dmod - 2.5*alog10(6.283185*(SBlim_size/ang_scale)^2))
	appmaglim = 17.7 - dmod
	
	scalearr=fillarr(0.01,0,10)
	ang = zang(scalearr,zmean,/wmap7,/silent)
	arcsec_scale = scalearr[closeto(ang,1.0)]
	
	;ps_start, filename='/Users/willettk/Astronomy/Research/GalaxyZoo/plots/gz2_bias_ratio_z'+string(zbins[0,k],format='(f4.2)')+'.ps', /color, /encap

	;cgplot, indgen(10), back='black', axescolor='white', /nodata, $
	;	charsize=cs, $
	;	position=[0.08,0.10,0.8,0.95], $
;	;	position=[0.08,0.10+(0.85/3.-0.025)*k,0.8,0.10+(0.85/3.-0.025)*(k+1)], $
	;	xtitle='M!IR!N', $
	;	ytitle='R!I50!N [kpc]', $
	;	xr=[minmag, maxmag], /xstyle, $
	;	yr=[minr50, maxr50], /ystyle
	
	;nm = 2
	;nr = 2
	
	;for i=0,nm-1 do begin
	;	for j=0,nr-1 do begin
	;
	;		ind = where(absr gt magarr[i] and absr le magarr[i]+delmag and $
	;			gz2main.petror50_r_kpc gt r50arr[j] and gz2main.petror50_r le r50arr[j]+delr50 and $
	;			gz2main.redshift gt zbins[0,k] and gz2main.redshift le zbins[1,k] and $
	;			gz2main.(wfind[3]+2) ge mag_r50_votelim,indcount)
	;		weights = gz2main[ind].(wfind[3]+2)
	;		el_raw = total(gz2main[ind].(wfind[1]) * weights) / total(weights)
	;		sp_raw = total(gz2main[ind].(wfind[2]) * weights) / total(weights)
	;
	;	if indcount ge mag_r50_binlim then $
	;		n_el[i,j,k] = el_raw
	;		n_sp[i,j,k] = sp_raw
	;		sp_el_ratio[i,j,k] = alog10(el_raw/sp_raw)
	;		;cgcolorfill, [0,0,delmag,delmag,0]+magarr[i], [0,delr50,delr50,0,0]+r50arr[j], $
	;		;	color=fix((alog10(el_raw/sp_raw) - ratio_color_min) / (ratio_color_max - ratio_color_min) * ncolors)
	;	
	;	endfor
	;endfor
	
	;print, string(zbins[0,k],format='(f5.2)')+'<z<'+string(zbins[1,k],format='(f5.2)')

	;cgplots, !x.crange, [arcsec_scale, arcsec_scale], color='white', linestyle=1, thick=3
	;cgplots, [appmaglim, appmaglim], !y.crange, color='white', linestyle=1, thick=3
	;cgplot, /over, SBlim_mag, SBlim_size, color='white', linestyle=1, thick=3
	;
	;cgtext, -19, 9.5, string(zbins[0,k],format='(f5.2)')+'<z<'+string(zbins[1,k],format='(f5.2)'), color='white', charsize=cs
	;
	;cgcolorbar, /right, range=[ratio_color_min, ratio_color_max], /vertical, title='log(n!Iel!N/n!Isp!N)', charsize=3.0

	;ps_end, /png

endfor

;save, n_el, n_sp, sp_el_ratio, filename=gz2dir+'sp_el_ratio'

; Figure A5

;sp_el_ratio_lowz = fltarr(nm,nr)
;for i=0,nm-1 do begin
;	for j=0,nr-1 do begin
;		zrow= sp_el_ratio[i,j,*]
;		sp_el_ratio_lowz[i,j] = sp_el_ratio[i,j,(where(zrow ne 0))[0]]
;	endfor
;endfor

;dummy = fltarr(nm,nr)
;for i=0,nm-1 do begin & for j=0,nr-1 do begin & zrow= sp_el_ratio[i,j,*] & dummy[i,j] = sp_el_ratio[i,j,(where(zrow ne 0 and zrow ne 9999.))[0]] & endfor & endfor &

restore,'~/Astronomy/Research/GalaxyZoo/sp_el_ratio_lowz.sav'

	cgplot, indgen(10), back='black', axescolor='white', /nodata, $
		charsize=cs, $
		position=[0.08,0.10,0.8,0.95], $
		xtitle='M!IR!N', $
		ytitle='R!I50!N [kpc]', $
		xr=[minmag, maxmag], /xstyle, $
		yr=[minr50, maxr50], /ystyle
	
	for i=0,nm-1 do begin
		for j=0,nr-1 do begin
	
			if sp_el_ratio_lowz[i,j] ne 0 then cgcolorfill, [0,0,delmag,delmag,0]+magarr[i], [0,delr50,delr50,0,0]+r50arr[j], $
				color=fix((sp_el_ratio_lowz[i,j] - ratio_color_min) / (ratio_color_max - ratio_color_min) * ncolors)
		
		endfor
	endfor
	
	zmean=0.085
	dl = lumdist(zmean, /wmap7, /silent)
	dmod = 5*alog10(dl*1e6) - 5.
	ang_scale = 1./zang(1.,zmean,/wmap7,/silent)
	SBapp = 23.0
	SBlim_size = indgen(200) / 10.0
	SBlim_mag = (SBapp - dmod - 2.5*alog10(6.283185d*(SBlim_size/ang_scale)^2))
	appmaglim = 17.7 - dmod
	
	scalearr=fillarr(0.01,0,10)
	ang = zang(scalearr,zmean,/wmap7,/silent)
	arcsec_scale = scalearr[closeto(ang,1.0)]
	
	cgplots, !x.crange, [arcsec_scale, arcsec_scale], color='white', linestyle=1, thick=3
	cgplots, [appmaglim, appmaglim], !y.crange, color='white', linestyle=1, thick=3
	cgplot, /over, SBlim_mag, SBlim_size, color='white', linestyle=1, thick=3
	
	cgtext, -19, 9.5, 'z = '+string(zmean,format='(f5.3)'), color='white', charsize=cs
	
	cgcolorbar, /right, range=[ratio_color_min, ratio_color_max], /vertical, title='log(n!Iel!N/n!Isp!N)', charsize=3.0

; Figure A6 - fit the data with a smooth function

x = magarr # replicate(1.0, nr)
y = replicate(1.0, nm) # r50arr
z = sp_el_ratio_lowz

; Steven's initial parameters, in the very weird way he organized them within the array

pinit = double([-0.3, 1.7, -22.8, 0.3, 0.3, -2.6, 1.5, 1.0, 1.0])
p = pinit[[5,6,1,0,7,8,2,3,4]]

result = gz2_task01_baseline_function(x,y,p)
nc = 50
cgloadct,5,ncolors=nc
cgcontour, result, magarr, r50arr, /fill, nlevels=nc, c_colors=indgen(nc), xtitle='M!IR!N', ytitle='R!I50!N [kpc]', charsize=2

; Assume Poissionian weights on each bin such that sig_ratio = (n_sp/n_el)^2 * (1/n_sp + 1/n_el)
task01_weights = sp_el_ratio_lowz * 0d + 0.05
nomask = where(sp_el_ratio_lowz ne 0.)
task01_weights[nomask] = 1.
szerr=sp_el_ratio_lowz*0

; Mask the bins without data
; Odd - mpfit2dfun isn't doing anything with the data. May be a problem with using weights as a mask?

	; try fixing p0 and p1 to semi-canonical values in parinfo

gz2_task01_baseline_fit = mpfit2dfun('gz2_task01_baseline_function', $
	magarr, r50arr, sp_el_ratio_lowz, szerr, p, weights=task01_weights,/quiet)

!p.multi=[0,1,2]

cgimage, $
;	result, $
	sp_el_ratio_lowz, $
	position=[0.15,0.1,0.9,0.9], $
	/axis, /scale, $
	;magarr, r50arr, $
	;/fill, nlevels=nc, c_colors=indgen(nc), $
	xtitle='M!IR!N', ytitle='R!I50!N [kpc]', charsize=2, $
	title='Starting parameters'

cgcontour, gz2_task01_baseline_function(x,y,gz2_task01_baseline_fit), magarr, r50arr, $
	position=[0.15,0.1,0.9,0.45], $
	/fill, nlevels=nc, c_colors=indgen(nc), $
	xtitle='M!IR!N', ytitle='R!I50!N [kpc]', charsize=2, title='Fit to data'


; End time

timeend = systime(1)

print,''
print,'Time elapsed: '+string(timeend-timestart,format='(f7.1)')
print,''

stop

end
