
;+
; NAME:
;       
;	STRIPE82_MRZ
;
; PURPOSE:
;
;	Plot the GZ2 Stripe 82 data in the M_R - R_50 - redshift plane.
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

pro stripe82_mrz, stop=stop, ps=ps, lowz=lowz, label=label, mingals = mingals, xbinsize = xbinsize, ybinsize = ybinsize

tstart = systime(1)

set_plot,'x'
device, retain=2

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
figsdir = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'

red, h100 = 0.704, omegalambda = 0.728
	
;if n_elements(lowz) gt 0 then begin
;	coadd_all = mrdfits(gz2dir+'fits/coadd1_coadd2_sample.fits',1,/silent)
;	coadd = coadd_all[where(coadd_all.redshift gt 0. and coadd_all.redshift le lowz)]
;endif else coadd = mrdfits(gz2dir+'fits/coadd1_coadd2.fits',1,/silent)

coadd = mrdfits(gz2dir+'fits/coadd1_coadd2_sample.fits',1,/silent)

tagnames=tag_names(coadd)
nt = n_elements(tagnames)

wfind1_temp=intarr(nt) 
wcind1_temp=intarr(nt) 
wfind2_temp=intarr(nt) 
wcind2_temp=intarr(nt) 

for i=0,nt-1 do begin
	if (strlen(tagnames[i]) ge 20) and (strmid(tagnames[i],strlen(tagnames[i])-19,19) eq 'WEIGHTED_FRACTION_1') then wfind1_temp[i]=1 
	if (strlen(tagnames[i]) ge 15) and (strmid(tagnames[i],strlen(tagnames[i])-14,14) eq 'TOTAL_WEIGHT_1') then wcind1_temp[i]=1 
	if (strlen(tagnames[i]) ge 20) and (strmid(tagnames[i],strlen(tagnames[i])-19,19) eq 'WEIGHTED_FRACTION_2') then wfind2_temp[i]=1 
	if (strlen(tagnames[i]) ge 15) and (strmid(tagnames[i],strlen(tagnames[i])-14,14) eq 'TOTAL_WEIGHT_2') then wcind2_temp[i]=1 
endfor

wfind1 = where(wfind1_temp, nw1)	; Weighted fraction
wcind1 = where(wcind1_temp, nc1)	; Total weight (count)
wfind2 = where(wfind2_temp, nw2)	; Weighted fraction
wcind2 = where(wcind2_temp, nc2)	; Total weight (count)

wfind = wfind2
wcind = wcind2

; Bin by M_R and R_50

zbin = where(coadd.redshift ge 0.03 and coadd.redshift lt 0.04)

z_binsize = 0.01
zarr = fillarr(z_binsize,0.00,0.15)
nz = n_elements(zarr)

xmin = -24
xmax = -17
ymin = 0
ymax = 10

if n_elements(xbinsize) eq 0 then xbinsize = 0.25
if n_elements(ybinsize) eq 0 then ybinsize = 0.50

xarr = fillarr(xbinsize,xmin,xmax)
yarr = fillarr(ybinsize,ymin,ymax)

nx = n_elements(xarr)
ny = n_elements(yarr)

hmin = 0
hmax = 68

cgloadct, 3


h_allz = fltarr(nx,ny)
if n_elements(mingals) eq 0 then mingals = 30



; Make an overlaid plot showing the redshift for the lowest z in each bin w/at least N galaxies. 




	
for i=0,nz - 1 do begin

	zbin = where(coadd.redshift ge zarr[i] and coadd.redshift lt zarr[i]+z_binsize, zcount)
	
	if zcount gt 0 then begin

		if keyword_set(ps) then begin
			ps_start, filename=figsdir+'stripe82_mrz/stripe82_mrz_'+string(i,format='(i02)')+'.ps', /color,/quiet
			textcolor = 'black'
		endif else begin
			textcolor = 'white'
		endelse

		h = hist_2d(coadd[zbin].petromag_mr, coadd[zbin].petror50_r_kpc, $
			bin1 = xbinsize, $
			min1 = xmin, $
			max1 = xmax, $
			bin2 = ybinsize, $
			min2 = ymin, $
			max2 = ymax)

		newind = where(h gt mingals and h_allz eq 0)
		h_allz[newind] = zarr[i]

		cgimage,h, position=[0.10,0.15,0.80,0.95], /axis, $
			minvalue=hmin, maxvalue=hmax, $
			axkeywords={xrange:[xmin,xmax], $
				yrange:[ymin,ymax], $
				xtitle:'M!IR!N [mag]', $
				ytitle:'R50 [kpc]'}
		
		cgcolorbar, /vert, position=[0.88,0.15,0.95,0.85], range=[hmin,hmax]

		cgtext, 0.83, 0.90, string(zarr[i],format='(f4.2)')+'<z<'+string(zarr[i]+z_binsize,format='(f4.2)'), charsize=1.5, color=textcolor, /normal, align=0

		if keyword_set(ps) then ps_end else wait, 0.3

	endif
	
endfor

cgloadct, 5, ncolors=nz-1
cgimage, h_allz, /axis, /scale, $
		position=[0.10,0.15,0.80,0.95], $
			axkeywords={xrange:[xmin,xmax], $
				yrange:[ymin,ymax], $
				xtitle:'M!IR!N [mag]', $
				ytitle:'R50 [kpc]'}

cgcolorbar, ncolors=nz-1,  /vert, position=[0.88,0.15,0.95,0.85], range=[min(zarr),max(zarr)]

; Note: my sizes are slightly different from the column PETROR_50_KPC in the GZ2 table since I use the WMAP7 cosmology,
; instead of the default 0.3, 0.7, 70. 

	dmod = dmodulus(median(zarr))
	ang_scale = dangular(median(zarr),/kpc) / 206265.0d
	SBapp = 23.0
	SBlim_size = (indgen(200)+1) / 10.0
	SBlim_mag = SBapp - dmod - 2.5*alog10(6.283185d * (SBlim_size/ang_scale)^2)

	cgplot, SBlim_mag, SBlim_size, /overplot, color='white', linestyle=2

	print,''
	print,'Area covered: ', n_elements(where(h_allz gt 0.)) * xbinsize * ybinsize
	print,''

if keyword_set(stop) then stop

end
