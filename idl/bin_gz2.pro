
;+
; NAME:
;       
;	BIN_GZ2
;
; PURPOSE:
;
;	Use HIST_ND to bin the GZ2 galaxies in MR, R50, and redshift.
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
;       Written by K. Willett                Dec 12
;-

pro bin_gz2, stop=stop

device,retain=2

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = gz2dir+'fits/'
savdir = gz2dir+'sav/'
figsdir = gz2dir+'gz2dropbox/figures/'
gz2tablefile = fitsdir+'gz2main_table_sample.fits'

	;gz2 = mrdfits(gz2tablefile, 1, /silent)
	;save,gz2,filename=savdir+'gz2main_table_sample.sav'

restore,savdir+'gz2main_table_sample.sav'

tagnames=tag_names(gz2)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [where(temparr)]
nr = n_elements(wfind)

task_count_ind_arr = [3,5,7,9,13,15,18,25,28,31,37]
wcind = ([0,wfind])[task_count_ind_arr] + 1
nr_arr = [3,2,2,2,4,2,3,7,3,3,6]
tc_all = [0]
for i=0,n_elements(nr_arr) - 1 do tc_all = [tc_all,intarr(nr_arr[i])+wcind[i]]
tc_all = tc_all[1:n_elements(tc_all)-1]

ngz = n_elements(gz2)

r50_kpc = gz2.petror50_r_kpc
mr = gz2.petromag_mr    
redshift = gz2.redshift
r90 = gz2.PETROR90_R

rad_in_arcsec = 180.*60.*60./!pi
dmod_arr = fltarr(ngz)
ang_scale_arr = fltarr(ngz)
min_classifications = 30

;for i=0,ngz - 1 do begin
;    if finite(redshift[i]) then begin
;        red, z=redshift[i],/silent
;        dmod_arr[i] = dmodulus(redshift[i])
;        ang_scale_arr[i] = dangular(redshift[i],/kpc) / rad_in_arcsec
;    endif
;endfor
;save, redshift, dmod_arr, ang_scale_arr, savdir+'gz2_main_cosmology.sav'
restore,savdir+'gz2_main_cosmology.sav'

mag_padding = 1.0
SB_padding = 1.0
surfacebrightness_app_lim = 23.0
absmag_lim = 17.0 - dmod_arr - mag_padding
r90_kpc = r90 * ang_scale_arr

maskind = where($
            (strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra') and $
            finite(redshift) and $
            mr le (absmag_lim) and $
            r90_kpc ge (3.0 * ang_scale_arr) and $
            (mr + dmod_arr + 2.5*alog10(6.283185d*(r50_kpc/ang_scale_arr)^2)) $
                le (surfacebrightness_app_lim - SB_padding) and $
            gz2.t01_smooth_or_features_total_weight ge min_classifications $
            ,count_om)

om = gz2[maskind]

vo = transpose([[redshift[maskind]],[mr[maskind]],[r50_kpc[maskind]]])
mins    =[0.01,-24.  , 0. ]
maxs    =[0.26,-16.  ,15. ]
binsizes=[0.01,  0.25, 0.5]

timestart = systime(1)

ho = hist_nd(vo,binsizes,min=mins,max=maxs-2*binsizes,reverse_indices=ri)

el_arr = ho*0d
sp_arr = ho*0d
size_ho = size(ho)

for i=0,size_ho[-1] - 1 do $
    if ri[i] ne ri[i+1] then begin
        el_arr[i] = total(om[ri[ri[i]:ri[i+1]-1]].t01_smooth_or_features_a01_smooth_weighted_fraction)
        sp_arr[i] = total(om[ri[ri[i]:ri[i+1]-1]].t01_smooth_or_features_a02_features_or_disk_weighted_fraction)
    endif

!p.multi=[0,5,6]
cgloadct,39

SBapp = 23.0
appmag_lim = 17.0
SBlim_size = findgen(150)/10.
zarr = fillarr(binsizes[0],mins[0],maxs[0])

; Plot the el/sp ratio at each redshift slice

for i=0,size_ho[1]-2 do begin
    zcenter = zarr[i] + binsizes[0]
    red, z = zcenter, /silent
    SBlim_mag = (SBapp - dmodulus(zcenter)- 2.5*alog10(6.283185*(SBlim_size/(dangular(zcenter,/kpc)/rad_in_arcsec))^2))
    absmag_lim = appmag_lim - dmodulus(zcenter)
    size_1arcsec = dangular(zcenter,/kpc)/rad_in_arcsec

    cgimage, alog10(el_arr[i,*,*] / sp_arr[i,*,*]), /axes,/scale,/interpolate,position=[0.1,0.1,0.9,0.9], xrange=[mins[1],maxs[1]],yrange=[mins[2],maxs[2]], charsize=3, minvalue=-2, maxvalue=2
    cgplots, SBlim_mag,SBlim_size, color='white', linestyle=2
    cgplots, replicate(absmag_lim,2), !y.crange, color='white', linestyle=2
    cgplots, !x.crange, replicate(size_1arcsec,2), color='white', linestyle=2

endfor
cgcolorbar,position=[0.15,0.05,0.95,0.10],range=[-2,2],color='black',charsize=3,charthick=2

;ind=[i+nx*(j+ny*k)]
;ri[ri[ind]:ri[ind+1]-1]

timeend = systime(1)

;print,timeend-timestart


; Load in the cube from the Python program

pycube = mrdfits(fitsdir+'task01_binned_rawlikelihood.fits',0)    

if keyword_set(stop) then stop

end
