
;+
; NAME:
;       
;	BIN_GZ2_FROM_PYTHON
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
;	Intended to be run by the galaxyzoo2 Python module, taking advantage of the 
;   REVERSE_INDICES keyword in the HIST_ND IDL function (critical for efficient summing of
;   the likelihood values in each z,M,R50 bin). Doing it in straight Python takes hours per task, 
;   vs. minutes in IDL. 
;
;   Python module takes this script and creates a temporary dummy script named TEMP_BIN_GZ2_FROM_PYTHON, 
;   adding lines at the top that define the task parameters and bin sizes for the script. This temporary
;   script is what is actually run and produces the FITS file; executing this file will result in an error
;   since necessary variables are not defined. 
;
; REVISION HISTORY
;       Written by K. Willett                Dec 12
;-

gz2dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/'
fitsdir = gz2dir+'fits/'

; Read in the masked data created by Python module

gz2 = mrdfits(fitsdir+vardef+'_data_for_idl.fits',1,/silent)
tagnames=tag_names(gz2)

; Create array of bins in redshift, abs mag, and size space

vo = transpose([[gz2.redshift],[gz2.petromag_mr],[gz2.petror50_r_kpc]])

; Define bin size
mins     = float([zmin, magmin, sizemin])
maxs     = float([zmax, magmax, sizemax])
binsizes = float([zstep, magstep, sizestep])

; Bin the GZ2 data in three dimensions
ho = hist_nd(vo,binsizes,min=mins,max=maxs,reverse_indices=ri)

task_names_wf_arr = strsplit(strcompress(task_names_wf,/remove_all),',',/ex)
nvar = n_elements(task_names_wf_arr)

; Create 4-D array, adding dimension for the number of responses per GZ2 task
size_ho = size(ho)
hcube = dblarr(nvar,size_ho[1],size_ho[2],size_ho[3])

; Use the reverse indices in HIST_ND to sum the vote fractions in each bin, 
; yielding the raw likelihoods of each task for every value of z, M_R, R_50

for j=0,nvar - 1 do begin
    hslice_temp = ho*0d
    temp_wf = where(tagnames eq strupcase(task_names_wf_arr[j]))
    for i=0,size_ho[-1] - 1 do $
        if ri[i] ne ri[i+1] then begin
            hslice_temp[i] = total(gz2[ri[ri[i]:ri[i+1]-1]].(temp_wf))
        endif
    hcube[j,*,*,*] = hslice_temp
endfor

; Define the edges and centers for the 3-D binned data
edges_redshift = dindgen((zmax + zstep - zmin)/zstep + 1)*zstep + zmin
edges_mag = dindgen((magmax + magstep - magmin)/magstep + 1)*magstep + magmin
edges_size = dindgen((sizemax + sizestep - sizemin)/sizestep + 1)*sizestep + sizemin

centers_redshift = edges_redshift[0:-2] + (edges_redshift[1]-edges_redshift[0])/2.
centers_mag = edges_mag[0:-2] + (edges_mag[1]-edges_mag[0])/2.
centers_size = edges_size[0:-2] + (edges_size[1]-edges_size[0])/2.

; Create structure of the edges and centers
a = {edges_redshift:edges_redshift, $
    edges_mag:edges_mag, $
    edges_size:edges_size, $
    centers_redshift:centers_redshift, $
    centers_mag:centers_mag, $
    centers_size:centers_size}

; Write the binned likelihoods to FITS file as primary HDU, 
; add the bin edges and centers as binary tables in secondary HDU. 
; Data can now be read back into Python.

mwrfits, transpose(hcube), fitsdir+vardef+'_idlbinned.fits', /create,/silent
mwrfits, a, fitsdir+vardef+'_idlbinned.fits',/silent

; Create temp file to indicate that process has finished

spawn, 'touch '+fitsdir+'/idlfilecreated'

end
