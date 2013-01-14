
;+
; NAME:
;       
;	GZ2_STRIPE82
;
; PURPOSE:
;
;	Plot the morphological attributes for single- and coadd-depth images from Stripe 82 in SDSS
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

pro gz2_stripe82_coadd1, ps=ps

; Import Stripe 82 data

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
gz2tablefile = gz2dir+'gz2table.fits'
gz2samplefile = gz2dir+'gz2sample.fits'

gz2_n_c1_coadd_file = gz2dir+'n_c1_coadd_willettk.fit'
gz2_n_c1_normal_file = gz2dir+'n_c1_normal_willettk.fit'
gz2_n_c1_spec_file = gz2dir+'n_c1_redshifts_willettk.fit'
gz2_n_c1_sample_file = gz2dir+'n_c1_sample_willettk.fit'

;gz2_n_c2_coadd_file = gz2dir+'n_c2_coadd_willettk.fit'
;gz2_n_c2_normal_file = gz2dir+'n_c2_normal_willettk.fit'
;gz2_n_c2_spec_file = gz2dir+'n_c2_redshifts_willettk.fit'
;gz2_n_c2_sample_file = gz2dir+'n_c2_sample_willettk.fit'

gz2_n_c2_coadd_file = gz2dir+'TEST_n_c2_coadd_willettk.fit'
gz2_n_c2_normal_file = gz2dir+'TEST_n_c2_normal_willettk.fit'
gz2_n_c2_spec_file = gz2dir+'TEST_n_c2_spec_willettk.fit'
gz2_n_c2_sample_file = gz2dir+'TEST_n_c2_sample_willettk.fit'

nc1_coadd = mrdfits(gz2_n_c1_coadd_file ,1)
nc1_normal = mrdfits(gz2_n_c1_normal_file ,1)
nc1_spec = mrdfits(gz2_n_c1_spec_file ,1)
nc1_sample = mrdfits(gz2_n_c1_sample_file ,1)

nc2_coadd = mrdfits(gz2_n_c2_coadd_file ,1)
nc2_normal = mrdfits(gz2_n_c2_normal_file ,1)
nc2_spec = mrdfits(gz2_n_c2_spec_file ,1)
nc2_sample = mrdfits(gz2_n_c2_sample_file ,1)

; ObjIDs for the two coadd samples are the same (although they will have different asset_IDs), since they
; are the same object in the Sloan database. However, I will have to do separate SQL queries in CasJobs
; to pull out the GZ2 data since they are separate assets. 

; DO NOT try to match anything in the full GZ2 table in IDL. Way too large for the memory. 

; Hangs up at 1853. What a random number. Everything scales fine before then (average 1.3 sec/match)
; Not a function of that specific integer - same problem if I start at higher index and try to match more than 1852 elements. Puzzling. Guess I could do it 1000 at a time. 

tagnames=tag_names(nc1_coadd)
nt = n_elements(tagnames)

temparr=intarr(nt) & for i=0,nt-1 do begin & temp = tagnames[i] & templen = strlen(temp) & if templen ge 18 then begin & temp2 = strmid(temp,templen-17,17) & if temp2 eq 'WEIGHTED_FRACTION' then temparr[i]=1 & endif & endfor

nbyn = 6
!p.multi=[0,nbyn,nbyn]

count=5 

; Scatter plot: coadd vs. normal

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd1_normal_scatter.ps',/color, /quiet
	cs=1.2
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

for i=0,nbyn^2-1 do begin 

	tagind_f = (where(temparr))[i]
	tagind_c = (where(temparr))[i] - 2

	countind = where(nc1_coadd.(tagind_c) ge count and nc1_normal.(tagind_c) ge count)
	coadd1_val = nc1_coadd[countind].(tagind_f)
	normal_val = nc1_normal[countind].(tagind_f)
	z = nc1_spec[countind].z

	cgplot, coadd1_val, normal_val, $
		psym=3, $
		xr=[0,1], /xstyle, $
		yr=[0,1], /ystyle, $
		title='Task '+strmid(tagnames[tagind_f],1,2)+', Answer '+string(i+1,format='(i02)'), $
		xtitle='coadd', ytitle='normal'
	cgplot, indgen(11)/10., indgen(11)/10., color='red' ,/overplot
endfor

if keyword_set(ps) then ps_end

; Histogram: coadd - normal

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd1_normal_histogram.ps',/color, /quiet
	cs=1.2
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

for i=0,nbyn^2-1 do begin 

	tagind_f = (where(temparr))[i]
	tagind_c = (where(temparr))[i] - 2

	countind = where(nc1_coadd.(tagind_c) ge count and nc1_normal.(tagind_c) ge count)
	coadd1_val = nc1_coadd[countind].(tagind_f)
	normal_val = nc1_normal[countind].(tagind_f)
	z = nc1_spec[countind].z

	cghistoplot, coadd1_val - normal_val, $
		title='Task '+strmid(tagnames[tagind_f],1,2)+$
		', Answer '+string(i+1,format='(i02)'), $
;		+' '+greek('mu')+'='+string(mean(coadd1_val - normal_val),format='(f5.2)')+', '+$
;		greek('sigma')+'='+string(stddev(coadd1_val - normal_val),format='(f5.2)'),$
		datacolor='blue',$
		/freq,$
		xr=[-1,1],/xstyle,$
		xtitle='coadd1 - normal'
	cgplots, [0,0], !y.crange, color='red' 

endfor

if keyword_set(ps) then ps_end

; Scatter plot: redshift vs. (coadd - normal)

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_redshift_deltaf_coadd1.ps',/color, /quiet
	cs=1.2
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

for i=0,nbyn^2-1 do begin 
	tagind_f = (where(temparr))[i]
	tagind_c = (where(temparr))[i] - 2

	countind = where(nc1_coadd.(tagind_c) ge count and nc1_normal.(tagind_c) gt count)
	coadd1_val = nc1_coadd[countind].(tagind_f)
	normal_val = nc1_normal[countind].(tagind_f)
	z = nc1_spec[countind].z

	cgplot, z, coadd1_val - normal_val, $
		title='Task '+strmid(tagnames[tagind_f],1,2)+$
		', Answer '+string(i+1,format='(i02)'), $
		;+' '+greek('mu')+'='+string(mean(coadd1_val - normal_val),format='(f5.2)')+', '+$
		;greek('sigma')+'='+string(stddev(coadd1_val - normal_val),format='(f5.2)'),$
		xtitle='z',$
		ytitle='coadd1 - normal',$
		psym=3,$
		xr=[-0.05,0.30], $
		yr=[-1.1,1.1], $
		/xstyle, /ystyle

	medarr=fltarr(26) 
	for j=0,25 do begin 
		zlower = j/100.
		zupper = (j+1)/100.
		arr = (coadd1_val - normal_val)[where(z gt zlower and z le zupper,zcount)]

		if zcount eq 0 then begin
			if arr[0] eq -1 then medarr[j] = !values.f_nan 
		endif else if zcount eq 1 then begin
			medarr[j] = arr[0] 
		endif else medarr[j] = mean(arr)

	endfor 
	cgplot,indgen(26)/100.,medarr,color='red',/overplot
endfor

if keyword_set(ps) then ps_end

; Plot the weighted fractions for Task 01 for normal, coadd1, and coadd2 in Stripe 82

!p.multi=[0,2,1]

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_redshift_bothcoadd_task01.ps',/color, /quiet
	cs=1.2
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

a01_ind_f = (where(temparr))[0]
a01_ind_c = (where(temparr))[0]-2
a02_ind_f = (where(temparr))[1]
a02_ind_c = (where(temparr))[1]-2
a03_ind_f = (where(temparr))[2]
a03_ind_c = (where(temparr))[2]-2

z = nc1_spec.z

; Thresholded likelihoods

deltaz=0.25
nbinz = 25

threshfrac_el_coadd1=fltarr(nbinz) 
threshfrac_sp_coadd1=fltarr(nbinz) 
threshfrac_ar_coadd1=fltarr(nbinz) 
threshfrac_uc_coadd1=fltarr(nbinz) 
threshfrac_el_coadd2=fltarr(nbinz) 
threshfrac_sp_coadd2=fltarr(nbinz) 
threshfrac_ar_coadd2=fltarr(nbinz) 
threshfrac_uc_coadd2=fltarr(nbinz) 
threshfrac_el_normal=fltarr(nbinz) 
threshfrac_sp_normal=fltarr(nbinz) 
threshfrac_ar_normal=fltarr(nbinz) 
threshfrac_uc_normal=fltarr(nbinz) 

	for j=0,nbinz-1 do begin 
		zlower = j/(nbinz/deltaz)
		zupper = (j+1)/(nbinz/deltaz)
		zind = where(z gt zlower and z le zupper, zcount)

		; Coadd 1

		temp = where(nc1_coadd[zind].(a01_ind_f) gt 0.8,ct_el_coadd1)
		ca_el_coadd1 = n_elements(nc1_coadd[zind].(a01_ind_f))
		if ca_el_coadd1 eq 0 then threshfrac_el_coadd1[j] = !values.f_nan else $
			threshfrac_el_coadd1[j] = float(ct_el_coadd1)/float(ca_el_coadd1)
		temp = where(nc1_coadd[zind].(a02_ind_f) gt 0.8,ct_sp_coadd1)
		ca_sp_coadd1 = n_elements(nc1_coadd[zind].(a02_ind_f))
		if ca_sp_coadd1 eq 0 then threshfrac_sp_coadd1[j] = !values.f_nan else $
			threshfrac_sp_coadd1[j] = float(ct_sp_coadd1)/float(ca_sp_coadd1)
		temp = where(nc1_coadd[zind].(a03_ind_f) gt 0.8,ct_ar_coadd1)
		ca_ar_coadd1 = n_elements(nc1_coadd[zind].(a03_ind_f))
		if ca_ar_coadd1 eq 0 then threshfrac_ar_coadd1[j] = !values.f_nan else $
			threshfrac_ar_coadd1[j] = float(ct_ar_coadd1)/float(ca_ar_coadd1)
		temp = where($
			nc1_coadd[zind].(a01_ind_f) lt 0.8 and $
			nc1_coadd[zind].(a02_ind_f) lt 0.8 and $
			nc1_coadd[zind].(a03_ind_f) lt 0.8, $
			ct_uc_coadd1)
		ca_uc_coadd1 = n_elements(nc1_coadd[zind])
		if ca_uc_coadd1 eq 0 then threshfrac_uc_coadd1[j] = !values.f_nan else $
			threshfrac_uc_coadd1[j] = float(ct_uc_coadd1)/float(ca_uc_coadd1)

		; Coadd 2

		temp = where(nc2_coadd[zind].(a01_ind_f) gt 0.8,ct_el_coadd2)
		ca_el_coadd2 = n_elements(nc2_coadd[zind].(a01_ind_f))
		if ca_el_coadd2 eq 0 then threshfrac_el_coadd2[j] = !values.f_nan else $
			threshfrac_el_coadd2[j] = float(ct_el_coadd2)/float(ca_el_coadd2)
		temp = where(nc2_coadd[zind].(a02_ind_f) gt 0.8,ct_sp_coadd2)
		ca_sp_coadd2 = n_elements(nc2_coadd[zind].(a02_ind_f))
		if ca_sp_coadd2 eq 0 then threshfrac_sp_coadd2[j] = !values.f_nan else $
			threshfrac_sp_coadd2[j] = float(ct_sp_coadd2)/float(ca_sp_coadd2)
		temp = where(nc2_coadd[zind].(a03_ind_f) gt 0.8,ct_ar_coadd2)
		ca_ar_coadd2 = n_elements(nc2_coadd[zind].(a03_ind_f))
		if ca_ar_coadd2 eq 0 then threshfrac_ar_coadd2[j] = !values.f_nan else $
			threshfrac_ar_coadd2[j] = float(ct_ar_coadd2)/float(ca_ar_coadd2)
		temp = where($
			nc2_coadd[zind].(a01_ind_f) lt 0.8 and $
			nc2_coadd[zind].(a02_ind_f) lt 0.8 and $
			nc2_coadd[zind].(a03_ind_f) lt 0.8, $
			ct_uc_coadd2)
		ca_uc_coadd2 = n_elements(nc2_coadd[zind])
		if ca_uc_coadd2 eq 0 then threshfrac_uc_coadd2[j] = !values.f_nan else $
			threshfrac_uc_coadd2[j] = float(ct_uc_coadd2)/float(ca_uc_coadd2)

		; Normal

		temp = where(nc1_normal[zind].(a01_ind_f) gt 0.8,ct_el_normal)
		ca_el_normal = n_elements(nc1_normal[zind].(a01_ind_f))
		if ca_el_normal eq 0 then threshfrac_el_normal[j] = !values.f_nan else $
			threshfrac_el_normal[j] = float(ct_el_normal)/float(ca_el_normal)
		temp = where(nc1_normal[zind].(a02_ind_f) gt 0.8,ct_sp_normal)
		ca_sp_normal = n_elements(nc1_normal[zind].(a02_ind_f))
		if ca_sp_normal eq 0 then threshfrac_sp_normal[j] = !values.f_nan else $
			threshfrac_sp_normal[j] = float(ct_sp_normal)/float(ca_sp_normal)
		temp = where(nc1_normal[zind].(a03_ind_f) gt 0.8,ct_ar_normal)
		ca_ar_normal = n_elements(nc1_normal[zind].(a03_ind_f))
		if ca_ar_normal eq 0 then threshfrac_ar_normal[j] = !values.f_nan else $
			threshfrac_ar_normal[j] = float(ct_ar_normal)/float(ca_ar_normal)
		temp = where($
			nc1_normal[zind].(a01_ind_f) lt 0.8 and $
			nc1_normal[zind].(a02_ind_f) lt 0.8 and $
			nc1_normal[zind].(a03_ind_f) lt 0.8, $
			ct_uc_normal)
		ca_uc_normal = n_elements(nc1_normal[zind])
		if ca_uc_normal eq 0 then threshfrac_uc_normal[j] = !values.f_nan else $
			threshfrac_uc_normal[j] = float(ct_uc_normal)/float(ca_uc_normal)


	endfor 

zarr = indgen(nbinz)/(nbinz/deltaz)

cs = 1
labelsize=1
nthick=5
cgplot, zarr, threshfrac_el_normal, $
	thick=nthick, $
	linestyle=0, $
	charsize=cs, $
	xr=[0,0.25], $
	yr=[0,1], $
	color='red', $
	xtitle='Redshift', $
	ytitle='Thresholded fraction (p>0.8)', $
	title='Task 01 - GZ2 Stripe 82'

cgplot, zarr, threshfrac_sp_normal, /overplot, thick = nthick, linestyle=0, color='blue'
cgplot, zarr, threshfrac_ar_normal, /overplot, thick = nthick, linestyle=0, color='green'
cgplot, zarr, threshfrac_uc_normal, /overplot, thick = nthick, linestyle=0, color='orange'

cgplot, zarr, threshfrac_el_coadd1, /overplot, thick = 2, linestyle=1, color='red'
cgplot, zarr, threshfrac_sp_coadd1, /overplot, thick = 2, linestyle=1, color='blue'
cgplot, zarr, threshfrac_ar_coadd1, /overplot, thick = 2, linestyle=1, color='green'
cgplot, zarr, threshfrac_uc_coadd1, /overplot, thick = 2, linestyle=1, color='orange'

cgplot, zarr, threshfrac_el_coadd2, /overplot, thick = 2, linestyle=2, color='red'
cgplot, zarr, threshfrac_sp_coadd2, /overplot, thick = 2, linestyle=2, color='blue'
cgplot, zarr, threshfrac_ar_coadd2, /overplot, thick = 2, linestyle=2, color='green'
cgplot, zarr, threshfrac_uc_coadd2, /overplot, thick = 2, linestyle=2, color='orange'

;cgerrplot, zarr, threshfrac_el_coadd1 - threshfrac_el_coadd1_err, threshfrac_el_coadd1 + threshfrac_el_coadd1_err, color='red'
;cgplot, z, nc1_coadd.(a01_ind_f), psym=3,color='red',/overplot
;cgerrplot, zarr, threshfrac_sp_coadd1 - threshfrac_sp_coadd1_err, threshfrac_sp_coadd1 + threshfrac_sp_coadd1_err, color='blue'
;cgplot, z, nc1_coadd.(a02_ind_f), psym=3,color='blue',/overplot
;cgerrplot, zarr, threshfrac_ar_coadd1 - threshfrac_ar_coadd1_err, threshfrac_ar_coadd1 + threshfrac_ar_coadd1_err, color='green'
;cgplot, z, nc1_coadd.(a03_ind_f), psym=3,color='green',/overplot

al_legend, /top, /left, ['Smooth','Features/disk','Artifact/star','Unclassified'],color=['red','blue','green','orange'],psym=28,charsize=labelsize
;al_legend, /top, /right, ['Coadd1','Normal'],color=['black','black'],linestyle=[0,2], thick=[3,1],charsize=labelsize

; Raw likelihoods

	frac_el_coadd1=fltarr(nbinz) 
	frac_el_coadd1_err=fltarr(nbinz) 
	frac_sp_coadd1=fltarr(nbinz) 
	frac_sp_coadd1_err=fltarr(nbinz) 
	frac_ar_coadd1=fltarr(nbinz) 
	frac_ar_coadd1_err=fltarr(nbinz) 
	frac_el_coadd2=fltarr(nbinz) 
	frac_el_coadd2_err=fltarr(nbinz) 
	frac_sp_coadd2=fltarr(nbinz) 
	frac_sp_coadd2_err=fltarr(nbinz) 
	frac_ar_coadd2=fltarr(nbinz) 
	frac_ar_coadd2_err=fltarr(nbinz) 
	frac_el_normal=fltarr(nbinz) 
	frac_el_normal_err=fltarr(nbinz) 
	frac_sp_normal=fltarr(nbinz) 
	frac_sp_normal_err=fltarr(nbinz) 
	frac_ar_normal=fltarr(nbinz) 
	frac_ar_normal_err=fltarr(nbinz) 

	for j=0,nbinz-1 do begin 

		zlower = j/(nbinz/deltaz)
		zupper = (j+1)/(nbinz/deltaz)
		zind = where(z gt zlower and z le zupper, zcount)

		;arr = nc1_coadd[where(z gt zlower and z le zupper,zcount)].(a01_ind_f)
		arr_el_coadd1 = nc1_coadd[zind].(a01_ind_f)
		arr_sp_coadd1 = nc1_coadd[zind].(a02_ind_f)
		arr_ar_coadd1 = nc1_coadd[zind].(a03_ind_f)
		arr_el_coadd2 = nc2_coadd[zind].(a01_ind_f)
		arr_sp_coadd2 = nc2_coadd[zind].(a02_ind_f)
		arr_ar_coadd2 = nc2_coadd[zind].(a03_ind_f)
		arr_el_normal = nc1_normal[zind].(a01_ind_f)
		arr_sp_normal = nc1_normal[zind].(a02_ind_f)
		arr_ar_normal = nc1_normal[zind].(a03_ind_f)

		if zcount eq 0 then begin
			frac_el_coadd1[j] = !values.f_nan 
			frac_el_coadd1_err[j] = !values.f_nan 
			frac_sp_coadd1[j] = !values.f_nan 
			frac_sp_coadd1_err[j] = !values.f_nan 
			frac_ar_coadd1[j] = !values.f_nan 
			frac_ar_coadd1_err[j] = !values.f_nan 
			frac_el_coadd2[j] = !values.f_nan 
			frac_el_coadd2_err[j] = !values.f_nan 
			frac_sp_coadd2[j] = !values.f_nan 
			frac_sp_coadd2_err[j] = !values.f_nan 
			frac_ar_coadd2[j] = !values.f_nan 
			frac_ar_coadd2_err[j] = !values.f_nan 
			frac_el_normal[j] = !values.f_nan 
			frac_el_normal_err[j] = !values.f_nan 
			frac_sp_normal[j] = !values.f_nan 
			frac_sp_normal_err[j] = !values.f_nan 
			frac_ar_normal[j] = !values.f_nan 
			frac_ar_normal_err[j] = !values.f_nan 
		endif else if zcount eq 1 then begin
			frac_el_coadd1[j] = arr_el_coadd1[0] 
			frac_el_coadd1_err[j] = arr_el_coadd1[0] 
			frac_sp_coadd1[j] = arr_sp_coadd1[0]
			frac_sp_coadd1_err[j] = arr_sp_coadd1[0]
			frac_ar_coadd1[j] = arr_ar_coadd1[0]
			frac_ar_coadd1_err[j] = arr_ar_coadd1[0]
			frac_el_coadd2[j] = arr_el_coadd2[0] 
			frac_el_coadd2_err[j] = arr_el_coadd2[0] 
			frac_sp_coadd2[j] = arr_sp_coadd2[0]
			frac_sp_coadd2_err[j] = arr_sp_coadd2[0]
			frac_ar_coadd2[j] = arr_ar_coadd2[0]
			frac_ar_coadd2_err[j] = arr_ar_coadd2[0]
			frac_el_normal[j] = arr_el_normal[0] 
			frac_el_normal_err[j] = arr_el_normal[0] 
			frac_sp_normal[j] = arr_sp_normal[0]
			frac_sp_normal_err[j] = arr_sp_normal[0]
			frac_ar_normal[j] = arr_ar_normal[0]
			frac_ar_normal_err[j] = arr_ar_normal[0]
		endif else begin
			frac_el_coadd1[j] = mean(arr_el_coadd1)
			frac_el_coadd1_err[j] = stddev(arr_el_coadd1)
			frac_sp_coadd1[j] = mean(arr_sp_coadd1)
			frac_sp_coadd1_err[j] = stddev(arr_sp_coadd1)
			frac_ar_coadd1[j] = mean(arr_ar_coadd1)
			frac_ar_coadd1_err[j] = stddev(arr_ar_coadd1)
			frac_el_coadd2[j] = mean(arr_el_coadd2)
			frac_el_coadd2_err[j] = stddev(arr_el_coadd2)
			frac_sp_coadd2[j] = mean(arr_sp_coadd2)
			frac_sp_coadd2_err[j] = stddev(arr_sp_coadd2)
			frac_ar_coadd2[j] = mean(arr_ar_coadd2)
			frac_ar_coadd2_err[j] = stddev(arr_ar_coadd2)
			frac_el_normal[j] = mean(arr_el_normal)
			frac_el_normal_err[j] = stddev(arr_el_normal)
			frac_sp_normal[j] = mean(arr_sp_normal)
			frac_sp_normal_err[j] = stddev(arr_sp_normal)
			frac_ar_normal[j] = mean(arr_ar_normal)
			frac_ar_normal_err[j] = stddev(arr_ar_normal)
		endelse

	endfor 

zarr = indgen(nbinz)/(nbinz/deltaz)

cgplot, zarr, frac_el_normal, $
	thick=nthick, $
	linestyle=0, $
	charsize=cs, $
	xr=[0,0.25], $
	yr=[0,1], $
	color='red', $
	xtitle='Redshift', $
	ytitle='Raw likelihoods', $
	title='Task 01 - GZ2 Stripe 82'

cgplot, /overplot, zarr, frac_sp_normal, linestyle=0, thick=nthick, color='blue'
cgplot, /overplot, zarr, frac_ar_normal, linestyle=0, thick=nthick, color='green'

cgplot, zarr, frac_el_coadd1, /overplot, linestyle=1, thick=2, color='red'
cgplot, zarr, frac_sp_coadd1, /overplot, linestyle=1, thick=2, color='blue'
cgplot, zarr, frac_ar_coadd1, /overplot, linestyle=1, thick=2, color='green'
                                                      
cgplot, zarr, frac_el_coadd2, /overplot, linestyle=2, thick=2, color='red'
cgplot, zarr, frac_sp_coadd2, /overplot, linestyle=2, thick=2, color='blue'
cgplot, zarr, frac_ar_coadd2, /overplot, linestyle=2, thick=2, color='green'

;al_legend, /top, /left, ['Smooth','Features/disk','Artifact/star'],color=['red','blue','green'],psym=28,charsize=labelsize
al_legend, /top, /right, ['Normal','Coadd1','Coadd2'],color=['black','black','black'],linestyle=[0,1,2], thick=[3,2,2],charsize=labelsize

if keyword_set(ps) then ps_end

; Plot difference in coadd2 - normal depth data as function of r-band surface brightness

llim = 10

nel_coadd2 = nc2_coadd.((where(temparr))[0]-3)
nsp_coadd2 = nc2_coadd.((where(temparr))[1]-3)
nratio_coadd2 = float(nel_coadd2) / float(nsp_coadd2)

nel_normal = nc2_normal.((where(temparr))[0]-3)
nsp_normal = nc2_normal.((where(temparr))[1]-3)
nratio_normal = float(nel_normal) / float(nsp_normal)

fel_coadd2 = nc2_coadd.((where(temparr))[0])
fel_normal = nc2_normal.((where(temparr))[0])
fsp_coadd2 = nc2_coadd.((where(temparr))[1])
fsp_normal = nc2_normal.((where(temparr))[1])

ntot_coadd2 = float(nel_coadd2 + nsp_coadd2)
ntot_normal = float(nel_normal + nsp_normal)

linexp = 'p[0] + p[1]*x'

nrind = where(finite(alog10(nratio_coadd2)) and finite(alog10(nratio_normal)))

nratio_alog_rmag_err = sqrt((2. * nsp_coadd2[nrind] * nsp_normal[nrind] + nsp_coadd2[nrind] + nsp_normal[nrind])/(nsp_coadd2[nrind] + nsp_normal[nrind]))

fsp_coadd2_err = nsp_coadd2/ntot_coadd2^2 * sqrt(nel_coadd2 + nsp_coadd2/ntot_coadd2^2)
fsp_normal_err = nsp_normal/ntot_normal^2 * sqrt(nel_normal + nsp_normal/ntot_normal^2)
deltafsp_err = sqrt(fsp_coadd2_err^2 + fsp_normal_err^2)
deltafsp_err[where(deltafsp_err eq 0.)] = 1.

nratio_coadd2_variance = nratio_coadd2 * (1. + nratio_coadd2)
nratio_normal_variance = nratio_normal * (1. + nratio_normal)
nratio_err = sqrt(nratio_coadd2_variance + nratio_normal_variance)

fel_coadd2_err = nel_coadd2/ntot_coadd2^2 * sqrt(nsp_coadd2 + nel_coadd2/ntot_coadd2^2)
fel_normal_err = nel_normal/ntot_normal^2 * sqrt(nsp_normal + nel_normal/ntot_normal^2)
deltafel_err = sqrt(fel_coadd2_err^2 + fel_normal_err^2)
deltafel_err[where(deltafel_err eq 0.)] = 1.

mu50_r = nc2_sample.mu50_r

	delmu = 0.2
	muarr = fillarr(delmu,18,23)
	nm = n_elements(muarr)

	; Compute mean, median, and linear fits for the samples as functions of redshift

	nratio_alog_arr=fltarr(nm) 
	for j=0,nm-1 do begin 
		mu50_r_ind = where(mu50_r ge muarr[j] and mu50_r lt muarr[j]+delmu,mucount)
		arr = (alog10(nratio_coadd2) - alog10(nratio_normal))[mu50_r_ind]
		if mucount eq 0 then begin
			if arr[0] eq -1 then nratio_alog_arr[j] = !values.f_nan 
		endif else if mucount eq 1 then begin
			nratio_alog_arr[j] = arr[0] 
		endif else nratio_alog_arr[j] = mean(arr,/nan)
	endfor 
	nratio_alog_mu50_r_linfit = mpfitexpr(linexp,$
		mu50_r[nrind], $
		alog10(nratio_coadd2[nrind]) - alog10(nratio_normal[nrind]), $
		sqrt((2. * nsp_coadd2[nrind] * nsp_normal[nrind] + nsp_coadd2[nrind] + nsp_normal[nrind])/(nsp_coadd2[nrind] + nsp_normal[nrind])), $
		perr = nratio_alog_mu50_r_linfit_err)


	fsp_arr=fltarr(nm) 
	for j=0,nm-1 do begin 
		mu50_r_ind = where(mu50_r ge muarr[j] and mu50_r lt muarr[j]+delmu,mucount)
		arr = (fsp_coadd2 - fsp_normal)[mu50_r_ind]
		if mucount eq 0 then begin
			if arr[0] eq -1 then fsp_arr[j] = !values.f_nan 
		endif else if mucount eq 1 then begin
			fsp_arr[j] = arr[0] 
		endif else fsp_arr[j] = mean(arr)
	endfor 
	fsp_mu50_r_linfit = mpfitexpr(linexp,$
		mu50_r, $
		fsp_coadd2 - fsp_normal, $
		fsp_err, $
		perr = fsp_mu50_r_linfit_err)

	nratio_arr=fltarr(nm) 
	for j=0,nm-1 do begin 
		mu50_r_ind = where(mu50_r ge muarr[j] and mu50_r lt muarr[j]+delmu,mucount)
		arr = (nratio_coadd2 - nratio_normal)[mu50_r_ind]
		if mucount eq 0 then begin
			if arr[0] eq -1 then nratio_arr[j] = !values.f_nan 
		endif else if mucount eq 1 then begin
			nratio_arr[j] = arr[0] 
		endif else nratio_arr[j] = mean(arr,/nan)
	endfor 
	nratio_mu50_r_linfit = mpfitexpr(linexp,$
		mu50_r[nrind], $
		nratio_coadd2[nrind] - nratio_normal[nrind], $
		nratio_err[nrind], $
		perr = fel_mu50_r_linfit_err)

	fel_arr=fltarr(nm) 
	for j=0,nm-1 do begin 
		mu50_r_ind = where(mu50_r ge muarr[j] and mu50_r lt muarr[j]+delmu,mucount)
		arr = (fel_coadd2 - fel_normal)[mu50_r_ind]
		if mucount eq 0 then begin
			if arr[0] eq -1 then fel_arr[j] = !values.f_nan 
		endif else if mucount eq 1 then begin
			fel_arr[j] = arr[0] 
		endif else fel_arr[j] = mean(arr)
	endfor 
	fel_mu50_r_linfit = mpfitexpr(linexp,$
		mu50_r, $
		fel_coadd2 - fel_normal, $
		fel_err, $
		perr = fel_mu50_r_linfit_err)

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd1_mu50.eps',/color, /quiet, xs=8, ys=8, /encap
	cs=1.2
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

!p.multi=[0,2,2]
cgplot, mu50_r[nrind], alog10(nratio_coadd2[nrind]) - alog10(nratio_normal[nrind]), $
	psym=3, $
	xr=[18,23], /xstyle, $
	charsize=cs, $
	xtitle=greek('mu')+'!I50,R!N [mag arcsec!E-2!N]', $
	ytitle='log(n!Iel!N/n!Isp!N)!Icoadd2!N - log(n!Iel!N/n!Isp!N)!Inormal!N'

	cgplot,muarr,nratio_alog_arr,color='red',/overplot, thick=3
	cgplot, /overplot, muarr, nratio_alog_mu50_r_linfit[0] + muarr*nratio_alog_mu50_r_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(mu50_r[nrind],alog10(nratio_coadd2[nrind])-alog10(nratio_normal[nrind])),format='(f7.4)'),charsize=cs, 18.5, 2.5, /data, color='blue'

cgplot, mu50_r, fsp_coadd2 - fsp_normal, $
	xr=[18,23], /xstyle, $
	psym=3, $
	charsize=cs, $
	xtitle=greek('mu')+'!I50,R!N [mag arcsec!E-2!N]', $
	ytitle='f!Isp,coadd2!N - f!Isp,normal!N'

	cgplot,muarr,fsp_arr,color='red',/overplot, thick=3
	cgplot, /overplot, muarr, fsp_mu50_r_linfit[0] + muarr*fsp_mu50_r_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(mu50_r, fsp_coadd2 - fsp_normal),format='(f7.4)'),charsize=cs, 18.5, 0.8, /data, color='blue'

cgplot, mu50_r, nratio_coadd2 - nratio_normal, $
	psym=3, $
	xr=[18,23], /xstyle, $
	charsize=cs, $
	xtitle=greek('mu')+'!I50,R!N [mag arcsec!E-2!N]', $
	ytitle='(n!Iel!N/n!Isp!N)!Icoadd2!N - (n!Iel!N/n!Isp!N)!Inormal!N'

	cgplot,muarr,nratio_arr,color='red',/overplot, thick=3
	cgplot, /overplot, muarr, nratio_mu50_r_linfit[0] + muarr*nratio_mu50_r_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(mu50_r[nrind], nratio_coadd2[nrind] - nratio_normal[nrind]),format='(f7.4)'),charsize=cs, 18.5, 30, /data, color='blue'

cgplot, mu50_r, fel_coadd2 - fel_normal, $
	xr=[18,23], /xstyle, $
	psym=3, $
	charsize=cs, $
	xtitle=greek('mu')+'!I50,R!N [mag arcsec!E-2!N]', $
	ytitle='f!Iel,coadd2!N - f!Iel,normal!N'

	cgplot,muarr,fel_arr,color='red',/overplot, thick=3
	cgplot, /overplot, muarr, fel_mu50_r_linfit[0] + muarr*fel_mu50_r_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(mu50_r, fel_coadd2 - fel_normal),format='(f7.4)'),charsize=cs, 18.5, 0.8, /data, color='blue'

if keyword_set(ps) then ps_end

; Plot difference in coadd2 - normal depth data as function of apparent r-band magnitude

rmag = nc2_sample.petromag_r

	delmag = 0.2
	rmagarr = fillarr(delmag,11,18)
	nm = n_elements(rmagarr)

	nratio_alog_rmag_mean=fltarr(nm) 
	nratio_alog_rmag_median=fltarr(nm) 
	for j=0,nm-1 do begin 
		rmag_ind = where(rmag ge rmagarr[j] and rmag lt rmagarr[j]+delmag,rmagcount)
		arr = (alog10(nratio_coadd2) - alog10(nratio_normal))[rmag_ind]
		if rmagcount eq 0 then begin
			nratio_alog_rmag_mean[j] = !values.f_nan 
			nratio_alog_rmag_median[j] = !values.f_nan 
		endif else if rmagcount eq 1 then begin
			nratio_alog_rmag_mean[j] = arr[0] 
			nratio_alog_rmag_median[j] = arr[0] 
		endif else begin
			nratio_alog_rmag_mean[j] = mean(arr,/nan)
			nratio_alog_rmag_median[j] = median(arr)
		endelse
	endfor 
	nratio_alog_rmag_linfit = mpfitexpr(linexp,$
		rmag[nrind], $
		alog10(nratio_coadd2[nrind]) - alog10(nratio_normal[nrind]), $
		sqrt((2. * nsp_coadd2[nrind] * nsp_normal[nrind] + nsp_coadd2[nrind] + nsp_normal[nrind])/(nsp_coadd2[nrind] + nsp_normal[nrind])), $
		perr = nratio_alog_rmag_linfit_err)

	fsp_rmag_mean=fltarr(nm) 
	fsp_rmag_median=fltarr(nm) 
	for j=0,nm-1 do begin 
		rmag_ind = where(rmag ge rmagarr[j] and rmag lt rmagarr[j]+delmag,rmagcount)
		arr = (fsp_coadd2 - fsp_normal)[rmag_ind]
		if rmagcount eq 0 then begin
			fsp_rmag_mean[j] = !values.f_nan 
			fsp_rmag_median[j] = !values.f_nan 
		endif else if rmagcount eq 1 then begin
			fsp_rmag_mean[j] = arr[0] 
			fsp_rmag_median[j] = arr[0] 
		endif else begin
			fsp_rmag_mean[j] = mean(arr)
			fsp_rmag_median[j] = median(arr)
		endelse
	endfor 
	fsp_rmag_linfit = mpfitexpr(linexp,$
		rmag, $
		fsp_coadd2 - fsp_normal, $
		fsp_err, $
		perr = fsp_rmag_linfit_err)


	nratio_rmag_mean=fltarr(nm) 
	nratio_rmag_median=fltarr(nm) 
	for j=0,nm-1 do begin 
		rmag_ind = where(rmag ge rmagarr[j] and rmag lt rmagarr[j]+delmag,rmagcount)
		arr = (nratio_coadd2 - nratio_normal)[rmag_ind]
		if rmagcount eq 0 then begin
			nratio_rmag_mean[j] = !values.f_nan 
			nratio_rmag_median[j] = !values.f_nan 
		endif else if rmagcount eq 1 then begin
			nratio_rmag_mean[j] = arr[0] 
			nratio_rmag_median[j] = arr[0] 
		endif else begin
			nratio_rmag_mean[j] = mean(arr,/nan)
			nratio_rmag_median[j] = median(arr)
		endelse
	endfor 
	nratio_rmag_linfit = mpfitexpr(linexp,$
		rmag[nrind], $
		nratio_coadd2[nrind] - nratio_normal[nrind], $
		nratio_err[nrind], $
		perr = fel_rmag_linfit_err)


	fel_rmag_mean=fltarr(nm) 
	fel_rmag_median=fltarr(nm) 
	for j=0,nm-1 do begin 
		rmag_ind = where(rmag ge rmagarr[j] and rmag lt rmagarr[j]+delmag,rmagcount)
		arr = (fel_coadd2 - fel_normal)[rmag_ind]
		if rmagcount eq 0 then begin
			fel_rmag_mean[j] = !values.f_nan 
			fel_rmag_median[j] = !values.f_nan 
		endif else if rmagcount eq 1 then begin
			fel_rmag_mean[j] = arr[0] 
			fel_rmag_median[j] = arr[0] 
		endif else begin
			fel_rmag_mean[j] = mean(arr)
			fel_rmag_median[j] = median(arr)
		endelse
	endfor 
	fel_rmag_linfit = mpfitexpr(linexp,$
		rmag, $
		fel_coadd2 - fel_normal, $
		fel_err, $
		perr = fel_rmag_linfit_err)

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd1_rmag.eps',/color, /quiet, xs=8, ys=8, /encap
	cs=1.2
	th=3
	thickline=3
	thinline=1
endif else begin
	cs=2
	th=1
	th=3
	thickline=1
	thinline=1
endelse

!p.multi=[0,2,2]
cgplot, rmag[nrind], alog10(nratio_coadd2[nrind]) - alog10(nratio_normal[nrind]), $
	psym=3, $
	xr=[13,18], /xstyle, $
	yr=[-0.5,0.2], /ystyle, $
	charsize=cs, $
	xtitle='m!Ipetro,R!N [mag]', $
	ytitle='log(n!Iel!N/n!Isp!N)!Icoadd2!N - log(n!Iel!N/n!Isp!N)!Inormal!N'

	cgplot,rmagarr,nratio_alog_rmag_mean,color='red',/overplot, thick=3
	cgplot,rmagarr,nratio_alog_rmag_median,color='green',/overplot, thick=3
	cgplot, /overplot, rmagarr, nratio_alog_rmag_linfit[0] + rmagarr*nratio_alog_rmag_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(rmag[nrind],alog10(nratio_coadd2[nrind])-alog10(nratio_normal[nrind])),format='(f7.4)'),charsize=cs, 13.5, 2.5, /data, color='blue'

cgplot, rmag, fsp_coadd2 - fsp_normal, $
	xr=[13,18], /xstyle, $
	yr=[-0.2,0.2], /ystyle, $
	psym=3, $
	charsize=cs, $
	xtitle='m!Ipetro,R!N [mag]', $
	ytitle='f!Isp,coadd2!N - f!Isp,normal!N'

	cgplot,rmagarr,fsp_rmag_mean,color='red',/overplot, thick=3
	cgplot,rmagarr,fsp_rmag_median,color='green',/overplot, thick=3
	cgplot, /overplot, rmagarr, fsp_rmag_linfit[0] + rmagarr*fsp_rmag_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(rmag, fsp_coadd2 - fsp_normal),format='(f7.4)'),charsize=cs, 13.5, 0.8, /data, color='blue'

cgplot, rmag, nratio_coadd2 - nratio_normal, $
	psym=3, $
	xr=[13,18], /xstyle, $
	yr=[-5,2], /ystyle, $
	charsize=cs, $
	xtitle='m!Ipetro,R!N [mag]', $
	ytitle='(n!Iel!N/n!Isp!N)!Icoadd2!N - (n!Iel!N/n!Isp!N)!Inormal!N'

	cgplot,rmagarr,nratio_rmag_mean,color='red',/overplot, thick=3
	cgplot,rmagarr,nratio_rmag_median,color='green',/overplot, thick=3
	cgplot, /overplot, rmagarr, nratio_rmag_linfit[0] + rmagarr*nratio_rmag_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(rmag[nrind], nratio_coadd2[nrind] - nratio_normal[nrind]),format='(f7.4)'),charsize=cs, 13.5, 30, /data, color='blue'

cgplot, rmag, fel_coadd2 - fel_normal, $
	xr=[13,18], /xstyle, $
	yr=[-0.2,0.2], /ystyle, $
	psym=3, $
	charsize=cs, $
	xtitle='m!Ipetro,R!N [mag]', $
	ytitle='f!Iel,coadd2!N - f!Iel,normal!N'

	cgplot,rmagarr,fel_rmag_mean,color='red',/overplot, thick=3
	cgplot,rmagarr,fel_rmag_median,color='green',/overplot, thick=3
	cgplot, /overplot, rmagarr, fel_rmag_linfit[0] + rmagarr*fel_rmag_linfit[1], $
		linestyle=0, $
		thick=3, $
		color='cyan'

	cgtext, greek('rho')+'='+string(correlate(rmag, fel_coadd2 - fel_normal),format='(f7.4)'),charsize=cs, 13.5, 0.8, /data, color='blue'

if keyword_set(ps) then ps_end

stop

end
