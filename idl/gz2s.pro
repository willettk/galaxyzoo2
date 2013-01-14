
;+
; NAME:
;       
;	GZ2_STRIPE82_COADD2
;
; PURPOSE:
;
;	Plot the morphological attributes for single- and coadd2-depth images from Stripe 82 in SDSS
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

pro gz2s, ps=ps

; Import Stripe 82 data

gz2dir = '~/Astronomy/Research/GalaxyZoo/fits/'
gz2 = mrdfits(gz2dir+'gz2_stripe82_normal_coadd2.fits',1)

tagnames=tag_names(gz2)
nt = n_elements(tagnames)

temparr_normal=intarr(nt) & for i=0,nt-1 do begin & temp = tagnames[i] & templen = strlen(temp) & if templen ge 20 then begin & temp2 = strmid(temp,templen-19,19) & if temp2 eq 'WEIGHTED_FRACTION_1' then temparr_normal[i]=1 & endif & endfor
temparr_coadd2=intarr(nt) & for i=0,nt-1 do begin & temp = tagnames[i] & templen = strlen(temp) & if templen ge 20 then begin & temp2 = strmid(temp,templen-19,19) & if temp2 eq 'WEIGHTED_FRACTION_2' then temparr_coadd2[i]=1 & endif & endfor

nbyn = 6

; Scatter plot: coadd vs. normal

!p.multi=[0,nbyn,nbyn]

count=5 

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd2_normal_scatter.ps',/color, /quiet
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

	tagind_normal_wf = (where(temparr_normal))[i]
	tagind_normal_ct = (where(temparr_normal))[i] - 2
	tagind_coadd2_wf = (where(temparr_coadd2))[i]
	tagind_coadd2_ct = (where(temparr_coadd2))[i] - 2
	countind = where(gz2.(tagind_normal_ct) ge count and gz2.(tagind_coadd2_ct) ge count)

	coadd2_val = gz2[countind].(tagind_coadd2_wf)
	normal_val = gz2[countind].(tagind_normal_wf)
	z = gz2[countind].redshift_1

	cgplot, coadd2_val, normal_val, $
		psym=3, $
		xr=[0,1], /xstyle, $
		yr=[0,1], /ystyle, $
		title='Task '+strmid(tagnames[tagind_normal_wf],1,2)+', Answer '+string(i+1,format='(i02)'), $
		xtitle='coadd2', ytitle='normal'
	cgplot, indgen(11)/10., indgen(11)/10., color='red' ,/overplot
endfor

if keyword_set(ps) then ps_end

; Histogram: coadd - normal

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd2_normal_histogram.ps',/color, /quiet
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

	tagind_normal_wf = (where(temparr_normal))[i]
	tagind_normal_ct = (where(temparr_normal))[i] - 2
	tagind_coadd2_wf = (where(temparr_coadd2))[i]
	tagind_coadd2_ct = (where(temparr_coadd2))[i] - 2
	countind = where(gz2.(tagind_normal_ct) ge count and gz2.(tagind_coadd2_ct) ge count)

	coadd2_val = gz2[countind].(tagind_coadd2_wf)
	normal_val = gz2[countind].(tagind_normal_wf)
	z = gz2[countind].redshift_1

	cghistoplot, coadd2_val - normal_val, $
		title='Task '+strmid(tagnames[tagind_normal_wf],1,2)+$
		', Answer '+string(i+1,format='(i02)'), $
;		+' '+greek('mu')+'='+string(mean(coadd2_val - normal_val),format='(f5.2)')+', '+$
;		greek('sigma')+'='+string(stddev(coadd2_val - normal_val),format='(f5.2)'),$
		datacolor='blue',$
		/freq,$
		xr=[-1,1],/xstyle,$
		xtitle='coadd2 - normal'
	cgplots, [0,0], !y.crange, color='red' 

endfor

if keyword_set(ps) then ps_end

; Scatter plot: redshift vs. (coadd - normal)

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_redshift_deltaf_coadd2.ps',/color, /quiet
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

	tagind_normal_wf = (where(temparr_normal))[i]
	tagind_normal_ct = (where(temparr_normal))[i] - 2
	tagind_coadd2_wf = (where(temparr_coadd2))[i]
	tagind_coadd2_ct = (where(temparr_coadd2))[i] - 2
	countind = where(gz2.(tagind_normal_ct) ge count and gz2.(tagind_coadd2_ct) ge count)

	coadd2_val = gz2[countind].(tagind_coadd2_wf)
	normal_val = gz2[countind].(tagind_normal_wf)
	z = gz2[countind].redshift_1

	cgplot, z, coadd2_val - normal_val, $
		title='Task '+strmid(tagnames[tagind_normal_wf],1,2)+$
		', Answer '+string(i+1,format='(i02)'), $
		;+' '+greek('mu')+'='+string(mean(coadd2_val - normal_val),format='(f5.2)')+', '+$
		;greek('sigma')+'='+string(stddev(coadd2_val - normal_val),format='(f5.2)'),$
		xtitle='z',$
		ytitle='coadd2 - normal',$
		psym=3,$
		xr=[-0.05,0.30], $
		yr=[-1.1,1.1], $
		/xstyle, /ystyle

	medarr=fltarr(26) 
	for j=0,25 do begin 
		zlower = j/100.
		zupper = (j+1)/100.
		arr = (coadd2_val - normal_val)[where(z gt zlower and z le zupper,zcount)]

		if zcount eq 0 then begin
			if arr[0] eq -1 then medarr[j] = !values.f_nan 
		endif else if zcount eq 1 then begin
			medarr[j] = arr[0] 
		endif else medarr[j] = mean(arr)

	endfor 
	cgplot,indgen(26)/100.,medarr,color='red',/overplot
endfor

if keyword_set(ps) then ps_end

; Plot difference in coadd2 - normal depth data as function of r-band surface brightness

llim = 10

nel_coadd2 = gz2.((where(temparr_coadd2))[0]-3)
nsp_coadd2 = gz2.((where(temparr_coadd2))[1]-3)
nratio_coadd2 = float(nel_coadd2) / float(nsp_coadd2)

nel_normal = gz2.((where(temparr_normal))[0]-3)
nsp_normal = gz2.((where(temparr_normal))[1]-3)
nratio_normal = float(nel_normal) / float(nsp_normal)

fel_coadd2 = gz2.((where(temparr_coadd2))[0])
fel_normal = gz2.((where(temparr_normal))[0])
fsp_coadd2 = gz2.((where(temparr_coadd2))[1])
fsp_normal = gz2.((where(temparr_normal))[1])

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

mu50_r = gz2.mu50_r_1

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
		/quiet, $
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
		/quiet, $
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
		/quiet, $
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
		/quiet, $
		mu50_r, $
		fel_coadd2 - fel_normal, $
		fel_err, $
		perr = fel_mu50_r_linfit_err)

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd2_mu50.eps',/color, /quiet, xs=8, ys=8, /encap
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

rmag = gz2.petromag_r_1

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
		/quiet, $
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
		/quiet, $
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
		/quiet, $
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
		/quiet, $
		rmag, $
		fel_coadd2 - fel_normal, $
		fel_err, $
		perr = fel_rmag_linfit_err)

if keyword_set(ps) then begin
	ps_start, filename='~/Astronomy/Research/GalaxyZoo/s82_coadd2_rmag.eps',/color, /quiet, xs=8, ys=8, /encap
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
	;yr=[-0.2,0.2], /ystyle, $
	yr=[-0.5,0.5], /ystyle, $
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

