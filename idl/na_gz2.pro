
;+
; NAME:
;       
;	NA_GZ2
;
; PURPOSE:
;
;	Compare classifications from Nair & Abraham vs. Galaxy Zoo 2
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
;       Written by K. Willett                May 12
;-

pro na_gz2, ps=ps, stop=stop, axial=axial, count=count

; Paths for files and plots

csv_path = '~/Astronomy/Research/GalaxyZoo/csv/'
fig_path = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'

; Bars in GZ2 and NA10

file=csv_path+'bars_na_gz2_axial_willettk.csv'
readcol, file, $
	objid, bar_na, bar_gz,$
	t02_edgeon_a05_no_count,t03_bar_a06_bar_count,t03_bar_a07_no_bar_count, $
	expAB_r,expAB_g,$
	format='a,l,f,i,i,i,f,f', $
	skip=1,/silent

	axialcut = where(alog10(1./expab_r) lt 0.3)
	if keyword_set(axial) then begin
		bar_na = bar_na[axialcut]
		bar_gz = bar_gz[axialcut]
		axialtitle='_axial'
	endif else axialtitle=''

	if n_elements(count) gt 0 then begin
		countind = where(t03_bar_a06_bar_count+t03_bar_a07_no_bar_count ge count)
		bar_na = bar_na[countind]
		bar_gz = bar_gz[countind]
		counttitle='_'+strtrim(count,2)
	endif else counttitle=''

	title=axialtitle+counttitle

	print,'Total number of overlap galaxies between samples: ',strtrim(n_elements(bar_na))

;!p.multi=[0,3,4]
;cs=2
;
;for i=0,9 do begin
;	cghistoplot, bar_na[where(bar_gz gt 0.1*i and bar_gz le 0.1*(i+1))],charsize=cs, /freq
;endfor

; Plot the binned GZ bar percentage vs. percentage of Nair galaxies in that bin range that have barflag > 0

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_bars'+axialtitle+counttitle+'.eps',/color, /quiet, xs=14, ys=8, /encap
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

nbins=25
nabinarr = fltarr(nbins)
for i=0,nbins-1 do begin
	gz_percent_ind = where(bar_gz ge i*1./(nbins) and bar_gz lt (i+1)*1./nbins)
	na_bar=where(bar_gz ge i*1./(nbins) and bar_gz lt (i+1)*1./nbins and bar_na gt 0,c_nabar)
	na_nobar=where(bar_gz ge i*1./(nbins) and bar_gz lt (i+1)*1./nbins and bar_na eq 0,c_nanobar)
	nabinarr[i] = float(c_nabar)/(c_nabar + c_nanobar)

endfor

; Plot 1

gzbinarr = fillarr(1./nbins,0,1.-1./nbins)
cgplot, gzbinarr,nabinarr, $
	charsize=cs, $
	xtitle='GZ2 bar fraction', $
	ytitle='Barred NA galaxies/All NA galaxies'

cgplot,/over,gzbinarr,gzbinarr,linestyle=1

corr_na_gz = correlate(gzbinarr,nabinarr)
print,''
print,'Spearman correlation factor: ',corr_na_gz

; Plot 2

cghistoplot,bar_gz[where(bar_na eq 0)], $
	charsize=cs,/outline,/freq,$
	datacolor='blue', $
	xtitle='GZ2 bar fraction'
cghistoplot,bar_gz[where(bar_na ge 1)], /outline, datacolor='red',/oplot,/freq
al_legend,charsize=cs,psym=28,color=['blue','red'],/top,/right,['NA no bar','NA bar']

; Plot 3

cghistoplot,bar_gz[where(bar_na gt 0)], $
	charsize=cs,/outline,$
	thick=th, $
;	/freq,$
	datacolor='Red', $
	xtitle='GZ2 bar fraction'
na_bar1 = where(na_binparse(bar_na,1),c1)
na_bar2 = where(na_binparse(bar_na,2),c2)
na_bar3 = where(na_binparse(bar_na,3),c3)
if c1 gt 4 then cghistoplot,bar_gz[na_bar1], /outline, datacolor='blue',/oplot, thick=th
if c2 gt 4 then cghistoplot,bar_gz[na_bar2], /outline, datacolor='green',/oplot, thick=th
if c3 gt 4 then cghistoplot,bar_gz[na_bar3], /outline, datacolor='purple',/oplot, thick=th

al_legend,charsize=cs,psym=28,color=['red','blue','green','purple'],/top,/left,['NA all bars','Strong bar','Int. bar','Weak bar']

; Plot 4

cghistoplot,bar_gz[where(bar_na gt 0)], $
	charsize=cs,/outline,$
	thick=th, $
;	/freq,$
	datacolor='Red', $
	xtitle='GZ2 bar fraction'

na_bar4 = where(na_binparse(bar_na,4),c4)
na_bar5 = where(na_binparse(bar_na,5),c5)
na_bar6 = where(na_binparse(bar_na,6),c6)
if c4 gt 4 then cghistoplot,bar_gz[na_bar4], /outline, datacolor='hot pink',/oplot, thick=th
if c5 gt 4 then cghistoplot,bar_gz[na_bar5], /outline, datacolor='orange',/oplot, thick=th
if c6 gt 4 then cghistoplot,bar_gz[na_bar6], /outline, datacolor='black',/oplot, thick=th

al_legend,charsize=cs,psym=28,color=['red','hot pink','orange','black'],/top,/left,['NA all bars','Ansae','Peanut','Nuclear bar']

; What type of bar (if any) do the sub-bar features correlate with?

;print, n_elements(setintersection(where(na_binparse(bar_na,4)),where(na_binparse(bar_na,1)))),c4
;print, n_elements(setintersection(where(na_binparse(bar_na,4)),where(na_binparse(bar_na,2)))),c4
;print, n_elements(setintersection(where(na_binparse(bar_na,4)),where(na_binparse(bar_na,3)))),c4
;print, n_elements(setintersection(where(na_binparse(bar_na,5)),where(na_binparse(bar_na,1)))),c5
;print, n_elements(setintersection(where(na_binparse(bar_na,5)),where(na_binparse(bar_na,2)))),c5
;print, n_elements(setintersection(where(na_binparse(bar_na,5)),where(na_binparse(bar_na,3)))),c5
;print, n_elements(setintersection(where(na_binparse(bar_na,6)),where(na_binparse(bar_na,1)))),c6
;print, n_elements(setintersection(where(na_binparse(bar_na,6)),where(na_binparse(bar_na,2)))),c6
;print, n_elements(setintersection(where(na_binparse(bar_na,6)),where(na_binparse(bar_na,3)))),c6

; Load the direct GZ2 data from the downloaded csv table; compare trends for the sample

;gz2file=csv_path+'gz2table_bars.csv'
;objid
;sample
;asset_id
;objid_str
;total_count
;total_weight
;t01_smooth_or_features_a01_smooth_count
;t01_smooth_or_features_a01_smooth_weight
;t01_smooth_or_features_a01_smooth_fraction
;t01_smooth_or_features_a01_smooth_weighted_fraction
;t01_smooth_or_features_a02_features_or_disk_count
;t01_smooth_or_features_a02_features_or_disk_weight
;t01_smooth_or_features_a02_features_or_disk_fraction
;t01_smooth_or_features_a02_features_or_disk_weighted_fraction
;t01_smooth_or_features_a03_star_or_artifact_count
;t01_smooth_or_features_a03_star_or_artifact_weight
;t01_smooth_or_features_a03_star_or_artifact_fraction
;t01_smooth_or_features_a03_star_or_artifact_weighted_fraction
;t01_smooth_or_features_total_count
;t01_smooth_or_features_total_weight
;t02_edgeon_a04_yes_count
;t02_edgeon_a04_yes_weight
;t02_edgeon_a04_yes_fraction
;t02_edgeon_a04_yes_weighted_fraction
;t02_edgeon_a05_no_count
;t02_edgeon_a05_no_weight
;t02_edgeon_a05_no_fraction
;t02_edgeon_a05_no_weighted_fraction
;t02_edgeon_total_count
;t02_edgeon_total_weight
;t03_bar_a06_bar_count
;t03_bar_a06_bar_weight
;t03_bar_a06_bar_fraction
;t03_bar_a06_bar_weighted_fraction
;t03_bar_a07_no_bar_count
;t03_bar_a07_no_bar_weight
;t03_bar_a07_no_bar_fraction
;t03_bar_a07_no_bar_weighted_fraction
;t03_bar_total_count
;t03_bar_total_weight
;readcol, gz2file, csv_objid,sample,asset_id,objid_str,total_count,total_weight,t03_bar_a06_bar_count,t03_bar_a06_bar_weight,t03_bar_a06_bar_fraction,t03_bar_a06_bar_weighted_fraction,t03_bar_a07_no_bar_count,t03_bar_a07_no_bar_weight,t03_bar_a07_no_bar_fraction,t03_bar_a07_no_bar_weighted_fraction,t03_bar_total_count,t03_bar_total_weight,$
;	format='a,a,l,a,i,f,i,f,f,f,i,f,f,f,i,f',/skip,/silent
;
;match, csv_objid, objid, aa, bb
;
;; Simple - what is the bar fraction (re KM11) in the overlap sample? What is it in the total original + extra samples?
;
;indall = where((sample eq 'original' or sample eq 'extra') and t03_bar_total_count gt 10,icount)
;indbar = where((sample eq 'original' or sample eq 'extra') and t03_bar_total_count gt 10 and t03_bar_a06_bar_weighted_fraction ge t03_bar_a07_no_bar_weighted_fraction,barcount)
;indnobar = where((sample eq 'original' or sample eq 'extra') and t03_bar_total_count gt 10 and t03_bar_a06_bar_weighted_fraction lt t03_bar_a07_no_bar_weighted_fraction,nobarcount)
;
;print,'Original + extra (all)',icount,float(barcount)/float(barcount + nobarcount)
;
;; Simple - what is the bar fraction (re KM11) in the overlap sample? What is it in the total original + extra samples?
;
;indall = where((sample[aa] eq 'original' or sample[aa] eq 'extra') and t03_bar_total_count[aa] gt 10,icount)
;indbar = where((sample[aa] eq 'original' or sample[aa] eq 'extra') and t03_bar_total_count[aa] gt 10 and t03_bar_a06_bar_weighted_fraction[aa] ge t03_bar_a07_no_bar_weighted_fraction[aa],barcount)
;indnobar = where((sample[aa] eq 'original' or sample[aa] eq 'extra') and t03_bar_total_count[aa] gt 10 and t03_bar_a06_bar_weighted_fraction[aa] lt t03_bar_a07_no_bar_weighted_fraction[aa],nobarcount)
;
;print,'Original + extra (overlap)',icount,float(barcount)/float(barcount + nobarcount)

;cgtext,/normal,0.52,0.95,alignment=0.5,charsize=cs,title

if keyword_set(ps) then ps_end


; RINGS




ringfile=csv_path+'rings_na_gz2_willettk.csv'
readcol, ringfile, $
	r_objid,ring_na,fring_na,$
	t08_odd_feature_total_count,$
	t08_odd_feature_a23_other_count,t08_odd_feature_a23_other_weighted_fraction,$
	ring_count,$
	ring_wfraction,$
	r_expAB_r,r_expAB_g, $
	format='a,i,i,i,i,f,i,f,f,f', $
	skip=1,/silent

r_axialcut = where(alog10(1./r_expab_r) lt 0.3)
if keyword_set(axial) then begin
	ring_na = ring_na[r_axialcut]
	ring_wfraction = ring_wfraction[r_axialcut]
	ring_count = ring_count[r_axialcut]
	r_axialtitle='_axial'
endif else axialtitle=''

title=axialtitle+counttitle

print,'Total number of overlap galaxies between ring samples: ',strtrim(n_elements(ring_na))

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_rings'+axialtitle+counttitle+'.eps',/color, /quiet, xs=14, ys=8, /encap
	cs=1.2
	th=3
	legsize=3
endif else begin
	cs=2
	th=1
	legsize = 1.2
endelse

!p.multi=[0,2,2]

nbins=25
nabinarr = fltarr(nbins)
for i=0,nbins-1 do begin
	gz_percent_ind = where(ring_wfraction ge i*1./(nbins) and ring_wfraction lt (i+1)*1./nbins)
	na_ring=where(ring_wfraction ge i*1./(nbins) and ring_wfraction lt (i+1)*1./nbins and ring_na gt 0,c_naring)
	na_noring=where(ring_wfraction ge i*1./(nbins) and ring_wfraction lt (i+1)*1./nbins and ring_na eq 0,c_nanoring)
	nabinarr[i] = float(c_naring)/(c_naring + c_nanoring)

endfor

; Plot 1

cghistoplot, ring_count, $
	binsize=2, $
	/log, $
	/outline, $
	datacolor='red', $
	xtitle='Number of votes for ring in GZ2', $
	charsize=cs

cghistoplot, /oplot, ring_count[where(ring_na gt 0)], $
	/outline, $
	datacolor='blue', $
	binsize=2

al_legend,charsize=cs,psym=28,color=['red','blue'],/top,/right,['All galaxies','NA rings']

;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,1))], $
;	/outline, $
;	datacolor='orange', $
;	binsize=2
;
;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,2))], $
;	/outline, $
;	datacolor='green', $
;	binsize=2
;
;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,3))], $
;	/outline, $
;	datacolor='purple', $
;	binsize=2

; Plot 2

rbinsize = 2
rbin = ceil(max(ring_count)/rbinsize)
rvotearr=fltarr(rbin)
rvotearr_inner=fltarr(rbin)
rvotearr_outer=fltarr(rbin)
rvotearr_nuclear=fltarr(rbin)
rbinarr = fillarr(rbinsize,min(ring_count),max(ring_count))
for i=0,rbin-1 do begin
	a = where(ring_count gt i*rbinsize)
	rvotearr[i] = float(n_elements(where(ring_na[a] gt 1)))/float(n_elements(a))
	rvotearr_nuclear[i] = float(n_elements(where(na_binparse(ring_na[a],1))))/float(n_elements(a))
	rvotearr_inner[i] = float(n_elements(where(na_binparse(ring_na[a],2))))/float(n_elements(a))
	rvotearr_outer[i] = float(n_elements(where(na_binparse(ring_na[a],3))))/float(n_elements(a))
endfor
cgplot, rbinarr[0:rbin-1],rvotearr, $
	thick = thickline, $
	xr=[0,60], $
	charsize=cs, $
	yr=[0,1.1],/ystyle, $
	xtitle='Number of votes for ring in GZ2',$ 
	ytitle='Ringed fraction (NA)'

cgplot, rbinarr[0:rbin-1],rvotearr_nuclear, $
	thick=thinline, $
	/over, $
	color='orange'

cgplot, rbinarr[0:rbin-1],rvotearr_inner, $
	thick=thinline, $
	/over, $
	color='green'

cgplot, rbinarr[0:rbin-1],rvotearr_outer, $
	thick=thinline, $
	/over, $
	color='purple'

al_legend,charsize=cs,psym=28,color=['black','orange','green','purple'],/top,/right,['NA all rings','Nuclear ring','Inner ring','Outer ring']
; Plot 3

bin3 = 0.05
cghistoplot, ring_wfraction, $
	/log, $
	binsize=bin3, $
	/outline, $
	datacolor='Red', $
	xtitle='GZ2 vote fraction', $
	charsize=cs

cghistoplot, /oplot, ring_wfraction[where(ring_na gt 0)], $
	/outline, $
	binsize=bin3, $
	datacolor='blue'

al_legend,charsize=cs,psym=28,color=['red','blue'],/top,/right,['All galaxies','NA rings']

; Plot 4

rbinsize = bin3
rbin = ceil(max(ring_wfraction)/rbinsize)
rvotearr=fltarr(rbin)
rvotearr_inner=fltarr(rbin)
rvotearr_outer=fltarr(rbin)
rvotearr_nuclear=fltarr(rbin)
rbinarr = fillarr(rbinsize,min(ring_wfraction),max(ring_wfraction))
for i=0,rbin-1 do begin
	a = where(ring_wfraction gt i*rbinsize)
	rvotearr[i] = float(n_elements(where(ring_na[a] gt 1)))/float(n_elements(a))
	rvotearr_nuclear[i] = float(n_elements(where(na_binparse(ring_na[a],1))))/float(n_elements(a))
	rvotearr_inner[i] = float(n_elements(where(na_binparse(ring_na[a],2))))/float(n_elements(a))
	rvotearr_outer[i] = float(n_elements(where(na_binparse(ring_na[a],3))))/float(n_elements(a))
endfor
cgplot, rbinarr[0:rbin-1],rvotearr, $
	thick=thickline, $
	charsize=cs, $
	yr=[0,1.1],/ystyle, $
	xtitle='GZ2 vote fraction',$ 
	ytitle='Ringed fraction (NA)'

cgplot, rbinarr[0:rbin-1],rvotearr_nuclear, $
	thick=thinline, $
	/over, $
	color='orange'

cgplot, rbinarr[0:rbin-1],rvotearr_inner, $
	thick=thinline, $
	/over, $
	color='green'

cgplot, rbinarr[0:rbin-1],rvotearr_outer, $
	thick=thinline, $
	/over, $
	color='purple'

;binsize=0.1
;cghistoplot,ring_wfraction[where(ring_na gt 0)], $
;	charsize=cs,/outline,$
;	binsize=binsize, $
;	thick=th, $
;	yr=[0,800], $
;;	/freq,$
;	datacolor='Red', $
;	xtitle='GZ2 vote fraction'
;na_ring1 = where(na_binparse(ring_na,1),c1)
;na_ring2 = where(na_binparse(ring_na,2),c2)
;na_ring3 = where(na_binparse(ring_na,3),c3)
;if c1 gt 4 then cghistoplot,ring_wfraction[na_ring1], /outline, datacolor='orange',/oplot, thick=th, binsize=binsize
;if c2 gt 4 then cghistoplot,ring_wfraction[na_ring2], /outline, datacolor='green',/oplot, thick=th, binsize=binsize
;if c3 gt 4 then cghistoplot,ring_wfraction[na_ring3], /outline, datacolor='purple',/oplot, thick=th, binsize=binsize
;
al_legend,charsize=cs,psym=28,color=['black','orange','green','purple'],/top,/right,['NA all rings','Nuclear ring','Inner ring','Outer ring']


; Plot 4

; Load the direct GZ2 data from the downloaded csv table; compare trends for the sample

;match, csv_objid, r_objid, ra, rb
;
;; Simple - what is the ring fraction in the overlap sample? What is it in the total original + extra samples?
;
;indall = where((sample eq 'original' or sample eq 'extra') and ring_count gt 10,icount)
;indbar = where((sample eq 'original' or sample eq 'extra') and ring_count gt 10 and t03_bar_a06_bar_weighted_fraction ge t03_bar_a07_no_bar_weighted_fraction,barcount)
;indnobar = where((sample eq 'original' or sample eq 'extra') and ring_count gt 10 and t03_bar_a06_bar_weighted_fraction lt t03_bar_a07_no_bar_weighted_fraction,nobarcount)
;
;print,'Original + extra (all) rings:',icount,float(barcount)/float(barcount + nobarcount)
;
;indall = where((sample[ra] eq 'original' or sample[ra] eq 'extra') and ring_count[ra] gt 10,icount)
;indbar = where((sample[ra] eq 'original' or sample[ra] eq 'extra') and ring_count[ra] gt 10 and t03_bar_a06_bar_weighted_fraction[ra] ge t03_bar_a07_no_bar_weighted_fraction[ra],barcount)
;indnobar = where((sample[ra] eq 'original' or sample[ra] eq 'extra') and ring_count[ra] gt 10 and t03_bar_a06_bar_weighted_fraction[ra] lt t03_bar_a07_no_bar_weighted_fraction[ra],nobarcount)
;
;print,'Original + extra (overlap) rings:',icount,float(barcount)/float(barcount + nobarcount)



if keyword_set(ps) then ps_end

; PAIRS/INTERACTING



pairfile=csv_path+'pairs_interacting_na_gz2_willettk.csv'
readcol, pairfile, $
	objid,pair_na,fpair_na,interaction_na,tails_na,$
	odd_count,$
	disturbed_count,disturbed_wfraction,$
	other_count,other_wfraction,$
	merger_count,merger_wfraction,$
	pair_expAB_r,pair_expAB_g, $
	format='a,i,i,i,i,i,i,f,i,f,i,f,f,f', $
	skip=1,/silent

p_axialcut = where(alog10(1./pair_expab_r) lt 0.3)
if keyword_set(axial) then begin
	pair_na = pair_na[p_axialcut]
	fpair_na = fpair_na[p_axialcut]
	interaction_na = interaction_na[p_axialcut]
	tails_na = tails_na[p_axialcut]
	p_axialtitle='_axial'
endif else axialtitle=''

title=axialtitle+counttitle

print,'Total number of overlap galaxies between NA/GZ2 samples: ',strtrim(n_elements(pair_na))

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_pairs'+axialtitle+counttitle+'.eps',/color, /quiet, xs=14, ys=8, /encap
	cs=1.2
	th=3
endif else begin
	cs=2
	th=1
endelse

!p.multi=[0,2,2]

; Plot 1

cghistoplot, merger_count, $
	binsize=2, $
	/log, $
	datacolor='black', $
	/outline, $
	xtitle='Number of votes for merger in GZ2', $
	charsize=cs

cghistoplot, /oplot, merger_count[where(pair_na gt 0)], $
	/outline, $
	datacolor='blue', $
	binsize=2

cghistoplot, /oplot, merger_count[where(interaction_na gt 2)], $
	/outline, $
	datacolor='red', $
	binsize=2

cghistoplot, /oplot, merger_count, $
	/outline, $
	datacolor='black', $
	binsize=2

al_legend,charsize=cs,psym=28,color=['black','blue','red'],/top,/right,['All galaxies','NA pairs','NA interacting']

;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,1))], $
;	/outline, $
;	datacolor='orange', $
;	binsize=2
;
;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,2))], $
;	/outline, $
;	datacolor='green', $
;	binsize=2
;
;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,3))], $
;	/outline, $
;	datacolor='purple', $
;	binsize=2

; Plot 2

mbinsize = 2
mbin = ceil(max(merger_count)/mbinsize)
mvotearr=fltarr(mbin)
pvotearr=fltarr(mbin)
mbinarr = fillarr(mbinsize,min(merger_count),max(merger_count))
for i=0,mbin-1 do begin
	a = where(merger_count gt i*mbinsize)
	mvotearr[i] = float(n_elements(where(interaction_na[a] gt 2)))/float(n_elements(a))
	pvotearr[i] = float(n_elements(where(pair_na[a] gt 1)))/float(n_elements(a))
endfor

cgplot, mbinarr[0:mbin-1],mvotearr, $
	charsize=cs, $
	color='red', $
	yr=[0,1.1],/ystyle, $
	xtitle='Number of votes for merger in GZ2',$ 
	ytitle='Interaction/pair fraction (NA)'

cgplot, mbinarr[0:mbin-1],pvotearr, $
	/over, $
	color='blue'

; Plot 3

bin3 = 0.05
cghistoplot, merger_wfraction, $
	/log, $
	binsize=bin3, $
	/outline, $
	datacolor='black', $
	xtitle='GZ2 merger vote fraction', $
	charsize=cs

cghistoplot, /oplot, merger_wfraction[where(pair_na gt 0)], $
	binsize=bin3, $
	/outline, $
	datacolor='blue'

cghistoplot, /oplot, merger_wfraction[where(interaction_na gt 2)], $
	binsize=bin3, $
	/outline, $
	datacolor='red'

al_legend,charsize=cs,psym=28,color=['black','blue','red'],/top,/right,['All galaxies','NA pairs','NA interacting']

;kstwo, merger_wfraction, merger_wfraction[where(pair_na gt 0)], d, p & print,d,p
;kstwo, merger_wfraction, merger_wfraction[where(interaction_na gt 0)], d, p & print,d,p
;kstwo, merger_wfraction[where(pair_na gt 0)], merger_wfraction[where(interaction_na gt 0)], d, p & print,d,p

; Plot 4

mbinsize = bin3
mbin = ceil(max(merger_wfraction)/mbinsize)
mvotearr=fltarr(mbin)
pvotearr=fltarr(mbin)
mbinarr = fillarr(mbinsize,min(merger_wfraction),max(merger_wfraction))
for i=0,mbin-1 do begin
	a = where(merger_wfraction gt i*mbinsize)
	mvotearr[i] = float(n_elements(where(interaction_na[a] gt 2)))/float(n_elements(a))
	pvotearr[i] = float(n_elements(where(pair_na[a] gt 1)))/float(n_elements(a))
endfor

cgplot, mbinarr[0:mbin-1],mvotearr, $
	charsize=cs, $
	color='red', $
	yr=[0,1.1],/ystyle, $
	xtitle='GZ2 merger vote fraction',$ 
	ytitle='Interaction/pair fraction (NA)'

cgplot, mbinarr[0:mbin-1],pvotearr, $
	/over, $
	color='blue'

if keyword_set(ps) then ps_end

; Plot 5

!p.multi=[0,2,2]
binsize=2

cghistoplot,merger_count[where(pair_na gt 0)], $
	charsize=cs,/outline,$
	binsize=binsize, $
	thick=3, $
	/log, $
	datacolor='black', $
	xtitle='GZ2 merger votes'
na_pair1 = where(na_binparse(pair_na,1),c1)
na_pair2 = where(na_binparse(pair_na,2),c2)
na_pair3 = where(na_binparse(pair_na,3),c3)
na_pair4 = where(na_binparse(pair_na,4),c4)
if c1 gt 4 then cghistoplot,merger_count[na_pair1], /outline, datacolor='blue',/oplot, thick=3, binsize=binsize
if c2 gt 4 then cghistoplot,merger_count[na_pair2], /outline, datacolor='green',/oplot, thick=3, binsize=binsize
if c3 gt 4 then cghistoplot,merger_count[na_pair3], /outline, datacolor='purple',/oplot, thick=3, binsize=binsize
if c4 gt 4 then cghistoplot,merger_count[na_pair4], /outline, datacolor='red',/oplot, thick=3, binsize=binsize

al_legend,charsize=cs/2.,psym=28,color=['black','blue','green','purple','red'],/top,/right,['All pairs','Close','Projected','Adjacent','Overlapping']

; Plot 6

cghistoplot,merger_count[where(interaction_na gt 0)], $
	charsize=cs,/outline,$
	binsize=binsize, $
	thick=3, $
	/log, $
	datacolor='black', $
	xtitle='GZ2 merger votes'
na_pair1 = where(na_binparse(interaction_na,1),c1)
na_pair2 = where(na_binparse(interaction_na,2),c2)
na_pair3 = where(na_binparse(interaction_na,3),c3)
na_pair4 = where(na_binparse(interaction_na,4),c4)
na_pair5 = where(na_binparse(interaction_na,5),c5)
na_pair6 = where(na_binparse(interaction_na,6),c6)
na_pair7 = where(na_binparse(interaction_na,7),c7)
na_pair8 = where(na_binparse(interaction_na,8),c8)
if c1 gt 4 then cghistoplot,merger_count[na_pair1], /outline, datacolor='blue',/oplot, thick=3, binsize=binsize
if c2 gt 4 then cghistoplot,merger_count[na_pair2], /outline, datacolor='green',/oplot, thick=3, binsize=binsize
if c3 gt 4 then cghistoplot,merger_count[na_pair3], /outline, datacolor='purple',/oplot, thick=3, binsize=binsize
if c4 gt 4 then cghistoplot,merger_count[na_pair4], /outline, datacolor='red',/oplot, thick=3, binsize=binsize
if c5 gt 4 then cghistoplot,merger_count[na_pair5], /outline, datacolor='yellow',/oplot, thick=3, binsize=binsize
if c6 gt 4 then cghistoplot,merger_count[na_pair6], /outline, datacolor='orange',/oplot, thick=3, binsize=binsize
if c7 gt 4 then cghistoplot,merger_count[na_pair7], /outline, datacolor='hotpink',/oplot, thick=3, binsize=binsize
if c8 gt 4 then cghistoplot,merger_count[na_pair8], /outline, datacolor='black',/oplot, thick=3, binsize=binsize

al_legend,charsize=cs/2.,psym=28,color=['black','blue','green','purple','red','yellow','orange','hotpink','black'],/top,/right,['Interacting','None','Disturbed','Warp','Shells','Short tail','Medium tail','Long tail','Bridge']

; Plot 7

binsize=0.1

cghistoplot,merger_wfraction[where(pair_na gt 0)], $
	charsize=cs,/outline,$
	binsize=binsize, $
	thick=3, $
;	yr=[0,800], $
;	/freq,$
	/log, $
	datacolor='black', $
	xtitle='GZ2 merger vote fraction'
na_pair1 = where(na_binparse(pair_na,1),c1)
na_pair2 = where(na_binparse(pair_na,2),c2)
na_pair3 = where(na_binparse(pair_na,3),c3)
na_pair4 = where(na_binparse(pair_na,4),c4)
if c1 gt 4 then cghistoplot,merger_wfraction[na_pair1], /outline, datacolor='blue',/oplot, thick=3, binsize=binsize
if c2 gt 4 then cghistoplot,merger_wfraction[na_pair2], /outline, datacolor='green',/oplot, thick=3, binsize=binsize
if c3 gt 4 then cghistoplot,merger_wfraction[na_pair3], /outline, datacolor='purple',/oplot, thick=3, binsize=binsize
if c4 gt 4 then cghistoplot,merger_wfraction[na_pair4], /outline, datacolor='red',/oplot, thick=3, binsize=binsize

al_legend,charsize=cs/2.,psym=28,color=['black','blue','green','purple','red'],/top,/right,['All pairs','Close','Projected','Adjacent','Overlapping']

; Plot 8

cghistoplot,merger_wfraction[where(interaction_na gt 0)], $
	charsize=cs,/outline,$
	binsize=binsize, $
	thick=3, $
;	yr=[0,800], $
;	/freq,$
	/log, $
	datacolor='black', $
	xtitle='GZ2 merger vote fraction'
na_pair1 = where(na_binparse(interaction_na,1),c1)
na_pair2 = where(na_binparse(interaction_na,2),c2)
na_pair3 = where(na_binparse(interaction_na,3),c3)
na_pair4 = where(na_binparse(interaction_na,4),c4)
na_pair5 = where(na_binparse(interaction_na,5),c5)
na_pair6 = where(na_binparse(interaction_na,6),c6)
na_pair7 = where(na_binparse(interaction_na,7),c7)
na_pair8 = where(na_binparse(interaction_na,8),c8)
if c1 gt 4 then cghistoplot,merger_wfraction[na_pair1], /outline, datacolor='blue',/oplot, thick=3, binsize=binsize
if c2 gt 4 then cghistoplot,merger_wfraction[na_pair2], /outline, datacolor='green',/oplot, thick=3, binsize=binsize
if c3 gt 4 then cghistoplot,merger_wfraction[na_pair3], /outline, datacolor='purple',/oplot, thick=3, binsize=binsize
if c4 gt 4 then cghistoplot,merger_wfraction[na_pair4], /outline, datacolor='red',/oplot, thick=3, binsize=binsize
if c5 gt 4 then cghistoplot,merger_wfraction[na_pair5], /outline, datacolor='yellow',/oplot, thick=3, binsize=binsize
if c6 gt 4 then cghistoplot,merger_wfraction[na_pair6], /outline, datacolor='orange',/oplot, thick=3, binsize=binsize
if c7 gt 4 then cghistoplot,merger_wfraction[na_pair7], /outline, datacolor='hotpink',/oplot, thick=3, binsize=binsize
if c8 gt 4 then cghistoplot,merger_wfraction[na_pair8], /outline, datacolor='black',/oplot, thick=3, binsize=binsize

al_legend,charsize=cs/2.,psym=28,color=['black','blue','green','purple','red','yellow','orange','hotpink','black'],/top,/right,['Interacting','None','Disturbed','Warp','Shells','Short tail','Medium tail','Long tail','Bridge']

; Check number of tails

!p.multi=[0,2,1]

cghistoplot,merger_count, $
	charsize=cs,/outline,$
	binsize=2, $
	thick=3, $
;	/log, $
	yr=[0,100], $
	datacolor='black', $
	xtitle='GZ2 merger votes'
na_tails1 = where(na_binparse(tails_na,1),c1)
na_tails2 = where(na_binparse(tails_na,2),c2)
na_tails3 = where(na_binparse(tails_na,3),c3)
na_tails4 = where(na_binparse(tails_na,4),c4)
if c1 gt 4 then cghistoplot,merger_count[na_tails1], /outline, datacolor='blue',/oplot, thick=3, binsize=2
if c2 gt 4 then cghistoplot,merger_count[na_tails2], /outline, datacolor='green',/oplot, thick=3, binsize=2
if c3 gt 4 then cghistoplot,merger_count[na_tails3], /outline, datacolor='purple',/oplot, thick=3, binsize=2
if c4 gt 4 then cghistoplot,merger_count[na_tails4], /outline, datacolor='red',/oplot, thick=3, binsize=2

al_legend,charsize=cs/2.,psym=28,color=['black','blue','green','purple','red'],/top,/right,['No tails','1 tail','2 tails','3+ tails','Bunny']

cghistoplot,merger_wfraction, $
	charsize=cs,/outline,$
	binsize=0.1, $
	thick=3, $
	;/log, $
	yr=[0,100], $
	datacolor='black', $
	xtitle='GZ2 merger vote fraction'
if c1 gt 4 then cghistoplot,merger_wfraction[na_tails1], /outline, datacolor='blue',/oplot, thick=3, binsize=0.1
if c2 gt 4 then cghistoplot,merger_wfraction[na_tails2], /outline, datacolor='green',/oplot, thick=3, binsize=0.1
if c3 gt 4 then cghistoplot,merger_wfraction[na_tails3], /outline, datacolor='purple',/oplot, thick=3, binsize=0.1
if c4 gt 4 then cghistoplot,merger_wfraction[na_tails4], /outline, datacolor='red',/oplot, thick=3, binsize=0.1

al_legend,charsize=cs/2.,psym=28,color=['black','blue','green','purple','red'],/top,/right,['No tails','1 tail','2 tails','3+ tails','Bunny']


; T-Types


ttypefile=csv_path+'ttypes_na_gz2_willettk.csv'
readcol, ttypefile, $
	t_objid,ttype,$
	gz_smooth_count,gz_smooth_wfraction,$
	gz_features_or_disk_count,gz_features_or_disk_wfraction,$
	gz_spiral_count,gz_spiral_wfraction,$
	gz_no_spiral_count,gz_no_spiral_wfraction,$
	gz_tight_count,gz_tight_wfraction,$
	gz_medium_count,gz_medium_wfraction,$
	gz_loose_count,gz_loose_wfraction, $
	format='a,i,i,f,i,f,i,f,i,f,i,f,i,f,i,f', $
	skip=1,/silent

match, t_objid, r_objid, tind, rind

t_axialcut = where(alog10(1./r_expab_r[tind]) lt 0.3)
if keyword_set(axial) then begin
	ttype = ttype[t_axialcut]
	gz_smooth_count=gz_smooth_count[t_axialcut]
	gz_smooth_wfraction=gz_smooth_wfraction[t_axialcut]
	gz_features_or_disk_count=gz_features_or_disk_count[t_axialcut]
	gz_features_or_disk_wfraction=gz_features_or_disk_wfraction[t_axialcut]
	gz_spiral_count=gz_spiral_count[t_axialcut]
	gz_spiral_wfraction=gz_spiral_wfraction[t_axialcut]
	gz_no_spiral_count=gz_no_spiral_count[t_axialcut]
	gz_no_spiral_wfraction=gz_no_spiral_wfraction[t_axialcut]
	gz_tight_count=gz_tight_count[t_axialcut]
	gz_tight_wfraction=gz_tight_wfraction[t_axialcut]
	gz_medium_count=gz_medium_count[t_axialcut]
	gz_medium_wfraction=gz_medium_wfraction[t_axialcut]
	gz_loose_count=gz_loose_count[t_axialcut]
	gz_loose_wfraction=gz_loose_wfraction[t_axialcut]
	t_axialtitle='_axial'
endif else axialtitle=''

!p.multi=[0,2,1]

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_ttype'+axialtitle+counttitle+'.eps',/color, /quiet,/encap,xs=15.0,ys=7.0
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

bin = 0.05

cghistoplot, gz_features_or_disk_wfraction[where(ttype eq -5)], $
	xtitle='Features/disk vote fraction (GZ2)', $
	thick=thickline, $
	xr=[-0.1,1.1], $
	yr=[0,2e3], $
;	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='red'
cghistoplot, gz_features_or_disk_wfraction[where(ttype ge -3 and ttype le 0)], $
	thick=thickline, $
;	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='lime green'

cghistoplot, gz_features_or_disk_wfraction[where(ttype gt 1 and ttype le 8)], $
	thick=thinline, $
;	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['E','S0','Sa-Sd'],color=['red','lime green','blue'], psym=28, /top, /left, charsize=cs*1.5

spiralcount = 10
totalsc = gz_spiral_count + gz_no_spiral_count

cghistoplot, gz_spiral_wfraction[where(ttype eq -5 and totalsc ge spiralcount)], $
	xtitle='Spiral vote fraction (GZ2)', $
	thick=thickline, $
	yr=[0,2e3], $
	xr=[-0.1,1.1], $
;	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='red'
cghistoplot, gz_spiral_wfraction[where(ttype ge -3 and ttype le 0 and totalsc ge spiralcount)], $
	thick=thickline, $
;	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='lime green'

cghistoplot, gz_spiral_wfraction[where(ttype gt 0 and ttype le 8 and totalsc ge spiralcount)], $
	thick=thinline, $
;	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['E','S0','Sa-Sd'],color=['red','lime green','blue'], psym=28, /top, /left, charsize=cs*1.5

if keyword_set(ps) then ps_end

	; Check brightness of the NA late-type galaxies with low GZ2 spiral WFs

	readcol, csv_path+'mags_na_gz2_willettk.csv', $
		mag_objid, junk, gmag, rmag, gmag_abs, $
		format='a,a,f,f,f', skip=1, /silent

	match, mag_objid, t_objid, mm, tt

	gmag = gmag[mm]
	rmag = rmag[mm]
	gmag_abs = gmag_abs[mm]

	spiral_cutoff = 0.8

	;for i=1,9 do begin

	;	spiral_cutoff = 0.1*i

	;	print,'Spiral cutoff: ',spiral_cutoff

	;	print,'N gals'
	;	print, n_elements(gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)])
	;	print, n_elements(gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)])
	;	print,'g'
	;	print, mean(gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)])
	;	print, mean(gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)])
	;	kstwo, gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)],gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)],d,p & print,d,p
	;	print,'r'
	;	print, mean(rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)])
	;	print, mean(rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)])
	;	kstwo, rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)],rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)],d,p & print,d,p
	;	print,'g abs'
	;	print, mean(gmag_abs[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)])
	;	print, mean(gmag_abs[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)])
	;	kstwo, gmag_abs[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)],gmag_abs[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)],d,p & print,d,p
	;	print,'g-r'
	;	print, mean(gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)] - $
	;		rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)])
	;	print, mean(gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)] - $
	;		rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)])
	;	kstwo, gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)] - rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction ge spiral_cutoff)],gmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)] - rmag[where(ttype ge 1 and ttype le 8 and totalsc gt spiralcount and gz_spiral_wfraction lt spiral_cutoff)], d, p & print,d,p

	;endfor


; Tightness of spirals

!p.multi=[0,2,3]

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_spiraltightness'+axialtitle+counttitle+'.ps',/color, /quiet
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

bin=0.05
cghistoplot, gz_tight_wfraction[where(totalsc ge spiralcount)], $
	thick=3, $
	xtitle='Tight spiral weighted fraction (GZ2)', $
	xr=[-0.1,1.1], $
	yr=[0,0.5],$ 
	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='black'

cghistoplot, gz_tight_wfraction[where(ttype ge -3 and ttype le 0 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='red'

cghistoplot, gz_tight_wfraction[where(ttype ge 1 and ttype le 2 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='orange'

cghistoplot, gz_tight_wfraction[where(ttype ge 3 and ttype le 4 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='yellow'

cghistoplot, gz_tight_wfraction[where(ttype ge 5 and ttype le 6 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='green'

cghistoplot, gz_tight_wfraction[where(ttype ge 7 and ttype le 8 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['Spirals','S0','Sa','Sb','Sc','Sd'], color=['black','red','orange','yellow','green','blue'], psym=28, /top, /right, charsize=cs/2.

bin=2
cghistoplot, gz_tight_count[where(totalsc ge spiralcount)], $
	thick=3, $
	xtitle='Tight spiral number of votes (GZ2)', $
	yr=[0,0.5],$ 
	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='black'

cghistoplot, gz_tight_count[where(ttype ge -3 and ttype le 0 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='red'

cghistoplot, gz_tight_count[where(ttype ge 1 and ttype le 2 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='orange'

cghistoplot, gz_tight_count[where(ttype ge 3 and ttype le 4 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='yellow'

cghistoplot, gz_tight_count[where(ttype ge 5 and ttype le 6 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='green'

cghistoplot, gz_tight_count[where(ttype ge 7 and ttype le 8 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['Spirals','S0','Sa','Sb','Sc','Sd'], color=['black','red','orange','yellow','green','blue'], psym=28, /top, /right, charsize=cs/2.

bin=0.05
cghistoplot, gz_medium_wfraction[where(totalsc ge spiralcount)], $
	thick=3, $
	xtitle='Medium spiral weighted fraction (GZ2)', $
	xr=[-0.1,1.1], $
	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='black'

cghistoplot, gz_medium_wfraction[where(ttype ge -3 and ttype le 0 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='red'

cghistoplot, gz_medium_wfraction[where(ttype ge 1 and ttype le 2 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='orange'

cghistoplot, gz_medium_wfraction[where(ttype ge 3 and ttype le 4 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='yellow'

cghistoplot, gz_medium_wfraction[where(ttype ge 5 and ttype le 6 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='green'

cghistoplot, gz_medium_wfraction[where(ttype ge 7 and ttype le 8 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['Spirals','S0','Sa','Sb','Sc','Sd'], color=['black','red','orange','yellow','green','blue'], psym=28, /top, /right, charsize=cs/2.

bin=2
cghistoplot, gz_medium_count[where(totalsc ge spiralcount)], $
	thick=3, $
	xtitle='Medium spiral number of votes (GZ2)', $
	yr=[0,0.5],$ 
	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='black'

cghistoplot, gz_medium_count[where(ttype ge -3 and ttype le 0 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='red'

cghistoplot, gz_medium_count[where(ttype ge 1 and ttype le 2 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='orange'

cghistoplot, gz_medium_count[where(ttype ge 3 and ttype le 4 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='yellow'

cghistoplot, gz_medium_count[where(ttype ge 5 and ttype le 6 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='green'

cghistoplot, gz_medium_count[where(ttype ge 7 and ttype le 8 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['Spirals','S0','Sa','Sb','Sc','Sd'], color=['black','red','orange','yellow','green','blue'], psym=28, /top, /right, charsize=cs/2.

bin=0.05
cghistoplot, gz_loose_wfraction[where(totalsc ge spiralcount)], $
	thick=3, $
	xtitle='Loose spiral weighted fraction (GZ2)', $
	yr=[0,0.5],$ 
	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='black'

cghistoplot, gz_loose_wfraction[where(ttype ge -3 and ttype le 0 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='red'

cghistoplot, gz_loose_wfraction[where(ttype ge 1 and ttype le 2 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='orange'

cghistoplot, gz_loose_wfraction[where(ttype ge 3 and ttype le 4 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='yellow'

cghistoplot, gz_loose_wfraction[where(ttype ge 5 and ttype le 6 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='green'

cghistoplot, gz_loose_wfraction[where(ttype ge 7 and ttype le 8 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['Spirals','S0','Sa','Sb','Sc','Sd'], color=['black','red','orange','yellow','green','blue'], psym=28, /top, /right, charsize=cs/2.

bin=2
cghistoplot, gz_loose_count[where(totalsc ge spiralcount)], $
	thick=3, $
	xtitle='Loose spiral number of votes (GZ2)', $
	yr=[0,0.5],$ 
	/freq, $
	/outline, $
	binsize=bin, $
	datacolor='black'

cghistoplot, gz_loose_count[where(ttype ge -3 and ttype le 0 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='red'

cghistoplot, gz_loose_count[where(ttype ge 1 and ttype le 2 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='orange'

cghistoplot, gz_loose_count[where(ttype ge 3 and ttype le 4 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='yellow'

cghistoplot, gz_loose_count[where(ttype ge 5 and ttype le 6 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='green'

cghistoplot, gz_loose_count[where(ttype ge 7 and ttype le 8 and totalsc ge spiralcount)], $
	thick=1, $
	/freq, $
	/outline, $
	binsize=bin, $
	/oplot, $
	datacolor='blue'

al_legend, ['Spirals','S0','Sa','Sb','Sc','Sd'], color=['black','red','orange','yellow','green','blue'], psym=28, /top, /right, charsize=cs/2.

if keyword_set(ps) then ps_end

; Colorful histogram for spiral tightness

!p.multi=[0,2,3]

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_spiraltightness_color'+axialtitle+counttitle+'.ps',/color, /quiet
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

sd = where(ttype ge 7  and ttype le 10)
sc = where(ttype ge 5  and ttype le 6)
sb = where(ttype ge 3  and ttype le 4)
sa = where(ttype ge 1  and ttype le 2)
s0 = where(ttype ge -3 and ttype le 0)
es = where(ttype eq -5)

bs = 0.1

; Tight spirals

tight_countlim = gz_tight_wfraction[where(totalsc ge spiralcount)]
ttype_countlim = ttype[where(totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='Tight spiral weighted fraction (GZ2)', $
	ytitle='Number of galaxies', $
	xr=[-0.1,1.1], $
	yr=[-200,1000],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	;print,csd,csc,csb,csa,cs0,ces

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='Tight spiral weighted fraction (GZ2)', $
	ytitle='Fraction per bin', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(tight_countlim gt x0 and tight_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	junk =  where(tight_countlim gt x0 and tight_countlim le x1,c_all)

	;print,csd,csc,csb,csa,cs0,ces

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

; Medium spirals

medium_countlim = gz_medium_wfraction[where(totalsc ge spiralcount)]
ttype_countlim = ttype[where(totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='Medium spiral weighted fraction (GZ2)', $
	ytitle='Number of galaxies', $
	xr=[-0.1,1.1], $
	yr=[-200,1500],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	;print,csd,csc,csb,csa,cs0,ces

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='Medium spiral weighted fraction (GZ2)', $
	ytitle='Fraction per bin', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(medium_countlim gt x0 and medium_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	junk =  where(medium_countlim gt x0 and medium_countlim le x1,c_all)

	;print,csd,csc,csb,csa,cs0,ces

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor

; Loose spirals

loose_countlim = gz_loose_wfraction[where(totalsc ge spiralcount)]
ttype_countlim = ttype[where(totalsc ge spiralcount)]

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='Loose spiral weighted fraction (GZ2)', $
	ytitle='Number of galaxies', $
	xr=[-0.1,1.1], $
	yr=[-200,4300],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	;print,csd,csc,csb,csa,cs0,ces

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif


endfor

al_legend, charsize=cs/1.5, ['E','S0','Sa','Sb','Sc','Sd'],color=['black','tomato','yellow','green','dark green','blue'], psym=28, /top, /right

cgplot, indgen(10), $
	/nodata, $
;	xtickinterval=0.01, $
	charsize=cs, $
	xtitle='Loose spiral weighted fraction (GZ2)', $
	ytitle='Fraction per bin', $
	xr=[-0.1,1.1], /xstyle, $
	yr=[-0.2,1.2],/ystyle

for i=0,9 do begin
	x0 = i*bs
	x1 = (i+1)*bs
	if i eq 0 then x0 = -0.01
	y0 = 0.
	y1 = 0.
	z_sd = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 7,csd)
	z_sc = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 5 and ttype_countlim le 6,csc)
	z_sb = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 3 and ttype_countlim le 4,csb)
	z_sa = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge 1 and ttype_countlim le 2,csa)
	z_s0 = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge -3 and ttype_countlim le 0,cs0)
	z_es = where(loose_countlim gt x0 and loose_countlim le x1 and ttype_countlim ge -5 and ttype_countlim le -4,ces)

	junk =  where(loose_countlim gt x0 and loose_countlim le x1,c_all)

	;print,csd,csc,csb,csa,cs0,ces

	if ces gt 0 then begin
		y0 = y1
		y1 = y0 + ces/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='black'
	endif

	if cs0 gt 0 then begin
		y0 = y1
		y1 = y0 + cs0/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='tomato'
	endif

	if csa gt 0 then begin
		y0 = y1
		y1 = y0 + csa/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='yellow'
	endif

	if csb gt 0 then begin
		y0 = y1
		y1 = y0 + csb/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='green'
	endif

	if csc gt 0 then begin
		y0 = y1
		y1 = y0 + csc/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='dark green'
	endif

	if csd gt 0 then begin
		y0 = y1
		y1 = y0 + csd/float(c_all)
		cgcolorfill, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], color='blue'
	endif

	cgtext, (x0+x1)/2., 1.05, /data, strtrim(c_all,2), alignment=0.5, charsize=cs/2.


endfor


if keyword_set(ps) then ps_end

if keyword_set(stop) then stop

print,''
end
