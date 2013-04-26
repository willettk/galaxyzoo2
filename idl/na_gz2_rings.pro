;+
; NAME:
;       
;	NA_GZ2_RINGS
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

pro na_gz2_rings, ps=ps, stop=stop, axial=axial, count=count

; Paths for files and plots

csv_path = '~/Astronomy/Research/GalaxyZoo/csv/'
fig_path = '~/Astronomy/Research/GalaxyZoo/datapaper/figures/'
fits_path = '~/Astronomy/Research/GalaxyZoo/fits/'

nagz2 = mrdfits(fits_path+'na_gz2_willettk.fit',1,/silent)

; Limits on previous questions

NED = where(nagz2.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge 0.227 and nagz2.t02_edgeon_a05_no_weighted_fraction ge 0.519)
NED_odd = where(nagz2.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge 0.227 and nagz2.t02_edgeon_a05_no_weighted_fraction ge 0.519 and nagz2.t06_odd_a14_yes_weighted_fraction ge 0.223)
NED_notodd = where(nagz2.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge 0.227 and nagz2.t02_edgeon_a05_no_weighted_fraction ge 0.519 and nagz2.t06_odd_a14_yes_weighted_fraction lt 0.223)

; How does the number of detected ringed galaxies depend on cutting on Task 06?

print,'no odd cut'
print, 'nuclear', n_elements(where(na_binparse(nagz2[NED].ring,1))),n_elements(where(na_binparse(nagz2[NED].ring,1)))/float(n_elements(where(nagz2[NED].ring gt 0)))
print, 'inner',   n_elements(where(na_binparse(nagz2[NED].ring,2))),n_elements(where(na_binparse(nagz2[NED].ring,2)))/float(n_elements(where(nagz2[NED].ring gt 0)))
print, 'outer',   n_elements(where(na_binparse(nagz2[NED].ring,3))),n_elements(where(na_binparse(nagz2[NED].ring,3)))/float(n_elements(where(nagz2[NED].ring gt 0)))
print,''
print,'odd'
print, 'nuclear', n_elements(where(na_binparse(nagz2[NED_odd].ring,1))),n_elements(where(na_binparse(nagz2[NED_odd].ring,1)))/float(n_elements(where(nagz2[NED_odd].ring gt 0)))
print, 'inner',   n_elements(where(na_binparse(nagz2[NED_odd].ring,2))),n_elements(where(na_binparse(nagz2[NED_odd].ring,2)))/float(n_elements(where(nagz2[NED_odd].ring gt 0)))
print, 'outer',   n_elements(where(na_binparse(nagz2[NED_odd].ring,3))),n_elements(where(na_binparse(nagz2[NED_odd].ring,3)))/float(n_elements(where(nagz2[NED_odd].ring gt 0)))
print,''
print,'not odd'
print, 'nuclear', n_elements(where(na_binparse(nagz2[NED_notodd].ring,1))),n_elements(where(na_binparse(nagz2[NED_notodd].ring,1)))/float(n_elements(where(nagz2[NED_notodd].ring gt 0)))
print, 'inner',   n_elements(where(na_binparse(nagz2[NED_notodd].ring,2))),n_elements(where(na_binparse(nagz2[NED_notodd].ring,2)))/float(n_elements(where(nagz2[NED_notodd].ring gt 0)))
print, 'outer',   n_elements(where(na_binparse(nagz2[NED_notodd].ring,3))),n_elements(where(na_binparse(nagz2[NED_notodd].ring,3)))/float(n_elements(where(nagz2[NED_notodd].ring gt 0)))
print,''

nagz2_NED = nagz2[NED]
nagz2_odd = nagz2[NED_odd]

ring_na_NED = nagz2_NED.Ring
ring_count_NED = nagz2_NED.t08_odd_feature_a19_ring_weight
ring_wfraction_NED = nagz2_NED.t08_odd_feature_a19_ring_weighted_fraction

ring_na_odd = nagz2_odd.Ring
ring_count_odd = nagz2_odd.t08_odd_feature_a19_ring_weight
ring_wfraction_odd = nagz2_odd.t08_odd_feature_a19_ring_weighted_fraction

;; RINGS
;
;ringfile=csv_path+'rings_na_gz2_willettk.csv'
;readcol, ringfile, $
;	r_objid,ring_na,fring_na,$
;	t08_odd_feature_total_count,$
;	t08_odd_feature_a23_other_count,t08_odd_feature_a23_other_weighted_fraction,$
;	ring_count,$
;	ring_wfraction,$
;	r_expAB_r,r_expAB_g, $
;	format='a,i,i,i,i,f,i,f,f,f', $
;	skip=1,/silent
;
;; This kluge is so that histogram bins end at vote fractions of 1.0 - KW, 27 Mar 2013
;;gzbar1 = where(bar_gz eq 1.0)
;;bar_gz[gzbar1] = 1d - 1d-3
;
;axialtitle=''
;counttitle=''
;if n_elements(count) eq 0 then count = 0
;
;if keyword_set(axial) then begin
;    cut = where(alog10(1./r_expab_r) lt 0.3 and t08_odd_feature_total_count ge count)
;	ring_na = ring_na[cut]
;	ring_wfraction = ring_wfraction[cut]
;	ring_count = ring_count[cut]
;	axialtitle='_axial'
;endif else begin
;    cut = where(t08_odd_feature_total_count ge count)
;    ring_na = ring_na[cut]
;    ring_wfraction = ring_wfraction[cut]
;    ring_count = ring_count[cut]
;    counttitle='_'+strtrim(count,2)
;endelse
;
;title=axialtitle+counttitle

;print,'Total number of overlap galaxies between ring samples: ',strtrim(n_elements(ring_na))

;if keyword_set(ps) then begin
;	ps_start, filename=fig_path+'na_rings'+axialtitle+counttitle+'.eps',/color, /quiet, xs=14, ys=8, /encap
;	cs=1.2
;	th=3
;	legsize=3
;	thickline = 7
;endif else begin
;	cs=2
;	th=1
;	legsize = 1.2
;	thickline = 1
;endelse
;
;!p.multi=[0,2,2]
;
;nbins=25
;nabinarr = fltarr(nbins)
;for i=0,nbins-1 do begin
;	gz_percent_ind = where(ring_wfraction ge i*1./(nbins) and ring_wfraction lt (i+1)*1./nbins)
;	na_ring=where(ring_wfraction ge i*1./(nbins) and ring_wfraction lt (i+1)*1./nbins and ring_na gt 0,c_naring)
;	na_noring=where(ring_wfraction ge i*1./(nbins) and ring_wfraction lt (i+1)*1./nbins and ring_na eq 0,c_nanoring)
;	nabinarr[i] = float(c_naring)/(c_naring + c_nanoring)
;
;endfor
;
;; Plot 1
;
;cghistoplot, ring_count, $
;	binsize=2, $
;	thick=thickline, $
;	/log, $
;	/outline, $
;	datacolor='red', $
;	xtitle='Number of votes for "ring" (GZ2)', $
;	ytitle='Count', $
;	charsize=cs
;
;cghistoplot, /oplot, ring_count[where(ring_na gt 0)], $
;	/outline, $
;	thick=thickline, $
;	datacolor='blue', $
;	binsize=2
;
;al_legend,charsize=cs,psym=28,color=['red','blue'],/top,/right,['All galaxies','Rings (NA10)']
;
;;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,1))], $
;;	/outline, $
;;	datacolor='orange', $
;;	binsize=2
;;
;;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,2))], $
;;	/outline, $
;;	datacolor='green', $
;;	binsize=2
;;
;;cghistoplot, /oplot, ring_count[where(na_binparse(ring_na,3))], $
;;	/outline, $
;;	datacolor='purple', $
;;	binsize=2
;
;; Plot 2
;
;rbinsize = 2
;rbin = ceil(max(ring_count)/rbinsize)
;rvotearr=fltarr(rbin)
;rvotearr_inner=fltarr(rbin)
;rvotearr_outer=fltarr(rbin)
;rvotearr_nuclear=fltarr(rbin)
;rbinarr = fillarr(rbinsize,min(ring_count),max(ring_count))
;for i=0,rbin-1 do begin
;	a = where(ring_count gt i*rbinsize)
;	rvotearr[i] = float(n_elements(where(ring_na[a] gt 1)))/float(n_elements(a))
;	rvotearr_nuclear[i] = float(n_elements(where(na_binparse(ring_na[a],1))))/float(n_elements(a))
;	rvotearr_inner[i] = float(n_elements(where(na_binparse(ring_na[a],2))))/float(n_elements(a))
;	rvotearr_outer[i] = float(n_elements(where(na_binparse(ring_na[a],3))))/float(n_elements(a))
;endfor
;cgplot, rbinarr[0:rbin-1],rvotearr, $
;	thick = thickline, $
;	xthick=thickline, $
;	ythick=thickline, $
;	xr=[0,60], $
;	charsize=cs, $
;	yr=[0,1.1],/ystyle, $
;	xtitle='Number of votes for "ring" (GZ2)', $
;	ytitle='Ringed fraction (NA10)'
;
;cgplot, rbinarr[0:rbin-1],rvotearr_nuclear, $
;	thick=thickline, $
;	/over, $
;	color='orange'
;
;cgplot, rbinarr[0:rbin-1],rvotearr_inner, $
;	thick=thickline, $
;	/over, $
;	color='green'
;
;cgplot, rbinarr[0:rbin-1],rvotearr_outer, $
;	thick=thickline, $
;	/over, $
;	color='purple'
;
;al_legend,charsize=cs,psym=28,color=['black','orange','green','purple'],/top,/right,['All rings','Nuclear ring','Inner ring','Outer ring']
;; Plot 3
;
;bin3 = 0.05
;cghistoplot, ring_wfraction, $
;	thick = thickline, $
;	/log, $
;	binsize=bin3, $
;	/outline, $
;	datacolor='Red', $
;	xtitle='Vote fraction for "ring" (GZ2)', $
;	ytitle='Count', $
;	charsize=cs
;
;cghistoplot, /oplot, ring_wfraction[where(ring_na gt 0)], $
;	thick = thickline, $
;	/outline, $
;	binsize=bin3, $
;	datacolor='blue'
;
;al_legend,charsize=cs,psym=28,color=['red','blue'],/top,/right,['All galaxies','Rings (NA10)']
;
;; Plot 4
;
;rbinsize = bin3
;rbin = ceil(max(ring_wfraction)/rbinsize)
;rvotearr=fltarr(rbin)
;rvotearr_inner=fltarr(rbin)
;rvotearr_outer=fltarr(rbin)
;rvotearr_nuclear=fltarr(rbin)
;rbinarr = fillarr(rbinsize,min(ring_wfraction),max(ring_wfraction))
;for i=0,rbin-1 do begin
;	a = where(ring_wfraction gt i*rbinsize)
;	rvotearr[i] = float(n_elements(where(ring_na[a] gt 1)))/float(n_elements(a))
;	rvotearr_nuclear[i] = float(n_elements(where(na_binparse(ring_na[a],1))))/float(n_elements(a))
;	rvotearr_inner[i] = float(n_elements(where(na_binparse(ring_na[a],2))))/float(n_elements(a))
;	rvotearr_outer[i] = float(n_elements(where(na_binparse(ring_na[a],3))))/float(n_elements(a))
;endfor
;cgplot, rbinarr[0:rbin-1],rvotearr, $
;	thick=thickline, $
;	xthick=thickline, $
;	ythick=thickline, $
;	charsize=cs, $
;	yr=[0,1.1],/ystyle, $
;	xtitle='Vote fraction for "ring" (GZ2)', $
;	ytitle='Ringed fraction (NA10)'
;
;cgplot, rbinarr[0:rbin-1],rvotearr_nuclear, $
;	thick=thickline, $
;	/over, $
;	color='orange'
;
;cgplot, rbinarr[0:rbin-1],rvotearr_inner, $
;	thick=thickline, $
;	/over, $
;	color='green'
;
;cgplot, rbinarr[0:rbin-1],rvotearr_outer, $
;	thick=thickline, $
;	/over, $
;	color='purple'
;
;al_legend,charsize=cs,psym=28,color=['black','orange','green','purple'],/top,/right,['All rings','Nuclear ring','Inner ring','Outer ring']
;
;if keyword_set(ps) then ps_end


; Second ring plot

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_gz2_rings.eps',/color, /quiet, xs=5, ys=8, /encap
	cs=1.2
	th=7
	thickline=3
	thinline=1
endif else begin
	cs=1.5
	th=1
	thickline=1
	thinline=1
endelse

!p.multi=[0,1,2]

; Rings vs. no rings

bs = 0.05
cghistoplot,ring_wfraction_NED[where(ring_na_NED eq 0)], $
	thick = th, $
	charsize=cs,/outline,/freq,$
	binsize = bs, $
	datacolor='blue', $
	xr=[-0.05, 1.05], /xstyle, $
	ytitle='fraction of galaxies', $
    ytickformat='(f4.2)',$
	xtitle='GZ2 ring fraction'
cghistoplot,ring_wfraction_NED[where(ring_na_NED ge 1)], binsize = bs, /outline, datacolor='red',/oplot,/freq, thick=th
;cghistoplot,ring_wfraction_odd[where(ring_na_odd eq 0)], binsize = bs, /outline, datacolor='blue',/oplot,/freq, thick=th+2
;cghistoplot,ring_wfraction_odd[where(ring_na_odd ge 1)], binsize = bs, /outline, datacolor='red',/oplot,/freq, thick=th+2
al_legend,charsize=cs,psym=28,color=['blue','red'],/top,/right,['NA10 no ring','NA10 ring']


; Rings by NA10 type

bs2=0.10
cghistoplot,ring_wfraction_NED[where(ring_na_NED gt 0)], $
	charsize=cs,/outline,$
	thick=th, $
;	/freq,$
	binsize = bs2, $
	datacolor='Red', $
	ytitle='number of galaxies', $
	xr=[-0.05, 1.05], /xstyle, $
	xtitle='GZ2 ring fraction'
na_ring1_NED = where(na_binparse(ring_na_NED,1),c1)
na_ring2_NED = where(na_binparse(ring_na_NED,2),c2)
na_ring3_NED = where(na_binparse(ring_na_NED,3),c3)
if c1 gt 4 then cghistoplot,ring_wfraction_NED[na_ring1_NED], /outline, binsize=bs2, datacolor='black',/oplot, thick=th
if c2 gt 4 then cghistoplot,ring_wfraction_NED[na_ring2_NED], /outline, binsize=bs2, datacolor='dark green',/oplot, thick=th
if c3 gt 4 then cghistoplot,ring_wfraction_NED[na_ring3_NED], /outline, binsize=bs2, datacolor='purple',/oplot, thick=th
;cghistoplot,ring_wfraction_odd[where(ring_na_odd gt 0)], binsize = bs2, /outline, datacolor='red',/oplot, thick=th+2
;na_ring1_odd = where(na_binparse(ring_na_odd,1),c1)
;na_ring2_odd = where(na_binparse(ring_na_odd,2),c2)
;na_ring3_odd = where(na_binparse(ring_na_odd,3),c3)
;if c1 gt 4 then cghistoplot,ring_wfraction_odd[na_ring1_odd], /outline, binsize=bs2, datacolor='black',/oplot, thick=th+1
;if c2 gt 4 then cghistoplot,ring_wfraction_odd[na_ring2_odd], /outline, binsize=bs2, datacolor='dark green',/oplot, thick=th+1
;if c3 gt 4 then cghistoplot,ring_wfraction_odd[na_ring3_odd], /outline, binsize=bs2, datacolor='purple',/oplot, thick=th+1

al_legend,charsize=cs,psym=28,color=['red','black','dark green','purple'],/top,/left,['NA10 all rings','Nuclear','Inner','Outer']


if keyword_set(ps) then ps_end








if keyword_set(stop) then stop

print,''
end

