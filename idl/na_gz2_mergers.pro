;+
; NAME:
;       
;	NA_GZ2_MERGERS
;
; PURPOSE:
;
;	Compare merger classifications from Nair & Abraham vs. Galaxy Zoo 2
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

pro na_gz2_mergers, ps=ps, stop=stop

; Paths for files and plots

csv_path = '~/Astronomy/Research/GalaxyZoo/csv/'
fig_path = '~/Astronomy/Research/GalaxyZoo/datapaper/figures/'
fits_path = '~/Astronomy/Research/GalaxyZoo/fits/'

nagz2 = mrdfits(fits_path+'na_gz2_willettk.fit',1,/silent)

; Limits on previous questions

galaxy = where(nagz2.t01_smooth_or_features_a03_star_or_artifact_weighted_fraction le 0.05)
galaxy_odd = where(nagz2.t01_smooth_or_features_a03_star_or_artifact_weighted_fraction le 0.05 and nagz2.t06_odd_a14_yes_weighted_fraction ge 0.223)
galaxy_notodd = where(nagz2.t01_smooth_or_features_a03_star_or_artifact_weighted_fraction le 0.05 and nagz2.t06_odd_a14_yes_weighted_fraction lt 0.223)

; How does the number of detected merging galaxies depend on cutting on Task 06?

print,'no odd cut'
print, string('interacting',format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,1))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,1)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, string('none',       format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,2))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,2)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, string('disturbed',  format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,3))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,3)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, string('warp',       format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,4))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,4)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, string('shells',     format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,5))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,5)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, string('tail short', format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,6))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,6)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, string('tail med',   format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,7))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,7)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, string('tail long',  format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy].dist,8))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy].dist,8)))/float(n_elements(where(nagz2[galaxy].dist gt 0)))
print, n_elements(where(nagz2[galaxy].dist gt 0))
print,''
print,'odd'
print, string('interacting',format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,1))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,1)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, string('none',       format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,2))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,2)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, string('disturbed',  format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,3))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,3)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, string('warp',       format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,4))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,4)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, string('shells',     format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,5))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,5)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, string('tail short', format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,6))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,6)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, string('tail med',   format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,7))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,7)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, string('tail long',  format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_odd].dist,8))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_odd].dist,8)))/float(n_elements(where(nagz2[galaxy_odd].dist gt 0)))
print, n_elements(where(nagz2[galaxy_odd].dist gt 0))
print,''
print,'not odd'
print, string('interacting',format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,1))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,1)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, string('none',       format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,2))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,2)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, string('disturbed',  format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,3))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,3)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, string('warp',       format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,4))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,4)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, string('shells',     format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,5))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,5)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, string('tail short', format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,6))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,6)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, string('tail med',   format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,7))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,7)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, string('tail long',  format='(a12)'), string(n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,8))),format='(i5)'),n_elements(where(na_binparse(nagz2[galaxy_notodd].dist,8)))/float(n_elements(where(nagz2[galaxy_notodd].dist gt 0)))
print, n_elements(where(nagz2[galaxy_notodd].dist gt 0))
print,''

nagz2_galaxy = nagz2[galaxy]
nagz2_odd = nagz2[galaxy_odd]

dist_na_galaxy = nagz2_galaxy.dist
merger_count_galaxy = nagz2_galaxy.t08_odd_feature_a24_merger_weight
merger_wfraction_galaxy = nagz2_galaxy.t08_odd_feature_a24_merger_weighted_fraction

dist_na_odd = nagz2_odd.dist
merger_count_odd = nagz2_odd.t08_odd_feature_a24_merger_weight
merger_wfraction_odd = nagz2_odd.t08_odd_feature_a24_merger_weighted_fraction

pair_na_galaxy = nagz2_galaxy.pair
pair_na_odd = nagz2_odd.pair

; Second merger plot

if keyword_set(ps) then begin
	ps_start, filename=fig_path+'na_gz2_dist.eps',/color, /quiet, xs=5, ys=4, /encap
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

!p.multi=[0,1,1]

; Disturbed galaxies

bs = 0.05
cghistoplot,merger_wfraction_odd[where(dist_na_odd eq 0)], $
    /nodata, $
	thick = th, $
	charsize=cs,/outline,$
    ;/freq,$
	binsize = bs, $
	datacolor='blue', $
;    yr=[0,0.1], $
    yr=[0,200], $
	xr=[-0.05, 1.05], /xstyle, $
	ytitle='number of galaxies', $
;    ytickformat='(f4.2)',$
	xtitle='GZ2 vote fraction for merger'
cghistoplot,merger_wfraction_odd[where(dist_na_odd ge 1)], binsize = bs, /outline, datacolor='red',/oplot, thick=th;,/freq, thick=th
cghistoplot,merger_wfraction_odd[where(pair_na_odd ge 1)], binsize = bs, /outline, datacolor='black',/oplot, thick=th;,/freq, thick=th
;cghistoplot,merger_wfraction_galaxy[where(dist_na_galaxy eq 0)], binsize = bs, /outline, datacolor='blue',/oplot,/freq, thick=th+2
;cghistoplot,merger_wfraction_galaxy[where(dist_na_galaxy ge 1)], binsize = bs, /outline, datacolor='red',/oplot,/freq, thick=th+2
;cghistoplot,merger_wfraction_galaxy[where(pair_na_galaxy ge 1)], binsize = bs, /outline, datacolor='black',/oplot,/freq, thick=th+2
al_legend,charsize=cs,psym=28,color=['blue','red','black'],/top,/right,['NA10 undisturbed','NA10 disturbed', 'NA10 pairs']

; NA10 interacting flags are screwed up (must have at least 11 categories, but only include 8). Also don't know difference between dist = 0 and dist = 1. 
; Can ask Nair in person sometime - for now, simply compare the paired and interacting fractions.

;; Disturbed galaxies by NA10 type
;
;bs2=0.10
;cghistoplot,merger_wfraction_galaxy[where(dist_na_galaxy gt 0)], $
;	charsize=cs,/outline,$
;	thick=th, $
;;	/freq,$
;	binsize = bs2, $
;	datacolor='Red', $
;	ytitle='number of galaxies', $
;	xr=[-0.05, 1.05], /xstyle, $
;	xtitle='GZ2 merger fraction'
;dist_na_pair1 = where(na_binparse(dist_na_galaxy,1),c1)
;dist_na_pair2 = where(na_binparse(dist_na_galaxy,2),c2)
;dist_na_pair3 = where(na_binparse(dist_na_galaxy,3),c3)
;dist_na_pair4 = where(na_binparse(dist_na_galaxy,4),c4)
;dist_na_pair5 = where(na_binparse(dist_na_galaxy,5),c5)
;dist_na_pair6 = where(na_binparse(dist_na_galaxy,6),c6)
;dist_na_pair7 = where(na_binparse(dist_na_galaxy,7),c7)
;dist_na_pair8 = where(na_binparse(dist_na_galaxy,8),c8)
;if c1 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair1], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='blue'
;if c2 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair2], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='green'
;if c3 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair3], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='purple'
;if c4 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair4], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='red'
;if c5 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair5], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='yellow'
;if c6 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair6], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='orange'
;if c7 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair7], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='hotpink'
;if c8 gt 4 then cghistoplot,merger_wfraction_galaxy[dist_na_pair8], /outline, /oplot, thick=th, binsize=bs2 ,  datacolor='black'
;
;al_legend,charsize=cs/2.,psym=28,color=['black','blue','green','purple','red','yellow','orange','hotpink','black'],/top,/right,['Interacting','None','Disturbed','Warp','Shells','Short tail','Medium tail','Long tail','Bridge']


if keyword_set(ps) then ps_end








if keyword_set(stop) then stop

print,''
end


