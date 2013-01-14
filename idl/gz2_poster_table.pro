
;+
; NAME:
;       
;	GZ2_POSTER_TABLE
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

pro gz2_poster_table, ps=ps, stop=stop

device,retain=2

timestart = systime(1)

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = gz2dir+'fits/'
savdir = gz2dir+'sav/'
figsdir = gz2dir+'gz2dropbox/figures/'
gz2tablefile = fitsdir+'gz2table.fits'

	;gz2 = mrdfits(gz2tablefile, 1, /silent)
	;save,gz2,filename=savdir+'gz2table.sav'

restore,savdir+'gz2table.sav'

tagnames=tag_names(gz2)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24

gz2_main            = gz2[where((strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra'),count_main)]
gz2_original        = gz2[where((strtrim(gz2.sample,2) eq 'original'),count_original)]
gz2_extra           = gz2[where((strtrim(gz2.sample,2) eq 'extra'),count_extra)]


taskprob = 0.5
threshprobarr = [0.5,0.8]
hilim = 25
lolim = 10

for i=0,n_elements(threshprobarr)-1 do begin

    threshprob = threshprobarr[i]
    print, ''
    print,'Threshprob = ',threshprob
    
    print,'Task 01'
    a01 = where(gz2_main.t01_smooth_or_features_total_weight ge hilim and gz2_main.t01_smooth_or_features_a01_smooth_weighted_fraction ge threshprob,na01)
    a02 = where(gz2_main.t01_smooth_or_features_total_weight ge hilim and gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge threshprob,na02)
    a03 = where(gz2_main.t01_smooth_or_features_total_weight ge hilim and gz2_main.t01_smooth_or_features_a03_star_or_artifact_weighted_fraction ge threshprob,na03)
    t01 = where(gz2_main.t01_smooth_or_features_total_weight ge hilim,nt01)
    a01_wf = gz2_main[t01].t01_smooth_or_features_a01_smooth_weighted_fraction
    a02_wf = gz2_main[t01].t01_smooth_or_features_a02_features_or_disk_weighted_fraction
    a03_wf = gz2_main[t01].t01_smooth_or_features_a03_star_or_artifact_weighted_fraction
    t01_w = gz2_main[t01].t01_smooth_or_features_total_weight
               
    print, na01, float(na01)/float(nt01) * 100., total(a01_wf), total(a01_wf)/n_elements(t01) * 100.
    print, na02, float(na02)/float(nt01) * 100., total(a02_wf), total(a02_wf)/n_elements(t01) * 100.
    print, na03, float(na03)/float(nt01) * 100., total(a03_wf), total(a03_wf)/n_elements(t01) * 100.


    print,'Task 02'
    a04 = where(gz2_main.t02_edgeon_total_weight ge lolim and gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge taskprob and gz2_main.t02_edgeon_a04_yes_weighted_fraction ge threshprob,na04)
    a05 = where(gz2_main.t02_edgeon_total_weight ge lolim and gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge taskprob and gz2_main.t02_edgeon_a05_no_weighted_fraction ge threshprob,na05)
    t02 = where(gz2_main.t02_edgeon_total_weight ge lolim,nt02)
    a04_wf = gz2_main[t02].t02_edgeon_a04_yes_weighted_fraction
    a05_wf = gz2_main[t02].t02_edgeon_a05_no_weighted_fraction
    t02_w = gz2_main[t02].t02_edgeon_total_weight

    print, na04, float(na04)/float(nt02) * 100., total(a04_wf), total(a04_wf)/n_elements(t02) * 100.
    print, na05, float(na05)/float(nt02) * 100., total(a05_wf), total(a05_wf)/n_elements(t02) * 100.
    
    print,'Task 03'
    a06 = where(gz2_main.t03_bar_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t03_bar_a06_bar_weighted_fraction ge threshprob,na06)
    a07 = where(gz2_main.t03_bar_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t03_bar_a07_no_bar_weighted_fraction ge threshprob,na07)
    t03 = where(gz2_main.t03_bar_total_weight ge lolim,nt03)
    a06_wf = gz2_main[t03].t03_bar_a06_bar_weighted_fraction
    a07_wf = gz2_main[t03].t03_bar_a07_no_bar_weighted_fraction
    t03_w = gz2_main[t03].t03_bar_total_weight

    print, na06, float(na06)/float(nt03) * 100., total(a06_wf), total(a06_wf)/n_elements(t03) * 100.
    print, na07, float(na07)/float(nt03) * 100., total(a07_wf), total(a07_wf)/n_elements(t03) * 100.
    
    print,'Task 04'
    a08 = where(gz2_main.t04_spiral_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge threshprob,na08)
    a09 = where(gz2_main.t04_spiral_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t04_spiral_a09_no_spiral_weighted_fraction ge threshprob,na09)
    t04 = where(gz2_main.t04_spiral_total_weight ge lolim,nt04)
    a08_wf = gz2_main[t04].t04_spiral_a08_spiral_weighted_fraction
    a09_wf = gz2_main[t04].t04_spiral_a09_no_spiral_weighted_fraction
    t04_w = gz2_main[t04].t04_spiral_total_weight

    print, na08, float(na08)/float(nt04) * 100., total(a08_wf), total(a08_wf)/n_elements(t04) * 100.
    print, na09, float(na09)/float(nt04) * 100., total(a09_wf), total(a09_wf)/n_elements(t04) * 100.
    
    
    print,'Task 05'
    a10 = where(gz2_main.t05_bulge_prominence_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t05_bulge_prominence_a10_no_bulge_weighted_fraction ge threshprob,na10)
    a11 = where(gz2_main.t05_bulge_prominence_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t05_bulge_prominence_a11_just_noticeable_weighted_fraction ge threshprob,na11)
    a12 = where(gz2_main.t05_bulge_prominence_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t05_bulge_prominence_a12_obvious_weighted_fraction ge threshprob,na12)
    a13 = where(gz2_main.t05_bulge_prominence_total_weight ge lolim and gz2_main.t02_edgeon_a05_no_weighted_fraction ge taskprob and gz2_main.t05_bulge_prominence_a13_dominant_weighted_fraction ge threshprob,na13)
    t05 = where(gz2_main.t05_bulge_prominence_total_weight ge lolim,nt05)
    a10_wf = gz2_main[t05].t05_bulge_prominence_a10_no_bulge_weighted_fraction
    a11_wf = gz2_main[t05].t05_bulge_prominence_a11_just_noticeable_weighted_fraction
    a12_wf = gz2_main[t05].t05_bulge_prominence_a12_obvious_weighted_fraction
    a13_wf = gz2_main[t05].t05_bulge_prominence_a13_dominant_weighted_fraction
    t05_w = gz2_main[t05].t05_bulge_prominence_total_weight

    print, na10, float(na10)/float(nt05) * 100., total(a10_wf), total(a10_wf)/n_elements(t05) * 100.
    print, na11, float(na11)/float(nt05) * 100., total(a11_wf), total(a11_wf)/n_elements(t05) * 100.
    print, na12, float(na12)/float(nt05) * 100., total(a12_wf), total(a12_wf)/n_elements(t05) * 100.
    print, na13, float(na13)/float(nt05) * 100., total(a13_wf), total(a13_wf)/n_elements(t05) * 100.
    
    
    print,'Task 06'
    a14 = where(gz2_main.t06_odd_total_weight ge hilim and gz2_main.t06_odd_a14_yes_weighted_fraction ge threshprob,na14)
    a15 = where(gz2_main.t06_odd_total_weight ge hilim and gz2_main.t06_odd_a15_no_weighted_fraction ge threshprob,na15)
    t06 = where(gz2_main.t06_odd_total_weight ge hilim,nt06)
    a14_wf = gz2_main[t06].t06_odd_a14_yes_weighted_fraction
    a15_wf = gz2_main[t06].t06_odd_a15_no_weighted_fraction
    t06_w = gz2_main[t06].t06_odd_total_weight

    print, na14, float(na14)/float(nt06) * 100., total(a14_wf), total(a14_wf)/n_elements(t06) * 100.
    print, na15, float(na15)/float(nt06) * 100., total(a15_wf), total(a15_wf)/n_elements(t06) * 100.
    
    
    print,'Task 07'
    a16 = where(gz2_main.t07_rounded_total_weight ge lolim and gz2_main.t01_smooth_or_features_a01_smooth_weighted_fraction ge taskprob and gz2_main.t07_rounded_a16_completely_round_weighted_fraction ge threshprob,na16)
    a17 = where(gz2_main.t07_rounded_total_weight ge lolim and gz2_main.t01_smooth_or_features_a01_smooth_weighted_fraction ge taskprob and gz2_main.t07_rounded_a17_in_between_weighted_fraction ge threshprob,na17)
    a18 = where(gz2_main.t07_rounded_total_weight ge lolim and gz2_main.t01_smooth_or_features_a01_smooth_weighted_fraction ge taskprob and gz2_main.t07_rounded_a18_cigar_shaped_weighted_fraction ge threshprob,na18)
    t07 = where(gz2_main.t07_rounded_total_weight ge lolim,nt07)
    a16_wf = gz2_main[t07].t07_rounded_a16_completely_round_weighted_fraction
    a17_wf = gz2_main[t07].t07_rounded_a17_in_between_weighted_fraction
    a18_wf = gz2_main[t07].t07_rounded_a18_cigar_shaped_weighted_fraction
    t07_w = gz2_main[t07].t07_rounded_total_weight

    print, na16, float(na16)/float(nt07) * 100., total(a16_wf), total(a16_wf)/n_elements(t07) * 100.
    print, na17, float(na17)/float(nt07) * 100., total(a17_wf), total(a17_wf)/n_elements(t07) * 100.
    print, na18, float(na18)/float(nt07) * 100., total(a18_wf), total(a18_wf)/n_elements(t07) * 100.
    
    print,'Task 08'
    a19 = where(gz2_main.t08_odd_feature_total_weight ge lolim and gz2_main.t06_odd_a14_yes_weighted_fraction ge taskprob and gz2_main.t08_odd_feature_a19_ring_weighted_fraction ge threshprob,na19)
    a20 = where(gz2_main.t08_odd_feature_total_weight ge lolim and gz2_main.t06_odd_a14_yes_weighted_fraction ge taskprob and gz2_main.t08_odd_feature_a20_lens_or_arc_weighted_fraction ge threshprob,na20)
    a21 = where(gz2_main.t08_odd_feature_total_weight ge lolim and gz2_main.t06_odd_a14_yes_weighted_fraction ge taskprob and gz2_main.t08_odd_feature_a21_disturbed_weighted_fraction ge threshprob,na21)
    a22 = where(gz2_main.t08_odd_feature_total_weight ge lolim and gz2_main.t06_odd_a14_yes_weighted_fraction ge taskprob and gz2_main.t08_odd_feature_a22_irregular_weighted_fraction ge threshprob,na22)
    a23 = where(gz2_main.t08_odd_feature_total_weight ge lolim and gz2_main.t06_odd_a14_yes_weighted_fraction ge taskprob and gz2_main.t08_odd_feature_a23_other_weighted_fraction ge threshprob,na23)
    a24 = where(gz2_main.t08_odd_feature_total_weight ge lolim and gz2_main.t06_odd_a14_yes_weighted_fraction ge taskprob and gz2_main.t08_odd_feature_a24_merger_weighted_fraction ge threshprob,na24)
    a28 = where(gz2_main.t08_odd_feature_total_weight ge lolim and gz2_main.t06_odd_a14_yes_weighted_fraction ge taskprob and gz2_main.t08_odd_feature_a38_dust_lane_weighted_fraction ge threshprob,na38)
    t08 = where(gz2_main.t08_odd_feature_total_weight ge lolim,nt08)
    a19_wf = gz2_main[t08].t08_odd_feature_a19_ring_weighted_fraction
    a20_wf = gz2_main[t08].t08_odd_feature_a20_lens_or_arc_weighted_fraction
    a21_wf = gz2_main[t08].t08_odd_feature_a21_disturbed_weighted_fraction
    a22_wf = gz2_main[t08].t08_odd_feature_a22_irregular_weighted_fraction
    a23_wf = gz2_main[t08].t08_odd_feature_a23_other_weighted_fraction
    a24_wf = gz2_main[t08].t08_odd_feature_a24_merger_weighted_fraction
    a38_wf = gz2_main[t08].t08_odd_feature_a38_dust_lane_weighted_fraction
    t08_w = gz2_main[t08].t08_odd_feature_total_weight

    print, na19, float(na19)/float(nt08) * 100., total(a19_wf), total(a19_wf)/n_elements(t08) * 100.
    print, na20, float(na20)/float(nt08) * 100., total(a20_wf), total(a20_wf)/n_elements(t08) * 100.
    print, na21, float(na21)/float(nt08) * 100., total(a21_wf), total(a21_wf)/n_elements(t08) * 100.
    print, na22, float(na22)/float(nt08) * 100., total(a22_wf), total(a22_wf)/n_elements(t08) * 100.
    print, na23, float(na23)/float(nt08) * 100., total(a23_wf), total(a23_wf)/n_elements(t08) * 100.
    print, na24, float(na24)/float(nt08) * 100., total(a24_wf), total(a24_wf)/n_elements(t08) * 100.
    print, na38, float(na38)/float(nt08) * 100., total(a38_wf), total(a38_wf)/n_elements(t08) * 100.
    
    
    print,'Task 09'
    a25 = where(gz2_main.t09_bulge_shape_total_weight ge lolim and gz2_main.t02_edgeon_a04_yes_weighted_fraction ge taskprob and gz2_main.t09_bulge_shape_a25_rounded_weighted_fraction ge threshprob,na25)
    a26 = where(gz2_main.t09_bulge_shape_total_weight ge lolim and gz2_main.t02_edgeon_a04_yes_weighted_fraction ge taskprob and gz2_main.t09_bulge_shape_a26_boxy_weighted_fraction ge threshprob,na26)
    a27 = where(gz2_main.t09_bulge_shape_total_weight ge lolim and gz2_main.t02_edgeon_a04_yes_weighted_fraction ge taskprob and gz2_main.t09_bulge_shape_a27_no_bulge_weighted_fraction ge threshprob,na27)
    t09 = where(gz2_main.t09_bulge_shape_total_weight ge lolim,nt09)
    a25_wf = gz2_main[t09].t09_bulge_shape_a25_rounded_weighted_fraction
    a26_wf = gz2_main[t09].t09_bulge_shape_a26_boxy_weighted_fraction
    a27_wf = gz2_main[t09].t09_bulge_shape_a27_no_bulge_weighted_fraction
    t09_w = gz2_main[t09].t09_bulge_shape_total_weight

    print, na25, float(na25)/float(nt09) * 100., total(a25_wf), total(a25_wf)/n_elements(t09) * 100.
    print, na26, float(na26)/float(nt09) * 100., total(a26_wf), total(a26_wf)/n_elements(t09) * 100.
    print, na27, float(na27)/float(nt09) * 100., total(a27_wf), total(a27_wf)/n_elements(t09) * 100.
    
    print,'Task 10'
    a28 = where(gz2_main.t10_arms_winding_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t10_arms_winding_a28_tight_weighted_fraction ge threshprob,na28)
    a29 = where(gz2_main.t10_arms_winding_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t10_arms_winding_a29_medium_weighted_fraction ge threshprob,na29)
    a30 = where(gz2_main.t10_arms_winding_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t10_arms_winding_a30_loose_weighted_fraction ge threshprob,na30)
    t10 = where(gz2_main.t10_arms_winding_total_weight ge lolim,nt10)
    a28_wf = gz2_main[t10].t10_arms_winding_a28_tight_weighted_fraction
    a29_wf = gz2_main[t10].t10_arms_winding_a29_medium_weighted_fraction
    a30_wf = gz2_main[t10].t10_arms_winding_a30_loose_weighted_fraction
    t10_w = gz2_main[t10].t10_arms_winding_total_weight

    print, na28, float(na28)/float(nt10) * 100., total(a28_wf), total(a28_wf)/n_elements(t10) * 100.
    print, na29, float(na29)/float(nt10) * 100., total(a29_wf), total(a29_wf)/n_elements(t10) * 100.
    print, na30, float(na30)/float(nt10) * 100., total(a30_wf), total(a30_wf)/n_elements(t10) * 100.
    
    print,'Task 11'
    a31 = where(gz2_main.t11_arms_number_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t11_arms_number_a31_1_weighted_fraction ge threshprob,na31)
    a32 = where(gz2_main.t11_arms_number_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t11_arms_number_a32_2_weighted_fraction ge threshprob,na32)
    a33 = where(gz2_main.t11_arms_number_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t11_arms_number_a33_3_weighted_fraction ge threshprob,na33)
    a34 = where(gz2_main.t11_arms_number_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t11_arms_number_a34_4_weighted_fraction ge threshprob,na34)
    a36 = where(gz2_main.t11_arms_number_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t11_arms_number_a36_more_than_4_weighted_fraction ge threshprob,na36)
    a37 = where(gz2_main.t11_arms_number_total_weight ge lolim and gz2_main.t04_spiral_a08_spiral_weighted_fraction ge taskprob and gz2_main.t11_arms_number_a37_cant_tell_weighted_fraction ge threshprob,na37)
    t11 = where(gz2_main.t11_arms_number_total_weight ge lolim,nt11)
    a31_wf = gz2_main[t11].t11_arms_number_a31_1_weighted_fraction
    a32_wf = gz2_main[t11].t11_arms_number_a32_2_weighted_fraction
    a33_wf = gz2_main[t11].t11_arms_number_a33_3_weighted_fraction
    a34_wf = gz2_main[t11].t11_arms_number_a34_4_weighted_fraction
    a36_wf = gz2_main[t11].t11_arms_number_a36_more_than_4_weighted_fraction
    a37_wf = gz2_main[t11].t11_arms_number_a37_cant_tell_weighted_fraction
    t11_w = gz2_main[t11].t11_arms_number_total_weight

    print, na31, float(na31)/float(nt11) * 100., total(a31_wf), total(a31_wf)/n_elements(t11) * 100.
    print, na32, float(na32)/float(nt11) * 100., total(a32_wf), total(a32_wf)/n_elements(t11) * 100.
    print, na33, float(na33)/float(nt11) * 100., total(a33_wf), total(a33_wf)/n_elements(t11) * 100.
    print, na34, float(na34)/float(nt11) * 100., total(a34_wf), total(a34_wf)/n_elements(t11) * 100.
    print, na36, float(na36)/float(nt11) * 100., total(a36_wf), total(a36_wf)/n_elements(t11) * 100.
    print, na37, float(na37)/float(nt11) * 100., total(a37_wf), total(a37_wf)/n_elements(t11) * 100.
    
    print,''

endfor





if keyword_set(stop) then stop

end

