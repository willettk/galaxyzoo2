
;+
; NAME:
;       
;	GZ2_TASK_HISTOGRAMS
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

pro gz2_task_histograms, ps=ps, stop=stop

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
gz2_stripe82_normal = gz2[where((strtrim(gz2.sample,2) eq 'stripe82'),count_stripe82_normal)]
gz2_stripe82_coadd1 = gz2[where((strtrim(gz2.sample,2) eq 'stripe82_coadd_1'),count_stripe82_coadd1)]
gz2_stripe82_coadd2 = gz2[where((strtrim(gz2.sample,2) eq 'stripe82_coadd_2'),count_stripe82_coadd2)]

if keyword_set(ps) then begin
	ps_start, filename=figsdir+'gz2_task_histograms.ps', /color, /quiet
	cs=1.2
	legendcs = 0.9
	th=3
	thickline=8
	thinline=3
endif else begin
	cs=2
	legendcs = cs
	th=1
	thickline=2
	thinline=1
endelse

color_main = 'black'
color_original = 'red'
color_extra = 'orange'
color_stripe82 = 'blue'
color_stripe82_coadd1 = 'green'
color_stripe82_coadd2 = 'dark green'

task_count_ind_arr = [3,5,7,9,13,15,18,25,28,31,37]
title_arr = ['Task 01 (smooth or features)', $
             'Task 02 (edge-on disk)', $
             'Task 03 (bar)', $
             'Task 04 (spiral structure)', $
             'Task 05 (bulge dominance)', $
             'Task 06 (anything odd)', $
             'Task 07 (rounded)', $
             'Task 08 (odd feature)', $
             'Task 09 (bulge shape)', $
             'Task 10 (arms winding)', $
             'Task 11 (arms number)']
ntask = n_elements(task_count_ind_arr)

;for i=0,ntask-1 do begin
;
;    !p.multi=[0,1,1]
;
;    task_count_ind = task_count_ind_arr[i]
;    
;    data = gz2_main.(wfind[task_count_ind]+2)
;    sidx = sort(data)
;    ndata = n_elements(data)
;    pc1  = data[sidx[1*ndata / 100]]
;    pc5  = data[sidx[5*ndata / 100]]
;    pc10 = data[sidx[10*ndata / 100]]
;
;    ;print,i
;    ;print,'1%', pc1
;    ;print,'5%', pc5
;    ;print,'10%', pc10
;    
;    bs=1
;    
;    cghistoplot, gz2_main.(wfind[task_count_ind]+2), $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_main,           $
;        thick=thickline,                $
;        title=title_arr[i], $
;        xtitle='Number of classifications per object', $
;        xr=[0,75], $
;        yr = [0,0.2]
;    
;    cgplots, replicate(pc1 ,2), !y.crange, thick=thinline, color=color_main, linestyle=3
;    cgplots, replicate(pc5 ,2), !y.crange, thick=thinline, color=color_main, linestyle=1
;    cgplots, replicate(pc10,2), !y.crange, thick=thinline, color=color_main, linestyle=2
;
;    cghistoplot, gz2_original.(wfind[task_count_ind]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_original, $
;        thick=thinline
;    cghistoplot, gz2_extra.(wfind[task_count_ind]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_extra, $
;        thick=thinline
;    cghistoplot, gz2_stripe82_normal.(wfind[task_count_ind]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_stripe82, $
;        thick=thinline
;    cghistoplot, gz2_stripe82_coadd1.(wfind[task_count_ind]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_stripe82_coadd1, $
;        thick=thinline
;    cghistoplot, gz2_stripe82_coadd2.(wfind[task_count_ind]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_stripe82_coadd2, $
;        thick=thinline
;    
;    al_legend,/top,/right,$
;        ['Original + Extra','Original','Extra','Stripe 82 normal','Stripe 82 coadd 1','Stripe 82 coadd 2'],$ 
;        color=[color_main, color_original, color_extra, color_stripe82, color_stripe82_coadd1, color_stripe82_coadd2], $
;        psym=28, $
;        symsize=2, $
;        charsize = legendcs
;
;endfor
	
if keyword_set(ps) then ps_end



    ; Look at thresholds and previous probs
    !p.multi=[0,3,3]
    wfarr = (indgen(9))/10.

    t01lim=10
    for i=0,n_elements(wfarr)-1 do begin
        cghistoplot,gz2_main[where(gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge wfarr[i] and gz2_main.t01_smooth_or_features_total_weight ge t01lim)].t02_edgeon_total_weight,$
        binsize=1,/outline,yr=[0,25000],title='Task 01 threshold on Task 02: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]

        cgtext, 5,20000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge wfarr[i] and gz2_main.t01_smooth_or_features_total_weight ge t01lim)].t02_edgeon_total_weight le t01lim)) / float(n_elements(where(gz2_main[where(gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge wfarr[0] and gz2_main.t01_smooth_or_features_total_weight ge t01lim)].t02_edgeon_total_weight le t01lim))),2)
        cgtext, 40,20000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge wfarr[i] and gz2_main.t01_smooth_or_features_total_weight ge t01lim)].t02_edgeon_total_weight gt t01lim)) / float(n_elements(where(gz2_main[where(gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction ge wfarr[0] and gz2_main.t01_smooth_or_features_total_weight ge t01lim)].t02_edgeon_total_weight gt t01lim))),2)

        endfor

    print,'Task 02'
    print,find99(gz2_main.t02_edgeon_total_weight,gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction,gz2_main.t01_smooth_or_features_total_weight,10)
    print,find99(gz2_main.t02_edgeon_total_weight,gz2_main.t01_smooth_or_features_a02_features_or_disk_weighted_fraction,gz2_main.t01_smooth_or_features_total_weight,20)
    print,''


; count number of galaxies lost to the left of the peak

    ; Keep same limit for Tasks 03 and 04 - tree doesn't branch
    t02lim=20
    for i=0,n_elements(wfarr)-1 do begin
        cghistoplot,gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[i] $ 
        and gz2_main.t02_edgeon_total_weight ge t02lim)].t03_bar_total_weight,$
        binsize=1,/outline,yr=[0,10000],title='Task 02 threshold on Task 03: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]
        cgtext, 5,8000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[i] and gz2_main.t02_edgeon_total_weight ge t02lim)].t03_bar_total_weight le t02lim)) / float(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[0] and gz2_main.t02_edgeon_total_weight ge t02lim)].t03_bar_total_weight le t02lim))),2)
        cgtext, 40,8000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[i] and gz2_main.t02_edgeon_total_weight ge t02lim)].t03_bar_total_weight gt t02lim)) / float(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[0] and gz2_main.t02_edgeon_total_weight ge t02lim)].t03_bar_total_weight gt t02lim))),2)

    endfor

    print,'Task 03'
    print,find99(gz2_main.t03_bar_total_weight,gz2_main.t02_edgeon_a05_no_weighted_fraction,gz2_main.t02_edgeon_total_weight,10)
    print,find99(gz2_main.t03_bar_total_weight,gz2_main.t02_edgeon_a05_no_weighted_fraction,gz2_main.t02_edgeon_total_weight,20)
    print,''


    for i=0,n_elements(wfarr)-1 do $
        cghistoplot,gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[i] $ 
        and gz2_main.t02_edgeon_total_weight ge t02lim)].t04_spiral_total_weight,$
        binsize=1,/outline,yr=[0,10000],title='Task 02 threshold on Task 04: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]

    print,'Task 04'
    print,find99(gz2_main.t04_spiral_total_weight,gz2_main.t02_edgeon_a05_no_weighted_fraction,gz2_main.t02_edgeon_total_weight,10)
    print,find99(gz2_main.t04_spiral_total_weight,gz2_main.t02_edgeon_a05_no_weighted_fraction,gz2_main.t02_edgeon_total_weight,20)
    print,''


    for i=0,n_elements(wfarr)-1 do begin
        cghistoplot,gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[i] $ 
        and gz2_main.t02_edgeon_total_weight ge t02lim)].t05_bulge_prominence_total_weight,$
        binsize=1,/outline,yr=[0,10000],title='Task 02 threshold on Task 03: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]
        cgtext, 5,8000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[i] and gz2_main.t02_edgeon_total_weight ge t02lim)].t05_bulge_prominence_total_weight le t02lim)) / float(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[0] and gz2_main.t02_edgeon_total_weight ge t02lim)].t05_bulge_prominence_total_weight le t02lim))),2)
        cgtext, 40,8000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[i] and gz2_main.t02_edgeon_total_weight ge t02lim)].t05_bulge_prominence_total_weight gt t02lim)) / float(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a05_no_weighted_fraction ge wfarr[0] and gz2_main.t02_edgeon_total_weight ge t02lim)].t05_bulge_prominence_total_weight gt t02lim))),2)

    endfor

    print,'Task 05'
    print,find99(gz2_main.t05_bulge_prominence_total_weight,gz2_main.t02_edgeon_a05_no_weighted_fraction,gz2_main.t02_edgeon_total_weight,10)
    print,find99(gz2_main.t05_bulge_prominence_total_weight,gz2_main.t02_edgeon_a05_no_weighted_fraction,gz2_main.t02_edgeon_total_weight,20)
    print,''


    ; Keep same limit for Tasks 10 and 11 - tree doesn't branch
    t04lim=20
    for i=0,n_elements(wfarr)-1 do begin
        cghistoplot,gz2_main[where(gz2_main.t04_spiral_a08_spiral_weighted_fraction ge wfarr[i] $ 
        and gz2_main.t04_spiral_total_weight ge t04lim)].t10_arms_winding_total_weight,$
        binsize=1,/outline,yr=[0,5000],title='Task 04 threshold on Task 10: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]
        cgtext, 5,4000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t04_spiral_a08_spiral_weighted_fraction ge wfarr[i] and gz2_main.t04_spiral_total_weight ge t04lim)].t10_arms_winding_total_weight le t04lim)) / float(n_elements(where(gz2_main[where(gz2_main.t04_spiral_a08_spiral_weighted_fraction ge wfarr[0] and gz2_main.t04_spiral_total_weight ge t04lim)].t10_arms_winding_total_weight le t04lim))),2)
        cgtext, 40,4000,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t04_spiral_a08_spiral_weighted_fraction ge wfarr[i] and gz2_main.t04_spiral_total_weight ge t04lim)].t10_arms_winding_total_weight gt t04lim)) / float(n_elements(where(gz2_main[where(gz2_main.t04_spiral_a08_spiral_weighted_fraction ge wfarr[0] and gz2_main.t04_spiral_total_weight ge t04lim)].t10_arms_winding_total_weight gt t04lim))),2)
    endfor

    print,'Task 10'
    print,find99(gz2_main.t10_arms_winding_total_weight,gz2_main.t04_spiral_a08_spiral_weighted_fraction,gz2_main.t04_spiral_total_weight,10)
    print,find99(gz2_main.t10_arms_winding_total_weight,gz2_main.t04_spiral_a08_spiral_weighted_fraction,gz2_main.t04_spiral_total_weight,20)
    print,''



    for i=0,n_elements(wfarr)-1 do $
        cghistoplot,gz2_main[where(gz2_main.t04_spiral_a08_spiral_weighted_fraction ge wfarr[i] $ 
        and gz2_main.t04_spiral_total_weight ge t04lim)].t11_arms_number_total_weight,$
        binsize=1,/outline,yr=[0,10000],title='Task 04 threshold on Task 11: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]

    print,'Task 10'
    print,find99(gz2_main.t11_arms_number_total_weight,gz2_main.t04_spiral_a08_spiral_weighted_fraction,gz2_main.t04_spiral_total_weight,10)
    print,find99(gz2_main.t11_arms_number_total_weight,gz2_main.t04_spiral_a08_spiral_weighted_fraction,gz2_main.t04_spiral_total_weight,20)
    print,''


    for i=0,n_elements(wfarr)-1 do begin
        cghistoplot,gz2_main[where(gz2_main.t02_edgeon_a04_yes_weighted_fraction ge wfarr[i] $ 
        and gz2_main.t02_edgeon_total_weight ge t02lim)].t09_bulge_shape_total_weight,$
        binsize=1,/outline,yr=[0,2000],title='Task 02 threshold on Task 09: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]
        cgtext, 5,1800,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a04_yes_weighted_fraction ge wfarr[i] and gz2_main.t02_edgeon_total_weight ge t02lim)].t09_bulge_shape_total_weight le t02lim)) / float(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a04_yes_weighted_fraction ge wfarr[0] and gz2_main.t02_edgeon_total_weight ge t02lim)].t09_bulge_shape_total_weight le t02lim))),2)
        cgtext, 40,1800,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a04_yes_weighted_fraction ge wfarr[i] and gz2_main.t02_edgeon_total_weight ge t02lim)].t09_bulge_shape_total_weight gt t02lim)) / float(n_elements(where(gz2_main[where(gz2_main.t02_edgeon_a04_yes_weighted_fraction ge wfarr[0] and gz2_main.t02_edgeon_total_weight ge t02lim)].t09_bulge_shape_total_weight gt t02lim))),2)
        endfor

    print,'Task 09'
    print,find99(gz2_main.t09_bulge_shape_total_weight,gz2_main.t02_edgeon_a04_yes_fraction,gz2_main.t02_edgeon_total_weight,10)
    print,find99(gz2_main.t09_bulge_shape_total_weight,gz2_main.t02_edgeon_a04_yes_fraction,gz2_main.t02_edgeon_total_weight,20)
    print,''

    for i=0,n_elements(wfarr)-1 do $
        cghistoplot,gz2_main[where(gz2_main.t01_smooth_or_features_a01_smooth_weighted_fraction ge wfarr[i] $ 
        and gz2_main.t01_smooth_or_features_total_weight ge t01lim)].t07_rounded_total_weight,$
        binsize=1,/outline,yr=[0,15000],title='Task 01 threshold on Task 07: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]

    print,'Task 07'
    print,find99(gz2_main.t07_rounded_total_weight,gz2_main.t01_smooth_or_features_a01_smooth_fraction,gz2_main.t01_smooth_or_features_total_weight,10)
    print,find99(gz2_main.t07_rounded_total_weight,gz2_main.t01_smooth_or_features_a01_smooth_fraction,gz2_main.t01_smooth_or_features_total_weight,20)
    print,''


    ;for i=0,n_elements(wfarr)-1 do $
    ;    cghistoplot,gz2_main[where(gz2_main.t01_smooth_or_features_total_weight ge t01lim)].t06_odd_total_weight,$
    ;    binsize=1,/outline,yr=[0,25000],title='Task 01 threshold on Task 06: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]

    t08lim = 10
    for i=0,n_elements(wfarr)-1 do begin
        cghistoplot,gz2_main[where(gz2_main.t06_odd_a14_yes_weighted_fraction ge wfarr[i] $
        and gz2_main.t06_odd_total_weight ge t08lim)].t08_odd_feature_total_weight,$
        binsize=1,/outline,yr=[0,35000],title='Task 06 threshold on Task 08: '+strtrim(wfarr[i],2),charsize=2,xr=[0,60]
        cgtext, 5,1800,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t06_odd_a14_yes_weighted_fraction ge wfarr[i] and gz2_main.t06_odd_total_weight ge t08lim)].t08_odd_feature_total_weight le t08lim)) / float(n_elements(where(gz2_main[where(gz2_main.t06_odd_a14_yes_weighted_fraction ge wfarr[0] and gz2_main.t06_odd_total_weight ge t08lim)].t08_odd_feature_total_weight le t08lim))),2)
        cgtext, 40,1800,/data, strtrim(n_elements(where(gz2_main[where(gz2_main.t06_odd_a14_yes_weighted_fraction ge wfarr[i] and gz2_main.t06_odd_total_weight ge t08lim)].t08_odd_feature_total_weight gt t08lim)) / float(n_elements(where(gz2_main[where(gz2_main.t06_odd_a14_yes_weighted_fraction ge wfarr[0] and gz2_main.t06_odd_total_weight ge t08lim)].t08_odd_feature_total_weight gt t08lim))),2)
        endfor

    print,'Task 08'
    print,find99(gz2_main.t08_odd_feature_total_weight,gz2_main.t06_odd_a14_yes_weighted_fraction,gz2_main.t06_odd_total_weight,10)
    print,find99(gz2_main.t08_odd_feature_total_weight,gz2_main.t06_odd_a14_yes_weighted_fraction,gz2_main.t06_odd_total_weight,20)

    stop


; Plot the difference in the total counts between Task 01 and Task 06. They should be identical,
; since these are the only tasks for both classifications. Odd that they are not - ask Steven
; for clarification. 

;ps_start, filename=figsdir+'gz2_task01_task06_diff.eps', /color, /quiet, xs=6, ys=6, /encap
;
;    cghistoplot, gz2_main.(wfind[task_count_ind_arr[0]]+2) - gz2_main.(wfind[task_count_ind_arr[5]]+2), $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_main,           $
;        thick=thickline,                $
;        xr=[0,30], $
;        yr = [0,0.45], $
;        ytickformat='(f3.1)', $
;        charsize = 1.0, $
;        title='Task 01 (smooth) - Task 06 (anything odd)', $
;        xtitle='difference in total counts'
;
;    cghistoplot, gz2_original.(wfind[task_count_ind_arr[0]]+2) - gz2_original.(wfind[task_count_ind_arr[5]]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_original,$
;        thick=thickline                
;
;    cghistoplot, gz2_extra.(wfind[task_count_ind_arr[0]]+2) - gz2_extra.(wfind[task_count_ind_arr[5]]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_extra,$
;        thick=thickline                
;
;    cghistoplot, gz2_stripe82_normal.(wfind[task_count_ind_arr[0]]+2) - gz2_stripe82_normal.(wfind[task_count_ind_arr[5]]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_stripe82,$
;        thick=thickline                
;
;    cghistoplot, gz2_stripe82_coadd1.(wfind[task_count_ind_arr[0]]+2) - gz2_stripe82_coadd1.(wfind[task_count_ind_arr[5]]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_stripe82_coadd1,$
;        thick=thickline                
;
;    cghistoplot, gz2_stripe82_coadd2.(wfind[task_count_ind_arr[0]]+2) - gz2_stripe82_coadd2.(wfind[task_count_ind_arr[5]]+2), $
;        /oplot, $
;        /outline, $
;        /freq, $
;        binsize = bs, $
;        datacolor=color_stripe82_coadd2,$
;        thick=thickline                
;
;ps_end
;
;


; End time

timeend = systime(1)

print,''
print,'Time elapsed: '+string(timeend-timestart,format='(f5.1)')+' sec.'
print,''

if keyword_set(stop) then stop

end
