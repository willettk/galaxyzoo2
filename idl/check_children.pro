
;+
; NAME:
;       
;	CHECK_CHILDREN
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
;	All hail J.D. Smith and HIST_ND.pro
;
;   Can also change the input file to compare any two sets of data with matched PDFs (such as Stripe 82)
;
; REVISION HISTORY
;       Written by K. Willett                Dec 2012
;-

function drawrand_uniform, nnum, xrange, keepseed = keepseed

    ;;; Initialize seed variable. 

    if ~keyword_set(keepseed) then begin
            dum = randomu(seed)
            defsysv,'!myseed',seed
    endif else seed=!myseed
    ret=randomu(seed,nnum)*(max(xrange)-min(xrange))+min(xrange)
    !myseed=seed
    return,ret
end

pro check_children, votemin, ps=ps, stop=stop

device,retain=2

timestart = systime(1)

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


gz2_original        = gz2[where((strtrim(gz2.sample,2) eq 'original'),count_original)]
gz2_extra           = gz2[where((strtrim(gz2.sample,2) eq 'extra'),count_extra)]

color_original = 'red'
color_extra = 'blue'
color_selected = 'dark green'

if keyword_set(ps) then begin
    ps_start, filename='~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/check_children.eps',$
        /encap,/color,xs=16,ys=10, /quiet
    cs = 3.2
    thickline=7
    thinline=2
endif else begin
    cs = 3
    thickline=2
    thinline=1
endelse
!p.multi=[0,3,2]

    r50_minval = 0.
    r50_maxval = 15.
    r50_binsize = 0.50

    ur_minval = 0.
    ur_maxval = 6.
    ur_binsize = 0.15

    rmag_minval = 10.
    rmag_maxval = 18.
    rmag_binsize = 0.25

    ur_original_histdata = histogram((gz2_original.petromag_u - gz2_original.petromag_r), binsize=ur_binsize, min=ur_minval,max=ur_maxval)
    ur_extra_histdata = histogram((gz2_extra.petromag_u - gz2_extra.petromag_r), binsize=ur_binsize, min=ur_minval,max=ur_maxval)
    r50_original_histdata = histogram(gz2_original.petror50_r, binsize=r50_binsize, min=r50_minval,max=r50_maxval)
    r50_extra_histdata    = histogram(gz2_extra.petror50_r   , binsize=r50_binsize, min=r50_minval,max=r50_maxval)

; Last thing to try. Create a 2D histogram of the original sample galaxies in (u-r) and R50. 
; Locate the closest match for an "extra" galaxy in the "original" sample based on color and size

data0 = gz2_original.petromag_u - gz2_original.petromag_r
data1 = gz2_original.petror50_r
data2 = gz2_original.petromag_r
vo = transpose([[data0],[data1],[data2]])
data0 = gz2_extra.petromag_u - gz2_extra.petromag_r
data1 = gz2_extra.petror50_r
data2 = gz2_extra.petromag_r
ve = transpose([[data0],[data1],[data2]])
binsizes = [ur_binsize,r50_binsize,rmag_binsize]
mins = [ur_minval,r50_minval,rmag_minval]
maxs = [ur_maxval,r50_maxval,rmag_maxval]

; Bin the data along the observational variables

hist_original = hist_nd(vo,binsizes,min=mins,max=maxs,reverse_indices=ri_original)
hist_extra    = hist_nd(ve,binsizes,min=mins,max=maxs,reverse_indices=ri_extra)

; Locate bins where the original sample has enough galaxies to draw from and recreate the extra sample demographics

oe = where(hist_original ge hist_extra and hist_extra gt 0,oecount)
sind=0
newinds= lonarr(total(hist_extra[oe]))

for i=0,oecount-1 do begin

    nnew = hist_extra[oe[i]]
    inds=ri_original[ri_original[oe[i]]:ri_original[oe[i]+1]-1]
    indind = fix(drawrand_uniform(nnew,[0,nnew]))
    newinds[sind:sind+nnew-1] = inds[indind] 
    sind += nnew

endfor

out = unique(newinds,ucount,/sort)
print,'Unique galaxies in selected sub-sample from "original":',n_elements(out)

    cghistoplot, gz2_original.petromag_u - gz2_original.petromag_r, $
        /outline, $
        /freq, $
        bin=ur_binsize, $
        mininput = ur_minval, $
        maxinput = ur_maxval, $
        charsize=cs, $
        ytickformat='(f4.2)', $
        xr = [ur_minval,ur_maxval], $
        yr=[0,0.15], $
        xtitle='(u-r)', $
        datacolor=color_original, $
        /oprobability, probcolor=color_original, $
        thick=thickline

    cghistoplot, gz2_extra.petromag_u - gz2_extra.petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=ur_binsize, $
        mininput = ur_minval, $
        maxinput = ur_maxval, $
        datacolor=color_extra, $
        /oprobability, probcolor=color_extra, $
        thick=thickline

    cghistoplot, gz2_original[newinds].petromag_u - gz2_original[newinds].petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=ur_binsize, $
        datacolor=color_selected, $
        mininput = ur_minval, $
        maxinput = ur_maxval, $
        /oprobability, probcolor=color_selected, $
        thick=thickline

    cghistoplot, gz2_original.petromag_r, $
        /outline, $
        /freq, $
        bin=rmag_binsize, $
        charsize=cs, $
        ytickformat='(f4.2)', $
        ytitle= ' ', $
        xtitle='m!Ir!N', $
        xr=[rmag_minval,rmag_maxval], $
        yr=[0,0.30], $
        mininput = rmag_minval, $
        maxinput = rmag_maxval, $
        datacolor=color_original, $
        /oprobability, probcolor=color_original, $
        thick=thickline

    cghistoplot, gz2_extra.petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=rmag_binsize, $
        datacolor=color_extra, $
        mininput = rmag_minval, $
        maxinput = rmag_maxval, $
        /oprobability, probcolor=color_extra, $
        thick=thickline

    cghistoplot, gz2_original[newinds].petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=rmag_binsize, $
        datacolor=color_selected, $
        mininput = rmag_minval, $
        maxinput = rmag_maxval, $
        /oprobability, probcolor=color_selected, $
        thick=thickline

    cghistoplot, gz2_original.petror50_r, $
        /outline, $
        /freq, $
        bin=r50_binsize, $
        charsize=cs, $
        ytickformat='(f4.2)', $
        ytitle= ' ', $
        xr=[r50_minval,r50_maxval], $
        yr=[0,0.25], $
        xtitle='r-band R!I50!N [arcsec]', $
        mininput = r50_minval, $
        maxinput = r50_maxval, $
        datacolor=color_original, $
        /oprobability, probcolor=color_original, $
        thick=thickline

    cghistoplot, gz2_extra.petror50_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=r50_binsize, $
        datacolor=color_extra, $
        mininput = r50_minval, $
        maxinput = r50_maxval, $
        /oprobability, probcolor=color_extra, $
        thick=thickline

    cghistoplot, gz2_original[newinds].petror50_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=r50_binsize, $
        datacolor=color_selected, $
        mininput = r50_minval, $
        maxinput = r50_maxval, $
        /oprobability, probcolor=color_selected, $
        thick=thickline

    al_legend, ['Original', 'Extra','Original','subsample'], $
        position=[6.0,0.175], $
        linestyle=0, $
        thick = thickline+2, $
        charsize = 1.6, $
        linsize = 0.15, $
        color=[color_original,color_extra,color_selected,'white']

    ; Add plots for three example tasks - 01, 02, 03

    taskarr = [0,3,5]
    ymaxarr = [0.10,0.65,0.50]
    xtitlearr = ['smooth vote fraction', 'edge-on vote fraction', 'bar vote fraction']
    ytitlearr = ['Relative frequency', ' ', ' ']

    for k = 0,2 do begin
        if n_elements(votemin) eq 0 then votemin=0

        original_selected = gz2_original[newinds]

        ind_original = where(gz2_original.(tc_all[k]) ge votemin,n_original)
        ind_extra = where(gz2_extra.(tc_all[k]) ge votemin,n_extra)
        ind_original_selected = where(original_selected.(tc_all[k]) ge votemin,n_selected)

        cghistoplot, gz2_extra[ind_extra].(wfind[taskarr[k]]), $
;            title=strmid(tagnames[wfind[k]],0,strlen(tagnames[wfind[k]])-18), $
            /outline, $
            /freq, $
            bin=0.05, $
            charsize=cs, $
            histdata = extra_task_histdata, $
            ytickformat='(f4.2)', $
            yr=[0,ymaxarr[k]], $
            xtitle=xtitlearr[k], $
            ytitle=ytitlearr[k], $
            datacolor=color_extra, $
            /oprobability, probcolor=color_extra, $
            thick=thickline

        cghistoplot, original_selected[ind_original_selected].(wfind[taskarr[k]]), $
            /oplot, $
            /outline, $
            /freq, $
            bin=0.05, $
            datacolor=color_selected, $
            histdata = selected_task_histdata, $
            /oprobability, probcolor=color_selected, $
            thick=thickline

        ;cgplots, replicate(mean(gz2_extra[ind_extra].(wfind[k])),2), !y.crange, $
        ;    color=color_extra, $
        ;    thick = 1, $
        ;    linestyle=2

        ;cgplots, replicate(mean(original_selected[ind_original_selected].(wfind[k])),2), !y.crange, $
        ;    color=color_selected, $
        ;    thick = 1, $
        ;    linestyle=2

    endfor

    al_legend, ['Extra','Original','subsample'], $
        position=[0.25,0.35], $
        linestyle=0, $
        thick = thickline+2, $
        charsize = 1.8, $
        linsize = 0.15, $
        color=[color_extra,color_selected,'white']

if keyword_set(ps) then ps_end

; Try the new sample from the original galaxies and compare it to the extra

print,''

!p.multi=[0,8,5]

for i=0,nr-1 do begin

    if n_elements(votemin) eq 0 then votemin=0
    original_selected = gz2_original[newinds]

    ind_original = where(gz2_original.(tc_all[i]) ge votemin,n_original)
    ind_extra = where(gz2_extra.(tc_all[i]) ge votemin,n_extra)
    ind_original_selected = where(original_selected.(tc_all[i]) ge votemin,n_selected)

    cghistoplot, gz2_original[ind_original].(wfind[i]), $
        title=strmid(tagnames[wfind[i]],0,strlen(tagnames[wfind[i]])-18), $
        /outline, $
        /freq, $
        bin=0.05, $
        charsize=1, $
        histdata = original_task_histdata, $
        ytickformat='(f3.1)', $
        datacolor=color_original, $
        /oprobability, probcolor=color_original, $
        thick=2

    cghistoplot, gz2_extra[ind_extra].(wfind[i]), $
        /oplot, $
        /outline, $
        /freq, $
        bin=0.05, $
        datacolor=color_extra, $
        histdata = extra_task_histdata, $
        /oprobability, probcolor=color_extra, $
        thick=2

    cghistoplot, original_selected[ind_original_selected].(wfind[i]), $
        /oplot, $
        /outline, $
        /freq, $
        bin=0.05, $
        datacolor=color_selected, $
        histdata = selected_task_histdata, $
        /oprobability, probcolor=color_selected, $
        thick=2

    ;cgplots, replicate(mean(gz2_original[ind_original].(wfind[i])),2), !y.crange, $
    ;    color=color_original, $
    ;    thick = 1, $
    ;    linestyle=2

    cgplots, replicate(mean(gz2_extra[ind_extra].(wfind[i])),2), !y.crange, $
        color=color_extra, $
        thick = 1, $
        linestyle=2

    cgplots, replicate(mean(original_selected[ind_original_selected].(wfind[i])),2), !y.crange, $
        color=color_selected, $
        thick = 1, $
        linestyle=2

;kstwo, gz2_original[ind_original].(wfind[i]), gz2_extra[ind_extra].(wfind[i]), d, p 
;kstwo, original_selected[ind_original_selected].(wfind[i]), gz2_extra[ind_extra].(wfind[i]), d, p 
kstwo, selected_task_histdata, extra_task_histdata, d, p 
print,string(strmid(tagnames[wfind[i]],0,strlen(tagnames[wfind[i]])-18),format='(a-45)'),$
    string(d,format='(f5.3)'), string(p,format='(f7.4)'), $
    string(abs(mean(original_selected[ind_original_selected].(wfind[i])) - mean(gz2_extra[ind_extra].(wfind[i]))),format='(f5.2)'), $
    string(n_original, format='(i8)'),string(n_extra, format='(i8)'),string(float(n_original)/n_extra,format='(f7.3)')

endfor
print,''

    cgplot, indgen(10), /nodata, xr=[0,1], yr=[0,1]
    al_legend, ['Original','Extra','Original (sub-selected'], $
        position=[0.1,0.7], $
        linestyle=0, $
        thick = 2, $
        charsize = 1, $
        linsize = 0.2, $
        color=[color_original,color_extra,color_selected]



if keyword_set(stop) then stop

end
