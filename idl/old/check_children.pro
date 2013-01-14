
;+
; NAME:
;       
;	
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

pro check_children, votemin, sortvariable=sortvariable, ps=ps, stop=stop

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
        /encap,/color,xs=16,ys=8, /quiet
    cs = 1
endif else begin
    cs = 3
endelse
wset,0
!p.multi=[0,3,1]

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

    if n_elements(sortvariable) eq 0 then sortvariable = 'ur'

    if sortvariable eq 'ur' then begin
        original_histdata = ur_original_histdata
        extra_histdata = ur_extra_histdata
        minval = ur_minval
        maxval = ur_maxval
        binsize = ur_binsize
        old_data_val = gz2_original.petromag_u - gz2_original.petromag_r
        sortstr = '(u-r) color'
    endif

    if sortvariable eq 'r50' then begin
        original_histdata = r50_original_histdata
        extra_histdata = r50_extra_histdata
        minval = r50_minval
        maxval = r50_maxval
        binsize = r50_binsize
        old_data_val = gz2_original.petror50_r
        sortstr = 'R50'
    endif

; Draw randomly from the PDF determined from the extra sample

np = n_elements(extra_histdata)
extra_histdata_norm = float(extra_histdata) / max(extra_histdata)
xrange = [minval,maxval]  ; Range of PDF
yrange = [0.,1.]  ; Range of random numbers and maximum prob

nnum = n_elements(gz2_extra)

;;; Set up some useful arrays
bad = lindgen(nnum)
pxi = fltarr(nnum)
xi = fltarr(nnum)
yi = fltarr(nnum)
nbad = nnum
cool = 0b                       ;loop control variable
ct = 0

xbinned = indgen(np) * (maxval - minval) / np + minval
xstep = xbinned[1] - xbinned[0]

while not cool do begin
    ;;; Get set of x values for values not yet settled
    xi[bad] = drawrand_uniform(nbad, xrange)
    ;;; Get set of y values for values not yet settled
    yi[bad] = drawrand_uniform(nbad, yrange)
    ;;; Get set of comparison values for indices not yet settled

    for i=0,nbad-1 do pxi[bad[i]] = extra_histdata_norm[closeto(xbinned,xi[bad[i]])]

;    pxi[bad] = call_function(userfunc, xi[bad])
    ;;; Is pxi > yi?  If not, mark the elements which require
    ;;; changing and rinse, repeat. Dry once there are no "bad"
    ;;; values 
    bad = where(yi gt pxi, nbad, complement = good)
    ngood = n_elements(good)
    ct += 1.
    if nbad eq 0 then cool = 1
endwhile

    ; Not at all sure this is valid, but it makes the u-r color histograms match much better. 
    xi += xstep

h_new = histogram(xi,binsize=binsize,min=minval,max=maxval)
h_old = histogram(old_data_val,min=minval,max=maxval,binsize=binsize, reverse_indices=ri)
newinds = lonarr(nnum)
sind = 0

for i=0,n_elements(h_new)-1 do begin

    ; For each bin in the self-index vector i, find the new colors that match 
    ; and randomly select an appropriate sub-sample from amongst them. 

    nb = h_new[i]
    if nb gt 0 then begin
        inds = ri[ri[i]:ri[i+1]-1]
        indind = fix(drawrand_uniform(nb,[0,nb]))
        newinds[sind:sind+nb-1] = inds[indind]

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; Trying this method of matching both color and size did not work 
        ; - only yielded 149 unique galaxies out of 22000+ expected. 

        ;; Now match the randomly-selected (u-r) color matches against the R50 sizes in that color bin

        ; inds = ri[ri[i]:ri[i+1]-1]
        ;;;; Set up some useful arrays
        ;bad_r50 = lindgen(nb)
        ;pxi_r50 = fltarr(nb)
        ;xi_r50 = fltarr(nb)
        ;yi_r50 = fltarr(nb)
        ;nbad_r50 = nb
        ;cool = 0b                       ;loop control variable
        ;ct = 0L
        ;xrange_r50=[0,n_elements(inds)-1]

        ;r50_original_histdata = histogram(gz2_original[inds].petror50_r,binsize=r50_binsize,min=r50_minval,max=r50_maxval)
        ;r50_original_histdata_norm = float(r50_original_histdata) / max(r50_original_histdata)
        ;nr50 = n_elements(r50_original_histdata)

        ;r50_bins = indgen(nr50) * float(r50_maxval - r50_minval) / nr50 + r50_minval
        ;r50_step = r50_bins[1] - r50_bins[0]

        ;while not cool do begin
        ;    xi_r50[bad_r50] = drawrand_uniform(nbad_r50, xrange_r50)
        ;    yi_r50[bad_r50] = drawrand_uniform(nbad_r50, yrange)
        ;    pxi_r50[bad_r50] = pdf_binned_norm(r50_original_histdata_norm,r50_bins,xi_r50[bad_r50])
        ;    bad_r50 = where(yi_r50 gt pxi_r50, nbad_r50, complement = good_r50)
        ;    ngood_r50 = n_elements(good_r50)
        ;    ct += 1
        ;    if nbad_r50 eq 0 then cool = 1
        ;endwhile

        ; Not at all sure this is valid, but it makes the u-r color histograms match much better. 
        ; xi_r50 += r50_step

        ; newinds[sind:(sind + nb) - 1] = inds[xi_r50]

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    endif 
    sind += nb

endfor

; Last thing to try. Create a 2D histogram of the original sample galaxies in (u-r) and R50. For every galaxy
; in the extra sample (just use a FOR loop - only 20K galaxies), locate the closest match on the histogram. 

data0 = gz2_original.petromag_u - gz2_original.petromag_r
data1 = gz2_original.petror50_r
vo = transpose([[data0],[data1]])
data0 = gz2_extra.petromag_u - gz2_extra.petromag_r
data1 = gz2_extra.petror50_r
ve = transpose([[data0],[data1]])
binsizes = [ur_binsize,r50_binsize]
mins = [ur_minval,r50_minval]
maxs = [ur_maxval,r50_maxval]

hist_original = hist_nd(vo,binsizes,min=mins,max=maxs,reverse_indices=ri_original)
hist_extra    = hist_nd(ve,binsizes,min=mins,max=maxs,reverse_indices=ri_extra)

oe = where(hist_original ge hist_extra and hist_extra gt 0,oecount)
sind=0
newinds_histnd = lonarr(total(hist_extra[oe]))

for i=0,oecount-1 do begin

    nnew = hist_extra[oe[i]]
    inds=ri_original[ri_original[oe[i]]:ri_original[oe[i]+1]-1]
    indind = fix(drawrand_uniform(nnew,[0,nnew]))
    newinds_histnd[sind:sind+nnew-1] = inds[indind] 
    sind += nnew

endfor

;ind=[i+nx*(j+ny*k)]
;ri[ri[ind]:ri[ind+1]-1]

out = unique(newinds,ucount,/sort)
print,'Unique galaxies in selected sub-sample from "original":',n_elements(out)

out = unique(newinds_histnd,ucount,/sort)
print,'Unique galaxies in selected sub-sample from "original":',n_elements(out)

newinds = newinds_histnd

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
        xtitle='petromag_u - petromag_r', $
        datacolor=color_original, $
        /oprobability, probcolor=color_original, $
        thick=2

    cghistoplot, gz2_extra.petromag_u - gz2_extra.petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=ur_binsize, $
        mininput = ur_minval, $
        maxinput = ur_maxval, $
        datacolor=color_extra, $
        /oprobability, probcolor=color_extra, $
        thick=2

    cghistoplot, gz2_original[newinds].petromag_u - gz2_original[newinds].petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=ur_binsize, $
        datacolor=color_selected, $
        mininput = ur_minval, $
        maxinput = ur_maxval, $
        /oprobability, probcolor=color_selected, $
        thick=2

    cghistoplot, gz2_original.petromag_r, $
        /outline, $
        /freq, $
        bin=rmag_binsize, $
        charsize=cs, $
        ytickformat='(f4.2)', $
        xtitle='petromag_r', $
        xr=[rmag_minval,rmag_maxval], $
        yr=[0,0.30], $
        mininput = rmag_minval, $
        maxinput = rmag_maxval, $
        datacolor=color_original, $
        /oprobability, probcolor=color_original, $
        thick=2

    cghistoplot, gz2_extra.petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=rmag_binsize, $
        datacolor=color_extra, $
        mininput = rmag_minval, $
        maxinput = rmag_maxval, $
        /oprobability, probcolor=color_extra, $
        thick=2

    cghistoplot, gz2_original[newinds].petromag_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=rmag_binsize, $
        datacolor=color_selected, $
        mininput = rmag_minval, $
        maxinput = rmag_maxval, $
        /oprobability, probcolor=color_selected, $
        thick=2

    cghistoplot, gz2_original.petror50_r, $
        /outline, $
        /freq, $
        bin=r50_binsize, $
        charsize=cs, $
        ytickformat='(f4.2)', $
        xr=[r50_minval,r50_maxval], $
        yr=[0,0.25], $
        xtitle='petroR50_r', $
        mininput = r50_minval, $
        maxinput = r50_maxval, $
        datacolor=color_original, $
        /oprobability, probcolor=color_original, $
        thick=2

    cghistoplot, gz2_extra.petror50_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=r50_binsize, $
        datacolor=color_extra, $
        mininput = r50_minval, $
        maxinput = r50_maxval, $
        /oprobability, probcolor=color_extra, $
        thick=2

    cghistoplot, gz2_original[newinds].petror50_r, $
        /oplot, $
        /outline, $
        /freq, $
        bin=r50_binsize, $
        datacolor=color_selected, $
        mininput = r50_minval, $
        maxinput = r50_maxval, $
        /oprobability, probcolor=color_selected, $
        thick=2

if keyword_set(ps) then ps_end

; Try the new sample from the original galaxies and compare it to the extra

print,''

wset,4
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
        /oprobability, probcolor=color_extra, $
        thick=2

    cghistoplot, original_selected[ind_original_selected].(wfind[i]), $
        /oplot, $
        /outline, $
        /freq, $
        bin=0.05, $
        datacolor=color_selected, $
        /oprobability, probcolor=color_selected, $
        thick=2

    ;cgplots, replicate(mean(gz2_original[ind_original].(wfind[i])),2), !y.crange, $
    ;    color=color_original, $
    ;    thick = 1, $
    ;    linestyle=2

    ;cgplots, replicate(mean(gz2_extra[ind_extra].(wfind[i])),2), !y.crange, $
    ;    color=color_extra, $
    ;    thick = 1, $
    ;    linestyle=2

;kstwo, gz2_original[ind_original].(wfind[i]), gz2_extra[ind_extra].(wfind[i]), d, p 
kstwo, original_selected[ind_original_selected].(wfind[i]), gz2_extra[ind_extra].(wfind[i]), d, p 
;print,string(strmid(tagnames[wfind[i]],0,strlen(tagnames[wfind[i]])-18),format='(a-45)'),$
;    string(d,format='(f5.3)'), string(p,format='(f7.3)'), $
;    string(n_original, format='(i8)'),string(n_extra, format='(i8)'),string(float(n_original)/n_extra,format='(f7.3)')

endfor
print,''

    cgplot, indgen(10), /nodata, xr=[0,1], yr=[0,1]
    al_legend, ['Original','Extra','Original '+sortstr+'-selected'], $
        position=[0.1,0.7], $
        linestyle=0, $
        thick = 2, $
        charsize = 1, $
        linsize = 0.2, $
        color=[color_original,color_extra,color_selected]


if keyword_set(stop) then stop

end
