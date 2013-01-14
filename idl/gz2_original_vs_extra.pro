
;+
; NAME:
;       
;	GZ2_ORIGINAL_VS_EXTRA
;
; PURPOSE:
;
;	Compare the "original" and "extra" samples in GZ2. Look to see if there are any obvious time-dependent biases. 
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


pro gz2_original_vs_extra, votemin, ps=ps, stop=stop

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
;gz2_main            = gz2[where((strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra'),count_main)]

color_main = 'black'
color_original = 'magenta'
color_extra = 'blue'

ps_start, filename='~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/gz2_original_vs_extra.eps',$
    /encap,/color,xs=16,ys=8, /quiet
!p.multi=[0,8,5]

print,''

; All objects

for i=0,nr-1 do begin

    if n_elements(votemin) eq 0 then votemin=0
    ind_original = where(gz2_original.(tc_all[i]) ge votemin,n_original)
    ind_extra = where(gz2_extra.(tc_all[i]) ge votemin,n_extra)

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

    ;cgplots, replicate(mean(gz2_original[ind_original].(wfind[i])),2), !y.crange, $
    ;    color=color_original, $
    ;    thick = 1, $
    ;    linestyle=2

    ;cgplots, replicate(mean(gz2_extra[ind_extra].(wfind[i])),2), !y.crange, $
    ;    color=color_extra, $
    ;    thick = 1, $
    ;    linestyle=2

kstwo, gz2_original[ind_original].(wfind[i]), gz2_extra[ind_extra].(wfind[i]), d, p 
;print,string(strmid(tagnames[wfind[i]],0,strlen(tagnames[wfind[i]])-18),format='(a-45)'),$
;    string(d,format='(f5.3)'), string(p,format='(f7.3)'), $
;    string(n_original, format='(i8)'),string(n_extra, format='(i8)'),string(float(n_original)/n_extra,format='(f7.3)')

endfor
print,''

    cgplot, indgen(10), /nodata, xr=[0,1], yr=[0,1]
    al_legend, ['Original','Extra'], $
        position=[0.1,0.7], $
        linestyle=0, $
        thick = 2, $
        charsize = 1, $
        linsize = 0.2, $
        color=[color_original,color_extra]


ps_end

; K-S tests for the weighted fractions in all tasks

; Histograms of the weighted fractions for the vote responses for all tasks

; Does the derived bias change if I only use the "original" GZ1 data?

if keyword_set(stop) then stop

end
