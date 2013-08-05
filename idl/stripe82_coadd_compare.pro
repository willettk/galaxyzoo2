
;+
; NAME:
;       
;	STRIPE82_COADD_COMPARE
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

pro stripe82_coadd_compare, stop=stop, ps=ps, lowz=lowz, label=label

tstart = systime(1)

set_plot,'x'
device, retain=2

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
figsdir = '~/Astronomy/Research/GalaxyZoo/datapaper/figures/'

if n_elements(lowz) gt 0 then begin
	coadd_all = mrdfits(gz2dir+'fits/coadd1_coadd2_sample.fits',1,/silent)
	coadd = coadd_all[where(coadd_all.redshift gt 0. and coadd_all.redshift le lowz)]
endif else coadd = mrdfits(gz2dir+'fits/coadd1_coadd2.fits',1,/silent)

tagnames=tag_names(coadd)
nt = n_elements(tagnames)

wfind1_temp=intarr(nt) 
wcind1_temp=intarr(nt) 
wfind2_temp=intarr(nt) 
wcind2_temp=intarr(nt) 

for i=0,nt-1 do begin
	if (strlen(tagnames[i]) ge 20) and (strmid(tagnames[i],strlen(tagnames[i])-19,19) eq 'WEIGHTED_FRACTION_1') then wfind1_temp[i]=1 
	if (strlen(tagnames[i]) ge 15) and (strmid(tagnames[i],strlen(tagnames[i])-14,14) eq 'TOTAL_WEIGHT_1') then wcind1_temp[i]=1 
	if (strlen(tagnames[i]) ge 20) and (strmid(tagnames[i],strlen(tagnames[i])-19,19) eq 'WEIGHTED_FRACTION_2') then wfind2_temp[i]=1 
	if (strlen(tagnames[i]) ge 15) and (strmid(tagnames[i],strlen(tagnames[i])-14,14) eq 'TOTAL_WEIGHT_2') then wcind2_temp[i]=1 
endfor

wfind1 = where(wfind1_temp, nw1)	; Weighted fraction
wcind1 = where(wcind1_temp, nc1)	; Total weight (count)
wfind2 = where(wfind2_temp, nw2)	; Weighted fraction
wcind2 = where(wcind2_temp, nc2)	; Total weight (count)

; Do the weighted vote fractions for the various GZ2 tasks agree between the coadd1 and coadd2 samples?

count = 10

if keyword_set(ps) then begin

	plotname='stripe82_coadd_compare_all.ps'

	ps_start, filename=figsdir+plotname, /color, /quiet
	cs=1.2
	legendcs = 0.9
	th=3
	thickline=10
	thinline=1
endif else begin
	cs=2
	legendcs = cs
	th=1
	th=3
	thickline=3
	thinline=1
endelse

meanarr = dblarr(nw1)
medianarr = dblarr(nw1)
stddevarr = dblarr(nw1)
tasknamearr = strarr(nw1)

;!p.multi=[0,2,2]
print,''
for i=0,nw1-1 do begin
	
	; Find which column the total weight for each response is in

	taskno1 = strmid(tagnames[wfind1[i]],1,2)
	task_weight_col1 = wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq taskno1)]
	taskno2 = strmid(tagnames[wfind2[i]],1,2)
	task_weight_col2 = wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq taskno2)]

	; Display data where the task for that galaxy had at least N responses

	taskind = where(coadd.(task_weight_col1) gt count and coadd.(task_weight_col2) gt count, taskcount)

	; Display data where the answer for that galaxy had at least N responses
	; This is only useful if it's the plurality of the answer, really.

	answerind = where(coadd.(wfind1[i]-2) gt count and coadd.(wfind2[i]-2) gt count,answercount)

	if answercount gt 5 then begin
		goodind = answerind 
		tag = ' '
	endif else begin
		goodind = taskind
		tag = 'T'
		print,tagnames[wfind1[i]-2],answercount
	endelse

	goodind = taskind

	delta_coadd = coadd[goodind].(wfind1[i]) - coadd[goodind].(wfind2[i])
	cghistoplot, delta_coadd, $
		bin = 0.05, $
		xr=[-1,1], /xstyle, $
		/outline,$
		backcolorname='white', axiscolor='black', datacolor='red', $
		charsize=1.2, $
		xtitle='coadd1 - coadd2', title='GZ2 Task '+string(taskno1,format='(i02)')+', N='+string(n_elements(goodind),format='(i6)') + ' '+tag

	cgplots, mean(delta_coadd), !y.crange, color='black', linestyle=1
	cgplots, median(delta_coadd), !y.crange, color='black', linestyle=2

	;cgplot, coadd[goodind].(wfind1[i]), coadd[goodind].(wfind2[i]), $
	;	xr=[0,1], /xstyle, $
	;	yr=[0,1], /ystyle, $
	;	psym = 3, $
	;	background='white', axiscolor='black', $
	;	charsize=cs, $
	;	xtitle='coadd1', $
	;	ytitle='coadd2', $
	;	title='GZ2 Task '+string(taskno1,format='(i02)')+', N='+string(n_elements(goodind),format='(i6)') + ' '+tag

	kstwo, coadd[goodind].(wfind1[i]), coadd[goodind].(wfind2[i]), d, p
	;print,string(tagnames[wfind1[i]],format='(a-70)'), string(mean(delta_coadd),format='(f10.3)'), string(median(delta_coadd),format='(f10.3)'), string(skewness(delta_coadd),format='(f10.3)'), string(probgauss(p),format='(f7.2)')

	; Weighted fractions of each response for both main and Stripe 82 samples. Note that a magnitude cut must be applied to Stripe 82, since it goes to r < 17.7

	;main_wf = gz2[where(gz2.(task_weight_col1[0]) gt count and (strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra'),maincount)].(wfind[i])
	;s82_normal_wf = gz2[where(gz2.(task_weight_col1[0]) gt count and strtrim(gz2.sample,2) eq 'stripe82' and gz2.petromag_r lt 17.0,s82count)].(wfind[i])

	; Print results to screen

	;print,string(strmid(tagnames[wfind[i]],0,strlen(tagnames[wfind[i]])-18),format='(a-70)'),' & ',string(mean(main_wf),format='(f10.3)'), ' & ',string(mean(s82_normal_wf),format='(f10.3)'),' & ',string(maincount,format='(i10)'),' & ',string(s82count,format='(i10)'),' & ',string(mean(main_wf)-mean(s82_normal_wf),format='(f11.3)'), ' & ',string(mean(main_wf)/mean(s82_normal_wf),format='(f11.3)') + ' \\'

	meanarr[i] = mean(delta_coadd)
	medianarr[i] = median(delta_coadd)
	stddevarr[i] = stddev(delta_coadd)
	tasknamearr[i] = string(taskno1,format='(i02)')

	;if i eq 27 then stop

endfor
print,''

if keyword_set(ps) then ps_end

; Make a truncated plot for the paper - display histogram for 2 or 3 major tasks?

; Alternate idea - compute the mean and standard deviation of the difference for each task, then plot all those on a single scatter plot

if keyword_set(ps) then begin

	plotname='stripe82_coadd_compare.eps'

	ps_start, filename=figsdir+plotname, /color, /quiet, xsize=10, ysize=5, /encap
	cs=1.2
endif else begin
	cs=2
	legendcs = cs
	th=1
	th=3
	thickline=3
	thinline=1
endelse

; Ellipticals

taskind = where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[0]],1,2))]) gt count and coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[0]],1,2))]) gt count, taskcount)
delta_coadd = coadd[taskind].(wfind1[0]) - coadd[taskind].(wfind2[0])

h = histogram(delta_coadd,binsize=0.05)
harr = fillarr(0.05,min(delta_coadd),max(delta_coadd))
cgplot, $
	harr, h/float(max(h)), $
	psym=10, $
	xr=[-1,1], /xstyle, $
	yr=[0,1], /ystyle, $
	position=[0.1,0.5,0.3,0.9], $
	/normal, $
	xticks=2, $
	xtickname=replicate(' ',3), $
	yticks=2, $
	;ytickname=replicate(' ',3), $
	ytitle=' ', $
	color='red', $
	charsize=1.0

cgtext, -0.9, 0.9, /data, 'Smooth', charsize=1.0
cgplots, replicate(median(delta_coadd),2), !y.crange, linestyle=2

; Spirals

taskind = where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[1]],1,2))]) gt count and coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[1]],1,2))]) gt count, taskcount)
delta_coadd = coadd[taskind].(wfind1[1]) - coadd[taskind].(wfind2[1])

h = histogram(delta_coadd,binsize=0.05)
harr = fillarr(0.05,min(delta_coadd),max(delta_coadd))
cgplot, $
	harr, h/float(max(h)), $
	psym=10, $
	xr=[-1,1], /xstyle, $
	yr=[0,1], /ystyle, $
	position=[0.3,0.5,0.5,0.9], $
	/noerase, $
	xticks=2, $
	xtickname=replicate(' ',3), $
	yticks=2, $
	ytickname=replicate(' ',3), $
	ytitle=' ', $
	color='red'

cgtext, 0.9, 0.9, /data, 'Features', charsize=1.0, align=1
cgtext, 0.9, 0.82, /data, 'or disk', charsize=1.0, align=1
cgplots, replicate(median(delta_coadd),2), !y.crange, linestyle=2

; Bar

taskind = where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[5]],1,2))]) gt count and coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[5]],1,2))]) gt count, taskcount)
delta_coadd = coadd[taskind].(wfind1[5]) - coadd[taskind].(wfind2[5])

h = histogram(delta_coadd,binsize=0.05)
harr = fillarr(0.05,min(delta_coadd),max(delta_coadd))
cgplot, $
	harr, h/float(max(h)), $
	psym=10, $
	xr=[-1,1], /xstyle, $
	yr=[0,1], /ystyle, $
	position=[0.1,0.1,0.3,0.5], $
	/noerase, $
	xticks=4, $
	xtickname=['-1.0','-0.5','0','0.5',' '], $
	yticks=2, $
	;ytickname=replicate(' ',3), $
	ytitle=' ', $
	color='red', $
	charsize=1.0

cgtext, -0.9, 0.9, /data, 'Bar', charsize=1.0
cgplots, replicate(median(delta_coadd),2), !y.crange, linestyle=2

; Noticeable bulge

taskind = where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[10]],1,2))]) gt count and coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[10]],1,2))]) gt count, taskcount)
delta_coadd = coadd[taskind].(wfind1[10]) - coadd[taskind].(wfind2[10])

h = histogram(delta_coadd,binsize=0.05)
harr = fillarr(0.05,min(delta_coadd),max(delta_coadd))
cgplot, $
	harr, h/float(max(h)), $
	psym=10, $
	xr=[-1,1], /xstyle, $
	yr=[0,1], /ystyle, $
	position=[0.3,0.1,0.5,0.5], $
	/noerase, $
	xticks=4, $
	xtickname=[' ','-0.5','0','0.5','1.0'], $
	yticks=2, $
	ytickname=replicate(' ',3), $
	ytitle=' ', $
	color='red', $
	charsize=1.0

cgtext, 0.9, 0.9, /data, 'Noticeable', charsize=1.0, align=1
cgtext, 0.9, 0.82, /data, 'bulge', charsize=1.0, align=1
cgplots, replicate(median(delta_coadd),2), !y.crange, linestyle=2

cgtext, 0.18, 0.02, '(coadd1 - coadd2)', charsize=1.3, /normal
cgtext, 0.04, 0.20, 'Histogram density (normalized)', charsize=1.3, /normal, orient=90

; Scatter plot
cgplot, meanarr, stddevarr, psym=16, $
	xtitle='Mean (coadd1 - coadd2)', $
	ytitle='Std. dev (coadd1 - coadd2)', $
	position=[0.6,0.1,0.95,0.9], $
	/noerase, $
	;xr=[-0.25,0.25], /xstyle, $
	;yr=[-0.25,0.25], /ystyle, $
	charsize=1.3

cgplots, !x.crange, [0,0], linestyle=1
cgplots, [0,0], !y.crange, linestyle=1

if keyword_set(label) then begin
	cgtext, meanarr[where(meanarr gt 0)], stddevarr[where(meanarr gt 0)], ' T'+tasknamearr[where(meanarr gt 0)], charsize=1.2, align=0
	cgtext, meanarr[where(meanarr lt 0)], stddevarr[where(meanarr lt 0)], 'T'+tasknamearr[where(meanarr lt 0)]+ ' ', charsize=1.2, align=1
endif

print,''
print,n_elements(where(abs(meanarr) gt 0.05))

if keyword_set(ps) then ps_end

; By how much is it different for just noticeable bulges?

taskind = where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[10]],1,2))]) gt count and coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[10]],1,2))]) gt count, taskcount)

print,'Answer 10'
val1 = coadd[taskind].(wfind1[10])
val2 = coadd[taskind].(wfind2[10])
print,mean(val1),stddev(val1)
print,mean(val2),stddev(val2)
print,'Coadd1 is '+string((mean(val2)/mean(val1) - 1d)*100,format='(f5.1)')+'% higher'
print,''

taskind = where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[11]],1,2))]) gt count and coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[11]],1,2))]) gt count, taskcount)

print,'Answer 11'
val1 = coadd[taskind].(wfind1[11])
val2 = coadd[taskind].(wfind2[11])
print,mean(val1),stddev(val1)
print,mean(val2),stddev(val2)
print,'Coadd2 is '+string((mean(val1)/mean(val2)-1d)*100,format='(f5.1)')+'% higher'

; Ellipticals

e_coadd1 = n_elements(where(coadd[where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[0]],1,2))]) ge count)].(wfind1[0]) ge 0.8))
e_coadd2 = n_elements(where(coadd[where(coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[0]],1,2))]) ge count)].(wfind2[0]) ge 0.8))
e_both = n_elements(where(coadd[where(coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[0]],1,2))]) ge count)].(wfind2[0]) ge 0.8 and coadd[where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[0]],1,2))]) ge count)].(wfind1[0]) ge 0.8))
e_either = n_elements(where(coadd[where(coadd.(wcind2[where(strmid(strtrim(tagnames[wcind2],2),1,2) eq strmid(tagnames[wfind2[0]],1,2))]) ge count)].(wfind2[0]) ge 0.8 or coadd[where(coadd.(wcind1[where(strmid(strtrim(tagnames[wcind1],2),1,2) eq strmid(tagnames[wfind1[0]],1,2))]) ge count)].(wfind1[0]) ge 0.8))
print,e_coadd1,double(e_coadd1)/n_elements(coadd)*100.,' ellipticals in coadd1'
print,e_coadd2,double(e_coadd2)/n_elements(coadd)*100.,' ellipticals in coadd2'
print,e_both,double(e_both)/n_elements(coadd)*100.,' ellipticals in both'
print,e_either,double(e_either)/n_elements(coadd)*100.,' ellipticals in either'

tend = systime(1)

print,''
print,'Time elapsed: ',string(tend - tstart,format='(f10.2)'),' sec'

if keyword_set(stop) then stop

end
