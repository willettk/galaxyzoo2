
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

tstart = systime(1)

gz2dir = '~/Astronomy/Research/GalaxyZoo/'

;gz2 = mrdfits(gz2dir+'fits/gz2_table_sample_match_coadd2.fits',1,/silent)
;restore, '~/Astronomy/Research/GalaxyZoo/sav/gz2_sample_table_coadd2.sav'

tagnames=tag_names(gz2)
nt = n_elements(tagnames)

temparr=intarr(nt) 
temparr2=intarr(nt) 
for i=0,nt-1 do begin
	if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 
	if (strlen(tagnames[i]) ge 13) and (strmid(tagnames[i],strlen(tagnames[i])-12,12) eq 'TOTAL_WEIGHT') then temparr2[i]=1 
endfor

wfind = where(temparr,  nw)
wcind = where(temparr2, nc)

; Do the weighted vote fractions for the various GZ2 tasks agree between the original + extra samples and the Stripe 82 normal depth sample?

count = 10

print,''
for i=0,nw-1 do begin
	
	; Find which column the total weight for each response is in

	taskno = strmid(tagnames[wfind[i]],1,2)
	task_weight_col = wcind[where(strmid(strtrim(tagnames[wcind],2),1,2) eq taskno)]

	; Weighted fractions of each response for both main and Stripe 82 samples. Note that a magnitude cut must be applied to Stripe 82, since it goes to r < 17.7

	main_wf = gz2[where(gz2.(task_weight_col) gt count and (strtrim(gz2.sample,2) eq 'original' or strtrim(gz2.sample,2) eq 'extra'),maincount)].(wfind[i])
	s82_normal_wf = gz2[where(gz2.(task_weight_col) gt count and strtrim(gz2.sample,2) eq 'stripe82' and gz2.petromag_r lt 17.0,s82count)].(wfind[i])

	; Print results to screen

	print,string(strmid(tagnames[wfind[i]],0,strlen(tagnames[wfind[i]])-18),format='(a-70)'),' & ',string(mean(main_wf),format='(f10.3)'), ' & ',string(mean(s82_normal_wf),format='(f10.3)'),' & ',string(maincount,format='(i10)'),' & ',string(s82count,format='(i10)'),' & ',string(mean(main_wf)-mean(s82_normal_wf),format='(f11.3)'), ' & ',string(mean(main_wf)/mean(s82_normal_wf),format='(f11.3)') + ' \\'

endfor
print,''

tend = systime(1)
print,string(round(tend - tstart))+' sec'
print,''

stop

end
