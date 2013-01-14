
;+
; NAME:
;       
;	GZ2_COADD_KS
;
; PURPOSE:
;
;	Run K-S tests on weighted fraction distributions for all the tasks in coadd1 and coadd2 data
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
;       Written by K. Willett                Aug 12
;-

gz2dir = '~/Astronomy/Research/GalaxyZoo/'
fitsdir = '~/Astronomy/Research/GalaxyZoo/fits/'
figsdir = '~/Astronomy/Research/GalaxyZoo/gz2dropbox/figures/'

gz2_coadd1 = mrdfits(fitsdir+'gz2_table_sample_match_coadd1.fits',1,/silent)
gz2_coadd2 = mrdfits(fitsdir+'gz2_table_sample_match_coadd2.fits',1,/silent)

tagnames=tag_names(gz2_coadd1)
nt = n_elements(tagnames)

temparr=intarr(nt) 
temparr2=intarr(nt) 
for i=0,nt-1 do begin
	if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 
	if (strlen(tagnames[i]) ge 13) and (strmid(tagnames[i],strlen(tagnames[i])-12,12) eq 'TOTAL_WEIGHT') then temparr2[i]=1 
endfor

wfind = where(temparr,  nw)
wcind = where(temparr2, nc)

gz2main = gz2_coadd2[where((strtrim(gz2_coadd2.sample,2) eq 'original' or strtrim(gz2_coadd2.sample,2) eq 'extra'),maincount)]
normal = gz2_coadd2[where(strtrim(gz2_coadd2.sample,2) eq 'stripe82'         and (gz2_coadd2.petromag_r - gz2_coadd2.extinction_r) le 17.0 and gz2_coadd2.petror90_r ge 3.,s82n_count)]
coadd1 = gz2_coadd1[where(strtrim(gz2_coadd1.sample,2) eq 'stripe82_coadd_1' and (gz2_coadd1.petromag_r - gz2_coadd1.extinction_r) le 17.0 and gz2_coadd1.petror90_r ge 3.,s82c1_count)]
coadd2 = gz2_coadd2[where(strtrim(gz2_coadd2.sample,2) eq 'stripe82_coadd_2' and (gz2_coadd2.petromag_r - gz2_coadd2.extinction_r) le 17.0 and gz2_coadd2.petror90_r ge 3.,s82c2_count)]

count = 10

!p.multi=[0,6,6]
for i=0,nw-2 do begin & task_weight_col = wcind[where(strmid(strtrim(tagnames[wcind],2),1,2) eq strmid(tagnames[wfind[i]],1,2))] & kstwo, coadd1[where(coadd1.(task_weight_col) gt count,c1)].(wfind[i]), coadd2[where(coadd2.(task_weight_col) gt count,c2)].(wfind[i]), d, p & print,i,' ',string(tagnames[wfind[i]],format='(a-70)'),d,p,probgauss(p) & cghistoplot, /freq, yr=[0,1],xr=[-0.1,1.1],/xstyle, coadd1[where(coadd1.(wcind[where(strmid(strtrim(tagnames[wcind],2),1,2) eq strmid(tagnames[wfind[i]],1,2))]) gt count)].(wfind[i]),/outline,datacolor='red',binsize=0.1,title=greek('sigma')+'='+string(probgauss(p),format='(f3.1)')+'  N!I1!N='+strtrim(c1,2)+', N!I2!N='+strtrim(c2,2),charsize=2.5 & cghistoplot,/freq,yr=[0,1],/oplot,coadd2[where(coadd2.(wcind[where(strmid(strtrim(tagnames[wcind],2),1,2) eq strmid(tagnames[wfind[i]],1,2))]) gt count)].(wfind[i]),/outline,datacolor='blue',binsize=0.1 & endfor

for i=0,0 do begin & task_weight_col = wcind[where(strmid(strtrim(tagnames[wcind],2),1,2) eq strmid(tagnames[wfind[i]],1,2))] & kstwo, coadd1[where(coadd1.(task_weight_col) gt count,c1)].(wfind[i]), coadd2[where(coadd2.(task_weight_col) gt count,c2)].(wfind[i]), d, p & print,i,' ',string(tagnames[wfind[i]],format='(a-70)'),d,p,probgauss(p) & cghistoplot, /freq,/oprob, probcolor='red', yr=[0,1],xr=[-0.1,1.1],/xstyle, coadd1[where(coadd1.(wcind[where(strmid(strtrim(tagnames[wcind],2),1,2) eq strmid(tagnames[wfind[i]],1,2))]) gt count)].(wfind[i]),/outline,datacolor='red',binsize=0.1,title=greek('sigma')+'='+string(probgauss(p),format='(f3.1)')+'  N!I1!N='+strtrim(c1,2)+', N!I2!N='+strtrim(c2,2),charsize=2.5 & cghistoplot,/freq,/oprob,probcolor='blue',/oplot,coadd2[where(coadd2.(wcind[where(strmid(strtrim(tagnames[wcind],2),1,2) eq strmid(tagnames[wfind[i]],1,2))]) gt count)].(wfind[i]),/outline,datacolor='blue',binsize=0.1 & endfor
	

stop

end
