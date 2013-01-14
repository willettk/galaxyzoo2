file='~/Astronomy/proposals/redspirals_chandra/masters2010/mnr_16503_sm_TableA1.txt'
file='~/Astronomy/proposals/redspirals_chandra/masters2010/mnr_16503_sm_TableA2.txt'

readcol, file, $
	objid, ra, dec, redshift, $
	absmag_r, absmag_r_err, color_gr, color_gr_err, redness, ellipticity, fracdeV, $
	format='a,f,f,f,f,f,f,f,f,f', $
	comment='%', $
	/silent, $
	skip=4
	
openw, lun, '~/Astronomy/proposals/redspirals_vla/first_bluespirals_upload.txt',/get_lun

for i=0,n_elements(ra) -1 do begin

	ratemp = sixty(ra[i]/15.)
	rastring = [string(ratemp[0:1],format='(i02)'),string(ratemp[2],format='(f04.1)')]
	dectemp = sixty(dec[i])
	if dectemp[0] eq 0 and dectemp[1] lt 0 then decsign = '-' else decsign = '+'
	if dectemp[0] gt 0 then decsign = '+'
	if dectemp[0] lt 0 then decsign = '-'
	decstring = [decsign+string(abs(dectemp[0]),format='(i02)'),string(abs(dectemp[1]),format='(i02)'),string(dectemp[2],format='(f04.1)')]
	printf,lun,[rastring,decstring]

endfor

close, lun

end
