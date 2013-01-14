;+
; NAME:
;       
;	VCREDSHIFT_OLD
;
; PURPOSE:
;
;	Plot spectra and determine redshifts for objects in the Violin Clef system
;
;	This routine uses data from my FIRST reduction, without spikes and better S/N. 
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
;       Written by K. Willett                Oct 11
;-

; Load wavelength-calibrated data from IRAF

mmtdir = '~/Astronomy/Research/GalaxyZoo/violinclef/mmt_data/'

file_exp1 = mmtdir+'dispcor.0064.ms.fits'
file_exp2 = mmtdir+'dispcor.0066.ms.fits'

fits_exp1 = readfits(file_exp1, hdr1,/silent)
fits_exp2 = readfits(file_exp2, hdr2,/silent)

; North to south, the four galaxies in the system are labeled A, B, C, D.

; A - northernmost galaxy in the system
; B - smaller galaxy to the south of A; possibly a satellite.
; C - center galaxy in the system
; D - southernmost galaxy in the system

smoothlevel = 5

spec_galA = smooth(fits_exp1[*,2,0],smoothlevel)
spec_galB = smooth(fits_exp1[*,1,0],smoothlevel)
spec_galC = smooth(fits_exp2[*,0,0],smoothlevel)
spec_galD = smooth(fits_exp1[*,0,0],smoothlevel)

;wave_galA = indgen(n_elements(spec_galA))*sxpar(hdr_galA,'CD1_1') + sxpar(hdr_galA,'CRVAL1')
;wave_galB = indgen(n_elements(spec_galB))*sxpar(hdr_galB,'CD1_1') + sxpar(hdr_galB,'CRVAL1')
wave_galC = indgen(n_elements(spec_galC))*sxpar(hdr_galC,'CD1_1') + sxpar(hdr_galC,'CRVAL1')
;wave_galD = indgen(n_elements(spec_galD))*sxpar(hdr_galD,'CD1_1') + sxpar(hdr_galD,'CRVAL1')

wat2_001 = strsplit(sxpar(hdr1,'WAT2_001'),/ex)
wat2_002 = strsplit(sxpar(hdr1,'WAT2_002'),/ex)
wat2_003 = strsplit(sxpar(hdr1,'WAT2_003'),/ex)
wat2_004 = strsplit(sxpar(hdr1,'WAT2_004'),/ex)

crval1_galD = wat2_001[6]
cd1_1_galD = wat2_001[7]
crval1_galB = wat2_002[8]
cd1_1_galB = wat2_002[9]
crval1_galA = wat2_003[9]
cd1_1_galA = wat2_003[10]

wave_galD = indgen(n_elements(spec_galD))*cd1_1_galD + crval1_galD
wave_galB = indgen(n_elements(spec_galB))*cd1_1_galB + crval1_galB
wave_galA = indgen(n_elements(spec_galA))*cd1_1_galA + crval1_galA

; Posited redshifts for the sample

z_galA = 0.094
z_galB = 0.095
z_galC = 0.094
z_galD = 0.094

; Plot spectra

ps_start, file='~/Astronomy/Research/GalaxyZoo/violinclef/mmt_data/vcredshift_old.ps', /color, /landscape,/quiet

!p.multi=[0,1,4]

angstrom = '!6!sA!r!u!9 %!6!n!3'

cleanA = where(abs(spec_galA) lt 150.)
cleanB = where(abs(spec_galB) lt 150.)
cleanC = where(abs(spec_galC) lt 150.)
cleanD = where(abs(spec_galD) lt 150.)

cs = 1.2

cgplot, wave_galA[cleanA], spec_galA[cleanA], psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,200], title='Violin Clef galaxy A'
vclines, z_galA

cgplot, wave_galB[cleanB], spec_galB[cleanB], psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,200], title='Violin Clef galaxy B'
vclines, z_galB

cgplot, wave_galC, spec_galC, psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,800], title='Violin Clef galaxy C'
vclines, z_galC

cgplot, wave_galD[cleanD], spec_galD[cleanD], psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,200], title='Violin Clef galaxy D'
vclines, z_galD

!p.multi=[0,1,1]

ps_end

end
