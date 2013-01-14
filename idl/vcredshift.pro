;+
; NAME:
;       
;	VCREDSHIFT
;
; PURPOSE:
;
;	Plot spectra and determine redshifts for objects in the Violin Clef system
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

file_galA = mmtdir+'disp.exp1.galA.ms.fits'
file_galB = mmtdir+'disp.exp1.galB.ms.fits'
file_galD = mmtdir+'disp.exp1.galD.ms.fits'
file_galC = mmtdir+'disp.exp2.galC.ms.fits'

fits_galA = readfits(file_galA, hdr_galA, /silent)
fits_galB = readfits(file_galB, hdr_galB, /silent)
fits_galC = readfits(file_galC, hdr_galC, /silent)
fits_galD = readfits(file_galD, hdr_galD, /silent)

; North to south, the four galaxies in the system are labeled A, B, C, D.

; A - northernmost galaxy in the system
; B - smaller galaxy to the south of A; possibly a satellite.
; C - center galaxy in the system
; D - southernmost galaxy in the system

smoothlevel = 5

spec_galA = smooth(fits_galA[*,0,0],smoothlevel)
spec_galB = smooth(fits_galB[*,0,0],smoothlevel)
spec_galC = smooth(fits_galC[*,0,0],smoothlevel)
spec_galD = smooth(fits_galD[*,0,0],smoothlevel)

wave_galA = indgen(n_elements(spec_galA))*sxpar(hdr_galA,'CD1_1') + sxpar(hdr_galA,'CRVAL1')
wave_galB = indgen(n_elements(spec_galB))*sxpar(hdr_galB,'CD1_1') + sxpar(hdr_galB,'CRVAL1')
wave_galC = indgen(n_elements(spec_galC))*sxpar(hdr_galC,'CD1_1') + sxpar(hdr_galC,'CRVAL1')
wave_galD = indgen(n_elements(spec_galD))*sxpar(hdr_galD,'CD1_1') + sxpar(hdr_galD,'CRVAL1')

; Posited redshifts for the sample

z_galA = 0.094
z_galB = 0.095
z_galC = 0.094
z_galD = 0.094

; Plot spectra

ps_start, file='~/Astronomy/Research/GalaxyZoo/violinclef/mmt_data/vcredshift_new.ps', /color, /landscape, /quiet

!p.multi=[0,1,4]

angstrom = '!6!sA!r!u!9 %!6!n!3'

cleanA = where(abs(spec_galA) lt 150.)
cleanB = where(abs(spec_galB) lt 150.)
cleanC = where(abs(spec_galC) lt 150.)
cleanD = where(abs(spec_galD) lt 150.)

cs = 1.2

plot, wave_galA[cleanA], spec_galA[cleanA], psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,200], title='Violin Clef galaxy A'
vclines, z_galA

plot, wave_galB[cleanB], spec_galB[cleanB], psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,200], title='Violin Clef galaxy B'
vclines, z_galB

plot, wave_galC, spec_galC, psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,800], title='Violin Clef galaxy C'
vclines, z_galC

plot, wave_galD[cleanD], spec_galD[cleanD], psym=10, charsize=cs, xtitle='Observed wavelength ['+angstrom+']',ytitle='Relative flux', /xstyle, xr=[3600,6900], /ystyle, yr=[-100,200], title='Violin Clef galaxy D'
vclines, z_galD

!p.multi=[0,1,1]

ps_end

end
