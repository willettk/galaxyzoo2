
;+
; NAME:
;       
;	VCLINES
;
; PURPOSE:
;
;	Overplot common lines in red SDSS galaxies to identify redshifts in the Violin Clef system
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

pro vclines, z

if n_elements(z) eq 0 then z = 0.

; Absorption features: Ca II K, Ca II H, G, MgI, NaI D

wavearr = [3934.8,3969.6,4305.6,5176.7,5895.6]
labelpos = [3800, 4000, 4320, 5200, 5920]
namearr = ['CaII K','CaII H','G', 'MgI','NaI D']
nw = n_elements(wavearr)

; Plot on current graph 

cs2 = 0.7

for i=0,nw-1 do begin
	cgplot, /over, $
		[wavearr[i],wavearr[i]] * (1. + z), !y.crange, $
		linestyle=2, color="Red"
	cgtext, labelpos[i]*(1. + z), $
		(!y.crange[1]-!y.crange[0]) * 0.8 + !y.crange[0], $
;		(i mod 2)*(!y.crange[1]-!y.crange[0])*0.1, $
		namearr[i], color="Red", /data, charsize=cs2
endfor


end
