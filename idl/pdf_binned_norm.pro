
;+
; NAME:
;       
;	PDF_BINNED_NORM
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

function pdf_binned_norm, h, bins, x

nb = n_elements(bins)
minval=min(bins)
maxval=max(bins)
binsize = (maxval - minval)/float(nb)
hnorm = h/float(max(h))

xh = histogram(float(x),binsize=binsize,min=minval,max=maxval,reverse_indices=ri_xh)
xnew = x*0.

for i=0,nb-1 do begin
    if xh[i] gt 0 then begin
        xtemp = ri_xh[ri_xh[i]:ri_xh[i+1]-1]
        xnew[xtemp] = hnorm[i]
    endif
endfor

return,xnew

end
