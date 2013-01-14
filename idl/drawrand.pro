;+
; NAME:
;         DrawRand
;
; PURPOSE:
;        Draw a random set of N values X following the probability
;        distribution function P(X).
;
; CATEGORY:
;
;        Statistical fun
;
; CALLING SEQUENCE:
;
;        result = DrawRand( USERFUNC, NNUM, XRANGE [, /BUILTIN, YMAX =] )
;
; INPUTS:
;
;        USERFUNC -  The name of a user-defined function that
;                    provides P(X). This function must be pre-
;                    compiled. The probablility function P(X)
;                    returned by USERFUNC should be "normalized"
;                    such that the most likely value X = Xo should
;                    produce P(Xo) = 1.0 with all other values of
;                    P(X) scaled accordingly and the minimum value
;                    P(X)_min = 0
;
;                    Alternatively, you may specify the range
;                    of possible P(X) values using the YMAX =
;                    keyword.
;
;                    USERFUNC is ignored if /BUILTIN is set
;        
;        NNUM     -  The number of random values to be returned
;
;        XRANGE   -  The limits on the range of values X can take.
;                    If this input is not supplied, then the
;                    function USERFUNC should return the xrange
;                    when no inputs are supplied to USERFUNC, i.e.
;                    the following, or similar, line should be in
;                    USERFUNC:
;
;                    if n_params() = 0 then return, xrange
;
; OUTPUTS:
;
;       result    -  A random set of NNUM values of X which follow
;                    the PDF defined in USERFUNC.
;
; KEYWORDS:
; 
;       YMAX      -  Specifies the maximum value returned by
;                    USERFUNC.
;       BUILTIN   -  Use the built-in function DRAWRAND_FUNC
;                    to return P(X).
;
; PROCEDURE:
;           For a single value (NNUM = 1):
;
;             1) Choose a uniform random number Xi between XRANGE[0] and
;                XRANGE[1] and a uniform random value Yi between 0 and 1
;             2) Evaluate P(Xi)
;                  a) if P(Xi) > Yi then return Xi
;                  b) else return to step 1
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
; Written by John "JohnJohn" Johnson 23-Apr-2004
; Heavily revised by Jason Wright 24-Apr-2004
;-

function drawrand_func, x
;;; Simple example of a linear pdf with a flat tail
if n_params() eq 0 then return, [0., 19.]
f = fltarr(n_elements(x))
lo = where(x le 19., ct)
if ct gt 0 then f[lo] = poly(x[lo], [1,-.05])
hi = where(x gt 19., ct)
if ct gt 0 then f[hi] = 0.05
return, f
end

function drawrand_uniform, nnum,xrange
;;; Initialize global seed variable
defsysv,'!myseed',exists=exists
if 1-keyword_set(exists) then begin
    dum = randomu(seed)
    defsysv,'!myseed',seed
endif
seed=!myseed
ret=randomu(seed,nnum)*(max(xrange)-min(xrange))+min(xrange)
!myseed=seed
return,ret
end

function drawrand, userfunc, nnum, xrange, ymax=ymax, builtin=builtin

;;; More initialization
if keyword_set(ymax) then yrange = [0, ymax] else yrange = [0,1]
if keyword_set(builtin) then userfunc = 'drawrand_func'
if n_elements(xrange) eq 0 then xrange = call_function(userfunc)
delx = max(xrange)-min(xrange) 

;;; Set up some useful arrays
bad = lindgen(nnum)
pxi = fltarr(nnum)
xi = fltarr(nnum)
yi = fltarr(nnum)
nbad = nnum
cool = 0b                       ;loop control variable
ct = 0
while not cool do begin
    ;;; Get set of x values for values not yet settled
    xi[bad] = drawrand_uniform(nbad, xrange)
    ;;; Get set of y values for values not yet settled
    yi[bad] = drawrand_uniform(nbad, yrange)
    ;;; Get set of comparison values for indices not yet settled
    pxi[bad] = call_function(userfunc, xi[bad])
    ;;; Is pxi > yi?  If not, mark the elements which require
    ;;; changing and rinse, repeat. Dry once there are no "bad"
    ;;; values 
    bad = where(yi gt pxi, nbad, complement = good)
    ngood = n_elements(good)
    ct += 1.
    if nbad eq 0 then cool = 1
endwhile
return, xi
end



