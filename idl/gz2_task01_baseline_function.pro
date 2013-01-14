
;+
; NAME:
;       
;	gz2_task01_baseline_function
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
;       Written by K. Willett                Jun 12
;-

function gz2_task01_baseline_s1, r50, p

	q1 = p[2]
	q2 = p[3]
	q3 = p[4]
	q4 = p[5]
	q5 = p[6]

	s1 = q1^((-1d)*(q2 + q3*(r50^q4))) + q5

	;print,'s1',s1

	return, s1

end

function gz2_task01_baseline_s2, r50, p

	r1 = p[7]
	r2 = p[8]
	q5 = p[6]

	s2 = r1 + r2*(gz2_task01_baseline_s1(r50,p) - q5)

	;print,'s2',s2

	return, s2

end

function gz2_task01_baseline_function, mr, r50, p

	mr = double(mr)
	r50 = double(r50)
	p = double(p)

	p1 = p[0]
	p2 = p[1]
	q1 = p[2]
	q2 = p[3]
	q3 = p[4]
	q4 = p[5]
	q5 = p[6]
	r1 = p[7]
	r2 = p[8]

	;func = p1 / (1d + exp((gz2_task01_baseline_s1(r50, p) - mr)/gz2_task01_baseline_s2(r50, p))) + p2
	func = p1 / (1d + exp((q1^((-1d)*(q2 + q3*(r50^q4))) + q5 - mr)/(r1 + r2*(q1^((-1d)*(q2 + q3*(r50^q4))))))) + p2
	

	;print,'func',func

	return, func

end

function gz2_task01_simplebaseline, mr, r50, p

	mr = double(mr)
	r50 = double(r50)
	p = double(p)

	func = p[0] * mr^2 + p[1] * r50^2

	return, func

end
