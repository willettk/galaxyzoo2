
;+
; NAME:
;       
;	FIND99
;
; PURPOSE:
;
;	Compute the GZ2 vote fraction threshold for the previous task at which 99% of responses
;       to the subsequent task will be kept with N >= vote limit for previous task. 
;
;       I know it's complicated. But I think it works. 
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
;       Usually called by GZ2_TASK_HISTOGRAMS.pro
;
; REVISION HISTORY
;       Written by K. Willett                27 Mar 2013
;-

function find99, thistask_votes, prev_response_wf, prev_task_votes, lim

    distribution_all = thistask_votes[where(prev_task_votes ge lim)]

    wflim = 0.0
    fraclim = 1.0
    while fraclim gt 0.99 do begin
        fraclim = n_elements(where(thistask_votes[where(prev_response_wf ge wflim and prev_task_votes ge lim)] gt lim)) / float(n_elements(where(distribution_all gt lim)))
        wflim += 1d-3
    endwhile

return, string(wflim,format='(f5.3)')

end
