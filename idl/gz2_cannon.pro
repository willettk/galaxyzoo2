
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

file='~/Astronomy/Research/GalaxyZoo/gz2_cannon.csv'
readcol, file, $
	objid, $
	asset_id, $
	total_count, $
	total_weight, $
	t01_smooth_or_features_a01_smooth_weighted_fraction, $
	t01_smooth_or_features_a03_star_or_artifact_weighted_fraction, $
	t05_bulge_prominence_a10_no_bulge_weighted_fraction, $
	t06_odd_a14_yes_weighted_fraction, $
	t08_odd_feature_a22_irregular_weighted_fraction, $
	/silent, $
	skipline=1, $
	format='a,a,i,f,f,f,f,f,f'

n = n_elements(objid)

shieldsim = where(total_count gt 20 and $
	total_weight gt 0.75*total_count and $
	t01_smooth_or_features_a01_smooth_weighted_fraction gt 0.5 and $
	t01_smooth_or_features_a03_star_or_artifact_weighted_fraction lt 0.1 and $
	t05_bulge_prominence_a10_no_bulge_weighted_fraction gt 0.75 and $
	t06_odd_a14_yes_weighted_fraction gt 0.5 and $
	t08_odd_feature_a22_irregular_weighted_fraction ge 0.4, $
	sscount)

print,sscount,n,float(sscount)/float(n)*100

stop

end
