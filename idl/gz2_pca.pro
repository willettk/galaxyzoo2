
;+
; NAME:
;       
;	GZ2_PCA
;
; PURPOSE:
;
;	Principal component analysis on GZ2 parameters
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

; Use NA10 classifications as "truth"

;na = mrdfits('~/Astronomy/Research/GalaxyZoo/fits/na_gz2_willettk.fit',1,/silent)

tagnames=tag_names(na)
nt = n_elements(tagnames)

temparr=intarr(nt) 
for i=0,nt-1 do if (strlen(tagnames[i]) ge 18) and (strmid(tagnames[i],strlen(tagnames[i])-17,17) eq 'WEIGHTED_FRACTION') then temparr[i]=1 

wfind = [0,where(temparr)]	; Important: index 1 corresponds to answer 01 ONLY for less than answer 24
wcind = [0,where(temparr)]	

;gz2_arr = [$
;	[na.(wfind[1])], $                  1 T01_SMOOTH_OR_FEATURES_A01_SMOOTH_WEIGHTED_FRACTION           
;	[na.(wfind[2])], $                  2 T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION
;	[na.(wfind[3])], $                  3 T01_SMOOTH_OR_FEATURES_A03_STAR_OR_ARTIFACT_WEIGHTED_FRACTION
;	[na.(wfind[4])], $                  4 T02_EDGEON_A04_YES_WEIGHTED_FRACTION
;	[na.(wfind[5])], $                  5 T02_EDGEON_A05_NO_WEIGHTED_FRACTION
;	[na.(wfind[6])], $                  6 T03_BAR_A06_BAR_WEIGHTED_FRACTION
;	[na.(wfind[7])], $                  7 T03_BAR_A07_NO_BAR_WEIGHTED_FRACTION
;	[na.(wfind[8])], $                  8 T04_SPIRAL_A08_SPIRAL_WEIGHTED_FRACTION
;	[na.(wfind[9])], $                  9 T04_SPIRAL_A09_NO_SPIRAL_WEIGHTED_FRACTION
;	[na.(wfind[10])], $                10 T05_BULGE_PROMINENCE_A10_NO_BULGE_WEIGHTED_FRACTION
;	[na.(wfind[11])], $                11 T05_BULGE_PROMINENCE_A11_JUST_NOTICEABLE_WEIGHTED_FRACTION
;	[na.(wfind[12])], $                12 T05_BULGE_PROMINENCE_A12_OBVIOUS_WEIGHTED_FRACTION
;	[na.(wfind[13])], $                13 T05_BULGE_PROMINENCE_A13_DOMINANT_WEIGHTED_FRACTION
;	[na.(wfind[14])], $                14 T06_ODD_A14_YES_WEIGHTED_FRACTION
;	[na.(wfind[15])], $                15 T06_ODD_A15_NO_WEIGHTED_FRACTION
;	[na.(wfind[16])], $                16 T07_ROUNDED_A16_COMPLETELY_ROUND_WEIGHTED_FRACTION
;	[na.(wfind[17])], $                17 T07_ROUNDED_A17_IN_BETWEEN_WEIGHTED_FRACTION
;	[na.(wfind[18])], $                18 T07_ROUNDED_A18_CIGAR_SHAPED_WEIGHTED_FRACTION
;	[na.(wfind[19])], $                19 T08_ODD_FEATURE_A19_RING_WEIGHTED_FRACTION
;	[na.(wfind[20])], $                20 T08_ODD_FEATURE_A20_LENS_OR_ARC_WEIGHTED_FRACTION
;	[na.(wfind[21])], $                21 T08_ODD_FEATURE_A21_DISTURBED_WEIGHTED_FRACTION
;	[na.(wfind[22])], $                22 T08_ODD_FEATURE_A22_IRREGULAR_WEIGHTED_FRACTION
;	[na.(wfind[23])], $                23 T08_ODD_FEATURE_A23_OTHER_WEIGHTED_FRACTION
;	[na.(wfind[24])], $                24 T08_ODD_FEATURE_A24_MERGER_WEIGHTED_FRACTION
;	[na.(wfind[25])], $                25 T08_ODD_FEATURE_A38_DUST_LANE_WEIGHTED_FRACTION
;	[na.(wfind[26])], $                26 T09_BULGE_SHAPE_A25_ROUNDED_WEIGHTED_FRACTION
;	[na.(wfind[27])], $                27 T09_BULGE_SHAPE_A26_BOXY_WEIGHTED_FRACTION
;	[na.(wfind[28])], $                28 T09_BULGE_SHAPE_A27_NO_BULGE_WEIGHTED_FRACTION
;	[na.(wfind[29])], $                29 T10_ARMS_WINDING_A28_TIGHT_WEIGHTED_FRACTION
;	[na.(wfind[30])], $                30 T10_ARMS_WINDING_A29_MEDIUM_WEIGHTED_FRACTION
;	[na.(wfind[31])], $                31 T10_ARMS_WINDING_A30_LOOSE_WEIGHTED_FRACTION
;	[na.(wfind[32])], $                32 T11_ARMS_NUMBER_A31_1_WEIGHTED_FRACTION
;	[na.(wfind[33])], $                33 T11_ARMS_NUMBER_A32_2_WEIGHTED_FRACTION
;	[na.(wfind[34])], $                34 T11_ARMS_NUMBER_A33_3_WEIGHTED_FRACTION
;	[na.(wfind[35])], $                35 T11_ARMS_NUMBER_A34_4_WEIGHTED_FRACTION
;	[na.(wfind[36])], $                36 T11_ARMS_NUMBER_A36_MORE_THAN_4_WEIGHTED_FRACTION
;	[na.(wfind[37])]]                  37 T11_ARMS_NUMBER_A37_CANT_TELL_WEIGHTED_FRACTION

gz2_arr = [$
	[na.(wfind[1])], $
	[na.(wfind[2])]]

etg = gz2_arr[where(na.ttype eq -5),*]

pca, etg, eigenvalues, eigenvectors, percentages, proj_obj, proj_atr, textout='~/Astronomy/Research/GalaxyZoo/gz2_pca',/silent

cgplot, etg[*,0] - mean(etg[*,0]), etg[*,1] - mean(etg[*,1]), psym=3, xr=[-1,1], yr=[-1,1]
cgplots, [0,eigenvectors[0,0]], [0,eigenvectors[0,1]], thick=3, color='red'
cgplots, [0,eigenvectors[1,0]], [0,eigenvectors[1,1]], thick=3, color='blue'

stop

end
