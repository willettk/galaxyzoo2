;pro dwarfagn_color_mass

;+
; NAME:
;       
;	
;
; PURPOSE:
;
;	Plot a quick color-mass diagram of the dwarf galaxy AGN + control sample in Ed Moran's sample
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
;	Stellar masses are from the MPA/JHU catalog, which had about a 71% retrieval rate for Ed's objects. Stellar masses are from DR7, though, and it's possible that the missing galaxies are from a newer release. 

; Stellar masses: http://home.strw.leidenuniv.nl/~jarle/SDSS/
;
; REVISION HISTORY
;       Written by K. Willett                28 Nov 2012
;-


dwarfagn_path = '~/Astronomy/Research/GalaxyZoo/dwarfagn/'

sm = mrdfits(dwarfagn_path+'totlgm_dr7_v5_2b.fit',1,/silent)
casid = mrdfits(dwarfagn_path+'CAS-IDs.fits',1,/silent)

moran_control = mrdfits(dwarfagn_path+'moran_control_crossmatch_willettk.fit',1,/silent)
moran_agn     = mrdfits(dwarfagn_path+'moran_agn_crossmatch_willettk.fit',1,/silent)

match,strtrim(casid.objid,2),strtrim(moran_control.dr7objid,2),starind_control,controlind
match,strtrim(casid.objid,2),strtrim(moran_agn.dr7objid,2),starind_agn,agnind

star_matched_control = sm[starind_control]
control_matched = moran_control[controlind]
star_matched_agn = sm[starind_agn]
agn_matched = moran_agn[agnind]

ulim = 21
ugood_control = where(control_matched.petromag_u lt ulim)
ugood_agn = where(agn_matched.petromag_u lt ulim)
star_matched_control = star_matched_control[ugood_control]
control_matched = control_matched[ugood_control]
star_matched_agn = star_matched_agn[ugood_agn]
agn_matched = agn_matched[ugood_agn]

device,retain=2

ps_start, filename=dwarfagn_path + 'colormass.eps',/color,/quiet,xs=5,ys=8,/encap

!p.multi=[0,1,2]
cs = 1
ss = 0.5
cgplot, star_matched_control.median, control_matched.petromag_u - control_matched.petromag_r, $
    psym=16, $
    color='dark green', $
    ytitle='(u - r)', $
    xtitle='Stellar mass [M!Isun!N]',$
    ;yr=[0,5],$
    charsize = cs, $
    symsize=ss

agn_ur_colorerr = sqrt(agn_matched.petromagerr_u^2 + agn_matched.petromagerr_r^2)
agn_starerr_lo = star_matched_agn.median - star_matched_agn.p16
agn_starerr_hi = star_matched_agn.p84 - star_matched_agn.median

oploterror, star_matched_agn.median, $
    agn_matched.petromag_u - agn_matched.petromag_r, $
    agn_starerr_hi, $
    agn_ur_colorerr, $
    /nohat, $
    /lobar, $
    psym=16, $
    color='red', $
    errcolor='red', $
    symsize=ss

oploterror, star_matched_agn.median, $
    agn_matched.petromag_u - agn_matched.petromag_r, $
    agn_starerr_hi, $
    agn_ur_colorerr, $
    /nohat, $
    /hibar, $
    psym=16, $
    color='red', $
    errcolor='red', $
    symsize=ss


al_legend,['Dwarf galaxies (control)','Dwarf AGN'], color=['dark green','red'], /bottom,/right, psym=16, charsize = cs

cgplot, star_matched_control.median, control_matched.petromag_g - control_matched.petromag_r, $
    psym=16, $
    color='dark green', $
    ytitle='(g - r)', $
    xtitle='Stellar mass [M!Isun!N]',$
    ;yr=[0,1.2],$
    charsize = cs, $
    symsize=ss

agn_gr_colorerr = sqrt(agn_matched.petromagerr_g^2 + agn_matched.petromagerr_r^2)

oploterror, star_matched_agn.median, $
    agn_matched.petromag_g - agn_matched.petromag_r, $
    agn_starerr_hi, $
    agn_gr_colorerr, $
    /nohat, $
    /lobar, $
    psym=16, $
    color='red', $
    errcolor='red', $
    symsize=ss

oploterror, star_matched_agn.median, $
    agn_matched.petromag_g - agn_matched.petromag_r, $
    agn_starerr_hi, $
    agn_gr_colorerr, $
    /nohat, $
    /hibar, $
    psym=16, $
    color='red', $
    errcolor='red', $
    symsize=ss

al_legend,['Dwarf galaxies (control)','Dwarf AGN'], color=['dark green','red'], /bottom,/right, psym=16, charsize = cs


; Current problem - which ones are the AGN? Find out and overplot in a different color. Also try Kevin's contours for the main-sequence GZ sample. 
; Outliers? Look at errors on u-band, r-band flux and see why the scatter exists. Need to rerun CasJobs query.

;cghistoplot, control_matched.petromag_u, datacolor='purple', /outline, thick=2, xr=[30,10], xtitle='Petrosian magnitude'
;cghistoplot, control_matched.petromag_g, datacolor='green', /outline, thick=2, /oplot
;cghistoplot, control_matched.petromag_r, datacolor='red', /outline, thick=2, /oplot
;cghistoplot, control_matched.petromag_i, datacolor='orange', /outline, thick=2, /oplot
;cghistoplot, control_matched.petromag_z, datacolor='black', /outline, thick=2, /oplot
;cgplots, [ulim,ulim], !y.crange, /data, linestyle=2
;al_legend, ['u','g','r','i','z'], psym=28, color=['purple','green','red','orange','black'], /top, /right
;

ps_end

end
