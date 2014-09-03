function fans_res,mono=mono,coll1=coll1,coll2=coll2

; to call: 
;	data = fans_res(mono='Cu', coll1=20,coll2=20)
;
;# 'mono' is either 'Cu' or 'Pg' for the Cu(220) or Pg(002) monochromators
;# coll1 is the collimation before the monochromator. Must be 10,20,40 or 60
;# coll2 is the collimation after the monochromator. Must be 10,20,40 or 60
; note using FANS the typical choice is 20'-20'

; returned array is (300,2). 
; column 0 is the energy value, column 1 the resolution in meV (=8.0655 cm-1)

compile_opt hidden

;# START RRESOLTION
d_cu=1.278                 ; Cu(220) d-spacing A
d_pg=3.354					; Pg(002) d-spacing A

;# an initial guesses
dmono=d_cu

case mono of
    'Cu': begin
	dmono=d_cu
	resol=fltarr(350,2)
    end
    'Pg': begin
	dmono=d_pg
	resol=fltarr(45,2)
    end
endcase


col1=coll1/60.0/57.269
col2=coll2/60.0/57.269
rmos=30./60./57.269
tmpd=col1^2+col2^2+4.*rmos^2
tmpn=(col1*col2)^2+(col1*rmos)^2+(col2*rmos)^2
dtheta=sqrt(tmpn/tmpd)

for i=0,(size(resol))[1]-1 do begin
    en=i+1
    resol(i,0)=en          ;the energy array
    wavel=sqrt(81.805/en)
    stheta=wavel/(2.*dmono)
    if (abs(stheta) lt 1) then begin
       cot_theta=sqrt(1-stheta*stheta)/stheta
    endif else  begin
       cot_theta=0.0
    endelse
    de=2.*en*cot_theta*dtheta
    resol(i,1)=sqrt(de*de+1.2*1.2)
endfor


return, resol
end
