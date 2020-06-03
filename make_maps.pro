function make_maps, skies, pix_no, pix_sort, xxs, yys, $
  gainmap, offmap, weights, quadmap=quadmap, varsky=varsky, wtmap,$
                    naive=naive
;
; This function constructs the weighted mean sky map
;
;INPUT:
; skies    = stack (n,n,m) of observations 
; pix_no   = stack of pixels numbers for each datum
; pix_sort = sort(pix_no)
; xxs      = x dimension of the array in which the pixel numbers are embedded
; yys      = y dimension of the array in which the pixel numbers are embedded
; gainmap  = gain map for the detector array
; offmap   = offset map for the detector array
; weights  = stack (n,n,m) of weights for each observation
;OPTIONAL INPUT:
; quadmap  = [2,2,M] array of variable offsets per image quadrant
;         OR [4,M] array of variable offsets per image column (mod 4)
;         OR [M] array of variable offsets per image
; varsky   = [M] array of variable sky intensities per image
;         OR [Nsub,M] array of variable sky intensities (by X rows) per image
;OPTIONAL OUTPUT:
; wtmap    = map of sum(weights) for each sky pixel.

; This program calibrates and averages data to produce a mean sky map
; on the specified pixels. 
;
; The code will recognize which type of variable offset is supplied (if any)
; by the dimensions of quadmap.
;
; In the 2001-10-30 distribution, this function was attached to the
; gain_only.pro, gain_off.pro, and gain_off_quad.pro functions. 
; Now, it is separate because the generalization for the variable offset 
; seems to be fairly simple and not too inefficient. 

naive=0


varoffsettype = (size(quadmap))[0]
varoffset = 0.0
varskytype = (size(varsky))[0]
vsky = 0.0

Px = n_elements(gainmap[*,0])
Py = n_elements(gainmap[0,*])
P = n_elements(gainmap)
Px2 = Px/2
P2 = P/2

print,'MAKING MAPS'
print,Px,Py
print,xxs
print,yys

nmaps = n_elements(xxs)
skymap = {region0: fltarr(xxs[0],yys[0])}
for i=1,nmaps-1 do skymap = create_struct(skymap,'region'+strtrim(i,2),$
  fltarr(xxs[i],yys[i]))
wtmap = skymap
starts = [0,uniq(pix_no[pix_sort])+1]
Gamma = n_elements(starts)-1

;ketron
;n_cumul = long(total(xxs*yys,/cumulative)) ;original
n_cumul = long64(total(double(xxs)*double(yys),/cumulative))
;stop
for i=0L,Gamma-1 do begin
  h = pix_sort[starts[i]:starts[i+1]-1]
  hmodP = h mod P
  case (varoffsettype) of
       1:  varoffset = quadmap[h/P]
       2:  varoffset = quadmap[h mod 4,h/P]
       3:  varoffset = quadmap[(h mod Px)/Px2,hmodP/P2,h/P]
    ELSE: ;varoffset = 0.0
  endcase
  case (varskytype) of
       1:  vsky = varsky[h/P]
       2:  begin
           Nsub = (size(varsky))[1]
           vsky = varsky[(h mod P)/Px/(Py/Nsub),h/P]
           end
       3:  message, 'Varskytype=3 was commented out';vsky = varsky[(h mod Px)/Px2,hmodP/P2,h/P]
    ELSE: ;vsky = 0.0
  endcase
  j = min(where(pix_no[pix_sort[starts[i]]]/n_cumul eq 0))
;  stop
  if (j eq 0) then begin 
    skymap.(j)[pix_no[pix_sort[starts[i]]]] = $
      total( ((skies[h]-offmap[hmodP]-varoffset)/gainmap[hmodP]-vsky)*weights[h] )
    wtmap.(j)[pix_no[pix_sort[starts[i]]]] = total(weights[h])
  endif else begin
    skymap.(j)[pix_no[pix_sort[starts[i]]] - n_cumul[j-1]] = $
      total( ((skies[h]-offmap[hmodP]-varoffset)/gainmap[hmodP]-vsky)*weights[h] )
    wtmap.(j)[pix_no[pix_sort[starts[i]]] - n_cumul[j-1]] = total(weights[h])
  endelse
endfor
for j=0,nmaps-1 do begin
  h = where(wtmap.(j) ne 0)
  skymap.(j)[h] = skymap.(j)[h]/wtmap.(j)[h]
endfor
if naive then begin
   ;naive_name = '/big_book/dnwilson/spitzer_nep_ch2/results_no_ep1_cbcd/naive'+string(randomn(seed)+10)+'.fits' ;add l=5 as attribute of string
   naive_name = '/big_book/dnwilson/spitzer_nep_ch1/transfer_function/results_temp3/naive'+string(randomn(seed)+10)+'.fits' ;add l=5 as attribute of string  
   mwrfits, skymap.region0, naive_name, /create
  ; pwd
   message, 'wrote '+ naive_name, /inf
;   stop
endif
return,skymap
end
