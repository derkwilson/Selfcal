pro make_sigs, skies, pix_no, pix_sort, map_no, $
              xxs, yys, gainmap, offmap, skymap, $
              weights, cr_flag, sigdet, sigsky, chi2, $
              quadmap=quadmap, varsky=varsky, n_sig=n_sig

; This function uses estimates of the rms errors wrt pixel and sky coord 
; to calculate new weights, new CR flags, and the chi-squared statistic.  

;INPUT:
; skies    = stack (n,n,m) of observations 
; pix_no   = stack of pixels numbers for each datum
; pix_sort = sort(pix_no)
; gainmap  = gain map for the detector array
; offmap   = offset map for the detector array
; skymap   = intensity map for the sky
;INPUT/OUTPUT
; weights  = stack (n,n,m) of weights for each observation
; cr_flag  = stack (Px,Py,M) of flags indicating CR hits (0=hit,1=OK)
;OUTPUT:
; sigdet   = gain rms deviation map for the detector array
; sigsky   = intensity rms deviation map for the sky
; chi2     = [Chi-squared statistic, number of degrees of freedom]
;OPTIONAL INPUT
; quadmap  = [2,2,M] array of variable offsets per image quadrant
;         OR [4,M] array of variable offsets per image column (mod 4)
;         OR [M] array of variable offsets per image
; n_sig    = number of sigma above test value to trigger CR flag (default=3)

; WARNING! 
; The uncertainites calculated here are not the formal uncertainites that 
; the algorithm can ultimately provide. The formal uncertainties are 
; optionally calculated in the main code as umap and psimap.

varoffsettype = (size(quadmap))[0]
varskytype = (size(varsky))[0]
if (n_elements(n_sig) ne 1) then n_sig = 3.
s = size(skies)
P = s[1]*s[2]
starts = [0,uniq(pix_no[pix_sort])+1]
Gamma = n_elements(starts)-1
n_cumul = long64(total(double(xxs)*double(yys),/cumulative))

deldet = skies
sigdet = fltarr(s[1],s[2])
sigsky = skymap
for i=0,n_tags(sigsky)-1 do sigsky.(i) = 0.0

for k=0L,s[3]-1 do begin
  case (varoffsettype) of
       1:  varoffset = quadmap[k]
       2:  varoffset = reform(quadmap[lindgen(s[1],s[2]) mod 4,k],s[1],s[2])
       3:  varoffset = congrid(quadmap[*,*,k],s[1],s[2])
    ELSE:  varoffset = 0.0
  endcase
  case (varskytype) of
       1:  vsky = varsky[k]
       2:  begin
           Nsub = (size(varsky))[1]
           vsky = varsky[lindgen(s[1],s[2])/(P/Nsub),k]
           end
;       3:  vsky = congrid(varsky[*,*,k],s[1],s[2])
    ELSE:  vsky = 0.0
  endcase
  if (map_no[k] eq 0) then begin
    deldet[*,*,k] = skies[*,*,k] - $
      gainmap*(skymap.(0)[pix_no[*,*,k]]+vsky) - offmap - varoffset
  endif else begin
    deldet[*,*,k] = skies[*,*,k] - offmap - varoffset - $
      gainmap*(skymap.(map_no[k])[pix_no[*,*,k]-n_cumul[map_no[k]-1]]+vsky)
  endelse
endfor

for i=0L,s[1]-1 do for j=0L,s[2]-1 do begin
  sigdet[i,j] = resistant_sigma(deldet[i,j,*],n_sig,mean,wt1)
  cr_flag[i,j,*] = wt1
endfor

for i=0L,Gamma-1 do begin
  h = pix_sort[starts[i]:starts[i+1]-1]
  if (n_elements(h) gt 1) then begin
    j = min(where(pix_no[h[0]]/n_cumul eq 0))
    if (j eq 0) then begin
      sigsky.(j)[pix_no[pix_sort[starts[i]]]] = $
        resistant_sigma(deldet[h],n_sig,mean,wt2)
    endif else begin
      sigsky.(j)[pix_no[pix_sort[starts[i]]]-n_cumul[j-1]] = $
        resistant_sigma(deldet[h],n_sig,mean,wt2)
    endelse
  endif else begin
    wt2 = 1
  endelse
  cr_flag[h] = cr_flag[h]+2*wt2
endfor

for k=0L,s[3]-1 do begin 
  if (map_no[k] eq 0) then layer = sigsky.(0)[pix_no[*,*,k]] else $
    layer = sigsky.(map_no[k])[pix_no[*,*,k]-n_cumul[map_no[k]-1]] 
  medsig = median(sigsky.(map_no[k])[where(sigsky.(map_no[k]) ne 0)])
  weights[*,*,k] = 1.0/(sigdet^2 + (layer + 2*medsig*(layer eq 0))^2)
endfor

;sigall = resistant_sigma(deldet,5.0,mean,wt3)
;cr_flag = cr_flag+4*wt3
minweights = min(weights)
medgain = median(gainmap)
siggain = resistant_sigma(gainmap,3.)
here = where(abs(gainmap-medgain) gt 10.*siggain,nhere)
if (nhere ge 1) then for i=0L,nhere-1 do $
  weights[here[i] mod s[1],here[i]/s[1],*] = minweights/1e8

cr_flag = cr_flag ne 0

; The chi^2 calculation is diagnostic only, 
; and has no impact on further iterations.

chi = deldet^2*weights
chi2 = total(chi*cr_flag)
nu = total(cr_flag)-Gamma-P-P*(total(abs(offmap)) ne 0)-n_elements(quadmap)

print,''
print,'  chi^2 = ',strtrim(chi2,2),'  chi^2/nu = ',strtrim(chi2/nu,2),$
  '  nu = ',strtrim(long64(nu),2),'  n_cr = ',strtrim(long64(total(cr_flag eq 0)),2)
print,''

chi2 = [chi2,nu]

end
