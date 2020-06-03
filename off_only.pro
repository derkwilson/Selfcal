; off_only.pro
;====================================================================
; Version: 2004-08-20
; 2004-08-20: NOT BACKWARD COMPATIBLE WITH EARLIER VERSIONS. 
; Uses new uniofm calling sequence and revised pixelization/mapping 
; procedures.
; Version: 2003-12-10
; 2002-06-12: Streamlined Section 3) by passing x,y calculations out of 
; get_pixel_numbers, which also means that the code is now independent of 
; changes internal to get_pixel_numbers.   
; 2002-06-20: Now uses new_pixel_numbers and astrometry structures, and 
; handles IRAC distortions.
; 2002-10-21: Changes for reduced memory and added use of scratch file.
; 2003-12-10: Modified gain_off.pro to omit gain part of solution, along
; the lines of gain_only.pro. A fixed and appropriate gain map is still 
; required as input.
;====================================================================
pro off_only, Di, astrom, astrom0, map_no, Wi, n_terms, Gp, Fp, Fq, $
          skymap, sigmadet, sigmasky, cr_flag, dark_flag, $
          n_sig=n_sig, scale=scale, term_limits=term_limits, $
          pix_no=pix_no, xxs=xxs, yys=yys, xo=xo, yo=yo, $
          umap=umap, psimap=psimap, ipsimap=ipsimap, channel=channel
;
; This procedure applies a least-squares algorithm (Fixsen, Moseley, & Arendt
; 1999) to extract sky intensities, and detector offsets from
; data sets of dithered images. Notation follows that of FM&A 1999.
;
;INPUT:
; Di         = stack (Px,Py,M) of observations
; astrom     = astrometry info for Di
; Wi         = stack (Px,Py,M) of sigma for each observation (aka weights)
; n_terms    = [min,max] numbers of terms to which the sum is carried out.
; Gp         = gain map for the detector array (aka gainmap)
;INPUT/OUTPUT:
; Fp         = offset map for the detector array (aka offmap)
; Fq         = (dummy parameter to maintain same call as other versions)
;OUTPUT:
; skymap     = intensity map for the sky
; sigmadet   = gain rms deviation map for the detector array
; sigmasky   = intensity rms deviation map for the sky
;INPUT/OUTPUT
; cr_flag    = stack (Px,Py,M) of flags indicating CR hits (0=hit,1=OK)
;OPTIONAL INPUT:
; dark_flag  = array of M flags indicating dark frames (0=dark,1=sky)
; n_sig      = number of sigma above test value to trigger CR flag
; scale      = the ratio of output to input pixel sizes
; term_limits= array of 2 lower limits on sizes of terms in the sum relative to
;              [previous term, current total]
;OPTIONAL INPUT/OUTPUT:
; pix_no     = stack of pixels numbers for each datum
; xxs        = x dimension of the array in which the pixel numbers are embedded
; yys        = y dimension of the array in which the pixel numbers are embedded
; xo         = x coordinate of the lower left corner
; yo         = y coordinate of the lower left corner
; umap       = input is the integer index one pixel of the detector array
;              output is the column of the U matix for the chosen detector
;              (not calculated if unspecified or input is out of range e.g. -1)
; psimap     = input is the integer index one observation in the Di stack
;              output is the column of the psi matix for the matching sky pixel
;              (not calculated if unspecified or input is out of range e.g. -1)
; ipsimap    = input is an integer flag
;              output is the diagonal of the Psi^(-1) matrix (a very slow calc)
;              (not calculated if unspecified or input le 0)
; channel    = IRAC channel number (1-4). If supplied, this triggers IRAC
;              array magnification, rotation, and distortion corrections.
;              Omit if no such corrections are desired.
;====================================================================

time0 = systime(1)

;1) Get dimensions of detector and data
;-----------------------------------------
sky_size = size(Di)
Px = sky_size[1]
Py = sky_size[2]
P = sky_size[1]*sky_size[2]
M = sky_size[3]
PxM = sky_size[5]
if (n_elements(dark_flag) ne M) then dark_flag = bytarr(M)+1B
do_umap = 0
if (n_elements(umap) eq 1) then $
  if ((umap ge 0) and (umap lt 2*P)) then do_umap=1
do_psimap = 0
if (n_elements(psimap) eq 1) then $
  if ((psimap ge 0) and (psimap lt PxM)) then do_psimap=1
do_ipsimap = 0
if (n_elements(ipsimap) eq 1) then $
  if (ipsimap gt 0) then do_ipsimap=1

; If n_terms is a single number, assume it is the maximum
if (n_elements(n_terms) eq 1) then maxterms=n_terms else maxterms=max(n_terms)
if (n_elements(n_terms) eq 1) then minterms=5       else minterms=min(n_terms)

sZ_limit = -1e20
sZ2_limit = -1e20
if (n_elements(term_limits) eq 2) then begin
 sZ_limit = term_limits[0]
 sZ2_limit = term_limits[1]
endif

; locate the directory where the compiled external C routines reside
search_path = expand_path(".:"+!path,/array,count=npath)
nfound = 0
i = 0
while (i lt npath) and (nfound ne 1) do begin
  c_dir = findfile(search_path[i]+"/spax.so",count=nfound)
  i = i+1
endwhile
c_dir = c_dir[0]
if (c_dir ne '') then c_dir=strmid(c_dir,0,strlen(c_dir)-strlen("spax.so")) $
  else message,'Compiled C routines (spax.so, etc.) cannot be located.'

help,/mem
print,'Finished 1)'

;2) Get pixels. Sort them. Identify unique pixels in sorted sequence
;-----------------------------------------
pix_no = pixel_numbers(astrom,astrom0,sky_size,map_no,dark_flag,$
         xxs,yys,xo,yo,xpix,ypix,scale=scale,distort_ch=channel)
pix_sort = sort(pix_no)
starts = [0,uniq(pix_no[pix_sort])+1]
Gamma = n_elements(starts)-1
h_sky = where(dark_flag eq 1)
xpix = temporary(xpix[*,*,h_sky])
ypix = temporary(ypix[*,*,h_sky])

help,/mem
print,'Finished 2)'

;3) Build estimated sky and weighted(data-model)
;-----------------------------------------
wG = Wi
for i=0L,M-1 do wG[*,*,i] = Wi[*,*,i]*Gp^2*cr_flag[*,*,i]
skymap = make_maps(Di,pix_no,pix_sort,xxs,yys,Gp,Fp,wG,wtmap)
wG = 0

;fill holes and interpolate skymap to fix data
kern=[[1,1,1],[1,0,1],[1,1,1]]
for i=0,n_tags(skymap)-1 do begin
  ws = convol(skymap.(i)*wtmap.(i),kern,/edge_wrap)
  w = convol(wtmap.(i),kern,/edge_wrap)
  h = where((wtmap.(i) eq 0) and (w gt 0),nh)
  if (nh ge 1) then skymap.(i)[h] = ws[h]/w[h]
endfor

Sa = fltarr(Px,Py,M)
for i=0,n_elements(h_sky)-1 do Sa[*,*,h_sky[i]] = $
  interpolate(skymap.(map_no[h_sky[i]]),xpix[*,*,h_sky[i]],ypix[*,*,h_sky[i]])
xpix=0B & ypix=0B

Di_H0i = Di
for i=0L,M-1 do begin
  Di_H0i[*,*,i] = Wi[*,*,i] * (Di[*,*,i]-Gp*Sa[*,*,i]-Fp) * cr_flag[*,*,i]
endfor

Sa = 0 ; save memory

help,/mem
print,'Finished 3)'

;4) build A^(-1/2) and detector subcomponents of "Y"
;-----------------------------------------
irootA = 1./sqrt(total(Wi,3,/double))
delHD = fltarr(P+Gamma)
delHD[0:P-1] = total(Di_H0i,3,/double)

help,/mem
print,'Finished 4)'

;5) build C^(-1/2) and sky subcomponents of "Y"
; and create pix_ind (0,...,Gamma-1)
;-----------------------------------------
irootC = fltarr(Gamma)
pix_ind = lonarr(Px,Py,M)
for i=0L,Gamma-1 do begin
  here = pix_sort[starts[i]:starts[i+1]-1]
  delHD[P+i] = total(Gp[here mod P]*Di_H0i[here],/double)
  pix_ind[here] = i
  irootC[i] = total(Wi[here]*Gp[here mod P]^2,/double)
endfor
Di_H0i = 0 ; save memory
irootC = 1./sqrt(temporary(irootC))

scrfile='escargot_scratch_'+string(form='(i8.8)',1e8*randomu(s,/double))+'.xdr'
save,Di,Wi,pix_no,pix_sort,file=scrfile
Di = 0

help,/mem
print,'Finished 5)'

;6) build T and reform T, irootA, and build diagonal of weights
;-----------------------------------------
T = temporary(Wi)
for i=0L,M-1 do T[*,*,i] = irootA*T[*,*,i]*Gp
for i=0L,Gamma-1 do begin
  here = pix_sort[starts[i]:starts[i+1]-1]
  T[here] = T[here] * irootC[i]
endfor

; C^(-1/2) Psi^(-1) C^(-1/2) = I - (T^t T)
if (do_ipsimap) then begin
  ipsimap = skymap
  for j=0,n_tags(ipsimap)-1 do ipsimap.(j) = 0.0
  n_cumul = lonarr(n_tags(ipsimap))
  for j=0,n_elements(n_cumul)-1 do n_cumul[j]=n_elements(ipsimap.(j))
  n_cumul = [0L,long(total(n_cumul,/cumul))]
  for i=0L,Gamma-1 do begin
    here = pix_sort[starts[i]:starts[i+1]-1]
    comp1 = total(T[here]^2,/double)
    region = max(where(n_cumul le pix_no[here[0]]))
    pix = pix_no[here[0]] - n_cumul[region]
    ipsimap.(region)[pix] = (1-comp1)/irootC[i]^2
  endfor
endif

pix_no = 0 ; save memory
pix_sort = 0 ; save memory

irootA = reform(irootA,P)

help,/mem
print,'Finished 6)'

;7) calculate Z0, Zn = T(Tt*Zn) + Z0, and Xp
;-----------------------------------------
stime0 = systime(1)
Z0 = fltarr(P)
status = call_external(c_dir+'spax.so','spax',T,pix_ind,$
  irootC*delHD[P:*],Z0,P,PxM)
Z0 = -Z0 + irootA*delHD[0:P-1]
;NOT OFFSETS! Z0[P:*] = Z0[P:*] - total(Z0[P:*],/double)/P  ; normalization
Z0 = Z0 - total(Z0,/double)/P  ; normalization
Z0 = float(Z0)
Zn = Z0
dZ0 = Z0
sZ = fltarr(maxterms+1)
sZ2 = fltarr(maxterms+1)

if (do_umap) then begin
  upix = fltarr(P)
  upix[umap] = 1.0
  zu0 = upix
  zun = upix
endif
if (do_psimap) then begin
  ualpha = fltarr(Gamma)
  ualpha[pix_ind[psimap]] = 1.0
  zeta0 = ualpha
  zetan = ualpha
endif

i = 0L
while ( ((abs(sZ[i]) gt sZ_limit) and (abs(sZ2[i]) gt sZ2_limit) and $
  (i lt maxterms)) or (i lt minterms) ) do begin
  i=i+1
  TTtZn = fltarr(P)
  TtZn = fltarr(Gamma)
  status = call_external(c_dir+'sptx.so','sptx',T,pix_ind,$
                          Zn,TtZn,P,PxM)
  status = call_external(c_dir+'spax.so','spax',T,pix_ind,$
                          TtZn,TTtZn,P,PxM)
;NOT OFFSETS!  TTtZn[P:*] = TTtZn[P:*] - total(TTtZn[P:*],/double)/P ; norm
  TTtZn[*] = TTtZn - total(TTtZn,/double)/P ; normalization
  dZ1 = TTtZn+Z0-Zn
  dZ1 = dZ1 - total(dZ1,/double)/P
  sZ[i] = total(dZ0*dZ1,/double)/total(dZ0^2,/double)
  dZ0 = dZ1
  Zn = TTtZn+Z0
  sZ2[i] = total(Zn*dZ1,/double)/total(Zn^2,/double)
  print,format='($,a,i4,1x,e12.3,1x,e12.3,a)','(TT)^n, n=',$
         i,sZ[i],sZ2[i],string(13b)

  if (do_umap) then begin
    TTtzun = fltarr(P)
    Ttzun = fltarr(Gamma)
    status = call_external(c_dir+'sptx.so','sptx',T,pix_ind,$
                        zun,Ttzun,P,PxM)
    status = call_external(c_dir+'spax.so','spax',T,pix_ind,$
                        Ttzun,TTtzun,P,PxM)
    zun = TTtzun+zu0
  endif
  if (do_psimap) then begin
    TtTzetan = fltarr(Gamma)
    Tzetan = fltarr(P)
    status = call_external(c_dir+'spax.so','spax',T,pix_ind,$
                        zetan,Tzetan,P,PxM)
    status = call_external(c_dir+'sptx.so','sptx',T,pix_ind,$
                          Tzetan,TtTzetan,P,PxM)
    zetan = TtTzetan+zeta0
  endif
endwhile
T = 0 ; save memory
pix_ind = 0 ; save memory
restore,scrfile
file_delete,scrfile

print,''
!p.multi = [0,2,2]
plot,sZ,psym=-1,yran=[0,1],ytitle='sZ',xtitle='N'
plot_io,sZ2,psym=-1,ytitle='sZ2',xtitle='N',yran=[1e-6,1e2]
!p.multi = 0
Xp = irootA*Zn
if (do_umap) then begin
  umap = reform(zun,Px,Py)
  print,'max(zun) = ',max(zun)
  print,'min(zun) = ',min(zun)
endif
if (do_psimap) then begin
  psi = total(irootC*ualpha,/double)*zetan
  print,'irootc(alpha) = ',total(irootC*ualpha,/double)
  print,'max(zetan) = ',max(zetan)
  psimap = skymap
  for j=0,n_tags(psimap)-1 do psimap.(j) = 0.0
  n_cumul = lonarr(n_tags(psimap))
  for j=0,n_elements(n_cumul)-1 do n_cumul[j]=n_elements(psimap.(j))
  n_cumul = [0L,long(total(n_cumul,/cumul))]
  for i=0L,Gamma-1 do begin
    pix = pix_no[pix_sort[starts[i]]]
    region = max(where(n_cumul le pix))
    pix = pix - n_cumul[region]
    psimap.(region)[pix] = psi[i]
  endfor
endif
stime1 = systime(1)
print,'series time (s) = ',stime1-stime0

help,/mem
print,'Finished 7)'

;8) find final gain and offset (Gp, Fp)
;-----------------------------------------
Fp = Fp+reform(Xp[0:P-1],Px,Py)
skymap = make_maps(Di,pix_no,pix_sort,xxs,yys,Gp,Fp,Wi*cr_flag,wtmap)

help,/mem
print,'Finished 8)'

;9) find uncertainties and make final map
;-----------------------------------------
make_sigs,Di,pix_no,pix_sort,map_no,xxs,yys,Gp,Fp,skymap,Wi,cr_flag,$
  sigmadet,sigmasky,chi2,n_sig=n_sig
skymap = make_maps(Di,pix_no,pix_sort,xxs,yys,Gp,Fp,Wi*cr_flag,wtmap)

help,/mem
time1 = systime(1)
print,' Total time (s) = ',time1-time0

end
