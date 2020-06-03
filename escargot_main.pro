;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;escargot_irac.pro
;
; This is an example of a driver procedure for processing data from 
; a single IRAC channel using the least-squares calibration code.  
; It was originally developed for use with processing raw data from separate 
; hi and lo brightness fields, but over time has acquired capabilities for 
; processing BCD data as well. 
; 
; This particular version has been adjusted to process
; IRAC BCD to determine delta-corrections to the offsets.
;
; Most of the parts a user would typically want or need to change are the 
; control parameters and input/output file specification in the first
; section. 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read input file info, then images, then initial gain & offset
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; DEFAULT CONTROL PARAMETERS FOR THE PROCEDURE
; version of the least sqaures procedure to use.
if (n_elements(procedure) ne 1) then procedure = 'off_only_1ff'
; number of iterations of the least-squares procedure
if (n_elements(n_iter) ne 1) then n_iter = 2
; number of term in series expansion to use
if (n_elements(n_terms) ne 1) then n_terms= 30
; n*sigma outlier level at which data get flagged as cosmic rays
if (n_elements(n_sig_cr) ne 1) then n_sig_cr=7.0
; number of data per detector to initially flag as CR hits (>=1)
if (n_elements(n_cr_init) ne 1) then n_cr_init = 1
; ratio of output (skymap)/input (detector) pixel sizes.
if (n_elements(scale0) ne 1) then scale0 = 1.0
; index of X_p (GOQ) for calc. of covariance, -1=none 
if (n_elements(umap0) ne 1) then umap0 = -1
; index of D_i => S_alpha for calc. of covariance, -1=none 
if (n_elements(psimap0) ne 1) then psimap0 = umap0
; flag for calculation of diagonal of covariance matrix.
if (n_elements(ipsimap0) ne 1) then ipsimap0 = -1 
; flag to reuse initial(old) weights, instead of those from make_sig
if (n_elements(use_old_weights) ne 1) then use_old_weights = 1
; flag to map data in celestial coordinates (else array coords of the hi data)
if (n_elements(map_celest) ne 1) then map_celest = 0

; "tags" are partial dir/filenames used to construct the input and 
; output file names 
; All matching files will be read so bad data and odd frametimes need to 
; be moved elsewhere first.
; In this example we are only processing one field, so tag_lo is given
; a dummy name and subsequently is ignored.
; In this example (and most IRAC cases) we read no dark frames.
frametime =  12.         ; frame time of the data to be processed (dummy param.)
tag_hi = '*bcd.fits'  ; file names for input hi brightness field
tag_lo = 'nosuchfile'    ; file names for input lo brightness field
tag_dark = 'nosuchfile'  ; file fames for input dark frames



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;User input: Location of ouput files;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;base name for output files
tag_save = '/big_book/dnwilson/spitzer_nep_ch1/differenced_c1_cbcd/epoch1-epoch2_skymaps/result'  ; base name for output files

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find and read photometry data 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; "hi" files must exist
;hi_files = findfile(tag_hi,count=nhi)
;hi_files1 = sort_fits_irac('gs','ch1',/hi)
;hi_files1=hi_files1.list1
;hi_files2 = sort_fits_irac('gs','ch1',/lo)
;hi_files2=hi_files2.list1
;hi_files = [hi_files1,hi_files2]
;nhi = n_elements(hi_files)
;files = hi_files


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;User input: location of files to be mosaicked;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;gets files from folder, nhi = number of files
files = file_search('/big_book/dnwilson/spitzer_nep_ch1/differenced_ch1_cbcd/epoch1-epoch2/*.fits',count=nhi) ;gets files from folder, nhi = number of files
print,size(files,/n_elements)
hi_files = files

;;dimensions of data tiles
XPIX = 256
YPIX = 256

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; "lo" files may or may not exist.
lo_files = findfile(tag_lo,count=nlo)
nsky = nhi + nlo
if (nlo gt 0) then files = [files,lo_files] 

print,'nhi,nlo = ',nhi,nlo

; "dark" files may or may not exist.
dfiles = findfile(tag_dark,count=ndark)
nframes = nsky+ndark
dark_flag = bytarr(nframes)+1
if (ndark gt 0) then begin 
  files = [files,dfiles]
  dark_flag[nsky:*] = 0
endif

skies = fltarr(XPIX,YPIX,nframes)
for i=0L,nframes-1 do begin
  skies[*,*,i] = readfits(files[i])
endfor
files_used = files



; Preprocess if raw data, (without SIPPFLAG in header)
head = headfits(files[nsky-1])
sippflag = sxpar(head,'SIPPFLAG',count=sip)
bcdflag = sxpar(head,'RAWFILE',count=bcd)
chnlnum = sxpar(head,(sip eq 0) ? 'CHNLNUM': 'ACHANID')
print,'CHNLNUM = ',chnlnum
; Processing got raw data only
if (sip+bcd eq 0) then begin
 ;photometry for channels 1 & 2 raw data has to be negated.
  if (chnlnum le 2) then skies = -skies
  ;images for channels 1 & 2 need to be flipped in the y direction
  if (chnlnum le 2) then for i=0,nframes-1 do $
    skies[*,*,i] = rotate(skies[*,*,i],7)
  ;convert from signed to unsigned int & remove wrap-around (in one step)
  h = where(skies lt -8000, nh)
  if (nh ge 1) then skies[h] = skies[h]+XPIX*YPIX ;;;???? 256.^2 for irac
endif

; replace NaN with zeros, weights set to zero later
hnan = where(finite(skies,/nan),nh)
if (nh ge 1) then skies[hnan] = 0.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read and adjust astrometry
; Use BCD pointing-refined keywords (*RFND) if present
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



extast,headfits(hi_files[0]),astr_blank

; Ingest from header WCS, or else saved SIP astr file
astrhi = replicate(astr_blank,nhi)
extast,head,ao,nopar
if (nopar ne -1) then begin	; headers contain WCS info, use them
  astrhi = replicate(astr_blank,nhi)
  for i=0,nhi-1 do begin
    hd = headfits(hi_files[i])
    extast,hd,aaa
    astrhi[i] = aaa
  endfor
endif else begin                ; headers lack WCS info
  print,'No astrometry found in headers'
  STOP  
endelse
astr0 = astrhi[0]
astrsave = astrhi
map_no = lonarr(nhi)

if (nlo ge 1) then begin
  astrlo = replicate(astr_blank,nlo)
  if (nopar ne -1) then begin     ; headers contain WCS info, use them
    astrlo = replicate(astr_blank,nlo)
    for i=0,nlo-1 do begin
      hd = headfits(lo_files[i])
      extast,hd,aaa
      astrlo[i] = aaa
    endfor
  endif else begin                ; headers lack WCS info
    print,'No astrometry found in headers'
    STOP  
  endelse
  astr0 = [astr0,astrlo[0]]
  astrsave = [astrsave,astrlo]
  map_no = [map_no,lonarr(nlo)+1]
endif
astr0.distort.name = 'IGNORE'
astr = temporary(astrsave)

;make convenient astrometry for darks.
if (ndark ge 1) then astr = [astr,replicate(astr[0],ndark)]

if (keyword_set(map_celest)) then astr0.cd = [[-1,0],[0,1]]*1.2/3600.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Organize input
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



meds = fltarr(4,nframes)
medsky = median(skies[*,*,0:nsky-1])
if (ndark ge 1) then meddark = median(skies[*,*,nsky:*]) else meddark = 0
for i=0,nsky-1 do begin
  for k=0,3 do begin
    meds[k,i] = median(skies[lindgen(64)*4+k,*,i])-medsky
;    skies[lindgen(64)*4+k,*,i]=skies[lindgen(64)*4+k,*,i]-meds[k,i]
  endfor
endfor
if (ndark ge 1) then begin 
  for i=nsky,nframes-1 do begin
    for k=0,3 do begin
      meds[k,i] = median(skies[lindgen(64)*4+k,*,i])-meddark
;      skies[lindgen(64)*4+k,*,i]=skies[lindgen(64)*4+k,*,i]-meds[k,i]
    endfor
  endfor
endif

cr_flag = bytarr(XPIX,YPIX,nframes)+1B
for i=0,255 do for j = 0,255 do begin
  a = sort(skies[i,j,0:nsky-1])
  cr_flag(i,j,a[nsky-n_cr_init:nsky-1]) = 0B
endfor
here = where(skies lt 0,nhere)
if (nhere ge 1) then cr_flag[here] = 0B
;for i=0,nframes-1 do cr_flag[*,*,i] = cr_flag[*,*,i]*(badpix ne 32)

sm_skies = skies
kern = exp(-0.5*dist(7)^2/1.4^2)
kern = shift(kern/total(kern),3,3)
for i=0,nframes-1 do sm_skies[*,*,i] = convol(skies[*,*,i],kern,/edge_wrap,/center)-skies[*,*,i]
;for i=0,nframes-1 do sm_skies[*,*,i] = smooth(skies[*,*,i],7,/edge)-skies[*,*,i]
weights = 1./(abs(skies)+abs(sm_skies*4)+0.1)^2
sm_skies = 0

case ndark of
     0: offmap = fltarr(XPIX,YPIX)
     1: offmap = skies[*,*,nframes-1]
  else: offmap = total(skies[*,*,nsky:*],3)/ndark
endcase
gainmap = fltarr(XPIX,YPIX)+1
;;for i=0,255 do for j=0,255 do gainmap[i,j] = $
;;  (median(skies[i,j,*])-offmap[i,j])/(medsky-meddark)
case (procedure) of
  'gain_off_quad' : quadmap = fltarr(2,2,nframes)
  'gain_off_4ff'  : quadmap = fltarr(4,nframes)
  'gain_off_1ff'  : quadmap = fltarr(nframes)
  'off_only_1ff'  : quadmap = fltarr(nframes)
        else      : quadmap = [0]
endcase

; If frame medians differ mostly because of the variable offset use this line
; Don't use it if medians differ because of intrinsic sky brightness
;quadmap = meds

if (nh ge 1) then weights[hnan] = weights[hnan]/1e8
weights0 = weights ; save original weights
scale = scale0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Run the procedure with n_iter iterations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



umap = -1
psimap = -1
ipsimap = -1
for i=0,n_iter-1 do begin
  status = execute(procedure+','+ $
    'skies,astr,astr0,map_no,weights,n_terms,'+$
    'gainmap,offmap,quadmap,skymap,sigmadet,sigmasky,'+$
    'cr_flag,dark_flag,n_sig=n_sig_cr,scale=scale,'+$
    'pix_no=pix_no,xxs=xxs,yys=yys,xo=xo,yo=yo,'+$
    'umap=umap,psimap=psimap,ipsimap=ipsimap')
  scale = 1.00  ; use scale=1 after first iteration, or else it's scale^n_iter
  if (i eq n_iter-2) then umap = umap0
  if (i eq n_iter-2) then psimap = psimap0
  if (i eq n_iter-2) then ipsimap = ipsimap0
  if (!d.name eq 'X') then begin
    tv,bytscl(gainmap,0,2)
    tv,bytscl(offmap,-50,100),256,0
    window,4*(i+1),ys=640,xs=1024
    tv,bytscl(skymap.(0),0.9*median(skies[*,*,0]),1.2*median(skies[*,*,0]))
    wset,0
  endif

; Ideally the new weights calculated by the least-squares code should be
; used in succesive iterations. But in some past cases, I've found that
; reusing the initial weights works better.
  if (keyword_set(use_old_weights)) then  weights = weights0

endfor



; remake skymap with final weights
spix = sort(pix_no)
skymap = make_maps(skies,pix_no,spix,xxs,yys,$
  gainmap,offmap,quadmap=quadmap,weights*cr_flag)



; make the coverage map
dummy = make_maps(skies*0+1,pix_no,spix,xxs,yys,$
  gainmap*0+1,offmap*0,weights gt 0,covermap)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; First step is to reformat some results to fit desired output





if (n_tags(psimap) eq 0) then begin 
  psimap = {region0: reform([psimap],1,1)}
  for i=1,n_tags(skymap)-1 do psimap = $
    create_struct(psimap,'region'+strtrim(i,2),reform([psimap.(0)],1,1))
endif
if (n_tags(ipsimap) eq 0) then begin 
  ipsimap = {region0: reform([ipsimap],1,1)}
  for i=1,n_tags(skymap)-1 do ipsimap = $
    create_struct(ipsimap,'region'+strtrim(i,2),reform([ipsimap.(0)],1,1))
endif
if (n_elements(umap) eq 1) then begin
  umap1 = [umap]
  umap2 = [umap]
endif else begin
  gsize = size(gainmap)
  umap1 = reform(umap[0:2*gsize[4]-1],gsize[1],gsize[2],2)
  qsize = size(quadmap)
  case (qsize[0]) of
       2: umap2 = reform(umap[2*gsize[4]:*],qsize[1],qsize[2])
    else: umap2 = umap[2*gsize[4]:*]
  endcase
endelse



;rotate Channels 1 and 2 calibration back to raw orientation
if (chnlnum le 2) then begin
  gainmap = rotate(gainmap,7)
  offmap = rotate(offmap,7)
  sigmadet = rotate(sigmadet,7)
  if (n_elements(umap) ne 1) then $
    for i=0,1 do umap1[*,*,i] = rotate(umap1[*,*,i],7)
endif

;reformat quadmap to fixed [4,Nframes] format for SIP's use
if (procedure eq 'gain_off_1ff') then $
  quadmap = transpose([[quadmap],[quadmap],[quadmap],[quadmap]])
if (procedure eq 'gain_off') then quadmap = fltarr(4,nframes)
if (procedure eq 'gain_only') then quadmap = fltarr(4,nframes)
if (procedure eq 'off_only_1ff') then $
  quadmap = transpose([[quadmap],[quadmap],[quadmap],[quadmap]])

;collect filenames
names = {filename:'',sky_data:1B}
names = replicate(names,nframes)
names.filename = files_used
names.sky_data = dark_flag

; Second step is to construct headers for FITS files

print,'BEGINNING PART K'

mkhdr,head_sky,skymap.(0)
mkhdr,head_gain,gainmap
mkhdr,head_quad,quadmap
sxaddpar,head_sky,'TITLE','Sky Map',' Type of data in this file'
sxaddpar,head_gain,'TITLE','Gain Map',' Type of data in this file'
sxaddpar,head_quad,'TITLE','Var. Offset Map',' Type of data in this file'
for i=0,2 do begin
  case i of 
    0: fhead = head_sky
    1: fhead = head_gain
    2: fhead = head_quad
  endcase
  sxaddpar,fhead,'EXTEND','T',' Extensions are permitted'
  sxaddhist,'  ',fhead
  sxaddhist,'Processed with the least-squares calibration procedure',fhead
  sxaddhist,'using the following input parameters',fhead

  ; ESCARGOT PARAMETERS ONLY
  sxaddpar,fhead,'NITER',n_iter,$
    'Number of iterations for procedure'
  sxaddpar,fhead,'NCRINIT',n_cr_init,$
    'Initial number of data/pixel flagged as CRs'
;  sxaddpar,fhead,'THETA0',th0,$
;    'Rot (deg) of skymap wrt. raw data coords'
  sxaddpar,fhead,'CAL_PROC',procedure,$
    'Least-squares calibration procedure used'

  ; GAIN_*.PRO PARAMETERS
  sxaddpar,fhead,'NTERMS1',n_terms[0],$
    'Minimum number of terms in series expansion'
  if (n_elements(n_terms) eq 2) then $
    sxaddpar,fhead,'NTERMS2',n_terms[1],$
      'Maximum number of terms in series expansion'
  if (n_elements(term_limits) ge 1) then $
  sxaddpar,fhead,'TERMLIM1',term_limits[0],$
    'Lower limit on term[i]/term[i-1]'
  if (n_elements(term_limits) eq 2) then $
  sxaddpar,fhead,'TERMLIM2',term_limits[1],$
    'Lower limit on term[i]/total(term[0:i])'
  sxaddpar,fhead,'NCRSIG',n_sig_cr,$
    'N*sigma limit for identifying CR hits'
  sxaddpar,fhead,'PIXSCALE',scale0,$
    'Ratio of sky pixel size / detector pixel size'
    
  ; SIP Keywords for calibration files, from last raw data header.
  if (i eq 1) then begin
    sxaddpar,fhead,'FILENAME','','sip_keyword_filename'
    sxaddpar,fhead,'SIPCFID' ,'','SAO Pipeline ID'
    sxaddpar,fhead,'ACHANID' ,chnlnum
    sxaddpar,fhead,'OBSMODE' ,sxpar(head,'A0617D00') ? "SUBARRAY" : "FULL_ARRAY"
    sxaddpar,fhead,'ORIAORID',sxpar(head,'AORKEY')
    sxaddpar,fhead,'FRAMETIM',frametime
    sxaddpar,fhead,'DATE-OBS',sxpar(head,'DATEFILE')
    sxaddpar,fhead,'HDRMODE' ,0
    sxaddpar,fhead,'PRIORITY',2
    sxaddpar,fhead,'RAWNAME' ,''
    sxaddpar,fhead,'SIPCLDIR',''
    sxaddpar,fhead,'CFILEDIR',''
    sxaddpar,fhead,'CALTYPE' ,'SFLAT   '
    sxaddpar,fhead,'ADCEID'  ,0
    sxaddpar,fhead,'ANUMREPS',1
    sxaddpar,fhead,'CBEGTIME','2003-09-02T22:23:00'
    sxaddpar,fhead,'CENDTIME','2010-12-31T00:00:00'
  endif

  case i of 
    0: head_sky = fhead
    1: head_gain = fhead
    2: head_quad = fhead
  endcase
endfor

putast,head_sky,astr0[0]
head_sigsky = head_sky
sxaddpar,head_sigsky,'TITLE','Sigma Map',' Type of data in this file'
head_cov = head_sky
sxaddpar,head_cov,'TITLE','Coverage Map',' Type of data in this file'
head_psi = head_sky
sxaddpar,head_psi,'TITLE','Psi Map',' Type of data in this file'
sxaddpar,head_psi,'PSIMAP0',psimap0,$
  'Pixel selected for calculation of Psi (psimap)'
head_ipsi = head_sky
sxaddpar,head_ipsi,'TITLE','Psi^(-1) Map',' Type of data in this file'
sxaddpar,head_ipsi,'IPSIMAP0',ipsimap0,$
  'Flag for calculation of diag of Psi (psimap)'

head_off = head_gain
sxaddpar,head_off,'TITLE','Offset Map',' Type of data in this file'
sxaddpar,head_off,'CALTYPE' ,'SDARK   '
head_sigdet = head_off
sxaddpar,head_sigdet,'TITLE','Sigma Map',' Type of data in this file'
head_umap1 = head_off
sxaddpar,head_umap1,'TITLE','U Map (gain/offset)',' Type of data in this file'
sxaddpar,head_umap1,'UMAP0',umap0,'Pixel selected for calculation of U (umap)'
head_umap2 = head_quad
sxaddpar,head_umap2,'TITLE','U Map (var. offset)',' Type of data in this file'
sxaddpar,head_umap2,'UMAP0',umap0,'Pixel selected for calculation of U (umap)'

namehead=['']
sxaddpar,namehead,'TITLE','Input Files','Listing of input data files'
sxaddpar,namehead,'TTYPE1','FILENAME','Input file name'
sxaddpar,namehead,'TTYPE2','SKY_DATA','Sky frame = 1, Dark frame = 0'

; Third step is to write the FITS files



for i=0L,n_tags(skymap)-1 do begin
  putast,head_sky,astr0[i]
  writefits,tag_save+'.skymap.'+string(i,form='(i2.2)')+'.fits',$
    skymap.(i),head_sky
  mwrfits,names,tag_save+'.skymap.'+string(i,form='(i2.2)')+'.fits',namehead
  fits_open,tag_save+'.skymap.'+string(i,form='(i2.2)')+'.fits',fcb,/append
  fits_write,fcb,sigmasky.(i),head_sigsky
  fits_write,fcb,covermap.(i),head_cov
  fits_write,fcb,psimap.(i),head_psi
  fits_write,fcb,ipsimap.(i),head_ipsi
  fits_close,fcb
endfor

writefits,tag_save+'.skyflat.fits',gainmap,head_gain

writefits,tag_save+'.skydark.fits',offmap,head_off
fits_open,tag_save+'.skydark.fits',fcb,/append
fits_write,fcb,quadmap,head_quad
fits_close,fcb
mwrfits,names,tag_save+'.skydark.fits',namehead
fits_open,tag_save+'.skydark.fits',fcb,/append
fits_write,fcb,sigmadet,head_sigdet
fits_write,fcb,umap1,head_umap1
fits_write,fcb,umap2,head_umap2
fits_close,fcb

; Fourth, write an IDL save set for various additional parameters and results

save,file=tag_save+'.xdr',astr,astr0,map_no,weights,cr_flag,$
  meds,medsky,meddark  ; internal escargot driver info.

print,'FINISHED'

end
