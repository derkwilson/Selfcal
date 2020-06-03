function pixel_numbers, astrom, astrom0, sky_size, map_no, dark_flag, $
    xn, yn, xo, yo, xpix, ypix, $
    scale=scale, distort_ch=distort_ch, distort_coeff=distort_coeff,$
                        name_pix_no=name_pix_no

;
; This function generates sky pixels numbers for all data in a set of images.
;
;RETURNED RESULT:
; (pixels)  = sky pixel positions in a [Px,Py,M] array
;INPUT:
; astrom    = array of M of astrometry structures extracted 
;             from FITS headers, or constructed from other info
; astrom0   = the output map(s) astrometry
; sky_size  = size(skies), which contains the dimensions and number
;             of images for which pixel numbers are needed
; map_no    = array of M indicies indicating appropriate map number per frame
; dark_flag = array of M flags indicating dark frames (0=dark,1=sky)
;OUTPUT:
; xn        = x dimension of the array in which the pixel numbers are embedded
; yn        = y dimension of the array in which the pixel numbers are embedded
; xo        = x coordinate of the lower left corner
; yo        = y coordinate of the lower left corner
; xpix	    = exact x-pixel coordinates
; ypix	    = exact y-pixel coordinates
;OPTIONAL INPUT:
; scale      = ratio of output to input pixel sizes (<1 -> subpixelization)
;              WARNING: scale will alter the values of astrom[0].cd
; distort_ch = channel number (1-4) for IRAC distortion correction if specified
;             (OBSOLETE. Distortion parameters should be in the FITS astrometry)
; distort_coeff = distortion correction if specified: [2,3] or [2,7] array
;               This overrides the default coefficients invoked by distort_ch
;             (OBSOLETE. Distortion parameters should be in the FITS astrometry)
;
; This function returns pixel numbers for each observation of
; a detector, whose dimensions are given in sky_size, as
; embedded in a pixel grid that is large enough to contain all
; of the dithered images.
;
; This function uses W. Landsman's IDL astrometry library procedures at 
; http://idlastro.gsfc.nasa.gov
;
; The returned pixel numbers are not absolute, i.e. do not
; correspond to a fixed sky location, as do DIRBE pixel numbers.
;
; To rotate an entire dataset, transform the astrometry of the target 
; header before you start:
; IDL> astrom[0].cd = [[cos(th),sin(th)],[-sin(th),cos(th)]]##astrom[0].cd
; Rotating and/or shifting parts of data set requires more detailed 
; spherical coordinate transformations of the individual headers. 
;
; astrom[0].crpix may be shifted by an integer number of pixels to match the 
; xn by yn map defined by the output pixels. Fractional pixels shifts between 
; the input astrometry and output pixels should no longer occur.
; astrom[0].cd is updated to include the specified scale factor
; Potential confusion may be avoided by using the desired values of astrom[0]
; from the start, rather than relying on the scale parameter to change 
; the pixels scales within the calibration code. 
;
; The lower left sky pixel of the final map is reserved for the location of 
; all dark frame observations.
;
; The xscale parameter in get_pixel_numbers is dropped here, since differing
; x and y scales should be reflected in the FITS header astrometry 
; in the CD matrix or CDELT parameters.

;exist = file_test(name_pix_no)
exist=0
if exist then begin
   restore, name_pix_no 
   return, pixels
endif else begin

   if (n_elements(scale) ne 1) then scale = 1.0

   xs = sky_size[1] & ys = sky_size[2] & zs = sky_size[3]
   xpix = fltarr(xs,ys,zs)
   ypix = fltarr(xs,ys,zs)
;pixels = lonarr(xs,ys,zs) original
   pixels = lon64arr(xs,ys,zs)  ;ketron
   x0 = double(lindgen(xs*ys) mod xs)
   y0 = double(lindgen(xs*ys) / xs)

; apply distortion in the FITS (u,v) coordinate system
   if (keyword_set(distort_ch) or (n_elements(distort_coeff) ge 6)) then begin
      x0 = x0+1-astrom[1].crpix[0]
      y0 = y0+1-astrom[1].crpix[1]
      undistort_iracuv,x0,y0,distort_ch,x0,y0,dcoeff=distort_coeff
      x0 = x0-1+astrom[1].crpix[0]
      y0 = y0-1+astrom[1].crpix[1]
   endif

   astrom0.cd = astrom0.cd*scale
   for i=0L,zs-1 do begin
      xy2ad,x0,y0,astrom[i],a,d
      ad2xy,a,d,astrom0[map_no[i]],xp,yp
      xpix[*,*,i] = xp
      ypix[*,*,i] = yp
;  print, minmax(xp)
;  print, minmax(yp)
;  print, 'i=',i
;  print, '---------'
   endfor

   nmaps = n_elements(astrom0)
;;- Original
;xo = lonarr(nmaps)
;yo = lonarr(nmaps)
;xn = lonarr(nmaps)
;yn = lonarr(nmaps)

   xo = lon64arr(nmaps)         ;ketron
   yo = lon64arr(nmaps)
   xn = lon64arr(nmaps)
   yn = lon64arr(nmaps)

   for i=0,nmaps-1 do begin
      h = where(map_no eq i)
      xo[i] = round(min(xpix[*,*,h],max=maxx))-1
      yo[i] = round(min(ypix[*,*,h],max=maxy))
      xn[i] = round(maxx)-xo[i]+1
      yn[i] = round(maxy)-yo[i]+1
      xpix[*,*,h] = xpix[*,*,h]-xo[i]
      ypix[*,*,h] = ypix[*,*,h]-yo[i]
      ; check that XPIX and YPIX > -0.5, since -0.5 would get rounded to -1 and
      ; generate negative values in PIXELS. Should be a very rare condition, but 
      ; it has occured in real data sets.
      hneg = where(xpix le -0.5,nh)
      if (nh ge 1) then xpix[hneg] = 0
      hneg = where(ypix le -0.5,nh)
      if (nh ge 1) then ypix[hneg] = 0
      pixels[*,*,h] = round(xpix[*,*,h])+xn[i]*round(ypix[*,*,h])
      if (i gt 0) then pixels[*,*,h] = pixels[*,*,h] + total(xn[0:i-1]*yn[0:i-1])
      astrom0[i].crpix = astrom0[i].crpix-[xo[i],yo[i]]
;  stop
   endfor

   here = where(dark_flag eq 0,nhere)
   if (nhere gt 0) then begin
      pixels[*,*,here] = 0
      xpix[*,*,here] = 0
      ypix[*,*,here] = 0
   endif
   
;   if keyword_set(name_pix_no) then begin
;      save, pixels, filename = name_pix_no
;      message, 'Saving ' + name_pix_no, /inf
;   endif
   return,pixels
endelse


end
