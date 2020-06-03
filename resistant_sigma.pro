;+
;NAME:
; resistant_sigma
;
;PURPOSE:
; An outlier-resistant determination of the mean and its standard deviation.
; It trims away outliers using the median and the median absolute deviation.
;
;CALLING SEQUENCE:
; sigma = resistant_sigma(data, cut, mean, wts)
;
;INPUT ARGUMENTS:
; data  = data to average
; cut   = Data more than this number of standard deviation from the
;             median is ignored. Suggested values: 2 and up.
;             Lower limit of cut > 1.0 is imposed.
;
;OUTPUT ARGUMENTS:
; mean  = the mean
; sigma = the standard deviation of the mean
; wts   = 0 if rejected, 1 if kept
;
;AUTHOR:   H. Freudenreich, STX, 1989; Second iteration added 5/91.
;          (Original version called resistant_mean.pro)
;MODIFIED: BA Franz, ARC, Changed from procedure to function and replaced
;          number rejected tracking with wts tracking. 6/92
;          R. Arendt. 3/00 : adpated resistant_mean to give a robust estimate
;          of the RMS sigma.
;          R. Arendt. 1/02 : corrected bug in conditional expression for 
;          application of sigma truncation correction, and added HTF's new
;          3rd order correction for sigma truncation correction. 
;          Streamlined code slightly.
;          R. Arendt. 3/02: Removed dependence on avg.pro
;-
function resistant_sigma,data,cut,mean,wts

;
; Must have at least 2 points
;
npts = n_elements(data)
if (npts le 1) then begin
    mean  = data(0)
    wts   = [0]
    return,0.0
endif 

;
; Limit on cut for truncation compensation of standard deviation
;
sc = cut>1.0

;
; Estimate standard deviation from median absolute deviation
;
mean        = median(data,/even)
absdev      = abs(data-mean)
sigma       = median(absdev,/even)/.6745
if (sigma lt 1.0e-24) then sigma = total(absdev)/n_elements(absdev)/.8

;
; Iterate twice
;
for iter=0,1 do begin
  wts         = absdev le sc*sigma
  goodpts     = data(where(wts,num_good))
  mean        = total(goodpts)/num_good
  
  ;
  ; Calculate new standard deviation and compensate for truncation
  ;
  sigma     = sqrt( total((goodpts-mean)^2)/num_good )    ; ?? num_good-1 ??
  if (sc le 3.50) then $
    sigma = sigma/(-0.15405+sc*(0.90723+sc*(-0.23584+sc*0.020142)))
endfor

return,sigma
end
