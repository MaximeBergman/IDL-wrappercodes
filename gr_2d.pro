pro gr_2d


;+
;NAME:
;     gr_2d
;PURPOSE:
;     Macro to run g(r) programme of Eric Weeks quickly and edit parameters
;     The procedure uses the g(r) subroutine ericgr2d.pro
;CALLING SEQUENCE:
;     gr_2d
;INPUT:
;     The input is the coordinate file created with findcoordinates_2d
;     This contains the particle information in the follwing collumns
;      1. x, 2. y, 3. area, 4. feature radius, 5. eccentricity, 6. framenumber 
;OUTPUT:
;     data file '*_gr2d.dat' containing two columns
;     1. r (um) 2. g(r) 
;AUTHORS AND HISTORY:
;     Janne-Mieke Meijer, Lund University, June 2015
;     Maxime Bergman & Janne-Mieke Meijer, Lund University, May 2016
; 
;-



 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;INPUT PARAMETERS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Define the coordinate file you want to analyse - outputfiles are put in same folder
;without '_crd.gdf' 
 filename   = 'E:/IDLruns/KS12PM_05wt_1112_tscan1' 
 ; example for PC: 'E:\datafolder\Sample_T_low_01' 
 ; example for MAC: '/Users/Max/Documents/EXP4/S5_T15_4_01' 
 
 
tracked   = 'no'         ; if your data is tracked, than fill in 'yes', if not than type in anything else
rmin      = 0.0          ; minumum radius in micrometer to take into acount
rmax      = 5.1            ; maximum radius in micrometer to take into acount
deltar    = 0.03            ; stepsize in the g(r) plot, this should  be no less than roughly 1/5 of your pixel size in micrometer

save_parameters = 'yes'   ;if set to 'yes', inputparameters are saved in the same crd_info.dat file as the inputparameters for the coordinates

;Below this point no more user input!!

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;READ IN DATASET;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;

;;;;from .gdf file;;;;
data=read_gdf(filename + '_crd.gdf')

;defining the windows 
window,7,xpos=800,ypos=1,xsize=1000, ysize = 500


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;CREATE GR;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;routine if the data is tracked data
if tracked eq 'yes' then begin
gr = ericgr2d(data,/track,rmin=rmin,rmax=rmax,deltar=deltar)
plot,gr(0,*),gr(1,*)

;routine if the data is only coordinates and frame numbers
endif else gr = ericgr2d(data,rmin=rmin,rmax=rmax,deltar=deltar)
plot,gr(0,*),gr(1,*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;SAVE DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

close,1
openw,1,strcompress(filename + '_gr2d.dat', /remove_all)
printf,1,'r(um)', 'g(r)', FORMAT = '(2(A13))'
printf,1,gr
close,1

if save_parameters eq 'yes' then begin
  close,1
  openw,1, filename + '_gr_info.txt' 
  printf,1,'2d gr settings:'
  printf,1,'rmin             =', rmin;
  printf,1,'rmax             =', rmax;
  printf,1,'deltar           =', deltar;
  close,1
endif


end