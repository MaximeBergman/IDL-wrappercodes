pro msd_2d

;+
;NAME: 
;     msd_2d
;PURPOSE:
;     Easy to run code to make an output of the mean square displacement program (msd.pro)
;     written in IDL by John C. Crocker and Eric R. Weeks (http://www.physics.emory.edu/faculty/weeks//idl/)
;CALLING SEQUENCE:
;     msd_2d
;INPUT:
;     tracked coordinate file '*_tracked.gdf' obtained from the track2d.pro macro
;     this contained 7 collumns with:
;     1. x 2. y 3. eccentricity 4. radius 5. area 6. frame number 7. particle ID 
;     the particle ID is the track identifier. 
;OUTPUT:
;     the outpuf file contains the determined values 
;     result is a (7,n_times) array of the form
;     [time,<x>,<y>,<x^2>,<y^2>,<x^2 + y^2>,<N>], where N is the
;     approx. independent number of data points in the average.
;PARAMETERS TO BE SET:
;     'maxtime'  If
;      this is not set, MSD computes out in dt until fewer than
;     1000 independent measurements are being averaged.  Setting
;     'maxtime' or 'mydts' overrides this limitation, of course.
;     'outfile' will save the result to a gdf file named 'outfile', which
;     is useful for batch file operation.
;     'micperpix' is set to the pixel size to convert the input data
;     file to microns if that has not been done already.  If it is
;     single number, the magnification is assumed to be isotropic,
;     otherwise, provide a dim-dimensional vector.
;     'timestep' is the number of seconds between frames (e.g. 1/60.)
;     'dim' is the spatial dimensionality of the data, defaults to 2.
;     'mydts' allows the user to supply a vector of (frame number) dt's
;     if he/she doesn't want my log-spaced points in time.
;     'erode' drops all data within 'erode' time steps of an end or
;      gap in a trajectory, to eliminate mistaken track contributions.
;      Very useful-- usually erode=2 or 3 works pretty well.
;      'minN' sets minimum number of independent measurements; only useful
;      maxtime is not used, e.g. in batch processing of plates
;      noplot allows you to NOT plot anything at the end. Added by GLH 02/18/2011
;      keepdrift does not remove <x>^2  
; NOTES:
;       The variances are presumably what you are interested in.  The
;       means are there merely to give you an idea of the average 'drift'
;       in the data (and because <x>^2 has to be subtracted from the
;       variances).  'N' can be used for error estimation: the error
;       in the variance < x(dt)^2 > ~ < x(dt)^2 >/sqrt( N(dt) ).
;-


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;INPUT PARAMETERS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Define the coordinate file you want to analyse - outputfiles are put in same folder
  filename   = 'E:/IDLruns/KS12PM_05wt_1112_tscan1' 
  ; example for PC: 'E:\datafolder\Sample_T_low_01' 
  ; example for MAC: '/Users/Max/Documents/EXP4/S5_T15_4_01' 
  
   dedrifted        = 'no'         ; set to 'yes' if you have corrected your tracked data for drifting. Outputfile will be '_msd_corr.gdf'
   maxframe         =  20           ; is the maximum 'dt' to calculate, in frames (however program will not return a value when certain number for criteria is not met
   save_parameters  = 'yes'         ; if set to 'yes', inputparameters are saved in the same crd_info.dat file as the inputparameters for the coordinates

  ; additional parameters
  ;timestep         = 0.036         ; 'timestep' is the number of seconds between frames (e.g. 1/60.) turn on in case you do not have a *_tstep.dat file, this is normally automatically peduced by findcoord_2d

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;READ IN DATASET;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;read in tracked data
  if dedrifted eq 'yes' then $
  filetoread = strcompress(filename + '_tracked_corr.gdf', /remove_all) $
  else filetoread = strcompress(filename + '_tracked.gdf', /remove_all)
  
  tracked = read_gdf(filetoread)
  
  ;read in time step
  if not(keyword_set(timestep)) then begin
    time = fltarr(1,1)
       close,1
       openr,1, Strcompress(filename +'_tstep.dat', /remove_all)
       readf,1, time
       close,1
    
    timestep = time(0,0)
    
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;OPENING ALL WINDOWS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
 window,10,title="Mean square displacement",xpos=1200,ypos=1,xsize=512, ysize = 512

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;RUN MSD;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  msdout = msd(tracked, maxtime=maxframe*timestep, timestep=timestep)
  ;print, msdout
  
  wset,10
  plot,msdout(0,*),msdout(5,*),$
    psym=circ(),  $ 
    /xl, /yl, $
    XTITLE='time (s)', YTITLE='msd (um^2)'


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;SAVE DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; make an output file
  if dedrifted eq 'yes' then $
    outputfile = strcompress(filename + '_msd_corr', /remove_all) $
  else outputfile = strcompress(filename + '_msd', /remove_all)


  write_gdf, msdout, strcompress(outputfile+'.gdf',/remove_all)
  
  close,1
  openw,1, strcompress(outputfile+'.dat',/remove_all), width = 250
  printf,1, 'time','<x>','<y>','<x^2>','<y^2>','<x^2+y^2>','<N>',  FORMAT ='(7(A13))'
  printf,1, msdout  
  close,1
  
  
  if save_parameters eq 'yes' then begin
    close,1
    openw,1, filename+ '_msd_info.txt' 
    printf,1,'msd settings'
    printf,1,'dedrifted            =', dedrifted
    printf,1,'timestep             =', timestep;
    printf,1,'maxframe             =', maxframe;
    close,1
  endif


end





