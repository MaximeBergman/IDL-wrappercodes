pro track_2d

;+
;NAME:
;     track_2d
;PURPOSE:
;     Macro to run the tracking procedure of Eric Weeks quickly and edit parameters
;CALLING SEQUENCE:
;     track_2d
;INPUT:
;      Input will be the coordinate *_crd.gdf file created by findcoordinates_2d
;      filename witouth the _crd.gdf extension.  
;OUTPUT:
;     2 data files '*_track.dat' and '*_track.gdf' containing for each particle 
;     1. x(um), 2. y(um), 3. area, 4. radius, 5. eccentricity, 6. framenumber 7. particle ID 
;     the coordinates are now sorted on their particle ID. 
;AUTHORS AND HISTORY:
;     Janne-Mieke Meijer, Lund University, June 2015
;     Maxime Bergman & Janne-Mieke Meijer, Lund University, May 2016
;-

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;INPUT PARAMETERS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Define the coordinate file you want to analyse - outputfiles are put in same folder
  filename   = 'E:/IDLruns/KS12PM_05wt_1112_tscan1' 
  ; example for PC: 'E:\datafolder\Sample_T_low_01' 
  ; example for MAC: '/Users/Max/Documents/EXP4/S5_T15_4_01' 
 
  maxdist         = 0.8         ;an estimate of the maximum distance that a particle would move in a single time interval;  if graph 4 goes to zero before the maxdist value than your doing great
  goodenough      = 2           ;set this keyword to eliminate all trajectories with fewer than goodenough valid positions.  This is useful for eliminating very short, mostly 'lost' trajectories
                                ;   due to blinking 'noise' particles in the data stream.
  memory          = 0           ;this is the number of time steps that a particle can be 'lost' and then recovered again.  If the particle reappears  after this number of frames has elapsed, it will be
                                ;   tracked as a new particle. The default setting is zero.this is useful if particles occasionally 'drop out' of the data.
  quiet           =  'no'       ;/quiet               ;use /quiet to not print any messages 
  save_parameters = 'yes'       ;if set to 'yes', inputparameters are saved in the same crd_info.dat file as the inputparameters for the coordinates
  
  ;optional parameters  
  ;savetracklength = 'no'       ;allows you to save the length of each track of a particle
  ;verbose         = 1          ;define with a variable to get more informational messages, place behind toggle command to turn off


  ;Below this point no more user input!!

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;READ IN DATASET;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;
  ; data in a gdf file? then use this reading string  
  inlees = strcompress(filename + '_crd.gdf', /remove_all)
  data = read_gdf(inlees)
  tres = data
  
  ;; data in a text file? then use this reading string
  ;  inlees = strcompress(filename + '_crd.dat', /remove_all)
  ;  data = read_ASCII(inlees)
  ;  tres = data.field1
  
  ;;data in a .dat file? then use this reading string
  ;     ;reading the number of particles from the seperate file
  ;     npt = fltarr(1,1)
  ;     close,1
  ;     openr,1, Strcompress(filename +'_npt.dat', /remove_all)
  ;     readf,1, npt
  ;     close,1
  ;  
  ;     header = StrArr(1)
  ;  
  ;     ;reading data
  ;     data = fltarr(6,npt)
  ;     inlees = strcompress(filename + '_crd.dat', /remove_all)
  ;     close, 1
  ;     openr, 1, inlees
  ;     readf, 1, data
  ;     close, 1


  ;Below this point no more user input!!


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;INITIALISE ARRAYS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tres = data
  npt = n_elements(data(0,*))
  maxfr = max(data(5,*))
  xmax = max(data(0,*))
  ymax = max(data(1,*))
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;OPENING ALL WINDOWS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  ; Opening the proper windows and arranging them nicely on the desktop
  window,8,xpos=1,ypos=1,xsize=512, ysize = 512
  window,9,xpos=540,ypos=1,xsize=512, ysize = 512 


  ;defining the placement of graphs in the windows
  !P.position=[0.15,0.15,0.85,0.85]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;TRACKING PARTICLES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;run the tracking code
  print, 'starting tracking code'
  t = track(tres,maxdist,goodenough=goodenough,memory=memory, verbose=verbose)
  print, 'done tracking'
  wset,9
  plottr,t,goodenough=goodenough;, x=[0,xmax],y=[0,ymax]
  
   ; plot the histogram of the the displacements between succesive frames
  wset,8
  blub=plot_hist1(mkpdf(t,1),/log)
  ;if this graph goes to zero before the maxdist value than your doing great
 
  ;make result array
  tracked = t(*,1:*)  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;SAVE DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
  outputfile = strcompress(filename + '_tracked', /remove_all)

  write_gdf, tracked, strcompress(outputfile+'.gdf',/remove_all)
  
  close,1
  openw,1, strcompress(outputfile+'.dat',/remove_all), width = 250
  
  printf,1, tracked
  close,1

;;in case you want to save the histogram of the the displacements between succesive frames
;  openw,1, strcompress(outputfile+'_hist.dat',/remove_all), width = 250
;  printf,1, blub
;  close,1
  
  if save_parameters eq 'yes' then begin
    close,1
    openw,1, filename+ '_track_info.txt' 
    printf,1,'tracking settings'
    printf,1,'maxdist                =', maxdist;
    printf,1,'goodenough             =', goodenough;
    printf,1,'memory                 =', memory;
    close,1
  endif
  
PRINT, 'done'

end


