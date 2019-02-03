pro de_drift_2d

;+
;NAME:
;     de_drift_2d
;PURPOSE:
;     Macro to run the dedrift procedure of Eric Weeks quickly and edit parameters in a clear way
;CALLING SEQUENCE:
;    de_drift_2d
;INPUT:
;     Input will be the tracked data file *_tracked.gdf
;OUTPUT:
;     output is a x and y drift corrected tracked data file *_tracked_corr.gdf
;AUTHORS AND HISTORY:
;     Janne-Mieke Meijer, Lund University, Feb 2016
;     Maxime Bergman & Janne-Mieke Meijer, Lund University, May 2016
;-


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;INPUT PARAMETERS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Define the coordinate file you want to analyse - outputfiles are put in same folder
  filename   = 'D:\IDL\ionic_microgel_project\blub\MC17_6wt_02_clust0_mc14000_corrxy' 
  ; example for PC: 'E:\datafolder\Sample_T_low_01' 
  ; example for MAC: '/Users/Max/Documents/EXP4/S5_T15_4_01' 
  

;  sm_size = 10  ; size for smoothing
;  good1   = 20   ; maximum tracklength you want plotted to check dedrifting script
;
;;Below this point no more user input!!
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;OPENING ALL WINDOWS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;window,2,title="window 2: average drift in X and Y + smoothing function",xpos=2520,ypos=1,xsize=512, ysize = 512
;window,3,title="window 3: particle tracks before dedrift correction", xpos=2000,ypos=520,xsize=512, ysize = 512
;window,4,title="window 4: particle tracks after dedrift correction", xpos=2520,ypos=520,xsize=512, ysize = 512
;
;!P.position=[0.1,0.1,0.9,0.9]
;  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;READ IN DATASET;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  file = strcompress(filename, /remove_all)
;  
;  tra=read_gdf(strcompress(file + '_tracked.gdf', /remove_all))
;  
;  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  ;;;;;;;;;;;;;;;DETERMINE DRIFT AND DEDRIFTC;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  
;  ;determine average motion in X and Y
;  wset, 2
;  mot=motion(tra)
;  wset, 2
;  trb=rm_motion(tra,mot,smooth=sm_size)
;  
;  ;plot tracks before and after dedrift correction
;  wset, 3
;  plottr,tra,good=good1
;  wset, 4
;  plottr,trb,good=good1
;  
;  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  ;;;;;;;;;;;;;;SAVE DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  
;  outputfile = strcompress(file + '_tracked_corr.gdf', /remove_all)
;  
;  write_gdf,trb,outputfile
;
;  
;
;print, 'Performed track corrections! :D" 

end
