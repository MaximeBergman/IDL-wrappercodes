pro findcoord_2d

;+
;NAME: 
;     findcoord_2d
;PURPOSE: 
;     Macro for running the procedyres written in IDL by John C. Crocker and Eric R. Weeks 
;     (http://www.physics.emory.edu/faculty/weeks//idl/)
;     Removes noise and smooths frames via bandpass filter. 
;     Then finds isotropic features in dataset.
;CALLING SEQUENCE: 
;     findcoord_2d
;INPUT:
;     Confocal Microscopy images or stack processed via ImageJ
;     imagefile: name of input file '*.tif' 
;OUTPUT:
;     Coordinate file '*_crd.gdf' and '*_crd.dat' containing the features found in the following output:
;     1. x(um), 2. y(um), 3. area, 4. radius, 5. eccentricity, 6. framenumber
;     Output file with total number of particles found '*_npt.dat' 
;     Output file with timestep between frames '*_tstep.dat' 
;     During the analysis the procedure will show: 
;           Tif image being analyzed
;           Bandpass filtered image
;           Bandpass filtered image with positions of features found
;           Tif image with positions of features found
;           Histogram of distribution of decimals in pixels for each feature
;           Histogram of feature area's
;           Scatter plot of area vs radius 
;AUTHORS AND HISTORY:
;     First developed by Roel Dullen, Dirk Aarts, Utrecht University (2002-2005)
;     Followed by addaptions by Volker de Villeneuve and Jan Hilhorst, Utrecht University (2006-2012)
;     Janne-Mieke Meijer, Lund University (June 2015)
;     Maxime Bergman, Lund University (July 2015) 
;     Janne-Mieke Meijer and Maxime Bergman, Lund University (May 2016)
;-

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;INPUT DATASET;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ; Define the .tif file to analyse 
 imagefile = 'E:/IDLruns/KS12PM_05wt_1112_tscan1.tif' 
 ; example for PC: 'E:\datafolder\Sample_T_low_01.tif' 
 ; example for MAC: '/Users/Max/Documents/EXP4/S5_T15_4_01.tif' 
 
 ; Define the output name without extention, this will be added automatically for the different outputs
 outputname   = 'E:/IDLruns/KS12PM_05wt_1112_tscan1' 
 ; example for PC: 'E:\datafolder\Sample_T_low_01' 
 ; example for MAC: '/Users/Max/Documents/EXP4/S5_T15_4_01' 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;INPUT PARAMETERS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;parameters for bandpass and feature finding - ADJUST THESE ALWAYS
  noise     = 2               ; variation in sphere size (pixels)
  diam      = 27              ; diameter of the spheres (pixels)
  minarea   = 155000          ; Always set to 100-1000 the first try, than adjust based on the int.area vs feature radius plot.
  maxarea   = 100000000        ; Always set to very high the first try, than adjust based on the int vs area histogram.
  
  ;additional parameters for feature finding
  maxecc    = 0.5            ; how spherical the particles are. 0 = cirkel, 1 = line. Everything smaller than maxecc is taken for analysis
  mindiam   = 0              ; minimum diameter filter for particle cload in int.area vs feature diameter plot
  lowlimit  = 0              ; low limit of the intensity that is taken into acount in the analysis (any pixel intensity < lowlimit is made 0)
  
  
  ;parameters for selection of frames to analyse
  ibegin    = 0              ; first frame number for analysis
  iend      = 1000000              ; end frame number for analysis
  step      = 1              ; stepsize between frames
  frames    = iend-ibegin+1  ; defining the number of images that are used
  
  ;additional parameters
  invertimage = 'no'         ; invert the images, this means image intensity will be flipped using: image = 255 - image
  quiet       = 'yes'          ; if set to quiet printing of the data will occur, of set to anything else, no printing will occur
  corrx = 1.                  ; set the correction factor for x (in case of miscallibration)
  corry = 1.                  ; set the correction factor for y (in case of miscalibaration)
  
  ;parameters for saving
  saved      = 'yes'          ; want to save the tiffs with coordinates in them? saved = 'yes' they are saved; saved = 'anything else' they are not saved
  splitco    = 'no'           ; want to save the coordinates per frame? saved = 'yes' they are saved; saved = 'anything else' they are saved in a single file
  ballview   = 'no'         ; want to save the coordinates for the ballviewer program? ballviewer = 'yes'
  
  ;if you have clusters
  clustercheck  = 'no'       ; if set to 'yes', will count clusters.
  clustervalue  = 1.25       ; how many times the diameter is overlap?
  throwclusters = 'no'       ; if set to 'yes', the found clusters will be averaged.

  ;if you have bleaching
  intcorr       = 'yes'      ; set to'yes' if you observe bleaching over time
  intcorr_step  = 500       ; define the step size for a 5% correction in minarea
  intcorr_rep   = 5         ; number of times to repeat 5% correction step, after this the correction turns to 3%. For datasets in which the bleaching seems to be nonlinear over time

  ;IMPORTANT: Below this no more user input!!
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;INITIALISE ARRAYS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;additional parameters
  k = 1
  radius = diam/2 ; for fover2d plot
  sep = diam+1 ; separation between found features (pixels) DO NOT CHANGE
  totalcountduo = 0 ; DO NOT CHANGE

  ;coordinaten-output-array:
  res       = fltarr(6)
  
  ok = file_test(imagefile)
  if not (ok) then begin
  print, imagefile, ' !!File cannot be found! Check filename!!
  return
  endif
  
  
;read in imagestack
    print, 'IDL is importing the tiff stack file...(this might take some time..)'
    imagestack = readtiffstack(imagefile)
    print, 'Import is finished'

  ;for saved
  if (not keyword_set(saved)) then begin
    print, 'NOTE: output imagefiles  will not be saved'
    saved='no'
  endif
  
  
  ; determine image dimensions from the first frame wanted...
  image0 = imagestack(*,*,0)     ;ELSE2.3
  dims      = size(image0,/dimensions)
  ;print,'dimensions are',dims
  dims0     = dims(0)
  dims1     = dims(1)
  pic      = fltarr(dims0,dims1)

 
  ;getting the x,y,z size information automatically from the files
  ok = query_tiff(imagefile, info)
  if (ok) then begin
      sizex = info.resolution
      if sizex(0,0) eq 1 then message,'Error: imagefile contains no size information'
      pixxy = 1/sizex(0,0)
      xlength = pixxy * dims(0)
  endif
   
  ns = n_elements(imagestack(0,0,*))
  if ns gt 200 then rep = 25 else rep=1
  if frames lt ns then iend = iend else iend = ns-1

   
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;OPENING ALL WINDOWS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Opening the proper windows and arranging them nicely on the desktop
  window,1,xpos=1,ypos=1,xsize=dims0, ysize = dims1
  window,2,xpos=525,ypos=1,xsize=dims0, ysize = dims1
  window,3,xpos=1050,ypos=1,xsize=dims0, ysize = dims1
  window,4,xpos=1,ypos=540,xsize=256, ysize = 256
  window,5,xpos=400,ypos=540,xsize=700, ysize = 500
  window,6,xpos=1020,ypos=540,xsize=700, ysize = 500

  ;defining the placement of graphs in the windows
  !P.position=[0.1,0.1,0.9,0.9]
   
   
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;START FEATURE FINDING PER FRAME;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
  ; loop over all the frames
  for i = ibegin,iend,step do begin
  
    if intcorr eq 'yes' then begin
      j = round(i/intcorr_step)
         
      if (j gt 0 and i eq ibegin) then minarea2 = (1-(0.05*j))*minarea ;this is in case you only want to analyse a small part of the data and start in the already bleached part
      if (j gt 0 and (j-k) eq 1 and j lt intcorr_rep ) then begin
        minarea2 = (1-(0.05*j))*minarea 
        k = k + 1
      endif
      if (j gt 0 and (j-k) eq 1 and j ge intcorr_rep ) then begin
          k = k + 1
          minarea2 = (1-((0.05*(intcorr_rep -1))+(0.03*(j-intcorr_rep))))*minarea
      endif 
    endif
  
  
    if (((i+1) mod rep eq 0) and not keyword_set(quiet)) then begin
      print,'processing frame'+$
        strcompress(i+1)+' out of'+strcompress(iend)+'....'
    endif
    
    aa = imagestack(*,*,i)
    
    if invertimage eq 'yes' then aa = 255-aa

    pic(*,*) = aa(*,*)
    pica = pic
    
    wset, 1
    tvscl, pica
    
    ; apply bandpass filter
    im = bpass(aa,noise,diam)
    imbyt = bytscl(im)
    
    wset,2
    tv,imbyt
    
    ; determine the features in the image
    dd = feature(im,diam,sep,masscut=minarea,min=lowlimit)
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;FEATURE SELECTION;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    nblub=n_elements(dd(0,*))
    
    if nblub gt 1 then begin
      ww= where(dd(3,*) ge mindiam, count)
      if (count ne 0) then dd=dd(*,ww)
      
      w2max = where(dd(2,*) le maxarea)
      dd = dd(*,w2max)
      
      eccmax = where(dd(4,*) le maxecc,countecc)
      ;if eccmax(0,0) ne -1 then dd = dd(*,eccmax) else dd = dd(*,0) <-- why is this in here
      if (countecc ne 0) then dd=dd(*,eccmax) ;<-- suggested code
      
      ;CLUSTER REMOVAL ROUTINE
      ndd = n_elements(dd(0,*))
      indicestoremove=fltarr(1)
      if clustercheck eq 'yes' then begin
        print, 'particles before cluster check', ndd
        
        for n = 1, ndd-1, 1 do begin ;a loop to check for clusters
          wduo = where((dd(0,n-1)-dd(0,*))^2 + (dd(1,n-1)-dd(1,*))^2 lt (clustervalue*diam)^2 and (dd(0,n-1)-dd(0,*))^2 + (dd(1,n-1)-dd(1,*))^2 gt 0.001, countduo, COMPLEMENT=notwduo);where are the clusters located
          ;print, wduo
          if wduo(0,0) ne -1 and countduo ge 1 then begin ;if clusters are found than average the coordinates
            connectvectors = fltarr(5, countduo)
            connectvectors(0,*) = dd(0,wduo) ;xcoor van alle NNs in kolom 0
            connectvectors(1,*) = dd(1,wduo) ;ycoor van alle NNs in kolom 1
            connectvectors(2,*) = dd(2,wduo) ; area van alle NNs in kolom 2
            connectvectors(3,*) = dd(3,wduo) ; elipticity van alle NNs in kolom 3

            ;remove the cluster from features
            indicestoremove=[[indicestoremove],wduo, n-1];          
            
          endif
        endfor
        
        dd = dd[*, Where(~Histogram(indicestoremove, MIN=0, MAX=ndd), /NULL)] ;also throw away reference particle
        ndd = n_elements(dd(0,*))
        indicestoremove=fltarr(1)
        if ndd gt 1 then totalcountduo=totalcountduo+countduo
        print, 'clusters thrown out, reduced features found to',ndd
        
      endif 

      if (dd(0) ne -1) then res = [[res],[dd,fltarr(1,ndd)+i]]
    endif ;feature finding per frame if-loop
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;CREATE AND SAVE RESULT IMAGES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
    
    ;plotting the coordinate analysis on the bandpass filter and the raw data
    wset,3
    picb = pic
    ee=fover2d(picb,dd,/circle,radius=diam/2)

    wset,2
    ee=fover2d(imbyt,dd,/circle,radius=diam/2)

    if saved eq 'yes' then begin                                                                                

      filenametif = strcompress(outputname + '_out.tif', /remove_all)
      
      if frames lt ns then begin 
        outputfile_im = strcompress(outputname + '_fr' + string(ibegin) + '_fr' + string(iend) + '_out.tiff', /remove_all)
        write_tiff, outputfile_im, ee
      endif else begin    
        if i eq ibegin then write_tiff, filenametif ,ee $
        else write_tiff, filenametif ,ee, /APPEND
      endelse 

    endif
    
   endfor     ; end of loop over the frames 

   

   ;plot the analysis of the coordinate-output-array
   wset, 4
   plot_hist,res(0,*) mod 1            ; plot of the subpixel resolution
   wset, 5
   plot_hist, res(2,*)                 ; plot a histogram of area
   wset, 6
   plot, res(2,*), res(3,*), psym=1   ; plot a cloud of area vs radius

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;CONVERT AND CORRECT COORDINATES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  coord3d        = res(*,1:*)
  coord3d(0:1,*) = coord3d(0:1,*)*pixxy   ;pixel to micrometer  (turn on for 3d coord, turn off for psi6 gr etc.)
  coord3d(0,*)   = coord3d(0,*) * corrx   ;correction for x
  coord3d(1,*)   = coord3d(1,*) * corry   ;correction for y
  coord          = coord3d


  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;SAVING COORDINATES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; 
  ; Writing the coordinates into separate data files
  if splitco eq 'yes' then begin
    for i=ibegin,iend do begin
      w=where(coord(2,*)eq i)
      coordspl = coord(0:1,w)
      if (i lt 10)   then filename = strcompress(outputname + '_000' + string(i) + 'coord.dat', /remove_all) $
      else if (i lt 100)  then filename = strcompress(outputname +  '_00' + string(i) + 'coord.dat', /remove_all) $
      else if (i lt 1000) then filename = strcompress(outputname +   '_0' + string(i) + 'coord.dat', /remove_all) $
      else                     filename = strcompress(outputname +    '_' + string(i) + 'coord.dat', /remove_all)
      
      close,1
      openw,1,filename, width = 100
      printf,1,coordspl
      close,1
      
    endfor
  endif
  
  ;save as .dat
  if frames lt ns then $
    outputfile = strcompress(outputname + '_fr' + string(ibegin) + '_' + string(iend) + '_crd.dat', /remove_all) $
  else outputfile = strcompress(outputname + '_crd.dat', /remove_all)
  
  close,1
  openw,1, outputfile, width = 250
  printf,1, 'x(um)','y(um)', 'area', 'radius', 'eccentricity', 'fr_number' , FORMAT ='(6(A13))'
  printf,1, coord
  close,1
    
  ;save as .gdf
  if frames lt ns then $
    outputfile1 = strcompress(outputname + '_fr' + string(ibegin) + '_' + string(iend) + '_crd.gdf', /remove_all) $
  else outputfile1 = strcompress(outputname + '_crd.gdf', /remove_all)
  write_gdf, coord, outputfile1
  
  ;This part is used to generate the coordinate files *.xyz that can be viewed in the ballviewer
  if ballview eq 'yes' then begin
    coordxyz = [coord(0:1,*),coord(4,*)]
  
    close,1
    openw,1,strcompress(outputname + '_crd.xyz', /remove_all)
    printf,1,npt,1,0
    printf,1,coordxyz
    close,1
  endif
  
  ; Print the amount of coordinates that were found 
  npt = n_elements(coord(0,*))

  if frames lt ns then $
    outputfilenpt = strcompress(outputname + '_fr' + string(ibegin) + '_' + string(iend) + '_npt.dat', /remove_all) $
  else outputfilenpt = strcompress(outputname + '_npt.dat', /remove_all)
  
  close,1
  openw,1, outputfilenpt
  printf,1, npt
  close,1

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;SAVING PARAMETERS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  
if frames lt ns then $
  outputinfo = strcompress(outputname + '_fr' + string(ibegin) + '_' + string(iend) + '_crd_info.txt', /remove_all) $
else outputinfo = strcompress(outputname + '_crd_info.txt', /remove_all)

ok = query_tiff(imagefile, info)

if (ok) then begin                                                                                      ;IF1

  str = info.description
  var = STRSPLIT(str, STRING(10b), /EXTRACT) ; Special Characters in Regular Expressions by doing string(10b)
  var2 = var.replace('=', '   =   ')
  nvar = n_elements(var2)

  if nvar ge 5 then begin
    var3 = var2[4].split('   =   ')

    tstep = var3[1]

        if frames lt ns then $
          outputfilet = strcompress(outputname + '_fr' + string(ibegin) + '_' + string(iend) + '_tstep.dat', /remove_all) $
        else    outputfilet = strcompress(outputname + '_tstep.dat', /remove_all)

    close,1
    openw,1, outputfilet
    printf,1, tstep
    close,1
  endif

  close,1
  openw,1, outputinfo, width = 100
  printf,1, 'noise          = ', noise
  printf,1, 'diam           = ', diam 
  printf,1, 'xlength        = ', xlength 
  printf,1, 'saved          = ', saved 
  printf,1, 'splitco        = ', splitco 
  printf,1, 'lowlimit       = ', lowlimit
  printf,1, 'maxecc         = ', maxecc 
  printf,1, 'mindiam        = ', mindiam
  printf,1, 'minarea        = ', minarea
  printf,1, 'maxarea        = ', maxarea
  printf,1, 'ibegin         = ', ibegin
  printf,1, 'iend           = ', iend
  printf,1, 'step           = ', step
  printf,1, 'frames         = ', frames
  printf,1, 'invertimage    = ', invertimage
  printf,1, 'tstep          = ', tstep
  printf,1, 'clustercheck   = ', clustercheck
  printf,1, 'corrx          = ', corrx
  printf,1, 'corry          = ', corry
  printf,1, '----------------------------------------------------------------------------------------'
  printf,1, 'Info from ImageJ:'
  for i = 0,nvar-1,1 do begin
    printf,1, var2[i]
  endfor
  close,1
  
endif else begin
  close,1
  openw,1, outputinfo, width = 100
  printf,1, 'noise          = ', noise
  printf,1, 'diam           = ', diam 
  printf,1, 'xlength        = ', xlength 
  printf,1, 'saved          = ', saved 
  printf,1, 'splitco        = ', splitco 
  printf,1, 'lowlimit       = ', lowlimit
  printf,1, 'maxecc         = ', maxecc 
  printf,1, 'mindiam        = ', mindiam
  printf,1, 'minarea        = ', minarea
  printf,1, 'maxarea        = ', maxarea
  printf,1, 'ibegin         = ', ibegin
  printf,1, 'iend           = ', iend
  printf,1, 'step           = ', step
  printf,1, 'frames         = ', frames
  printf,1, 'invertimage    = ', invertimage
  printf,1, 'tstep          = ', tstep
  printf,1, 'clustercheck   = ', clustercheck
  printf,1, 'corrx          = ', corrx
  printf,1, 'corry          = ', corry
  close,1
endelse

print, 'done'

end