clear all ; close all
HideCursor

portn = 57424 ; 

trigsON = false ; 
if trigsON
    lptwrite(portn,0) %Initiate triggers
    lptwrite(portn,0)
    lptwrite(portn,1) % start trials
end

stimtime = 3 ; % stimulus time in seconds 
resttime = 1 ; % post stimulus time 
focustime = 0 ; % minimum focus time (+ 1-2s jitter)
nparams = 8 ; 
ntrials = 4 ;

    rtrigs = 1:nparams ; rtrigs = repmat(rtrigs,[1,ntrials]) ; 
    randinds = randperm(length(rtrigs)) ; 
    stimtrigs = rtrigs(randinds) ; 
    
    HideCursor
    Screen('Preference', 'SkipSyncTests', 1);
	AssertOpenGL;
	screens=Screen('Screens');
	screenNumber=max(screens);
    
	white=WhiteIndex(screenNumber);
	black=BlackIndex(screenNumber);
	gray=(white+black)/2;
        if round(gray)==white
            gray=black;
        end
	inc=white-gray ;	
   
	% Open a double buffered fullscreen window and draw a gray background 
	% to front and back buffers:
	[w ,screenRect]=Screen('OpenWindow',screenNumber, gray) ;
    w2 = screenRect(3)/2 ; h2 = screenRect(4)/2 ; 
    % Enable alpha-blending:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
    % Show initial gray screen:
    Screen('Flip', w) ;
    f=0.05 ;
    cyclespersecond=6 ;
    angle=0 ;
    texsize=250 ; % Half-Size of the grating image.
    % Calculate parameters of the grating:
    p=ceil(1/f) ;  % pixels/cycle    
    fr=f*2*pi ;
    visiblesize=800 ; %2*texsize+1 ;
   
    orientations = 0:23:180 ; 
    [xg,yg] = meshgrid(-800:.5:800,-800:.5:800) ; 
    sg = sin(xg) ; 
    sgtex = Screen('MakeTexture',w,uint8(mat2gray(sg)*255)) ; 
    
    % red crosshair (radius 8 pixels, dark red)
    [x,y] = meshgrid(-700:700,-700:700) ; 
    xhair = zeros(size(x,1),size(x,2),4) ; 
    xhair(:,:,1) = gray ; 
    xhair(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 5) ; 
    xhairtex = Screen('MakeTexture',w,xhair) ; 
    
     % smaller brighter red crosshair
    xhairsmall = zeros(size(x,1),size(x,2),4) ;
    xhairsmall(:,:,1) = white ; 
    xhairsmall(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 3) ; 
    xhairsmalltex = Screen('MakeTexture',w,xhairsmall) ;
    
    % black crosshair (radius 8 pixels, dark red)
    xhair2 = zeros(size(x,1),size(x,2),4) ;
    xhair2(:,:,1) = black ; xhair2(:,:,2) = black ; xhair2(:,:,3) = black ; 
    xhair2(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 5) ; 
    xhair2tex = Screen('MakeTexture',w,xhair2) ;
    
    % circular aperture mask 
    circmask(:,:,1) = x.*0 + gray ; 
    circmask(:,:,2) = double(sqrt(x.^2 + y.^2) > 400)*white ;
    circmasktex = Screen('MakeTexture',w,circmask) ; 
   
    % make the polar annulus, should vary in both size and radius...
    [xp,yp] = meshgrid(-200:200,-200:200) ; 
    [theta,rho] = cart2pol(xp,yp) ; 
    % should move outwards more quickly, and also expand as it moves out.
    radinds = 5:1:200 ; clear pmasks_ext pmasks_int ; 
    maxsqr = floor(sqrt(length(radinds))) ; 
    indincr = maxsqr/length(radinds) ; 
    circinds = floor((1:indincr:maxsqr).^2) ;
    for i=1:length(circinds)
        if i < 10 
        pmasks_ext(:,:,i) = rho > radinds(circinds(i)) ; 
        pmasks_int(:,:,i) = rho > radinds(circinds(i)) + radinds(circinds(i)) ; 
        else
        pmasks_ext(:,:,i) = rho > radinds(circinds(i)) ; 
        pmasks_int(:,:,i) = rho > radinds(circinds(i)) + radinds(circinds(i)) ; 
        end
    end
    pbm = pmasks_ext + pmasks_int ; 
    disp3d(pbm==1) ;
    
    % smaller circular aperture mask
    smcircmask(:,:,1) = x.*0 + gray ; 
    smcircmask(:,:,2) = double(sqrt(x.^2 + y.^2) > 25)*white ;
    smcircmasktex = Screen('MakeTexture',w,smcircmask) ; 
    
    % Query duration of monitor refresh interval:
    ifi=Screen('GetFlipInterval', w) ;
    
    waitframes = 1;
    waitduration = waitframes * ifi ;
    
    % Translate requested speed of the grating (in cycles per second)
    % into a shift value in "pixels per frame", assuming given
    % waitduration: This is the amount of pixels to shift our "aperture" at
    % each redraw:
    shiftperframe = cyclespersecond * p * waitduration ;
    
    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp for our "WaitBlanking" emulation:
    vbl=Screen('Flip', w) ;
    
    frameRate=Screen('FrameRate',screenNumber) ;
    % Convert movieDuration in seconds to duration in frames to draw:
    stimframes = stimtime * frameRate ; 
    restframes = resttime * frameRate ; 
    focusframes = focustime * frameRate ; 
    % Use realtime priority for better timing precision:
    priorityLevel=MaxPriority(w) ;
    Priority(priorityLevel) ;    
    rotincr = 365/stimframes ; 
    rots = round(0:rotincr:365) ;     
    tic    
    currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
    for presentation = 1:size(presentorder,2) 
        restcount = 1 ; 
        focuscount = 1 ;
        stimcount = 1 ;
        poststimcount = 1 ;
        jitter = rand*frameRate*4 ; %frameRate/2 + rand*frameRate ; % jitter .5-1.5s long
        DONE = false ; 
        trigsent = false ; 
        while ~DONE             
            if currentstate == 0 % present the black crosshair (rest period)
                restcount = restcount + 1 ;   
                Screen('DrawTexture',w,xhair2tex) ;  
                if restcount > restframes
                    currentstate = 1 ; 
                end
            end         
            if currentstate == 1 % present the dark red and bright red crosshairs (focus period)
                focuscount = focuscount + 1 ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
                if focuscount > focusframes + jitter
                    if trigsON
                        lptwrite(portn,stimtrigs(presentation)) ;
                        WaitSecs(0.004) ;
                        lptwrite(portn,0) ; 
                    end
                    currentstate = 2 ;
                end
            end            
            if currentstate == 2 % present the stimulus               
                xoffset = mod(stimcount*shiftperframe/2,900);  % Shift the grating by "shiftperframe" pixels per frame:
                srcRect=[xoffset 0 xoffset + visiblesize visiblesize];           
                Screen('DrawTexture',w,sgtex,srcRect,[],orientations(stimtrigs(presentation))) ;   
                Screen('DrawTexture',w,masktex(ceil(stimcount/4))) ;   
                %Screen('DrawTexture',w,quadtex(trig1s(presentation)),[],[],rots(stimcount)) ; 
                %Screen('DrawTexture',w,circmasktex) ;                 
                Screen('DrawTexture',w,xhairtex) ; 
                Screen('DrawTexture',w,xhairsmalltex) ;  
                stimcount = stimcount + 1 ; 
                if stimcount > stimframes 
                    currentstate = 3 ;                
                end
            end              
            if currentstate == 3 % present the post stimulus fixation
                poststimcount = poststimcount + 1 ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
                if poststimcount > frameRate/2 ;% .5 seconds post stimulus fixation
                   currentstate = 0 ; % restart the loop  
                   DONE = true ;
                end    
            end         
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.           
            if KbCheck% Abort demo if any key is pressed:
               Screen('CloseAll') ;
               break
            end         
        end            
    end
    %The same commands wich close onscreen and offscreen windows also close
	%textures.
	Screen('CloseAll') ;
    toc