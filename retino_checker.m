clear all ; close all

%%%% designing a 2d parameter space randomization
%%%% you want the quadrants to repeat linearly ie 1,2,3,4,1,2,3,4,1,2,3,4
%%% but within each you want to permute the orientations






%cd c:/mscripts/finalstims/

stimtime = 30 ; % stimulus time in seconds
resttime = .75 ; % post stimulus time
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

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
    visiblesize=3000 ; %2*texsize+1 ;
    
    % generate the sequence in which the stimulus will be presented
    % (pseudorandom)
    nparams = 4 ; 
    ntrials = 2 ; 
    params = 1:nparams ; 
    allparams = repmat(params,[1,ntrials]) ; 
    perms = randperm(length(allparams)) ; 
    randparams = allparams(perms) ; 
   
    
     % create the checkerboard
    [x,y] = meshgrid(-800:800,-800:800) ; 
    [theta,rho] = cart2pol(x,y) ; 
    radincrs = 10:10:200 ; radincrs = cumsum(radincrs) ; 
    evods = mod(1:length(radincrs),2) ; clear rhos ; 
    for i=1:length(radincrs)
        if i==1 ; 
            rhos(:,:,i) = (rho <= radincrs(i)).*evods(i) ; 
        else
            rhos(:,:,i) = (rho <= radincrs(i) & rho > radincrs(i-1)).*evods(i) ; 
        end
    end

    thetaincrs = -pi:.31:pi ; 
    evods = mod(1:length(thetaincrs),2) ; 
    clear thetas ; 
    for i=1:length(thetaincrs)
        if i==1
            thetas(:,:,i) = (theta <= thetaincrs(i)).*evods(i) ;
        else
            thetas(:,:,i) = (theta <= thetaincrs(i) & theta > thetaincrs(i-1)).*evods(i) ; 
        end
    end

    srhos = sum(rhos,3) ; 
    sthetas = sum(thetas,3) ; 
    sthetas(srhos==1) = ~sthetas(srhos==1) ; 
    sthetas = uint8(mat2gray(sthetas)*255) ; 

    sthetas2 = sum(thetas,3) ; 
    sthetas2(srhos==0) = ~sthetas(srhos==0) ; 
    sthetas2 = uint8(mat2gray(sthetas2)*255) ; 

    checktex1 = Screen('MakeTexture',w,sthetas) ; 
    checktex2 = Screen('MakeTexture',w,sthetas2) ; 
   
    % put the stimuli texture IDs in an array for easy indexing
    trigs = [1,2,3,4] ;
    
    blackb = zeros(3000,3000) ; 
    whiteb = ones(3000,3000)*255 ; 
    blacktex = Screen('MakeTexture',w,blackb) ;
    whitetex = Screen('MakeTexture',w,whiteb) ;

    % red crosshair (radius 8 pixels, dark red)
    [x,y] = meshgrid(-600:600,-600:600) ; 
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
        
    fnbhd = 15 ; 
    % quadrants 
    [x,y] = meshgrid(-1200:1200,-1200:1200) ; 
    botleft(:,:,2) = (1-((x<0) .* (y>0)))*255 ; botleft(:,:,1) = x.*0 + gray ; 
    botleft(:,:,2) = imfilter(squeeze(botleft(:,:,2)),fspecial('gaussian',fnbhd,fnbhd)) ; 
    quadtex(1) = Screen('MakeTexture',w,botleft) ; 
    
    topleft(:,:,2) = (1-((x<0) .* (y<0)))*255 ; topleft(:,:,1) = x.*0 + gray ; 
    topleft(:,:,2) = imfilter(squeeze(topleft(:,:,2)),fspecial('gaussian',fnbhd,fnbhd)) ; 
    quadtex(2) = Screen('MakeTexture',w,topleft) ; 
    
    topright(:,:,2) = (1-((x>0) .* (y<0)))*255 ; topright(:,:,1) = x.*0 + gray ; 
    topright(:,:,2) = imfilter(squeeze(topright(:,:,2)),fspecial('gaussian',fnbhd,fnbhd)) ; 
    quadtex(3) = Screen('MakeTexture',w,topright) ; 
    
    botright(:,:,2) = (1-((x>0) .* (y>0)))*255 ; botright(:,:,1) = x.*0 + gray ; 
    botright(:,:,2) = imfilter(squeeze(botright(:,:,2)),fspecial('gaussian',fnbhd,fnbhd)) ; 
    quadtex(4) = Screen('MakeTexture',w,botright) ; 
    
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
    
    hzs = 8 ; 
    for h=1:length(hzs) ; 
        stimarr = zeros(1,stimframes) ; 
        hz = hzs(h) ; fs = frameRate/hz ; 
        icount = 0 ; 
        for i=1:fs:length(stimarr)
            if mod(icount,2) == 0 
                stimarr(i:i+fs) = 1 ; 
            end
            icount = icount + 1 ; 
        end
        stimarrs(h,:) = stimarr ; 
    end
    
    rotincr = 365/stimframes ; 
    rots = round(0:rotincr:365) ;     
        
    tic
    currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
    for presentation = 1:length(randparams)
        restcount = 1 ; 
        focuscount = 1 ;
        stimcount = 1 ;
        poststimcount = 1 ;
        jitter = frameRate/3 + rand*frameRate ; %frameRate/2 + rand*frameRate ; % jitter .5-1.5s long
        DONE = false ; 
        trigsent = false ; 
        stimarr = stimarrs(1,:) ; 
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
                    currentstate = 2 ;
                end
            end            
            if currentstate == 2 % present the stimulus               
                xoffset = 0 ;%mod(stimcount*shiftperframe,360);  % Shift the grating by "shiftperframe" pixels per frame:
                srcRect=[xoffset 0 xoffset + visiblesize visiblesize];    
                if stimarr(stimcount) == 0
                    Screen('DrawTexture',w,checktex1,[],[]) ;                
                else
                    Screen('DrawTexture',w,checktex2,[],[]) ;
                end
                Screen('DrawTexture',w,quadtex(randparams(presentation)),[],[],rots(stimcount)) ; 
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
              %  Screen('DrawTexture',w,circmasktex) ;                                        
                stimcount = stimcount + 1 ; 
                if stimcount > stimframes 
                    currentstate = 3 ;                
                end
            end              
            if currentstate == 3 % present the post stimulus fixation
                poststimcount = poststimcount + 1 ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
                if poststimcount > frameRate/2 % .5 seconds post stimulus fixation
                   currentstate = 0 ; % restart the loop  
                   DONE = true ;
                end    
            end         
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.           
            if KbCheck% Abort demo if any key is pressed:
               Screen('CloseAll');
               break
            end         
        end            
    end
    %The same commands wich close onscreen and offscreen windows also close
	%textures.
	Screen('CloseAll');
    toc