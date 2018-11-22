clear all ; close all
HideCursor

portn = 57424 ; 

stimtypes = {'right pink','left pinky','right thumb','left thumb','right index','left index','right toes','left toes','tongue'} ; 


trigsON = false ; 
Screen('Preference', 'DefaultFontSize', [100]);
if trigsON
    lptwrite(portn,0) %Initiate triggers
    lptwrite(portn,0)
    lptwrite(portn,1) % start trials
end
stimtime = 1 ; % stimulus time in seconds
resttime = 1 ; % post stimulus time
focustime = 1 ; % minimum focus time (+ 1-2s jitter)

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
    
    % generate the sequence in which the stimulus will be presented
    % (pseudorandom)
    nparams = 9 ; 
    npresents = nparams*3 ;
    combos = zeros(1,nparams) ; 
    presentorder = zeros(1,nparams)  ;
    presentcount = 1; 
    while sum(sum(combos)) ~= npresents
        rq = ceil(rand*nparams) ;      
        if combos(rq) < npresents/nparams
            combos(rq) = combos(rq) + 1 ;
            presentorder(presentcount) = rq ; 
            presentcount = presentcount + 1 ; 
        end  
    end  

    orientations = 0:16:360 ; 
    
    [xg,yg] = meshgrid(-300:.5:300,-300:.5:300) ; 
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
   
    % quadrants 
    [x,y] = meshgrid(-900:900,-900:900) ; 
    botright(:,:,2) = (1-((x>0) .* (y>0)))*255 ; botright(:,:,1) = x.*0 + gray ; quadtex(1) = Screen('MakeTexture',w,botright) ; 
    topright(:,:,2) = (1-((x>0) .* (y<0)))*255 ; topright(:,:,1) = x.*0 + gray ; quadtex(2) = Screen('MakeTexture',w,topright) ; 
    botleft(:,:,2) = (1-((x<0) .* (y>0)))*255 ; botleft(:,:,1) = x.*0 + gray ; quadtex(3) = Screen('MakeTexture',w,botleft) ; 
    topleft(:,:,2) = (1-((x<0) .* (y<0)))*255 ; topleft(:,:,1) = x.*0 + gray ; quadtex(4) = Screen('MakeTexture',w,topleft) ; 
    top(:,:,2) = (1-((y<0)))*255 ; top(:,:,1) = x.*0 + gray ; quadtex(5) = Screen('MakeTexture',w,top) ; 
    bot(:,:,2) = (1-((y>0)))*255 ; bot(:,:,1) = x.*0 + gray ; quadtex(6) = Screen('MakeTexture',w,bot) ; 
    right(:,:,2) = (1-((x>0)))*255 ; right(:,:,1) = x.*0 + gray ; quadtex(7) = Screen('MakeTexture',w,right) ; 
    left(:,:,2) = (1-((x<0)))*255 ; left(:,:,1) = x.*0 + gray ; quadtex(8) = Screen('MakeTexture',w,left) ; 
    fovea(:,:,2) = (1-(sqrt(x.^2+y.^2)<300))*255 ; fovea(:,:,1) = x.*0 + gray ; quadtex(9) = Screen('MakeTexture',w,fovea) ; 
    periphery(:,:,2) = (1-(sqrt(x.^2+y.^2)>300))*255 ; periphery(:,:,1) = x.*0 + gray ; quadtex(10) = Screen('MakeTexture',w,periphery) ; 
     
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
    tic
    currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
    for presentation = 1:size(presentorder,2) 
        restcount = 1 ; 
        focuscount = 1 ;
        stimcount = 1 ;
        poststimcount = 1 ;
        jitter = 0 ; 
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
                Screen('DrawTexture',w,xhair2tex) ;               
                Screen('DrawText',w,stimtypes{presentorder(presentation)},w2-400,h2-300) ; 
                if focuscount > focusframes + jitter
                    cs = ceil(presentorder(presentation)/10) ; 
                    rs = mod(presentorder(presentation)-1,10)+1 ; 
                    if trigsON
                        lptwrite(portn,str2num([num2str(cs),num2str(ct)])) ;
                        WaitSecs(0.004) ;
                        lptwrite(portn,0) ; 
                    end
                    disp(['cs = ',num2str(cs),' rs = ',num2str(rs)]) ; 
                    currentstate = 2 ;
                end
            end            
            if currentstate == 2 % present the stimulus                           
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
                if poststimcount > 0 ;% .5 seconds post stimulus fixation
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
