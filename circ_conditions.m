clear all ; close all

stimtime = 5 ; % stimulus time in seconds
resttime = 2.5 ; % post stimulus time
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
    visiblesize=3000 ; %2*texsize+1 ;
    
    % generate the sequence in which the stimulus will be presented
    % (pseudorandom)
    nparams = 15 ; 
    npresents = nparams*1 ;
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

    clear allstims ; 
    sfreqcount = 1 ;
    for sfreq=0.25:0.25:1.25
        tfreqcount = 1 ; 
        for tfreq=0.25:.25:.75
            lim = sfreq*200 ; 
            [xg,yg] = meshgrid(-lim:sfreq:lim,-lim:sfreq:lim) ; 
            [th,rho] = cart2pol(xg,yg) ;
            for i=1:50
                allstims{sfreqcount}(tfreqcount,i,:,:) = sin(rho+mod(i,50)*tfreq) ; 
            end   
            tfreqcount = tfreqcount + 1 ; 
        end
        sfreqcount = sfreqcount + 1 ; 
    end
    
    for i=1:length(allstims)
        for j=1:size(allstims{i},1)
            for k=1:size(allstims{i},2)
                sintex(i,j,k) = Screen('MakeTexture',w,uint8(mat2gray(squeeze(allstims{i}(j,k,:,:)))*255)) ;  
            end
        end
    end
    
    % red crosshair (radius 8 pixels, dark red)
    [x,y] = meshgrid(-200:200,-200:200) ; 
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
    circmask(:,:,2) = double(sqrt(x.^2 + y.^2) > 500)*white ;
    circmasktex = Screen('MakeTexture',w,circmask) ; 
   
    %{
    [x,y] = meshgrid(-900:900,-900:900) ; 
    botright(:,:,2) = (1-((x>0) .* (y>0)))*255 ; botright(:,:,1) = x.*0 + gray ; quadtex(1) = Screen('MakeTexture',w,botright) ; 
    topright(:,:,2) = (1-((x>0) .* (y<0)))*255 ; topright(:,:,1) = x.*0 + gray ; quadtex(2) = Screen('MakeTexture',w,topright) ; 
    botleft(:,:,2) = (1-((x<0) .* (y>0)))*255 ; botleft(:,:,1) = x.*0 + gray ; quadtex(3) = Screen('MakeTexture',w,botleft) ; 
    topleft(:,:,2) = (1-((x<0) .* (y<0)))*255 ; topleft(:,:,1) = x.*0 + gray ; quadtex(4) = Screen('MakeTexture',w,topleft) ; 
    %}
      
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
        jitter = frameRate/3 + rand*frameRate ; %frameRate/2 + rand*frameRate ; % jitter .5-1.5s long
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
                    cs = floor((presentorder(presentation)-1)/size(sintex,1))+1 ; 
                    ct = mod(presentorder(presentation)-1,3)+1
                    disp(['cs = ',num2str(cs),' ct = ',num2str(ct)]) ; 
                    currentstate = 2 ;
                end
            end            
            if currentstate == 2 % present the stimulus               
                xoffset = 0 ;               
                wg = size(allstims{cs},3) ; hg = size(allstims{cs},4) ; 
                srcRect=[xoffset 0 xoffset + wg hg] ;    
                Screen('DrawTexture',w,sintex(cs,ct,mod(ceil(stimcount),50)+1),srcRect,[-wg*1.5+w2-1,-hg*1.5+h2-1,wg*1.5+w2-1,hg*1.5+h2-1]) ;   
                %Screen('DrawTexture',w,quadtex(ct)) ; 
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
