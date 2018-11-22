clear all ; close all
HideCursor

portn = 57424 ; 

trigsON = false ; 

if trigsON
    lptwrite(portn,0) %Initiate triggers
    lptwrite(portn,0)
    lptwrite(portn,1) % start trials
end
stimtime = 1 ; % stimulus time in seconds
resttime = 1 ; % post stimulus time
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
    nparams = 20 ; 
    npresents = nparams*2 ;
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

    sfreq = 1.35 ; tfreq = 0.2 ; 
    lim = sfreq*100 ;
    [xg,yg] = meshgrid(-lim:sfreq:lim,-lim:sfreq:lim) ; 
    [~,rho] = cart2pol(xg,yg) ;
    for i=1:250
        allstims{1}(i,:,:) = sin(rho+mod(i,250)*tfreq) ; 
    end   
    sfreq = 0.15 ; tfreq = 0.85 ; 
    lim = sfreq*100 ;
    [xg,yg] = meshgrid(-lim:sfreq:lim,-lim:sfreq:lim) ; 
    [~,rho] = cart2pol(xg,yg) ;
    for i=1:250
        allstims{2}(i,:,:) = sin(rho+mod(i,250)*tfreq) ; 
    end   
    for i=1:length(allstims)
        for k=1:size(allstims{i},1)
            sintex(i,k) = Screen('MakeTexture',w,uint8(mat2gray(squeeze(allstims{i}(k,:,:)))*255)) ;  
        end        
    end
    
    %{
    clear allstims ; 
    sfreqcount = 1 ;
    for sfreq=0.25:0.25:1.25
        tfreqcount = 1 ; 
        for tfreq=0.25:.25:.75
            lim = sfreq*150 ; 
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
    %}
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
    circmask(:,:,2) = double(sqrt(x.^2 + y.^2) > 500)*white ;
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
                    %cs = floor((presentorder(presentation)-1)/size(sintex,2))+1 ; 
                    %ct = mod(presentorder(presentation)-1,3)+1 ; 
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
                xoffset = 0 ;               
                wg = size(allstims{cs},2) ; hg = size(allstims{cs},3) ; 
                srcRect=[xoffset 0 xoffset + wg hg] ; scfact = 2.5 ; 
                Screen('DrawTexture',w,sintex(cs,mod(ceil(stimcount),250)+1),srcRect,[-wg*scfact+w2-1,-hg*scfact+h2-1,wg*scfact+w2-1,hg*scfact+h2-1]) ;   
                Screen('DrawTexture',w,quadtex(rs)) ; 
                Screen('DrawTexture',w,circmasktex) ;                 
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
