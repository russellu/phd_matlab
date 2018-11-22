clear all ; close all

%%%%% ALWAYS CHANGE TRIAL COUNT
trialCount = 1 ; 

number_of_presentations_per_condition = 7 ; % how many times to present each stimulus per run

stimTimes = {} ;
rTimes = {} ; 

cd c:/mscripts/finalfinalstims/

stimtime = 2 ; % stimulus time in seconds (time that stimulus is on)
resttime = .5 ; % rest time in seconds (time for black crosshair)
focustime = .5 ; % focus time in seconds (pre stimulus red crosshair)
posttime = .0 ; % post stimulus red crosshair time in seconds

totaltime =  stimtime + resttime + focustime + posttime ; 

TR = 2 ; % TR in seconds (for 3 TR pause at beginning)

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
    visiblesize=2*texsize+1 ;
    
    % create the masks
    smallsz = 75; largesz = 200 ; 
    [x,y] = meshgrid(-300:300,-300:300) ; 
    smask = sqrt(x.^2 + y.^2) < smallsz ; 
    smask2 = (sqrt(x.^2 + y.^2) < largesz & sqrt(x.^2 + y.^2) >smallsz) ;
    q1 = y>0&x>0 ; q2 = y>0&x<0 ; q3 = y<0&x>0 ; q4 = y<0&x<0 ; 
    s_br = smask.*q1 ; s_bl = smask.*q2 ; s_tr = smask.*q3 ; s_tl = smask.*q4 ; 
    l_br = smask2.*q1 ; l_bl = smask2.*q2 ; l_tr = smask2.*q3 ; l_tl = smask2.*q4 ; 
    alls = {s_tr,s_tl,l_tr,l_tl,s_br,s_bl,l_br,l_bl} ; 
    stims{1} = s_bl ; 
    stims{2} = s_tl+s_br ; 
    stims{3} = s_tr ; 
    stims{4} = l_br ;
    stims{5} = l_tr+l_bl ; 
    stims{6} = l_tl ;
    stims{7} = s_br ; 
    stims{8} = s_tr+s_bl ; 
    stims{9} = s_tl ; 
    stims{10} = l_bl ;
    stims{11} = l_tl+l_br ; 
    stims{12} = l_tr ;
    
    %for i=1:max(size(stims)) ; subplot(3,4,i) ; imshow(stims{i}) ; end
    % create the compatibility matrix 
    for i=1:size(stims,2) ; for j=1:size(stims,2) ; cmat(i,j) = max(max(stims{i}+stims{j})) ; end ; end 
    bcount = 1;  for i=1:size(cmat,1) ; for j=1:size(cmat,2) ; if j<i ; if cmat(i,j) ==2 ; badinds(bcount,:) = [i,j] ; bcount = bcount + 1 ; end ; end ; end ; end
    nparams = 12 ;  
    % generate the sequence in which the stimulus will be presented
    % (pseudorandom)
    % randperm approach
    clear randps ; for i=1:1000 ; randps(i,:) = randperm(nparams) ; end
    first = randps(1,:) ; presentorder = zeros(1,120) ; presentorder(1:nparams) = first ; 
    icount = 1 ; 
    while icount < max(size(presentorder))/nparams 
        for i=1:max(size(randps))
            index = ceil(rand*max(size(randps))) ; 
            inds1 = randps(index,1:2) ; inds2 = first(11:12) ; s = [inds1',inds2'] ; 
            if ~isempty(intersect(inds1,inds2)) && isempty(intersect(badinds,s,'rows')) 
                first = randps(index,:) ; presentorder(icount*nparams+1:icount*nparams+nparams) = first ; icount = icount + 1 ;                 
                break ; 
            end
        end
    end
    count = 1 ; for i=2:max(size(presentorder)) ; if cmat(presentorder(i-1),presentorder(i))==2 ; count = count + 1 ; disp('error') ; end ; end ; count
    % make mask textures
    for i=1:max(size(stims)) ; 
        circmask(:,:,1) = stims{i}.*0 + gray ; 
        circmask(:,:,2) = (1-stims{i}).*white ;
        circmasktex = Screen('MakeTexture',w,circmask) ; 
        allmasks(i) = circmasktex ; 
    end 
    % Create one single static grating image:
    grating = load('grating.mat') ; grating = grating.grating ; 
    gratingtex = Screen('MakeTexture', w, grating) ; 
    
    % put the stimuli texture IDs in an array for easy indexing
    stimarr = [gratingtex] ; 
    trigs = [10] ;
        
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
    
     %%% blue crosshair, for startup screen:
    xhairb = zeros(size(x,1),size(x,2),4) ;
    xhairb(:,:,3) = gray ; 
    xhairb(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 5) ; 
    xhairbtex = Screen('MakeTexture',w,xhairb) ;
      
    % smaller circular aperture mask
    smcircmask(:,:,1) = x.*0 + gray ; 
    smcircmask(:,:,2) = double(sqrt(x.^2 + y.^2) > 25)*white ;
    smcircmasktex = Screen('MakeTexture',w,smcircmask) ; 
    
    % make the multiple crosshair textures for the countdown
    mx1 = zeros(size(x,1),size(x,2),4) ;
    mx1(:,:,1) = black ; mx1(:,:,2) = black ; mx1(:,:,3) = black ; 
    mx1(:,:,4) = white*double(sqrt((x+20).^2 + (y).^2) < 5) ; 
    mx1tex = Screen('MakeTexture',w,mx1) ;
    
    mx2 = zeros(size(x,1),size(x,2),4) ;
    mx2(:,:,1) = black ; mx2(:,:,2) = black ; mx2(:,:,3) = black ; 
    mx2(:,:,4) = white*double(sqrt((x+40).^2 + (y).^2) < 5) ; 
    mx2tex = Screen('MakeTexture',w,mx2) ;
    
    mx3 = zeros(size(x,1),size(x,2),4) ;
    mx3(:,:,1) = black ; mx3(:,:,2) = black ; mx3(:,:,3) = black ; 
    mx3(:,:,4) = white*double(sqrt((x+60).^2 + (y).^2) < 5) ; 
    mx3tex = Screen('MakeTexture',w,mx3) ;
    
    mx4 = zeros(size(x,1),size(x,2),4) ;
    mx4(:,:,1) = black ; mx4(:,:,2) = black ; mx4(:,:,3) = black ; 
    mx4(:,:,4) = white*double(sqrt((x+80).^2 + (y).^2) < 5) ; 
    mx4tex = Screen('MakeTexture',w,mx4) ;
    
    mx5 = zeros(size(x,1),size(x,2),4) ;
    mx5(:,:,1) = black ; mx5(:,:,2) = black ; mx5(:,:,3) = black ; 
    mx5(:,:,4) = white*double(sqrt((x+100).^2 + (y).^2) < 5) ; 
    mx5tex = Screen('MakeTexture',w,mx5) ;

    % Query duration of monitor refresh interval:
    ifi=Screen('GetFlipInterval', w) ;
    
    waitframes = 1;
    waitduration = waitframes * ifi ;
    
    % Translate requested speed of the grating (in cycles per second)
    % into a shift value in "pixels per frame", assuming given
    % waitduration: This is the amount of pixels to shift our "aperture" at
    % each redraw:
    shiftperframe= cyclespersecond * p * waitduration ;
    
    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp for our "WaitBlanking" emulation:
    vbl=Screen('Flip', w) ;
    
    frameRate=Screen('FrameRate',screenNumber) ;
    % Convert movieDuration in seconds to duration in frames to draw:
    stimframes = stimtime * frameRate ; 
    restframes = resttime * frameRate ;
    focusframes = focustime * frameRate ;
    postframes = posttime * frameRate ;
    % Use realtime priority for better timing precision:
    priorityLevel=MaxPriority(w) ;
    Priority(priorityLevel) ;
    
    currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on
    smallblockcounter = 1 ;
    juststarting = true ; 
        
    
       %%% INITIATE COUNTDOWN
        
                Screen('DrawTexture',w,xhair2tex) ;Screen('DrawTexture',w,mx5tex) ;Screen('DrawTexture',w,mx4tex) ;Screen('DrawTexture',w,mx3tex) ;Screen('DrawTexture',w,mx2tex) ;Screen('DrawTexture',w,mx1tex) ; 
                Screen('Flip', w, vbl + (waitframes - 0.5) * ifi); pause(1) ; 
                Screen('DrawTexture',w,xhair2tex) ;Screen('DrawTexture',w,mx4tex) ;Screen('DrawTexture',w,mx3tex) ;Screen('DrawTexture',w,mx2tex) ;Screen('DrawTexture',w,mx1tex) ; 
                Screen('Flip', w, vbl + (waitframes - 0.5) * ifi); pause(1) ; 
                Screen('DrawTexture',w,xhair2tex) ;Screen('DrawTexture',w,mx3tex) ;Screen('DrawTexture',w,mx2tex) ;Screen('DrawTexture',w,mx1tex) ; 
                Screen('Flip', w, vbl + (waitframes - 0.5) * ifi); pause(1) ; 
                Screen('DrawTexture',w,xhair2tex) ;Screen('DrawTexture',w,mx2tex) ;Screen('DrawTexture',w,mx1tex) ; 
                Screen('Flip', w, vbl + (waitframes - 0.5) * ifi); pause(1) ; 
                Screen('DrawTexture',w,xhair2tex) ;Screen('DrawTexture',w,mx1tex) ; 
                Screen('Flip', w, vbl + (waitframes - 0.5) * ifi); pause(1) ;      
                Screen('DrawTexture',w,xhair2tex) ; 
                Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);     
                Screen('DrawTexture',w,xhairbtex) ;  
                Screen('DrawTexture',w,xhairbtex) ;  
                Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);       
                                               
                pause(3*TR)
            
                startTime = GetSecs ;
    
    for presentation = 1:size(presentorder,2) 
        save(['stimTimes',num2str(trialCount)],'stimTimes') ; 
        restcount = 1 ; 
        focuscount = 1 ;
        stimcount = 1 ;
        poststimcount = 1 ;
        jitter = 0 ;
        DONE = false ; 
        trigsent = false ; 
        rTimes{size(rTimes,2)+1} = GetSecs ; 
        
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
                if focuscount > focusframes 
                    currentstate = 2 ;                    
                    stimTimes{size(stimTimes,2)+1} = [GetSecs-startTime,presentorder(presentation)] ;                     
                end
            end            
            if currentstate == 2 % present the stimulus               
                xoffset = mod(stimcount*shiftperframe,360*8);  % Shift the grating by "shiftperframe" pixels per frame:
                srcRect=[xoffset 0 xoffset + visiblesize visiblesize];                     
                Screen('DrawTexture',w,stimarr(1),srcRect,[]) ;                  
                Screen('DrawTexture',w,allmasks(presentorder(presentation))) ;  
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
                stimcount = stimcount +1 ; 
                if stimcount > stimframes 
                    currentstate = 3 ;                
                end
            end              
            if currentstate == 3 % present the post stimulus fixation
                poststimcount = poststimcount + 1 ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
                if poststimcount > postframes 
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
   vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.           

    pause(5)
    FEND = GetSecs ;   
    stend = [startTime,FEND] ;     
    save(['stimTimes',num2str(trialCount)],'stimTimes') ; 
	Screen('CloseAll');
    diff(stend)   
    for i=1:size(stimTimes,2)
       timess(i) = stimTimes{i}(1) ;      
    end
    