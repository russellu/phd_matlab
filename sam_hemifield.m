clear all ; close all

cd c:/mscripts/finalstims/

stimtime = 25 ; % stimulus time in seconds
resttime = .75 ; % post stimulus time
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

HideCursor
Screen('Preference', 'SkipSyncTests', 1);
AssertOpenGL;
screens=Screen('Screens') ;
screenNumber=max(screens) ;

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

% generate the sequence in which the stimulus will be presented
% (pseudorandom)
nparams = 4 ; 
npresents = nparams*15 ;
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

% Create one single static grating image:
grating = sin(meshgrid(-200:.25:200,-50:.25:50)) ; 
grating = uint8(mat2gray(grating)*255) ; 
gratingtex = Screen('MakeTexture', w, grating) ;

% put the stimuli texture IDs in an array for easy indexing
stimarr = [gratingtex,gratingtex,gratingtex,gratingtex] ; 
trigs = [10,11,12,13] ;

% red crosshair (radius 8 pixels, dark red)
[x,y] = meshgrid(-600:600,-600:600) ; 
xhair = zeros(size(x,1),size(x,2),4) ;
xhair(:,:,1) = gray ; 
xhair(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 10) ; 
xhairtex = Screen('MakeTexture',w,xhair) ;

 % smaller brighter red crosshair
xhairsmall = zeros(size(x,1),size(x,2),4) ;
xhairsmall(:,:,1) = white ; 
xhairsmall(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 6) ; 
xhairsmalltex = Screen('MakeTexture',w,xhairsmall) ;

% black crosshair (radius 8 pixels, dark red)
xhair2 = zeros(size(x,1),size(x,2),4) ;
xhair2(:,:,1) = black ; xhair2(:,:,2) = black ; xhair2(:,:,3) = black ; 
xhair2(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 10) ; 
xhair2tex = Screen('MakeTexture',w,xhair2) ;

% circular aperture mask (left)
circmask(:,:,1) = x.*0 + gray ; 
circmask(:,:,2) = double(sqrt((x+300).^2 + (y).^2) > 240)*white ;
circmasktexright = Screen('MakeTexture',w,circmask) ; 
% circular aperture mask (right)
circmask(:,:,1) = x.*0 + gray ; 
circmask(:,:,2) = double(sqrt((x-300).^2 + (y).^2) > 240)*white ;
circmasktexleft = Screen('MakeTexture',w,circmask) ; 

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
shiftperframe= cyclespersecond * p * waitduration ;
shiftperframe = shiftperframe/2 ; 
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

currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on
smallblockcounter = 1 ;
juststarting = true ; 

jits = zeros(1,size(presentorder,2)) ; 
tic
for presentation = 1:size(presentorder,2) 
    restcount = 1 ; 
    focuscount = 1 ;
    stimcount = 1 ;
    poststimcount = 1 ;
    jitter = frameRate/3 + rand*frameRate ; %frameRate/2 + rand*frameRate ; % jitter .5-1.5s long
    jits(presentation) = jitter ; 
    DONE = false ; 
    trigsent = false ; 
    rotangle = 90 ; 
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
            xoffset = mod(stimcount*shiftperframe,360);  % Shift the grating by "shiftperframe" pixels per frame:
            srcRect=[xoffset 0 xoffset + visiblesize visiblesize] ;
            Screen('DrawTexture',w,stimarr(presentorder(presentation)),srcRect,[0,0,screenRect(3)*.5,screenRect(4)],rotangle) ;
            Screen('DrawTexture',w,circmasktexright,[],[0,0,screenRect(3),screenRect(4)]) ;
            Screen('DrawTexture',w,xhairtex) ;
            Screen('DrawTexture',w,xhairsmalltex) ;
            stimcount = stimcount + 1 ;
            
            if mod(stimcount,200)==0
                rotangle = rand*360 ; 
            end
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
save('jits','jits') ; 
Screen('CloseAll');
toc

