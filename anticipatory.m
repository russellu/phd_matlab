
clear all ; close all

cd c:/mscripts/finalstims/

stimtime = 2 ; % stimulus time in seconds
resttime = 2 ; % post stimulus time
focustime = 2 ; % minimum focus time (+ 1-2s jitter)

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
grating = load('grating.mat') ; grating = grating.grating ; 
gratingtex = Screen('MakeTexture', w, grating) ;

% create the low freq grating:
lowgrating = load('lowgrating.mat') ; lowgrating = lowgrating.lowgrating ;
lowgratingtex = Screen('MakeTexture',w,lowgrating) ; 

% create the three contrast levels
c1 = load('c1.mat') ; c1 = c1.c1 ; 
contrastlow = Screen('MakeTexture',w,c1) ;
c2 = load('c2.mat') ; c2 = c2.c2 ;  
contrastmed = Screen('MakeTexture',w,c2) ;

% create the plaid to overlay on the grating
plaid = load('plaid.mat') ; plaid = plaid.plaid ; 
plaidtex = Screen('MakeTexture', w, plaid) ;

% create the two randomization levels
randarrs = load('randarrs.mat') ; randarrs = randarrs.randarrs ; 
rnd10 = squeeze(randarrs(:,:,2)) ; 
rnd10tex = Screen('MakeTexture', w, rnd10) ;

rnd50 = uint8(mat2gray(squeeze(randarrs(:,:,3)))*255) ;
rnd50tex = Screen('MakeTexture', w, rnd50) ;

% put the stimuli texture IDs in an array for easy indexing
stimarr = [gratingtex,gratingtex,gratingtex,gratingtex] ; 
trigs = [10,11,12,13] ;

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

% circular aperture mask
circmask(:,:,1) = x.*0 + gray ; 
circmask(:,:,2) = double(sqrt(x.^2 + y.^2) > 150)*white ;
circmasktex = Screen('MakeTexture',w,circmask) ; 

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
shiftperframe = shiftperframe/10 ; 
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
            srcRect=[xoffset 0 xoffset + visiblesize visiblesize];   
            if presentorder(presentation) ~= 4
                Screen('DrawTexture',w,stimarr(presentorder(presentation)),srcRect,[]) ;            
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
                Screen('DrawTexture',w,circmasktex) ;            
            else
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
                Screen('DrawTexture',w,circmasktex) ;  
            end
            
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
save('jits','jits') ; 
Screen('CloseAll');
toc

