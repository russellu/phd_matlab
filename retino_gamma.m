clear all ; close all
HideCursor

cd c:/mscripts2/finalstims ; 
rands = load('randarrs') ; randarrs = rands.randarrs ; 

stimtime = 5 ; % stimulus time in seconds 
resttime = 3 ; % post stimulus time 
focustime = 1 ; % minimum focus time (+ 1-2s jitter)
nparams = 3 ; % 6 orientations for the starting angle, and the 2 annulus types
ntrials = 4 ;

HideCursor
Screen('Preference', 'SkipSyncTests', 1) ;
AssertOpenGL ;
screens=Screen('Screens') ;
screenNumber=max(screens) ;

white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);
gray=(white+black)/2;
    if round(gray)==white
        gray=black;
    end
inc=white-gray ;	

stimtrigs = repmat(1:nparams,[1,nparams*ntrials]) ; 
randinds = randperm(nparams*ntrials) ; 
stimtrigs = stimtrigs(randinds) ; 


% Open a double buffered fullscreen window and draw a gray background 
% to front and back buffers:
[w ,screenRect]=Screen('OpenWindow',screenNumber, gray) ;
w2 = screenRect(3)/2 ; h2 = screenRect(4)/2 ; 
screenWidth = screenRect(3) ; screenHeight = screenRect(4) ; 
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
visiblesize=screenWidth ; %2*texsize+1 ;

% for the different orientations:
[xg,yg] = meshgrid(-800:1:800,-800:1:800) ; 
sg = sin(xg) ; 
sgtex = Screen('MakeTexture',w,uint8(mat2gray(sg)*255)) ; 

% red crosshair (radius 8 pixels, dark red)
[x,y] = meshgrid(-50:50,-50:50) ; 
xhair = zeros(size(x,1),size(x,2),4) ; 
xhair(:,:,1) = gray/2 ; 
xhair(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 12) ; 
xhairtex = Screen('MakeTexture',w,xhair) ; 

 % smaller brighter red crosshair
xhairsmall = zeros(size(x,1),size(x,2),4) ;
xhairsmall(:,:,1) = gray ; 
xhairsmall(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 8) ; 
xhairsmalltex = Screen('MakeTexture',w,xhairsmall) ;

xhairsmallest = zeros(size(x,1),size(x,2),4) ;
xhairsmallest(:,:,1) = white ; 
xhairsmallest(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 4) ; 
xhairsmallesttex = Screen('MakeTexture',w,xhairsmallest) ;

% black crosshair (radius 8 pixels, dark red)
xhair2 = zeros(size(x,1),size(x,2),4) ;
xhair2(:,:,1) = black ; xhair2(:,:,2) = black ; xhair2(:,:,3) = black ; 
xhair2(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 12) ; 
xhair2tex = Screen('MakeTexture',w,xhair2) ;


randinds = [1,3,13] ; 
for i=1:length(randinds)
    randvals(:,:,1) = randarrs(:,:,randinds(i)) ; 
    randvals(:,:,2) = white ; 
    randtex(i) = Screen('MakeTexture',w,randvals) ;  
end



% Query duration of monitor refresh interval :
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
for presentation = 1:length(stimtrigs) ; 
    restcount = 1 ; 
    focuscount = 1 ;
    stimcount = 1 ;
    poststimcount = 1 ;
    jitter = 0 ; %rand*frameRate*4 ; %frameRate/2 + rand*frameRate ; % jitter .5-1.5s long
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
            Screen('DrawTexture',w,xhairsmallesttex) ;  
            if focuscount > focusframes + jitter
                currentstate = 2 ;
            end
        end            
        if currentstate == 2 % present the stimulus  
            xoffset = mod(stimcount*shiftperframe/2,900);  % Shift the grating by "shiftperframe" pixels per frame:
            srcRect=[xoffset 0 xoffset + 500 500];          
            Screen('DrawTexture',w,randtex(stimtrigs(presentation)),srcRect,[0,0,screenWidth,screenHeight]) ; 
            Screen('DrawTexture',w,xhairtex) ; 
            Screen('DrawTexture',w,xhairsmalltex) ;             
            Screen('DrawTexture',w,xhairsmallesttex) ;  

            stimcount = stimcount + 1 ; 
            if stimcount > stimframes 
                currentstate = 3 ;                
            end
        end              
        if currentstate == 3 % present the post stimulus fixation
            poststimcount = poststimcount + 1 ;
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;  
            Screen('DrawTexture',w,xhairsmallesttex) ;  
            if poststimcount > 0% .5 seconds post stimulus fixation
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

% don't forget the pause at the END!   
Screen('CloseAll') ;
toc