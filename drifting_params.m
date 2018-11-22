clear all ; close all
stimtime = 2 ; % stimulus time in seconds
resttime = 1 ; % post stimulus time
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

HideCursor
Screen('Preference', 'SkipSyncTests', 0);
AssertOpenGL;
screens=Screen('Screens');
screenNumber=max(screens);

white=WhiteIndex(screenNumber) ;
black=BlackIndex(screenNumber) ;
gray=(white+black)/2 ;
    if round(gray)==white
        gray=black;
    end
inc=white-gray ;	

% Open a double buffered fullscreen window and draw a gray background 
% to front and back buffers:
[w ,screenRect]=Screen('OpenWindow',screenNumber, gray) ;
w2 = screenRect(3)/2 ; h2 = screenRect(4)/2 ; 
screenWidth = screenRect(3) ; screenHeight = screenRect(4) ; 
squareRect = [0,0,screenRect(4),screenRect(4)] ; 
% Enable alpha-blending:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
% Show initial gray screen:
Screen('Flip', w) ;
f=0.05 ;
cyclespersecond=1.5 ;
angle=0 ;
texsize=250 ; % Half-Size of the grating image.
% Calculate parameters of the grating:
p=ceil(1/f) ;  % pixels/cycle    
fr=f*2*pi ;
visiblesize=200 ; %2*texsize+1 ;

randparams =[1,2,3,4,5]; 

% red crosshair (radius 8 pixels, dark red)
[x,y] = meshgrid(-600:600,-600:600) ; 
xhair = zeros(size(x,1),size(x,2),4) ;
xhair(:,:,1) = gray ; 
xhair(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 10) ; 
xhairtex = Screen('MakeTexture',w,xhair) ;

% smaller brighter red crosshair
xhairsmall = zeros(size(x,1),size(x,2),4) ;
xhairsmall(:,:,1) = white ; 
xhairsmall(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 5) ; 
xhairsmalltex = Screen('MakeTexture',w,xhairsmall) ;

% black crosshair (radius 8 pixels, dark red)
xhair2 = zeros(size(x,1),size(x,2),4) ;
xhair2(:,:,1) = black ; xhair2(:,:,2) = black ; xhair2(:,:,3) = black ; 
xhair2(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 10) ; 
xhair2tex = Screen('MakeTexture',w,xhair2) ;

% circle mask 
circmask = zeros(size(x,1),size(x,2),2) ;
circmask(:,:,1) = gray ; 
circmask(:,:,2) = white*double(sqrt(x.^2 + y.^2) > 350) ; 
circmasktex = Screen('MakeTexture',w,circmask) ;

contarrs = zeros(400,400,5); 
[xg,yg] = meshgrid(-200:1:199,-200:1:199); 
sinx = sin(xg); 
contlvls = [.05,.15,.25,.5,1]; 
for i=1:length(contlvls)
    contarrs(:,:,i) = 127.5+sinx*127.5*contlvls(i); 
end
for i=1:size(contarrs,3); 
contmask(:,:,1) = uint8((contarrs(:,:,i))) ; 
contmask(:,:,2) = white ; 
contmasktex(i) = Screen('MakeTexture',w,contmask) ;
end

inds = randperm(numel(sinx));
prevnpix = 1; prevsinx = sinx; 
rndlvls = [0.05,0.1,0.25,0.5,1]; 
for i=1:length(rndlvls)
   npix = rndlvls(i)*numel(sinx); 
   newsinx = sinx; 
   newsinx(inds(1:npix/2)) = sinx(inds(npix/2:npix-1)) ; 
   newsinx(inds(npix/2:npix-1)) = sinx(inds(1:npix/2)); 
   randarrs(:,:,i) = newsinx;
   prevsinx = newsinx; 
   prevnpix = npix; 
end

randmask = zeros(size(xg,1),size(xg,2),2) ;
for i=1:size(randarrs,3); 
randmask(:,:,1) = uint8(mat2gray(randarrs(:,:,i))*255) ; 
randmask(:,:,2) = white ; 
randmasktex(i) = Screen('MakeTexture',w,randmask) ;
end

circrads = [25,50,100,150,200]; 
circmask = zeros(size(xg)); 
for i=1:size(randarrs,3); 
circmask = zeros(size(xg,1),size(xg,2),2) ;
circmask(:,:,1) = gray ; 
circmask(:,:,2) = white*double(sqrt(xg.^2 + yg.^2) > circrads(i)) ; 
circmasktexes(i) = Screen('MakeTexture',w,circmask) ;
end

% Query duration of monitor refresh interval:
ifi=Screen('GetFlipInterval', w) ;

waitframes = 1;
waitduration = waitframes * ifi ;

% Translate requested speed of the grating (in cycles per second)
% into a shift value in "pixels per frame", assuming given
% waitduration: This is the amount of pixels to shift our "aperture" at
% each redraw:

% Perform initial Flip to sync us to the VBL and for getting an initial
% VBL-Timestamp for our "WaitBlanking" emulation:
vbl=Screen('Flip', w) ;
frameRate=Screen('FrameRate',screenNumber) ;
restframes = resttime*frameRate ; 
% Convert movieDuration in seconds to duration in frames to draw:
stimframes = stimtime * frameRate ; 
focusframes = focustime * frameRate ;

ntexes = length(contmasktex); 
ncycles = 10; % 10 complete cycles of the randomization mask. 
frameMultFactor = (ntexes*ncycles)/stimframes;

% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w) ;
Priority(priorityLevel) ;
startT = GetSecs ; 
currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
for presentation = 1:length(randparams)
    restcount = 1 ; 
    focuscount = 1 ;
    stimcount = 1 ;
    poststimcount = 1 ;
    jitter = 0 ; 
    DONE = false ; 
    trigsent = false ; 
    currentround = 0 ; 
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
            srcRect=[stimcount/10 0 stimcount/10 + 200 200];   
            hwdiff = screenRect(3)-screenRect(4); 
            dispRect = [hwdiff/2,0,screenRect(3)-hwdiff/2,screenRect(4)]; 
            Screen('DrawTexture',w,contmasktex(5),srcRect,dispRect); 
            Screen('DrawTexture',w,circmasktexes(presentation),[],dispRect); 
            Screen('DrawTexture',w,xhairtex);  
            Screen('DrawTexture',w,xhairsmalltex);
            Screen('DrawTexture',w,circmasktexes(3),[],dispRect); 
            stimcount = stimcount + 1 ; 
            if stimcount > stimtime*frameRate  
                currentstate = 3 ;         
            end
        end              
        if currentstate == 3 % present the post stimulus fixation
            poststimcount = poststimcount + 1 ;
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;  
            if poststimcount >= 0 
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
vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.      
%spause(10) ; 
Screen('CloseAll');


