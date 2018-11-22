clear all ; close all
stimtime = 300 ; % stimulus time in seconds
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

randparams = 1; 

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

sizearrs = zeros(400,400,100); 
[xg,yg] = meshgrid(-100:.5:99.5,-100:.5:99.5); sinx = sin(xg); icount =1 ; 
maxcirc = 80; 
for i=0.001:0.01:1
    sizearrs(:,:,icount) = (sqrt((xg.^2 + yg.^2)) < i*maxcirc+10).*sinx; 
    icount = icount + 1; 
end
doublesizearrs = zeros(size(sizearrs,1),size(sizearrs,2),size(sizearrs,3)*2); 
doublesizearrs(:,:,1:end/2) = sizearrs; 
for i=1:size(sizearrs,3) ; doublesizearrs(:,:,end/2+i) = sizearrs(:,:,end-i+1) ; end
sizemask = zeros(size(xg,1),size(xg,2),2) ;
for i=1:size(doublesizearrs,3); 
sizemask(:,:,1) = uint8(mat2gray(doublesizearrs(:,:,i))*255) ; 
sizemask(:,:,2) = white ; 
sizemasktex(i) = Screen('MakeTexture',w,sizemask) ;
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

ntexes = length(sizemasktex); 
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
            %srcRect=[xoffset 0 xoffset + 400 400];   
            hwdiff = screenRect(3)-screenRect(4); 
            dispRect = [hwdiff/2,0,screenRect(3)-hwdiff/2,screenRect(4)]; 
            Screen('DrawTexture',w,sizemasktex(mod(round(stimcount*frameMultFactor),199)+1),[],dispRect); 
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;
            %Screen('DrawTexture',w,circmasktex,[],[]); 
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


