clear all ; close all
stimtime = 120 ; % stimulus time in seconds
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
frameRate=Screen('FrameRate',screenNumber) ;

randparams = 1:10; 

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

[xg,yg] = meshgrid(-800:.5:800,-100:.5:100) ; xg = xg(:,1:end-1); yg = yg(:,1:end-1); 
sinx = sin(xg)<-.35; newsin = zeros(size(sinx)); 
sizesin = zeros(size(sinx,1),size(sinx,2),2) ;
sizesin(:,:,2) = white ; 
sizesin(:,:,1) = uint8(mat2gray(sinx)*255) ; 
sizesintex = Screen('MakeTexture',w,sizesin) ;

ncycles = 4; 
radcount = 400; 
radmin = 20; radmax=450; 
radstep = (radmax-radmin)/radcount; icount=0;
for i=radmin:radstep:radmax-radstep
    icount = icount +1; 
    circmask = zeros(size(x,1),size(x,2),2) ;
    circmask(:,:,1) = gray ; 
    circmask(:,:,2) = white*double(sqrt(x.^2 + y.^2) > radmin+radstep*icount) ; 
    circularmasks(icount) = Screen('MakeTexture',w,circmask) ;
end
circincr = (icount*ncycles*2)/(frameRate*stimtime);
circinds = ceil(0:circincr:icount); circinds(circinds==0) = 1 ; circinds(circinds>icount) = icount; 
newcircinds = zeros(1,length(circinds)*2); 
newcircinds(1:end/2) = circinds; newcircinds(end/2+1:end) = fliplr(circinds); 
repcircinds = repmat(newcircinds,[1,4]); 


% Query duration of monitor refresh interval:
ifi=Screen('GetFlipInterval', w) ;

waitframes = 1;
waitduration = waitframes * ifi ;

% Translate requested speed of the grating (in cycles per second)
% into a shift value in "pixels per frame", assuming given
% waitduration: This is the amount of pixels to shift our "aperture" at
% each redraw:
%shiftperframe = cyclespersecond * p * waitduration ;

% Perform initial Flip to sync us to the VBL and for getting an initial
% VBL-Timestamp for our "WaitBlanking" emulation:
vbl=Screen('Flip', w) ;
restframes = resttime*frameRate ; 
shiftperframe = ((frameRate/stimtime)*(size(sinx,2)/screenRect(3)))/(2.1); 
% Convert movieDuration in seconds to duration in frames to draw:
stimframes = stimtime * frameRate ; 
focusframes = focustime * frameRate ;
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
            xoffset = stimcount*shiftperframe;  % Shift the grating by "shiftperframe" pixels per frame:
            srcRect=[xoffset 0 xoffset + 400 400];   
            hwdiff = screenRect(3)-screenRect(4); 
            dispRect = [hwdiff/2,0,screenRect(3)-hwdiff/2,screenRect(4)]; 
            Screen('DrawTexture',w,sizesintex,srcRect,dispRect); 
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;
            Screen('DrawTexture',w,circularmasks(repcircinds(stimcount)),[],[]); 
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


 