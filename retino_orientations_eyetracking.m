clear all ; close all
sub = 'greg' ; 
%name = [sub,'_stim_left'] ; 
name = [sub,'retinotopy'] ; 
% cryocooler => sequence dependent
%handle = CedrusResponseBox('Open','COM3') ; 

% Eye-tracker variables
hostname = '192.168.1.222';
order = 1;

% Connect to eye tracker recorder
%eye_connect = getEyeTrackerConnect(hostname);

%outlet = create_marker_outlet();

%{
lptwrite(57424,0) %Initiate triggers
lptwrite(57424,0)
lptwrite(57424,98) % start trials
%}
% Start triggers


%eyeTrigger('str', order, eye_connect);
%order = order + 1;

stimtime = 18 ; % stimulus time in seconds
resttime = 2 ; % post stimulus time
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

HideCursor
Screen('Preference', 'SkipSyncTests', 1);
AssertOpenGL;
screens=Screen('Screens');
screenNumber=max(screens);

trigsON = false ; 

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
cyclespersecond=6 ;
angle=0 ;
texsize=250 ; % Half-Size of the grating image.
% Calculate parameters of the grating:
p=ceil(1/f) ;  % pixels/cycle    
fr=f*2*pi ;
visiblesize=3000 ; 

randparams = ([1,2,3,4]) ; 

% create the checkerboard
[x,y] = meshgrid(-150:150,-150:150) ; 
[theta,rho] = cart2pol(x,y) ; 
radincrs = 10:10:200 ; radincrs = cumsum(radincrs) ; 
evods = mod(1:length(radincrs),2) ; clear rhos ; 
for i=1:length(radincrs)
    if i==1 ; 
        rhos(:,:,i) = (rho <= radincrs(i)).*evods(i) ; 
    else
        rhos(:,:,i) = (rho <= radincrs(i) & rho > radincrs(i-1)).*evods(i) ; 
    end
end

thetaincrs = -pi:.31:pi ; 
evods = mod(1:length(thetaincrs),2) ; 
clear thetas ; 
for i=1:length(thetaincrs)
    if i==1
        thetas(:,:,i) = (theta <= thetaincrs(i)).*evods(i) ;
    else
        thetas(:,:,i) = (theta <= thetaincrs(i) & theta > thetaincrs(i-1)).*evods(i) ; 
    end
end

srhos = sum(rhos,3) ; 
sthetas = sum(thetas,3) ; 
sthetas(srhos==1) = ~sthetas(srhos==1) ; 
sthetas = uint8(mat2gray(sthetas)*255) ; 
sthetas = uint8(ones(size(sthetas))*255) ; 

sthetas2 = sum(thetas,3) ; 
sthetas2(srhos==0) = ~sthetas(srhos==0) ; 
sthetas2 = uint8(mat2gray(sthetas2)*255) ; 
sthetas2 = uint8(zeros(size(sthetas2))*255) ; 

checktex1 = Screen('MakeTexture',w,sthetas) ; 
checktex2 = Screen('MakeTexture',w,sthetas2) ; 

thetaonly1 = sum(thetas,3) ; 
thetaonly1 = uint8(mat2gray(thetaonly1)*255) ; 
thetaonly2 = ~sum(thetas,3) ; 
thetaonly2 = uint8(mat2gray(thetaonly2)*255) ; 

thetatex1 = Screen('MakeTexture',w,thetaonly1) ; 
thetatex2 = Screen('MakeTexture',w,thetaonly2) ; 

% put the stimuli texture IDs in an array for easy indexing
trigs = [1] ;

blackb = zeros(150,150) ; 
whiteb = ones(150,150)*255 ; 
blacktex = Screen('MakeTexture',w,blackb) ;
whitetex = Screen('MakeTexture',w,whiteb) ;

% red crosshair (radius 8 pixels, dark red)
[x,y] = meshgrid(-150:150,-150:150) ; 
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
restframes = resttime*frameRate ; 
% Convert movieDuration in seconds to duration in frames to draw:
stimframes = stimtime * frameRate ; 

focusframes = focustime * frameRate ;
% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w) ;
Priority(priorityLevel) ;
hzs = 16 ; 
for h=1:length(hzs) ; 
    stimarr = zeros(1,stimframes) ; 
    hz = hzs(h) ; fs = frameRate/hz ; 
    icount = 0 ; 
    for i=1:fs:length(stimarr)
        if mod(icount,2) == 0 
            stimarr(i:i+fs) = 1 ; 
        end
        icount = icount + 1 ; 
    end
    stimarrs(h,:) = stimarr ; 
end    

% make the stimuli:
[xg,yg] = meshgrid(-150:150,-150:150) ; 
[th,rh] = cart2pol(xg,yg) ; 
quad = th>0 & th<pi/2 ; 
wedgemat(:,:,2) = uint8(mat2gray(1-quad)*255) ; 
wedgemat(:,:,1) = gray ; 
wedgetex = Screen('MakeTexture',w,wedgemat) ; 

[xg2,yg2] = meshgrid(-180:.5:180,-180:.5:180) ; 
sinu(:,:,2) = uint8(ones(size(xg2))*255) ; 
sinu(:,:,1) = uint8(sin(xg2)*255) ; 
sinutext = Screen('MakeTexture',w,sinu) ; 

circmask = double(sqrt(xg.^2 + yg.^2) < 90);
texmat(:,:,2) = uint8((mat2gray(1-circmask))*255) ; 
texmat(:,:,1) = gray ; 
circtex = Screen('MakeTexture',w,texmat) ; 
clear rhosin
for i=1:52 
    [xg,yg] = meshgrid(-250:250,-250:250) ; 
    [th,rh] = cart2pol(xg,yg) ;
    rhosin(:,:,i) = sin(rh+i/4) ; clear tmat ; 
    tmat(:,:,2) = uint8(zeros(size(xg))+255) ; 
    tmat(:,:,1) = uint8((rhosin(:,:,i))*255) ; 
    gratetex(i) = Screen('MakeTexture',w,tmat) ; 
end

%circinds = find(wedgemat(:,:,2)==255) ; circinds2 = find(circmask(:,:,1)==0) ; 
%a = (double(tmat(:,:,1)).*double(wedgemat(:,:,2)==0)) ; a=a.*circmask(:,:,1) ; a(circinds) = 122.5  ; a(circinds2) = 122.5 ; 

% Eye tracker connection start
%{
xhairwhite = zeros(size(x,1),size(x,2),4) ;
xhairwhite(:,:,2) = white; % Green
xhairwhite(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 10) ; % Alpha
xhairwhitetex = Screen('MakeTexture',w,xhairwhite) ;
Screen('DrawTexture',w,xhairtex) ;  
Screen('DrawTexture',w,xhairwhitetex) ; 
vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.           
        
while ~outlet.have_consumers()
    if KbCheck% Abort demo if any key is pressed:
       Screen('CloseAll');
       break
    end
    pause(5);
end
Screen('DrawTexture',w,xhair2tex);
vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
pause(5);
outlet.push_sample({'str'});
%}

rotstep = 400/(stimframes+10) ; 
rotangles = 1:rotstep:399.5 ; 
lrotangles = length(rotangles) ; 
startT = GetSecs ; 
tic
currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
for presentation = 1:length(randparams)
    restcount = 1 ; 
    focuscount = 1 ;
    stimcount = 1 ;
    poststimcount = 1 ;
    jitter = 0 ; %frameRate/3 + rand*frameRate ; %frameRate/2 + rand*frameRate ; % jitter .5-1.5s long
    DONE = false ; 
    trigsent = false ; 
    stimarr = stimarrs(1,:) ; 
    currentround = 0 ; 
    while ~DONE             
        if currentstate == 0 % present the black crosshair (rest period)
            restcount = restcount + 1 ;   
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;            
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
                stimtimes(1,presentation) = GetSecs - startT ;
                order = order + 1 ;
            end
        end            
        if currentstate == 2 % present the stimulus               
            xoffset = 0 ; % mod(stimcount*shiftperframe,360);  % Shift the grating by "shiftperframe" pixels per frame:
            srcRect=[0 0 screenRect(3) screenRect(4)];   
            sqrRect = [(screenRect(3)-screenRect(4))/2,0,screenRect(4)+(screenRect(3)-screenRect(4))/2,screenRect(4)] ; 

            if randparams(presentation)==1
                Screen('DrawTexture',w,gratetex(floor(mod(stimcount/2,51)+1)),[],sqrRect) ;
                Screen('DrawTexture',w,wedgetex,[],sqrRect,rotangles(stimcount)) ; 
                fwdrots(stimcount) = rotangles(stimcount) ; 
            elseif randparams(presentation)==2 
                Screen('DrawTexture',w,sinutext,[],[],rotangles(stimcount)) ; 
            elseif randparams(presentation)==3
                Screen('DrawTexture',w,gratetex(floor(mod(stimcount/2,51)+1)),[],sqrRect) ;
                Screen('DrawTexture',w,wedgetex,[],sqrRect,lrotangles - rotangles(stimcount)) ; 
                revrots(stimcount) = rotangles(lrotangles-stimcount) ; 
            elseif randparams(presentation)==4
                Screen('DrawTexture',w,sinutext,[],[], rotangles(lrotangles-stimcount)) ; 
            end
            Screen('DrawTexture',w,circtex,[],sqrRect) ; 
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;
            stimcount = stimcount + 1 ; 
            if stimcount > stimframes 
                currentstate = 3 ;         
                stimtimes(2,presentation) = GetSecs - startT ; 
            end
        end              
        if currentstate == 3 % present the post stimulus fixation
            poststimcount = poststimcount + 1 ;
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;  
            if poststimcount > 0 %frameRate/2 % .5 seconds post stimulus fixation
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
pause(10) ; 
Screen('CloseAll');
save('stimtimes','stimtimes') ; 
toc
