clear all ; close all

% cryocooler => sequence dependent

stimtime = 2 ; % stimulus time in seconds
resttime = 1.75 ; % post stimulus time
tetanictime = 20; 
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

nstims_pre = 20;
nstims_post = 20; 

HideCursor
Screen('Preference', 'SkipSyncTests', 1) ;
AssertOpenGL ;
screens=Screen('Screens') ;
screenNumber=max(screens) ;

trigsON = false ; 

%lptwrite(57424,0) %Initiate triggers
%lptwrite(57424,0)
%lptwrite(57424,98) %Start trials

white=WhiteIndex(screenNumber) ;
black=BlackIndex(screenNumber) ;
gray=(white+black)/2 ;
if round(gray)==white
    gray=black ;
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

% 1+2+2+4+8 = 17 retinotopic, 16 orientations and 16 plaids 

trigs_pre = repmat([1,2],[1,nstims_pre]); trigs_pre = trigs_pre(randperm(length(trigs_pre))); 
trigs_post = repmat([11,12],[1,nstims_post]); trigs_post = trigs_post(randperm(length(trigs_post))); 

trigs = [trigs_pre,3,trigs_post] ; 
randparams = trigs;

% create the checkerboard
[x,y] = meshgrid(-150:150,-150:150) ; 
[theta,rho] = cart2pol(x,y) ; 
radincrs = 10:10:200 ; radincrs = cumsum(radincrs) ; 
evods = mod(1:length(radincrs),2) ; clear rhos ; 
for i=1:length(radincrs)
    if i==1 
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
trigs = [1,2,3,4] ;

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
tetanicframes = tetanictime*frameRate; 

focusframes = focustime * frameRate ;
% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w) ;
Priority(priorityLevel) ;
hzs = 13 ; 
for h=1:length(hzs) 
    stimarr = zeros(1,tetanicframes) ; 
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

% make the retinotopic stimuli:
[xg,yg] = meshgrid(-150:150,-150:150) ; 
width = length(xg) ; height = length(yg) ; 
circmask = double(sqrt(xg.^2 + yg.^2) < 350);
fullfield = sqrt(xg.^2 + yg.^2) < 150 ; fullfield = fullfield.*circmask ;
left = (xg<0) .* circmask ; right = (xg>=0) .* circmask ;
allstims(:,:,1) = fullfield ; allstims(:,:,2) = left ; allstims(:,:,3) = right ; 

[xgg2,ygg2] = meshgrid(-100:1:100,-400:1:400) ; 
xgg2(:) = 0; 
for i=2:4:size(xgg2,1)-4; xgg2(:,i-1:i) = 1; end
sinu(:,:,2) = uint8(ones(size(xgg2))*255) ; 
sinu(:,:,1) = uint8((((xgg2))*255)) ; 
sinutext = Screen('MakeTexture',w,sinu) ; 

xgg3 = zeros(size(xgg2)); 
for i=4:8:size(xgg3,1)-4; xgg3(:,i-3:i) = 1; end
sinulow(:,:,2) = uint8(ones(size(xgg3))*255) ; 
sinulow(:,:,1) = uint8((((xgg3))*255)) ; 
sinutextlow = Screen('MakeTexture',w,sinulow) ; 

invsinu(:,:,2) = uint8(ones(size(xgg2))*255) ; 
invsinu(:,:,1) = uint8((((1-xgg2))*255)) ; 
invsinutext = Screen('MakeTexture',w,invsinu) ; 

for i=1:size(allstims,3)
    texmat(:,:,2) = uint8((1-allstims(:,:,i))*255) ; 
    texmat(:,:,1) = gray ; 
    alltex(i) = Screen('MakeTexture',w,texmat) ; 
end

% make the grating: 
circmask = double(sqrt(xg.^2 + yg.^2) < 60);
texmat(:,:,2) = uint8(mat2gray(1-circmask)*255) ; 
texmat(:,:,1) = gray ; 
circtex = Screen('MakeTexture',w,texmat) ; 

startT = GetSecs ; 
tic
currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
for presentation = 1:length(randparams)
    restcount = 1 ; 
    focuscount = 1 ;
    stimcount = 1 ;
    poststimcount = 1 ;
    jitter = (1.5*frameRate) + (rand-.5)*frameRate ; 
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
                %lptwrite(57424,randparams(presentation)) ;
                %WaitSecs(0.004) ;
                %lptwrite(57424,0) ;
            end
        end            
        if currentstate == 2 % present the stimulus       
            
            srcRect=[0 0 size(xgg2,1) size(xgg2,2)];                      
            sqrRect = [(screenRect(3)-screenRect(4))/2, 0, screenRect(4)+(screenRect(3)-screenRect(4))/2, screenRect(4)] ;
            sqrRect(1:2) = sqrRect(1:2)+100; sqrRect(3:4) = sqrRect(3:4)-100; 
            largeRect = [sqrRect(1)-500,sqrRect(2)-500,sqrRect(3)+500,sqrRect(4)+500] ; 

            %Screen('DrawTexture',w,sinutext,srcRect,sqrRect) ;      
            
            if randparams(presentation)==3
                if stimarrs(stimcount)==1 
                    Screen('DrawTexture',w,sinutext,[0 0 size(xg,1) size(xg,2)],sqrRect) ; 
                else
                    Screen('DrawTexture',w,invsinutext,[0 0 size(xg,1) size(xg,2)],sqrRect) ; 
                end                              
            elseif randparams(presentation)==1 || randparams(presentation)==11
                if stimcount > frameRate*.5
                    xoffset = ((stimcount - frameRate*.5)*shiftperframe)/2;  % Shift the grating by "shiftperframe" pixels per frame:
                    mvsrcRect=[0+xoffset/4 0 size(xg,1)+xoffset/4 size(xg,2)];   
                    Screen('DrawTexture',w,sinutext,mvsrcRect,sqrRect) ;    
                else
                    Screen('DrawTexture',w,sinutext,[0 0 size(xg,1) size(xg,2)],sqrRect) ;  
                end
            elseif randparams(presentation)==2 || randparams(presentation)==12
                if stimcount > frameRate*.5
                    xoffset = ((stimcount - frameRate*.5)*shiftperframe)/2;  % Shift the grating by "shiftperframe" pixels per frame:
                    mvsrcRect=[0+xoffset/4 0 size(xg,1)+xoffset/4 size(xg,2)];   
                    Screen('DrawTexture',w,sinutextlow,mvsrcRect,sqrRect) ;    
                else
                    Screen('DrawTexture',w,sinutextlow,[0 0 size(xg,1) size(xg,2)],sqrRect) ;  
                end
            end
             
            Screen('DrawTexture',w,circtex,[],largeRect) ; 
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;
            stimcount = stimcount + 1 ; 
            if randparams(presentation) == 3
                if stimcount > tetanicframes 
                    currentstate = 3 ;         
                    stimtimes(2,presentation) = GetSecs - startT ; 
                end
            elseif stimcount > stimframes 
                currentstate = 3 ;         
                stimtimes(2,presentation) = GetSecs - startT ; 
            end
        end              
        if currentstate == 3 % present the post stimulus fixation
            poststimcount = poststimcount + 1 ;
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;  
            if poststimcount > 0 %0 %frameRate/2 % .5 seconds post stimulus fixation
               currentstate = 0 ; % restart the loop  
               DONE = true ;
            end    
        end         
        vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.           
        if KbCheck % Abort demo if any key is pressed:
           Screen('CloseAll');
           break
        end
    end
end
%save(name,'stimtimes') ; 
vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw. 
%lptwrite(57424,0) ;%Initiate triggers
%lptwrite(57424,0) ;
%lptwrite(57424,99) ;% start trials
pause(10) ;
Screen('CloseAll') ;
toc