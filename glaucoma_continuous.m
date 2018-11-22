clear all ; close all
sub = 'laurent' ; 
%name = [sub,'_stim_left'] ; 
name = [sub,'_stim_left'] ; 

%handle = CedrusResponseBox('Open','COM3') ; 

stimtime = 20 ; % stimulus time in seconds
resttime = 6 ; % post stimulus time
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

HideCursor
Screen('Preference', 'SkipSyncTests', 1);
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

%randparams = repmat([1,3,2,4],[1,4]) ; 
randparams = [13,13,1,2,3,4,5,6,7,8,9,10,11,12,13,13]; 
resttimes = [6,6,6,6,6,6,3,3,3,3,3,3,3,6,6,6] ; 
allrounds = [1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1] ; 
%randparams = ones(1,10)*3 ; 
stimtimes = zeros(2,length(randparams)) ; 
stimtypes = randparams ; 

% create the checkerboard
[x,y] = meshgrid(-800:800,-800:800) ; 
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

sthetas2 = sum(thetas,3) ; 
sthetas2(srhos==0) = ~sthetas(srhos==0) ; 
sthetas2 = uint8(mat2gray(sthetas2)*255) ; 

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

blackb = zeros(3000,3000) ; 
whiteb = ones(3000,3000)*255 ; 
blacktex = Screen('MakeTexture',w,blackb) ;
whitetex = Screen('MakeTexture',w,whiteb) ;

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

fnbhd = 15 ; 
% quadrants 
[x,y] = meshgrid(-1200:1200,-1200:1200) ; 
[th,rh] = cart2pol(x,y) ; 
botleft(:,:,2) = (1-(th<pi/4 & th>0))*255 ; botleft(:,:,1) = x.*0 + gray ; 
botleft(:,:,2) = imfilter(squeeze(botleft(:,:,2)),fspecial('gaussian',fnbhd,fnbhd)) ; 
quadtex(1) = Screen('MakeTexture',w,botleft) ; 

% circle mask
circ = zeros(size(x,1),size(x,2),2) ;
circ(:,:,2) = (sqrt(x.^2 + y.^2) > 600).*white ; 
circ(:,:,1) = gray ; 
circtex = Screen('MakeTexture',w,circ) ;


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

for i=1:length(resttimes)
    restframes(i) = resttimes(i) * frameRate ;
end

focusframes = focustime * frameRate ;
% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w) ;
Priority(priorityLevel) ;

hzs = 8 ; 
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

% construct the circular annulus MASKS for radial mapping
[x,y] = meshgrid(-screenWidth/4:2:screenWidth/4,-screenHeight/4:2:screenHeight/4) ; 
[theta,rho] = cart2pol(x,y) ; 
incr = (log(stimframes))/(stimframes) ; 
outer = exp(1:incr:(log(stimframes)+1)) ;

startind = round(length(outer))/4 ; 
outer2 = outer(startind:end) ; 
outer=  (imresize(outer2,[1,length(outer)])) ;
% if the annulus is too big for the screen
%outer = outer/2 ;
minouter = min([screenWidth,screenHeight]) ; 
ind = find(abs(outer-minouter) == min(abs(outer-minouter))) ; 
newouter = outer(1:ind) ; 
outer = imresize(newouter,[1,length(outer)]) ; 
outer = outer / 2 ; 
fovthresh = 80 ; 
rhos = zeros(size(x,1),size(x,2),stimframes) ; 
for i=1:length(outer)
    d4 = sqrt(outer(i))*2 ;
    rhos(:,:,i) = rho < outer(i) & rho > outer(i)-d4 ;     
end
for i=1:size(rhos,3)
    rhoarr(:,:,2) = (1-squeeze(rhos(:,:,i)))*255 ; rhoarr(:,:,1) = x.*0 + gray ; 
    rhotex(i) = Screen('MakeTexture',w,rhoarr) ; 
end
rotincr = 365/stimframes ; 
rots = round(0:rotincr:365) ;    

% moving bar stimulus
[xg,yg] = meshgrid(-300:300) ; 
barframes = stimframes/4 ; 
barincr = size(xg,1)./barframes ; 
barcount = 1 ; 
for i=-150:barincr/2:150
    width = (abs(i))/4+5 ; 
    bars(:,:,barcount) = xg>i & xg <i+width ;  
    barcount = barcount + 1 ; 
end
for i=1:size(bars,3)
    bararr(:,:,2) = (1-squeeze(bars(:,:,i)))*255 ; bararr(:,:,1) = xg.*0 + gray ; 
    bartex(i) = Screen('MakeTexture',w,bararr) ; 
end




stimtypeframes = [stimframes/8,stimframes/8,...
                  stimframes,stimframes,stimframes,stimframes,...
                  stimframes/4,stimframes/4,stimframes/4,stimframes/4,stimframes/4,stimframes/4,stimframes/4,stimframes/4,...
                  stimframes/8,stimframes/8] ; 
%{
CedrusResponseBox('FlushEvents', handle);
evt = CedrusResponseBox('GetButtons', handle) ;
disp('waiting...') ;
while isempty(evt) & ~KbCheck
    evt = CedrusResponseBox('GetButtons', handle) ;
    resp = CedrusResponseBox('FlushEvents', handle)  ;
    if resp==-1
        break ;
    end
end       
%}
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
            Screen('DrawTexture',w,xhair2tex) ;  
            if restcount > restframes(presentation)
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
            end
        end            
        if currentstate == 2 % present the stimulus               
            xoffset = 0 ;%mod(stimcount*shiftperframe,360);  % Shift the grating by "shiftperframe" pixels per frame:
            srcRect=[xoffset 0 xoffset + visiblesize visiblesize];   
            if randparams(presentation)==1 || randparams(presentation) == 2 || randparams(presentation) >= 5
                if stimarr(stimcount) == 0
                    Screen('DrawTexture',w,checktex1,[],[]) ;                
                else
                    Screen('DrawTexture',w,checktex2,[],[]) ;
                end
            elseif randparams(presentation)==3 || randparams(presentation) == 4 
                if stimarr(stimcount) == 0
                    Screen('DrawTexture',w,thetatex1,[],[]) ;                
                else
                    Screen('DrawTexture',w,thetatex2,[],[]) ;
                end
            end          
            if randparams(presentation)==1 
                Screen('DrawTexture',w,quadtex(1),[],[],rots(stimcount)) ; 
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ;  
                Screen('DrawTexture',w,circtex) ;                                        
            elseif randparams(presentation)==2
                Screen('DrawTexture',w,quadtex(1),[],[],rots(length(rots)-stimcount)) ; 
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
                Screen('DrawTexture',w,circtex) ;                                        
            elseif randparams(presentation)==3
                Screen('DrawTexture',w,rhotex(length(rhotex)-stimcount),[],[0,0,screenWidth,screenHeight]) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
                Screen('DrawTexture',w,circtex) ;                                        
            elseif randparams(presentation)==4
                Screen('DrawTexture',w,rhotex(stimcount),[],[0,0,screenWidth,screenHeight]) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
                Screen('DrawTexture',w,circtex) ;                                        
            elseif randparams(presentation)==5
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],0) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
            elseif randparams(presentation)==6
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],180) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
            elseif randparams(presentation)==7
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],90) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
            elseif randparams(presentation)==8
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],270) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
            elseif randparams(presentation)==9
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],45) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
            elseif randparams(presentation)==10
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],135) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
            elseif randparams(presentation)==11
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],315) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
            elseif randparams(presentation)==12
                Screen('DrawTexture',w,bartex(stimcount),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],225) ;
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
           elseif randparams(presentation)==13
                Screen('DrawTexture',w,xhairtex) ;  
                Screen('DrawTexture',w,xhairsmalltex) ; 
                Screen('DrawTexture',w,circtex) ;                                        
            end 
            stimcount = stimcount + 1 ; 
            if stimcount > stimtypeframes(presentation) 
                currentround = currentround + 1 ; 
                stimcount = 1 ; 
                if currentround >= allrounds(presentation) 
                    currentstate = 3 ;         
                    stimtimes(2,presentation) = GetSecs - startT ; 
                end
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
%The same commands wich close onscreen and offscreen windows also close
%textures.
save(name,'stimtimes') ; 
vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.      
pause(10) ; 
Screen('CloseAll');
toc


