clear all ; close all

% cryocooler => sequence dependent

stimtime = 1.75 ; % stimulus time in seconds
resttime = 1.75 ; % post stimulus time
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

HideCursor
Screen('Preference', 'SkipSyncTests', 1) ;
AssertOpenGL ;
screens=Screen('Screens') ;
screenNumber=max(screens) ;

trigsON = false ; 

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
trigs = [31,32,33,34,35,1,2,3,4,5,14,15,16,17,19,20] ; 
prots = [90,180] ; % arots(end) = [] ; 
randparams = repmat(trigs,[1,10]) ; 
perm = randperm(length(randparams)) ; 
randparams = randparams(perm) ; 

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

% make the retinotopic stimuli:
[xg,yg] = meshgrid(-150:150,-150:150) ; 
width = length(xg) ; height = length(yg) ; 
circmask = double(sqrt(xg.^2 + yg.^2) < 350);
fullfield = sqrt(xg.^2 + yg.^2) < 150 ; fullfield = fullfield.*circmask ;
left = (xg<0) .* circmask ; right = (xg>=0) .* circmask ; top = (yg<=0).*circmask ; bot = (yg>0).*circmask ; 
[th,rh] = cart2pol(xg,yg) ; 
wedge1 = (th>0 & th < pi/4).*circmask ; wedge2 = (th >= pi/4 & th < pi/2).*circmask ; 
wedge3 = (th>=pi/2 & th < 3*pi/4).*circmask ; wedge4 = (th>=3*pi/4 & th <= pi).*circmask ; 
wedge5 = (th<=0 & th < -(pi-pi/4)).*circmask ; wedge6 = (th >= -(pi-pi/4) & th < -(pi-pi/2)).*circmask ;  
wedge7 = (th >= -(pi-pi/2) & th < -(pi-3*pi/4)).*circmask ; wedge8 = (th >= -(pi-3*pi/4) & th <= 0).*circmask ; 
allwedge(:,:,1) = wedge1 ; allwedge(:,:,2) = wedge2 ; allwedge(:,:,3) = wedge3 ; allwedge(:,:,4) = wedge4 ; 
allwedge(:,:,5) = wedge5 ; allwedge(:,:,6) = wedge6 ; allwedge(:,:,7) = wedge7 ; allwedge(:,:,8) = wedge8 ; 
ring1 = double(sqrt(xg.^2 + yg.^2) < 400) & double(sqrt(xg.^2 + yg.^2) >=200) ;
ring2 = double(sqrt(xg.^2 + yg.^2) < 200) & double(sqrt(xg.^2 + yg.^2) >=50) ;
ring3 = double(sqrt(xg.^2 + yg.^2) < 50) & double(sqrt(xg.^2 + yg.^2) >=0) ;
allrings(:,:,1) = ring1 ; allrings(:,:,2) = ring2 ; allrings(:,:,3) = ring3 ; 
quad1 = (th>0 & th<pi/2).*circmask ; quad2 = (th>=pi/2 & th<=pi).*circmask ; 
quad3 = (th>=-pi & th<-pi/2).*circmask ; quad4 = (th>=-pi/2 & th<=0).*circmask ; 
allquads(:,:,1) = quad1  ; allquads(:,:,2) = quad2 ; allquads(:,:,3) = quad3 ; allquads(:,:,4) = quad4 ; 
rots = 22.5:22.5:361 ; 
bar1 = (yg<10 & yg>-10).*circmask ; bar2 = imrotate(bar1,rots(1),'nearest','crop').*circmask ; 
bar3 = imrotate(bar1,rots(2),'nearest','crop').*circmask ; bar4 = imrotate(bar1,rots(3),'nearest','crop').*circmask ;
bar5 = imrotate(bar1,rots(4),'nearest','crop').*circmask ; bar6 = imrotate(bar1,rots(5),'nearest','crop').*circmask ;
bar7 = imrotate(bar1,rots(6),'nearest','crop').*circmask ; bar8 = imrotate(bar1,rots(7),'nearest','crop').*circmask ; 
allbars(:,:,1) = bar1 ; allbars(:,:,2) = bar2 ; allbars(:,:,3) = bar3 ; allbars(:,:,4) = bar4 ; 
allbars(:,:,5) = bar5 ; allbars(:,:,6) = bar6 ; allbars(:,:,7) = bar7 ; allbars(:,:,8) = bar8 ; 
clear allstims ; 
allstims(:,:,1) = fullfield ; allstims(:,:,2) = left ; allstims(:,:,3) = right ; allstims(:,:,4) = top ; allstims(:,:,5) = bot ; 
allstims(:,:,6:13) = allwedge ; allstims(:,:,14:17) = allquads ; allstims(:,:,18:20) = allrings ; 
for i=1:8
   allstims(:,:,i+20) = imrotate(allwedge(:,:,i),22.5,'crop') ;  
end

[xg2,y2g] = meshgrid(-20:.5:20,-20:.5:20) ; 
conts = [0.01,0.05,0.15,0.35,1];
for i=1:5
    sinu(:,:,2) = uint8(ones(size(xg2))*255) ; 
    sinu(:,:,1) = uint8(((sin(xg2))*127.5)*conts(i) + 127.5) ; 
    allsinu(:,:,i) = sinu(:,:,1); 
    %sinutext(i) = Screen('MakeTexture',w,sinu) ; 
end

for i=1:5 ; subplot(1,5,i) ; mask = imresize(allstims(:,:,1),[81,81]); img = double(allsinu(:,:,6-i)).*round(imresize(allstims(:,:,1),[81,81])); img(mask<.9) = 127.5; imagesc(img,[0,255]) ; set(gca,'XTickLabel',[],'YTickLabel',[]); colormap gray; end


%{
thplus = th+pi ; 
rightedge = thplus<pi+(pi/2-pi/8) & thplus > pi-(pi/2)+pi/8 ; 
leftedge = fliplr(rightedge) ; 
allstims(:,:,29) = rightedge ; allstims(:,:,30) = leftedge ; 
%}
for i=1:size(allstims,3)
    texmat(:,:,2) = uint8((1-allstims(:,:,i))*255) ; 
    texmat(:,:,1) = gray ; 
    alltex(i) = Screen('MakeTexture',w,texmat) ; 
end

% make the grating: 
circmask = double(sqrt(xg.^2 + yg.^2) < 70);
%texmat(:,:,2) = uint8(mat2gray(1-circmask)*255) ; 
%texmat(:,:,1) = gray ; 
%circtex = Screen('MakeTexture',w,texmat) ; 
clear rhosin
for i=1:200 
    [xg,yg] = meshgrid(-50:50,-50:50) ; 
    [th,rh] = cart2pol(xg,yg) ;
    rhosin = sin(rh+i/16) ;  
    %tmat(:,:,2) = uint8(zeros(size(xg))+255) ; 
    %tmat(:,:,1) = uint8((rhosin)*255) ; 
    %gratetex(i) = Screen('MakeTexture',w,tmat) ; 
end


titles = {'full','left','right','top','bottom','bottom right','bottom left','top left','top right','periphery','fovea'};
retstims = [1,2,3,4,5,14,15,16,17,19,20]; 
for i=1:length(retstims) ; subplot(1,11,i) ; maski = imresize(squeeze(allstims(:,:,retstims(i))).*allstims(:,:,1),[101,101]) ; imagesc(rhosin.*maski); colormap gray; set(gca,'XTickLabel',[],'YTickLabel',[]);  end


% make the plaids: 
[xg,yg] = meshgrid(-350:.25:350,-50:.25:50) ; 
sinxg = sin(xg) ; 
binsin = sinxg <0 ; 
clear tmat ; 
tmat(:,:,2) = uint8(zeros(size(xg))+255) ; 
tmat(:,:,1) = uint8(mat2gray(binsin)*85) ; 
plaidtex = Screen('MakeTexture',w,tmat) ; 
% make the second plaid
binsin = sinxg<0 ; 
clear tmat ; 
tmat(:,:,2) = uint8(zeros(size(xg))+255) ; 
tm2 = tmat(:,:,2) ; tm2(binsin==0) = 0 ; tmat(:,:,2) = tm2 ; 
tmat(:,:,1) = uint8(mat2gray(binsin)*190) ; 
plaidtex2 = Screen('MakeTexture',w,tmat) ; 

% make the grating itself
[xg,yg] = meshgrid(-350:.25:350,-50:.25:50) ; 
sinxg = sin(xg) ; 
clear tmat ; 
tmat(:,:,2) = uint8(ones(size(xg))*255) ; 
tmat(:,:,1) = uint8((sinxg)*255) ; 
sintex = Screen('MakeTexture',w,tmat) ; 
newrots = 0:11.25:180 ; newrots(end) = [] ; 

startT = GetSecs ; 
tic
currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
for presentation = 1:length(randparams)
    restcount = 1 ; 
    focuscount = 1 ;
    stimcount = 1 ;
    poststimcount = 1 ;
    jitter = 0 ; %(resttime*frameRate) + (rand-.5)*frameRate ; 
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
            if focuscount > focusframes %+ jitter
                currentstate = 2 ;
                stimtimes(1,presentation) = GetSecs - startT ;
                lptwrite(57424,randparams(presentation)) ;
                WaitSecs(0.004) ;
                lptwrite(57424,0) ;
            end
        end            
        if currentstate == 2 % present the stimulus               
            xoffset = mod(stimcount*shiftperframe,360)/2;  % Shift the grating by "shiftperframe" pixels per frame:
            srcRect=[0+xoffset/4 0 size(xg,1)+xoffset/4 size(xg,2)];   
            sqrRect = [(screenRect(3)-screenRect(4))/2, 0, screenRect(4)+(screenRect(3)-screenRect(4))/2, screenRect(4)] ;
            largeRect = [sqrRect(1)-500,sqrRect(2)-500,sqrRect(3)+500,sqrRect(4)+500] ; 

            if randparams(presentation)<=30
                %Screen('DrawTexture',w,gratetex(mod(stimcount,51)+1),[],sqrRect) ;
                Screen('DrawTexture',w,gratetex(stimcount),[],sqrRect) ;
                Screen('DrawTexture',w,alltex(randparams(presentation)),[],sqrRect) ; 
            %elseif randparams(presentation) < 40
            %    Screen('DrawTexture',w,plaidtex,srcRect,sqrRect) ; 
            %    Screen('DrawTexture',w,plaidtex2,srcRect,sqrRect,prots(randparams(presentation)-20)) ;        
            else
                %Screen('DrawTexture',w,sintex,srcRect,sqrRect,newrots(randparams(presentation)-40)) ;                
                Screen('DrawTexture',w,sinutext(randparams(presentation)-30),srcRect,sqrRect) ;      
            end
            
            Screen('DrawTexture',w,circtex,[],largeRect) ; 
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
        if KbCheck % Abort demo if any key is pressed:
           Screen('CloseAll');
           break
        end
    end
end
%save(name,'stimtimes') ; 
vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw. 
lptwrite(57424,0) ;%Initiate triggers
lptwrite(57424,0) ;
lptwrite(57424,99) ;% start trials
pause(10) ;
Screen('CloseAll') ;
toc