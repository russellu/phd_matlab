clear all ; close all
%HideCursor

stimtime = 3 ; % stimulus time in seconds 
resttime = 0 ; % post stimulus time 
focustime = 0 ; % minimum focus time (+ 1-2s jitter)
nparams = 8 ; % 6 orientations for the starting angle, and the 2 annulus types
ntrials = 4 ;

params = 1:nparams ;
randparams = [1,1,1,1] ; % initialization vector
%while (min(abs(diff(randparams))) < 2) | (min(abs(diff(randparams(1:2:end)))) < 2) | (min(abs(diff(randparams(2:2:end)))) < 2)
for i=1:10000
paramrep = repmat(params,[1,ntrials]) ; 
randinds = randperm(length(paramrep)) ; 
randparams = paramrep(randinds) ; 
allrandparams(i,:) = randparams ; 
end
sumparams = (sum(abs(diff(allrandparams,1,2)),2)) ; maxparam = find(sumparams==max(sumparams),1) ; 
randparams = allrandparams(maxparam,:) ; 
% 10 different retinotopic configurations...40 trials
startangles = [0,60,120,180,240,300] ;
orientations = [0,30,60,90] ; 
orientationInds = 1:length(randparams) ; 

for i=1:nparams 
    indsi = find(randparams==i) ; 
    orients(indsi) = randperm(4) ; 
end

for i=1:length(randparams) ; triggers(i) = str2num([num2str(randparams(i)),num2str(orients(i))]) ; end



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

% circular aperture mask 
circmask(:,:,1) = x.*0 + gray ; 
circmask(:,:,2) = double(sqrt(x.^2 + y.^2) > 400)*white ;
circmasktex = Screen('MakeTexture',w,circmask) ; 

fnbhd = 15 ; 
% quadrants 
[x,y] = meshgrid(-800:800,-800:800) ; 
botleft(:,:,2) = (1-((x<0) .* (y>0)))*255 ; botleft(:,:,1) = x.*0 + gray ; 
botleft(:,:,2) = imfilter(squeeze(botleft(:,:,2)),fspecial('gaussian',fnbhd,fnbhd)) ; 
quadtex(1) = Screen('MakeTexture',w,botleft) ; 

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
restframes = resttime * frameRate ; 
focusframes = focustime * frameRate ; 

% construct the circular annulus
[x,y] = meshgrid(-screenWidth/4:2:screenWidth/4,-screenHeight/4:2:screenHeight/4) ; 
[theta,rho] = cart2pol(x,y) ; 
incr = (log(stimframes))/stimframes ; 
outer = exp(1:incr:(log(stimframes)+1)) ;
fovthresh = 30 ; 
rhos = zeros(size(x,1),size(x,2),stimframes) ; 
for i=1:length(outer)
    if i < fovthresh+1
        rhos(:,:,i) = rho < outer(i) ; 
    else
        rhos(:,:,i) = rho < outer(i) & rho > outer(i-fovthresh) ;     
    end
end

for i=1:size(rhos,3)
    rhoarr(:,:,2) = (1-squeeze(rhos(:,:,i)))*255 ; rhoarr(:,:,1) = x.*0 + gray ; 
    rhotex(i) = Screen('MakeTexture',w,rhoarr) ; 
end


startT = getSecs ; 

% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w) ;
Priority(priorityLevel) ;    
rotincr = 365/stimframes ; 
rots = round(0:rotincr:365) ;     
tic    
currentstate = 0 ; % currentstate = 0 => rest, currentstate = 1 => focus fixation, currentstate = 2 => stimulus on    
for presentation = 1:length(triggers)
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
                disp([num2str(randparams(presentation)),' ',num2str(orientations(orients(presentation))),' ',num2str(orients(presentation))]) ; 
            end
        end            
        if currentstate == 2 % present the stimulus  
            if randparams(presentation) < 7
                srcRect=[0 0 0 + 800 800];           
                Screen('DrawTexture',w,sgtex,srcRect,[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],orientations(orients(presentation))) ;   
                Screen('DrawTexture',w,quadtex(1),[],[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],rots(stimcount)+startangles(randparams(presentation))) ; 
            else
                if randparams(presentation) == 7
                    srcRect=[0 0 0 + 800 800];           
                    Screen('DrawTexture',w,sgtex,srcRect,[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],orientations(orients(presentation))) ;   
                    Screen('DrawTexture',w,rhotex(stimcount),[],[0,0,screenWidth,screenHeight]) ;
                elseif randparams(presentation) == 8 
                    srcRect=[0 0 0 + 800 800] ;           
                    Screen('DrawTexture',w,sgtex,srcRect,[-screenWidth,-screenHeight,screenWidth*2,screenHeight*2],orientations(orients(presentation))) ;   
                    Screen('DrawTexture',w,rhotex(length(rhotex)-stimcount),[],[0,0,screenWidth,screenHeight]) ;
                end
            end
            %Screen('DrawTexture',w,circmasktex) ;                 
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