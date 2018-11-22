clear all ; close all

[tetanic,freq1] = audioread('C:\shared\sounds\Tetanic_13Hz_1000Hz_-6dBFS_0.05s.aiff');  
[hz750,~] = audioread('C:\shared\sounds\750hz_sinwave_50ms.aiff');  
[hz1000,~] = audioread('C:\shared\sounds\Single_sound_1000Hz_-6dBFS_0.05s.aiff');  

tetanic(end:end+1100) = 0; 

nrchannels = 2; 
sound_tetanic = [tetanic,tetanic]'; 
sound_750 = [hz750,hz750]'; 
sound_1000 = [hz1000,hz1000]'; 

nstims_pre = 3; 
nstims_post = 3;

%soundsc(wavedata,freq)

stimtime = 2 ; % stimulus time in seconds
resttime = 1.75 ; % post stimulus time
tetanictime = 15; 
focustime = 0 ; % minimum focus time (+ 1-2s jitter)

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
                %lptwrite(57424,randparams(presentation)) ;
                %WaitSecs(0.004) ;
                %lptwrite(57424,0) ;
            end
        end            
        if currentstate == 2 % present the stimulus       
            

            if (randparams(presentation)==1 || randparams(presentation)==11) && stimcount==1
                soundsc(sound_750,freq1)
            elseif (randparams(presentation)==2 || randparams(presentation)==12) && stimcount==1
                soundsc(sound_1000,freq1)
            elseif randparams(presentation)==3 && (stimcount==1 || mod(stimcount,round(0.9746*frameRate))==0)
                soundsc(sound_tetanic,freq1)
            end
            
            Screen('DrawTexture',w,xhairtex) ;  
            Screen('DrawTexture',w,xhairsmalltex) ;
            stimcount = stimcount + 1 ; 
            if randparams(presentation)==3
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
%lptwrite(57424,0) ;%Initiate triggers
%lptwrite(57424,0) ;
%lptwrite(57424,99) ;% start trials
pause(10) ;
Screen('CloseAll') ;
toc