% simple script to display two stimulus types in a block design
% focus_time and post_stim_time are not needed in such a simple paradigm,
% so they are both set to 0. 
% the script cycles through the stimulus blocks defined by the variable
% 'randparams', and displays stimuli in a block design sequence
% anything involving 'cedrusresponsebox' is specific to the MRI scanner
% stimulus computer, and can remain commented-out for now.

clear all ; close all
%handle = CedrusResponseBox('Open','COM3') ; % for MRI scanner

stim_time = 1 ; % stimulus time in seconds
rest_time = 1 ; % post stimulus time
focus_time = 0 ; % minimum focus time (+ 1-2s jitter)
post_stim_time = 0 ; % post-stimulus fixation time 

endpausetime = 1  ; % time to display gray screen at end of script (sec) 

HideCursor % do not show the mouse cursor during stimulus presentation
Screen('Preference', 'SkipSyncTests', 1) ; % skip synchronization tests
AssertOpenGL;
screens=Screen('Screens'); % get all screens on current setup
screenNumber=max(screens); % get the primary screen

white=WhiteIndex(screenNumber) ; % get max luminosity value for screen
black=BlackIndex(screenNumber) ; % get min luminosity value for screen
gray=(white+black)/2 ; % calculate gray as middle luminosity value of screen
if round(gray)==white
    gray=black;
end
inc=white-gray ;	

% Open a double buffered fullscreen window and draw a gray background 
% to front and back buffers:
[w ,screenRect]=Screen('OpenWindow',screenNumber, gray) ;
w2 = screenRect(3)/2 ; h2 = screenRect(4)/2 ; 
screenWidth = screenRect(3) ; screenHeight = screenRect(4) ; % screen dimensions
squareRect = [0,0,screenRect(4),screenRect(4)] ; 
% Enable alpha-blending:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
% Show initial gray screen:
Screen('Flip', w) ;

% randparams is a vector that determines the presentation sequence of the
% different stimuli that you have constructed (for this script, there are
% only two different types). you can also randomize the presentation order
% by using the function 'randperm' 
randparams = [1,2,1,2] ; % (1=high contrast grating, 2=low contrast grating)

% make the stimuli:
[xg,yg] = meshgrid(-150:150,-150:150) ; 

% build the high contrast grating texture
grate_tex_mat(:,:,1) = uint8(mat2gray(sin(xg))*255) ; % texture pattern channel (all values between 0-255)
grate_tex_mat(:,:,2) = 255 ; % opacity channel
grate_tex = Screen('MakeTexture',w,grate_tex_mat) ; % generate the texture

% build the low contrast grating texture
low_grate_tex_mat(:,:,1) = uint8(mat2gray(sin(xg))*255/3 + 85) ; % texture pattern channel (all values between 0-255)
low_grate_tex_mat(:,:,2) = 255 ; % opacity channel
low_grate_tex = Screen('MakeTexture',w,low_grate_tex_mat) ; % generate the texture

% build the circular mask texture
circmask = double(sqrt(xg.^2 + yg.^2) < 150) ; % define the radius of the circular mask 
circtexmat(:,:,2) = uint8(mat2gray(1-circmask)*255) ; % opacity (anything within circular mask is see-through)
circtexmat(:,:,1) = gray ; % texture pattern channel (simple gray background)
circtex = Screen('MakeTexture',w,circtexmat) ; % generate the texture

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

% Perform initial Flip to sync us to the VBL and for getting an initial
% VBL-Timestamp for our "WaitBlanking" emulation:
vbl=Screen('Flip', w) ;

frameRate=Screen('FrameRate',screenNumber) ; % get the frameRate of the monitor
% calculate the number of frames for each 'state' (rest,stim,focus,post)
rest_frames = rest_time * frameRate ; 
stim_frames = stim_time * frameRate ; 
focus_frames = focus_time * frameRate ;
post_stim_frames = post_stim_time * frameRate ; 

% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w) ;
Priority(priorityLevel) ;

%{
% for MRI scanner
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

startT = GetSecs ; % time at which the stimulus started (starts timing after all the texture generation and other setup is complete)
tic
currentstate = 0 ; % currentstate=0 => rest, currentstate=1 => focus fixation, currentstate=2 => stimulus on, currentstate=3 => post-stimulus fixation    
% start the main loop that will cycle through all blocks (number of blocks defined by the length of the 'randparams' variable)
for presentation = 1:length(randparams) % for all stimulus blocks 
    restcount = 1 ; % count (in frames) of the current block's resting period
    focuscount = 1 ; % count (in frames) of the current block's focus period
    stimcount = 1 ; % count (in frames) of the current block's stimulus period
    poststimcount = 1 ; % count (in frames) of the current block's post-stimulus period
    jitter = 0 ; % random jitter (in frames) to add to the current block's focus period
    DONE = false ; 
    while ~DONE % while not DONE, update the screen every frame, and wait at the end if the update took less time than the frame duration            
        if currentstate == 0 % present the black crosshair (rest period)
            restcount = restcount + 1 ; % increment rest count once per frame (will switch to 1 once all rest frames have passed)   
            Screen('DrawTexture',w,xhair2tex) ;  % draw the black crosshair
            if restcount > rest_frames % if the current frame count is greater than the number of rest frames, switch states
                currentstate = 1 ; % switch to the pre-stimulus focus state
            end
        end         
        if currentstate == 1 % present the dark red and bright red crosshairs (focus period)
            focuscount = focuscount + 1 ; % increment focus count once per frame (will switch to 2 once all focus frames have passed)
            Screen('DrawTexture',w,xhairtex) ;  % draw dark red crosshair
            Screen('DrawTexture',w,xhairsmalltex) ; % draw smaller crosshair (bright red)
            if focuscount > focus_frames + jitter % if current frame count is greater than the number of focus frames, switch states
                currentstate = 2 ; % switch to the stimulus state
            end
        end            
        if currentstate == 2 % present the stimulus    
            
            % sqrRect is the portion of the screen on which to display the
            % textures, you can change this to be full screen, or any
            % sub-set of the screen that you wish to display your texture
            % on (the texture will be automatically re-sized depending on
            % the dimensions of 'sqrRect'). right now, it is set to display
            % a square texture in the middle of the screen, which i
            % recommend you leave as is
            sqrRect = [(screenRect(3)-screenRect(4))/2,0,screenRect(4)+(screenRect(3)-screenRect(4))/2,screenRect(4)] ; 
           
            % draw the sinusoidal grating texture in the middle of the
            % screen
            if randparams(presentation)==1 % present the first stimulus type (high contrast grating)
                Screen('DrawTexture',w,grate_tex,[],sqrRect) ; 
                
            elseif randparams(presentation)==2 % present the second stimulus type (low contrast grating) 
                Screen('DrawTexture',w,low_grate_tex,[],sqrRect) ; 

            elseif randparams(presentation)==3 % present the third stimulus type
                % add another stimulus type here
            end
            
            
            
            % draw the circular mask texture in the middle of the screen
            Screen('DrawTexture',w,circtex,[],sqrRect) ; 
            
            % draw the crosshair
            Screen('DrawTexture',w,xhairtex) ; % dark red
            Screen('DrawTexture',w,xhairsmalltex) ; % bright red
            stimcount = stimcount + 1 ; % increment the stimulus frame counter
            if stimcount > stim_frames % if current frame count is greater than the number of stimulus frames, switch states
                currentstate = 3 ; % switch to post-stimulus fixation period        
            end
        end              
        if currentstate == 3 % present the post stimulus fixation
            poststimcount = poststimcount + 1 ; % increment the post-stimulus frame counter
            Screen('DrawTexture',w,xhairtex) ; % dark red crosshair
            Screen('DrawTexture',w,xhairsmalltex) ; % bright red crosshair 
            if poststimcount > post_stim_frames % if current frame count is greater than the number of post-stimulus fixation frames
               currentstate = 0 ; % reset the state
               DONE = true ; % move to the next block 
            end    
        end         
        vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.           
        if KbCheck % Abort demo if any key is pressed:
           Screen('CloseAll');
           break
        end
    end
end
vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);% Flip 'waitframes' monitor refresh intervals after last redraw.      
pause(endpausetime) ; % show a gray screen for endpausetime seconds
Screen('CloseAll'); % end the onScreenWindow process
toc


 