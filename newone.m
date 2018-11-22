clear all 
close all

%HideCursor
save_name='S3_seq2';

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
    visiblesize=2*texsize+1 ;

    % Create one single static grating image:
    [x,y]=meshgrid(-2*texsize:2*texsize + p, -texsize:texsize);   
    gb(:,:,1) = gray + inc*cos(fr*x) ; 
    gb(:,:,2) = white*ones(size(x)) ; 
    gf(:,:,1) = gray + inc*cos(fr*x) ; 
    gf(:,:,2) = gray*ones(size(x)) ; 
    
     % red crosshair (radius 8 pixels, dark red)
    xhair = zeros(size(x,1),size(x,2),4) ;
    xhair(:,:,1) = gray ; 
    xhair(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 10) ; 
    xhairtex = Screen('MakeTexture',w,xhair) ;
    
    % smaller brighter red crosshair
    %xhairsmall = zeros(size(x,1),size(x,2),4) ;
    %xhairsmall(:,:,1) = white ; 
    %xhairsmall(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 25) ; 
    %xhairsmalltex = Screen('MakeTexture',w,xhairsmall) ;
    
    % black crosshair (radius 8 pixels, dark red)
    xhair2 = zeros(size(x,1),size(x,2),4) ;
    xhair2(:,:,1) = black ; xhair2(:,:,2) = black ; xhair2(:,:,3) = black ; 
    xhair2(:,:,4) = white*double(sqrt(x.^2 + y.^2) < 10) ; 
    xhair2tex = Screen('MakeTexture',w,xhair2) ;
    
    % Store grating in texture:
    gtexb=Screen('MakeTexture', w, gb) ;
    gtexf=Screen('MakeTexture', w, gf) ;

    % Query duration of monitor refresh interval:
    ifi=Screen('GetFlipInterval', w) ;
    
    waitframes = 1;
    waitduration = waitframes * ifi ;
    
    % Translate requested speed of the grating (in cycles per second)
    % into a shift value in "pixels per frame", assuming given
    % waitduration: This is the amount of pixels to shift our "aperture" at
    % each redraw:
    shiftperframe= cyclespersecond * p * waitduration ;
    
    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp for our "WaitBlanking" emulation:
    disp('flipping1') ; 
    vbl=Screen('Flip', w) ;
    
    frameRate=Screen('FrameRate',screenNumber) ;
    % Convert movieDuration in seconds to duration in frames to draw:
   
    % Use realtime priority for better timing precision:
    priorityLevel=MaxPriority(w) ;
    Priority(priorityLevel) ;
    blockcounter = 1 ; 
    task = false ;
    rest = true ; 

    % Tell matlab command window to stop listening to KB input
    %ListenChar(2);
allowedkeys = {KbName('space'),KbName('a'),KbName('b'),KbName('c'),KbName('d'),KbName('e'),KbName('f'),KbName('g'), KbName('h'), KbName('i'),KbName('j'),KbName('k'), KbName('l'), KbName('m'), KbName('n'), KbName('o'), KbName('p'),KbName('q'), KbName('r'), KbName('s'),KbName('t'),KbName('u'),KbName('v'),KbName('w'),KbName('x'),KbName('y'),KbName('z')} ;
kc_end = KbName('return'); % Quit on space or return
pressed = {} ;
nblocks = 2 ;
blockseqs = 2 ;

%seq = {'a','d','f','s','d','a','s','f'} ; %4 2 1 3 2 4 3 1

seq = {'s','d','a','f','s','f','a','d'} ; % 3 2 4 1 3 1 4 2
keysWereDown = ones(1,256);
starttime = GetSecs ;
nblocks = 2 ;
blockseqs = 2 ;
% init
pressed_key=cell(0);
time_key=cell(0);
time_correct_sequence=cell(0);
at_what_time_correct_sequence=cell(0);
Fs=2000;
t=(0:2000)'/Fs;
high=sin(200*2*pi*t);
pause(5) ;
for m=1:nblocks
    ncorrects = 0 ;
    sound(high)
    seqStart = GetSecs ;
    faults = 0 ;
    %seqcount =0 ;
    
    while ncorrects < blockseqs
        % We want keysAreDown, which contains the state of every key
        [keyIsDown ctime keysAreDown] = KbCheck();
        % The newly released keys are ones that were down but now are not
        keysReleased = (keysAreDown==0) & (keysWereDown==1);
        % The newly pressed keys are ones that were not down but now are
        keysPressed  = (keysAreDown==1) & (keysWereDown==0);
        % Check any of the target keys were pressed
        for i=1:size(allowedkeys,2)
            if any(keysPressed(allowedkeys{i}))
                % lptwrite(57424,0) ;
                %lptwrite(57424,allowedkeys{i}) ;
                disp(strcat(KbName(allowedkeys{i}),' _ pressed')) ;
                pressed{pcount+1} = KbName(allowedkeys{i}) ;
                pressed_key{end+1,m} = KbName(allowedkeys{i}) ;
                time_key{end+1,m}=GetSecs-seqStart;
                pcount = pcount + 1 ;
            end
           vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
           Screen('DrawTexture', w, xhair2tex);        
        end
         vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
        Screen('DrawTexture', w, xhairtex);
        % scan the current sequence
        if size(pressed,2) >= 1
          %disp('checking sequence') ;
            for x=1:size(pressed,2)
                
                if pressed{x} ~= seq{x} % if wrong
                    beep
                    pressed = {} ;
                    pcount = 0 ;
                    disp('incorrect') ;
                    faults = faults +  1 ;
                    break ;
                    
                end
            end
             if any(keysPressed(kc_end))
            disp('return pressed: end')
            breaktime = true ;
            break
        end
        keysWereDown = keysAreDown;
    end
   
    allfaults(m) = faults ;
    timex(m) = GetSecs - seqStart ;    
    sound(high)
    vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
     Screen('DrawTexture', w, xhairtex); 
      pause(12) ; % pause at the end of each block
     
   
end
end
figure;bar(timex) ;
figure;bar(cell2mat(time_correct_sequence)')
figure;bar(diff(cell2mat(at_what_time_correct_sequence))')

save(save_name,'pressed_key','time_key','time_correct_sequence','timex','at_what_time_correct_sequence')

    