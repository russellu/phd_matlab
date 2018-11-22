%clean
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %  Define experimental parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Screen angular size : TBD from experimental setup)
% sizeDegX = 24; % 2*atand(halfScreenSizeX/screenDist
sizeDegX = 24; % 2*atand(halfScreenSizeX/screenDist
sizeDegY = 18; % 2*atand(halfScreenSizeY/screenDist

% DEBUG : scales stimDuration, dot speed, etc. because image display is slow in matlab figures
scalingFactor = 1;
HideCursor
% Dot definition
nDots = 1000;
nBlackDots = nDots/2;
nWhiteDots = nDots - nBlackDots;
dotSize = 4;
dotSpeedDeg = 16; % Dot angular speed (deg/s)
dotLifetime = 0.3; % in (s)
dotAtCenterFraction = 1/20; % Fraction of dots outside the screen regenerated close to the center

dotDegSpeed = dotSpeedDeg/scalingFactor;
dotLifetime = dotLifetime*scalingFactor;

% Fization dot definition
fixDotRadiusDeg = 0.5;
fixDotColor = [1 0 0];

% Stimulus/rest definition
stimDur=8;
nReps = 2;
stimDur=stimDur *scalingFactor;
restDur=5;

stimFreq = 1/(stimDur);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %  Main program
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try % Avoid getting stuck with an uncloseable Psychtoolbox window
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %  Setup Pschtoolbox variables
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sca
    
    % Here we call some default settings for setting up Psychtoolbox
    PsychDefaultSetup(2);
    
    % Get the screen numbers. This gives us a number for each of the screens
    % attached to our computer.
    screens = Screen('Screens');
    
    % To draw we select the maximum of these numbers. So in a situation where we
    % have two screens attached to our monitor we will draw to the external
    % screen.
    screenNumber = max(screens);
    
    % Define black, white and grey
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    grey = white / 2;
    
    % dotSize x dotSize rectangle for display
    dotTemplate = [0 0 dotSize dotSize];
    
    % Open an on screen window
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
    
    % Get the size of the on screen window
    [sizePixX, sizePixY] = Screen('WindowSize', window);
    
    % Get the centre coordinate of the window
    [centerPixX, centerPixY] = RectCenter(windowRect);
    
    % Measure the vertical refresh rate of the monitorc
    dt = Screen('GetFlipInterval', window);
    
    % Numer of frames to wait when specifying good timing
    waitframes = 1;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %  Setup computation
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Screen
    pixPerDeg = sizePixX/sizeDegX;
    
    % Fixation dot
    fixDotRadiusPix = fixDotRadiusDeg*pixPerDeg;
    fixDotRadiusPix = floor(fixDotRadiusPix/2);
    fixDotPosCart = [centerPixX, centerPixY];
    
    % Time and radius increment at each frame
    
    %     dt = dt*scalingFactor;
    dr = pixPerDeg * dotSpeedDeg * dt; % pix/deg * deg/s * s/increment = pix/increment
    
    tVect = 0:dt:(stimDur-dt);
    nFramesStim = numel(tVect);
    nFramesRest = round(restDur/dt);
    
    % Define useful functions
    genDotPosCart = @(N) [((sizePixX-1)*rand(N,1) + 1 - centerPixX) , ...
        ((sizePixY-1)*rand(N,1) + 1 - centerPixY)];
    maskCount1D = @(x) sum(logical(x)); % Count number of 1's in 1D mask
    
    % Same as cart2pol and pol2cart, but input and output are in single vectors [x, y] and [th, r]
    cartesian2polar = @(x_y) [atan2(x_y(:,2),x_y(:,1)), hypot(x_y(:,1),x_y(:,2))];
    polar2cartesian = @(th_r) [th_r(:,2).*cos(th_r(:,1)), th_r(:,2).*sin(th_r(:,1))];
    moveDots = @(x, dx) bsxfun(@plus, x, dx);
    
    % Pre-compute dot color to speed-up image formation
    greySinMod = grey*(-0.5*cos(2*pi*tVect*stimFreq)+0.5);
    blackSinMod = grey - greySinMod;
    whiteSinMod = grey + greySinMod;
    dotColor = [repmat(blackSinMod,[nBlackDots 1]); repmat(whiteSinMod,[nWhiteDots 1])];
    % tmp = rand(size(dotColor,1)); % This could make indexing faster... (rnd
    % color, vs rand age --> could improve indexing speed with deadMask)
    % [~, rndInd] = sort(tmp);
    % dotColor= dotColor(rndInd,:);
    

    % Retreive the maximum priority number
    topPriorityLevel = MaxPriority(window);
    
    % Finally we do the same as the last example except now we additionally
    % tell PTB that no more drawing commands will be given between coloring the
    % screen and the flip command. This, under some circumstances, can help
    % acheive good timing.
    Priority(topPriorityLevel);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %  Main loop : Stimulus presentation
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tStimEnd = zeros(nReps,1);
    tRestEnd = zeros(nReps,1);
    
    vbl = Screen('Flip', window);
    tStart = tic;
    Screen('DrawDots', window, fixDotPosCart, fixDotRadiusPix, fixDotColor, [], 2);
    for rep = 1:nReps
        % Initialize points
        dotPosCart = genDotPosCart(nDots);
        dotPosPol = cartesian2polar(dotPosCart);
        % age = ((0:nDots-1)/nDots*dotLifetime)' + dotLifetime;
        age = ((randperm(nDots)-1)/nDots*dotLifetime)'; % Dot age is randomly initialized between 0 and dotLifetime-dt
        
        %     figure;
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %  Sinusoidal stimulus presentation
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for frame = 1:nFramesStim
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% %  Compute dot position at each frame
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Main position is kept in polar coordinates, since radial movement is simpler. New dot
            % positions are generated in cartesian coordinates (to be uniformly distributed in the screen)
            % and then converted to polar coordinates. Before display, coordinates are converted to
            % cartesian and rounded.
            
            agePrev = age;
            
            % Find dead dots
            deadMask = age >= dotLifetime;
            nDeadDots = maskCount1D(deadMask);
            
            % Update age
            age(deadMask) = agePrev(deadMask) - dotLifetime;
            
            %     % Regenerate dead dots randomly in the image (x,y)
            dotPosCart(deadMask,:) = genDotPosCart(nDeadDots);
            
            % Convert to polar coordinates (theta, r)
            dotPosPol(deadMask,:) = cartesian2polar(dotPosCart(deadMask,:));
            
            % Move dots a distance of dr along the direction they were headed in
            if nDeadDots < nDots
                dotPosPol(~deadMask,:) = moveDots(dotPosPol(~deadMask,:), [0, dr]);
            end
            
            %     % New dots only move a fraction of dr
            %     if nDeadDots > 0
            %         moveVect = bsxfun(@times, [0, dr], age(deadMask)/dt);
            %         dotPosPol(deadMask,:) = moveDots(dotPosPol(deadMask,:), moveVect);
            %     end
            
            % Convert to x-y coordinates for display
            dotPosCart = polar2cartesian(dotPosPol);            
            dotPosCartScreen = bsxfun(@plus, dotPosCart, [centerPixX centerPixY]);
            
            % Find pts outside the image
            outMask = any(bsxfun(@gt, dotPosCartScreen, [sizePixX sizePixY]), 2) | ...
                any(dotPosCartScreen<1, 2);
            nOutVox = maskCount1D(outMask);
            
            % Regenerate dots outside of screen (forget movement... could become a recursive problem).
            % Divide x and y coordinates by 10 for dotAtCenterFraction of new dots
            dotPosCart(outMask,:) = genDotPosCart(nOutVox);
%             ind = randperm(nOutVox, round(nOutVox*dotAtCenterFraction));
            ind = randperm(nOutVox);
            ind = ind(1:round(nOutVox*dotAtCenterFraction));
            tmp = dotPosCart(outMask,:);
            tmp(ind,:) = tmp(ind,:)./10;
            dotPosCart(outMask,:) = tmp;
            dotPosPol(outMask,:) = cartesian2polar(dotPosCart(outMask,:));
            dotPosCartScreen(outMask,:) = bsxfun(@plus, dotPosCart(outMask,:), [centerPixX centerPixY]);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% %  Display frame
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            allRects = [dotPosCartScreen' - dotSize/2 ; dotPosCartScreen' + dotSize/2];
            
            % Draw the rect to the screen
            Screen('FillRect', window, repmat(dotColor(:,frame)', [3 1]), allRects);
            
            Screen('DrawDots', window, fixDotPosCart, fixDotRadiusPix, fixDotColor, [], 2);
            
            % Tell PTB no more drawing commands will be issued until the next flip
            Screen('DrawingFinished', window);
            
            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * dt);
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% %  Old
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         % Build image
            %         im = grey * ones(sizePixY,sizePixX);
            %         validInd = sub2ind([sizePixY, sizePixX], dotPosCartScreen(:,2), dotPosCartScreen(:,1));
            %         % Crude way to make larger (5x5) dots
            %         for kk = [-2*sizePixY-2 -2*sizePixY-1 -2*sizePixY -2*sizePixY+1 -2*sizePixY+2 ...
            %                 -sizePixY-2 -sizePixY-1 -sizePixY -sizePixY+1 -sizePixY+2 ...
            %                 -2 -1 0 +1 +1 ...
            %                 sizePixY-2 sizePixY-1 sizePixY sizePixY+1 sizePixY+2 ...
            %                 2*sizePixY-2 2*sizePixY-1 2*sizePixY 2*sizePixY+1 2*sizePixY+2]
            %             imInd = validInd + kk;
            %
            %             ind = imInd > 0 & imInd < sizePixX*sizePixY;
            %             im(imInd(ind)) = dotColor(ind,frame);
            %         end
            %         % Fixation dot
            %         centerInd = sub2ind([sizePixY, sizePixX], round(centerPixY),round(centerPixX));
            %
            %         for kk = -fixDotRadiusPix:fixDotRadiusPix
            %             ind = centerInd+kk*sizePixY;
            %             im(ind-fixDotRadiusPix:ind+fixDotRadiusPix) = white;
            %         end
            %
            %
            %         %   figure
            %         imshow(im);
            %         title(num2str(tVect(frame)));
            %         % Increment time
            %         age = age + dt;
            %         tCurr = toc(tStart);
            %
            %         pause(frame*dt-tCurr)
            
        end
        tStimEnd(rep) = toc(tStart);
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %  Rest presentation
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for frame = 1:nFramesRest
            Screen('DrawDots', window, fixDotPosCart, fixDotRadiusPix, [0 0 1], [], 2);
            
            % Tell PTB no more drawing commands will be issued until the next flip
            Screen('DrawingFinished', window);
            
            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * dt);
        end
        tRestEnd(rep) = toc(tStart);
    end

    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %  Finish & clean-up
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    repVect = (1:nReps)';
    tStimEndExpected = repVect*stimDur + (repVect-1)*restDur;
    tRestEndExpected = repVect*(stimDur + restDur);
    
    fprintf('tStimEnd tRestEnd tStimEndExpected tRestEndExpected\n')
    disp([tStimEnd tRestEnd tStimEndExpected tRestEndExpected])
    
    Priority(0);
    sca
    
catch
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end