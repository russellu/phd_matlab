clear all ; close all ;
[hz1000,freq] = audioread('C:\shared\sounds\Tetanic_13Hz_1000Hz_-6dBFS_0.05s.aiff');  
nrchannels = 2; 
wavedata = [hz1000,hz1000]'; 
repetitions = 100; 
soundsc(wavedata,freq)



% Perform basic initialization of the sound driver:
InitializePsychSound;

% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of freq and nrchannels sound channels.
% This returns a handle to the audio device:
try
    % Try with the 'freq'uency we wanted:
    pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
catch
    % Failed. Retry with default frequency as suggested by device:
    fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
    fprintf('Sound may sound a bit out of tune, ...\n\n');

    psychlasterror('reset');
    pahandle = PsychPortAudio('Open', [], [], 0, [], nrchannels);
end

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pahandle, wavedata);

% Start audio playback for 'repetitions' repetitions of the sound data,
% start it immediately (0) and wait for the playback to start, return onset
% timestamp.
t1 = PsychPortAudio('Start', pahandle, repetitions, 0, 1);


% Wait for release of all keys on keyboard:
KbReleaseWait;

fprintf('Audio playback started, press any key for about 1 second to quit.\n');

lastSample = 0;
lastTime = t1;

% Stay in a little loop until keypress:
while ~KbCheck
    % Wait a seconds...
    WaitSecs(1);
    
    % Query current playback status and print it to the Matlab window:
    s = PsychPortAudio('GetStatus', pahandle);
    % tHost = GetSecs;
    
    % Print it:
    fprintf('\n\nAudio playback started, press any key for about 1 second to quit.\n');
    fprintf('This is some status output of PsychPortAudio:\n');
    disp(s);
    
    realSampleRate = (s.ElapsedOutSamples - lastSample) / (s.CurrentStreamTime - lastTime);
    fprintf('Measured average samplerate Hz: %f\n', realSampleRate);
    
    tHost = s.CurrentStreamTime;
    clockDelta = (s.ElapsedOutSamples / s.SampleRate) - (tHost - t1);
    clockRatio = (s.ElapsedOutSamples / s.SampleRate) / (tHost - t1);
    fprintf('Delta between audio hw clock and host clock: %f msecs. Ratio %f.\n', 1000 * clockDelta, clockRatio);
end

% Stop playback:
PsychPortAudio('Stop', pahandle);

% Close the audio device:
PsychPortAudio('Close', pahandle);

% Done.
fprintf('Demo finished, bye!\n');
