% Build the path of the files to import
AnatDir   = fullfile(tutorial_dir, 'anatomy');
RawFile   = fullfile(tutorial_dir, 'data', 'tutorial_EEG.bin');
ElcFile   = fullfile(tutorial_dir, 'data', 'tutorial_electrodes.elc');
SpikeFile = fullfile(tutorial_dir, 'data', 'tutorial_spikes.txt');

% Start a new report
bst_report('Start');
% Subject name
SubjectName = 'sepi01';

% === ANATOMY ===
% Process: Import anatomy folder
bst_process('CallProcess', 'process_import_anatomy', [], [], ...
    'subjectname', SubjectName, ...
    'mrifile',     {AnatDir, 'FreeSurfer'}, ...
    'nvertices',   15000, ...
    'nas', [135, 222,  75], ...
    'lpa', [ 57, 118,  68], ...
    'rpa', [204, 119,  76], ...
    'ac',  [131, 145, 110], ...
    'pc',  [130, 119, 111], ...
    'ih',  [128, 134, 170]);

% Process: Generate BEM surfaces
bst_process('CallProcess', 'process_generate_bem', [], [], ...
    'subjectname', SubjectName, ...
    'nscalp',      1922, ...
    'nouter',      1922, ...
    'ninner',      1922, ...
    'thickness',   4);

% === LINK CONTINUOUS FILE ===
% Process: Create link to raw file
sFilesRaw = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
    'subjectname',    SubjectName, ...
    'datafile',       {RawFile, 'EEG-DELTAMED'}, ...
    'channelreplace', 1, ...
    'channelalign',   0);

% Process: Set channel file
sFilesRaw = bst_process('CallProcess', 'process_import_channel', sFilesRaw, [], ...
    'channelfile',  {ElcFile, 'XENSOR'}, ...
    'usedefault',   1, ...
    'channelalign', 1);

% Process: Set channels type
sFilesRaw = bst_process('CallProcess', 'process_channel_settype', sFilesRaw, [], ...
    'sensortypes', 'SP1, SP2, RS, PHO, DELR, DELL, QR, QL', ...
    'newtype',     'MISC');

% Process: Project electrodes on scalp
sFilesRaw = bst_process('CallProcess', 'process_channel_project', sFilesRaw, []);

% Process: Events: Import from file
sFilesRaw = bst_process('CallProcess', 'process_evt_import', sFilesRaw, [], ...
    'evtfile', {SpikeFile, 'ARRAY-TIMES'}, ...
    'evtname', 'SPIKE');

% Process: Snapshot: Sensors/MRI registration
bst_process('CallProcess', 'process_snapshot', sFilesRaw, [], ...
    'target',   1, ...  % Sensors/MRI registration
    'modality', 4, ...  % EEG
    'orient',   1, ...  % left
    'comment',  'MEG/MRI Registration');

% === EVALUATION 50 Hz ===
% Process: Power spectrum density (Welch)
sFilesPsd = bst_process('CallProcess', 'process_psd', sFilesRaw, [], ...
    'timewindow',  [], ...
    'win_length',  10, ...
    'win_overlap', 50, ...
    'clusters',    {}, ...
    'isvolume',    0, ...
    'sensortypes', 'EEG', ...
    'edit', struct(...
         'Comment',         'Power', ...
         'TimeBands',       [], ...
         'Freqs',           [], ...
         'ClusterFuncTime', 'none', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0));

% Process: Snapshot: Frequency spectrum
bst_process('CallProcess', 'process_snapshot', sFilesPsd, [], ...
    'target',   10, ...  % Frequency spectrum
    'modality', 4, ...   % EEG
    'comment',  'Power spectrum density');

% === HIGH-PASS FILTER ===
% Process: Import MEG/EEG: Time
sFiles = bst_process('CallProcess', 'process_import_data_time', sFilesRaw, [], ...
    'subjectname', SubjectName, ...
    'condition',   '', ...
    'timewindow',  [], ...
    'split',       0, ...
    'usectfcomp',  1, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);

% Process: High-pass:0.5Hz
sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
    'highpass',    0.5, ...
    'lowpass',     0, ...
    'mirror',      1, ...
    'sensortypes', '', ...
    'overwrite',   1);

% Process: Compute noise covariance
sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
    'baseline', [120, 130], ...
    'dcoffset', 1, ...
    'method',   1, ...  % Full noise covariance matrix
    'copycond', 0, ...
    'copysubj', 0);

% === IMPORT EVENTS ===
% Process: Import MEG/EEG: Events
sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
    'subjectname', SubjectName, ...
    'condition',   '', ...
    'eventname',   'SPIKE', ...
    'timewindow',  [], ...
    'epochtime',   [-0.1, 0.3], ...
    'createcond',  1, ...
    'ignoreshort', 1, ...
    'usectfcomp',  1, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);

% Process: Average: By condition (subject average)
sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype',    3, ...
    'avg_func',   1, ...  % Arithmetic average: mean(x)
    'keepevents', 0);

% Process: Snapshot: Recordings time series
sFiles = bst_process('CallProcess', 'process_snapshot', sFiles, [], ...
    'target',   5, ...  % Recordings time series
    'modality', 4, ...  % EEG
    'comment',  'Average spike');

% Process: Snapshot: Recordings topography (contact sheet)
sFiles = bst_process('CallProcess', 'process_snapshot', sFiles, [], ...
    'target',   7, ...  % Recordings topography (contact sheet)
    'modality', 4, ...  % EEG
    'orient',   1, ...  % left
    'contact_time',   [-40, 110], ...
    'contact_nimage', 16, ...
    'comment',  'Average spike');

% === SOURCE MODELING ===
% Process: Compute head model
sFiles = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
    'sourcespace', 1, ...
    'eeg',         3, ...  % OpenMEEG BEM
    'openmeeg',    struct(...
         'BemSelect', [0, 0, 1], ...
         'BemCond', [1, 0.0125, 1], ...
         'BemNames', {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemFiles', {{}}, ...
         'isAdjoint', 0, ...
         'isAdaptative', 1, ...
         'isSplit', 0, ...
         'SplitLength', 4000));

% Process: Compute sources
sFiles = bst_process('CallProcess', 'process_inverse', sFiles, [], ...
    'method', 1, ...  % Minimum norm estimates (wMNE)
    'wmne', struct(...
         'NoiseCov',      [], ...
         'InverseMethod', 'wmne', ...
         'ChannelTypes',  {{}}, ...
         'SNR',           3, ...
         'diagnoise',     0, ...
         'SourceOrient',  {{'fixed'}}, ...
         'loose',         0.2, ...
         'depth',         1, ...
         'weightexp',     0.5, ...
         'weightlimit',   10, ...
         'regnoise',      1, ...
         'magreg',        0.1, ...
         'gradreg',       0.1, ...
         'eegreg',        0.1, ...
         'ecogreg',       0.1, ...
         'seegreg',       0.1, ...
         'fMRI',          [], ...
         'fMRIthresh',    [], ...
         'fMRIoff',       0.1, ...
         'pca',           1), ...
    'sensortypes', 'EEG', ...
    'output', 1);  % Kernel only: shared

% Process: Snapshot: Sources (one time)
sFiles = bst_process(CallProcess', 'process_snapshot', sFiles, [], ...
    'target',   8, ...  % Sources (one time)
    'modality', 1, ...  % MEG (All)
    'orient',   3, ...  % top
    'time',     0, ...
    'comment', 'Average spike');

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);