clear all ; close all
stimnames{1} = 'unperturbed'  ;
stimnames{2} = 'contrast_5%'  ;
stimnames{3} = 'contrast_33%'  ;
stimnames{4} = 'plaid'  ;
stimnames{5} = 'rnd_10%'  ;
stimnames{6} = 'rnd_60%'  ;
subjects = {
    'alex' %no
    'alexandra3'%yes
    'audrey' %yes
    'charest' %no
    'esteban'%yes
    'fabio'%yes
    'gabriella'%yes
    'gab' %yes
    'genevieve'
    'julie'%yes
    'jeremie'
    'katrine'%yes
    'leila'%no
    'lisa'%yes
    'marc'% yes
    'marie'%yes
    'mathieu'%yes
    'maxime'
    'menglu'%yes
    'mingham'%yes
    'olga'%no
    'patricia' % yes
    'po'%yes
    'russell'%yes
    'suhan2'%yes
    'sunachakan' % yes
    'tah'
    'tegan2'%yes
    'vincent'%no
    'ychele'%yes    
} ;
rootDir = 'c:/shared/allres/' ; 
rtrigs{1} = 'S 11' ;% grating
rtrigs{2} = 'S 12' ;% 5%
rtrigs{3} = 'S 13' ;% 33%
rtrigs{4} = 'S 14' ;% plaid
rtrigs{5} = 'S 15' ;% rnd 10%
rtrigs{6} = 'S 16' ;% rnd 60%
nparams = size(rtrigs,2) ; 
for sub=1:size(subjects,1)
cd([rootDir,subjects{sub}]) ; 
 resamps=dir('resamp*set') ; 
    disp(['PROCESSING SUBJECT: ',subjects{sub}])
    for res=1:size(resamps,1)
        EEG = pop_loadset(resamps(res).name) ; 
        if sub==1 && res==1 ; chanlocs = EEG.chanlocs ; elseif sub==size(subjects,1) ; EEG.chanlocs = chanlocs ; end ; 
        EEG = pop_loadset(resamps(res).name) ; 
        EEG = pop_resample(EEG,256) ; 
        EEG = pop_runica(EEG,'icatype','runica') ; 
        EEG = eeg_checkset(pop_saveset(EEG,['ica_',num2str(res)])) ; 
    end
end