clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','tegan','valerie'} ; 
gamma_names = {'gamma_01','gamma_02'} ;
allstim_names = {'allstim_01','allstim_02'} ; 

for s=1:length(subs)
    cd(['C:\shared\badger_eeg\',subs{s}]) ; 

    gammas=dir('*gamma*Pulse*set') ; 
    allstims=dir('*allstim*Pulse*set') ; 
    rests=dir('*rest*Pulse*set') ; 
    movies=dir('*movie*Pulse*set') ; 
    %{
    for g=1:length(rests)
        EEG = pop_loadbv('.',allstims(g).name) ; 
        if g==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    end
    %}
    merged = pop_loadset(rests(1).name) ; 
    morig = merged ; 

    remt = 0 ;    
    for i=1:length(merged.clustlist);
       removetimes = merged.clustlist{i} - remt ; 
       remt = remt + length(merged.clustlist{i}) ; 
       merged = pop_select(merged,'nopoint',[removetimes']) ; 
    end 
    lats = {merged.urevent.latency} ; 
    types = {merged.urevent.type} ; 
    r128s = find(strcmp('R128',types)) ; 
    rlats = cell2mat(lats(r128s)) ; 
    start = rlats(1)+merged.srate ; stop = rlats(length(rlats))-merged.srate ; 
    merged = pop_select(merged,'point',[start,stop]) ; 
    
    
    merged2 = pop_loadset(movies(1).name) ; 
    morig2 = merged2 ; 
    
    remt = 0 ; 
    for i=1:length(merged2.clustlist);
       removetimes = merged2.clustlist{i} - remt ; 
       remt = remt + length(merged2.clustlist{i}) ; 
       merged2 = pop_select(merged2,'nopoint',[removetimes']) ; 
    end 
    lats = {merged2.urevent.latency} ;
    types = {merged2.urevent.type} ; 
    r128s = find(strcmp('R128',types)) ; 
    rlats = cell2mat(lats(r128s)) ; 
    start = rlats(1)+merged2.srate ; stop = rlats(length(rlats))-merged2.srate ; 
    merged2 = pop_select(merged2,'point',[start,stop]) ; 
    
    merged = pop_mergeset(merged,merged2) ; 
    
    trigs = {'S 11','S 12','S 13','S 14','S 21','S 22','S 23','S 24','S 31','S 32','S 33','S 34','S 41','S 42','S 43','S 44','S 51','S 52','S 53','S 54','S 61','S 62','S 63','S 64'} ; 
    
    epochmerged = merged ; %pop_epoch(merged,trigs,[-2,12]) ; 
    filt_broad = epochmerged ; filt_broad.data = eegfiltfft(filt_broad.data,merged.srate,1,120) ; broad_ica = pop_runica(filt_broad,'runica') ; 
    filt_narrow = epochmerged ; filt_narrow.data = eegfiltfft(filt_narrow.data,merged.srate,40,60) ; narrow_ica = pop_runica(filt_narrow,'runica') ; 
    
    broad_merged = ica_applyweights(merged,broad_ica) ;  
    narrow_merged = ica_applyweights(merged,narrow_ica) ; 
    
    %{
    gamma_sets{1} = broad_merged ; gamma_sets{2} = narrow_merged ; 
    stims = trigs ;
    for i=1:length(gamma_sets)
        for stim=1:length(stims)
            epochs = pop_epoch(gamma_sets{i},stims(stim),[-2,7]) ; 
            for c=1:size(epochs.icaact,1) 
                [ersp(i,stim,c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.icaact(c,:,:)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
                    'plotersp','off','plotitc','off','timesout',100,'freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0) ;
            end
        end
    end
    for i=1:2 ; figure ; for j=1:64  ; subplot(5,13,j) ; imagesc(squeeze(mean(ersp(i,1:end,j,:,:),2)),[-8,8]) ; end ; suptitle(subs{s}) ; end 
    %}
    %pop_saveset(broad_merged,'allstim_broad_merged.set') ; 
    %pop_saveset(narrow_merged,'allstim_narrow_merged.set') ;   
    pop_saveset(broad_merged,'motion_merged.set') ; 
end


