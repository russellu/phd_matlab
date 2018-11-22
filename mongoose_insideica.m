clear all ; close all ; 
subs = {'MONG_01_RB','MONG_02_DP','MONG_03_CG','MONG_04_AB','MONG_05_SG','MONG_06_TS'} ; 
gamma_names = {'gamma_01','gamma_02'} ;
allstim_names = {'allstim_01','allstim_02'} ; 

for s=6:length(subs)
    cd(['c:\shared\mongoose\PROJET_MONGOOSE\raw\',subs{s}]) ; 

    %gammas=dir('*gamma*Pulse*set') ; 
    %allstims=dir('*allstim*Pulse*set') ; 
    rests=dir('*REST*Pulse*vhdr') ; 
    movies=dir('*FIX*Pulse*vhdr') ; 
    %{
    for g=1:length(rests)
        EEG = pop_loadbv('.',allstims(g).name) ; 
        if g==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end
    end
    %}
    merged = pop_loadbv('.',rests(1).name) ; 
     
    merged2 = pop_loadbv('.',movies(1).name) ; 
   
    merged = pop_mergeset(merged,merged2) ; 
    
    
    filt_broad = merged ; filt_broad.data = eegfiltfft(filt_broad.data,merged.srate,1,120) ; broad_ica = pop_runica(filt_broad,'runica') ; 
    filt_narrow = merged ; filt_narrow.data = eegfiltfft(filt_narrow.data,merged.srate,40,60) ; narrow_ica = pop_runica(filt_narrow,'runica') ; 
    
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
    pop_saveset(broad_merged,'broad_merged.set') ; 
end


