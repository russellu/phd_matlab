clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;


for sub=1:length(subs)
    cd(['c:/shared/badger_eeg2/',subs{sub}]) ; ls 
    bcgicas = dir('bcgica*gamma*set') ; 
    for i=1:length(bcgicas) ; EEG = pop_loadset(bcgicas(i).name) ; EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end ; end
    filt = eegfiltfft(merged.data,250,10,15) + eegfiltfft(merged.data,250,40,70)*3 ; 
    filteeg = merged ; filteeg.data = filt ; 
    ep = pop_epoch(filteeg,{'S  1','S  2','S  3'},[-2,7]) ; 
    ica = pop_runica(ep,'runica') ; 
    newica = ica_applyweights(merged,ica) ; 
    ep = pop_epoch(newica,{'S  1','S  2','S  3'},[-2,7]) ;     
    clear ersp 
    for i=1:size(ep.icaact,1)
        [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(i,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                                'plotersp','off','plotitc','off','baseline',0,'timesout',100,'freqs',[1,100],'nfreqs',50,'winsize',64) ; 
    end
    figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-6,6]) ; title(i) ; end ; suptitle(subs{sub}) ; 
    tp(ica) ; suptitle(subs{sub}) ; 
    highfreqs{1} = ica.icaweights ; highfreqs{2} = ica.icasphere ; 
    save('highfreqs','highfreqs') ; 
    
    
end

hcomps = {[32,6,48],[9,18,20],[12,10,4],[24,9,16],[8,15,52],[11,21,33],[46,33,15],[25,16,21],[17,7,21]} ;
